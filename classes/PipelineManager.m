classdef PipelineManager < handle
    % PIPELINEMANAGER manages data processing pipelines.
    % This class allows the creation of analysis pipeline and manages the execution of
    % functions. In addition, it controls for failed / completed  steps in the pipeline.
    
    properties               
        b_ignoreLoggedFiles logical % If true, PIPELINEMANAGER will ignore identical jobs previously run.        
    end
    properties (SetAccess = private)    
        
        % Structure array containing steps of the pipeline. 
        pipe = struct('className', '','argsIn', {},'argsOut',{},'outFileName','',...
            'inputFileName', '','lvl', [], 'b_save2Dat', logical.empty, 'datFileName',...
            '', 'opts',[],'name','');% !!If the fields are changed, please apply 
                                        % them to the set method of this property.
                                        
        fcnDir char % Directory of the analysis functions.        
        funcList struct % structure containing the info about each function in the "fcnDir".
        ProtocolObj Protocol % Protocol Object.
        b_taskState logical % Boolean indicating if the task in pipeline was successful (TRUE) or not (FALSE).
        tmp_LogBook % Temporarily stores the table from PROTOCOL.LOGBOOKFILE
        tmp_BranchPipeline % Temporarily stores LogBook from a Hierarchical branch.
        PipelineSummary % Shows the jobs run in the current Pipeline
        tmp_TargetObj % % Temporarily stores an object (TARGEROBJ).
    end
    properties (Access = private)
        current_task % Task structure currently running.
        current_pipe % Pipeline currently running.
    end
    methods
        % Constructor
        function obj = PipelineManager(ProtocolObj)
            % PIPELINEMANAGER Construct an instance of this class
            %  Input:
            %   ProtocolObj(Protocol): Protocol object. This input is
            %   needed to have access to the protocol's hierarchy of
            %   objects.
            
            p = inputParser;
            addRequired(p,'ProtocolObj', @(x) isa(x, 'Protocol'));
            parse(p,ProtocolObj);                        
            obj.ProtocolObj = p.Results.ProtocolObj;
            root = getenv('Umitoolbox');
            obj.fcnDir = fullfile(root, 'Analysis');
            obj.createFcnList;            
            obj.b_ignoreLoggedFiles = false;
        end
        % SETTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.pipe(obj, pipe)
            % Pipeline structure setter. If pipe is empty, create an empty
            % structure containing tasks fields.
                        
            if isempty(fieldnames(pipe))
                pipe = struct('className', '','argsIn', {},'argsOut',{},'outFileName','',...
                    'inputFileName', '','lvl', [], 'b_save2Dat', logical.empty, 'datFileName',...
                    '', 'opts',[],'name','');
            end
            % Check if all fields exist:
            if ~all(ismember(fieldnames(pipe),fieldnames(obj.pipe)))
                error('umIToolbox:PipelineManager:InvalidInput',...
                    'The pipeline structure provided is invalid!');
            end
            obj.pipe = pipe;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setOpts(obj, func)
            % SETOPTS opens an INPUTDLG for entry of optional variables
            % (OPTS) of methods in the Pipeline.
            
            idx = obj.check_funcName(func);
            if isempty(idx)
                return
            end
            
            S = obj.funcList(idx).info.opts;
            fields = fieldnames(S);
            b_isNum = structfun(@(x) isnumeric(x), S);
            b_isLogic = structfun(@(x) islogical(x), S);
            for i = 1:length(fields)
                if b_isNum(i)
                    fields{i} = [fields{i} ' (numeric)'];
                elseif b_isLogic(i)
                    fields{i} = [fields{i} ' (logical: 0 or 1)'];
                end
            end
            prompt = fields;
            dlgtitle = ['Set optional parameters for ' obj.funcList(idx).name];
            dims = [1 35];
            definput = structfun(@(x) {num2str(x)}, S);
            opts.Resize = 'on';
            answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
            if isempty(answer)
                disp('Operation cancelled by User');
                return
            end
            fields = fieldnames(S);
            for i = 1:length(answer)
                if b_isNum(i)
                    obj.funcList(idx).info.opts.(fields{i}) = str2double(answer{i});
                elseif b_isLogic(i)
                    obj.funcList(idx).info.opts.(fields{i}) = logical(str2double(answer{i}));
                else
                    obj.funcList(idx).info.opts.(fields{i}) = answer{i};
                end
            end
            disp(['Optional Parameters set for function : ' obj.funcList(idx).name]);
        end
        function addTask(obj,className, func, varargin)
            % This method adds an analysis function, or task, to the
            % pipeline. Here, we can choose to save the output of a given
            % task as a .DAT file. ADDTASK will create a string containing 
            % the task that will be evaluated during pipeline execution.
            % Inputs:
            %   className (str): name of the object (Subject, Acquitision,
            %       Modality, etc) that the task will run.
            %   func (str || char):  name or index of the analysis function
            %       contained in obj.funcList property.
            %   b_save2Dat (bool): Optional. True means that the output data 
            %       from the function will be saved as a .DAT file.
            %   datFileName(char): Optional. Name of the file to be saved.
            %       If not provided, the analysis function's default
            %       filename will be used.
            
            p = inputParser;
            addRequired(p,'className', @(x) ischar(x) || isstring(x));
            addRequired(p, 'func', @(x) ischar(x) || isnumeric(x));
            addOptional(p, 'b_save2Dat', false, @islogical);
            addOptional(p, 'datFileName', '', @ischar);
            parse(p,className, func, varargin{:});
            
            % Check if the function name is valid:
            idx = obj.check_funcName(p.Results.func);
            if isempty(idx)
                warning('Operation cancelled! The function "%s" does not exist!',...
                    p.Results.func);
                return
            end
            % Create "task" structure. This is the one that will be added
            % to the pipeline:
            task = obj.funcList(idx).info;
            task.inputFileName = '';
            task.className = p.Results.className;
            task.name = obj.funcList(idx).name;
            task.b_save2Dat = p.Results.b_save2Dat;
            task.datFileName = p.Results.datFileName;
            % Determine Level of the task in the Protocol's hierarchy.
            switch p.Results.className
                %                 case 'Protocol' % Disabled for now.
                %                     lvl = 4;
                case 'Subject'
                    task.lvl = 3;
                case 'Acquisition'
                    task.lvl = 2;
                otherwise
                    task.lvl = 1;
            end
            % Control for steps already in the pipeline:
            if any(strcmp(task.name, {obj.pipe.name}))
                warning('Operation cancelled! The function "%s" already exists in the Pipeline!',....
                    task.name);
                return
            end
            
            % Case where this task is the 1st step of the pipeline:
            if isempty(obj.pipe)
                task.inputFileName = obj.getFirstInputFile(task);                
            end
            % Control for multiple outputs from the previous step:
            % Here, we assume that functions with multiple outputs
            % create only "Files" and not "data".
            % Therefore, we will update the function string to load
            % one of the "Files" before running the task.
            
            % Look from bottom to top of the pipeline for tasks with files
            % as outputs. This is necessary because not all analysis
            % functions have outputs.
            
            if any(strcmp(task.argsIn, 'data'))
                 for idxOutFile = length(obj.pipe):-1:1
                     if any(strcmp('outFile', obj.pipe(idxOutFile).argsOut))
                         break
                     end
                 end
                 if iscell(obj.pipe(idxOutFile).outFileName)
                     % Ask user to select a file:
                     disp('Controlling for multiple outputs')
                     w = warndlg({'Previous step has multiple output files!',...
                         'Please, select one to be analysed!'});
                     waitfor(w);
                     [indxFile, tf] = listdlg('ListString', obj.pipe(idxOutFile).outFileName,...
                         'SelectionMode','single');
                     if ~tf
                         disp('Operation cancelled by User')
                         return
                     end
                     task.inputFileName = obj.pipe(idxOutFile).outFileName{indxFile};
                 else
                     task.inputFileName = obj.pipe(idxOutFile).outFileName;
                 end
            end
            
            % Save to Pipeline:
            obj.pipe = [obj.pipe; task];
            disp(['Added "' task.name '" to pipeline.']);
           
            % Control for data to be saved as .DAT files for task:
            if ~task.b_save2Dat
                return
            end
            if ~any(strcmp('outData', task.argsOut)) 
                warning(['Cannot save output to .DAT file for the function'...
                    ' "%s" \nbecause it doesn''t have any data as output!'], task.name);
                return
            end
            % Save datFileName as default output name from task's function:
            if isempty(task.datFileName)
                obj.pipe(end).datFileName = obj.pipe(end).outFileName;
            % OR update datFileName to add file extension:
            else
                obj.pipe(end).datFileName = [obj.pipe(end).datFileName, '.dat'];
            end
            
                
        end
        function varargout = showPipeSummary(obj)
            % This method creates a summary of the current pipeline.
            
            % Output (optional): if an output variable exists, it creates a
            % character array, if not, the information is displayed in the
            % command window.
            
            if isempty(obj.pipe)
                disp('Pipeline is empty!')
                if nargout == 1
                    varargout{1} = '';
                end
                return
            end
            
            str = sprintf('%s\n', 'Pipeline Summary:');
            for i = 1:length(obj.pipe)
                str =  [str sprintf('--->> Step # %d <<---\n', i)];
                if isempty(obj.pipe(i).opts)
                    opts = {'none'; 'none'};
                else
                    opts = [fieldnames(obj.pipe(i).opts)';...
                        cellfun(@(x) num2str(x), struct2cell(obj.pipe(i).opts), 'UniformOutput', false)'];
                end
                txt = sprintf('Function name : %s\nOptional Parameters:\n',...
                    obj.pipe(i).name);
                str = [str, txt, sprintf('\t%s : %s\n', opts{:})];
                if obj.pipe(i).b_save2Dat
                    str = [str, sprintf('Data to be saved as : "%s"\n', obj.pipe(i).datFileName)];
                end
                if ~isempty(obj.pipe(i).inputFileName)
                    str = [str, sprintf('Input File Name : "%s"\n', obj.pipe(i).inputFileName)];
                end
                str = [str, sprintf('--------------------\n')];
            end
            if nargout == 0
                disp(str)
            else
                varargout{1} = str;
            end
        end        
        function showFuncList(obj)
            % Displays a list of analysis function from "obj.funcList" in
            % the command window.
            disp('List of available functions (index : Name) :');
            for i = 1:length(obj.funcList)
                fprintf('%d : %s\n', i, obj.funcList(i).name);
            end
        end
        
        
        % TO BE CHANGED... %%
        function run_pipeline(obj)
            % RUN_PIPELINE runs the tasks in OBJ.PIPE
            lbf = matfile(obj.ProtocolObj.LogBookFile);
            obj.tmp_LogBook = lbf.LogBook;
            obj.PipelineSummary = obj.ProtocolObj.createEmptyTable;
            % Identify the maximum level in the hierarchy to run the
            % pipeline.
            top_lvl = max([obj.Pipe.level]);
            % Get indexes of Filtered Objects from OBJ.PROTOCOLOBJ.QUERYFILTER function.
            idxList = obj.ProtocolObj.Idx_Filtered;
            % Identify branches in the hierarchy that will be processed in
            % a single pipeline run.
            ppIdx = zeros(size(idxList,1),1);
            switch top_lvl
                % The pipeline will run at the modality level.
                case 1
                    ppIdx = (1:size(idxList,1))';
                    % The pipeline will run at the acquisition level.
                case 2
                    a = 1;
                    uniqA = unique(idxList(:,[1 2]),'rows');
                    for i = 1:size(uniqA,1)
                        idx = all(idxList(:,[1 2]) == uniqA(i,:),2);
                        ppIdx(idx) = a;
                        a = a+1;
                    end
                    % The pipeline will run at the subject level.
                case 3
                    a = 1;
                    uniqS = unique(idxList(:,1),'rows');
                    for i = 1:length(uniqS)
                        idx = all(idxList(:,1) == uniqS(i,:),2);
                        ppIdx(idx) = a;
                        a = a+1;
                    end
            end
            f = waitbar(0,'Analysing data...', 'Name', 'Pipeline progress');
            % run pipeline at each branch of the hierarchy
            uniq_branch = unique(ppIdx);
            for i = 1:length(uniq_branch)
%                 obj.State = true;
                waitbar(i/length(uniq_branch),f);
                branch = idxList(ppIdx == uniq_branch(i),:);
                obj.run_tasksOnBranch(branch);
%                 if ~obj.State
%                     continue
%                 end
            end
            
            obj.PipelineSummary(1,:) = []; % remove dummy row.
            
            LogBook = obj.PipelineSummary;
            save(obj.ProtocolObj.LogBookFile, 'LogBook');
            cd(obj.ProtocolObj.SaveDir);
            
            disp(obj.PipelineSummary)
            % Save Protocol Object:
            protocol = obj.ProtocolObj;
            save([obj.ProtocolObj.SaveDir obj.ProtocolObj.Name '.mat'], 'protocol');
            disp('Protocol object Saved!');
            waitbar(1,f,'Finished!');
            delete(f)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Pipeline Management methods:       
        function savePipe(obj, filename)
            % SAVEPIPE saves the structure OBJ.PIPE in a .MAT file in the
            % folder PIPELINECONFIGFiles inside MAINDIR of OBJ.PROTOCOLOBJ.            
            targetDir = fullfile(obj.ProtocolObj.SaveDir, 'PipeLineConfigFiles');
            [~,~] = mkdir(targetDir);            
            pipeStruct = obj.pipe;
            txt = jsonencode(pipeStruct);
            fid = fopen(fullfile(targetDir,[filename '.json']), 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            disp(['Pipeline saved as "' filename '" in ' targetDir]);               
        end
        function loadPipe(obj, filename)
            % LOADPIPE loads the structure PIPE inside FILENAME and assigns
            % it to OBJ.PIPE property.
            if ~isfile(filename)                
                filename = fullfile(obj.ProtocolObj.SaveDir, 'PipeLineConfigFiles',...
                    [filename, '.json']);
            end
            txt = fileread(filename);
            obj.pipe = jsondecode(txt);
            disp('Pipeline Loaded!');
            obj.showPipeSummary;
        end        
        function reset_pipe(obj)
            % This function erases the pipe property and resets the funcList
            % property to default parameter values.            
            obj.pipe = struct();
            obj.funcList = struct.empty;
            obj.createFcnList;
        end
    end
    
    methods (Access = private)        
        
        function subTasks = pipeSplitter(obj)
            % PIPESPLITTER segments the pipeline in sub-pipelines that are run
            % at each level of the Hierarchy.
            
            % Find consecutive levels
            lvls = [obj.Pipe.level];
            idx = ones(1,length(lvls));
            a = 1;
            for i = 1:length(lvls)-1
                lvl = lvls(i);
                next_lvl = lvls(i+1);
                if lvl ~= next_lvl
                    idx(i) = a;
                    a = a+1;
                    idx(i+1) = a;
                else
                    idx(i) = a;
                    idx(i+1) = a;
                end
            end
            uniq_idx = unique(idx);
            subTasks = cell(1,numel(uniq_idx));
            for i = 1:numel(uniq_idx)
                b_idx = ( idx == uniq_idx(i) );
                subTasks{i} = obj.pipe(b_idx);
            end
        end
        function run_tasksOnBranch(obj, branch)
            % RUN_TASKSONBRANCH finds the object to run the tasks in the
            % pipeline.
            
            % Split pipeline if there is more than one level.
            ppLine = obj.pipeSplitter;
            obj.tmp_Branchpipeline = obj.ProtocolObj.createEmptyTable;
            for i = 1:length(ppLine)
                subtasks = ppLine{i};
                obj.current_pipe = subtasks;
                lvl = subtasks.level;
                switch lvl
                    case 1
                        targetIdxArr = unique(branch,'rows');
                    case 2
                        targetIdxArr = unique(branch(:,[1 2]), 'rows');
                    case 3
                        targetIdxArr = unique(branch(:,1), 'rows');
                end
                
                for j = 1:size(targetIdxArr,1)
                    obj.getTargetObj(targetIdxArr(j,:));
                    obj.tmp_TargetObj.LastLog = obj.ProtocolObj.createEmptyTable;
                    obj.readFilePtr;
                    for k = 1:length(subtasks)
                        obj.current_task = subtasks(k);
                        obj.run_taskOnTarget;
                        if ~obj.State
                            return
                        end
                    end
                    obj.tmp_Branchpipeline = [obj.tmp_Branchpipeline; obj.tmp_TargetObj.LastLog];
                end
            end
            obj.eraseIntermediateFiles();
        end
        function run_taskOnTarget(obj)
            % RUN_TASKONTARGET runs a task in the pipeline structure array in
            % TASK.
            %   It checks if the command TASK.FUNCNAME was already
            %   sucessfully performed by comparing it to the LOGBOOK from
            %   the PROTOCOL object. Also, it appends the information of
            %   the task on the object's LASTLOG.
            
            cd(obj.tmp_TargetObj.SaveFolder);
            LastLog = obj.ProtocolObj.createEmptyTable;
            % Fill out LOGTABLE:
            % Parse TARGETOBJ.SAVEFOLDER to get Subject/Acquisition/Recording
            % names:
            str = erase(obj.tmp_TargetObj.SaveFolder, obj.ProtocolObj.SaveDir);
            str = split(str, filesep); str = str(cellfun(@(x) ~isempty(x), str));
            for i = 1:length(str)
                LastLog(:,i) = str(i);
            end
            task = obj.current_task;
            task.InputFile_UUID = 'None';
            task = populateFuncStr(obj, task);
            if ~strcmp(task.InputFile_UUID, 'None')
                LastLog.InputFile_UUID = {task.InputFile_UUID};
                LastLog.InputFile_Path = {task.Input};
            end
            LastLog.ClassName = {class(obj.tmp_TargetObj)};
            LastLog.Job = {task.funcStr};
            % Check if Job was already performed:
            b_isLogged = obj.checkInFilePtr(task);
            if ~b_isLogged || (b_isLogged && obj.IgnoreLoggedFiles)
                % Run the step:
                try
                    if strcmp(task.Input, 'missing')
                        errID = 'MATLAB:Umitoolbox:pipelineManager:FileNotFound';
                        errmsg = ['Input File for function ' task.Name ' not found!'];
                        error(errID,errmsg);
                    end
                    disp(['Running ' task.Name '...']);
                    % Load options structure in the workspace.
                    opts = task.opts;
                    % Evaluate function string:
                    eval(task.funcStr);
                    %
                    state = true;
                    LastLog.Messages = 'No Errors';
                catch ME
                    state = false;
                    LastLog.Messages = {getReport(ME)};
                    LastLog.Messages_short = {getReport(ME, 'basic','hyperlinks','off')};
                end
                LastLog.Completed = state;
                LastLog.RunDateTime = datetime('now');
                obj.tmp_TargetObj.LastLog = [obj.tmp_TargetObj.LastLog; LastLog];
                obj.pipelineSummary = [obj.pipelineSummary; LastLog];
                obj.State = state;
                if LastLog.Completed
                    disp('Task Completed!')
                    if exist('out', 'var')
                        if ischar(out)
                            out = {out};
                        end
                        for i = 1:length(out)
                            SaveFolder = task.SaveIn;
                            if endsWith(out{i}, '.dat')
                                mDt_file = matfile(strrep(fullfile(SaveFolder, out{i}), '.dat', '_info.mat'),'Writable', true);
                            else
                                mDt_file = matfile(out{i},'Writable', true);
                            end
                            % Inheritance of MetaData from last File created (Different ones ONLY by different function).
                            lastFile = task.Input;
                            if isfile(lastFile)
                                lastMetaData = matfile(strrep(lastFile, '.dat', '_info.mat'));
                                props = setdiff(properties(lastMetaData), properties(mDt_file));
                                for k = 1:length(props)
                                    eval(['mDt_file.' props{k} '= lastMetaData.' props{k} ';'])
                                end
                            end
                            fileUUID = mDt_file.fileUUID;
                            if iscell(fileUUID)
                                fileUUID = [fileUUID{:}];
                            end
                            task.File_UUID = fileUUID;
                            task.FileName = out{i};
                            obj.current_task = task;
                            obj.write2FilePtr(task);
                        end
                    end
                else
                    disp('Failed!')
                end
            else
                disp([task.Name ' Skipped!'])
                obj.State = true;
                return
            end
        end
        function getTargetObj(obj, targetIdx)
            % GETTARGETOBJ finds the object TARGETOBJ indicated by the
            % index TARGETIDX inside PROTOCOL.
            
            tgtSz = size(targetIdx,2);
            switch tgtSz
                case 1
                    targetObj = obj.ProtocolObj.Array.ObjList(targetIdx);
                case 2
                    targetObj = obj.ProtocolObj.Array.ObjList(targetIdx(1)).Array.ObjList(targetIdx(2));
                case 3
                    targetObj = obj.ProtocolObj.Array.ObjList(targetIdx(1)).Array.ObjList(targetIdx(2)).Array.ObjList(targetIdx(3));
            end
            obj.tmp_TargetObj = targetObj;
        end        
        
        function out = getFirstInputFile(obj, funcInfo)
            % This method verifies is the function has any data as input.
            % If yes, then a creates an dialog box containing a list of .DAT files
            % to choose as input. This method is called by ADDTASK only
            % when the first task of a pipeline is created.
            % Input:
            %   funcInfo (struct): structure containing the task's function info.
            % Output:
            %   out (char): name of input file. Empty if file does not exist.
            
            % Initialize output structure:
            out = '';
                        
            % Control for tasks that do not have any data as input:
            if ~any(strcmp('data', funcInfo.argsIn))
                return
            end
            
            % Get all existing objects from the selected items in Protocol
            % object:
            idx = obj.ProtocolObj.Idx_Filtered;
            idx = unique(idx,'rows');
            classes = {};
            for i = 1:size(idx,1)
                classes{i} = class(obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3)));
            end
            classes = [{'Subject', 'Acquisition'}  classes]; classes = unique(classes);
            % Select Object containing input file:
            [indx,tf] = listdlg('PromptString', {'Select the Object containing', 'the input for the function :',  FuncName},'ListString',classes, 'SelectionMode', 'single');
            if ~tf
                return
            end
            switch classes{indx}
                case 'Subject'
                    targetObj = obj.ProtocolObj.Array.ObjList(idx(1,1));
                case 'Acquisition'
                    targetObj = obj.ProtocolObj.Array.ObjList(idx(1,1)).Array.ObjList(idx(1,2));
                otherwise
                    for i = 1:size(idx,1)
                        b_isMod = strcmp(class(obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3))), classes{indx});
                        if b_isMod
                            break
                        end
                    end
                    targetObj = obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3));
            end
            % Display a list of .DAT files from the selected
            % object(targetObj):
            FileList = dir(fullfile(targetObj.SaveFolder, '*.dat'));
            if isempty(FileList)
                warndlg(['No Files found in ' targetObj.SaveFolder], 'pipeline warning!', 'modal')
                return
            else
                [indx,tf] = listdlg('PromptString', {'Select the File from ' classes{indx},...
                    'as input for the function :',  FuncName},'ListString',{FileList.name}, 'SelectionMode',...
                    'single');
                if ~tf
                    return
                end
            end
            
            out = FileList{indx}.name;
            
        end
        
%         function readFilePtr(obj)
%             % READFILEPTR loads the content of FILEPTR.JSON in a structure
%             % stored in OBJ.TMP_FILEPTR.
%             txt = fileread(obj.tmp_TargetObj.FilePtr);
%             a = jsondecode(txt);
%             for i = 1:numel(a.Files)
%                 a.Files(i).Folder = tokenizePath(a.Files(i).Folder, obj.tmp_TargetObj, 'detokenize');
%                 a.Files(i).InputFile_Path = tokenizePath(a.Files(i).InputFile_Path, obj.tmp_TargetObj, 'detokenize');
%             end
%                 obj.tmp_FilePtr = a;
%         end
%         function write2FilePtr(obj, task)
%             % WRITE2FILEPTR writes the File information stored in structure FILEINFO
%             % in OBJ.TMP_TARGETOBJ.FILEPTR.
%             
%             %Initialize
%             FileInfo = struct('Name', task.FileName, 'UUID', task.File_UUID, 'Folder', task.SaveIn, 'InputFile_Path', task.Input,...
%                 'InputFile_UUID', task.InputFile_UUID, 'creationDateTime', datestr(now), 'FunctionInfo', ...
%                 struct('Name', task.Name, 'DateNum', task.DateNum, 'Job', task.funcStr, 'opts', task.opts));
%             
%             FileList = obj.tmp_FilePtr.Files;
%             % Check for Files already logged on FilePtr
%             idx = false(length(FileList),2);
%             for i = 1:length(FileList)
%                 idx(i,1) = strcmp(FileInfo.Name, FileList(i).Name);
%                 idx(i,2) = strcmp(FileInfo.FunctionInfo.Name, FileList(i).FunctionInfo.Name);
%             end
%             idx = all(idx,2);
%             % If there are no Files 
%             if isempty(FileList)
%                 obj.tmp_FilePtr.Files = FileInfo;
%             % If there are files and one identical, replace it.
%             elseif ~isempty(FileList) && any(idx)
%                 obj.tmp_FilePtr.Files(idx) = FileInfo;
%             % If there are files and none identical: Append
%             else
%                 obj.tmp_FilePtr.Files = [FileList; FileInfo];
%             end
%             for i = 1:numel(obj.tmp_FilePtr.Files)
%                 obj.tmp_FilePtr.Files(i).Folder = tokenizePath(obj.tmp_FilePtr.Files(i).Folder, obj.tmp_TargetObj);
%                 obj.tmp_FilePtr.Files(i).InputFile_Path = tokenizePath(obj.tmp_FilePtr.Files(i).InputFile_Path, obj.tmp_TargetObj);
%             end
%             txt = jsonencode(obj.tmp_FilePtr);
%             fid = fopen(obj.tmp_TargetObj.FilePtr, 'w');
%             fprintf(fid, '%s', txt);
%             fclose(fid);
%         end
        function createFcnList(obj)
            % This function creates a structure containing all information
            % about the analysis functions inside the "Analysis" folder.
            % This information is stored in the "funcList" property of
            % pipelineManager.
            
            % !!For now, it will read only folders directly below the
            % "Analysis" folder. Subfolders inside these folders will not
            % be read.
            
            % Set Defaults:
            default_Output = '';
            default_opts = struct();
            
            disp('Creating Fcn list...');
            list = dir(fullfile(obj.fcnDir, '\*\*.m'));
            for i = 1:length(list)
                out = parseFuncFile(list(i));
                % Validate if all input arguments from the function are
                % "valid" inputs keywords:
                kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts', 'object'};
                kwrds_out = {'outFile', 'outData', 'metaData', 'outDataStat'};
                if all(ismember(out.argsIn, kwrds_args)) && all(ismember(out.argsOut, kwrds_out))
%                     disp(list(i).name);
                    [~,list(i).name, ~] = fileparts(list(i).name);
                    list(i).info = out;
                    
                    obj.funcList = [obj.funcList ; list(i)];
                end
                
            end
            disp('Function list created!');
            function info = parseFuncFile(fcnStruct)
                % This function reads the code inside a function listed in
                % "fcnStruct" to retrieve input arguments, outputs and
                % default options. 
                % Input:
                %   fcnStruct (struct): structure containing the file info
                %   from a given analysis function.
                % Output:
                %   info (struct): structure containing list of arguments
                %   (inputs and outputs) and optional parameters from the
                %   analysis function.
                
                info = struct('argsIn', {},'argsOut', {}, 'outFileName', '', 'opts', []);
                txt = fileread(fullfile(fcnStruct.folder, fcnStruct.name));
                funcStr = regexp(txt, '(?<=function\s*).*?(?=\n)', 'match', 'once');
                funcStr = erase(funcStr(1:end-1), ' '); % remove white spaces and an extra line from the function string.
                outStr = regexp(funcStr,'.*(?=\=)', 'match', 'once');
                out_args = regexp(outStr, '\[*(\w*)\,*(\w*)\]*', 'tokens', 'once');
                idx_empty = cellfun(@isempty, out_args);
                info(1).argsOut = out_args(~idx_empty);
                [~,funcName,~] = fileparts(fcnStruct.name);
                expInput = ['(?<=' funcName '\s*\().*?(?=\))'];
                str = regexp(funcStr, expInput, 'match', 'once');
                str = strip(split(str, ','));
                idx_varargin = strcmp(str, 'varargin');
                info.argsIn = str(~idx_varargin);
                expOutput = 'default_Output\s*=.*?(?=\n)';
                str = regexp(txt, expOutput, 'match', 'once');
                if isempty(str)
                    default_Output = '';
                else
                    eval(str)
                end
                info.outFileName = default_Output;
                expOpts = 'default_opts\s*=.*?(?=\n)';
                str = regexp(txt, expOpts, 'match', 'once');
                if ~isempty(str)
                    eval(str)
                    info.opts = default_opts;
                    info.argsIn{end+1} = 'opts';
                end
            end
        end
        function task = populateFuncStr(obj, task)
            % POPULATEFUNCSTR replaces keywords in TASK.FUNCSTR with the
            % info contained in OBJ.TMP_TARGETOBJ. It is used in
            % OBJ.RUN_TASKONTARGET.

            % Replace Input string:
            switch task.Input
                case 'RawFolder'
                    folder = obj.tmp_TargetObj.RawFolder;
                    task.Input = folder; 
                case 'SaveFolder'
                    folder = obj.tmp_TargetObj.SaveFolder;
                    task.Input = folder;
                case 'object'
                    task.Input = 'obj.tmp_TargetObj';
                     otherwise
                    % ALL THIS SECTION NEEDS TO BE CHANGED!!
                    if isempty(obj.tmp_FilePtr.Files)
                        task.Input = 'missing';
                        task.funcStr = '';
                        return
                    else
                        idx = strcmp(task.Input, {obj.tmp_FilePtr.Files.Name});
                        if sum(idx) == 0
                            try
                                % Try to find file with name different from
                                % default:
                                idx_pipe = strcmp(task.Input, {obj.current_pipe.Output});
                                prev_task = obj.current_pipe(idx_pipe);
                                % Find file in FilePtr from function in prev_task:
                                fcn_info = arrayfun(@(x) x.FunctionInfo, obj.tmp_FilePtr.Files);
                                idxFcnName = strcmp({fcn_info.Name}, prev_task.Name);
                                idxFcnDate = [fcn_info.DateNum] == prev_task.DateNum;
                                idx = idxFcnName & idxFcnDate;
                                if sum(idx) == 1
                                    filePath = fullfile(obj.tmp_TargetObj.SaveFolder, obj.tmp_FilePtr.Files(idx).Name);
                                else
                                    task.Input = 'missing';
                                    task.funcStr = '';
                                    return
                                end
                            catch
                                task.Input = 'missing';
                                task.funcStr = '';
                                return
                            end
                        end
                        if strcmp(obj.tmp_FilePtr.Files(idx).Folder, 'RawFolder')
                            filePath = fullfile(obj.tmp_TargetObj.RawFolder, obj.tmp_FilePtr.Files(idx).Name);
                        else
                            filePath = fullfile(obj.tmp_TargetObj.SaveFolder, obj.tmp_FilePtr.Files(idx).Name);
                        end
                    end
                    task.Input = filePath;
                    inputMetaData = strrep(filePath, '.dat', '_info.mat');
                    mDt_input = matfile(inputMetaData); fileUUID = mDt_input.fileUUID;
                    if iscell(fileUUID)
                        fileUUID = [fileUUID{:}];
                    end
                    task.InputFile_UUID = fileUUID;
            end
            % Replace SaveFolder string:
            switch task.SaveIn
                case 'RawFolder'
                    folder = obj.tmp_TargetObj.RawFolder;
                    task.SaveIn = folder;
                case 'SaveFolder'
                    folder = obj.tmp_TargetObj.SaveFolder;
                    task.SaveIn = folder;                   
            end
            if ~strcmp(task.Input, 'obj.tmp_TargetObj')
                funcStr = [task.Name '(''' task.Input ''',''' task.SaveIn ''''];
            else
                funcStr = [task.Name '(' task.Input ',''' task.SaveIn ''''];
            end
            % Fix empty input character "~":
            funcStr = strrep(funcStr, '''~''', '~');
                
            % Add optionals: 
            if ~isempty(task.opts)  
                funcStr = [funcStr ', opts'];
            end
            if ~isempty(task.Output)
                funcStr = ['out = ' funcStr];
            end
            funcStr = [funcStr ');'];
            task.funcStr = funcStr;
        end
        %%%%%%%%%% NEW METHODS %%%%%%%%
        
         function idx_fcn = check_funcName(obj, func)
            % This function is used by "setOpts" and "addTask" methods to
            % validate if the input "func" is valid.
            % Input
            %   func (numeric OR char) : index OR name of a function from
            %   obj.funcList.
            % Output:
            %   idx_fcn(bool): index of "func" in "obj.funcList". Returns
            %   empty if the function was not found in the list.
            idx_fcn = [];
            if isnumeric(func)
                idx_fcn = func == 1:length(obj.funcList);
                msg = ['Function with index # ' num2str(func)];
            else
                idx_fcn = strcmp(func, {obj.funcList.name});
                msg = ['Function "' func '"'];
            end
            if ~any(idx_fcn)
                disp([msg ' not found in the function list!']);
                idx_fcn = [];
            end
        end
    end
end

