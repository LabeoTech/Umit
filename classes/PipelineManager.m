classdef PipelineManager < handle
    % PIPELINEMANAGER creates and manages pre-processing and
    % analysis pipelines.
    %  Creates a structure with the info necessary to run a simple pipeline
    %  as a sequence of steps(i.e. step 1, step 2, step3, ... step N).
    %  Controls for failed / completed  steps in the pipeline.
    properties
        EraseIntermediate % Boolean to indicate if PIPELINEMANAGER will delete intermediate files in the pipeline. It keeps only the first and last file in the pipeline.
        IgnoreLoggedFiles % Boolean. If true, PIPELINEMANAGER will ignore identical jobs previously run (logged in OBJ.PROTOCOLOBJ.LOGBOOKFILE).
        FuncRootDir % Directory where all functions are located.
        Pipe % Array containing steps of the pipeline.
    end
    properties (SetAccess = private, SetObservable)
        ProtocolObj % Protocol Object.
        State % Boolean indicating if the task in pipeline was successful (TRUE) or not (FALSE).
        tmp_LogBook % Temporarily stores the table in PROTOCOL.LOGBOOKFILE
        tmp_BranchPipeline % Temporarily stores LogBook from a Hierarchical branch.
        PipelineSummary % Shows the jobs run in the current Pipeline
        FunctionList % Structure containing the Public methods of all Classes in the ISAtoolbox with its inputs, optional parameters and outputs.
        tmp_FilePtr % Temporarily stores content of FILEPTR.JSON file from TARGETOBJ.
        tmp_TargetObj % % Temporarily stores an object (TARGEROBJ).
        block_pipe = false % Boolean that, if true, the pipeline won't accept any other steps.
    end
    properties (Access = private)
        current_task % Task structure currently running.
        current_pipe % Pipeline currently running.
    end
    methods
        % Constructor
        function obj = PipelineManager(Pipe, ProtocolObj, FuncRootDir)
            % PIPELINEMANAGER Construct an instance of this class
            %   PIPE is a structure containing all the information of the pipeline.
            if nargin > 0
                obj.Pipe = Pipe;
                obj.ProtocolObj = ProtocolObj;
                obj.FuncRootDir = FuncRootDir;
            end
            obj.createFunctionList;
            obj.EraseIntermediate = false;
            obj.IgnoreLoggedFiles = false;
        end
        % Set functions
        function set.ProtocolObj(obj, ProtocolObj)
            validateattributes(ProtocolObj, {'Protocol'}, {'nonempty'});
            obj.ProtocolObj = ProtocolObj;
        end
        function set.FuncRootDir(obj, FuncRootDir)
            obj.validate_path(FuncRootDir)
            addpath(genpath(FuncRootDir));
            obj.FuncRootDir = FuncRootDir;
        end
        %%% Validators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function validate_path(~, input)
            if ~isfolder(input)
                errID = 'IsaToolbox:InvalidInput';
                msg = 'Input is not a valid folder or it is not in MATLAB''s path.';
                error(errID, msg);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function addTask(obj,className, funcName, varargin)
            p = inputParser;
            addRequired(p,'className', @(x) ischar(x) || isstring(x));
            addRequired(p,'funcName', @(x) ischar(x) || isstring(x));
            addOptional(p, 'params', struct, @isstruct);
            parse(p,className, funcName, varargin{:});

            idx = strcmp(p.Results.funcName, {obj.FunctionList.Name});
            funcInfo = obj.FunctionList(idx);
            task = obj.validateParams(p.Results.params,funcInfo);
            if obj.block_pipe
               disp('No more steps allowed. Last step has no output.')
                return
            end
            if isempty(task.Output)
                warndlg(['The function ' task.Name ' has no outputs. '...
                    'No more steps can be added to the pipeline after this one.'],...
                    'Pipeline warning')
                obj.block_pipe = true;
            end
                
            % Determine Level of the step in the hierarchy.
            switch p.Results.className
                %                 case 'Protocol' % Disabled for now.
                %                     lvl = 4;
                case 'Subject'
                    lvl = 3;
                case 'Acquisition'
                    lvl = 2;
                otherwise
                    lvl = 1;
            end
            % Adds a step to the pipeline:
            b_valid = obj.validateTask(p.Results.className, task);
            if ~b_valid
                return
            end
            % Determine the task INPUT:                        
            b_inputIsFile = strcmp(task.Input, 'File');
            if isempty(obj.Pipe) && b_inputIsFile
                task.Input = obj.askForFirstInput(p.Results.funcName);
                if isempty(task.Input)
                    disp('Operation cancelled by user!')
                    return
                end
            elseif ~isempty(obj.Pipe) && b_inputIsFile
                if ~contains(obj.Pipe(end).Output, '.dat')
                    prevStep = obj.Pipe(end);
                    errordlg({['Error adding ' task.Name ' to pipeline.'],...
                        ['The output of the previous step "' prevStep.Name ...
                        '" (' prevStep.Output ') is incompatible with the current function input (' task.Input ').'],...
                        'Try rearranging the pipeline to fix this.'}, 'Pipeline order error')
                    return
                end    
                task.Input = obj.Pipe(end).Output;
            end

            % Choose input if multiple options are available:
            if iscell(task.Input) && numel(task.Input) > 1
                task.Input = obj.pickAnInput(task);
            end

            % Extra info in task structure:
            task.className = p.Results.className;
            task.level = lvl;
            obj.Pipe = [obj.Pipe task];
        end
        function opts = setOpts(obj, funcName)
            % SETOPTS opens an INPUTDLG for entry of optional variables
            % (OPTS) of methods in the Pipeline. Output: Structure
            % containing the variables and values.
            idx = strcmp(funcName, {obj.FunctionList.Name});
            S = obj.FunctionList(idx).opts;
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
            dlgtitle = ['Set optional parameters for ' funcName];
            dims = [1 35];
            definput = structfun(@(x) {num2str(x)}, S);
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            fields = fieldnames(S);
            opts = struct;
            for i = 1:length(answer)
                if b_isNum(i)
                    opts.(fields{i}) = str2double(answer{i});
                elseif b_isLogic(i)
                    opts.(fields{i}) = logical(str2double(answer{i}));
                else
                    opts.(fields{i}) = answer{i};
                end
            end
        end        
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
        function showPipeSummary(obj)
            % SHOWPIPESUMMARY displays a summary of the pipeline steps in
            % the command window.
            
            disp(repmat('-', 2,100));
            fields = setdiff(fieldnames(obj.Pipe),{'opts', 'Output'});
            for i = 1:length(obj.Pipe)
                disp(['Step # ' num2str(i)])
                for j = 1:length(fields)
                   fprintf('%s : %s \n', fields{j}, num2str(obj.Pipe(i).(fields{j})))
                end
                disp('##Optional Parameters##')
                if ~isempty(obj.Pipe(i).opts)
                Optsfields = fieldnames(obj.Pipe(i).opts);
                for j = 1:length(Optsfields)
                   fprintf('%s : %s \n', Optsfields{j}, num2str(obj.Pipe(i).opts.(Optsfields{j})))
                end
                else
                    disp(' No custom parameters selected.')
                end
                disp(repmat('-', 1,100));           
            end
            disp(repmat('-', 1,100));           
        end
        function savePipe(obj, filename)
            % SAVEPIPE saves the structure OBJ.PIPE in a .MAT file in the
            % folder PIPELINECONFIGFiles inside MAINDIR of OBJ.PROTOCOLOBJ.
            savedir = obj.ProtocolObj.SaveDir;
            cd(savedir)
            [~,~] = mkdir('PipeLineConfigFiles');
            cd('PipeLineConfigFiles');
            pipeStruct = obj.Pipe;
            txt = jsonencode(pipeStruct);
            fid = fopen([filename '.json'], 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
        end
        function loadPipe(obj, filename)
            % LOADPIPE loads the structure PIPE inside FILENAME and assigns
            % it to OBJ.PIPE property.
            if ~isfile(filename)
                path = which([obj.MainDir PipeLineConfigFiles]);
                filename = fullfile(path, filename);
            end
            txt = fileread(filename);
            a = jsondecode(txt);
            obj.Pipe = a;
        end
        
    end
    
    methods (Access = private)
        function newParams = validateParams(~, params,funcInfo)
            % VALIDATEPARAMS checks for user-input parameters in PARAMS and
            % changes the FUNCINFO fields accordingly. Outputs NEWPARAMS with
            % the updated info.
            newParams = funcInfo; newParams.opts = '';
            if isempty(fieldnames(params))
                return
            end
            paramsFields = fieldnames(params);
            for i = 1:length(paramsFields)
                newParams.opts.(paramsFields{i}) = params.(paramsFields{i});
            end
        end
        function b_valid = validateTask(obj, className, task)
            % VALIDATETASK verifies if the task is already listed in Pipeline.
            b_valid = true;
            
            if isempty(obj.Pipe)
                return
            else
                if any(all([strcmp(className, {obj.Pipe.className}); strcmp(task.Name, {obj.Pipe.Name})],1))
                    disp(['Function with name ' task.Name ' in class ' className ' already in the Pipeline and was ignored'])
                    b_valid = false;
                end
            end
        end
        function b_isLogged = checkInFilePtr(obj, task)
           % CHECKINFILEPTR compares the information in LASTLOG with the JSON file content
           % of OBJ.TMP_TARGETOBJ. This function is used in
           % OBJ.RUN_TASKONTARGET
           %Initialize:
           b_isLogged = false;
           idx = [];
           FileInfo = obj.tmp_FilePtr.Files;
           
           if isempty(FileInfo)
               return
           end
           idx = false(length(FileInfo),5);
           for i = 1:length(FileInfo)
               idx(i,1) = strcmp(FileInfo(i).FunctionInfo.Job, task.funcStr);
               idx(i,2) = strcmp(FileInfo(i).FunctionInfo.Name, task.Name);
               idx(i,3) = isequaln(FileInfo(i).FunctionInfo.DateNum, task.DateNum);
               idx(i,4) = strcmp(FileInfo(i).InputFile_UUID, task.InputFile_UUID);
               idx(i,5) = isequaln(FileInfo(i).FunctionInfo.opts, task.opts);
           end
           idx = all(idx,2);
           b_isLogged = any(idx);    
        end
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
                subTasks{i} = obj.Pipe(b_idx);
            end
        end
        function run_tasksOnBranch(obj, branch)
            % RUN_TASKSONBRANCH finds the object to run the tasks in the
            % pipeline.
            
            % Split pipeline if there is more than one level.
            ppLine = obj.pipeSplitter;
            obj.tmp_BranchPipeline = obj.ProtocolObj.createEmptyTable;
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
                    obj.tmp_BranchPipeline = [obj.tmp_BranchPipeline; obj.tmp_TargetObj.LastLog];
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
                        errID = 'MATLAB:Umitoolbox:PipelineManager:FileNotFound';
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
                end
                LastLog.Completed = state;
                LastLog.RunDateTime = datetime('now');
                obj.tmp_TargetObj.LastLog = [obj.tmp_TargetObj.LastLog; LastLog];
                obj.PipelineSummary = [obj.PipelineSummary; LastLog];
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
        function newInput = pickAnInput(~,task)
            %PICKANINPUT prompts a LISTDLG to user-input of a function's
            % input when more than one output is possible.
            [indx,tf] = listdlg('PromptString', {'Select one of the inputs from', ['the function ' task.Name]}, 'ListString',task.Input, 'SelectionMode', 'single');
            if tf
                newInput = task.Input{indx};
            else
                newInput = '';
            end
        end
        function out = askForFirstInput(obj, FuncName)
            % ASKFORFIRSTINPUT creates an input prompt to get user to give the INPUT
            % of the first task in the pipeline.
            out = '';
            idx = obj.ProtocolObj.Idx_Filtered;
            idx = unique(idx,'rows');
            classes = {};
            for i = 1:size(idx,1)
                classes{i} = class(obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3)));
            end
            classes = [{'Subject', 'Acquisition'}  classes]; classes = unique(classes);
            
            [indx,tf] = listdlg('PromptString', {'Select the Object containing', 'the input for the function :',  FuncName},'ListString',classes, 'SelectionMode', 'single');
            if ~tf
                return
            end
            switch classes{indx}
                case 'Subject'
                    obj.tmp_TargetObj = obj.ProtocolObj.Array.ObjList(idx(1,1));
                case 'Acquisition'
                    obj.tmp_TargetObj = obj.ProtocolObj.Array.ObjList(idx(1,1)).Array.ObjList(idx(1,2));
                otherwise
                    for i = 1:size(idx,1)
                        b_isMod = strcmp(class(obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3))), classes{indx});
                        if b_isMod
                            break
                        end
                    end
                    obj.tmp_TargetObj = obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3));
            end
            obj.readFilePtr;
            FileList = {obj.tmp_FilePtr.Files.Name};
            if isempty(FileList)
                warndlg(['No Files found in ' obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3)).SaveFolder], 'Pipeline warning!', 'modal')
                return
            else
                [indx,tf] = listdlg('PromptString', {'Select the File from ' classes{indx},  'as input for the function :',  FuncName},'ListString',FileList, 'SelectionMode', 'single');
                if ~tf
                    return
                end
            end
            out = FileList{indx};
        end
        function eraseIntermediateFiles(obj)
            % ERASEINTERMEDIATEFILES deletes the file generated from the previous
            % step in the pipeline if OBJ.ERASEINTERMEDIATE is True, except
            % the first and last files in the pipeline.
            
            % Temporary solution for removing dummy rows from table OBJ.TMP_BRANCHPIPELINE:
            idx = strcmp(obj.tmp_BranchPipeline.Subject(:,1), 'None');
            obj.tmp_BranchPipeline(idx,:) = [];
            
            if ~obj.EraseIntermediate || height(obj.tmp_BranchPipeline) < 2
                return
            end
            
            for i = 3:height(obj.tmp_BranchPipeline)
                step = obj.tmp_BranchPipeline(i,:);
                
                if strcmp(step.InputFile_Path, 'None')
                    continue
                end
                
                delete(step.InputFile_Path{:});
                delete(strrep(step.InputFile_Path{:}, '.dat', '_info.mat'));
                disp(['File ' step.InputFile_Path{:} ' deleted!'])
                
                [~,fileName,ext] = fileparts(step.InputFile_Path{:});
                idx = strcmp([fileName ext], {obj.tmp_FilePtr.Files.Name});
                obj.tmp_FilePtr.Files(idx) = [];
                obj.tmp_BranchPipeline.InputFile_Path(i) = {'Deleted'};
                obj.tmp_BranchPipeline.InputFile_UUID(i) = {'Deleted'};
            end
            
            % Update OBJ.PIPELINESUMMARY
            for i = 1:height(obj.tmp_BranchPipeline)
                if strcmp(obj.tmp_BranchPipeline.InputFile_UUID(i), 'None')
                    continue
                end
                idx = strcmp(obj.PipelineSummary.Job(:), obj.tmp_BranchPipeline.Job(i));
                obj.PipelineSummary(idx,:) = obj.tmp_BranchPipeline(i,:);
            end
            % Update objects' FilePtrs.
            txt = jsonencode(obj.tmp_FilePtr);
            fid = fopen(obj.tmp_TargetObj.FilePtr, 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
        end
        function readFilePtr(obj)
            % READFILEPTR loads the content of FILEPTR.JSON in a structure
            % stored in OBJ.TMP_FILEPTR.
            txt = fileread(obj.tmp_TargetObj.FilePtr);
            a = jsondecode(txt);
            for i = 1:numel(a.Files)
                a.Files(i).Folder = tokenizePath(a.Files(i).Folder, obj.tmp_TargetObj, 'detokenize');
                a.Files(i).InputFile_Path = tokenizePath(a.Files(i).InputFile_Path, obj.tmp_TargetObj, 'detokenize');
            end
                obj.tmp_FilePtr = a;
        end
        function write2FilePtr(obj, task)
            % WRITE2FILEPTR writes the File information stored in structure FILEINFO
            % in OBJ.TMP_TARGETOBJ.FILEPTR.
            
            %Initialize
            FileInfo = struct('Name', task.FileName, 'UUID', task.File_UUID, 'Folder', task.SaveIn, 'InputFile_Path', task.Input,...
                'InputFile_UUID', task.InputFile_UUID, 'creationDateTime', datestr(now), 'FunctionInfo', ...
                struct('Name', task.Name, 'DateNum', task.DateNum, 'Job', task.funcStr, 'opts', task.opts));
            
            FileList = obj.tmp_FilePtr.Files;
            % Check for Files already logged on FilePtr
            idx = false(length(FileList),2);
            for i = 1:length(FileList)
                idx(i,1) = strcmp(FileInfo.Name, FileList(i).Name);
                idx(i,2) = strcmp(FileInfo.FunctionInfo.Name, FileList(i).FunctionInfo.Name);
            end
            idx = all(idx,2);
            % If there are no Files 
            if isempty(FileList)
                obj.tmp_FilePtr.Files = FileInfo;
            % If there are files and one identical, replace it.
            elseif ~isempty(FileList) && any(idx)
                obj.tmp_FilePtr.Files(idx) = FileInfo;
            % If there are files and none identical: Append
            else
                obj.tmp_FilePtr.Files = [FileList; FileInfo];
            end
            for i = 1:numel(obj.tmp_FilePtr.Files)
                obj.tmp_FilePtr.Files(i).Folder = tokenizePath(obj.tmp_FilePtr.Files(i).Folder, obj.tmp_TargetObj);
                obj.tmp_FilePtr.Files(i).InputFile_Path = tokenizePath(obj.tmp_FilePtr.Files(i).InputFile_Path, obj.tmp_TargetObj);
            end
            txt = jsonencode(obj.tmp_FilePtr);
            fid = fopen(obj.tmp_TargetObj.FilePtr, 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
        end
        function createFunctionList(obj)
            % CREATEFUNCTIONLIST parses all functions inside
            % OBJ.FUNCROOTDIR to get function information.
            
            % NEEDS TO BE IMPROVED OR CHANGED TO A BETTER METHOD TO GET
            % FUNCTIONS INPUTS/OUTPUTS. It is too dependent on how the
            % functions are coded.
            
            %Initialize:
            List = struct();
            default_Output = '';
            default_opts = struct();
            %
            cd(obj.FuncRootDir);
            funcList = dir('*/*.m');
            for i = 1:length(funcList)
                cd(funcList(i).folder);
                List(i).Name = funcList(i).name(1:end-2); % Remove ".m" from function name.
                List(i).FunctionFolder = funcList(i).folder;
                List(i).DateNum = funcList(i).datenum;
                info = parseFuncFile(funcList(i).name);
                fn = fieldnames(info.args);
                for k = 1:length(fn)
                    List(i).(fn{k}) = info.args.(fn{k});
                end
                List(i).opts = info.opts;
            end
            obj.FunctionList = List;
            function info = parseFuncFile(funcFile)
                info.args = [];
                info.opts = [];
                txt = fileread(funcFile);
                funcStr = regexp(txt, '(?<=function\s*).*?(?=\n)', 'match', 'once');
                expInput = ['(?<=' funcFile(1:end-2) '\s*\().*?(?=\))'];
                str = regexp(funcStr, expInput, 'match', 'once');
                str = split(str, ',');
                info.args.Input = strtrim(str{1});
                info.args.SaveIn = strtrim(str{2});
                expOutput = 'default_Output\s*=.*?(?=\n)';
                str = regexp(txt, expOutput, 'match', 'once');
                if isempty(str)
                    default_Output = '';
                else
                    eval(str)
                end
                info.args.Output = default_Output;
                expOpts = 'default_opts\s*=.*?(?=\n)';
                str = regexp(txt, expOpts, 'match', 'once');
                if ~isempty(str)
                    eval(str)
                    info.opts = default_opts;
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
                    idx = strcmp(task.Input, {obj.tmp_FilePtr.Files.Name});
                    if isempty(obj.tmp_FilePtr.Files) || sum(idx) == 0
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
                    elseif strcmp(obj.tmp_FilePtr.Files(idx).Folder, 'RawFolder')
                        filePath = fullfile(obj.tmp_TargetObj.RawFolder, obj.tmp_FilePtr.Files(idx).Name);
                    else
                        filePath = fullfile(obj.tmp_TargetObj.SaveFolder, obj.tmp_FilePtr.Files(idx).Name);
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
    end
end

