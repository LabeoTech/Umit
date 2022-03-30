classdef PipelineManager < handle
    % PIPELINEMANAGER manages data processing pipelines.
    % This class allows the creation of analysis pipeline and manages the execution of
    % functions. In addition, it controls for failed / completed  steps in the pipeline.
    
    properties
        b_ignoreLoggedFiles  = false %(bool) If true, PIPELINEMANAGER will ignore identical jobs previously run.
        b_saveDataBeforeFail = false %(bool) If true, the more recent data ("current_data") in the pipeline will be...
        % saved to a file when an error occurs.
    end
    properties (SetAccess = private)
        ClassName char % Name of the class that the pipeline analysis functions will run.
        ClassLevel int16 % Level of the class in protocol's hierarchy (1 = Modality, 2 = Acquisition, 3= Subject);
        % Structure array containing steps of the pipeline.
        pipe = struct('className', '','argsIn', {},'argsOut',{},'outFileName','',...
            'inputFileName', '','lvl', [], 'b_save2File', logical.empty, 'datFileName',...
            '', 'opts',[],'name','');% !!If the fields are changed, please apply
        % the same changes to the property's set method.
        
        fcnDir char % Directory of the analysis functions.
        funcList struct % structure containing the info about each function in the "fcnDir".
        ProtocolObj Protocol % Protocol Object.
        b_taskState = false % Boolean indicating if the task in pipeline was successful (TRUE) or not (FALSE).
        tmp_LogBook % Temporarily stores the table from PROTOCOL.LOGBOOKFILE
        tmp_BranchPipeline % Temporarily stores LogBook from a Hierarchical branch.
        PipelineSummary % Shows the jobs run in the current Pipeline
        tmp_TargetObj % % Temporarily stores an object (TARGEROBJ).
        current_task % Task structure currently running.
        current_pipe % Pipeline currently running.
        current_data % Data available in the workspace during pipeline.
        current_metaData % MetaData associated with "current_data".
        current_outFile cell % List of file names created as output from some of the analysis functions.
        b_state logical % True if a task of a pipeline was successfully executed.
        pipeFirstInput = '' % Name of the first data to be used by the Pipeline.
        % It can be the name of an existing file, or
        % "outFile" for a function that creates a file
        % such as "run_ImagesClassification".
        
        h_wbItem % Handle of waitbar dialog showing the progress of the pipeline across protocol's objects.
        h_wbTask % Handle of waitbar dialog showing the progress of the pipeline in the current object.
        targetObjFullID char % Full ID of the current object.
        % For example, if the current object is a
        % modality the full ID is:
        % SubjID--AcqID--ModID.
    end
    
    methods
        % Constructor
        function obj = PipelineManager(ProtocolObj, ClassName)
            % PIPELINEMANAGER Constructs an instance of this class
            %  Input:
            %   ProtocolObj(Protocol): Protocol object. This input is
            %   needed to have access to the protocol's hierarchy of
            %   objects.
            %   ClassName (str): Name of an existing class from
            %   ProtocolObj. The pipeline analysis functions will run on objects of
            %   this class only.
            
            p = inputParser;
            addRequired(p,'ProtocolObj', @(x) isa(x, 'Protocol'));
            addRequired(p,'ClassName', @ischar);
            parse(p,ProtocolObj, ClassName);
            obj.ProtocolObj = p.Results.ProtocolObj;
            
            % Set level of selected class:
            switch p.Results.ClassName
                case 'Subject'
                    obj.ClassLevel = 3;
                case 'Acquisition'
                    obj.ClassLevel = 2;
                otherwise
                    obj.ClassLevel = 1;
                    % Check if the selected class exists in the Filtered objects
                    % from Protocol:
                    if isempty(obj.ProtocolObj.Idx_Filtered)
                        obj.ProtocolObj.clearFilterStruct;
                        obj.ProtocolObj.queryFilter;
                    end
                    
                    % Find the class from the list of Modality objects from protocol's filtered
                    % objects:
                    items = obj.ProtocolObj.extractFilteredObjects(3);
                    idx_class = cellfun(@(x) isa(x,p.Results.ClassName), items);
                    % Check if all acquisitions contain the selected class:
                    idx_acq = ismember(obj.ProtocolObj.Idx_Filtered(:,2),...
                        unique(obj.ProtocolObj.Idx_Filtered(idx_class,2)));
                    % Control for missing modalities with the selected class
                    % name and remode acquisitions that do not have the
                    % modality.
                    if ~any(idx_acq)
                        % Throw error if there are no modalities with the name
                        % of the selected class.
                        errID = 'umIToolbox:PipelineManager:MissingInput';
                        errMsg = ['The protocol object with name "' p.Results.ClassName ...
                            '" does not exist in the filtered list of objects!'];
                        error(errID, errMsg)
                    elseif ~all(idx_acq)
                        % Remove acquisitions that lack the selected class from
                        % the protocol "Idx_Filtered" property:
                        items_errList = items(~idx_acq);
                        str = {};
                        for i = 1:length(items_errList)
                            str{i} = strjoin({items_errList{i}.MyParent.MyParent.ID, ...
                                items_errList{i}.MyParent.ID, items_errList{i}.ID},' -- ');
                        end
                        warn_msg = [{['The following Acquisition(s) lack objects with name "' ...
                            p.Results.ClassName '" and will be removed from this pipeline session:']}; str];
                        warndlg(warn_msg, 'Warning! Missing Objects in protocol')
                        % Update list of selected objects from protocol:
                        obj.ProtocolObj.Idx_Filtered = obj.ProtocolObj.Idx_Filtered(idx_acq, :);
                    end
            end
            
            obj.ClassName = p.Results.ClassName;
            root = getenv('Umitoolbox');
            obj.fcnDir = fullfile(root, 'Analysis');
            obj.createFcnList;
            obj.b_ignoreLoggedFiles = false;
            obj.b_saveDataBeforeFail = false;
        end
        % SETTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.pipe(obj, pipe)
            % Pipeline structure setter. If pipe is empty, create an empty
            % structure containing tasks fields.
            
            if isempty(fieldnames(pipe))
                pipe = struct('className', '','argsIn', {},'argsOut',{},'outFileName','',...
                    'inputFileName', '','lvl', [], 'b_save2File', logical.empty, 'datFileName',...
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
            if isempty(S)
                disp(['The function ' obj.funcList(idx).name ' does not have any optional parameters.']);
                return
            end
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
        
        function addTask(obj,func, varargin)
            % This method adds an analysis function, or task, to the
            % pipeline. Here, we can choose to save the output of a given
            % task as a .DAT file. ADDTASK will create a string containing
            % the task that will be evaluated during pipeline execution.
            % Inputs:
            %   func (str || char):  name or index of the analysis function
            %       contained in obj.funcList property.
            %   b_save2File (bool): Optional. True means that the output data
            %       from the function will be saved as a .DAT file.
            %   datFileName(char): Optional. Name of the file to be saved.
            %       If not provided, the analysis function's default
            %       filename will be used.
            
            p = inputParser;
            addRequired(p, 'func', @(x) ischar(x) || isnumeric(x));
            addOptional(p, 'b_save2File', false, @islogical);
            addOptional(p, 'datFileName', '', @ischar);
            parse(p,func, varargin{:});
            
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
            task.className = obj.ClassName;
            task.lvl = obj.ClassLevel;
            task.name = obj.funcList(idx).name;
            % Add default values for fields used to save the data: These
            % values will be updated later:
            task.b_save2File = false;
            task.datFileName = '';
            
            % Control for steps IDENTICAL to the task that are already in the pipeline:
            idx_equal = arrayfun(@(x) isequaln(task,x), obj.pipe);
            if any(idx_equal)
                warning('Operation cancelled! The function "%s" already exists in the Pipeline!',....
                    task.name);
                return
            end
            
            % Look for first input file to the pipeline;
            if isempty(obj.pipeFirstInput)
                task.inputFileName = obj.getFirstInputFile(task);
                if task.inputFileName == 0
                    disp('Operation Cancelled by User')
                    return
                end
            end
            
            % Control for multiple outputs from the previous step:
            % Here, we assume that functions with multiple outputs
            % create only "Files" and not "data".
            % Therefore, we will update the function string to load
            % one of the "Files" before running the task.
            
            % Look from bottom to top of the pipeline for tasks with files
            % as outputs. This is necessary because not all analysis
            % functions have outputs.
            
            for i = length(obj.pipe):-1:1
                if ismember('outData', obj.pipe(i).argsOut)
                    break
                elseif any(strcmp(task.argsIn, 'data')) && any(strcmp('outFile', obj.pipe(i).argsOut))
                    if iscell(obj.pipe(i).outFileName)
                        % Ask user to select a file:
                        disp('Controlling for multiple outputs')
                        w = warndlg({[obj.pipe(i).name ' has multiple output files!'],...
                            'Please, select one to be analysed!'});
                        waitfor(w);
                        [indxFile, tf] = listdlg('ListString', obj.pipe(i).outFileName,...
                            'SelectionMode','single');
                        if ~tf
                            disp('Operation cancelled by User')
                            return
                        end
                        task.inputFileName = obj.pipe(i).outFileName{indxFile};
                    else
                        task.inputFileName = obj.pipe(i).outFileName;
                    end
                end
            end
            
            % Save to Pipeline:
            obj.pipe = [obj.pipe; task];
            disp(['Added "' task.name '" to pipeline.']);
            
            % Control for data to be saved as .DAT files for task:
            if ~p.Results.b_save2File
                return
            elseif ~any(ismember({'outData', 'outDataStat'}, task.argsOut))
                warning(['Cannot save output to .DAT file for the function'...
                    ' "%s" \nbecause it doesn''t have any data as output!'], task.name);
                obj.pipe(end).b_save2File = false;
                obj.pipe(end).datFileName = '';      
                return
            else
                % switch the b_save2File variable to TRUE:
                obj.pipe(end).b_save2File = true;
            end
            
            % Save datFileName as default output name from task's function
            % if it wasn't previously defined by the User:
            if isempty(p.Results.datFileName)
                obj.pipe(end).datFileName = obj.pipe(end).outFileName;
            else
                % OR update datFileName to add file extension:
                [~,~,ext]= fileparts(p.Results.datFileName);
                if isempty(ext)
                    [~,~,ext_def] = fileparts(obj.pipe(end).outFileName);
                    obj.pipe(end).datFileName = [p.Results.datFileName, ext_def];
                end                
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
            
            str = sprintf('Pipeline Summary:\nRun on object: "%s"\n\n', obj.ClassName);
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
                if obj.pipe(i).b_save2File
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
        
        function saveLastData(obj)
            % This method ensures that the last step of the pipeline that
            % ouputs data ("outData") is saved. It is called from "run_pipeline".
            % This method updates the "datFileName" and "b_save2File" variables
            % from the last "saveble" step of the pipeline.
            %
            
            % Look for last step that ouputs a "data":
            for i = length(obj.pipe):-1:1
                if ~isempty(obj.pipe(i).outFileName)
                    break
                end
            end
            
            if obj.pipe(i).b_save2File & ~isempty(obj.pipe(i).datFileName)
                return
            end
            
            
            % Save file of last step if the function outputs Data or
            % StatsData:
            if isempty(obj.pipe(i).datFileName) && any(ismember({'outData', 'outDataStat'},...
                    obj.pipe(i).argsOut))
                obj.pipe(i).datFileName = obj.pipe(i).outFileName;
                obj.pipe(i).b_save2File = true;
            end
            
        end
                
        function run_pipeline(obj)
            % RUN_PIPELINE runs the tasks in OBJ.PIPE
            lbf = matfile(obj.ProtocolObj.LogBookFile);
            obj.tmp_LogBook = lbf.LogBook;
            obj.PipelineSummary = obj.ProtocolObj.createEmptyTable;
            
            % Check if last data will be saved:
            obj.saveLastData;
            % Initialize waitbars:
            obj.setWaitBar('Initialize')
            
            % Get indexes of Filtered Objects from OBJ.PROTOCOLOBJ.QUERYFILTER function.
            switch obj.ClassLevel
                case 1
                    % Modality
                    targetIdxArr = unique(obj.ProtocolObj.Idx_Filtered,'rows');
                case 2
                    % Acquisition
                    targetIdxArr = unique(obj.ProtocolObj.Idx_Filtered(:,[1 2]), 'rows');
                case 3
                    % Subject
                    targetIdxArr = unique(obj.ProtocolObj.Idx_Filtered(:,1), 'rows');
            end
            
            for i = 1:size(targetIdxArr,1)
                % Clear current data,  metaData and File List before starting the pipeline:
                obj.current_data = []; obj.current_metaData = [];obj.current_outFile = {};
                obj.b_state = true;
                % Get handle of current object:
                obj.getTargetObj(targetIdxArr(i,:));
                % Initialize Log table:
                obj.tmp_TargetObj.LastLog = obj.ProtocolObj.createEmptyTable;
                % Create Full ID of object:
                myParent = obj.tmp_TargetObj.MyParent;
                myFullName = {obj.tmp_TargetObj.ID};
                while ~isa(myParent, 'Protocol')
                    myFullName = [myFullName, {myParent.ID}];
                    myParent = myParent.MyParent;
                end
                myFullName = fliplr(myFullName);
                obj.targetObjFullID = strjoin(myFullName, ' -- ');
                % Update waitbars:
                obj.setWaitBar('UpdateItem', i, size(targetIdxArr,1));
                fprintf([repmat('-',1,50),'\n']);
                fprintf('Object Name: %s\n\n', obj.targetObjFullID)
                % Run pipeline in each target object:
                for j = 1:length(obj.pipe)
                    obj.current_task = obj.pipe(j);
                    obj.setWaitBar('UpdateTask', j/length(obj.pipe));
                    fprintf('Running task # %d/%d ----->>>>>\n',j,length(obj.pipe));
                    obj.run_taskOnTarget;
                    
                    if ~obj.b_state
                        break
                    end
                    % Control for Pipeline cancelling by User:
                    if getappdata(obj.h_wbItem, 'b_abortPipe')
                        % Delete waitbars and abort inner loop:
                        delete([obj.h_wbItem, obj.h_wbTask])
                        break
                    end
                    % This pause is here to allow the WaitBar to update
                    % during the execution of this method.
                    pause(.001);
                end
                %                 % When the pipeline reaches the last step, force to save the current
                %                 % data to a file:
                %                 if j == length(obj.pipe)
                %                     obj.saveDataToFile(obj.pipe(end));
                %                 end
                % Update Pipeline summary table:
                obj.PipelineSummary = [obj.PipelineSummary; obj.tmp_TargetObj.LastLog];
                fprintf([repmat('-',1,50),'\n']);
                % Abort outer loop if user cancels pipeline:
                if ~ishandle(obj.h_wbItem)
                    break
                end
            end
            
            % Remove "empty" rows from the Pipeline Summary Log table:
            idx_emptyRow = all(strcmp('None',table2cell(obj.PipelineSummary(:,1:5))),2);
            obj.PipelineSummary(idx_emptyRow,:) = [];
            % Update LogBook with Pipeline Summary table:
            LogBook = [obj.tmp_LogBook; obj.PipelineSummary];
            % Save Log Book to file:
            save(obj.ProtocolObj.LogBookFile, 'LogBook');
            % Show Pipeline Summary in command window:
            disp(obj.PipelineSummary)
            
            % Save Protocol Object:
            protocol = obj.ProtocolObj;
            save([obj.ProtocolObj.SaveDir obj.ProtocolObj.Name '.mat'], 'protocol');
            disp('Protocol object Saved!');
            delete([obj.h_wbItem, obj.h_wbTask]);
            
        end
        % Pipeline Management methods:
        function savePipe(obj, filename)
            % SAVEPIPE saves the structure OBJ.PIPE in a .JSON file in the
            % folder PIPELINECONFIGFILES inside the SAVEDIR of OBJ.PROTOCOLOBJ.
            
            targetDir = fullfile(obj.ProtocolObj.SaveDir, 'PipeLineConfigFiles');
            [~,~] = mkdir(targetDir);
            pipeStruct = obj.pipe;
            txt = jsonencode(pipeStruct);
            fid = fopen(fullfile(targetDir,[filename '.json']), 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            disp(['Pipeline saved as "' filename '" in ' targetDir]);
        end
        
        function loadPipe(obj, pipeFile)
            % LOADPIPE loads the structure PIPE inside FILENAME and assigns
            % it to OBJ.PIPE property.
            % Input:
            %   pipeFile(char): full path to the .JSON file containing the
            %   pipeline config.
            txt = fileread(pipeFile);
            new_pipe = jsondecode(txt);
            
            % erase current pipeline:
            obj.reset_pipe;
            
            % Add new tasks:
            for i = 1:length(new_pipe)
                indx_name = find(strcmp(new_pipe(i).name, {obj.funcList.name}));
                % Update funcList with custom opts settings:
                if ~isequaln(new_pipe(i).opts,obj.funcList(indx_name).info.opts)
                    obj.funcList(indx_name).info.opts = new_pipe(i).opts;
                end
                % add tasks to pipeline:
                if new_pipe(i).b_save2File
                    obj.addTask(indx_name, true, new_pipe(i).datFileName);
                else
                    obj.addTask(indx_name);
                end
            end
        end
        
        function reset_pipe(obj, varargin)
            % This function erases the pipe property and resets the funcList
            % property to default parameter values.
            % Input:
            %   flag(char): (optional) type 'all' to reset function list in
            %   addition to the pipeline.
            flag = '';
            if nargin > 1
                flag = varargin{:};
            end
            obj.pipe = struct();
            obj.pipeFirstInput = '';
            if strcmp(flag, 'all')
                obj.funcList = struct.empty;
                obj.createFcnList;
            end
            % Clear current data,  metaData and File List:
            obj.current_data = []; obj.current_metaData = [];obj.current_outFile = {};
        end
    end
    
    methods (Access = private)
        
        function run_taskOnTarget(obj)
            % RUN_TASKONTARGET runs a task in the pipeline structure array in
            % TASK.
            %   It checks if the command TASK.FUNCNAME was already
            %   sucessfully performed by comparing it to the LOGBOOK from
            %   the PROTOCOL object. Also, it appends the information of
            %   the task on the object's LASTLOG.
            
            % Current task:
            task = obj.current_task;
            % Initialize empty Log for current object:
            LastLog = obj.ProtocolObj.createEmptyTable;
            % Fill out Log with Subject/Acquisition/Modality IDs :
            cnt = 2;
            tmpObj = obj.tmp_TargetObj;
            ID_list = {tmpObj.ID};
            while ~isa(tmpObj.MyParent, 'Protocol')
                tmpObj = tmpObj.MyParent;
                ID_list{cnt} = tmpObj.ID;
                cnt = cnt+1;
            end
            LastLog(:,1:3) = fliplr(ID_list);
            clear tmpObj
            % Add class name to table:
            LastLog(:,4) = {task.className};
            LastLog(:,5) = {task.name};
            %%%
            
            % Create function string and update log table:
            task.funcStr = createFcnString(obj, task);
            LastLog.Job = {task.funcStr};
            %  Execute the task:
            try
                % Control for missing input files:
                if task.inputFileName
                    errID = 'MATLAB:Umitoolbox:PipelineManager:FileNotFound';
                    errmsg = ['Input File for function ' task.name ' not found!'];
                    assert(isfile(fullfile(obj.tmp_TargetObj.SaveFolder, task.inputFileName)),...
                        errID,errmsg);
                    obj.loadInputFile(task);
                end
                % Check for data already run and skip step if so:
                b_skipStep = obj.checkDataHistory(task);
                
                if b_skipStep && ~obj.b_ignoreLoggedFiles
                    disp(['Skipped function : "' task.name '"!']);
                    LastLog.Messages_short = 'Skipped';
                    LastLog.Completed = true;
                    return
                end
                fprintf('\tFunction Name: %s \n\n',task.name);
                % Load options structure in the workspace.
                opts = task.opts; %#ok the "opts" structure is used in the EVAL function.
                
                % Evaluate function string:
                eval(task.funcStr);
                % Update log table and tell other methods that the function
                % was successfully run:
                obj.b_state = true;
                LastLog.Messages = 'No Errors';
                % Update data history of current data with task:
                obj.updateDataHistory(task);
            catch ME
                obj.b_state = false;
                LastLog.Messages = {getReport(ME)};
                LastLog.Messages_short = {getReport(ME, 'basic','hyperlinks','off')};
                disp('FAILED!');                
            end
            % Save data to file:
            if task.b_save2File && obj.b_state                
                % Look for tasks with output data from the current step and
                % save the data to a .DAT or .MAT file:
                obj.saveDataToFile(task,false)
            elseif obj.b_saveDataBeforeFail && ~obj.b_state
                % Look for tasks from the previous steps given that the
                % current one failed and save it to a file:
                indx = find(strcmp(task.name, {obj.pipe.name}));
                if indx > 1
                    obj.saveDataToFile(obj.pipe(indx-1),true);
                end
            end
            
            % Update log table of target object:
            LastLog.Completed = obj.b_state;
            LastLog.RunDateTime = datetime('now');
            obj.tmp_TargetObj.LastLog = [obj.tmp_TargetObj.LastLog; LastLog];
            % Remove "empty" rows from the target Object Log table:
            idx_emptyRow = all(strcmp('None',table2cell(obj.tmp_TargetObj.LastLog(:,1:5))),2);
            obj.tmp_TargetObj.LastLog(idx_emptyRow,:) = [];
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
            % This method verifies if the function has any data as input.
            % If yes, then a creates an dialog box containing a list of .DAT files
            % to choose as input. This method is called by ADDTASK only
            % when the first task of a pipeline is created.
            % Input:
            %   funcInfo (struct): structure containing the task's function info.
            % Output:
            %   out (char): name of input file. Empty if file does not exist.
            
            out = '';
            % Control for function that creates the first input:
            if any(strcmp('outFile', funcInfo.argsOut))
                obj.pipeFirstInput = 'outFile';
                return            
            elseif any(strcmp('outData', funcInfo.argsOut)) && ...
                    ~any(ismember({'data', 'dataStat'}, funcInfo.argsIn))
                obj.pipeFirstInput = 'outData';
                return
            end
            
            % Control for tasks that do not have any data as input:
            if ~any(ismember({'data', 'dataStat'}, funcInfo.argsIn))
                return
            end
            
            
            % Get all existing objects from the selected items in Protocol
            % object:
            idx = unique(obj.ProtocolObj.Idx_Filtered,'rows');
            classes = {};
            for i = 1:size(idx,1)
                classes{i} = class(obj.ProtocolObj.Array.ObjList(idx(i,1)).Array.ObjList(idx(i,2)).Array.ObjList(idx(i,3)));
            end
            classes = [{'Subject', 'Acquisition'}  classes]; classes = unique(classes);
            % Select Object containing input file:
            [indx,tf] = listdlg('PromptString', {'Select the Object containing', 'the input for the function :',...
                funcInfo.name},'ListString',classes, 'SelectionMode', 'single');
            if ~tf
                out = 0;
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
            
            % Display a list of files from the selected
            % object(targetObj):
            datFileList = dir(fullfile(targetObj.SaveFolder, '*.dat'));
            if any(strcmp('dataStat', funcInfo.argsIn))
                % Select only .MAT files containing "dataHistory" variable.
                % This was a way to exclude other .MAT files in the folder
                % that are not Stats files.
                matFileList = dir(fullfile(targetObj.SaveFolder, '*.mat'));
                matFilesMap = arrayfun(@(x) matfile(fullfile(x.folder,x.name)), matFileList, 'UniformOutput',false);
                b_validMat = cellfun(@(x) isprop(x, 'dataHistory'), matFilesMap);
                [~,datFileNames,~] = arrayfun(@(x) fileparts(x.name), datFileList, 'UniformOutput', false);
                [~,matFileNames,~] = arrayfun(@(x) fileparts(x.name), matFileList, 'UniformOutput', false);
                b_statMat = ~ismember(matFileNames, datFileNames);
                % Update index of valid .mat files:
                b_validMat = b_validMat & b_statMat;
                FileList = matFileList(b_validMat);
            else
                FileList = datFileList;
            end
            if isempty(FileList)
                warndlg(['No valid Files found in ' targetObj.SaveFolder], 'pipeline warning!', 'modal')
                out = 0;
                return
            else
                [indx,tf] = listdlg('PromptString', {'Select the File from ' classes{indx},...
                    'as input for the function :',  funcInfo.name},'ListString',{FileList.name}, 'SelectionMode',...
                    'single');
                if ~tf
                    out = 0;
                    return
                end
            end
            
            % Save first input
            out = FileList(indx).name;
            obj.pipeFirstInput = out;
        end
        
        function createFcnList(obj)
            % This function creates a structure containing all information
            % about the analysis functions inside the "Analysis" folder.
            % This information is stored in the "funcList" property of
            % pipelineManager.
            
            % !!For now, it will read only folders directly below the
            % "Analysis" folder. Subfolders inside these folders will not
            % be read!
            
            % Set Defaults:
            default_Output = '';
            default_opts = struct();
            
            disp('Creating Fcn list...');
            list = dir(fullfile(obj.fcnDir, '\*\*.m'));
            for i = 1:length(list)
                out = parseFuncFile(list(i));
                % Validate if all input arguments from the function are
                % "valid" inputs keywords:
                kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts', 'object', 'dataStat'};
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
                info = struct('argsIn', {},'argsOut', {}, 'outFileName', '', 'opts', []);
                txt = fileread(fullfile(fcnStruct.folder, fcnStruct.name));
                funcStr = erase(regexp(txt, '(?<=function\s*).*?(?=\r*\n)', 'match', 'once'),' ');
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
                % SPECIAL CASE. Look for "object"s as optional arguments in
                % function first lines:
                expObj = 'default_object\s*=';
                str = regexp(txt, expObj, 'match', 'once');
                if ~isempty(str)                    
                    info.argsIn{end+1} = 'object';
                end
            end
        end
        
        function idx_fcn = check_funcName(obj, func)
            % This function is used by "setOpts" and "addTask" methods to
            % validate if the input "func" is valid.
            % Input
            %   func (numeric OR char) : index OR name of a function from
            %   obj.funcList.
            % Output:
            %   idx_fcn(bool): index of "func" in "obj.funcList". Returns
            %   empty if the function was not found in the list.
            
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
        
        function fcnStr = createFcnString(obj, task)
            % This method creates a string containing the function to be
            % called during the current task of the pipeline by an EVAL
            % statement. This method is called from "run_taskOnTarget"
            % method.
            % !! If new input or output arguments are created, meaning new
            % argument keywords, this function has to be updated in order the
            % function string to work.
            
            % Input:
            %   task (struct): info of current function to be run on the
            %       current object.
            % Output:
            %   fcnStr (char): string containing call to analysis function
            %   in the current task.
            
            % Create analysis function string:
            fcnStr = '';
            % Replace input argument names:
            % Important! Put "dataStat" and "outDataStat" first in the string list, otherwise 
            % the replace function will create a non-existant property name:
            
            argsIn = replace(task.argsIn, ["RawFolder", "SaveFolder", "dataStat","metaData", "object", "data"],...
                {['''' obj.tmp_TargetObj.RawFolder ''''],['''' obj.tmp_TargetObj.SaveFolder ''''], 'obj.current_data',...
                'obj.current_metaData', 'obj.tmp_TargetObj', 'obj.current_data'});
            
            argsOut = replace(task.argsOut, ["outDataStat", "metaData", "outData", "outFile"],...
                {'obj.current_data', 'obj.current_metaData', 'obj.current_data', 'obj.current_outFile'}); 
            if isempty(argsOut)
                fcnStr = [fcnStr ';' task.name '(' strjoin(argsIn,',') ');'];
            else
                fcnStr = [fcnStr ';' '[' strjoin(argsOut, ',') ']=' task.name '(' strjoin(argsIn,',') ');'];
            end
            fcnStr = strip(fcnStr,'left', ';');
            
        end
        
        function updateDataHistory(obj, step)
            % This function creates or  updates the "dataHistory" structure
            % and saves the information to the metaData structure/matfile.
            % The dataHistory contains all information about the functions'
            % parameters used to create the current "data" and when it was run.
            %
            % Input:
            %    step(struct) : current step of the pipeline;
            
            funcInfo = obj.funcList(strcmp(step.name, {obj.funcList.name}));
            % Create a local structure with the function's info:            
            curr_dtHist = genDataHistory(funcInfo, step.funcStr, step.opts,'none');           
            % First, we need to know if the output is a "data", a .DAT file or a .MAT file:
            if any(strcmp(step.argsOut, 'outFile'))
                % In case the step ouput is .DAT file(s):
                
                % Get only filename instead of full path:
                [~, filenames, ext] = cellfun(@(x) fileparts(x), obj.current_outFile,...
                    'UniformOutput', false);
                curr_dtHist.outputFile_list = join([filenames',ext'],'');
                
                for i = 1:length(obj.current_outFile)
                    % Map existing metaData file to memory:
                    mtD = matfile(strrep(obj.current_outFile{i}, '.dat', '.mat'));
                    mtD.Properties.Writable = true;
                    % Create or update "dataHistory" structure:
                    if isprop(mtD, 'dataHistory')
                        mtD.dataHistory = [mtD.dataHistory; curr_dtHist];
                    else
                        mtD.dataHistory = curr_dtHist;
                    end
                end
            elseif any(strcmp(step.argsOut, 'outDataStat'))
                % In case of step output is .MAT file(s):
                if isfield(obj.current_data, 'dataHistory')
                    obj.current_data.dataHistory = [obj.current_data.dataHistory; curr_dtHist];
                else
                    obj.current_data.dataHistory = curr_dtHist;
                end
                
            else
                % In case of step output is a data array:
                if isfield(obj.current_metaData, 'dataHistory')
                    obj.current_metaData.dataHistory = [obj.current_metaData.dataHistory; curr_dtHist];
                else
                    obj.current_metaData.dataHistory = curr_dtHist;
                end
            end
        end
        
        function b_skip = checkDataHistory(obj,step)
            % This method verifies if the function to be run in "step" was
            % already performed or not on the current data.
            % If so, the pipeline step will be skipped.
            % Input:
            %   step (struct): structure containing the function information that
            %   will run on the data.
            % Output:
            %   b_skip (bool): True if the step was already run on the
            %   data and should be skipped.
            
            b_skip = false;
            % Find function info in Function List:
            fcnInfo = obj.funcList(strcmp(step.name, {obj.funcList.name}));
            % Find step info in object's dataHistory:
            
            % For retro-compatibility with data created in previous
            % versions of umIT:
            if ~isfield(obj.current_metaData, 'dataHistory')
                return
            end
            
            dH = obj.current_metaData.dataHistory(strcmp(step.name,...
                {obj.current_metaData.dataHistory.name}));
            % If the function's creation date AND the function string AND optional parameters are
            % the same, we consider that the current step was already run.
            if isempty(dH)
                return
            elseif ( isequal(datetime(fcnInfo.date), dH.creationDatetime) &&...
                    strcmp(step.funcStr, dH.funcStr) ) && isequaln(step.opts, dH.opts)
                b_skip = true;
            end
        end
        
        function loadInputFile(obj,step)
            % This function loads the data and metaData (if applicable)
            % from a .DAT or .MAT file indicated by the inputFileName field
            % in the "step" structure.
            % Input:
            %   step (struct): structure containing the info of the current
            %   task that will be executed inside method
            %   "run_taskOnTarget".
            
            
            if endsWith(step.inputFileName, '.dat')
                %If the InputFile is a .DAT file:
                [obj.current_data, obj.current_metaData] = ...
                    loadDatFile(fullfile(obj.tmp_TargetObj.SaveFolder, step.inputFileName));
            else
                % If the InputFile is a .MAT file
                obj.current_data = load(fullfile(obj.tmp_TargetObj.SaveFolder,...
                    step.inputFileName));
                % Erase current metaData, since it will not be associated
                % with the curren_data anymore:
                obj.current_metaData = [];
            end
            
        end
               
        function saveDataToFile(obj, step, b_failed)
            % This methods looks back in the pipeline from "step" for tasks
            % with "data" or "stats data" as output and saves the current data to a
            % .DAT or .MAT file.
            % Input:
            %    step(struct) : info of the current task in the pipeline.
            %    b_failed (bool): If TRUE, this function will ignore the
            %    "b_save2File" from "step" and save the data.
            % Get the pipeline until the task in "step":
            indx = find(strcmp(step.name, {obj.pipe.name}));
            subPipe = obj.pipe(1:indx);
            % Look back in pipeline for steps with "data" or "stats data"
            % as output and save the current data using the task's info:            
            for i = length(subPipe):-1:1
                task = subPipe(i);                   
                % If the pipeline failed, and the datFileName was not set,
                % use the default file name to save the data:
                if b_failed & isempty(task.datFileName) & ischar(task.outFileName)
                    task.datFileName = task.outFileName;
                end
                
                if endsWith(task.datFileName, '.dat')
                    save2Dat(fullfile(obj.tmp_TargetObj.SaveFolder,task.datFileName),...
                        obj.current_data, obj.current_metaData);
                    return
                elseif endsWith(task.datFileName, '.mat')
                    S = obj.current_data;
                    save(fullfile(obj.tmp_TargetObj.SaveFolder,task.datFileName),...
                        '-struct', 'S', '-v7.3');
                    return
                end
            end
        end
        %%%%%%%%%%% WAITBAR METHODS: %%%%%%%%%%%%%%%
        function setWaitBar(obj, tag, varargin)
            % This method creates two "waitbar" dialogs.
            % The first shows the progress of the pipeline runs across objects
            % while the second one shows the progress of tasks in a given
            % object.
            % Inputs:
            %   tag (char): "Initialize" : creates the 2 waitbars.
            %               "UpdateItem" : updates bar1.
            %               "UpdateTask" : updates bar2.
            %   barVal (float): (Optional) fractional value of the bar.
            
            % Control for invalid waitbar handles:
            if ~strcmp(tag, 'Initialize')
                b_HandleExist = (ishandle(obj.h_wbItem) & ishandle(obj.h_wbTask));
                if ~b_HandleExist
                    return
                end
            end
            
            
            switch tag
                case 'Initialize'
                    obj.h_wbItem = waitbar(0,'Initializing Pipeline...',...
                        'Name','Pipeline Progress',...
                        'CreateCancelBtn', @obj.wb_cancelBtn);
                    setappdata(obj.h_wbItem, 'b_abortPipe', 0);
                    obj.h_wbItem.Children(2).Title.Interpreter = 'none';
                    
                    obj.h_wbTask = waitbar(0,'Initializing Task...','Name','',...
                        'CloseRequestFcn',@DoNothing);
                    obj.h_wbTask.Children(1).Title.Interpreter = 'none';
                    
                case 'UpdateItem'
                    waitbar(varargin{1}/varargin{2}, obj.h_wbItem, ['Item ' num2str(varargin{1})...
                        '/' num2str(varargin{2})]);
                    obj.h_wbTask.Name = obj.targetObjFullID;
                case 'UpdateTask'
                    waitbar(varargin{1}, obj.h_wbTask, ['Running "' obj.current_task.name '"']);
            end
            
            function DoNothing(~,~)
                % Empty callback to avoid closing Waitbar #2
            end
        end
        
        function wb_cancelBtn(obj,src,evnt)
            % Callback of cancel button in waitbar1 ("h_wbItem").
            % This callback triggers the cancellation of the current
            % pipeline when the user clicks on the cancel button.
            
            fprintf('>>>>>>>>>>>>>>>>>>Cancelling Pipeline...>>>>>>>>>>>\n');
            
            src.String = 'Wait!';
            if strcmp(evnt.EventName, 'Action')
                set(src.Parent.Children(2).Title, 'String', 'Please Wait. Stopping Pipeline...');
                setappdata(obj.h_wbItem, 'b_abortPipe', 1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end
