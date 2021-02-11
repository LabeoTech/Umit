classdef PipelineManager < handle
    % PIPELINEMANAGER Class that creates and manages pre-processing and
    % analysis pipelines.
    %  Creates a structure with the info necessary to run a simple pipeline
    %  as a sequence of steps(i.e. step 1, step 2, step3, ... step N).
    %  Controls for run / failed / completed / aborted steps in the
    %  pipeline.
    properties
        EraseIntermediate % Boolean to indicate if PIPELINEMANAGER will delete intermediate files in the pipeline. It keeps only the first and last file in the pipeline.
        IgnoreLogBook % Boolean. If true, PIPELINEMANAGER will ignore identical jobs previously run (logged in OBJ.PROTOCOLOBJ.LOGBOOKFILE).
    end
    properties (SetAccess = private, SetObservable)
        Pipe % Array containing steps of the pipeline.
        ProtocolObj % Protocol Object.
        State % Boolean indicating if the task in pipeline was successful (TRUE) or not (FALSE).
        tmp_LogBook % Temporary stores the table in PROTOCOL.LOGBOOKFILE
        tmp_BranchPipeline % Temporary stores LogBook from a Hierarchical branch.
        PipelineSummary % Shows the jobs run in the current Pipeline
        FunctionList % Structure containing the Public methods of all Classes in the ISAtoolbox with its inputs, optional parameters and outputs.
    end
    
    methods
        % Constructor
        function obj = PipelineManager(Pipe, ProtocolObj)
            % PIPELINEMANAGER Construct an instance of this class
            %   PIPE is a structure containing all the information of the pipeline.
            if nargin > 0
                obj.Pipe = Pipe;
                obj.ProtocolObj = ProtocolObj;
            end
            obj.readJSON;
            obj.EraseIntermediate = false;
            obj.IgnoreLogBook = false;
        end
        % Set functions
        function set.ProtocolObj(obj, ProtocolObj)
            mustBeA(ProtocolObj, 'Protocol');
            obj.ProtocolObj = ProtocolObj;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
        function addTask(obj,className, funcName, params)
            arguments
                obj
                className char
                funcName char
                params.input char
                params.opts struct
                params.output char
            end
            funcInfo = getfield(obj.FunctionList, className, funcName);
            params = obj.validateParams(params,funcInfo);
            
            % Determine Level of the step in the hierarchy.
            switch className
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
            % Validates if FUNCNAME exists in CLASSNAME:
            b_valid = obj.validateTask(className, funcName);
            if ~b_valid
                return
            end
            task.className = className;
            task.funcName = funcName;
            % Determine the task INPUT:
            if isempty(obj.Pipe) && ~isempty(params.input)
                task.input = obj.askForFirstInput(funcName);
                if isempty(task.input)
                    disp('Operation cancelled by user!')
                    return
                end
            elseif isempty(obj.Pipe) && isempty(params.input)
                task.input = params.input;
            else
                task.input = obj.Pipe(end).output;
            end
            if iscell(task.input) && numel(task.input) > 1
                task.input = obj.pickAnInput(task);
            end
            task.output = params.output;
            task.opts = params.opts;
            task.level = lvl;
            funcStr = obj.createFuncStr(task,funcName);
            task.funcStr = funcStr;
            obj.Pipe = [obj.Pipe;task];
        end
        
        function opts = setOpts(obj, className, funcName)
            % SETOPTS opens an INPUTDLG for entry of optional variables
            % (OPTS) of methods in the Pipeline. Output: Structure
            % containing the variables and values.
            S = getfield(obj.FunctionList, className, funcName, 'opts');
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
                    for i = 1:length(uniqA)
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
                obj.State = true;
                waitbar(i/length(uniq_branch),f);
                branch = idxList(ppIdx == uniq_branch(i),:);
                obj.run_tasksOnBranch(branch);
                if ~obj.State
                    continue
                end
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
        
        
        function savePipe(obj, filename)
            % SAVEPIPE saves the structure OBJ.PIPE in a .MAT file in the
            % folder PIPELINECONFIGFiles inside MAINDIR of OBJ.PROTOCOLOBJ.
            maindir = obj.ProtocolObj.MainDir;
            cd(maindir)
            [~]= mkdir('PipeLineConfigFiles');
            cd('PipeLineConfigFiles');
            pipeStruct = obj.Pipe;
            save([filename '.mat'], 'pipeStruct');
        end
        
        function loadPipe(obj, filename)
            % LOADPIPE loads the structure PIPE inside FILENAME and assigns
            % it to OBJ.PIPE property.
            if ~isfile(filename)
                path = which([obj.MainDir PipeLineConfigFiles]);
                filename = fullfile(path, filename);
            end
            a = load(filename);
            obj.Pipe = a.pipeStruct;
            clear a
        end
        
    end
    
    methods (Access = private)
        function readJSON(obj)
            % Creates FUNCTIONLIST from JSON file in the ISAtoolbox subFunc
            % folder.
            JSONfile = which('ISAtoolboxFunctionInventory.json');
            a = fileread(JSONfile);
            b = jsondecode(a);
            obj.FunctionList = b.ISAtoolbox_methodsInventory;
        end
        
        function b_valid = validateTask(obj, className, funcName)
            % VALIDATETASK verifies if the class named CLASSNAME contains
            % the function named FUNCNAME. In addition, it ignores
            % repetitions. Outputs TRUE if the task is valid and FALSE if
            % it is not.
            b_valid = true;
            funcList = methods(className); funcList = setdiff(funcList, [{className}; methods('handle')]); % lists all public methods from class of type CLASSNAME.
            if ~any(strcmp(funcList, funcName))
                disp(['Cannot find function with name ''' funcName ''' in object of type ''' className '''. Check the function name and try again.']);
                b_valid = false;
                return
            end
            
            if isempty(obj.Pipe)
                return
            else
                if any(all([strcmp(className, {obj.Pipe.className}); strcmp(funcName, {obj.Pipe.funcName})],1))
                    disp(['Function with name ' funcName ' in class ' className ' already in the Pipeline and was ignored'])
                    b_valid = false;
                end
            end
        end
        function isLogged = checkInLogBook(~,LogBook, LogTable)
            % This function checks if the job in the pipeline has already
            % been successfully processed.
            isLogged = false;
            if isempty(LogBook)
                return
            end
            c1 = table2array(LogBook(:,1:6));
            c2 = table2array(LogTable(:,1:6));
            idx = false(size(c1));
            for i = 1:size(c1,2)
                idx(:,i) = strcmp(c2(i), c1(:,i));
            end
            idx = all(idx,2);
            a = any(idx);
            b = any(LogBook.Completed(idx));
            
            if a && b
                isLogged = true;
            end
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
                    targetObj = obj.getTargetObj(targetIdxArr(j,:));
                    LastLog = obj.ProtocolObj.createEmptyTable;
                    targetObj.LastLog = LastLog;
                    for k = 1:length(subtasks)
                        task = subtasks(k);
                        obj.run_taskOnTarget(targetObj, task);
                        if ~obj.State
                            return
                        end
                    end
                    obj.tmp_BranchPipeline = [obj.tmp_BranchPipeline; targetObj.LastLog];
                end
            end
            obj.eraseIntermediateFiles();
        end
        
        function run_taskOnTarget(obj, targetObj, task)
            % RUN_TASKONTARGET runs a task in the pipeline structure array in
            % TASK.
            %   It checks if the command TASK.FUNCNAME was already
            %   sucessfully performed by comparing it to the LOGBOOK from
            %   the PROTOCOL object. Also, it appends the information of
            %   the task on the object's LASTLOG .
            cd(targetObj.SaveFolder);
            LastLog = obj.ProtocolObj.createEmptyTable;
            % Fill out LOGTABLE:
            % Parse TARGETOBJ.SAVEFOLDER to get Subject/Acquisition/Recording
            % names:
            str = erase(targetObj.SaveFolder, obj.ProtocolObj.SaveDir);
            str = split(str, filesep); str = str(cellfun(@(x) ~isempty(x), str));
            for i = 1:length(str)
                LastLog(:,i) = str(i);
            end
            % Create evalStr:
            evalStr = ['targetObj.' task.funcStr];
            % Get input file info:
            if ~isempty(task.input)
                evalStr = strrep(evalStr, task.input , ['targetObj.' task.input]);
                task.funcStr = strrep(task.funcStr, task.input, ['targetObj.' task.input]);
                task.input = targetObj.(task.input);
                LastLog.InputFile_Path = {task.input};
                inputMetaData = strrep(task.input, '.dat', '_info.mat');
                mDt = matfile(inputMetaData);
                LastLog.InputFile_UUID = {mDt.fileUUID};
            end
            LastLog.ClassName = {class(targetObj)};
            LastLog.Job = {task.funcStr};
            % Check if Job was already performed:
            if ~obj.IgnoreLogBook
                b_isLogged = obj.checkInLogBook(obj.tmp_LogBook, LastLog);
            else
                b_isLogged = false;
            end
            
            if ~b_isLogged
                % Run the step:
                try
                    disp(['Running ' task.funcName '...']);
                    % Create Eval String:
                    eval(evalStr);
                    state = true;
                    LastLog.Messages = 'No Errors';
                catch ME
                    state = false;
                    LastLog.Messages = {getReport(ME)};
                end
                LastLog.Completed = state;
                LastLog.RunDateTime = datetime('now');
                targetObj.LastLog = [targetObj.LastLog; LastLog];
                obj.PipelineSummary = [obj.PipelineSummary; LastLog];
                obj.State = state;
                if LastLog.Completed
                    disp('Done!')
                else
                    disp('Failed!')
                end
            else
                disp([task.funcName ' Skipped!'])
                return
            end
        end
        
        function targetObj = getTargetObj(obj, targetIdx)
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
        end
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
                newParams.(paramsFields{i}) = params.(paramsFields{i});
            end
        end
        
        function funcStr = createFuncStr(~, params, funcName)
            % CREATEFUNCSTR checks for user-input parameters in PARAMS and
            % changes the FUNCINFO fields accordingly. Outputs the string
            % that will be passed in the EVAL function in
            % OBJ.RUN_TASKONTARGET and the VARS structure.
            
            funcStr = [funcName '('];
            % add input variable to the string:
            if ~isempty(params.input)
                funcStr = [funcStr num2str(params.input) ','];
            end
            % add optional variables to the string:
            if ~isempty(params.opts)
                optsField = fieldnames(params.opts);
                optsStr = '';
                for i = 1:length(optsField)
                    if isnumeric(params.opts.(optsField{i})) || islogical(params.opts.(optsField{i}))
                        str = ['''' optsField{i} ''',' num2str(params.opts.(optsField{i})) ','];
                    else
                        str = ['''' optsField{i} ''',''' params.opts.(optsField{i}) ''','];
                    end
                    optsStr= [optsStr str];
                end
                funcStr = [funcStr optsStr];
            end
            funcStr = [funcStr ');'];
            funcStr = strrep(funcStr, ',)', ')');
        end
        
        function newInput = pickAnInput(~,task)
            %PICKANINPUT prompts a LISTDLG to user-input of a function's
            % input when more than one output is possible.
            [indx,tf] = listdlg('PromptString', {'Select one of the inputs from', ['the function ' task.funcName]}, 'ListString',task.input, 'SelectionMode', 'single');
            if tf
                newInput = task.input{indx};
            else
                newInput = '';
            end
        end
        
        function out = askForFirstInput(obj, funcName)
            % ASKFORFIRSTINPUT creates an input prompt to get user to give the INPUT
            % of the first task in the pipeline.
            out = '';
            classes = fieldnames(obj.FunctionList);
            [indx,tf] = listdlg('PromptString', {'Select the Object containing', 'the input for the function :',  funcName},'ListString',classes, 'SelectionMode', 'single');
            if ~tf
                return
            end
            PropList = properties(classes{indx});
            idx = contains(PropList, '_file');
            PropList = PropList(idx);
            [indx,tf] = listdlg('PromptString', {'Select the File from ' classes{indx},  'as input for the function :',  funcName},'ListString',PropList, 'SelectionMode', 'single');
            if ~tf
                return
            end
            out = PropList{indx};
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
            
            for i = 2:height(obj.tmp_BranchPipeline)
                step = obj.tmp_BranchPipeline(i,:); 
                delete(step.InputFile_Path{:});
                delete(strrep(step.InputFile_Path{:}, '.dat', '_info.mat'));
                obj.tmp_BranchPipeline.InputFile_Path(i) = {'Deleted'};
                disp(['File ' step.InputFile_Path{:} ' deleted!'])
            end 
            
            % Update OBJ.PIPELINESUMMARY
            for i = 1:height(obj.tmp_BranchPipeline)
                idx = strcmp(obj.PipelineSummary.InputFile_UUID(:), obj.tmp_BranchPipeline.InputFile_UUID(i));
                obj.PipelineSummary(idx,:) = obj.tmp_BranchPipeline(i,:);
            end
            
        end
    
    end
end

