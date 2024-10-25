classdef PipelineManager < handle
    % PIPELINEMANAGER manages data processing pipelines.
    % This class allows the creation of analysis pipeline and manages the execution of
    % functions. In addition, it controls for failed / completed  steps in the pipeline.
    
    properties
        b_ignoreLoggedFiles  = false %(bool) If true, PIPELINEMANAGER will ignore identical jobs previously run.
        b_saveDataBeforeFail = false %(bool) If true, the more recent data ("current_data") in the pipeline will be saved to a file when an error occurs.
        b_overwriteFiles = false %(bool) If true, the pipeline will overwrite existing files with same names. If false, a "_n"umber will be appended to the new file name.
    end
    properties (SetAccess = {?DataViewer})
        % Structure array containing steps of the pipeline:
        pipe = struct('argsIn', {},'argsOut',{},'outFileName','',...
            'inputFileName', '', 'b_save2File', logical.empty, 'saveFileName',...
            '', 'opts',struct.empty,'opts_vals',struct.empty,...
            'opts_def',struct.empty ,'name','', 'inputFrom','', 'seq',[],'seqIndx',[],...
            'b_hasDataIn',false,'b_hasDataOut',false,'b_hasFileOut',false,'b_paramsSet',false);% !!If the fields are changed, please apply to the "set" method as well.
        funcList struct % structure containing the info about each function in the "fcnDir".
    end
    properties (SetAccess = private)
        % the same changes to the property's set method.
        fcnDir char % Directory of the analysis functions.
        PipelineSummary % Shows the jobs run in the current Pipeline
        current_seq = 0 % index of current sequence in pipeline.
        SaveFolderList % List of Save folders
        RawFolderList % List of Raw folders corresponding to the items in SaveFolderList.
    end
    properties (GetAccess = {?DataViewer})
        current_data % Data available in the workspace during pipeline.
        current_info % Structure with info about the functions ran on "current_data". This will be saved to the dataHistory file.
        current_saveFolder % Current path to save folder during pipeline execution.
        current_rawFolder % Current path to raw folder during pipeline execution.
        dataHistory % Local copy of the content from the dataHistory file in the data's SaveFolder.
        dv_inputFilename % Name of the input file (DataViewer only).
        b_inputFromDataViewer = false % TRUE if the first input from the pipeline comes from the DataViewer app.
        b_state logical % True if a task of a pipeline was successfully executed.
    end
    properties (Access = private)
        ProtocolSaveFolder char % Path to Protocols Save Folder. Used to locate PipelineConfig and Error Log folders.
        folderLog table % Pipeline LogTable stored in the "pipeLog.mat" file in the SaveFolder.
        maxLogRows = 100 % Maximum number of rows in the folderLog table. If exeeded, older entries will be erased.
        tmp_BranchPipeline % Temporarily stores LogBook from a Hierarchical branch.
        current_pipe % Pipeline currently running.
        current_outFile cell % List of file names created as output from some of the analysis functions.
        % It can be the name of an existing file, or
        % "outFile" for a function that creates a file
        % such as "run_ImagesClassification".
        h_wbItem % Handle of waitbar dialog showing the progress of the pipeline across protocol's objects.
        h_wbTask % Handle of waitbar dialog showing the progress of the pipeline in the current object.
        timeTag % timestamp used as "tag" for temporary files created during the pipeline execution.
        current_seqIndx = 0 % index of step in current sequence in pipeline.
        b_newSeq = false; % boolean to indicate if a new sequence was created in the previous step.
        b_pipeIsValid = false % TRUE, if the pipeline passed all validations and is ready to execution.        
    end
    
    methods
        % Constructor
        function obj = PipelineManager(SaveFolderList,varargin)
            % PIPELINEMANAGER Constructs an instance of this class
            %  Input:
            %   SaveFolderList(cell array of chars): List of Save folder
            %       paths containing .dat and .mat files.
            %   RawFolderList(cell array of chars): List of Raw folders
            %       associated with the save folders.
            %   Optionals:
            %   ProtocolSaveFolder (char): Path to "Protocol" folder. This
            %   will be used to save/load pipelines, and save pipeline
            %   error logs.
            
            p = inputParser;
            validationFun = @(x) (iscell(x) && ischar([x{:}])) || (ischar(x));
            addRequired(p,'SaveFolderList', validationFun);
            addOptional(p,'RawFolderList',{''}, validationFun);
            addOptional(p,'ProtocolSaveFolder',pwd,@ischar);
            parse(p,SaveFolderList, varargin{:});
            RawFolderList = p.Results.RawFolderList;
            obj.ProtocolSaveFolder = p.Results.ProtocolSaveFolder;
            
            % Further validate each Item from Save and Raw folder lists.
            if ischar(SaveFolderList);SaveFolderList = {SaveFolderList};end
            % criterion #1 - Both lists should have the same length:
            errID = 'umIToolbox:PipelineManager:MissingInput';
            errMsg = 'SaveFolderList and RawFolderList should have the same number of items!';
            assert(isequaln(length(SaveFolderList),length(RawFolderList)),errID,errMsg)
            % criterion #2 - All SaveFolders should exist:
            assert(all(isfolder(SaveFolderList)), 'One or more folders in SaveFolderList do not exist!');
            
            % For the Raw Folder list, replace missing folders with
            % "MISSING" string and raise a warning. Reason: sometimes, the raw data
            % is not available (different HD, PC etc.).
            b_missingFolders = ~isfolder(RawFolderList);
            if any(b_missingFolders)
                warning('A total of %d out of %d Raw Folder(s) do not exist! Beware! Functions that use this parameter will fail!',...
                    sum(b_missingFolders),length(b_missingFolders));
                RawFolderList(b_missingFolders) = {'MISSING'};
            end
            % Store folder lists
            obj.SaveFolderList = SaveFolderList;
            obj.RawFolderList = RawFolderList;
            % Create timestamp tag:
            obj.timeTag = datestr(datetime('now'),'_ddmmyyyyHHMMSS');

            if isdeployed
                [obj.fcnDir,~,~] = fileparts(which('funcTemplate.m'));
                a = load(fullfile(obj.fcnDir,'deployFcnList.mat'));
                obj.funcList = a.out; % Get the structure "out" created inside the function "umitFcnReader".
            else
                rootDir = erase(mfilename('fullpath'),['classes' filesep 'PipelineManager']);
                if isempty(rootDir)
                    error('Umitoolbox environment variable not found!')
                end
                obj.fcnDir = fullfile(rootDir, 'Analysis');
                obj.createFcnList;
            end
        end
        % SETTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.pipe(obj,pipe)
            % Pipeline structure setter. If pipe is empty, create an empty
            % structure containing tasks fields.
            
            if isempty(fieldnames(pipe))
                pipe = struct('argsIn', {},'argsOut',{},'outFileName','',...
                    'inputFileName', '','b_save2File', logical.empty, 'saveFileName',...
                    '', 'opts',struct.empty,'opts_vals',struct.empty,...
                    'opts_def',struct.empty ,'name','', 'inputFrom','',...
                    'seq',[],'seqIndx',[],'b_hasDataIn',false,'b_hasDataOut',false,...
                    'b_hasFileOut',false,'b_paramsSet',false);
            end
            % Check if all fields exist:
            if ~all(ismember(fieldnames(pipe),fieldnames(obj.pipe)))
                error('umIToolbox:PipelineManager:InvalidInput',...
                    'The pipeline structure provided is invalid!');
            end
            obj.pipe = pipe;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setOpts(obj,varargin)
            % SETOPTS opens an INPUTDLG for entry of optional variables
            % (OPTS) of methods in the Pipeline.
            
            p = inputParser;
            addRequired(p,'obj');
            addOptional(p,'funcName','',@ischar)
            addOptional(p,'seq',[],@isnumeric)
            parse(p,obj, varargin{:});
            funcName = p.Results.funcName;
            seq = p.Results.seq;
            
            % If pipeline is empty, abort:
            assert(~isempty(obj.pipe), 'Pipeline is empty. Create a pipeline first, then set the options!');
            % If no arguments are provided show the pipeline as a series of
            % sequences and ask input in command line:
            if any([isempty(funcName), isempty(seq)])
                fprintf('%s\nAvailable functions separated by sequence:\n',repmat('-',1,100))
                for ii = 1:obj.current_seq
                    thisSeq = obj.pipe(arrayfun(@(x) any(x.seq == ii) & numel(x.seqIndx) == 1,obj.pipe));
                    fprintf('----->> Sequence #%d :\n',ii);
                    fprintf('\t-->%s\n', thisSeq.name);
                end
            end
            if isempty(funcName)
                funcName = input('Type the name of the function to set parameters: ','s');
            end
            if isempty(seq)
                seq = input(['Type the sequence number of the function "' funcName '" :']);
            end
            % Check if the function exists in the selected pipeline sequence:
            idxName = strcmpi(funcName, {obj.pipe.name}); idxSeq = arrayfun(@(x) any(x.seq == seq),obj.pipe);
            if size(idxName,1)~= size(idxSeq,1)
                idxSeq = idxSeq';
            end
            idxFunc = idxName & idxSeq;
            if ~any(idxFunc)
                error(['The function "' funcName '" does not exist in pipeline sequence #' num2str(seq) '!']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Check if the function has optional parameters:
            if isempty(obj.pipe(idxFunc).opts)
                disp(['The function ' obj.pipe(idxFunc).name ' does not have any optional parameters.']);
                return
            end
            
            if isfield(obj.pipe(idxFunc), 'opts_vals')
                disp('Optional parameter values found!')
                currVals = cellfun(@(x) obj.pipe(idxFunc).opts.(x),fieldnames(obj.pipe(idxFunc).opts), 'UniformOutput',false);
                defVals = cellfun(@(x) obj.pipe(idxFunc).opts_def.(x),fieldnames(obj.pipe(idxFunc).opts_def), 'UniformOutput',false);
                listVals = cellfun(@(x) obj.pipe(idxFunc).opts_vals.(x),fieldnames(obj.pipe(idxFunc).opts_vals), 'UniformOutput',false);
                typeVals = {};
                for i = 1:length(listVals)
                    if isnumeric(listVals{i}) && numel(listVals{i}) > 2
                        typeVals{i} = 'numericArray'; % Array of numerical values.
                    elseif isnumeric(listVals{i}) && numel(listVals{i}) == 2
                        typeVals{i} = 'numericRange'; % 1x2 array of numerical values indicating lower and upper bounds.
                    elseif islogical(listVals{i})
                        typeVals{i} = 'logical'; % 1x2 array of logical values = [true;false];
                    elseif all(cellfun(@ischar, listVals{i})) && numel(listVals{i}) > 1 && size(listVals{i},1) < size(listVals{i},2)
                        typeVals{i} = 'charArray'; % Cell array of strings.
                    elseif all(cellfun(@ischar, listVals{i})) && numel(listVals{i}) > 1 && size(listVals{i},1) > size(listVals{i},2)
                        typeVals{i} = 'charArrayMultiSelect'; % Cell array of strings with multi-selection option.
                    else
                        typeVals{i} = 'mixArray'; % Cell array of strings and numbers.
                    end
                end
                out = buildInputDlg(obj.pipe(idxFunc).name,fieldnames(obj.pipe(idxFunc).opts),currVals,defVals,listVals,typeVals);
            end
            if isempty(out)
                return
            end
            % Save parameters to 'opts' structure:
            for i = 1:numel(fieldnames(obj.pipe(idxFunc).opts))
                obj.pipe(idxFunc).opts.(out{i,1}) = out{i,2};
            end
            obj.pipe(idxFunc).b_paramsSet = true; %
            disp(['Optional Parameters set for function : ' obj.pipe(idxFunc).name]);
        end
        
        function state = addTask(obj,func,varargin)
            % This method adds an analysis function, or task, to the
            % pipeline. Here, we can choose to save the output of a given
            % task as a .DAT file. ADDTASK will create a string containing
            % the task that will be evaluated during pipeline execution.
            % Inputs:
            %   func (str || char):  name or index of the analysis function
            %       contained in obj.funcList property.
            %   b_save2File (bool): Optional. True means that the output data
            %       from the function will be saved as a .DAT file.
            %   saveFileName(char): Optional. Name of the file to be saved.
            %       If not provided, the analysis function's default
            %       filename will be used.
            %   inputFrom (char): Optional. Name of the input function to
            %       the task. If you set this parameter, a new sequence of the
            %       pipeline will be created !
            %   inputFileName (char): Optional. Name of the input file to
            %       the task. If "inputFrom" is not provided, we assume that
            %       the input file comes from the hard drive.
            %   fromSeq (numeric): Optional. Sequence number of the input
            %       function set in "inputFrom" parameter. If not provided, we
            %       assume that the function comes from the sequence "1". This
            %       parameter is ignored if the input comes from the disk.
            % Output:
            %   state (bool): FALSE, if failed to add the task to the
            %   pipeline.
            
            p = inputParser;
            addRequired(p, 'func', @(x) ischar(x) || isnumeric(x));
            addOptional(p, 'b_save2File', false, @islogical);
            addOptional(p, 'saveFileName', '', @ischar);
            addParameter(p, 'inputFrom','', @ischar);
            addParameter(p,'inputFileName','',@ischar);
            addParameter(p,'fromSeq',1,@isPositiveIntegerValuedNumeric);
            parse(p,func, varargin{:});
            %
            currSeq = obj.current_seq;
            state = false;
            % Check if the function name is valid:
            idxFunc = strcmpi(p.Results.func, {obj.funcList.name});
            assert(any(idxFunc),['Operation aborted! The function "' p.Results.func '" does not exist!']);
            % Create "task" structure. This is the one that will be added
            % to the pipeline:
            task = obj.funcList(idxFunc).info;
            task.name = obj.funcList(idxFunc).name;
            % Add default values for fields used to save the data: These
            % values will be updated later:
            task.b_save2File = p.Results.b_save2File;
            task.saveFileName = p.Results.saveFileName;
            task.inputFileName = p.Results.inputFileName;
            task.b_paramsSet = false;
            task.inputFrom = '';
            % Control inputFrom:
            if ~isempty(p.Results.inputFrom)
                % Check if the task function needs data:
                assert(any(strcmpi(task.argsIn, 'data')), ['The function "' task.name '" does not have "data" as input. Operation aborted!'])
                if ~strcmpi(p.Results.inputFrom, '_LOCAL_')
                    % Find input function in the existing pipeline:
                    idxInputFcn = ( strcmpi(p.Results.inputFrom, {obj.pipe.name}) & arrayfun(@(x) any(x.seq == round(p.Results.fromSeq)),obj.pipe)' );
                    assert(any(idxInputFcn), 'Input function not found! Operation aborted!');
                    % Save index of input function to the current task:
                    task.inputFrom = find(idxInputFcn); % Use pipeline index instead of name to point to input function.
                    % Check if the input function generates data:
                    assert(any(ismember({'outData','outFile'},obj.pipe(idxInputFcn).argsOut)),...
                        ['The function "' task.inputFrom '" does not have any "data" as output. Operation aborted!']);
                    % Check if the inputFileName exists in the list of
                    % outFileNames from functions with "outFile":
                    if ~isempty(task.inputFileName) && any(strcmpi('outFile',{obj.pipe(idxInputFcn).argsOut}))
                        assert(any(strcmpi(task.inputFileName, {obj.pipe(idxInputFcn).argsOut})),'Input file name is invalid!');
                    end
                else
                    % Set "inputFrom" to zero when the input comes from the
                    % disk.
                    task.inputFrom = 0;
                end
            else
                % Control for cases where the User sets an input file name
                % without setting the "inputFrom" parameter. In this case,
                % we force set the "inputFrom" to the disk.
                if ~isempty(task.inputFileName)
                    task.inputFrom = 0;
                end
            end
            %%% SPECIAL CASE - DataViewer
            % For the first step of the pipeline avoid asking for an input
            % file:
            if obj.b_inputFromDataViewer&& obj.current_seq == 0
                if isempty(task.inputFrom) || task.inputFrom ~= 0
                    task.inputFrom = -1;% Set inputFrom to "-1" the data from DataViewer.
                else
                    % The task.inputFrom field will be zero when generating a script
                    % from DataViewer GUI.
                    task = obj.setInput(task);
                end
                if task.dependency
                    obj.current_seqIndx = 0;
                    task = obj.addDependency(task);
                    task.seq = obj.current_seq;
                    task.seqIndx = obj.current_seqIndx + 1;w
                    task.inputFrom = length(obj.pipe);
                else
                    obj.current_seq = 1;
                    obj.current_seqIndx = 1;
                    task.seq = obj.current_seq; task.seqIndx = obj.current_seqIndx;
                end
                if ~isempty(obj.current_data)
                    [~,file,ext] = fileparts(obj.dv_inputFilename);
                    task.inputFileName = [file ext]; % Get name of the input file in DataViewer.
                end
                % Remove dependency field from "task")
                task = rmfield(task, 'dependency');
                %
                % Add task to pipeline:
                obj.pipe = [obj.pipe; task];
                fprintf('Added "%s" to sequence #%d of the pipeline.\n',task.name, obj.current_seq)
                state = true;
                return
            end
            % Check if the new function needs another one to work:
            if task.dependency
                task = obj.addDependency(task);
            end
            % Remove dependency field from "task")
            task = rmfield(task, 'dependency');
            %%%------------------------------------------------------------
            % Add branch and sequence index to the function:
            task = obj.setInput(task);
            if isempty(task)
                % If failed, return
                obj.current_seq = currSeq; % Reset current sequence
                return
            end
            % Control for duplicate functions in the same sequence:
            assert(~any(strcmpi(task.name,{obj.pipe(arrayfun(@(x) any(x.seq == obj.current_seq),obj.pipe)).name})), ...
                ['Operation Aborted! The function "' task.name '" already exists in the pipeline!'])
            % Perform checks on save options:
            if task.b_save2File
                % Check if the saveFileName is correct:
                [~,~,ext] = fileparts(task.outFileName);
                if ~any(strcmpi('outData', task.argsOut))
                    warning(['Cannot save output to .DAT file for the function'...
                        ' "%s" \nbecause it doesn''t have any data as output!'], task.name);
                    task.b_save2File = false;
                    task.saveFileName= '';
                elseif isempty(task.saveFileName)
                    % Save as default:
                    task.saveFileName = task.outFileName;
                elseif ~endsWith(task.saveFileName, ext)
                    % Append file exension:
                    task.saveFileName = [task.saveFileName, ext];
                end
            end
            % Add task to pipeline:
            obj.pipe = [obj.pipe; task];
            fprintf('Added "%s" to sequence #%d of the pipeline.\n',task.name, obj.current_seq)
            state = true;
        end
        %%%%% TO BE DONE --------------------------------------------------
        %         function state = rmTask(obj,func,varargin)
        %             % RMTASK removes a given function from the existing pipeline and
        %             % updates the remaining steps in order to maintain the pipeline
        %             % workflow.
        %             % Inputs:
        %             %   func (str || char):  name or index of the analysis function
        %             %       contained in obj.funcList property.
        %             %   fromSeq (positive integer): Optional. Sequence number of the input
        %             %       function set in "inputFrom" parameter. If not provided, we
        %             %       assume that the function comes from the sequence "1". This
        %             %       parameter is ignored if the input comes from the disk.
        %             % Output:
        %             %   state (bool): FALSE, if failed to add the task to the
        %             %   pipeline.
        %
        %             p = inputParser;
        %             addRequired(p,'func',@(x) ischar(x) || isnumeric(x));
        %             addParameter(p,'fromSeq',1,@isPositiveIntegerValuedNumeric);
        %             parse(p,func,varargin{:});
        %             %
        %             state = false;
        %             % Check if the function exists in the pipeline sequence:
        %             idxFunc = ( strcmpi(p.Results.func,{obj.pipe.name}) && arrayfun(@(x) any(x.seq == p.Results.fromSeq),obj.pipe) );
        %             assert(any(idxFunc),['Operation aborted! The function "' p.Results.func '" does not exist!']);
        %
        %         end
        %%%%% -------------------------------------------------------------
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
            
            str = sprintf('Pipeline Summary:\n\n');
            for i = 1:length(obj.pipe)
                str =  [str sprintf('--->> Step # %d <<---\n', i)];
                if isempty(obj.pipe(i).opts)
                    opts = {'none'; 'none'};
                else
                    fn = fieldnames(obj.pipe(i).opts)';
                    vals = {};
                    for j = 1:length(fn)
                        val = obj.pipe(i).opts.(fn{j});
                        if isnumeric(val) || islogical(val)
                            vals{j} = num2str(val);
                        elseif iscell(val) && ischar(val{1})
                            vals{j} = strjoin(val, ', ');
                        else
                            vals{j} = val;
                        end
                    end
                    opts = [fn;vals];
                end
                txt = sprintf('Function name : %s\nOptional Parameters:\n',...
                    obj.pipe(i).name);
                str = [str, txt, sprintf('\t%s : %s\n', opts{:})];
                if obj.pipe(i).b_save2File
                    str = [str, sprintf('Data to be saved as : "%s"\n', obj.pipe(i).saveFileName)];
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
        
        function run_pipeline(obj)
            % RUN_PIPELINE runs the tasks in OBJ.PIPE
            
            % Reset PipelineSummary table:
            obj.PipelineSummary = table();
            
            if ~obj.b_pipeIsValid
                % Force saving sequences' last steps (Except for DataViewer)
                obj.validatePipeline; % Run pipeline validation
            end
            % Recheck:
            if ~obj.b_pipeIsValid
                warning('Pipeline execution aborted! Pipeline is invalid.')
                return
            end

            % Initialize waitbars:
            obj.setWaitBar('Initialize')
            
            for ii = 1:length(obj.SaveFolderList)
                % Clear current data, metaData and File List for each item in te pipeline:
                if ~obj.b_inputFromDataViewer
                    obj.current_data = []; obj.current_info = []; obj.current_outFile = {};
                end
                obj.b_state = true;
                % Update current save and raw folders:
                obj.current_saveFolder = obj.SaveFolderList{ii};
                obj.current_rawFolder = obj.RawFolderList{ii};
                obj.loadFolderLog(obj.SaveFolderList{ii});
                % Update Data History from current save folder:
                obj.dataHistory = obj.loadDataHistory(obj.current_saveFolder);
               
                % Update waitbars:
                obj.setWaitBar('UpdateItem', ii, length(obj.SaveFolderList));
                fprintf([repmat('-',1,50),'\n']);
                fprintf('Save Folder : %s\n', obj.current_saveFolder);
                fprintf('Raw Folder : %s\n', obj.current_rawFolder);
                % Run one sequence at a time:
                for jj = 1:max([obj.pipe.seq],[],'all')
                    obj.current_seq = jj;
                    % Get current pipeline sequence:
                    thisSeq = obj.pipe(arrayfun(@(x) any(x.seq == obj.current_seq), obj.pipe));
                    % Skip pipeline steps:
                    if ~obj.b_ignoreLoggedFiles
                        [thisSeq, skippedFcns,selFile] = obj.skipSteps(thisSeq);
                        if isempty(thisSeq)
                            % When all sequence is skipped:
                            fprintf('All steps skipped from the current sequence!\n')
                            % Initialize empty Log for current object:
                            % Add class name to table:
                            txt = sprintf('The file %s contains the following steps:',selFile);
                            txt = [txt,sprintf('\n%s', skippedFcns{:})];%#ok
                            dummyTask = struct('name',['Skipped sequence # ' num2str(obj.current_seq)],...
                                'inputFileName',selFile,'saveFileName',selFile);
                            obj.updateFolderLog(dummyTask,datetime('now'),true,txt);                             
                        elseif ~isempty(skippedFcns)
                            % When some steps are skipped:
                            fprintf('The following steps will be skipped:\n')
                            fprintf('\t"%s"\n',skippedFcns{:});
                        end
                    end
                    % Run pipeline sequence in each target object:
                    for kk = 1:length(thisSeq)
                        task = thisSeq(kk);
                        obj.setWaitBar('UpdateTask', kk/length(thisSeq),'taskName',task.name);
                        fprintf('Running task # %d/%d ----->>>>>\n',kk,length(thisSeq));
                        % Execute step on target object:
                        obj.run_taskOnTarget(task);
                        % Control for Pipeline cancelling by User:
                        if getappdata(obj.h_wbItem, 'b_abortPipe')
                            % Delete waitbars and abort inner loop:
                            delete([obj.h_wbItem, obj.h_wbTask])
                            break
                        end
                        % If execution failed, abort.
                        if ~obj.b_state
                            break
                        end
                        % This pause is here to allow the WaitBar to update
                        % during the execution of this method.
                        pause(.001);
                    end
                    fprintf([repmat('-',1,50),'\n']);
                    % Remove temporary files with appended with the "timeTag":
                    obj.deleteTemporaryFiles(obj.current_saveFolder);
                    % Abort outer loop if user cancels pipeline or a pipeline execution failed:
                    if ~ishandle(obj.h_wbItem) || ~obj.b_state
                        break
                    end
                end
                % Abort outer loop if user cancels pipeline:
                if ~ishandle(obj.h_wbItem)
                    break
                end
                % Save folderLog file:
                obj.saveFolderLog (obj.current_saveFolder);                                
            end
                       
            % Show Pipeline Summary in command window:
            disp(obj.PipelineSummary);
            % Delete progress bars:
            delete([obj.h_wbItem, obj.h_wbTask]);
            % Request User permission to save error report:
            obj.genErrorReport;
        end
        
        function validatePipeline(obj, varargin)
            % VALIDATEPIPELINE checks if the data generated in the pipeline
            % is saved to the disk. If the user forgets to save the end
            % points of the pipeline sequences, this function will force the
            % saving of the data using the function's default ouput
            % filename.
            if nargin > 1
                % Change function behaviour when using it for script generation
                b_isGenScript = varargin{:};
            else
                b_isGenScript = false;
            end
            
            obj.b_pipeIsValid = false;
            % Check each sequence from bottom to top until a "data" output
            % is found, then force it to be saved:
            
            lastFileNameList = {''};
            for ii = 1:obj.current_seq
                % Isolate sequence:
                idxSeq = arrayfun(@(x) any(x.seq == ii),obj.pipe);
                thisSeq = obj.pipe(idxSeq);
                for jj = length(thisSeq):-1:1
                    if ~any(strcmpi(thisSeq(jj).argsOut, 'outData'))
                        % Go to the previous step if the last one doesn't
                        % generate data as output.
                        continue
                    end
                    % Force saving last step of the sequence (Except for DataViewer):
                    if ~thisSeq(jj).b_save2File && ( ~obj.b_inputFromDataViewer || b_isGenScript )
                        thisSeq(jj).b_save2File = true;
                        thisSeq(jj).saveFileName = thisSeq(jj).outFileName;
                    end
                    
                    if obj.b_inputFromDataViewer && ~thisSeq(jj).b_save2File
                        % DATAVIEWER: Skip duplicate save file names for
                        % steps that will not be saved.
                        break
                    end
                    % Check for output files with the same name
                    idxDuplicate = strcmpi(thisSeq(jj).saveFileName,lastFileNameList);
                    if any(idxDuplicate)
                        % If a filename from a sequence already exists,
                        % ask user to change it. Otherwise, the data
                        % from one sequence will be OVERWRITTEN!!
                        [~,defName,ext] = fileparts(lastFileNameList{idxDuplicate});
                        newFileName  = [defName '_' num2str(sum(idxDuplicate)+1), ext]; % Append sequence index to new filename
                        qd = changeFileName(thisSeq(jj).name, lastFileNameList{idxDuplicate});
                        waitfor(qd);
                        thisSeq(jj).saveFileName = newFileName;
                    end
                    lastFileNameList = [lastFileNameList,{thisSeq(jj).saveFileName}];%#ok
                    break
                end
                % Save changes to pipeline:
                obj.pipe(idxSeq) = thisSeq;
                obj.b_pipeIsValid = true;
            end
            disp('Pipeline is ready to run.');
            %%%%%--Local function ------------------------------------------
            function uf = changeFileName(fcnName, fileName)
                % This function creates an input dialog so the User can type in
                % the name of the file to be saved
                uf = uifigure('Name','Duplicate file names found!','NumberTitle',...
                    'off', 'CreateFcn',{@movegui,'center'}, 'Position',[0 0 470 150],...
                    'Visible','off', 'CloseRequestFcn',@okBtnPushed);
                msg = sprintf(['A file with name "%s" is already set to be saved at the end of the pipeline.\n' ...
                    'Please, type a different file name for the function "%s":'], fileName,fcnName);
                uilabel(uf, 'Position',[10 80 460 70], 'Text', msg);
                ef = uieditfield(uf,'Position',[10 60 300 30],'Value',newFileName);
                uibutton(uf,'Position',[10,10,80,34],'Text','Ok', 'ButtonPushedFcn',@okBtnPushed);
                uf.Visible = 'on';
                % OK button Callback
                function okBtnPushed(source,~)
                    newFileName = ef.Value;
                    [~,name,extension] = fileparts(newFileName);
                    if isempty(extension)
                        % Force original extension:
                        newFileName = [name,ext];
                    end
                    % Close figure;
                    if isa(source,'matlab.ui.Figure')
                        delete(source);
                    else
                        delete(source.Parent)
                    end
                end
            end
        end
        
        function savePipe(obj,pipeFile)
            % SAVEPIPE saves the structure OBJ.PIPE in a .PIPE file in the
            % folder PIPELINECONFIGFILES inside the SAVEDIR of OBJ.PROTOCOLOBJ.
            % If this function is called from DataViewer, it saves in
            % DataViewer's folder set by the user.
            
            if ~obj.b_pipeIsValid
                obj.validatePipeline;
            end
            if ~obj.b_pipeIsValid
                error('Failed to validate pipeline. Save Pipeline aborted!')
            end
            [path, filename,ext] = fileparts(pipeFile);
            if ~strcmpi(filename,'.pipe') || isempty(ext)
                ext = '.pipe'; % Force extension.
            end
            % If the user did not enforce a path, use ProtocolSaveFolder:
            if isempty(path)
                if ~obj.b_inputFromDataViewer
                    path = fullfile(obj.ProtocolSaveFolder, 'PipeLineConfigFiles');
                    [~,~] = mkdir(path);
                else
                    path = obj.ProtocolSaveFolder;
                end
            end
            pipeStruct = obj.pipe;
            save(fullfile(path,[filename ext]), 'pipeStruct','-mat');
            disp(['Pipeline saved as "' filename ext '" in ' path]);
        end
        
        function loadPipe(obj,pipeFile)
            % LOADPIPE loads the structure PIPE inside FILENAME and assigns
            % it to OBJ.PIPE property.
            % Input:
            %   pipeFile(char): full path to the .JSON file containing the
            %   pipeline config.
            
            % Read Pipeline Config file:
            a = load(pipeFile,'-mat');
            new_pipe = a.pipeStruct;
            % Reset current pipeline:
            obj.reset_pipe;
            % Check if all functions listed in "new_pipe" exist:
            idx = ismember({new_pipe.name},{obj.funcList.name});
            if ~all(idx)
                h = errordlg(['The following functions do not exist:' sprintf('\n"%s"',new_pipe(~idx).name)], ...
                    'Failed to load pipeline!');
                waitfor(h);
                return
            end            
            obj.pipe = new_pipe;
            % Update current sequence index:
            obj.current_seq = max([obj.pipe.seq],[],'all');
            disp('Pipeline loaded!')
            % Validate loaded pipeline:
            obj.validatePipeline;
        end
        
        function reset_pipe(obj)
            % This function erases the pipe property and resets the funcList
            % property to default parameter values.
            obj.pipe = struct();
            % Reset some properties:
            if ~obj.b_inputFromDataViewer
                obj.current_data = [];
            end
            obj.current_outFile = {};
            obj.current_seq = 0; obj.current_seqIndx = 0;
            obj.b_pipeIsValid = false;
            %
            disp('Pipeline erased!')
        end
        
        function generateScript(obj, filename)
            %   GENERATESCRIPT Generates a MATLAB script file from the pipeline structure.
            %   Generates a script file that executes the steps defined in the pipeline.
            %
            %   Input:
            %   - filename: The name of the script file to be generated.
            %
            %   This method iterates through the pipeline steps, generates MATLAB script
            %   code for each step, and appends it to the output script file specified
            %   by 'filename'. The script includes loading input data (if required),
            %   setting optional parameters, executing the pipeline step function,
            %   and saving output data to files (if required). Additionally, it may
            %   include a section to delete temporary files.
            %
            
            % Check if the pipeline is empty, and if so, return.
            if isempty(obj.pipe)
                return
            end
            
            % Check if the pipeline is properly built:
            obj.validatePipeline(true);
            
            % Create the initial text for the script.
            txt = sprintf(['%%%% Pipeline script.\n%% Script generated by PipelineManager class @'...
                '%s\n\n%% clearvars;%% Clear workspace (commented by default)\n%%'...
                'Here is a summary of the pipeline:\n'], datestr(now(), 'HH:MM:ss dd-mm-yyyy'));
            % Add pipeline summary to script docstring:
            summaryTxt = '';
            tmp = obj.showPipeSummary;
            idxNL = 1+regexp(tmp,'\n');
            idxNL = [1 idxNL(1:end-1)];
            for ii = 2:length(idxNL)
                summaryTxt = [summaryTxt  sprintf('%%\t\t') tmp(idxNL(ii-1):idxNL(ii)-1)];%#ok
            end
            txt = [txt summaryTxt];
            % Initialize working folder to the current directory.
            txt = [txt sprintf('%%%% Set variables:\nFolder = pwd; %% By Default, set the current folder as the working folder.\n%%%% Pipeline execution\n\n')];
            
            % Loop through each step in the pipeline.
            curSeq = 0;
            for ii = 1:length(obj.pipe)
                if obj.pipe(ii).seq(1) ~= curSeq
                    txt = [txt sprintf('%% Sequence %d:\n',obj.pipe(ii).seq(1))];
                    curSeq = obj.pipe(ii).seq(1);
                end
                % Generate comments for the current step.
                txt = [txt sprintf('%% Execute "%s":\ndisp(''[Running...] %s'');\n', obj.pipe(ii).name,obj.pipe(ii).name)];%#ok
                
                % Add data loading string if input data file is specified.
                if ~isempty(obj.pipe(ii).inputFileName) && ~strcmpi(obj.pipe(ii).inputFileName, 'data')
                    % If the function requires a file as input, create a loading string.
                    if endsWith(obj.pipe(ii).inputFileName, '.dat')
                        txt = [txt sprintf('[data, metaData] = loadDat(''%s''); %% Load input data\n', obj.pipe(ii).inputFileName)];%#ok
                    else
                        txt = [txt sprintf('data = load(%s); %% Load input data\n', obj.pipe(ii).inputFileName)];%#ok
                    end
                end
                
                % Set input and output argument names.
                argsIn = obj.pipe(ii).argsIn;
                argsIn = replace(argsIn, {'RawFolder', 'SaveFolder'}, {'Folder', 'Folder'});
                argsIn(strcmpi(argsIn,'object')) = [];
                
                % Add optional parameters if applicable.
                optsStr = '';
                fcnStr = '';
                if obj.pipe(ii).b_paramsSet
                    fn = fieldnames(obj.pipe(ii).opts);
                    for jj = 1:length(fn)
                        if ischar(obj.pipe(ii).opts.(fn{jj}))
                            val = ['''' obj.pipe(ii).opts.(fn{jj}) ''''];
                        elseif iscell(obj.pipe(ii).opts.(fn{jj})) && ischar([obj.pipe(ii).opts.(fn{jj}){:}])
                            val = '{{';
                            for kk = 1:length(obj.pipe(ii).opts.(fn{jj}))
                                val = [val '''' obj.pipe(ii).opts.(fn{jj}){kk} '''', ','];
                            end
                            val(end:end+1) = '}}';
                        else
                            val = ['[' strjoin(arrayfun(@num2str, obj.pipe(ii).opts.(fn{jj}), 'UniformOutput', false), ';') ']'];
                        end
                        optsStr = [optsStr ',''' fn{jj} ''',' val];
                    end
                    optsStr = ['opts = struct(' strip(optsStr, 'left', ',') '); % Optional parameters'];
                else
                    argsIn(strcmpi(argsIn,'opts')) = [];
                end
                
                argsOut = obj.pipe(ii).argsOut;
                argsOut = replace(argsOut, 'outData', 'data');
                
                % Create function string for the current step.
                if isempty(argsOut)
                    fcnStr = [fcnStr ';' obj.pipe(ii).name '(' strjoin(argsIn, ',') ');'];%#ok
                elseif numel(argsOut) == 1
                    fcnStr = [fcnStr ';' strjoin(argsOut, ',') ' = ' obj.pipe(ii).name '(' strjoin(argsIn, ',') ');'];%#ok
                else
                    fcnStr = [fcnStr ';' '[' strjoin(argsOut, ',') '] = ' obj.pipe(ii).name '(' strjoin(argsIn, ',') ');'];%#ok
                end
                fcnStr = strip(fcnStr, 'left', ';');
                txt = [txt sprintf('%s\n%s\n\n', optsStr, fcnStr)];%#ok
                
                % File Saving section:
                if obj.pipe(ii).b_save2File
                    % Save data to file.
                    if endsWith(obj.pipe(ii).saveFileName, '.dat')
                        txt = [txt sprintf('saveDat(fullfile(Folder, ''%s''), data); %% Save data to .DAT file "%s" \n', obj.pipe(ii).saveFileName, obj.pipe(ii).saveFileName)];%#ok
                    else
                        txt = [txt sprintf('save(fullfile(Folder, %s), ''-struct'', ''data'', ''-v7.3''); %% Save data to .MAT file "%s" \n', obj.pipe(ii).saveFileName, obj.pipe(ii).saveFileName)];%#ok
                    end
                end
                txt = [txt sprintf('disp(''[Completed] %s'');\n',obj.pipe(ii).name)];
            end
            
            % Special case:
            % Create a section to delete temporary files created inside the script.
            idx = contains({obj.pipe.saveFileName}, num2str(obj.timeTag));
            if any(idx)
                % Create deletion section.
                delStr = sprintf('%%%% Delete temporary files\n%% Delete .dat/mat files created during the execution of this script.\n');
                delStr  = [delStr sprintf('disp(''Deleting temporary files...'');\n')];
                fNames = {obj.pipe.saveFileName};
                fNames(~idx) = [];
                for ii = 1:length(fNames)
                    if endsWith(fNames{ii}, '.dat')
                        delStr = [delStr sprintf('delete(''%s'');  delete(''%s'');\n', fNames{ii}, strrep(fNames{ii}, '.dat', '.mat'))];%#ok
                    else
                        delStr = [delStr sprintf('delete(''%s'');\n', fNames{ii})];%#ok
                    end
                end
                txt = [txt delStr];
            end
            
            % Add a closing comment to the script.
            txt = [txt sprintf('%%%%\ndisp(''Script execution completed!'');\n%%%%%%%%%%%%%%%%%%%%%%%% END OF FILE %%%%%%%%%%%%%%%%%%%%%%%%')];
            
            % Save the generated script text to the specified file.
            fid = fopen(filename, 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            
            % Get the folder containing the generated script.
            folder = fileparts(filename);
            
            % Display a message indicating that the script has been generated.
            disp('Script generated!')
            
            % If the folder is empty, set it to the current directory.
            if isempty(folder)
                folder = pwd;
            end
            
            % Open the folder containing the generated script.
            openFolder(folder);
        end
        
        %%%%%%--Pipeline Visualization  -----------------------------------
        function fH = drawPipe(obj,varargin)
            % DRAWPIPE creates a Directed Acyclic Graph (DAG) representation
            % of the pipeline.
            % Inputs:
            %   fH (handle, optional): handle to the figure where the
            %       DAG will be plotted. The figure must be created with the
            %       command "figure". "UIfigure" not supported!
            % Output:
            %   fH (handle): handle to the figure with plotted DAG.
            p = inputParser;
            addRequired(p,'obj')
            addOptional(p,'fH',[],@(x) isa(x,'matlab.ui.Figure') | isempty(x))
            parse(p,obj, varargin{:});
            if isempty(p.Results.fH)
                % Create figure when no axis is provided:
                fH = figure('Name','Pipeline Visualization', 'MenuBar','none', 'Toolbar','none','NumberTitle','off');
            else
                fH = p.Results.fH;
                clf(fH)
            end
            % Create local copy of the pipeline. We will change some info
            % just for the sake of plotting.
            pp = obj.pipe;
            % Abort if the pipeline is empty.
            if isempty(pp)
                return
            end
            if obj.b_inputFromDataViewer
                dskName = 'DataViewer';
            else
                dskName = 'Disk';
            end
            % Get original figure position:
            origPos = fH.Position;
            % GUI elements' paramters:
            bounds = [.2 .1 .8 .9]; % Normalized to facilitate editing.
            xSpacing = 50; % Minimal distances between edges of buttons on X;
            ySpacing = 60; % Idem in Y.
            btnHeight = 35; % Button Height in points.
            btnFontSize = 11; % Button Font size
            arrowHeadSz = 7; % Arrow head size in points;
            % Button colors:
            myRed = [.9 0 0];
            myGreen = [0 .85 0];
            myGray = [.92 .92 .92];
            % Create UIContext Menu for extra options:
            cm = uicontextmenu(fH);
            uimenu('Parent',cm,'Label','Save to file','Callback',{@chooseSaveFileName,obj}); % Menu to change save filename.
            % Create panel to be able to lock the figure during PushButton
            % Callback execution:
            pan = uipanel('Parent', fH, 'Position',[0 0 1 1], 'Title','Setting parameter...','Visible','off');
            % Create "Disk" button at the middle of the figure:
            diskBtn = uicontrol(fH,'Style','pushbutton','String',dskName,'Enable','off', 'FontSize',btnFontSize, 'Tag','origin');
            diskBtn.Position([3 4]) = [diskBtn.Extent(3) btnHeight];
            %             diskBtn.Position(4) = btnHeight;
            % Re-calculate sequence indices for 2+ sequences for a better
            % display:
            if any([pp.seq]>1)
                for ii = 2:max([pp.seq])
                    seqPos = find(arrayfun(@(x) any(x.seq == ii),pp));
                    if pp(seqPos(1)).inputFrom == 0
                        continue
                    end
                    for jj = 2:length(seqPos)
                        seqIndxPos = pp(seqPos(jj)).seq == ii;
                        pp(seqPos(jj)).seqIndx(seqIndxPos) = max([pp(seqPos).seqIndx]) + 1;
                    end
                end
            end
            % Create buttons for each function:
            for ii = 1:length(pp)
                btnArr(ii) = uicontrol(fH,'Style','pushbutton','String',pp(ii).name, 'FontSize',btnFontSize);
                % Add Push button callback:
                btnArr(ii).UserData = ii; % Store pipeline index in UserData;
                btnArr(ii).Callback = {@callSetOpts,obj,pan}; % Call setOpts to set function's paramerets
                btnArr(ii).Position([3,4]) = [btnArr(ii).Extent(3)+10 btnHeight];% Avoid word wrapping
                btnArr(ii).BackgroundColor = fH.Color; % Make button "invisible" the color will be given by the CData property.
                if pp(ii).b_hasDataOut & any(~contains(pp(ii).argsOut, 'outFile','IgnoreCase',true))
                    % Add context menu for functions with output data.
                    btnArr(ii).UIContextMenu = cm;
                end
                btnArr(ii).Tooltip = obj.genToolTipTxt(pp(ii));
            end
            % Make all buttons the same width
            maxW = max(arrayfun(@(x) x.Position(3), btnArr));
            arrayfun(@(x) set(x, 'Position',[x.Position(1), x.Position(2), maxW, x.Position(4)]),btnArr);
            % Improve buttons appearance:
            myGreen = obj.prettyfyBtn(myGreen,btnArr(1).Position([3 4]), fH.Color(1)); % Here, we assume that the figure's color is a tone of gray.
            myRed = obj.prettyfyBtn(myRed,btnArr(1).Position([3 4]), fH.Color(1));
            myGray = obj.prettyfyBtn(myGray,btnArr(1).Position([3 4]), fH.Color(1));
            for ii = 1:length(pp)
                % Change button background color if parameters were set:
                if isempty(obj.pipe(ii).opts)
                    btnArr(ii).CData = myGray; % Gray. No parameters.
                elseif obj.pipe(ii).b_paramsSet
                    btnArr(ii).CData = myGreen; % Green. Parameters already set.
                else
                    btnArr(ii).CData = myRed; % Red. Parameters not set.
                end
            end
            % Check if the figure is sufficiently large to accomodate all buttons
            % without overlap:
            btnX = arrayfun(@(x) x.Position(3), btnArr);
            seqIndxList = arrayfun(@(x) x.seqIndx(1,1),pp);
            newFigX = fH.Position(3);
            for ii = 1:numel(unique(seqIndxList))
                idx = seqIndxList == ii;
                sumX = ceil(sum(btnX(idx)) + xSpacing*(sum(idx)+1));
                if sumX > newFigX
                    newFigX = sumX;
                end
            end
            fH.Position(3) = newFigX;
            seqList = unique([pp.seq]);
            % Do the same in Y.
            nLvls = max(seqIndxList) + 1;
            minYsize = nLvls*btnHeight + (nLvls-1)*ySpacing;
            if minYsize > fH.Position(4)
                fH.Position(4) = minYsize;
            end
            % Transform to pixels:
            bounds = [bounds([1 3])*fH.Position(3) bounds([2 4])*fH.Position(4)];
            bounds = bounds([1 3 2 4]);
            items = [diskBtn, btnArr];
            % Get list of sequences and sequence indices for each function:
            
            % Distribute horizontally:
            ctrX = zeros(1,length(items));
            ctrX(1) = bounds(1)+(bounds(3)-bounds(1))/2; % Put "Disk" on the middle;
            if max(seqList) == 1
                % If there is only one sequence, put it in the middle of
                % the figure
                ctrX(2:end) = ctrX(1);
            else
                % Spread sequences evenly across the figure width.
                seqPairs = zeros(length(pp),2);
                for ii = 1:length(pp)
                    % Check if the step input is from a different sequence,
                    % if so, ensure that there are no arrow crossings by
                    % rearranging the x positions:
                    if isempty(pp(ii).inputFrom) || pp(ii).inputFrom == 0 || pp(ii).inputFrom == -1
                        seqPairs(ii,:) = [pp(ii).seq(1,1), pp(ii).seq(1,1)];
                    else
                        seqPairs(ii,:) = [pp(pp(ii).inputFrom).seq(1,1), pp(ii).seq(1,1)];
                    end
                end
                % Minimize the number of arrows crossing sequences:
                seqPairs = unique(seqPairs,'rows');
                crossSeqs = seqPairs(seqPairs(:,1)~=seqPairs(:,2),:);
                permIndx = seqList;
                if ~isempty(crossSeqs)
                    allPerms = perms(seqList);
                    crossSum = zeros(size(allPerms,1),1);
                    for ii = 1:size(allPerms,1)
                        [~,idx_target] = ismember(crossSeqs(:,2), allPerms(ii,:));
                        [~,idx_source] = ismember(crossSeqs(:,1), allPerms(ii,:));
                        crossSum(ii) = sum(abs(idx_target - idx_source));
                    end
                    % Find sequence with minimal number of crossings:
                    permIndx = allPerms(find(crossSum == min(crossSum),1,'last'),:);
                end
                % Get sequence X positions:
                xPos = round(linspace(bounds(1), bounds(3),max(seqList)));
                xPos = xPos(permIndx);
                Xdict = containers.Map(seqList,xPos);
                % Update ctrX, except "Disk":
                for ii = 1:length(pp)
                    ctrX(ii+1) = Xdict(pp(ii).seq(1,1));
                end
            end
            % Distribute vertically:
            nLvl = max([pp.seqIndx]);
            Y = fliplr(round(linspace(bounds(2), bounds(4),nLvl+1)));
            Ydict = containers.Map(1:nLvl,Y(2:end));
            ctrY = zeros(1,length(items));
            ctrY(1) = Y(1);
            for ii = 1:length(pp)
                indx = pp(ii).seqIndx(1,1);
                ctrY(ii+1) = Ydict(indx);
            end
            % Fix X and Y positions so the buttons will be centered:
            for ii = 1:length(items)
                items(ii).Position(1) = ctrX(ii) - (items(ii).Position(3)/2);
                items(ii).Position(2) = ctrY(ii) - (items(ii).Position(4)/2);
            end
            % Add arrows between buttons:
            arrYsource = arrayfun(@(x) x.Position(2), items);
            arrYtarget = arrayfun(@(x,y) x + y.Position(4)/2,ctrY,items);
            % Draw arrows and texts:
            for ii = seqList
                seq = find(arrayfun(@(x) any(x.seq == ii),pp));
                if pp(seq(1)).inputFrom == 0 |  pp(seq(1)).inputFrom == -1 %#ok
                    % For when the input comes from the "Disk" or "DataViewer":
                    an = annotation('arrow',[0,0],[1,1],'Units','pixels', 'HeadWidth',arrowHeadSz);
                    an.X = [ctrX(1) ctrX(seq(1)+1)];
                    if ~isequal(arrYsource(1), arrYtarget(seq(1)+1))
                        % Create arrows between columns
                        half_distY = (arrYsource(1) - arrYtarget(seq(1)+1))/2;
                        ln1 = annotation('line',[0,0],[1,1],'Units','pixels');
                        ln1.X = [an.X(1), an.X(1)];
                        ln1.Y = [arrYsource(1), arrYsource(1) - half_distY];
                        ln2 = annotation('line',[0,0],[1,1],'Units','pixels');
                        ln2.X = [an.X(1), an.X(2)];
                        ln2.Y = [arrYsource(1) - half_distY, arrYsource(1) - half_distY];
                        an.X = [an.X(2) an.X(2)];
                        an.Y = [arrYsource(1) - half_distY, arrYtarget(seq(1)+1)];
                    else
                        an.Y = [arrYsource(1) arrYtarget(seq(1)+1)];
                    end
                end
                for jj = 1:length(seq)-1
                    an = annotation('arrow',[0,0],[1,1],'Units','pixels', 'HeadWidth',arrowHeadSz);
                    an.X = [ctrX(seq(jj)+1) ctrX(seq(jj+1)+1)];
                    if ~isequal(arrYsource(seq(jj)+1), arrYtarget(seq(jj+1)+1))
                        % Create arrows between columns
                        half_distY = (arrYsource(seq(jj)+1) - arrYtarget(seq(jj+1)+1))/2;
                        ln1 = annotation('line',[0,0],[1,1],'Units','pixels');
                        ln1.X = [an.X(1), an.X(1)];
                        ln1.Y = [arrYsource(seq(jj)+1), arrYsource(seq(jj)+1) - half_distY];
                        ln2 = annotation('line',[0,0],[1,1],'Units','pixels');
                        ln2.X = [an.X(1), an.X(2)];
                        ln2.Y = [arrYsource(seq(jj)+1) - half_distY, arrYsource(seq(jj)+1) - half_distY];
                        an.X = [an.X(2) an.X(2)];
                        an.Y = [arrYsource(seq(jj)+1) - half_distY, arrYtarget(seq(jj+1)+1)];
                    else
                        an.Y = [arrYsource(seq(jj)+1) arrYtarget(seq(jj+1)+1)];
                    end
                end
            end
            drawnow;
            % Put panel on top:
            uistack(pan, 'top');
            % Reposition the figure so the top-left corner is in the same
            % position as the original figure:
            fH.Position(1) = origPos(1);
            fH.Position(2) = origPos(4) + origPos(2) - fH.Position(4);
            %%%%% Local functions -----------------------------------------
            
            %%%%% Figure Callbacks ----------------------------------------
            function callSetOpts(src,~,obj,panel)
                % This callback calls the method "setOpts" and changes the
                % color of the button when parameters were changed.
                
                % Get pipeline index:
                ppIndx = src.UserData;
                fcnInfo = obj.pipe(ppIndx);
                if isempty(fcnInfo.opts)
                    % Abort, if no parameters exist
                    return
                end
                % Block figure interaction by turning uipanel visible
                panel.Visible = 'on';
                % Call "setOpts"
                obj.setOpts(fcnInfo.name, fcnInfo.seq(1,1));
                % Update button's tooltip text:
                src.Tooltip = obj.genToolTipTxt(obj.pipe(ppIndx));
                % Change the button color to green:
                src.CData = myGreen;
                jiggleFig(ancestor(src, 'figure')) % Update tooltips
                % Show figure content
                panel.Visible = 'off';
            end
            
            function chooseSaveFileName(src,~,obj)
                % Create input dialog to select/change SaveFileName
                btn = gco;
                idxFcn = btn.UserData;
                step = obj.pipe(idxFcn);
                % Create input dialog box:
                if isempty(step.saveFileName)
                    defName = step.outFileName;
                else
                    defName = step.saveFileName;
                end
                [~,~,ext] = fileparts(defName);
                answer = inputdlg('Type file name:','Save step as',[1 50],{defName});
                if isempty(answer)
                    disp('Operation cancelled by User')
                    return
                end
                [~,newName,~] = fileparts(answer{1});
                if isempty(newName)
                    obj.pipe(idxFcn).b_save2File = false;
                    obj.pipe(idxFcn).saveFileName = '';
                else
                    % Update saveFileName in pipeline:
                    obj.pipe(idxFcn).b_save2File = true;
                    obj.pipe(idxFcn).saveFileName = [newName, ext];
                end
                % Update tooltips:
                figH = ancestor(src,'figure');
                btnList = findobj(figH,'Type','uicontrol');
                btn = btnList([btnList.UserData] == idxFcn);
                btn.Tooltip = obj.genToolTipTxt(obj.pipe(idxFcn));
                jiggleFig(ancestor(src, 'figure'))
            end
            
            function jiggleFig(figH)
                % For some reason, the tooltip doesn't update until we
                % resize the figure...
                figH.Position(3) = figH.Position(3) + 1;
                figH.Position(3) = figH.Position(3) - 1;
            end
            
        end
        %%%%%%-------------------------------------------------------------
    end
    methods (Access = {?DataViewer})
        %%%%%-- Methods for interfacing with DataViewer -------------------
        function loadDataFromDummyProtocol(obj,data,filename,dataInfo)
            % This methods imports the imaging data from DataViewer to this
            % class.
            if ~obj.b_inputFromDataViewer
                % Exclusive to "dummy" protocol instance.
                return
            end
            % Update current data:
            obj.current_data = data;
            % Store filename:
            obj.dv_inputFilename = filename;
            
            % Load DataHistory
            obj.dataHistory = obj.loadDataHistory(obj.ProtocolObj.SaveDir);
            if isempty(obj.dataHistory)
                return
            end
            % Update current info for data stored in RAM
            if exist('dataInfo','var')
                obj.current_info = dataInfo;
            else
                % Update input file name of first step in the pipeline:
                obj.pipe(1).inputFileName = obj.dv_inputFilename;
                % Update current data history:
                idxFile = strcmp({obj.dataHistory.filename},filename);
                if any(idxFile)
                    obj.current_info = obj.dataHistory(idxFile).info;
                else
                    obj.current_info = [];
                end
            end
        end
        
        function dataHistory = loadDataHistory(~,folder)
            % LOADDATAHISTORY loads the dataHistory structure from the "dataHistory.mat"
            % file in the protocol's saveFolder.
            % This function also cleans-up the dataHistory structure and
            % file if there are missing files in the folder.
            
            dataHistory = struct.empty(0,1);
            if ~isfile(fullfile(folder,'dataHistory.mat')) 
                return
            end
            load(fullfile(folder,'dataHistory.mat'));%#ok
            % Check if all files listed in dataHistory still exist in the
            % saveFolder:
            fileList = getFileList(folder,'all');            
            % Clean-up! Remove missing files from dataHistory:
            dataHistory(~ismember({dataHistory.filename}, fileList)) = [];
            if isempty(dataHistory)
                delete(fullfile(folder,'dataHistory.mat'));
            else
                % Overwrite dataHistory.mat file:
                save(fullfile(folder,'dataHistory.mat'),'dataHistory');
            end
        end
        
        function saveDataHistory(obj,SaveFolder, filename)
            % SAVEDATAHISTORY creates or overwrites the data history for a
            % .dat or .mat file (filename) inside the dataHistory.mat file.
            
            % Create a new dataHistory.mat file, if it doesn't exists in
            % the SaveFolder
            if isempty(obj.current_info)
                return
            end
            thisDataHistory = struct('filename',filename,'info',obj.current_info);
            if ~isfile(fullfile(SaveFolder,'dataHistory.mat'))
                dataHistory = thisDataHistory;%#ok
                save(fullfile(SaveFolder,'dataHistory.mat'),'dataHistory');
                return
            end
            % Reload dataHistory to ensure that the local one is
            % up-to-date (dataHistory clean-up):
            obj.dataHistory = obj.loadDataHistory(SaveFolder);
            
            % Append info to existing data history file:
            idxFile = strcmp({obj.dataHistory.filename},filename);
            if any(idxFile)
                % Overwrite the existing file with the new data history:
                obj.dataHistory(idxFile) = thisDataHistory;
            else
                % Append DataHistory file with new entry:
                obj.dataHistory(end+1) = thisDataHistory;
            end
            % Overwrite dataHistory.mat file:
            dataHistory = obj.dataHistory;%#ok
            save(fullfile(SaveFolder,'dataHistory.mat'),'dataHistory');
        end
        
    end
    
    methods (Access = private)
        
        function saveFolderLog(obj,SaveFolder)
            % SAVEFOLDERLOG saves the obj.folderLog table to the
            % "pipeLog.mat" file in the "SaveFolder".
            % This function checks for the height limit of the table and
            % up to "maxLogRows" and removes older entries before saving.
            
            if height(obj.folderLog) > obj.maxLogRows
                % Remove older entries if log exceeds maximum number of
                % rows (oler entries are on top).
                obj.folderLog(1:height(obj.folderLog) - obj.maxLogRows,:) = [];
            end
            % Save Folder log:
            pipeLog = obj.folderLog;
            save(fullfile(SaveFolder,'pipeLog.mat'),'pipeLog');
        end
        
        function loadFolderLog(obj,SaveFolder)
            % LOADFOLDERLOG loads the log table from the SaveFolder in the
            % obj.folderLog property.
            
            if ~isfile(fullfile(SaveFolder,'pipeLog.mat'))
                obj.folderLog = table();
                return
            end
            load(fullfile(SaveFolder,'pipeLog.mat'));%#ok
            obj.folderLog =  pipeLog;
        end
        
        function updateFolderLog(obj,task,runDateTime,b_completed,errMsg)
            % UPDATEFOLDERLOG creates/appends task entries in the
            % folderLog and PipelineSummary tables.
                        
            if isa(errMsg,'MException')
                MessageLong = getReport(errMsg,'extended','hyperlinks','off');
                MessageShort = getReport(errMsg,'basic','hyperlinks','off');
            elseif isempty(errMsg)
                MessageLong = 'No Errors';MessageShort = 'No Errors';
            else
                MessageLong = errMsg; MessageShort = errMsg;
            end
            thisLog = table({task.name},{task.inputFileName},{task.saveFileName}, b_completed, runDateTime,{MessageShort},{MessageLong},...
                'VariableNames', {'TaskName','InputFile','OutputFile', 'Completed', 'RunDateTime', 'Messages_short','Messages'});          
            % Update Pipeline summary table:
            obj.PipelineSummary = [obj.PipelineSummary; thisLog];
            % Update folder log:
            obj.folderLog = [obj.folderLog; thisLog];
        end
        
        function run_taskOnTarget(obj, task)
            % RUN_TASKONTARGET runs a task in the pipeline structure array in
            % TASK.
            % Input
            %   task(struct): current step from "obj.pipe".
            
            % Create function string and update log table:
            funcStr = createFcnString(obj, task);
            Messages = '';
            %  Execute the task:
            try
                % Control for missing input files:
                if (~isempty(task.inputFileName) && ~strcmpi(task.inputFileName,'data') ) && task.inputFrom ~= -1
                    errID = 'MATLAB:Umitoolbox:PipelineManager:FileNotFound';
                    errmsg = ['Input File for function ' task.name ' not found!'];
                    assert(isfile(fullfile(obj.current_saveFolder, task.inputFileName)),...
                        errID,errmsg);
                    obj.loadInputFile(task);
                end
                fprintf('\tFunction Name: %s \n\n',task.name);
                % Load options structure in the workspace.
                opts = task.opts; %#ok the "opts" structure is used in the EVAL function.
                % Evaluate function string:
                eval(funcStr);
                % Update log table and tell other methods that the function
                % was successfully run:
                obj.b_state = true;                
                % Update data history of current data with task:
                obj.updateStepInfo(task);
            catch ME
                obj.b_state = false;
                Messages = ME;
                disp('FAILED!');
            end
            % Save data to file:
            if task.b_save2File && obj.b_state
                % Look for tasks with output data from the current step and
                % save the data to a .DAT or .MAT file:
                obj.saveDataToFile(task,false)
            elseif obj.b_saveDataBeforeFail && ~obj.b_state
                obj.saveDataToFile(task,true);
            end
            % Update log table of current folder:
            obj.updateFolderLog(task,datetime('now'),obj.b_state,Messages);
            
        end
        
        function [newSeq,skippedSteps, selFile] = skipSteps(obj,thisSeq)
            % SKIPSTEPS looks in the folder "folderName" for .dat/.mat files
            % that can potentially replace the N first steps of the pipeline
            % sequence "thisSeq".
            % The criteria to replace the pipeline steps are:
            %    1- Same creation date time of the functions.
            %    2- Same function name
            %    3- Same function "opts" parameters
            %    4- No intermediate steps from the file is used as input to
            %    other sequences
            %    Inputs:
            %        thisSeq (struct): current sequence of the pipeline.
            %    Outputs:
            %        newSeq (struct): updated sequence without the redundant
            %            steps.
            %        skippedSteps (cell): list of skipped function names from
            %            "thisSeq".
            
            if obj.b_inputFromDataViewer && ~isempty(obj.current_info)
                % For DataViewer, also, look at the dataHistory of the
                % current data.
                % Temporarily add the dataHistory of the data in RAM to the
                % "dataHistory" property. This will be removed at the end
                % of this function.
                obj.dataHistory =  [obj.dataHistory; struct('filename','self','info',obj.current_info)];
            end
            
            newSeq = thisSeq;
            skippedSteps = {};
            newSeqArr = cell(size(obj.dataHistory));
            skippedStepsArr = newSeqArr;
            selFile = '';
            % Abort,if the dataHistory file doesn't exists.
            if isempty(obj.dataHistory)
                return
            end
            % Compare dataHistory with pipeline sequence:
            for ii = 1:length(obj.dataHistory)
                [newSeqArr{ii}, skippedStepsArr{ii}] = compareDataHistory(obj,thisSeq,obj.dataHistory(ii));
            end
            % Check for matches:
            if all(cellfun(@isempty,skippedStepsArr))
                return
            end
            
            % Select file with the largest number of steps:
            
            % Check if there is a file that contains all the sequence:
            idxSkipAll = cellfun(@isempty,newSeqArr);
            if any(idxSkipAll)
                indxSkip = find(idxSkipAll,1,'first');
                selFile = obj.dataHistory(indxSkip).filename;
                % Copy file if the "saveFileName of the current sequence has a different "saveFileName"
                for ii = length(thisSeq):-1:1
                    if thisSeq(ii).saveFileName
                        break
                    end
                end
                if ~strcmpi(selFile,'self')
                    thisSeq(ii).inputFileName = selFile;
                    obj.loadInputFile(thisSeq(ii));% Load step;
                end
                if ~isempty(thisSeq(ii).saveFileName) && ~strcmpi(thisSeq(ii).saveFileName,selFile) && ~obj.b_inputFromDataViewer
                    obj.saveDataToFile(thisSeq(ii),false)
                end
                newSeq = {}; skippedSteps = skippedStepsArr{indxSkip};
                return
            end
            nSteps = cellfun(@numel,skippedStepsArr);
            % Find file(s) with the maximum number of steps to be skipped:
            idxMax = find(nSteps == max(nSteps));
            if length(idxMax) > 1
                % When more than one file exist, try to find the one with
                % the same name as the inputFileName(s) in the current
                % sequence. If none is found, just pick the first one that has the same dataHistory.
                for ii = length(idxMax):-1:1
                    if any(arrayfun(@(x) strcmp(obj.dataHistory(idxMax(ii)).filename,x.inputFileName),thisSeq)) ||...
                            strcmpi(obj.dataHistory(idxMax(ii)).filename,'self')
                        break
                    end
                end
                idxMax = idxMax(ii);
            end
            selFile = obj.dataHistory(idxMax).filename;
            newSeq = newSeqArr{idxMax};
            skippedSteps = skippedStepsArr{idxMax};
            % For DataViewer, when the data is already in RAM:
            if strcmpi(selFile, 'self')
                newSeq(1).inputFileName = '';
                newSeq(1).inputFrom = -1;
            end
            % Remove "self" entry from dataHistory:
            idxSelf = strcmp({obj.dataHistory.filename},'self');
            obj.dataHistory(idxSelf) = [];
            %%%%%--Local function ------------------------------------------
            function [outSeq,skipNames] = compareDataHistory(obj,seqIn, dhIn)
                % COMPAREDATAHISTORY compares the dataHistory of "fileIn"
                % with the pipeline sequence "seqIN" and outputs the
                % updated sequence "outSeq" and the list of skipped steps
                % "skipNames".
                
                outSeq = seqIn;
                skipNames = {};
                
                % Compare datetimes:
                [~,locB] = ismember({dhIn.info.name},{obj.funcList.name});
                locB(locB == 0) = []; % Remove non-existent functions.
                if isempty(locB)
                    % If none of the functions exist, abort.
                    return
                end
                
                % Compare function names and optional parameters:
                thisHistory = struct();
                for jj = 1:length(seqIn)
                    thisHistory(jj).name = seqIn(jj).name;
                    thisHistory(jj).inputFileName = seqIn(jj).inputFileName;
                    thisHistory(jj).opts = seqIn(jj).opts;
                end
                % IF there is an input file, prepend the dataHistory to the current
                % sequence:
                idxFromDisk = find(~cellfun(@isempty,{seqIn.inputFileName}) & [seqIn.inputFrom] == 0,1,'first');
                idxFromDataViewer = find([seqIn.inputFrom] == -1,1,'first');
                if any(idxFromDisk > 0)
                    idxFrom = idxFromDisk;
                    % Load the input file dataHistory:
                    [~,inputFile,ext] = fileparts(seqIn(idxFromDisk).inputFileName);                    
                    idxFile = strcmp({obj.dataHistory.filename},[inputFile,ext]);
                    if any(idxFile)
                        ppInfo = obj.dataHistory(idxFile).info;
                    end
                elseif  any(idxFromDataViewer > 0) && ~isempty(obj.current_info)
                    idxFrom = idxFromDataViewer;                    
                    ppInfo = obj.current_info;
                end
                
                if exist('ppInfo','var')
                    %%%% For retrocompatibility
                    fNames = {'name','inputFileName','opts'};
                    for kk = 1:length(fNames)
                        if ~isfield(ppInfo,fNames{kk})
                            ppInfo(1).(fNames{kk}) = [];
                        end
                    end
                    %%%%%
                    prepend_info = struct();
                    for jj = 1:length(ppInfo)
                        prepend_info(jj).name = ppInfo(jj).name;
                        prepend_info(jj).inputFileName = ppInfo(jj).inputFileName;
                        prepend_info(jj).opts = ppInfo(jj).opts;
                    end
                    thisHistory = horzcat(prepend_info,thisHistory(idxFrom:end));
                end
                %%%%%% ----------------------------------------------------
                % Check for the existence of consecutive equal steps:
                b_isEqual = false(size(dhIn.info));
                for jj = 1:length(dhIn.info)
                    % Compare name and opts:
                    b_isEqual(jj) = (strcmpi(thisHistory(jj).name, dhIn.info(jj).name) && ...
                        strcmpi(thisHistory(jj).inputFileName, dhIn.info(jj).inputFileName) && ...
                        isequaln(thisHistory(jj).opts,dhIn.info(jj).opts));
                    if jj == length(thisHistory) || ~b_isEqual(jj)
                        break
                    end
                end
                if ~all(b_isEqual)
                    return
                end
                % Get indices of sequence corresponding to the dataHistory:
                [~,indxEqual] = ismember({dhIn.info(b_isEqual).name},{seqIn.name});indxEqual(indxEqual == 0) = [];
                % If any consecutive steps were equal, check if the intermediate steps
                % from dataHistory are inputs to other sequences:
                b_isInputFcn = arrayfun(@(x) numel(x) > 1, seqIn(indxEqual(1:end-1)));
                if any(b_isInputFcn)
                    return
                end
                % If all steps are to be skipped, the outSeq will return
                % empty
                outSeq = seqIn;
                outSeq(indxEqual) = []; % Erase first N steps corresponding to the dataHistory.
                if ~isempty(outSeq)
                    % If the sequence is partially skipped, update the
                    % sequence:
                    indx = find(~cellfun(@isempty,({outSeq.inputFrom})), 1,'first');
                    if any(indx)
                        outSeq(indx).inputFrom = 0; % Input from Disk.
                        outSeq(indx).inputFileName = dhIn.filename;
                        % Reset sequence indices:
                        for jj = 1:length(outSeq)
                            outSeq(jj).seqIndx = jj;
                        end
                    end
                end
                % Return steps to be skipped:
                skipNames = {seqIn(indxEqual).name};
            end
        end
        %%%%%%--Helpers for "addTask" method -----------------------------
        function task = setInput(obj,task)
            % SETINPUT selects the input to the function in "task". It
            % controls for multiple outputs and for functions with no
            % input. It updates the fields of "task" with the input
            % information.
            % !If the User cancels any of the dialogs, the function returns
            % an empty array!
            
            orig_seq = obj.current_seq;
            
            % For when this is the first item in the pipeline:
            if obj.current_seq == 0
                obj.current_seq = 1;
                obj.current_seqIndx = 1;
                task.seq = obj.current_seq;
                task.seqIndx = obj.current_seqIndx;
                if task.b_hasDataIn
                    task.inputFrom = 0; % Set input index to zero when the data comes from the Hard Drive.
                    task.inputFileName = obj.selectInputFileName(task.inputFrom, task.name);
                elseif (task.b_hasDataOut || task.b_hasFileOut ) && ~task.b_hasDataIn
                    task.inputFrom = 0;
                end
                if task.inputFileName == 0 % Means that the operation was cancelled by User within a GUI.
                    task = [];
                end
                return
            end
            % Procedure for adding new items to an existing pipeline:
            if ~isempty(task.inputFrom) || ( (task.b_hasDataOut || task.b_hasFileOut) && ~task.b_hasDataIn )
                % If the User forces the creation of a new branch or
                %   if the task function generates data from the disk.
                obj.current_seq = obj.current_seq + 1;
                obj.current_seqIndx = 1; % Reset current sequence index.
                if isempty(task.inputFrom)
                    % Force input to be "_LOCAL_":
                    task.inputFrom = 0;
                elseif task.inputFrom == 0
                    % prompt to file selection:
                    if isempty(task.inputFileName)
                        task.inputFileName = obj.selectInputFileName(0,task.name);
                    end
                else
                    obj.current_seqIndx = 2;
                    % Append new sequence to input function:
                    idx = strcmpi(obj.pipe(task.inputFrom).name,{obj.pipe.name}) & arrayfun(@(x) any(ismember(x.seq,[obj.pipe(task.inputFrom).seq])),obj.pipe)';
                    if ~any(idx)
                        % If there is no input function in the Sequence #1,
                        % abort:
                        warning('Input function "%s"  to function "%s" not found! Operation aborted.', obj.pipe(task.inputFrom).name, task.name)
                        obj.current_seq = orig_seq;
                        task = [];
                        return
                    end
                    obj.pipe(idx).seq = [obj.pipe(idx).seq obj.current_seq]; % Append sequence number
                    obj.pipe(idx).seqIndx = [obj.pipe(idx).seqIndx 1]; % Append sequence index as "1".
                    % Check input file:
                    if ~isempty(task.inputFileName)
                        %skip
                    elseif any(strcmpi('outFile',obj.pipe(idx).argsOut))
                        % Get File from input function:
                        task.inputFileName = obj.selectInputFileName(find(idx),task.name);%#ok
                        
                    else
                        % Ensure that the input function saves the
                        % "outData" to the disk.
                        if isempty(obj.pipe(idx).saveFileName)
                            % Create temporary file:
                            [~,defName,ext] = fileparts(obj.pipe(idx).outFileName);
                            obj.pipe(idx).saveFileName = [defName, obj.timeTag, ext];
                            obj.pipe(idx).b_save2File = true;
                        end
                        task.inputFileName = obj.pipe(idx).saveFileName;
                    end
                end
                task.seq = obj.current_seq;
                task.seqIndx = obj.current_seqIndx;
                if task.inputFileName == 0
                    obj.current_seq = orig_seq;
                    task = [];
                end
                return
            end
            
            % For the rest of the steps:
            if task.b_hasDataOut && ~task.b_hasDataIn
                % For functions that do not have "data" input but create an
                % output
                task.inputFrom = 0;
            elseif ~task.b_hasDataOut && ~task.b_hasFileOut && ~task.b_hasDataIn
                % Control for functions that do not have any data input or
                % output. *Generally, these are functions that changes
                % auxiliary files such as meta Data*.
                task.inputFrom = [];
            else
                % If this is not the first step, look backwards on the pipeline
                % to find the first function with "outData" or "outFile":
                for ii = length(obj.pipe):-1:1
                    inputFileName = obj.selectInputFileName(ii, task.name);
                    if ~isempty(inputFileName) | inputFileName == 0%#ok
                        task.inputFrom = ii;
                        task.inputFileName = inputFileName;
                        break
                    end
                end
                % If all sequence was scanned, and no function with an output
                % was found, force input from the Disk:
                if ii == 1 && isempty(inputFileName)
                    task.inputFrom = 0; % Input from Disk;
                    task.inputFileName = obj.selectInputFileName(0,task.name);
                end
            end
            % Update step index:
            obj.current_seqIndx = obj.current_seqIndx + 1;
            % Add sequence info to task:
            task.seq = obj.current_seq;
            task.seqIndx = obj.current_seqIndx;
            if task.inputFileName == 0
                obj.current_seq = orig_seq;
                task = [];
            end
        end
        
        function filename = selectInputFileName(obj,source,target)
            % SELECTINPUTFILENAME creates a dialog box for the User to
            % select a file that will be the input to the "task". If
            % the "task" field "inputFrom" is a function with multiple
            % files, it will show the possible outputs from the given
            % function. However, if the "inputFrom" is "_LOCAL_", the
            % prompt will show a list of existing files from the first
            % element in the saveFolderList.
            % Input:
            %   source (int): index of function from "pipe" or 0 for "_LOCAL_".
            filename = '';
                       
            % For functions as "source":
            if source ~=0
                % Get source function info:
                funcInfo = obj.pipe(source);
                if isempty(funcInfo)
                    warning('Input function not found. To create a new branch, the input function must be in the first sequence')
                    return
                end
                % Check if the selected input function has any outputs:
                if ~any(ismember(funcInfo.argsOut, {'outData','outFile'}))
                    % If not outputs are found, abort!
                    return
                end
                % If the output is "outData", return, otherwise, create
                % dialog box for file selection:
                if any(ismember(funcInfo.argsOut, 'outData'))
                    filename = 'data';
                    return
                end
                % For functions with "outFile":
                if numel(funcInfo.outFileName) == 1
                    % For a single file:
                    filename = funcInfo.outFileName{1};
                else
                    % Create dialog box so the user selects the file:
                    [indxFile, tf] = listdlg('ListString', funcInfo.outFileName,...
                        'SelectionMode','single','PromptString',{'Select a file from:',...
                        ['"' funcInfo.name '" as input to :' ], ['"' target '":']}, 'ListSize',[250,280],'Name', 'Select file');
                    if ~tf
                        disp('Operation cancelled by User')
                        filename = 0;
                        return
                    end
                    filename = funcInfo.outFileName{indxFile};
                end
            else
                % For "_LOCAL_" source, display the .dat files from the
                % first item of the SaveFolderList:                                
                fileList = getFileList(obj.SaveFolderList{1},'dat');
                               
                % If not files exist in the item's save folder, raise a
                % warning and abort:
                if isempty(fileList)
                    w = warndlg(['No valid data files found in folder: ' item.SaveFolder], 'Operation aborted!');
                    waitfor(w);
                    filename = 0;
                    return
                end
                % Create dialog box
                [indxFile, tf] = listdlg('ListString', fileList,'SelectionMode','single',...
                    'PromptString',{'Select a file as input to:', ['"' target '"']},'ListSize',[250,280],'Name','Select file');
                if ~tf
                    disp('Operation cancelled by User');
                    filename = 0;
                    return
                end
                filename = fileList{indxFile};
            end
        end
        
        function task = addDependency(obj,task)
            % ADDDEPENDENCY adds a step to the pipeline before the new "task"
            % with the function necessary to run the "task".
            
            % Check if the dependency function exists:
            if ~any(strcmpi(task.dependency,{obj.funcList.name}))
                % Abort, if the dependency fcn doesn't exist.
                error(['Failed to add ' task.name 'to the pipeline. The dependency function ' task.dependency ' does not exist!'])
            end
            % Check if it already exist in the pipeline:
            if any(strcmpi(task.dependency, {obj.pipe.name}))
                % Cancel operation if the dependency function already exists
                % in the pipeline.
                return
            end
            % Add dependency function to the current sequence
            if obj.current_seq > 0
                obj.addTask(task.dependency,'fromSeq',obj.current_seq);
            else
                obj.addTask(task.dependency);
            end
            
        end
        %%%%%%-------------------------------------------------------------
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
            opts_values = struct();
            
            disp('Creating Fcn list...');
            list = dir(fullfile(obj.fcnDir, '*','*.m'));
            for i = 1:length(list)
                out = parseFuncFile(list(i));
                % Validate if all input arguments from the function are
                % "valid" inputs keywords:
                kwrds_args = {'data', 'SaveFolder', 'RawFolder', 'opts'};
                kwrds_out = {'outFile', 'outData'};
                if all(ismember(out.argsIn, kwrds_args)) && all(ismember(out.argsOut, kwrds_out))
                    [~,list(i).name, ~] = fileparts(list(i).name);
                    list(i).info = out;
                    list(i).info.opts_def = list(i).info.opts; % Duplicate default params.
                    % add boolean to indicate if the function as inputs and
                    % outputs (This will be used by other methods)
                    list(i).info.b_hasDataIn = any(contains(list(i).info.argsIn,'data','IgnoreCase',true));
                    list(i).info.b_hasDataOut = any(contains(list(i).info.argsOut,'outData','IgnoreCase',true));
                    list(i).info.b_hasFileOut = any(contains(list(i).info.argsOut,'outFile','IgnoreCase',true));
                    obj.funcList = [obj.funcList ; list(i)];
                end
            end
            disp('Function list created!');
            function info = parseFuncFile(fcnStruct)
                info = struct('argsIn', {},'argsOut', {}, 'outFileName', '', 'opts', [],...
                    'opts_def',[],'opts_vals',[], 'dependency','');
                % Read the '.m' file content and exclude comments (lines
                % starting with the "%" character):
                fid = fopen(fullfile(fcnStruct.folder, fcnStruct.name),'r');
                txt = '';
                while 1
                    tline = fgetl(fid);
                    if tline == -1
                        break
                    end
                    if ~startsWith(strip(tline),'%')
                        txt=[txt,sprintf('%s\n',tline)];%#ok
                    end
                end
                fclose(fid);
                clear fid tline
                % Parse function header to get input and output variables:
                funcStr = erase(regexp(txt, '(?<=function\s*).*?(?=\r*\n)', 'match', 'once'),' ');
                outStr = regexpi(funcStr,'.*(?=\=)', 'match', 'once');
                out_args = regexpi(outStr, '\[*(\w*)\,*(\w*)\]*', 'tokens', 'once');
                info(1).argsOut = out_args(~cellfun(@isempty, out_args));
                [~,funcName,~] = fileparts(fcnStruct.name);
                expInput = ['(?<=' funcName '\s*\().*?(?=\))'];
                str = regexpi(funcStr, expInput, 'match', 'once');
                str = strip(split(str, ','));
                info.argsIn = str(~strcmp(str, 'varargin'));
                % Get Default outputs:
                expOutput = 'default_Output\s*=.*?(?=\n)';
                str = regexpi(txt, expOutput, 'match', 'once');
                if isempty(str)
                    default_Output = '';
                else
                    eval(str)
                end
                info.outFileName = default_Output;
                % Get dependent function:
                dependency = '';
                depStr = regexpi(txt, 'dependency\s*=.*?(?=\n)', 'match', 'once');
                if ~isempty(depStr)
                    eval(depStr)
                end
                info.dependency = dependency;
                % Parse default optional params struct:
                expOpts = 'default_opts\s*=.*?(?=\n)';
                str = regexpi(txt, expOpts, 'match', 'once');
                if ~isempty(str)
                    eval(str)
                    info.opts = default_opts;
                    info.argsIn{end+1} = 'opts';
                    % Parse optional params values struct:
                    optsVals = 'opts_values\s*=.*?(?=\n)';
                    str_opts = regexpi(txt, optsVals, 'match', 'once');
                    if ~isempty(str_opts)
                        eval(str_opts)
                        info.opts_vals = opts_values;
                    end
                end
                % SPECIAL CASE. Look for "object"s as optional arguments in
                % function first lines:
                expObj = 'default_object\s*=';
                str = regexpi(txt, expObj, 'match', 'once');
                if ~isempty(str)
                    info.argsIn{end+1} = 'object';
                end
            end
        end
        
        function fcnStr = createFcnString(obj,task)
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
            
            argsIn = replace(task.argsIn, ["RawFolder", "SaveFolder", "data"],...
                {['''' obj.current_rawFolder ''''],['''' obj.current_saveFolder ''''],...
                'obj.current_data'});
            argsOut = replace(task.argsOut, ["outData", "outFile"],...
                {'obj.current_data', 'obj.current_outFile'});
            if isempty(argsOut)
                fcnStr = [fcnStr ';' task.name '(' strjoin(argsIn,',') ');'];
            else
                fcnStr = [fcnStr ';' '[' strjoin(argsOut, ',') ']=' task.name '(' strjoin(argsIn,',') ');'];
            end
            fcnStr = strip(fcnStr,'left', ';');
            
        end
        
        function updateStepInfo(obj,step)
            % This function creates or  updates the "current_info" structure
            %
            % The dataHistory contains all information about the functions'
            % parameters used to create the current "data" and when it was run.
            %
            % Input:
            %    step(struct) : current step of the pipeline;
            
            funcInfo = obj.funcList(strcmp(step.name, {obj.funcList.name}));
            % Create a local structure with the function's info:
            thisStep = obj.getStepInfo(funcInfo, step.opts, step.inputFileName, '');
            % Append to current dataHistory:
            obj.current_info = [obj.current_info;thisStep];
            
            % If the data was already saved as a .dat or .mat file by the
            % function, update the function's info directly in the
            % dataHistory.mat file:
            if any(strcmp(step.argsOut, 'outFile'))
                obj.current_info(end).outFileName = obj.current_outFile; % Update outFileName list with the actual files generated by the function in "step".
                % In case the step ouput is a .DAT or .MAT file, update the
                % dataHistory for each one:
                for jj = 1:length(obj.current_outFile)
                    obj.saveDataHistory(obj.current_saveFolder,obj.current_outFile{jj});
                end
            end
        end
        
        function out = getStepInfo(~,fcnInfo, optsStruct, inputFileName, outFileName)
            % This function creates a structure containing information about an
            % analysis function.
            
            % Inputs:
            %   fcnInfo (struct): structure containing the function's basic informations with
            %       fields:
            %           -name (char): name of the analysis function.
            %           -folder (char): path where the analysis function file is located.
            %           -creationDatetime(datetime): timestamp of the creation of the
            %               analysis function file.
            %   optsStruct (struct): structure containing optional parameters of the
            %       analysis function.
            %   inputFileName(cell|char): name of the input file(s) to the function.
            %   outFileName(cell|char): name of the output file(s) from the function.
            %   This field is used just by functions that create files already. Just to
            %   keep track of the files that were created.
            % Output:
            %   out (struct): structure with the information necessary for the
            %       dataHistory variable in the data's metaData.
            
            out = struct('runDatetime', datetime('now'), 'name', {fcnInfo.name},...
                'folder', {fcnInfo.folder}, 'creationDatetime',...
                datetime(fcnInfo.datenum, 'ConvertFrom', 'datenum'),...
                'opts', optsStruct, 'inputFileName',inputFileName, 'outFileName',outFileName);
        end
              
        function loadInputFile(obj,step)
            % This function loads the data and metaData (if applicable)
            % from a .DAT or .MAT file indicated by the inputFileName field
            % in the "step" structure.
            % Input:
            %   step (struct): structure containing the info of the current
            %   task that will be executed inside method
            %   "run_taskOnTarget".
            
            disp(['Loading input file : ' step.inputFileName])
            
            if endsWith(step.inputFileName, '.dat')
                %If the InputFile is a .DAT file:
                obj.current_data = loadDat(fullfile(obj.current_saveFolder, step.inputFileName));
            else
                % If the InputFile is a .MAT file
                obj.current_data = load(fullfile(obj.current_saveFolder,...
                    step.inputFileName));
            end
            
            obj.current_info = [];
            
            % Update Data History file handle:
            obj.dataHistory = obj.loadDataHistory(obj.current_saveFolder);
            if ~isempty(obj.dataHistory)
                % Update current info:
                idxFile = strcmp({obj.dataHistory.filename},step.inputFileName);
                if any(idxFile)
                    obj.current_info = obj.dataHistory(idxFile).info;
                end
            end
        end
        
        function saveDataToFile(obj,step,b_failed)
            % This methods looks back in the pipeline from "step" for tasks
            % with "data" or "stats data" as output and saves the current data to a
            % .DAT or .MAT file.
            % Input:
            %    step(struct) : info of the current task in the pipeline.
            %    b_failed (bool): If TRUE, this means that the "step"
            %       execution failed. In this case, we look for the previous
            %       step that generated the "current_data" and save it. If
            %       FALSE, we just save the data to a file with name
            %       "step.saveFileName".
            
            % Get the pipeline until the task in "step":
            indx = obj.findStep(step);
            %%%% FOR TESTING
            assert(~isempty(indx),'FAILED TO FIND STEP IN PIPELINE');
            assert(numel(indx) == 1,'DUPLICATE STEPS FOUND IN PIPELINE! DEBUG!');
            %%%%%
            % Look back in pipeline for steps with "data" or "stats data"
            % as output and save the current data using the task's info:
            if b_failed
                % Ignore the current step index:
                indx = indx-1;
            end
            if indx == 0
                warning('First step of the pipeline failed! No previous data exists to be saved. Operation aborted!')
                return
            end
            for ii = indx:-1:1
                if obj.pipe(ii).b_hasDataOut
                    % Skip step without data output
                    break
                end
            end
            % Change save file name if the user decides not to overwrite
            % the data:
            saveFileName = obj.pipe(ii).saveFileName;
            if ~obj.b_overwriteFiles
                [~,name, ext] = fileparts(saveFileName);
                fList = getFileList(obj.current_saveFolder,ext);
                cnt = 2;
                newName = name;
                while 1
                    if any(strcmpi([newName, ext],fList))
                        newName = [name '_' num2str(cnt)];
                    else
                        break
                    end
                    cnt = cnt + 1;
                end
                % Save new filename to pipeline:
                saveFileName = [newName ext];
            end
            
            if b_failed
                if isempty(obj.pipe(ii).saveFileName) && ischar(obj.pipe(ii).outFileName)
                    % If the previous step was a function with data output,
                    % use the default file name to save the data:
                    [~,name,ext] = fileparts(obj.pipe(ii).outFileName);
                    saveFileName = [name '_recovered' ext]; % string to append to files saved before an error.
                elseif contains(obj.pipe(ii).saveFileName, obj.timeTag,'IgnoreCase',true)
                    % If the previous step already saved a temporary file, just
                    % rename it.
                    saveFileName = replace(obj.pipe(ii).saveFileName,obj.timeTag,'_recovered');
                    movefile(fullfile(obj.current_saveFolder,obj.pipe(ii).saveFileName),...
                        fullfile(obj.current_saveFolder,saveFileName));
                    % Save data history:
                    obj.saveDataHistory(obj.current_saveFolder,saveFileName);
                    return
                end
            end
            
            % Save current data to file:
            if endsWith(saveFileName,'.dat')
                % For .dat files:
                saveDat(fullfile(obj.current_saveFolder,saveFileName),...
                    obj.current_data);
            else
                % For .mat files:
                disp('Writing data to .MAT file ...')
                S = obj.current_data;
                save(fullfile(obj.current_saveFolder,saveFileName),...
                    '-struct', 'S', '-v7.3');
                disp(['Data saved in : "' fullfile(obj.current_saveFolder,saveFileName) '"'])
            end
            % Save data history:
            obj.saveDataHistory(obj.current_saveFolder,saveFileName);
        end
        
        function deleteTemporaryFiles(obj,folder)
            % DELETETEMPORARYFILES removes .dat and .mat files from
            % "folder" that were automatically generated during the
            % pipeline due to the existence of branches. The files to be
            % deleted are appended with the "timeTag" of the current
            % pipeline.
            
            % Get list of files in folder:
            list = getFileList(folder,'all');
            % delete files with the "timeTag":
            idxTemp = contains(list, obj.timeTag,'IgnoreCase',true);
            if any(idxTemp)
                fprintf('Deleting pipeline temporary files in folder "%s"...\n',folder);
                cellfun(@(x) delete(fullfile(folder,x)),list(idxTemp));
                disp('Done.')
            end
        end
        
        function genErrorReport(obj)
            % GENERRORREPORT saves the error messages from the Pipeline
            % Summary table to a .txt file inside the folder
            % "PipelineErrorLogs" folder in protocol's save folder.
            
            % Check if any error occurred in the pipeline
            if all(obj.PipelineSummary.Completed)
                return
            end
            
            % Ask User if he/she wants to save the error log:
            answer = questdlg('One or more steps of the pipeline failed. Generate error log file?',...
                'Errors found!', 'Yes','No','Yes');
            waitfor(answer);
            if strcmpi(answer, 'no')
                return
            end
            % Repackage error messages into string:
            str = sprintf(['-------------------- Pipeline error log --------------------\n'...
                'Pipeline executed at:%s\nReport generated at: %s\nTotal number of failed tasks: %d\n%s\n'],...
                datestr(datetime(obj.timeTag,'InputFormat','_ddMMyyyyHHmmss')),...
                datestr(datetime('now')),sum(obj.PipelineSummary.Completed), repmat('-',1,60));
            % Get error messages in execution order:
            errorTab = obj.PipelineSummary(~obj.PipelineSummary.Completed,:);
            headers = errorTab.Properties.VariableNames;
            for ii = 1:height(errorTab)
                info = table2cell(errorTab(ii,:));
                id = [headers([1:5,9]);[info([1:5]) {datestr(info{9})}]];
                % Trim error messsage to remove references to
                % PipelineManager methods:
                idx = strfind(info{10},'Error in PipelineManager');
                msg = info{10}(1:idx(1)-1);
                str = [str, sprintf('Recording Info:\n\t%s : %s\n\t%s : %s\n\t%s : %s\n\t%s : %s\n\t%s : %s\n\t%s : %s\n',id{:})];%#ok
                str = [str, sprintf('Error Message:\n""\n%s\n""\n%s\n', msg,repmat('-',1,60))];%#ok
            end
            % Save to .txt file:
            if obj.b_inputFromDataViewer
                % For DataViewer:
                folder = obj.ProtocolSaveFolder;
            else
                folder = fullfile(obj.ProtocolSaveFolder,'PipelineErrorLogs');
            end
            if ~exist(folder,'dir')
                mkdir(folder);
            end
            filename = ['error_log_' datestr(datetime(obj.timeTag,'InputFormat','_ddMMyyyyHHmmss'),'dd_mm_yyyy_HH_MM'), '.txt'];
            fid = fopen(fullfile(folder,filename), 'w');
            fprintf(fid,'%s',str);
            fclose(fid);
            % Open system's file explorer:
            openFolder(folder);
        end
        %%%%%%-- WAITBAR methods-------------------------------------------
        function setWaitBar(obj,tag,varargin)
            % This method creates two "waitbar" dialogs.
            % The first shows the progress of the pipeline runs across objects
            % while the second one shows the progress of tasks in a given
            % object.
            % Inputs:
            %   tag (char): "Initialize" : creates the 2 waitbars.
            %               "UpdateItem" : updates bar1.
            %               "UpdateTask" : updates bar2.
            
            p = inputParser;
            addRequired(p,'obj')
            addRequired(p,'tag',@ischar)
            addOptional(p, 'currBarVal',0, @(x) isnumeric(x) & x >= 0)
            addOptional(p, 'totalBarVal',1, @(x) isnumeric(x) & x >= 0)
            addParameter(p,'taskName','',@ischar);
            parse(p,obj,tag,varargin{:})
            barVal = p.Results.currBarVal;
            totalBarVal = p.Results.totalBarVal;
            taskName = p.Results.taskName;
            clear p
            
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
                    waitbar(barVal/totalBarVal, obj.h_wbItem, ['Item ' num2str(barVal)...
                        '/' num2str(totalBarVal)]);
                    obj.h_wbTask.Name = obj.current_saveFolder;
                case 'UpdateTask'
                    waitbar(barVal, obj.h_wbTask, ['Running "' taskName '"']);
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
        %%%%%%---Pipeline Viz helper methods ------------------------------
        function txt = textifyOpts(~,opts)
            % TEXTIFYOPTS creates a formatted text version of the "opts"
            % structure.
            
            logicMap = containers.Map([true,false],{'Yes','No'});
            fn = fieldnames(opts);
            vals = cell(size(fn));
            for ii = 1:length(fn)
                val = opts.(fn{ii});
                if isnumeric(val)
                    vals{ii} = num2str(val);
                elseif ischar(val)
                    vals{ii} = val;
                elseif iscell(val)
                    vals{ii} = strjoin(val,', ');
                elseif islogical(val)
                    vals{ii} = logicMap(val);
                else
                    error('Unknown value data type')
                end
            end
            if size(fn,1)> size(fn,2)
                fn = fn'; vals = vals';
            end
            info = vertcat(fn,vals);
            txt = [sprintf('Parameters:\n'), sprintf('%s: "%s"\n',info{:})];
        end
        
        function tipTxt = genToolTipTxt(obj,step)
            % Add tooltips for each one:
            if ~isempty(step.opts)
                tipTxt = obj.textifyOpts(step.opts);
                if ~isempty(step.saveFileName)
                    tipTxt = [tipTxt, sprintf('\n%s',repmat('-',1,20)),  sprintf('\nSave to file: "%s"',step.saveFileName)];
                end
            else
                tipTxt = 'No Parameters';
            end
            if ~isempty(step.inputFrom)
                sourceFile = step.inputFileName;
                if step.inputFrom == 0
                    source = 'Disk';
                elseif step.inputFrom == -1
                    % SPECIAL CASE for DataViewer
                    source = 'DataViewer';
                else
                    source = obj.pipe(step.inputFrom).name;
                end
                % Add input source to data tip:
                tipTxt = [tipTxt, sprintf('\n%s',repmat('-',1,20)),sprintf('\nInput From: "%s"',source)];
                if ~isempty(sourceFile)
                    if contains(sourceFile, obj.timeTag, 'IgnoreCase',true)
                        % Replace temporary file names with "data"
                        sourceFile = 'data';
                    end
                    tipTxt = [tipTxt, sprintf('\nInput File: "%s"',sourceFile)];
                end
            end
        end
        
        function rgb = prettyfyBtn(~, color, btnSz, bcgCol)
            % This creates a CData with the given color (rgb triplet)
            % that mimics a button with rounded corners.
            w = round(btnSz(1)); % width
            h = round(btnSz(2)); % height
            % Choose the radius of the rounded corners
            r = round(.05*w); % radius
            % Define the x and y grids using meshgrid
            [x, y] = meshgrid(1:w, 1:h);
            % Define a mask for the rounded rectangle
            mask = ((x <= r) & (y <= r) & (sqrt((x - r).^2 + (y - r).^2)) <= r) | ... % top left corner
                ((x >= w - r) & (y <= r) & (sqrt((x - w + r).^2 + (y - r).^2)) <= r) | ... % top right corner
                ((x <= r) & (y >= h - r) & (sqrt((x - r).^2 + (y - h + r).^2)) <= r) | ... % bottom left corner
                ((x >= w - r) & (y >= h - r) & (sqrt((x - w + r).^2 + (y - h + r).^2)) <= r); % bottom right corner
            % Create the rounded rectangle
            rect = ones(size(x));
            rect(mask) = 0;
            msk1 = zeros(size(x));
            msk1(r:end-r,:) = 1;
            msk1(:,r:end-r) = 1;
            antiMsk = ~(~rect | msk1);
            rim = bwmorph(~antiMsk,'remove');
            r = repmat(color(1),h,w); r(antiMsk) = bcgCol; r(rim) = 0;
            g = repmat(color(2),h,w);g(antiMsk) = bcgCol; g(rim) = 0;
            b = repmat(color(3),h,w); b(antiMsk) = bcgCol; b(rim) = 0;
            rgb = cat(3,r,g,b);
        end
        %%%%%%-------------------------------------------------------------
        function indx = findStep(obj,step)
            % This function gives the index of the task "step" in the pipeline "pipe"
            % It compares all fields from the structure to identify the
            % step.
            
            fn = fieldnames(step);
            % Exclude fields that may change during pipeline execution:
            fn = setdiff(fn, {'seqIndx', 'inputFrom', 'b_save2File', 'inputFileName'});
            idx = false(length(obj.pipe), length(fn));
            for ii = 1:length(fn)
                idx(:,ii) = arrayfun(@(x) isequaln(step.(fn{ii}),x.(fn{ii})),obj.pipe);
            end
            indx = find(all(idx,2));
        end
        %%%%%%---Pipeline Viz helper methods ------------------------------
    end
end

%%%%%--Local functions-----------------------------------------------------

% These functions work with the "setOpts" method to create an input dialog
% box to set the optional parameters for a given function.

function out = buildInputDlg(fcnName,fields,currVals,defVals,listVals,typeVals)
% This function creates an input dialog for user data entry.

% Create output variable with current values:
out = {};
currOpts = {};
for i = 1:length(fields)
    currOpts(i,1) = fields(i);
    currOpts(i,2) = currVals(i);
end
dlg = uifigure('Name',['Set Parameters for: ' fcnName], 'NumberTitle','off','Position',[0,0,310,240],...
    'MenuBar','none', 'ToolBar', 'none','Visible','off', 'Resize', 'on', 'CloseRequestFcn', @figCloseRequest);
movegui(dlg, 'center');
myFontSize = 12;

% Create grid layout:
g = uigridlayout(dlg);
g.ColumnWidth = {max(cellfun(@(x) length(x), fields))*myFontSize,'1x'};
g.ColumnWidth = {'1x','1x'};
g.ColumnSpacing = 5;
% Set RowHeight of gridLayout depending on the type of variable:
idx = strcmpi('charArrayMultiSelect', typeVals);
rh = {};
for i = 1:length(idx)
    if idx(i)
        rh{i} = '1x';
    else
        rh{i} = 30;
    end
end
g.RowHeight = [rh, {60}];
% Update figure height:
dlg.InnerPosition(4) = sum([g.RowHeight{:}, 2*g.RowSpacing]);
for i = 1:length(fields)
    lb = uilabel(g,'Text', [fields{i} ': ']);
    lb.FontSize = myFontSize;
    lb.HorizontalAlignment = 'center';
    lb.Layout.Row = i;
    lb.Layout.Column = 1;
    % Create Value entry objects:
    switch typeVals{i}
        case{'numericArray','logical', 'charArray'}
            vo = uidropdown(g);
        case 'charArrayMultiSelect'
            vo = uipanel(g, 'Scrollable', 'on');
        case {'numericRange'}
            vo = uieditfield(g, 'numeric', 'ValueChangedFcn', @lockTextField);
        otherwise
            vo = uieditfield(g, 'ValueChangedFcn', @lockTextField);
    end
    % Set position of element in uigrid:
    vo.Layout.Row = i;
    vo.Layout.Column = 2;
    
    % Populate elements with current values:
    switch typeVals{i}
        case 'numericArray'
            vo.Items = arrayfun(@(x) num2str(x), listVals{i}, 'UniformOutput',false);
            vo.Value = vo.Items(listVals{i} == currVals{i});
        case 'charArray'
            vo.Items = listVals{i};
            vo.Value = vo.Items(strcmp(listVals{i},currVals{i}));
        case 'charArrayMultiSelect'
            % Create uipanel with series of checkboxes:
            idxDef = ismember(listVals{i},currVals{i});
            glChar = uigridlayout(vo,[length(listVals{i}),1]);
            glChar.RowHeight = repmat({20},size(listVals{i}));
            glChar.Scrollable = 'on';
            for jj = 1:length(listVals{i})
                c = uicheckbox('Parent',glChar, 'Text', listVals{i}{jj},...
                    'Value', idxDef(jj), 'ValueChangedFcn', @lockCheckBox);
                c.Layout.Row = jj;
            end
            % Resize figure to accomodate checkbox list:
            dlg.Position(4) = dlg.Position(4) + 20*length(listVals{i});
            movegui(dlg, 'center');
        case 'logical'
            keys = {'Yes','No'};
            [~,locB] = ismember([true, false], listVals{i});
            vo.Items = keys(locB);
            vo.Value = vo.Items(listVals{i} == currVals{i});
        case 'numericRange'
            vo.Value = currVals{i};
            % Set min and max range:
            vo.Limits = [listVals{i}];
        otherwise
            vo.Value = num2str(currVals{i});
            idx_char = cellfun(@ischar, listVals{i});
            name = strjoin(cellfun(@(x) ['"', x, '"'],listVals{i}(idx_char), 'UniformOutput', false),', ');
            if ~all(idx_char)
                % Check if numeric values are scalar or an array.
                if numel(listVals{i}{~idx_char}) > 1
                    vo.Tooltip = ['Type ' name ' or a range of numbers separated by a semicolon (ex. 0;2;3;4).'];
                else
                    vo.Tooltip = ['Type ' name ' or a single number.'];
                end
            else
                vo.Tooltip = 'Type a string';
            end
    end
end
% Add "Save button"
g2 = uigridlayout(g);
g2.Padding = [15 0 15 0];
btnSave = uibutton(g2, 'Text', 'Ok', 'BackgroundColor', [0 .7 0], 'ButtonPushedFcn',@saveParams);
btnSave.Layout.Row = 1;
btnSave.Layout.Column = 1;
btnReset = uibutton(g2, 'Text', 'Reset', 'BackgroundColor', [.7 .7 .7],...
    'ButtonPushedFcn',@change2Defs);
btnReset.Layout.Row = 1;
btnReset.Layout.Column = 2;
btnReset.Tooltip = 'Reset current values to function''s default';
dlg.Visible = 'on';
pause(0.2);
uiwait(dlg);
% Return current values if user closes the figure:
% if isempty(out)
%     out = currOpts;
%     return
% end
% Clean "out" before completion:
for i = 1:size(out,1)
    out{i,1} = currOpts{i,1};
    switch typeVals{i}
        case 'numericArray'
            out{i,2} = str2double(out{i,2});
        case {'charArray', 'charArrayMultiSelect'}
            %Do Nothing%
        case 'numericRange'
            %Do Nothing%
        case 'logical'
            if strcmp(out{i,2},'Yes')
                out{i,2} = true;
            else
                out{i,2} = false;
            end
        case 'mixArray'
            % Transform string digits into double:
            tmp = str2num(erase(out{i,2},' '));
            if ~isempty(tmp)
                out{i,2} = tmp;
            end
        otherwise
            % This should not be reached.
            error('Unknown data type')
    end
end
% Callbacks for pushbutton and uifigure:
    function change2Defs(src,~)
        % This function changes all fields back to the default values.
        disp('Changing to default...');
        gl = src.Parent.Parent;
        nRows = length(gl.RowHeight)-1;
        layout_info = arrayfun(@(x) get(x, 'Layout'),gl.Children);
        for ii = 1:nRows
            valueObj = gl.Children([layout_info.Row] == ii & [layout_info.Column] == 2);
            if isa(valueObj, 'matlab.ui.control.DropDown')
                switch typeVals{ii}
                    case {'numericArray', 'logical'}
                        idx_def = listVals{ii} == defVals{ii};
                    case 'charArray'
                        idx_def = strcmp(listVals{ii}, defVals{ii});
                end
                valueObj.Value = valueObj.Items(idx_def);
            elseif isa(valueObj,'matlab.ui.container.Panel')
                % For charArrayMultiSelect case:
                idxDefReset = strcmpi(defVals{ii}, listVals{ii});
                arrayfun(@(x,y) set(x,'Value', y), valueObj.Children.Children, idxDefReset);
            else
                valueObj.Value = defVals{ii};
            end
        end
    end
    function saveParams(src,~)
        % This callback saves the parameters selected by the user in the
        % GUI to the "out" variable and closes the figure;
        gl = src.Parent.Parent;
        nRows = length(gl.RowHeight)-1;
        layout_info = arrayfun(@(x) get(x, 'Layout'),gl.Children);
        for ii = 1:nRows
            valueObj = gl.Children([layout_info.Row] == ii & [layout_info.Column] == 2);
            if isa(valueObj,'matlab.ui.container.Panel')
                % Special case: get checkbox values inside uipanel:
                out{ii,2} = listVals{ii}(arrayfun(@(x) x.Value, valueObj.Children.Children));
            else
                out{ii,2} = valueObj.Value;
            end
        end
        delete(src.Parent.Parent.Parent); % Close figure.
    end
    function figCloseRequest(src,~)
        % Displays a message to user and erase output:
        out = {};
        src.Visible = 'off';
        w = warndlg('Operation cancelled by User! Changes won''t be applied!', 'Set options Cancelled', 'modal');
        waitfor(w);
        delete(src);
    end
    function lockCheckBox(src,~)
        % This callback avoids the unchecking of the last checked box.
        idxState = arrayfun(@(x) x.Value, src.Parent.Children);
        if ~any(idxState)
            src.Value = 1;
        end
    end
    function lockTextField(src,evnt)
        % This callback avoids leaving any input field EMPTY.
        if ( isempty(evnt.Value) ) || ( ~isnumeric(evnt.Value)  && isempty(strip(evnt.Value)) )
            beep
            src.Value = evnt.PreviousValue;
        end
    end
end

