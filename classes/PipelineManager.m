classdef PipelineManager < handle
    % PIPELINEMANAGER manages data processing pipelines.
    % This class allows the creation of analysis pipeline and manages the execution of
    % functions. In addition, it controls for failed / completed  steps in the pipeline.
    
    properties
        b_ignoreLoggedFiles  = false %(bool) If true, PIPELINEMANAGER will ignore identical jobs previously run.
        b_saveDataBeforeFail = false %(bool) If true, the more recent data ("current_data") in the pipeline will be...
        % saved to a file when an error occurs.
        PipelineSummary % Shows the jobs run in the current Pipeline
    end
    properties (SetAccess = ?PipelineConfig_dlgBox)
        % Structure array containing steps of the pipeline:
        pipe = struct('argsIn', {},'argsOut',{},'outFileName','',...
            'inputFileName', '', 'b_save2File', logical.empty, 'saveFileName',...
            '', 'opts',struct.empty,'opts_vals',struct.empty,...
            'opts_def',struct.empty ,'name','', 'inputFrom','', 'seq',[],'seqIndx',[]);% !!If the fields are changed, please apply
        % the same changes to the property's set method. 
    end
    properties (SetAccess = private)       
        funcList struct % structure containing the info about each function in the "fcnDir".
        ProtocolObj Protocol % Protocol Object.
        fcnDir char % Directory of the analysis functions.
        ClassName char % Name of the class that the pipeline analysis functions will run.
        ClassLevel int16 % Level of the class in protocol's hierarchy (1 = Modality, 2 = Acquisition, 3= Subject);
    end
    properties (Access = private)               
        b_taskState = false % Boolean indicating if the task in pipeline was successful (TRUE) or not (FALSE).
        tmp_LogBook % Temporarily stores the table from PROTOCOL.LOGBOOKFILE
        tmp_BranchPipeline % Temporarily stores LogBook from a Hierarchical branch.        
        tmp_TargetObj % % Temporarily stores an object (TARGEROBJ).
        current_task % Task structure currently running.
        current_pipe % Pipeline currently running.
        current_data % Data available in the workspace during pipeline.
        current_metaData % MetaData associated with "current_data".
        current_outFile cell % List of file names created as output from some of the analysis functions.
        b_state logical % True if a task of a pipeline was successfully executed.        
        % It can be the name of an existing file, or
        % "outFile" for a function that creates a file
        % such as "run_ImagesClassification".
        h_wbItem % Handle of waitbar dialog showing the progress of the pipeline across protocol's objects.
        h_wbTask % Handle of waitbar dialog showing the progress of the pipeline in the current object.
        targetObjFullID char % Full ID of the current object.
        % For example, if the current object is a
        % modality the full ID is:
        % SubjID--AcqID--ModID.
        timeTag % timestamp used as "tag" for temporary files created during the pipeline execution.
        current_seq = 0 % index of current sequence in pipeline .
        current_seqIndx = 0 % index of step in current sequence in pipeline.
        b_newSeq = false; % boolean to indicate if a new sequence was created in the previous step.
        b_pipeIsValid = false % TRUE, if the pipeline passed all validations and is ready to execution.
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
            % Create timestamp tag:
            obj.timeTag = datestr(datetime('now'),'_ddmmyyyyHHMMSS');
            % Set level of selected class:
            switch p.Results.ClassName
                case 'Subject'
                    obj.ClassLevel = 1;
                case 'Acquisition'
                    obj.ClassLevel = 2;
                otherwise
                    obj.ClassLevel = 3;
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
            obj.b_ignoreLoggedFiles = false;
            obj.b_saveDataBeforeFail = false;
        end
        % SETTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.pipe(obj, pipe)
            % Pipeline structure setter. If pipe is empty, create an empty
            % structure containing tasks fields.
            
            if isempty(fieldnames(pipe))
                pipe = struct('argsIn', {},'argsOut',{},'outFileName','',...
                    'inputFileName', '','b_save2File', logical.empty, 'saveFileName',...
                    '', 'opts',struct.empty,'opts_vals',struct.empty,...
                    'opts_def',struct.empty ,'name','', 'inputFrom','',...
                    'seq',[],'seqIndx',[]);
            end
            % Check if all fields exist:
            if ~all(ismember(fieldnames(pipe),fieldnames(obj.pipe)))
                error('umIToolbox:PipelineManager:InvalidInput',...
                    'The pipeline structure provided is invalid!');
            end
            obj.pipe = pipe;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setOpts(obj, varargin)
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
            idxFunc = ( strcmpi(funcName, {obj.pipe.name}) & arrayfun(@(x) any(x.seq == seq) & numel(x.seqIndx) == 1,obj.pipe)' );
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
            disp(['Optional Parameters set for function : ' obj.pipe(idxFunc).name]);
        end
        
        function state = addTask(obj,func, varargin)
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
            % Output:
            %   state (bool): FALSE, if failed to add the task to the
            %   pipeline.
            
            p = inputParser;
            addRequired(p, 'func', @(x) ischar(x) || isnumeric(x));
            addOptional(p, 'b_save2File', false, @islogical);
            addOptional(p, 'saveFileName', '', @ischar);
            addParameter(p, 'inputFrom','', @(x) ischar(x) || isnumeric(x));
            parse(p,func, varargin{:});
            
            %
            state = false;
            % Check if the function name is valid:
            idxFunc = strcmpi(p.Results.func, {obj.funcList.name});
            assert(any(idxFunc),['Operation aborted! The function "' p.Results.func '" does not exist!']);
            
            % Create "task" structure. This is the one that will be added
            % to the pipeline:
            task = obj.funcList(idxFunc).info;
            task.name = obj.funcList(idxFunc).name;
            task.inputFileName = '';
            % Add default values for fields used to save the data: These
            % values will be updated later:
            task.b_save2File = p.Results.b_save2File;
            task.saveFileName = p.Results.saveFileName;
            task.inputFrom = p.Results.inputFrom;
            
            % Control inputFrom:
            if ~isempty(task.inputFrom) && ~strcmpi(task.inputFrom, '_LOCAL_')
                % Check if the task function needs data:
                assert(any(strcmpi(task.argsIn, 'data')), ['The function "' task.name '" does not have "data" as input. Operation aborted!'])
                % Check if the input function generates data:
                assert(any(ismember({'outData','outFile'},obj.funcList(idxFunc).info.argsOut)),...
                    ['The function "' task.inputFrom '" does not have "data" as input. Operation aborted!'])
            end
            % Add branch and sequence index to the function:
            task = obj.setInput(task);
            if isempty(task)
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
                        
            obj.validatePipeline; % Run pipeline validation
            
            % Recheck:
            if ~obj.b_pipeIsValid
                warning('Pipeline execution aborted! Pipeline is invalid.')
                return
            end
            % Reset sequence counter:
            obj.current_seq = 1;
            lbf = matfile(obj.ProtocolObj.LogBookFile);
            obj.tmp_LogBook = lbf.LogBook;
            obj.PipelineSummary = obj.ProtocolObj.createEmptyTable;            
            % Initialize waitbars:
            obj.setWaitBar('Initialize')
            
            % Get indexes of Filtered Objects from OBJ.PROTOCOLOBJ.QUERYFILTER function.
            switch obj.ClassLevel
                case 3
                    % Modality
                    targetIdxArr = unique(obj.ProtocolObj.Idx_Filtered,'rows');
                case 2
                    % Acquisition
                    targetIdxArr = unique(obj.ProtocolObj.Idx_Filtered(:,[1 2]), 'rows');
                case 1
                    % Subject
                    targetIdxArr = unique(obj.ProtocolObj.Idx_Filtered(:,1), 'rows');
            end
            
            for ii = 1:size(targetIdxArr,1)
                % Clear current data,  metaData and File List before starting the pipeline:
                obj.current_data = []; obj.current_metaData = [];obj.current_outFile = {};
                obj.b_state = true;
                % Get handle of current object:
                obj.getTargetObj(targetIdxArr(ii,:));
                if isempty(obj.tmp_TargetObj.LastLog)
                    % Initialize Log table:
                    obj.tmp_TargetObj.LastLog = obj.ProtocolObj.createEmptyTable;
                end
                % Create Full ID of object:
                ID_list = repmat({'null'}, 1,3);
                if isa(obj.tmp_TargetObj, 'Modality')
                    ID_list = {obj.tmp_TargetObj.MyParent.MyParent.ID, obj.tmp_TargetObj.MyParent.ID, obj.tmp_TargetObj.ID};
                elseif isa(obj.tmp_TargetObj, 'Acquisition')
                    ID_list(1,[1,2]) = {obj.tmp_TargetObj.MyParent.ID, obj.tmp_TargetObj.ID};
                else
                    ID_list{1} = obj.tmp_TargetObj.ID;
                end
                obj.targetObjFullID = strjoin(ID_list, ' -- ');
                % Update waitbars:
                obj.setWaitBar('UpdateItem', ii, size(targetIdxArr,1));
                fprintf([repmat('-',1,50),'\n']);
                fprintf('Object Name: %s\n\n', obj.targetObjFullID)
                % Run one sequence at a time:
                for jj = 1:max([obj.pipe.seq],[],'all')
                    obj.current_seq = jj;
                    % Get current pipeline sequence:
                    thisSeq = obj.pipe(arrayfun(@(x) any(x.seq == obj.current_seq), obj.pipe));
                    % Skip pipeline steps:
                    if ~obj.b_ignoreLoggedFiles
                        [thisSeq, skippedFcns,selFile] = obj.skipSteps(thisSeq,obj.tmp_TargetObj.SaveFolder);
                        if isempty(thisSeq)
                            % When all sequence is skipped:
                            fprintf('All steps skipped from the current sequence\n')                            
                            % Initialize empty Log for current object:
                            LastLog = obj.ProtocolObj.createEmptyTable;
                            % Fill out Log with Subject/Acquisition/Modality IDs :
                            LastLog(1,1:3) = strsplit(obj.targetObjFullID, ' -- ');
                            % Add class name to table:
                            txt = sprintf('The file %s contains the following steps:',selFile);
                            txt = [txt,sprintf('\n%s', skippedFcns{:})];
                            LastLog(1,[4:8,11]) = {obj.ClassName, ['Skipped Sequence #' num2str(obj.current_seq)],...
                                'none',selFile,true, txt};  
                            
%                             obj.updateTargetObjLog(LastLog);
                            obj.PipelineSummary = [obj.PipelineSummary; LastLog];
                        elseif ~isempty(skippedFcns)
                            % When some steps are skipped:
                            fprintf('The following steps will be skipped:\n')
                            fprintf('\t"%s"\n',skippedFcns{:});
                        end
                    end
                    % Run pipeline sequence in each target object:
                    for kk = 1:length(thisSeq)
                        obj.current_task = thisSeq(kk);
                        obj.setWaitBar('UpdateTask', kk/length(thisSeq));
                        fprintf('Running task # %d/%d ----->>>>>\n',kk,length(thisSeq));
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
                    
                    % Remove temporary files with appended with the "timeTag":
                    obj.deleteTemporaryFiles(obj.tmp_TargetObj.SaveFolder);
                    % Abort outer loop if user cancels pipeline:
                    if ~ishandle(obj.h_wbItem)
                        break
                    end
                end
%                 % Update Pipeline summary table:
%                 obj.PipelineSummary = [obj.PipelineSummary; obj.tmp_TargetObj.LastLog];
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
            % Request User permission to save error report:
            obj.genErrorReport;
        end
        
         function validatePipeline(obj)
            % VALIDATEPIPELINE checks if the data generated in the pipeline
            % is saved to the disk. If the user forgets to save the end
            % points of the pipeline sequences, this function will force the
            % saving of the data using the function's default ouput
            % filename.
            
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
                    % Force saving last step of the sequence:
                    if ~thisSeq(jj).b_save2File
                        thisSeq(jj).b_save2File = true;
                        thisSeq(jj).saveFileName = thisSeq(jj).outFileName;
                    end
                    % Check for output files with the same name
                    idxDuplicate = strcmpi(thisSeq(jj).saveFileName,lastFileNameList);
                    if any(idxDuplicate)
                        % If a filename from a sequence already exists,
                        % ask user to change it. Otherwise, the data
                        % from one sequence will be OVERWRITTEN!!
                        [~,defName,ext] = fileparts(lastFileNameList{idxDuplicate});
                        newFileName  = [defName '_' num2str(ii), ext]; % Append sequence index to new filename
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
                % OKbutton Callback
                function okBtnPushed(source,~)
                    newFileName = ef.Value;
                    if ~endsWith(newFileName,'.dat','IgnoreCase',true) || ~endsWith(newFileName,',mat', 'IgnoreCase',true)
                        newFileName = [newFileName,ext];
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
        
        function savePipe(obj, filename)
            % SAVEPIPE saves the structure OBJ.PIPE in a .JSON file in the
            % folder PIPELINECONFIGFILES inside the SAVEDIR of OBJ.PROTOCOLOBJ.
            
            [path, file,ext] = fileparts(filename);
            if ~strcmpi(ext,'.json')
                ext = '.json';
            end
            if isempty(path)
                path = fullfile(obj.ProtocolObj.SaveDir, 'PipeLineConfigFiles');
            end
            if ~exist(path,'dir')
                [~,~] = mkdir(targetDir);
            end
            pipeStruct = obj.pipe;           
            txt = jsonencode(pipeStruct);
            fid = fopen(fullfile(path,[file ext]), 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            disp(['Pipeline saved as "' file '" in ' path]);
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
            % Check if all functions listed in "new_pipe" exist:
            idx = ismember({new_pipe.name},{obj.funcList.name});            
            if ~all(idx)
                h = errordlg(['The following functions do not exist:' sprintf('\n"%s"',new_pipe(~idx).name)], ...
                    'Failed to load pipeline!');
                waitfor(h);
                return
            end
            %%%----- RETROCOMPATIBILITY SECTION ---------------------------
            stale_fields = setdiff(fieldnames(new_pipe), fieldnames(obj.pipe));
            missing_fields = setdiff(fieldnames(obj.pipe), fieldnames(new_pipe));
            if ~isempty(stale_fields) || ~isempty(missing_fields)   
               w = warndlg('Deprecated pipeline file found. The file will be updated now!');                
               waitfor(w);
                
                for ii = 1:length(new_pipe)
                    % Add inputFrom:
                    if new_pipe(ii).b_save2File
                        try
                            state = obj.addTask(new_pipe(ii).name,true,new_pipe(ii).datFileName);
                        catch
                            state = obj.addTask(new_pipe(ii).name,true,new_pipe(ii).saveFileName);
                        end
                    else
                        state = obj.addTask(new_pipe(ii).name);
                    end
                    if ~state
                        error('Failed to load pipeline!');
                    end
                    % Update parameters
                    obj.pipe(end).opts = new_pipe(ii).opts;
                end
                % Overwrite old .JSON file:                
                obj.savePipe(filename);
                obj.validatePipeline;
                return
            end                                       
            %%%------------------------------------------------------------
            % Add new tasks:
            fn = fieldnames(obj.pipe);
            for i = 1:length(new_pipe)                
                for k = 1:numel(fn)
                    obj.pipe(i).(fn{k}) = new_pipe(i).(fn{k});
                end
            end
            disp('Pipeline loaded!')
            % Validate loaded pipeline:
            obj.validatePipeline;
        end
        
        function reset_pipe(obj)
            % This function erases the pipe property and resets the funcList
            % property to default parameter values.            
            obj.pipe = struct();                       
            % Clear current data,  metaData and File List:
            obj.current_data = []; obj.current_metaData = [];obj.current_outFile = {};
            % Reset sequence:
            obj.current_seq = 0; obj.current_seqIndx = 0;
            disp('Pipeline erased!')
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
            LastLog(:,1:3) = strsplit(obj.targetObjFullID, ' -- ');
            % Add class name to table:
            LastLog(:,4) = {obj.ClassName};
            LastLog(:,5) = {task.name};
            %%%
            % Create function string and update log table:
            task.funcStr = createFcnString(obj, task);
            LastLog.Job = {task.funcStr};
            % Add full path of input file, if applicable:
            b_hasInputFile = false;
            if strcmpi(task.inputFileName,'data')                
                % Do nothing
            elseif ~isempty(task.inputFileName)
                b_hasInputFile = true;
                LastLog.InputFile_Path = fullfile(obj.tmp_TargetObj.SaveFolder, task.inputFileName);
            end
            %  Execute the task:
            try
                % Control for missing input files:
                if b_hasInputFile
                    errID = 'MATLAB:Umitoolbox:PipelineManager:FileNotFound';
                    errmsg = ['Input File for function ' task.name ' not found!'];
                    assert(isfile(fullfile(obj.tmp_TargetObj.SaveFolder, task.inputFileName)),...
                        errID,errmsg);
                    obj.loadInputFile(task);
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
                LastLog.Messages = {getReport(ME,'extended', 'hyperlinks','off')};
                LastLog.Messages_short = {getReport(ME, 'basic','hyperlinks','off')};
                disp('FAILED!');
            end
            thisSeq = obj.pipe([obj.pipe.seq] == obj.current_seq);
            indx = find(strcmp(task.name, {thisSeq.name}));
            % Save data to file:
            if task.b_save2File && obj.b_state
                % Look for tasks with output data from the current step and
                % save the data to a .DAT or .MAT file:
                obj.saveDataToFile(thisSeq,indx,false)
            elseif obj.b_saveDataBeforeFail && ~obj.b_state
                % Look for tasks from the previous steps given that the
                % current one failed and save it to a file:                
                if indx > 1
                    obj.saveDataToFile(thisSeq,indx-1,true);
                end
            end            
            
            LastLog.Completed = obj.b_state;
            LastLog.RunDateTime = datetime('now');
            % Remove "empty" rows from LastLog:
            idx_emptyRow = all(strcmp('None',table2cell(LastLog(:,1:5))),2);
            LastLog(idx_emptyRow,:) = [];
            % Update log table of target object:
            obj.updateTargetObjLog(LastLog);
            % Update Pipeline summary table:
            obj.PipelineSummary = [obj.PipelineSummary; LastLog];
            
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
        
        function [newSeq,skippedSteps, selFile] = skipSteps(obj,thisSeq,folderName)
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
           %        folderName(char): full path to the folder containing
           %            files.
           %    Outputs:
           %        newSeq (struct): updated sequence without the redundant
           %            steps.
           %        skippedSteps (cell): list of skipped function names from
           %            "thisSeq".  
                                           
           % Get list of valid files in the folder:
           fileList = getFileList(folderName, 'all');
           %
           newSeq = thisSeq; 
           skippedSteps = {};
           newSeqArr = cell(size(fileList));
           selFile = '';
           skippedStepsArr = newSeqArr;
           % Compare dataHistory with pipeline sequence:
           for ii = 1:length(fileList)
               [newSeqArr{ii}, skippedStepsArr{ii}] = compareDataHistory(obj,thisSeq,fullfile(folderName,fileList{ii}));
           end
           % Check for matches:
           if all(cellfun(@isempty,skippedStepsArr))
               return
           end
           % Select file with the largest number of steps:
           % Check if there is a file that contains all sequence:
           idxSkipAll = cellfun(@isempty,newSeqArr);
           if any(idxSkipAll)
              indxSkip = find(idxSkipAll,1,'first');
              selFile = fileList{indxSkip};
              step.inputFileName = selFile;
              step.saveFileName = thisSeq(end).saveFileName;
              step.name = thisSeq(end).name;
              obj.loadInputFile(step);% Load step;
              if ~isempty(step.saveFileName) && ~strcmpi(fileList{indxSkip}, step.saveFileName)
                obj.saveDataToFile(thisSeq, length(thisSeq),false); %Save to file 
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
               % sequence.  IF none is found, just pick the first one that has the same dataHistory.
               for ii = length(idxMax):-1:1
                   if any(arrayfun(@(x) strcmp(fileList{idxMax(ii)},x.inputFileName),thisSeq))
                       break
                   end
               end
               idxMax = idxMax(ii);
           end
           selFile = fileList{idxMax};
           newSeq = newSeqArr{idxMax};
           skippedSteps = skippedStepsArr{idxMax};                      
           %%%%%--Local function ------------------------------------------
            function [outSeq,skipNames] = compareDataHistory(obj,seqIn, fileIn)
                % COMPAREDATAHISTORY compares the dataHistory of "fileIn"
                % with the pipeline sequence "seqIN" and outputs the
                % updated sequence "outSeq" and the list of skipped steps
                % "skipNames".                
                outSeq = seqIn;
                skipNames = {};
                % Load data history:
                [path,file,ext] = fileparts(fileIn);
                a = load(fullfile(path,[file '.mat']),'dataHistory');
                if isempty(fieldnames(a))
                    return
                end
                dataHistory = a.dataHistory; clear a
                %%% For retrocompatibility
                if ~strcmpi(fieldnames(dataHistory(1)), 'inputFileName')
                    dataHistory(1).inputFileName = '';
                end
                %%%
                % Compare datetimes:
                [~,locB] = ismember({dataHistory.name},{obj.funcList.name});
                if locB
                    idxSameDate = cellfun(@(x,y) strcmpi(datestr(x),y),{dataHistory.creationDatetime},...
                        {obj.funcList(locB).date});
                    
                    %                 if ~all(idxSameDate)
                    %                     return
                    %                 end
                end
                % Compare function names, inputFile and optional parameters :
                thisHistory = struct();
                for jj = 1:length(seqIn)
                    thisHistory(jj).name = seqIn(jj).name;
                    thisHistory(jj).inputFileName = seqIn(jj).inputFileName;
                    thisHistory(jj).opts = seqIn(jj).opts;
                end
                % If there is an input file in the first step of the pipeline sequence,
                % prepend the dataHistory to the current sequence:                
                idxFromFile = find(~cellfun(@isempty,{seqIn.inputFileName}),1,'first');                
                if idxFromFile == 1
                    % Load the input file dataHistory:
                    [~,inputFile,~] = fileparts(seqIn(idxFromFile).inputFileName);
                    try
                        dh = load(fullfile(path,[inputFile,'.mat']),'dataHistory'); dh = dh.dataHistory;
                        
                        prepend_info = struct();
                        for jj = 1:length(dh)
                            prepend_info(jj).name = dh(jj).name;
                            prepend_info(jj).inputFileName = dh(jj).inputFileName;
                            prepend_info(jj).opts = dh(jj).opts;
                            thisHistory = horzcat(prepend_info,thisHistory(idxFromFile:end));
                        end
                        clear inputFile prepend_info
                    catch
                        %
                    end
                    
                end
                % Check for the existence of consecutive equal steps:
                b_isEqual = false(size(dataHistory));                
                for jj = 1:length(dataHistory)
                    % Compare name and opts:
                    b_isEqual(jj) = (strcmpi(thisHistory(jj).name, dataHistory(jj).name) && ...
                        strcmpi(thisHistory(jj).inputFileName, dataHistory(jj).inputFileName) && ...
                        isequaln(thisHistory(jj).opts,dataHistory(jj).opts));
                    if jj == length(thisHistory) || ~b_isEqual(jj)
                        break
                    end
                end
                if ~all(b_isEqual)
                    return
                end                                                                              
                % Get indices of sequence corresponding to the dataHistory:
                [~,indxEqual] = ismember({dataHistory(b_isEqual).name},{seqIn.name});indxEqual(indxEqual == 0) = [];
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
                    outSeq(indx).inputFrom = '_LOCAL_';
                    outSeq(indx).inputFileName = [file ext];
                    % Reset sequence indices:
                    for jj = 1:length(outSeq)
                        outSeq(jj).seqIndx = jj;
                    end
                end
                % Return steps to be skipped:
                skipNames = {seqIn(indxEqual).name};                               
            end                                 
        end
        
        %%%%%%--Helpers for "addTask" method -----------------------------
        
        function task = setInput(obj, task)
            % SETINPUT selects the input to the function in "task". It
            % controls for multiple outputs and for functions with no
            % input. It updates the fields of "task" with the input
            % information.
            % !If the User cancels any of the dialogs, the function returns
            % an empty array!
            
            b_hasDataIn = any(strcmpi('data',task.argsIn));
            b_hasDataOut = any(ismember({'outData','outFile'},task.argsOut));
            
            % For the first step of the pipeline:
            if obj.current_seq == 0
                obj.current_seq = 1;
                obj.current_seqIndx = 1;
                task.seq = obj.current_seq;
                task.seqIndx = obj.current_seqIndx;
                if b_hasDataOut && b_hasDataIn
                    task.inputFrom = '_LOCAL_';
                    task.inputFileName = obj.selectInputFileName(task.inputFrom, task.name);
                elseif b_hasDataOut && ~b_hasDataIn
                    task.inputFrom = '_LOCAL_';
                end
                if task.inputFileName == 0
                    task = [];
                end
                return
            end
            % For when a new sequence needs to be created:
            
            if ~isempty(task.inputFrom) || (b_hasDataOut && ~b_hasDataIn)
                % If the User forces the creation of a new branch or
                %   if the task function generates data from the disk.
                obj.current_seq = obj.current_seq + 1;
                obj.current_seqIndx = 1; % Reset current sequence index.
                if isempty(task.inputFrom)
                    task.inputFrom = '_LOCAL_';
                elseif strcmpi(task.inputFrom, '_LOCAL_')
                    % prompt to file selection:
                    task.inputFileName = obj.selectInputFileName('_LOCAL_',task.name);
                else
                    obj.current_seqIndx = 2;
                    % Append new sequence to input function:
                    idx = strcmpi(task.inputFrom,{obj.pipe.name}) & arrayfun(@(x) any(x.seq == 1),obj.pipe)';
                    if ~any(idx)
                        % If there is no input function in the Sequence #1,
                        % abort:
                        task = [];
                        return
                    end
                    obj.pipe(idx).seq = [obj.pipe(idx).seq obj.current_seq]; % Append sequence number
                    obj.pipe(idx).seqIndx = [obj.pipe(idx).seqIndx 1]; % Append sequence index as "1".
                    % Check input file:
                    if any(strcmpi('outFile',obj.pipe(idx).argsOut))
                        % Get File from input function:
                        task.inputFileName = obj.selectInputFileName(obj.pipe(idx).name,task.name);
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
                    task = [];
                end
                return
            end
            
            % For the rest of the steps:
            if b_hasDataOut && ~b_hasDataIn
                % For functions that do not have "data" input but create an
                % output
                task.inputFrom = '_LOCAL_';
            elseif ~any(strcmpi('data',task.argsIn)) && ~any(strcmpi('outData',task.argsOut))
                % Control for functions that do not have any data input or
                % output. *Generally, these are functions that changes
                % auxiliary files such as meta Data*.
                task.inputFrom = '';
            else
                % If this is not the first step, look backwards on the pipeline
                % to find the first function with "outData" or "outFile":
                for ii = length(obj.pipe):-1:1
                    step = obj.pipe(ii);
                    inputFileName = obj.selectInputFileName(step.name, task.name);
                    if ~isempty(inputFileName) | inputFileName == 0%#ok
                        task.inputFrom = step.name;
                        task.inputFileName = inputFileName;
                        break
                    end
                end
            end
            % Update step index:
            obj.current_seqIndx = obj.current_seqIndx + 1;
            % Add sequence info to task:
            task.seq = obj.current_seq;
            task.seqIndx = obj.current_seqIndx;
            if task.inputFileName == 0
                task = [];
            end
        end
        
        function filename = selectInputFileName(obj,source, target)
            % SELECTINPUTFILENAME creates a dialog box for the User to
            % select a file that will be the input to the "task". If
            % the "task" field "inputFrom" is a function with multiple
            % files, it will show the possible outputs from the given
            % function. However, if the "inputFrom" is "_LOCAL_", the
            % prompt will show a list of existing files from the first
            % element in "protocol".
            % Input:
            %   source (char): function name from "pipe" propery or
            %   "_LOCAL_".
            filename = '';
            % For functions as "source":
            if ~strcmpi(source, '_LOCAL_')
                % Get main(first) sequence:
                seqList = obj.pipe(arrayfun(@(x) any(x.seq == obj.current_seq), obj.pipe));
                % Get function info from pipeline:
                funcInfo = seqList(strcmpi(source, {seqList.name}));
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
                        ['"' source '" as input to :' ], ['"' target '":']}, 'ListSize',[250,280],'Name', 'Select file');
                    if ~tf
                        disp('Operation cancelled by User')
                        filename = 0;
                        return
                    end
                    filename = funcInfo.outFileName{indxFile};
                end
            else
                % For "_LOCAL_" source, look inside the folder of the first
                % element from "protocol" and display the .dat and .mat
                % files:
                item = obj.ProtocolObj.extractFilteredObjects(obj.ClassLevel); item = item{1};
                % Create list of .dat and .mat files:
                fileList = getFileList(item.SaveFolder, 'all');
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
                kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts', 'object'};
                kwrds_out = {'outFile', 'outData', 'metaData'};
                if all(ismember(out.argsIn, kwrds_args)) && all(ismember(out.argsOut, kwrds_out))
                    [~,list(i).name, ~] = fileparts(list(i).name);
                    list(i).info = out;
                    list(i).info.opts_def = list(i).info.opts; % Duplicate default params.
                    obj.funcList = [obj.funcList ; list(i)];
                end
            end
            disp('Function list created!');
            function info = parseFuncFile(fcnStruct)
                info = struct('argsIn', {},'argsOut', {}, 'outFileName', '', 'opts', [],...
                    'opts_def',[],'opts_vals',[]);
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
                outStr = regexp(funcStr,'.*(?=\=)', 'match', 'once');
                out_args = regexp(outStr, '\[*(\w*)\,*(\w*)\]*', 'tokens', 'once');
                info(1).argsOut = out_args(~cellfun(@isempty, out_args));
                [~,funcName,~] = fileparts(fcnStruct.name);
                expInput = ['(?<=' funcName '\s*\().*?(?=\))'];
                str = regexp(funcStr, expInput, 'match', 'once');
                str = strip(split(str, ','));
                info.argsIn = str(~strcmp(str, 'varargin'));
                % Get Default outputs:
                expOutput = 'default_Output\s*=.*?(?=\n)';
                str = regexp(txt, expOutput, 'match', 'once');
                if isempty(str)
                    default_Output = '';
                else
                    eval(str)
                end
                info.outFileName = default_Output;
                % Parse default optional params struct:
                expOpts = 'default_opts\s*=.*?(?=\n)';
                str = regexp(txt, expOpts, 'match', 'once');
                if ~isempty(str)
                    eval(str)
                    info.opts = default_opts;
                    info.argsIn{end+1} = 'opts';
                    % Parse optional params values struct:
                    optsVals = 'opts_values\s*=.*?(?=\n)';
                    str_opts = regexp(txt, optsVals, 'match', 'once');
                    if ~isempty(str_opts)
                        eval(str_opts)
                        info.opts_vals = opts_values;
                    end
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
            if isprop(obj.tmp_TargetObj, 'RawFolder')
                argsIn = replace(task.argsIn, ["RawFolder", "SaveFolder", "data","metaData", "object"],...
                    {['''' obj.tmp_TargetObj.RawFolder ''''],['''' obj.tmp_TargetObj.SaveFolder ''''], 'obj.current_data',...
                    'obj.current_metaData', 'obj.tmp_TargetObj'});
            else
                argsIn = replace(task.argsIn, [ "SaveFolder", "data","metaData", "object"],...
                    {['''' obj.tmp_TargetObj.SaveFolder ''''], 'obj.current_data',...
                    'obj.current_metaData', 'obj.tmp_TargetObj'});
            end
            
            argsOut = replace(task.argsOut, ["outData", "metaData","outFile"],...
                {'obj.current_data', 'obj.current_metaData', 'obj.current_outFile'});
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
            curr_dtHist = genDataHistory(funcInfo, step.funcStr, step.opts,'none', step.inputFileName);
            % First, we need to know if the output is a "data", a .DAT file or a .MAT file:
            if any(strcmp(step.argsOut, 'outFile'))
                % In case the step ouput is .DAT file(s):
                
                % Get only filename instead of full path:
                [~, filenames, ext] = cellfun(@(x) fileparts(x), obj.current_outFile,...
                    'UniformOutput', false);
                if size(filenames,2) > size(filenames,1)
                    filenames = filenames';
                    ext = ext';
                end
                curr_dtHist.outputFile_list = join([filenames,ext],'');
                                             
                for i = 1:length(obj.current_outFile)
                    % Map existing metaData file to memory:
                    mtD = matfile(strrep(obj.current_outFile{i}, '.dat', '.mat'));
                    mtD.Properties.Writable = true;
                    % Create or update "dataHistory" structure:
                    if isprop(mtD, 'dataHistory')
                        mtD.dataHistory = appendDataHistory(mtD.dataHistory, curr_dtHist);
                    else
                        mtD.dataHistory = curr_dtHist;
                    end
                end
            elseif any(strcmp(step.argsOut, 'outData')) && endsWith(step.outFileName, '.mat')
                % In case of step output is .MAT file(s):
                if isfield(obj.current_data, 'dataHistory')
                    obj.current_data.dataHistory = appendDataHistory(obj.current_data.dataHistory, curr_dtHist);
                else
                    obj.current_data.dataHistory = curr_dtHist;
                end
                
            else
                % In case of step output is a data array:
                if isfield(obj.current_metaData, 'dataHistory')
                    obj.current_metaData.dataHistory = appendDataHistory(obj.current_metaData.dataHistory, curr_dtHist);
                else
                    obj.current_metaData.dataHistory = curr_dtHist;
                end
            end
            %%%%%--Local function -----------------------------------------
            function out = appendDataHistory(dh_original, new_dh)
                % This function appends the data history "dh" to
                % "current_dataHistory" property of obj.
                               
                % Account for missing fields (FOR RETROCOMPATIBILITY)
                fn = setdiff(fieldnames(new_dh),fieldnames(dh_original));
                for ii = 1:length(fn)
                    dh_original(1).(fn{ii}) = [];
                end
                % Append fields
                fn = setdiff(fieldnames(dh_original), fieldnames(new_dh));
                for ii = 1:length(fn)
                    new_dh(1).(fn{ii}) = [];
                end
                out = vertcat(dh_original,new_dh);
            end
        end
        
        function updateTargetObjLog(obj, LastLog)
           % UPDTETARGETOBJLOG appends the Log of the current pipeline to 
           % the object's "LastLog" property up to 500 rows. When the 500
           % rows are filled, each row added removes the last row.
           
           if height(obj.tmp_TargetObj.LastLog) + height(LastLog) < 500
               obj.tmp_TargetObj.LastLog = [obj.tmp_TargetObj.LastLog; LastLog];
               return
           end
           
           % For each new row, remove one from the object's "LastLog"
           % table:           
           obj.tmp_TargetObj.LastLog = [obj.tmp_TargetObj.LastLog; LastLog];
           obj.tmp_TargetObj.LastLog(1:height(obj.tmp_TargetObj.LastLog)-500,:) = [];
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
        
        function saveDataToFile(obj, seq, indxStep, b_failed)
            % This methods looks back in the pipeline sequence from "step" for tasks
            % with "data" or "stats data" as output and saves the current data to a
            % .DAT or .MAT file.
            % Input:
            %    seq(struct) : pipeline sequence with the current task in the pipeline.
            %    indxStep(scalar): index of the step of the sequence "seq"
            %    to be saved.
            %    b_failed (bool): If TRUE, this function will ignore the
            %    "b_save2File" from "step" and save the data.
            % Get the pipeline until the task in "step":
            indx = find(strcmp(seq(indxStep).name, {seq.name}));
            subPipe = seq(1:indx);
            % Look back in pipeline for steps with "data" or "stats data"
            % as output and save the current data using the task's info:
            for i = length(subPipe):-1:1
                task = subPipe(i);
                % If the pipeline failed, and the saveFileName was not set,
                % use the default file name to save the data:
                if b_failed & isempty(task.saveFileName) & ischar(task.outFileName)
                    task.saveFileName = task.outFileName;
                end
                
                if endsWith(task.saveFileName, '.dat')
                    save2Dat(fullfile(obj.tmp_TargetObj.SaveFolder,task.saveFileName),...
                        obj.current_data, obj.current_metaData);
                    return
                elseif endsWith(task.saveFileName, '.mat')
                    disp('Writing data to .MAT file ...');
                    S = obj.current_data;
                    save(fullfile(obj.tmp_TargetObj.SaveFolder,task.saveFileName),...
                        '-struct', 'S', '-v7.3');
                    disp(['Data saved in : "' fullfile(obj.tmp_TargetObj.SaveFolder,task.saveFileName) '"']);
                    return
                end
            end
        end
        
        function deleteTemporaryFiles(obj, folder)
            % DELETETEMPORARYFILES removes .dat and .mat files from
            % "folder" that were automatically generated during the
            % pipeline due to the existance of branches. The files to be
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
           folder = fullfile(obj.ProtocolObj.SaveDir,'PipelineErrorLogs');
           if ~exist(folder,'dir')
               mkdir(folder);
           end
           filename = ['error_log_' datestr(datetime(obj.timeTag,'InputFormat','_ddMMyyyyHHmmss'),'dd_mm_yyyy_HH_MM'), '.txt'];
           fid = fopen(fullfile(folder,filename), 'w');
           fprintf(fid,'%s',str);
           fclose(fid);          
           % Open system's file explorer:
           switch computer
               case 'PCWIN64'
                   % Windows
                   winopen(folder);
               case 'GLNXA64'
                   % Linux
                   system(['gnome-open ' folder]);
               case 'MACI64'
                   % MacOS
                   system(['open ' folder]);
               otherwise
                   disp(['Error log saved at: ' folder]);
           end            
        end        
        %%%%%%-- WAITBAR methods-------------------------------------------
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
        %%%%%%-------------------------------------------------------------
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
            if isnumeric(currVals{i})
                vo.Value = strjoin(arrayfun(@num2str,currVals{i},'UniformOutput',false),';');
            else
                vo.Value = num2str(currVals{i});
            end
            idx_char = cellfun(@ischar, listVals{i});
            name = strjoin(cellfun(@(x) ['"', x, '"'],listVals{i}(idx_char), 'UniformOutput', false),', ');
            if ~all(idx_char)
                vo.Tooltip = ['Type ' name ' or a number.'];
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
waitfor(dlg);
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
