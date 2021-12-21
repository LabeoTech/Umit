classdef DataViewer_pipelineMngr < handle
    % This class is a simpler version of a PipelineManager class from
    % umIToolbox. This class will create and manage a small analysis
    % pipeline for imaging datasets inside the GUI "DataViewer".
    %
    
    properties
        pipe struct
        fcnDir char
        funcList struct
    end
    
    properties (SetAccess = private)
        data % numerical array containing imaging data
        metaData % structure or matfile containing meta data associated with "data".
        SaveFolder % folder where data created will be stored (Save Directory).
        RawFolder % folder where the Raw data are stored (Pertinent for Data Import functions only!)
        outFile % list of file names created by some functions that generate .DAT files instead of data as outputs (e.g. run_ImagesClassification)
        dataHistory % structure containing the History of analysis functions applied to the current data.        
    end    
   
    methods
        function obj = DataViewer_pipelineMngr(data, metaData, SaveFolder, RawFolder)
            rootDir = getenv('Umitoolbox');
            obj.fcnDir = fullfile(rootDir, 'Analysis');
            obj.createFcnList;
            obj.data = data;
            obj.metaData = metaData;
            obj.SaveFolder = SaveFolder;
            obj.RawFolder = RawFolder; 
        end
        %%%%% SETTERS %%%%
        function set.metaData(obj, metaData)
            % This set function makes sure that relevant variables inside "metaData"
            % are kept intact throughout the analysis pipeline.
            
            % Transform matFiles in structure:
            if isa(metaData, 'matlab.io.MatFile')
                metaData = load(metaData.Properties.Source);
            end
            % For cases where "metaData" will be overwritten:
            if ~isempty(obj.metaData)
                disp('Updating metaData...');
                % Append new fields:
                newFields = setdiff(fieldnames(metaData), fieldnames(obj.metaData));
                for i = 1:numel(newFields)
                    obj.metaData.(newFields{i}) = metaData.(newFields{i});
                end
                % Update variables related to data dimension:
                obj.metaData.datSize = metaData.datSize;
                obj.metaData.datLength = metaData.datLength;
                obj.metaData.dim_names = metaData.dim_names;
                
                % Replace dataHistory field:
                if isfield(metaData, 'dataHistory')
                    obj.metaData.dataHistory = metaData.dataHistory;
                end
            else
                obj.metaData = metaData;
            end
            
        end
        %%%%%%%%%%%%%%%%%% 
        function setOpts(obj, func)
            % SETOPTS opens an INPUTDLG for entry of optional variables
            % (OPTS) of methods in the Pipeline. Output: Structure
            % containing the variables and values.
            [b_OK, idx] = obj.check_funcName(func);
            if ~b_OK
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
        
        function addTask(obj, func, varargin)
            % This function adds an analysis function to the pipeline.
            
            % Parse Inputs:
            p = inputParser;
            addRequired(p, 'func', @(x) ischar(x) || isnumeric(x));
            addOptional(p, 'b_save2Dat', false, @islogical);
            addOptional(p, 'datFileName', '', @ischar);
            parse(p, func, varargin{:});
            
            % Check if the function name is valid:
            [b_OK, idx] = obj.check_funcName(p.Results.func);
            if ~b_OK
                return
            end
            % Create temporary structure with function info:
            funcInfo = obj.funcList(idx).info;
            funcInfo.name = obj.funcList(idx).name;
            
            % Replace Input argument names:
            argsIn = replace(funcInfo.argsIn, {'RawFolder', 'SaveFolder', 'data', 'metaData'}, ...
                {['''' obj.RawFolder ''''], ['''' obj.SaveFolder ''''], 'obj.data', 'obj.metaData'});            
            % Replace output argument names:
            argsOut = replace(funcInfo.argsOut, {'outData', 'metaData', 'outDataStat', 'outFile'},...
                {'obj.data', 'obj.metaData', 'obj.data', 'obj.outFile'});
            % Build function string:
            funcInfo.funcStr = ['[' strjoin(argsOut, ',') ']=' funcInfo.name '('...
                strjoin(argsIn, ',') ');'];
            % Add optional fields to funcInfo structure:
            funcInfo.b_save2Dat = false; funcInfo.datFileName = '';
            % Add step to pipeline:
            if isempty(obj.pipe)
                % Check if all input arguments exist. If not, throw an error:
                props = setdiff(argsIn, {'opts'});
                idx_miss = false(1,length(props));
                for i = 1:length(idx_miss)
                    eval(['idx_miss(' num2str(i) ') = isempty(' props{i} ');'])
                end
                if any(idx_miss)                   
                    errID = 'umIToolbox:DataViewer_pipelineMngr:MissingInputs';
                    error(errID, 'Cannot add task! The following inputs are missing:\n"%s"', ...
                        strjoin(erase(props(idx_miss), 'obj.'), ' , '));
                end
                obj.pipe = [obj.pipe; funcInfo];
                
            elseif contains(funcInfo.name, {obj.pipe.name})
                warning('The function "%s" already exists in the Pipeline!', funcInfo.name);
                return
            else
                % Control for multiple outputs from the previous step:
                % Here, we assume that functions with multiple outputs
                % create only "Files" and not "data" in the workspace.
                % Given that, we will update the function string to load
                % one of the "Files" in the workspace before running the 
                % current step.
                
                if (numel(obj.pipe(end).outFileName) > 1) && any(strcmp(obj.pipe(end).argsOut, 'outFile'))
                    disp('Controlling for multiple outputs')
                    w = warndlg({'Previous step has multiple output files!',...
                        'Please, select one to be analysed!'});
                    waitfor(w);
                    [indx, tf] = listdlg('ListString', obj.pipe(end).outFileName,...
                        'SelectionMode','single');
                    if ~tf
                        disp('Operation cancelled by User')
                        return
                    end
        
                    % Update function string to load "data" and "metaData"
                    % from the selected file. Here, we assume that the
                    % given file is located in the "SaveFolder".
                    str = ['[obj.data, obj.metaData] = loadDatFile(fullfile(obj.SaveFolder, '''...
                        obj.pipe(end).outFileName{indx} '''));'];
                    funcInfo.funcStr = [str, funcInfo.funcStr];
                end
                % Save to Pipeline:
                obj.pipe = [obj.pipe; funcInfo];
            end
                       
            % Control for existing data to be saved for the current
            % function:
            if p.Results.b_save2Dat
                if any(strcmp('obj.data', argsOut))
                    obj.pipe(end).b_save2Dat = p.Results.b_save2Dat;
                    if isempty(p.Results.datFileName)
                        obj.pipe(end).datFileName = obj.pipe(end).outFileName;
                    else
                        obj.pipe(end).datFileName = p.Results.datFileName;
                    end
                else
                    warning(['Cannot save output to .DAT file for the function'...
                        ' "%s" \nbecause it doesn''t have any data as output!'], funcInfo.name);
                end
            else
                
            end
            disp(['Added "' obj.funcList(idx).name '" to pipeline.']); 
        end
        
        function showFuncList(obj)
            % Displays a list of analysis function from "obj.funcList" in
            % the command window.
            disp('List of available functions:');
            for i = 1:length(obj.funcList)
                fprintf('%d : %s\n', i, obj.funcList(i).name);
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
                    str = [str, sprintf('Save to file: "%s"\n', fullfile(obj.SaveFolder, obj.pipe(i).datFileName))];
                end
                str = [str, sprintf('--------------------\n')];
            end
            if nargout == 0
                disp(str)
            else
                varargout{1} = str;
            end
        end
        
        function varargout = run_pipeline(obj)
            % Very (I mean... very!) simple function to run pipeline.
            % This functions renames the input and output variables
            % contained in the function string and evaluates it
            % sequentially:            
            % Outputs (optional):
            %   outMsg (str) : Message containing the string "Pipeline
            %   Finished" or an Exception report if errors were found
            %   during execution.
                       
            disp('Running pipeline on Data...');

            h = waitbar(0, 'Initiating pipeline...');
            h.Children.Title.Interpreter = 'none';
            for i = 1:length(obj.pipe)
                waitbar(i/length(obj.pipe), h,...
                    ['Processing ' obj.pipe(i).name '... Step' num2str(i)...
                    '/' num2str(length(obj.pipe))]);
                opts = obj.pipe(i).opts; %#ok "opts" is used in an "eval" function.
                try
                    % Check if step was already run:
                    b_skip = obj.checkDataHistory(obj.pipe(i));
                    if b_skip
                        disp(['Step #' num2str(i) ' skipped!'])
                        continue
                    end
                    % Run the function:
                    eval(obj.pipe(i).funcStr);
                    % Update the metaData with the current function info:
                    obj.updateDataHistory(obj.pipe(i));
                    % If user wants, save the output to a file:
                    if obj.pipe(i).b_save2Dat
                        save2Dat(obj.pipe(i).datFileName, obj.data, obj.metaData);
                    end
                    outMsg = 'Pipeline Finished.';
                catch ME
                    outMsg = getReport(ME,'extended', 'hyperlinks', 'off');
                    break
                end
            end
            close(h);
            if nargout == 0
                disp(outMsg)
            else
                varargout{1} = outMsg;
            end
        end
        
        function reset_pipe(obj)
            % This function erases the pipe property and resets the funcList
            % property.
            obj.pipe = struct.empty;
            obj.funcList = struct.empty;
            obj.createFcnList;
        end
        
    end
    
    methods (Access = private)
        
        function [b_OK, idx_fcn]= check_funcName(obj, func)
            % This function is used by "setOpts" and "addTask" methods to
            % validate if the input "func" is valid.
            % Input
            %   func (numeric OR char) : index OR name of a function from
            %   obj.funcList.
            % Output:
            %   b_OK (bool): True if "func" exists in the list.
            %   idx_fcn(bool): index of "func" in "obj.funcList".
            b_OK = true;
            if isnumeric(func)
                idx_fcn = func == 1:length(obj.funcList);
                msg = ['Function with index # ' num2str(func)];
            else
                idx_fcn = strcmp(func, {obj.funcList.name});
                msg = ['Function "' func '"'];
            end
            if ~any(idx_fcn)
                disp([msg ' not found in the function list!']);
                b_OK = false;
            end
            
        end
        
        function createFcnList(obj)
            % This function creates a structure containing all information
            % about the analysis functions inside the "Analysis" folder.
            
            % Set Defaults:
            default_Output = '';
            default_opts = struct();
            
            disp('Creating Fcn list...');
            list = dir(fullfile(obj.fcnDir, '\*\*.m'));
            for i = 1:length(list)
                out = parseFuncFile(list(i));
                % Validate if all input arguments from the function are
                % "valid" inputs keywords:
                kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts'};
                kwrds_out = {'outFile', 'outData', 'metaData', 'outDataStat'};
                if all(ismember(out.argsIn, kwrds_args)) && all(ismember(out.argsOut, kwrds_out))
                    disp(list(i).name);
                    [~,list(i).name, ~] = fileparts(list(i).name);
                    list(i).info = out;
                    obj.funcList = [obj.funcList ; list(i)];
                end
                
            end
            disp('Function list created!');
            function info = parseFuncFile(fcnStruct)
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
        
        function updateDataHistory(obj, step)
           % This function creates or  updates the "dataHistory" structure
           % and saves the information to the metaData structure/matfile.
           % The dataHistory contains all information about the functions' 
           % parameters used to create the current "data" and when it was run.
           %
           % Input:
           %    step(struct) : current step of the pipeline;
           
           disp('Building Data History...')
           timestamp = datetime('now');
           funcInfo = obj.funcList(strcmp(step.name, {obj.funcList.name}));
           % Create a local structure with the function's info:
           curr_dtHist = struct('runDatetime', timestamp, 'name', {funcInfo.name},...
               'folder', {funcInfo.folder}, 'creationDatetime', datetime(funcInfo.date),...
               'opts', step.opts, 'funcStr', {step.funcStr}, 'outputFile_list', 'none');
           % First, we need to know if the output is a "data" or a file:
           if any(strcmp(step.argsOut, 'outFile'))
               curr_dtHist.outputFile_list = obj.outFile;
               for i = 1:length(obj.outFile)
                   % Map existing metaData file to memory:
                   mtD = matfile(strrep(obj.outFile{i}, '.dat', '.mat'));
                   mtD.Properties.Writable = true;
                   % Create or update "dataHistory" structure:
                   if isprop(mtD, 'dataHistory')
                       mtD.dataHistory = [mtD.dataHistory; curr_dtHist];
                   else
                       mtD.dataHistory = curr_dtHist;
                   end
               end
           else
               if isfield(obj.metaData, 'dataHistory')
                   obj.metaData.dataHistory = [obj.metaData.dataHistory; curr_dtHist];
               else
                   obj.metaData.dataHistory = curr_dtHist;
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
            
            disp('Checking step...');
            b_skip = false;
            % Find function info in Function List:
            fcnInfo = obj.funcList(strcmp(step.name, {obj.funcList.name}));
            % Find step info in object's dataHistory:
            dH = obj.metaData.dataHistory(strcmp(step.name, {obj.metaData.dataHistory.name}));
            % If the function's creation date AND the function string are
            % the same, we consider that the current step was already run.
            if ( isequal(datetime(fcnInfo.date), dH.creationDatetime) & strcmp(step.funcStr, dH.funcStr) )
                b_skip = true;
            end
        end
    end
end