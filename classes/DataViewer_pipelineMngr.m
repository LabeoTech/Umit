classdef DataViewer_pipelineMngr < handle
    % This class is a simpler version of a PipelineManager class from
    % umIToolbox. This class will create and manage a small analysis
    % pipeline for imaging datasets inside the GUI "DataViewer".
    %
    
    properties
        fcnDir char % Directory of the analysis functions.
        funcList struct % structure containing the info about each function in the "fcnDir"
    end
    
    properties (SetAccess = private)
        % structure containing the info of each step of the pipeline:
        pipe = struct('argsIn', {},'argsOut',{},'outFileName','','opts',struct.empty,...
            'opts_vals',struct.empty,'opts_def',struct.empty, 'name','',...
            'funcStr', '','b_save2File', logical.empty, 'datFileName', '');
        
        data % numerical array containing imaging data
        metaData % structure or matfile containing meta data associated with "data".
        SaveFolder % folder where data created will be stored (Save Directory).
        RawFolder % folder where the Raw data are stored (Pertinent for Data Import functions only!)
        outFile % list of file names created by some functions that generate .DAT files instead of data as outputs (e.g. run_ImagesClassification)
        dataHistory % structure containing the History of analysis functions applied to the current data.
    end
    
    methods
        function obj = DataViewer_pipelineMngr(data, metaData, SaveFolder, RawFolder)
            if isdeployed
                [obj.fcnDir,~,~] = fileparts(which('funcTemplate.m'));                                
                a = load(fullfile(obj.fcnDir,'deployFcnList.mat'));
                obj.funcList = a.out; % Get the structure "out" created inside the function "umitFcnReader".
            else
                rootDir = erase(mfilename('fullpath'), ['classes' filesep 'DataViewer_pipelineMngr']);
                if isempty(rootDir)
                    error('Umitoolbox environment variable not found!')
                end
                obj.fcnDir = fullfile(rootDir, 'Analysis');
                obj.createFcnList;
            end
            
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
        
        function set.pipe(obj, pipe)
            % Pipeline structure setter. If pipe is empty, create an empty
            % structure containing tasks fields.
            
            if isempty(fieldnames(pipe)) || isempty(pipe)
                pipe = struct('argsIn', {},'argsOut',{},'outFileName','','opts',struct.empty,...
                    'opts_vals',struct.empty,'opts_def',struct.empty, 'name','','funcStr', '','b_save2File', logical.empty, 'datFileName', '');
            end
            % Check if all fields exist:
            if ~all(ismember(fieldnames(pipe),fieldnames(obj.pipe)))
                error('umIToolbox:DataViewer_pipelineMngr:InvalidInput',...
                    'The pipeline structure provided is invalid!');
            end
            obj.pipe = pipe;
        end
        
        %%%%%%%%%%%%%%%%%%
        function setOpts(obj, func)
            % SETOPTS opens an INPUTDLG for entry of optional variables
            % (OPTS) of methods in the Pipeline.
            
            % Check if the function exists:
            idx = obj.check_funcName(func);
            if isempty(idx)
                return
            end
            % Check if the function has optional parameters:
            S = obj.funcList(idx).info;
            if isempty(S.opts)
                disp(['The function ' obj.funcList(idx).name ' does not have any optional parameters.']);
                return
            end
            
            % Check for the existence of "opts_values" structure:
            if isfield(S, 'opts_vals')
                disp('Optional parameter values found!')
                currVals = cellfun(@(x) S.opts.(x),fieldnames(S.opts), 'UniformOutput',false);
                defVals = cellfun(@(x) S.opts_def.(x),fieldnames(S.opts_def), 'UniformOutput',false);
                listVals = cellfun(@(x) S.opts_vals.(x),fieldnames(S.opts_vals), 'UniformOutput',false);
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
                out = buildInputDlg(obj.funcList(idx).name,fieldnames(S.opts),currVals,defVals,listVals,typeVals);               
            else
                % Older version of fcn params input dialog:  
                fields = fieldnames(S.opts);                             
                b_isNum = structfun(@(x) isnumeric(x), S.opts);
                b_isLogic = structfun(@(x) islogical(x), S.opts);
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
                definput = structfun(@(x) {num2str(x)}, S.opts);
                opts.Resize = 'on';
                out = inputdlg(prompt,dlgtitle,dims,definput,opts);
                if isempty(out)
                    disp('Operation cancelled by User');
                    return
                end                 
                for i = 1:length(out)
                    if b_isNum(i)
                        out{i} = str2double(out{i});
                    elseif b_isLogic(i)
                        out{i} = logical(str2double(out{i}));                    
                    end
                end 
                out = horzcat(fieldnames(S.opts),out);
            end
            % Save parameters to 'opts' structure:
            for i = 1:numel(fieldnames(obj.funcList(idx).info.opts))
                obj.funcList(idx).info.opts.(out{i,1}) = out{i,2};
            end
            disp(['Optional Parameters set for function : ' obj.funcList(idx).name]);
        end
        
        function addTask(obj, func, varargin)
            % This function adds an analysis function to the pipeline.
            
            % Parse Inputs:
            p = inputParser;
            addRequired(p, 'func', @(x) ischar(x) || isnumeric(x));
            addOptional(p, 'b_save2File', false, @islogical);
            addOptional(p, 'datFileName', '', @ischar);
            parse(p, func, varargin{:});
            
            % Check if the function name is valid:
            idx = obj.check_funcName(p.Results.func);
            if isempty(idx)
                return
            end
            % Create temporary structure with function info:
            funcInfo = obj.funcList(idx).info;
            funcInfo.name = obj.funcList(idx).name;
            
            % Replace Input argument names:
            argsIn = replace(funcInfo.argsIn, {'RawFolder', 'SaveFolder', 'data', 'metaData'}, ...
                {['''' obj.RawFolder ''''], ['''' obj.SaveFolder ''''], 'obj.data', 'obj.metaData'});
            % Replace output argument names:
            argsOut = replace(funcInfo.argsOut, {'metaData', 'outData', 'outFile'},...
                {'obj.metaData', 'obj.data', 'obj.outFile'});
            % Build function string:
            funcInfo.funcStr = ['[' strjoin(argsOut, ',') ']=' funcInfo.name '('...
                strjoin(argsIn, ',') ');']; funcInfo.funcStr = erase(funcInfo.funcStr, '[]=');
            % Add optional fields to funcInfo structure:
            funcInfo.b_save2File = p.Results.b_save2File; funcInfo.datFileName = p.Results.datFileName;
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
            elseif contains(funcInfo.name, {obj.pipe.name})
                % Abort, if the user tries to add the same function more
                % than once to the pipeline:
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
            end
            
            % Add step to pipeline:
            obj.pipe = [obj.pipe; funcInfo];
            
            % Control for existing data to be saved for the current
            % function:
            if obj.pipe(end).b_save2File && any(strcmp('obj.data', argsOut))
                %  If a custom filename was not provided, use the
                %  default one:
                if isempty(obj.pipe(end).datFileName)
                    obj.pipe(end).datFileName = obj.pipe(end).outFileName;
                else
                    % Add .DAT extension, if user just provided the
                    % filename:
                    [~,~,ext]= fileparts(obj.pipe(end).datFileName);
                    [~,~,ext_def] = fileparts(obj.pipe(end).outFileName);
                    if isempty(ext)
                        obj.pipe(end).datFileName = [obj.pipe(end).datFileName, ext_def];
                    end
                end
            elseif obj.pipe(end).b_save2File && ~any(strcmp('obj.data', argsOut))
                warning(['Cannot save output to .DAT file for the function'...
                    ' "%s" \nbecause it doesn''t have any data as output!'], funcInfo.name);
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
            outMsg = '';
            for i = 1:length(obj.pipe)
                waitbar(i/length(obj.pipe), h,...
                    ['Processing ' obj.pipe(i).name '... Step' num2str(i)...
                    '/' num2str(length(obj.pipe))]);
                
                try
                    % Check if step was already run:
                    b_skip = obj.checkDataHistory(obj.pipe(i));
                    if b_skip
                        disp(['Step #' num2str(i) ' skipped!'])
                        continue
                    end
                    % Load options structure in the workspace.
                    opts = obj.pipe(i).opts;%#ok the "opts" structure is used in the EVAL function.
                    % Run the function: 
                    eval(obj.pipe(i).funcStr);
                    % Update the metaData with the current function info:
                    obj.updateDataHistory(obj.pipe(i));
                    % If user wants, save the output to a file:
                    if obj.pipe(i).b_save2File
                        save2Dat(fullfile(obj.SaveFolder, obj.pipe(i).datFileName), obj.data, obj.metaData);
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
        
        function savePipe(obj)
            % SAVEPIPE saves the structure OBJ.PIPE to a .JSON file in the
            % Input:
            %   filename (char): Name of the file that will contain the
            %   pipeline structure "obj.pipe".
            
            [file, Path] = uiputfile('*.json', 'Save Pipeline as ...', 'pipeConfigFile');
            if file == 0
                disp('Operation Cancelled by User.')
                return
            end
            pipeStruct = obj.pipe;
            txt = jsonencode(pipeStruct);
            fid = fopen(fullfile(Path, file), 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            disp(['Pipeline saved as "' file '" in ' Path]);
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
                if ~isequaln(new_pipe(i).opts,obj.funcList(indx_name).info.opts)
                    obj.funcList(indx_name).info.opts = new_pipe(i).opts;
                end
                
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
            if strcmp(flag, 'all')
                obj.funcList = struct.empty;
                obj.createFcnList;
            end
        end
        
    end
    
    methods (Access = private)
        
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
        
        function createFcnList(obj)
            % This function creates a structure containing all information
            % about the analysis functions inside the "Analysis" folder.
            
            % Set Defaults:
            default_Output = '';
            default_opts = struct();
            opts_values = struct();
            list = dir(fullfile(obj.fcnDir, '\*\*.m'));            
            for i = 1:length(list)
                out = parseFuncFile(list(i));
                % Validate if all input arguments from the function are
                % "valid" inputs keywords:
                kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts'};
                kwrds_out = {'outFile', 'outData', 'metaData'};
                if all(ismember(out.argsIn, kwrds_args)) && all(ismember(out.argsOut, kwrds_out))
                    disp(list(i).name);
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
                    % Parse optional params values struct:
                    optsVals = 'opts_values\s*=.*?(?=\n)';
                    str_opts = regexp(txt, optsVals, 'match', 'once');
                    if ~isempty(str_opts)
                        eval(str_opts)
                        info.opts_vals = opts_values;
                    end
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
            funcInfo = obj.funcList(strcmp(step.name, {obj.funcList.name}));
            % Create a local structure with the function's info:
            curr_dtHist = genDataHistory(funcInfo, step.funcStr, step.opts,'none');
            
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
            
            % For retro-compatibility with data created in previous
            % versions of umIT:
            if ~isfield(obj.metaData, 'dataHistory')
                return
            end
            
            dH = obj.metaData.dataHistory(strcmp(step.name, {obj.metaData.dataHistory.name}));
            % If the function's creation date AND the function string AND optional parameters are
            % the same, we consider that the current step was already run.
            if isempty(dH)
                return
            elseif ( isequal(datetime(fcnInfo.date), dH.creationDatetime) &&...
                    strcmp(step.funcStr, dH.funcStr) ) && isequaln(step.opts, dH.opts)
                b_skip = true;
            end
        end
    end
end


% Local functions
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
    'MenuBar','none', 'ToolBar', 'none','Visible','off', 'CloseRequestFcn', @figCloseRequest);
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
            vo = uieditfield(g, 'numeric');
        otherwise
            vo = uieditfield(g);
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
            idxDef = strcmp(listVals{i},currVals{i});
            glChar = uigridlayout(vo,[length(listVals{i}),1]);
            glChar.RowHeight = repmat({20},size(listVals{i}));            
            for jj = 1:length(listVals{i})
                c = uicheckbox('Parent',glChar, 'Text', listVals{i}{jj}, 'Value', idxDef(jj));
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
if isempty(out)
    out = currOpts;
    return
end
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
        % Just displays a message to user.
        src.Visible = 'off';
        w = warndlg('Operation cancelled by User! Changes won''t be applied!', 'Set options Cancelled!', 'modal');
        waitfor(w);
        delete(src);
    end
end
