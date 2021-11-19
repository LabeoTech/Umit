classdef DataViewer_pipelineMngr < handle
    % This class is a simple version of a PipelineManager class from
    % umIToolbox. This class will create and manage a small analysis
    % pipeline for imaging datasets inside the GUI "funcImag_viz".
    
    properties
        pipe struct
        fcnDir char
        funcList struct
    end
    properties (Access = private)
        current_data % Path to a .dat file or a numeric matrix with imaging data.
        current_metaData % Matfile object with metadata associated with the imaging data in "current_data".
        current_folder % Path of 'SaveFolder'.
    end
    methods
        function obj = DataViewer_pipelineMngr(curr_data, curr_metaData, curr_folder)
            rootDir = getenv('Umitoolbox');
            obj.fcnDir = fullfile(rootDir, 'Analysis');
            obj.createFcnList;
            obj.current_data = curr_data;
            obj.current_metaData = curr_metaData;
            obj.current_folder = curr_folder;
        end
        
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
        
        function addTask(obj, func)
            % This function adds an analysis function to the pipeline.
            
            [b_OK, idx] = obj.check_funcName(func);
            if ~b_OK
                return
            end
            funcInfo = obj.funcList(idx).info;
            funcInfo.name = obj.funcList(idx).name;
            
            % Check if function already exists in current pipe and update function string:
            funcInfo.funcStr = strrep(funcInfo.funcStr, 'SaveFolder','obj.current_folder');
            funcInfo.funcStr = strrep(funcInfo.funcStr, 'RawFolder','obj.current_folder');
            % Replace "data", "outData" and "metaData" by obj's respective
            % properties:
            funcInfo.funcStr = strrep(funcInfo.funcStr, 'data,', 'obj.current_data,');
            funcInfo.funcStr = strrep(funcInfo.funcStr, 'metaData', 'obj.current_metaData');
            funcInfo.funcStr = strrep(funcInfo.funcStr, 'outData', 'obj.current_data');
            if isempty(obj.pipe)
                % Rebuild funcStr by replacing 'data' and 'metaData'
                % variables in the input arguments:
%                 funcInfo.inStr = strrep(funcInfo.inStr, 'data','obj.current_data');
%                 funcInfo.inStr = strrep(funcInfo.inStr, 'metaData','obj.current_metaData');
%                 funcInfo.funcStr = [funcInfo.outStr, funcInfo.inStr];
                obj.pipe = [obj.pipe; funcInfo];
            elseif contains(funcInfo.name, {obj.pipe.name})
                disp('Function already exists in Pipeline');
                return
            else
                % Change input names of function string:
                % Possible input-output matches:
%                 funcInfo.inStr = strrep(funcInfo.inStr, 'data','outData');
%                 funcInfo.inStr = strrep(funcInfo.inStr, 'file','outFile');
%                 funcInfo.funcStr = [funcInfo.outStr, funcInfo.inStr];
                % Save to Pipeline:
                obj.pipe = [obj.pipe; funcInfo];
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
                str =  [str sprintf('-->> Step # %d <<---\n', i)];
                if isempty(obj.pipe(i).opts)
                    opts = {'none'; 'none'};
                else
                    opts = [fieldnames(obj.pipe(i).opts)';...
                        cellfun(@(x) num2str(x), struct2cell(obj.pipe(i).opts), 'UniformOutput', false)'];
                end
                txt = sprintf('Function name : %s\nJob: "%s"\nOptional Parameters:\n',...
                    obj.pipe(i).name, obj.pipe(i).funcStr);
                str = [str, txt, sprintf('\t%s : %s\n', opts{:}), sprintf('-------------------\n')];                
            end
            if nargout == 0 
                disp(str)
            else 
                varargout{1} = str;
            end
        end
        
        function out = run_pipeline(obj)
            % Very (I mean... very!) simple function to run pipeline.
            
            % Output
            %   out(struct) : structure containing "data" and "metaData" of
            %   transformed data.
            
            h = waitbar(0, 'Initiating pipeline...');
            h.Children.Title.Interpreter = 'none';
            for i = 1:length(obj.pipe)
                waitbar(i/length(obj.pipe), h,...
                    ['Processing ' obj.pipe(i).name '... Step' num2str(i) '/' num2str(length(obj.pipe))]);
                opts = obj.pipe(i).opts;
                eval(obj.pipe(i).funcStr);
            end
            close(h);
%             for i = 1:numel(obj.pipe(end).output)
%                 eval(['out.' (obj.pipe(end).output{i}) ' = ' obj.pipe(end).output{i} ';'])
%             end
            out.data = obj.current_data;
            out.metaData = obj.current_metaData;
            disp('Finished pipeline!');
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
                idx_fcn = (func == 1:length(obj.funcList));
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
                disp(list(i).name);
                out = parseFuncFile(list(i));
                % Validate if all input arguments from the function are
                % "valid" inputs keywords:
                kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder'};
                kwrds_out = {'outFile', 'outData', 'metaData'};
                if all(ismember(out.args, kwrds_args)) && all(ismember(out.output, kwrds_out))
                    [~,list(i).name, ~] = fileparts(list(i).name);
                    list(i).info = out;
                    obj.funcList = [obj.funcList ; list(i)];
                    
                else
                    disp(['"' list(i).name '" isn''t valid umIT function!']);
                end
                
            end
            disp('Function list created!');
            function info = parseFuncFile(fcnStruct)
                info = struct('args', {},'outStr', '', 'inStr', '','output', {}, 'outFileName', '', 'opts', []);
                txt = fileread(fullfile(fcnStruct.folder, fcnStruct.name));
                funcStr = regexp(txt, '(?<=function\s*).*?(?=\n)', 'match', 'once');
                funcStr = erase(funcStr(1:end-1), ' '); % remove white spaces and an extra line from the function string.
                outStr = regexp(funcStr,'.*(?=\=)', 'match', 'once');
                out_args = regexp(outStr, '\[*(\w*)\,*(\w*)\]*', 'tokens', 'once');
                idx_empty = cellfun(@isempty, out_args);
                info(1).output = out_args(~idx_empty);
                info.outStr = outStr;
                [~,funcName,~] = fileparts(fcnStruct.name);
                expInput = ['(?<=' funcName '\s*\().*?(?=\))'];
                str = regexp(funcStr, expInput, 'match', 'once');
                str = strip(split(str, ','));
                info.args = setdiff(str, {'varargin'});
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
                    info.funcStr = [strrep(funcStr, 'varargin', 'opts'), ';'];
                else
                    funcStr = erase(funcStr, ',varargin');
                    info.funcStr = [funcStr, ';'];
                end
                info.inStr = erase(info.funcStr,outStr);
            end
        end
        
        
        
        
        
        
        
        
    end
end