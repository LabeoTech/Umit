function umitFcnReader(appName)
% This function is used prior to compiling an umIT application as an .EXE
% file.
% The information of analysis functions from umIT's "Analysis" folder will be
% read and saved in a .mat file.
% This file will be used by DataViewer_pipelineMngr and PipelineManager
% classes to be able to display the optional parameters of the analysis functions.

% Try to locate umIT's "Analysis" folder:
root = erase(mfilename('fullpath'),['Analysis' filesep 'umitFcnReader']);
if isempty(root)
    % Try to locate "Analysis" folder by searching for the funcTemplate.m
    % function inside it.
    [root,~,~] = fileparts(which('funcTemplate.m'));
else
    root = fullfile(root,'Analysis');
end

% Set Defaults:
default_Output = '';
default_opts = struct();
opts_values = struct();

disp('Creating Fcn list...');
list = dir(fullfile(root, '\*\*.m'));
out = [];
for i = 1:length(list)
    fcn_info = parseFuncFile(list(i));
    if strcmpi(appName,'umIT')
        % Validate if all input arguments from the function are
        % "valid" inputs keywords:
        kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts', 'object', 'dataStat'};
        kwrds_out = {'outFile', 'outData', 'metaData', 'outDataStat'};
    elseif strcmpi(appName,'DataViewer')
        kwrds_args = {'data', 'metaData', 'SaveFolder', 'RawFolder', 'opts'};
        kwrds_out = {'outFile', 'outData', 'metaData'};
    else
        error('Unknown app! It must be either "umIT" or "DataViewer"')
    end
    if all(ismember(fcn_info.argsIn, kwrds_args)) && all(ismember(fcn_info.argsOut, kwrds_out))
        [~,list(i).name, ~] = fileparts(list(i).name);
        [~,list(i).folder,~] = fileparts(list(i).folder);
        list(i).info = fcn_info;
        list(i).info.opts_def = list(i).info.opts; % Duplicate default params.
        
        out = [out;list(i)];
    end
end

save(fullfile(root,'deployFcnList.mat'),'out')
disp('Function list Saved to .MAT file!');

% Helper function:
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
