function outFile = run_HemoCompute(SaveFolder, varargin)
% RUN_HEMOCOMPUTE calls the function HEMOCOMPUTE from the IOI library (LabeoTech).
% For more info, read the HEMOCOMPUTE docstring.

% Defaults:
default_Output = {'HbO.dat', 'HbR.dat'}; % This line is here just for Pipeline management.
default_opts = struct('FilterSet', 'GCaMP', 'b_normalize', true, 'Red', true, 'Green', true, 'Amber', true);
opts_values = struct('FilterSet',{{'GCaMP', 'jRGECO' , 'none'}},'b_normalize',[true,false],'Red',[false, true], 'Green',[false, true],'Amber',[false, true]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputPafunction outFile = run_HemoCompute(SaveFolder, data, varargin)
% RUN_HEMOCOMPUTE Call HEMOCOMPUTE from the IOI library (LabeoTech).
%
% This wrapper dispatches to the RAM-safe or in-memory implementation
% depending on the type of "data", while preserving the master-branch
% physiological parametrization.
%
% Inputs
% ------
% SaveFolder : char
%   Folder containing the input reflectance files and where outputs will be
%   written.
%
% data : numeric array or char
%   Numeric input or a .dat filename used only to determine whether
%   RAM-safe mode should be enabled.
%
% opts (optional struct)
%   Structure with fields:
%       - FilterSet             : excitation/emission filter set
%       - b_normalize           : normalize input channels if needed
%       - Red                   : use red illumination
%       - Green                 : use green illumination
%       - Amber                 : use amber illumination
%       - HbT_concentration_uM  : total hemoglobin concentration (uM)
%       - StO2perc              : oxygen saturation (%)
%
% Output
% ------
% outFile : cell
%   Expected output filenames {'HbO.dat','HbR.dat'}.
%
% Notes
% -----
% This function preserves the master physiological parameters while using
% the Astrocyte RAM-safe execution path.

% Defaults:
default_Output = {'HbO.dat', 'HbR.dat'}; %#ok This line is here just for Pipeline management.
default_opts = struct( ...
    'FilterSet', genFieldsFromFile('FilterSets.mat',1), ...
    'b_normalize', true, ...
    'Red', true, ...
    'Green', true, ...
    'Amber', true, ...
    'HbT_concentration_uM', 100, ...
    'StO2perc', 60);

opts_values = struct( ... %#ok  This is here only as a reference for PIPELINEMANAGER.m.
    'FilterSet', {genFieldsFromFile('FilterSets.mat')'}, ...
    'b_normalize', [true,false], ...
    'Red', [false, true], ...
    'Green', [false, true], ...
    'Amber', [false, true], ...
    'HbT_concentration_uM', [eps, Inf], ...
    'StO2perc', [eps, 100]);

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addRequired(p, 'data', @(x) (isnumeric(x) && ndims(x) == 3) || ischar(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));

parse(p, SaveFolder, data, varargin{:});

SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = default_Output;
clear p

% Translate opts to illumination list:
fields = fieldnames(opts);
fields = setdiff(fields, {'FilterSet', 'b_normalize', 'HbT_concentration_uM', 'StO2perc'});
idx = cellfun(@(x) opts.(x), fields);
list = fields(idx)';

% Decide RAM-safe mode from input type:
b_lowRAMmode = ischar(data);

% Run HemoCompute function from IOI library:
try
    HemoCompute( ...
        SaveFolder, SaveFolder, ...
        opts.FilterSet, ...
        list, ...
        opts.b_normalize, ...
        opts.HbT_concentration_uM, ...
        opts.StO2perc, ...
        b_lowRAMmode); %#ok Result is used only to catch runtime errors.

catch ME
    ME = addCause(ME, MException( ...
        'umIToolbox:run_HemoCompute:UnkwnownError', ...
        'Function HemoCompute failed!'));
    throw(ME)
end

disp('Finished HemoCompute.')
endrser;
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, varargin{:});
%Initialize Variables:
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = default_Output;
clear p
% Translate opts to char cell array:
fields = fieldnames(opts);
fields = setdiff(fields, {'FilterSet', 'b_normalize'});
idx = cellfun(@(x) opts.(x), fields);
list = fields(idx)';

% Run HemoCompute function from IOI library:
try
    HbO = HemoCompute(SaveFolder,SaveFolder, opts.FilterSet, list, opts.b_normalize); %#ok The output is here just to catch assertion errors from HemoCompute.
catch ME
    ME = addCause(ME, MException('umIToolbox:run_HemoCompute:UnkwnownError',...
        'Function HemoCompute Failed! Check Matlab command window for messages!'));
    throw(ME)
end
disp('Finished HemoCompute.')
end