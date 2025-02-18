function outData = run_SpeckleMapping(data, varargin)
% RUN_SPECKLEMAPPING calls the function
% SPECKLEMAPPING from the IOI library (LabeoTech).
% In brief, this function calculates the spatial OR temporal standard
% deviation of a Laser Speckle Contrast Imaging dataset.
% Consult the documentation of "SpeckleMapping" function for details
%
% Inputs :
%   data = numeric array (e.g. Y,X,T image time series).
%
% Outputs:
%    out = structure containing the speckle map.


% Defaults:
default_Output = 'std_speckle.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('sType', 'Temporal', 'bSaveMap', false,'bLogScale', false);
opts_values = struct('sType',{{'Spatial', 'Temporal'}},'bSaveMap', [true,false], 'bLogScale', [true,false]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, varargin{:});
opts = p.Results.opts;
clear p
% Run SpeckleMapping function from IOI library:
disp(['Calculating ' opts.sType ' standard deviation in speckle data ...'])
outData = genDataStructure(SpeckleMapping(data, opts.sType, '', opts.bSaveMap, opts.bLogScale));
disp('Finished Speckle Mapping.')
end
