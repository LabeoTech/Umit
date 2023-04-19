function [outData, metaData] = run_SpeckleMapping(data, metaData, varargin)
% RUN_SPECKLEMAPPING calls the function
% SPECKLEMAPPING from the IOI library (LabeoTech).
% In brief, this function calculates the spatial OR temporal standard
% deviation of a Laser Speckle Contrast Imaging dataset.
% Consult the documentation of SpeckleMapping function for details on
% inputs.

% Defaults:
default_Output = 'std_speckle.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('sType', 'Temporal', 'bSaveMap', false,'bLogScale', false);
opts_values = struct('sType',{{'Spatial', 'Temporal'}},'bSaveMap', [true,false], 'bLogScale', [true,false]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
data = p.Results.data;
origMetaData = p.Results.metaData;
opts = p.Results.opts;
clear p

% Run SpeckleMapping function from IOI library:
disp(['Calculating ' opts.sType ' standard deviation in speckle data ...'])
outData = SpeckleMapping(data, opts.sType, '', opts.bSaveMap, opts.bLogScale);
% Create meta data:
% Change "Freq" property to zero:
origMetaData.Freq = 0;
% Use first and second input's dimensions:
new_dims = origMetaData.dim_names(1:2);
% Create new meta data file based on function's input metaData:
metaData = genMetaData(outData, new_dims, origMetaData);
disp('Finished Speckle Mapping.')
end
