function [outData, metaData] = run_SpeckleMapping(SaveFolder, varargin)
% RUN_SPECKLEMAPPING calls the function
% SPECKLEMAPPING from the IOI library (LabeoTech).
% In brief, this function calculates the spatial OR temporal standard
% deviation of a Laser Speckle Contrast Imaging dataset.
% Consult the documentation of SpeckleMapping function for details on
% inputs.

% Defaults:
default_Output = 'std_speckle.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('sType', 'Temporal', 'channel', 'speckle', 'bSaveMap', false,...
    'bLogScale', false);

%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, varargin{:});
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p

% Run SpeckleMapping function from IOI library:
disp(['Calculating ' opts.sType ' standard deviation in ' opts.channel '...'])
outData = SpeckleMapping(SaveFolder, opts.sType, opts.channel, opts.bSaveMap, opts.bLogScale);
origMetaData= load(fullfile(SaveFolder, [opts.channel '.mat']));
%%%%% Check the shape of data depending on the type of STD calculation:
% Create new meta data file based on function's input metaData:
% Use first and second input's dimensions:
metaData = genMetaData(outData, origMetaData.dim_names(1:2), speckleMetaData);
disp('Finished Speckle Mapping.')
end
