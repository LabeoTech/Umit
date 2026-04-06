function [outData, metaData] = run_SpeckleMapping(SaveFolder, data, varargin)
% RUN_SPECKLEMAPPING calls the function SPECKLEMAPPING from the IOI library.
%
% In brief, this function calculates the spatial OR temporal standard
% deviation of a Laser Speckle Contrast Imaging dataset.
%
% Inputs:
%   SaveFolder : folder containing the channel files
%   data       : numeric 3D array or a .dat filename (used to trigger RAM-safe mode)
%   opts       : optional structure with fields:
%       - sType     : 'Spatial' or 'Temporal'
%       - channel   : channel name (default: 'speckle')
%       - bSaveMap  : save TIFF map
%       - bLogScale : apply -log10 transform
%
% Outputs:
%   outData   : speckle map (numeric) or filename (RAM-safe mode)
%   metaData  : metadata associated with outData

% Defaults:
default_Output = 'std_speckle.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct('sType', 'Temporal', 'channel', 'speckle', 'bSaveMap', false, 'bLogScale', false);
opts_values = struct('sType', {{'Spatial', 'Temporal'}}, 'channel', {{'fluo_475','fluo', 'red', 'green', 'yellow', 'speckle'}}, 'bSaveMap', [true,false], 'bLogScale', [true,false]); %#ok

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addRequired(p, 'data', @(x) (isnumeric(x) && ndims(x) == 3) || ischar(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, SaveFolder, data, varargin{:});

SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p

% -------------------------------------------------------------------------
% Dispatch mode
% -------------------------------------------------------------------------
disp(['Calculating ' opts.sType ' standard deviation in ' opts.channel '...'])

b_lowRAMmode = ischar(data);

outData = SpeckleMapping(SaveFolder, opts.sType, opts.channel, ...
                         opts.bSaveMap, opts.bLogScale, b_lowRAMmode);

origMetaData = load(fullfile(SaveFolder, [opts.channel '.mat']));

% -------------------------------------------------------------------------
% Build output metadata
% -------------------------------------------------------------------------
origMetaData.Freq = 0;
new_dims = origMetaData.dim_names(1:2);

if ischar(outData)
    metaData = origMetaData;
    metaData.dim_names = new_dims;
    metaData.datLength = 1;
    metaData.Freq = 0;
    metaData.datSize = origMetaData.datSize(1:2);
    metaData.datFile = outData;
else
    metaData = genMetaData(outData, new_dims, origMetaData);
end

disp('Finished Speckle Mapping.')
end