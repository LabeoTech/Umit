function [outData, metaData] = SPCM(data, metaData, varargin)
% This function creates a Seed-Pixel correlation map (SPCM) by calculating
% the Pearson's zero-lag temporal correlation between pixels.
% Inputs:
%   data (3D numerical matrix): Image time series with dimensions {'Y','X','T}.
%   metaData (struct): structure containing smeta data associated with "data".
%   FisherZ_transform(bool, optinal): If true,the Fisher's Z tranformation is
%       applied to the data. Otherwise, the data is expressed as Pearson's
%       correlation values.
% Outputs:
%   outData (3D numerical matrix): Seed-Pixel correlation maps with
%   dimensions {'Y', 'X', 'S'}. Here, 'S' is the seed dimension.
%   metaData (struct): structure containing smeta data associated with "data".


% Defaults:
default_Output = 'SPCMap.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('b_FisherZ_transform', false);
opts_values = struct('b_FisherZ_transform',[false, true]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
% Check if data is a 3-D matrix with dimensions 'X', 'Y' and 'T':
[idx, locB] = ismember({'Y', 'X', 'T'}, metaData.dim_names);
if ~all(idx) || any(locB(1:2)>2)
    error('Umitoolbox:SPCM:WrongInput', ...
        'Input Data must be a a 3-D matrix with dimensions "Y", "X" and "T".');
end
% Permute data to have 'Y', 'X', 'T' dimensions (for retrocompatibility with older versions of ImagesClassification fcn):
data = permute(data, locB);
% Calculate SeedPixel Correlation:
sz_dat = size(data);
data = reshape(data, [], size(data,3))';
disp('Calculating Pearson''s correlation...');
% Calculate Pearson's correlation:
outData = corrcoef(data);
clear data
% Apply Z Fisher transformation to corr Data:
if opts.b_FisherZ_transform
    disp('Applying Fisher''s Z transform...');    
    % Z-Fisher Transform:
    outData = atanh(outData);   
end
outData = reshape(outData, [sz_dat(1) sz_dat(2) sz_dat(1)*sz_dat(2)]);

% Create MetaData structure:
dim_names = {'Y', 'X', 'S'};
metaData = genMetaData(outData,dim_names, metaData);
disp('Finished SPCM');
end