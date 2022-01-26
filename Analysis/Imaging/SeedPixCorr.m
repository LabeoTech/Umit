function [outData, metaData] = SeedPixCorr(data, metaData, varargin)
% This function calculates a pixel-wise temporal correlation. 



% Defaults:
default_Output = 'SeedPixCorr.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('imageResizeTo', -1, 'FisherZ_transform', false);

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

% Find NaNs and replace them with zeros:
% idx_nan = isnan(data);
% data(idx_nan) = 0;
% Check if data is a 3-D matrix with dimensions 'X', 'Y' and 'T':
[idx, locB] = ismember({'Y', 'X', 'T'}, metaData.dim_names);
if ~all(idx) || any(locB(1:2)>2)
    error('Umitoolbox:SeedPixCorr:WrongInput', ...
        'Input Data must be a a 3-D matrix with dimensions "Y", "X" and "T".');
end
    
% Permute data to have 'Y', 'X', 'T' dimensions:
% data = permute(data, locB); % disable for now until ImagesClassification outputs data as {"Y","X","T"}. BrunoO 26/01/2022

% Calculate SeedPixel Correlation:
% Preserve data Aspect Ratio:
data_size = size(data);
if opts.imageResizeTo == -1
%     dataAspectRatio = 1;
    xy_size = data_size([1 2]);
    A = data;
else
    dataAspectRatio = opts.imageResizeTo/max(data_size([1 2]));
    xy_size = round(data_size([1 2]).*dataAspectRatio);
    A = imresize3(data, [xy_size(1), xy_size(2), size(data,3)], 'nearest');
end
clear data;
B = reshape(A, [], size(A,3))';
clear A
% % Put NaNs back to data:
% B(idx_nan) = NaN;
% [CM, P] = corr(B); % Removed P-Value matrix for now...
disp('Calculating correlation...');
outData = corr(B);
clear B
% Apply Z Fisher transformation to corr Data:
if opts.FisherZ_transform
    outData = atanh(outData);
end
outData = reshape(outData, [xy_size(1) xy_size(2) xy_size(1)*xy_size(2)]);
% P = single(reshape(P, [xy_size(1) xy_size(2) xy_size(1)*xy_size(2)]));% Removed P-Value matrix for now...

% Create MetaData structure:
dim_names = {'Y', 'X', 'S'};
metaData = genMetaData(outData,dim_names, metaData);
disp('Finished SPCM');
end