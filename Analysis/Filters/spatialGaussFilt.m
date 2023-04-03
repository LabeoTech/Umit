function outData = spatialGaussFilt(data, metaData, varargin)
% SPATIALGAUSSFILT performs a spatial gaussian filter in an image time
% series.


% Limitations:
% The data must be an Image time series OR a Correlation map with
% dimensions: {Y,X,T}.

% Defaults:
default_Output = 'spatFilt.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('Sigma',1);
opts_values = struct('Sigma',[0,Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p,'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
outData = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
% Validate if "data" is an Image Time Series:
errID = 'umIToolbox:spatialGaussFilt:InvalidInput';
errMsg = 'Wrong Input Data type. Data must contain dimensions "X" and "Y"';
assert(all(ismember({'Y', 'X'}, metaData.dim_names)), errID, errMsg);
assert(opts.Sigma>0, errID,'Sigma must be a positive value!');
% Apply gaussian filter:
disp('Filtering frames...')
% Permute the data so the X and Y are the 1st and 2nd dimensions:
[~,xyIndx] = ismember({'Y','X'},metaData.dim_names);
permIndx = [xyIndx, setdiff(1:length(metaData.dim_names),xyIndx)];
outData = permute(outData,permIndx);
dataSz = size(outData);
outData = reshape(outData,dataSz(1), dataSz(2), []);
for ii = 1:size(outData,3)
    outData(:,:,ii) = imgaussfilt(outData(:,:,ii), opts.Sigma, 'FilterDomain', 'spatial'); % Forced FilterDomain to "spatial" to avoid problems when the data has Infs or NaNs.);
end
outData = reshape(outData,dataSz);
outData = ipermute(outData,permIndx);
disp('Finished with Spatial Gaussian filter.')      
end


