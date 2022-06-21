function outData = spatialGaussFilt(data, metaData, varargin)
% SPATIALGAUSSFILT performs a spatial gaussian filter in an image time
% series.


% Limitations:
% The data must be an Image time series OR a Correlation map with
% dimensions: {Y,X,T or S}.

% Defaults:
default_Output = 'spatFilt.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('KernelSize', 5,'Sigma',1);
opts_values = struct('KernelSize', [1:2:51], 'Sigma',[0,Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
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
errID = 'umIToolbox:GSR:InvalidInput';
errMsg = ['Wrong Input Data type. Data must be an Image time series with ' ...
    'dimensions "X", "Y", "T" or  a correlation map  with dimensions "X","Y","S".'];
assert(all(ismember(metaData.dim_names,{'Y', 'X', 'T'})) || ...
    all(ismember(metaData.dim_names,{'Y', 'X', 'S'})), errID, errMsg);

% Apply gaussian filter:
disp('Filtering frames...')
w = waitbar(0,'Filtering frames...');
for i = 1:size(outData,3)
    waitbar(i/size(outData,3),w);
    outData(:,:,i) = imgaussfilt(outData(:,:,i), opts.Sigma,'FilterSize', opts.KernelSize);
end
close(w)
disp('Finished with Spatial Gaussian filter.')      
end


