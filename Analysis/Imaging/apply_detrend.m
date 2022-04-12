function outData = apply_detrend(data, metaData)
% APPLY_DETREND applies a linear detrend to the time domain of image time
% series or image time series split by events.
% Inputs: 
%   data (3D or 4D numerical matrix): Image time series ('Y','X','T') or
%       image time series split by events ('E', 'Y', 'X', 'T').
%   metaData (struct/matfile): structure containing the meta data associated with
%   "data".
% Output:
%   outData: Detrended "data".
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
parse(p, data, metaData);
% Further validation of inputs:
errMsg = 'Invalid Input. Data must be an Image time series!';
errID = 'umIToolbox:apply_detrend:WrongInput';
assert(all(ismember({'Y','X','T'}, p.Results.metaData.dim_names)), errID, errMsg);
% Initialize Variables:
data = p.Results.data;
metaData = p.Results.metaData;
clear p
%
orig_sz = size(data);
idx_T = strcmp('T', metaData.dim_names);
data = reshape(data, prod(orig_sz(~idx_T)), orig_sz(idx_T));
% Calculate linear trend:

disp('Detrending...');
delta_y = median(data(:,end-7:end),2, 'omitnan') - median(data(:,1:7),2,'omitnan');
delta_x =(size(data,2)- 7);
M = delta_y./delta_x; clear delta_*

trend = bsxfun(@times,M,linspace(-2,size(data,2)-3,...
    size(data,2))) + median(data(:,1:7),2,'omitnan');
% Automatic selection of normalization/subtraction-only depending on the
% average value of the data:
if mean(data,'all','omitnan') > 1 
    % Normalize "raw" data
    disp('Normalizing data...');
    outData = bsxfun(@rdivide,bsxfun(@minus,data,trend),trend);
else
    % Just remove trend if data was already normalized
    outData = bsxfun(@minus,data,trend);
end
outData = reshape(outData, orig_sz);
disp('Finished detrend!');
end