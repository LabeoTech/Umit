function outData = apply_detrend(data, metaData)
% APPLY_DETREND applies a linear detrend to the time domain of image time
% series or image time series split by events. To calculate the linear
% trend, this function uses some frames at the start and at the end of the
% time series to calculate the slope. If the data is an image time series
% split by events, the number of frames corresponds to the baseline time
% stored in the "preEventTime_sec" variable from "metaData".
% Inputs: 
%   data (3D or 4D numerical matrix): Image time series ('Y','X','T') or
%       image time series split by events ('E', 'Y', 'X', 'T').
%   metaData (struct/matfile): structure containing the meta data associated with
%   "data".
% Output:
%   outData: Detrended "data".

% Defaults:
default_Output = 'data_detrended.dat'; %#ok This line is here just for Pipeline management.
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
% Check for baseline info:
if isfield(metaData, 'preEventTime_sec')
    frames = metaData.preEventTime_sec*metaData.Freq;
    if frames <=2
        frames = 3;
    end
    % use odd number of frames:
    if mod(frames,2) == 0
        frames = frames + 1;
    end
else
    frames = 7; % Force to 7, the number of frames to calculate linear slope.
end
delta_y = median(data(:,end-frames:end),2, 'omitnan') - median(data(:,1:frames),2,'omitnan');
delta_x =(size(data,2)- frames);
M = delta_y./delta_x; clear delta_*
b = median(data(:,1:frames),2,'omitnan');
trend = bsxfun(@times,M,linspace(-2,size(data,2)-3,...
    size(data,2))) + b;
% Remove trend if data was already normalized    
outData = data - trend + b;
outData = reshape(outData, orig_sz);
disp('Finished detrend!');
end