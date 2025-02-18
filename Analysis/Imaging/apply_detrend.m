function outData = apply_detrend(data)
% APPLY_DETREND applies a linear detrend to the time domain of image time
% series or image time series split by events. To calculate the linear
% trend, this function uses some frames at the start and at the end of the
% time series to calculate the slope. If the data is an image time series
% split by events, the number of frames corresponds to the baseline time
% stored in the "baselinePeriod" variable from the "events.mat" file.
%
% Inputs:
%   data (3D num array | struct): Image time series ('Y','X','T') or
%       image time series split by events ('E', 'Y', 'X', 'T').
% Output:
%   outData (3D num array | struct): Detrended "data".

% Defaults:
default_Output = 'data_detrended.dat'; %#ok This line is here just for Pipeline management.
%%% Arguments parsing and validation %%%
% Validator for data from structure:
validateDataStructure = @(x) isstruct(x) && isDatStat(x) && isDataImageTimeSeries(x) && x.b_hasEvents;
% Parse inputs:
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x) == 3) || validateDataStructure(x)); % Validate if the input is numerical or a structure
parse(p, data);
%
outData = data;clear data p;

if isstruct(outData)
    % Detrend data already split by events:
    fn = fieldnames(outData.data);
    for ii = 1:length(fn)
        outData.data.(fn{ii}) = detrendData(outData.data.(fn{ii}),outData.baselinePeriod,outData.FrameRateHz);
    end
else
    % Detrend image time series:
    outData = detrendData(outData,0,0);
end        

disp('Finished detrend!');
end
% Local function:
function data_detrended = detrendData(data,baselinePeriod,FrameRateHz)

orig_sz = size(data);
% Here, we assume that the Time dimension is the last one.
data = reshape(data, [], orig_sz(end));
% Calculate linear trend:
disp('Detrending...');
% Check for baseline info:
if baselinePeriod
    frames = round(baselinePeriod*FrameRateHz);
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
data_detrended = data - trend + b;
data_detrended = reshape(data_detrended, orig_sz);
end