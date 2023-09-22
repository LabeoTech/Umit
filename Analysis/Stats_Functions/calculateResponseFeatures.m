function outData = calculateResponseFeatures(data, varargin)
% CALCULATERESPONSEFEATURES calculates the response peak amplitude latencies and the area
% under the curve (AUC) for a given response signal. A "response" is defined as a
% signal that crosses a baseline level by a threshold (1 std by default).
%
% The input data for this function should be the output of the function
% "getDataFromROI," representing temporal response profiles of each ROI to a stimulus.
%
% Inputs:
%   data (struct): ROI data with temporal signals split by event.
%   Parameters:
%       - data (struct): ROI data with temporal signal split by event.
%       - opts (struct):
%           - ResponsePolarity (str): Polarity of the response ('positive' or 'negative').
%           - STD_threshold (scalar): Standard deviation factor used to set the threshold for
%             response onset and offset detection.
%           - TimeWindow_sec (scalar, [start, end], or 'all'): Post-event time window (in
%             seconds) to be used to calculate the response features. Use 'all' to consider
%             the entire response, or provide a range of values separated by a semicolon as "start;end".
%
% Output:
%   outData (struct): Structure containing extracted response features for ROIs.
%
% Defaults:
%   - ResponsePolarity: 'positive'
%   - STD_threshold: 1
%   - TimeWindow_sec: 'all'
%
% Example usage:
%   data = getDataFromROI(...); % Obtain ROI data
%   % Using the entire post-event time:
%   opts = struct('ResponsePolarity', 'positive', 'STD_threshold', 2, 'TimeWindow_sec', 'all');
%   responseFeatures = calculateResponseFeatures(data, opts);
%   % Using a specific time window (e.g., 2 seconds starting at 1.5 seconds post-event):
%   opts = struct('ResponsePolarity', 'negative', 'STD_threshold', 1.5, 'TimeWindow_sec', '1.5;3.5');
%   responseFeatures = calculateResponseFeatures(data, opts);

% Defaults:
dependency = 'getDataFromROI'; %#ok Dependent function that will be automatically added to the pipeline before this one.
default_Output = 'RespFeatures.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('STD_threshold', 1, 'ResponsePolarity', 'positive', 'TimeWindow_sec', 'all');
opts_values = struct('STD_threshold', [eps Inf], 'ResponsePolarity', {{'positive', 'negative'}}, 'TimeWindow_sec',{{'all'}}); %#ok This is here only as a reference for PIPELINEMANAGER.m.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isstruct(x)); % Validate if the input is a structure:
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x));
% Parse inputs:
parse(p,data, varargin{:});
% Initialize Variables:
outData = p.Results.data;
opts = p.Results.opts;
clear p
%%%%%%%%%%%%%%%%
% Check if the input data is "time vector split by events":
assert(all(ismember(upper(outData.dim_names), {'O','E','T'})),...
    'Wrong input data type. Data must be a time vector split by event(s).');
% For each ROI, calculate the average response and use the average time
% vector to calculate the peak stats:
ROIdata = struct();
evntList = unique(outData.eventID);
% Get frame list of time window:
evntFr = round(outData.preEventTime_sec*outData.Freq);
if strcmpi(opts.TimeWindow_sec,'all')
    frOn = evntFr+1;
    frOff = size(outData.data{1},2) - evntFr;
else       
    % Get frames values:
    frOn = round((outData.preEventTime_sec + opts.TimeWindow_sec(1))*outData.Freq);
    frOff = round((outData.preEventTime_sec + opts.TimeWindow_sec(2))*outData.Freq);
    % Reset to default if input values are out of range
    if frOn > size(outData.data{1},2) || frOff > size(outData.data{1},2) || frOn > frOff
        warning('TimeWindow onset is out of range! Reset to default ("all")')
        frOn = evntFr +1;
        frOff = size(outData.data{1},2);    
    end
end
% locate time and event dimensions:
idxEdim = find(strcmpi(outData.dim_names, 'E'));
% idxTdim = find(strcmpi(outData.dim_names), 'T');
for ii = 1:length(outData.data)
    % Instantiate output arrays:
    PeakAmp_arr = nan(size(evntList),'single');
    PeakLat_arr = PeakAmp_arr;
    onsetAmp_arr = PeakAmp_arr;
    onsetLat_arr = PeakAmp_arr;
    AUCamp_arr = PeakAmp_arr;
    avgAmp_arr = PeakAmp_arr;
    
    for jj = 1:length(evntList)
        idx = ( outData.eventID == evntList(jj) );
        % Calculate average response to the current event:
        avgResp = squeeze(mean(outData.data{ii}(:,idx,:),idxEdim,'omitnan'));
        if strcmpi(opts.ResponsePolarity, 'negative')
            % Flip the data around its mean to make the response peak
            % positive:
            avgResp = -1.*avgResp;
        end                
        % Calculate threshold as a multiple of the standard deviation of the
        %   baseline period from the average response:
        avgBsln = mean(avgResp(1:evntFr),'omitnan'); % Average baseline amplitude
        thr = avgBsln + opts.STD_threshold*std(avgResp(1:evntFr),0,1, 'omitnan');
        % Calculate response features inside the time window:
        
        % Find the maximum ('peak') response value:
        [PeakValue,indxPeak] = max(avgResp(frOn:frOff),[],1);
        indxPeak = frOn + indxPeak;
        % Calculate the Average amplitude value:
        avgAmp_arr(jj) = mean(avgResp(frOn:frOff),'omitnan') - avgBsln;
        % Calculate the Area Under the Curve amplitude:
        AUCamp_arr(jj) = trapz(avgResp(frOn:frOff)) - trapz(avgResp(1:evntFr));
        % Calculate Peak amplitude:
        PeakAmp_arr(jj) = PeakValue - avgBsln;
        if PeakValue > thr
            % Calculate the onset(threshold crossing point) amplitude and latencies in seconds:
            onsetIndx = frOn + find(avgResp(frOn:frOff) > thr,1,'first');
            % Calculate onset amplitude:
            onsetAmp_arr(jj) = avgResp(onsetIndx) - avgBsln;
            % Calculate onset latency:
            onsetLat_arr(jj) = (onsetIndx - evntFr)/outData.Freq;
            % Calculate the Peak latency:
            PeakLat_arr(jj) = (indxPeak - evntFr)/outData.Freq;
        end
    end
    % Put data inside "ROIdata" structure:
    ROIdata(ii).PeakAmplitude = PeakAmp_arr';
    ROIdata(ii).PeakLatency = PeakLat_arr';
    ROIdata(ii).TimeWindowAverageAmplitude = avgAmp_arr';
    ROIdata(ii).AUCamplitude = AUCamp_arr';
    ROIdata(ii).OnsetAmplitude = onsetAmp_arr';
    ROIdata(ii).OnsetLatency = onsetLat_arr';
end
% Update metaData:
outData.eventID = evntList;
outData.label = outData.eventNameList;
outData.data = ROIdata;
outData.dim_names = {'O','E'};
end
