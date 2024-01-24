function outData = calcResponsePeak(data, varargin)
% CALCRESPONSEPEAK calculates the amplitude, peak and onset latencies of a
% response. Here, a "response" means a signal that crosses the "baseline"
% level by a threshold (1 std by default). 
% The input data for this function must be the output of the function
% "getDataFromROI". Thus, the input data is a time vector split by events
% representing the temporal response profile of each ROI to a given stimulus.

% Inputs:
%   data: structure containing the data and meta data created by the
%       "getDataFromROI" function.
%   Parameters:
%       - ResponsePolarity (positive or negative). If positive, the
%       amplitude will be calculated using the maximum response. If
%       negative, the minimum will be used.
%       - TimeWindow_sec (all | positive pair of numbers): Post-event onset
%       Time window where to look for the peak. If "all", the whole
%       post-event time will be used. To set a narrow window, type the
%       start and end (in seconds) separated by a semicolon. Ex:
%       To look for the peak only between 0.5 and 1.2 seconds after the
%       onset of the stimulus, type: 0.5;1.2.
% Output:
%   outData: structure containing stats-ready data extracted from ROIs.

% Defaults:
default_Output = 'PeakAmp.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('ResponsePolarity','positive', 'TimeWindow_sec', 'all');
opts_values = struct('ResponsePolarity', {{'positive','negative'}}, 'TimeWindow_sec',{{'all',[0 Inf]}});%#ok This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isstruct(x)); % Validate if the input is a structure:
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality')); % Used by the umIToobox app ONLY.

% Parse inputs:
parse(p,data, varargin{:});
% Initialize Variables:
outData = p.Results.data;
opts = p.Results.opts;
% object = p.Results.object;
clear p
%%%%%%%%%%%%%%%%
% Check if the input data is "time vector split by events":
assert(all(ismember(upper(outData.dim_names), {'O','E','T'})),...
    'Wrong input data type. Data must be a time vector split by event(s).');
if numel(opts.TimeWindow_sec)~= 2 || diff(opts.TimeWindow_sec < 0)
    error('Invalid time window. The input must be a pair of increasing positive numbers.');
end
% locate time and event dimensions:
idxEdim = find(strcmpi(outData.dim_names,'E'));
idxTdim = find(strcmpi(outData.dim_names,'T'));

% For each ROI, calculate the average response and use the average time
% vector to calculate the peak stats:
ROIdata = cell(length(outData.data),1);
evntList = unique(outData.eventID);
% Get frame list of time window:
evntFr = round(outData.preEventTime_sec*outData.Freq);
if strcmpi(opts.TimeWindow_sec,'all')
    frOn = evntFr+1;
    frOff = size(outData.data{1},idxTdim) - evntFr;
else
    % Get frames values:
    frOn = round(opts.TimeWindow_sec(1)*outData.Freq) + evntFr;
    frOff = round(opts.TimeWindow_sec(2)*outData.Freq) + evntFr;
    % Reset to default if input values are out of range
    if frOn > size(outData.data{1},idxTdim) || frOff > size(outData.data{1},idxTdim) || frOn > frOff
        warning('TimeWindow onset is out of range! Reset to default ("all")')
        frOn = evntFr +1;
        frOff = size(outData.data{1},idxTdim);
    end
end

for ii = 1:length(outData.data)
        
    % Instantiate output arrays:
    PeakAmp_arr = nan(size(evntList),'single');    
    % Be sure that the data is arranged properly ('O','E','T'):
    data = outData.data{ii};
    data = permute(data,[setdiff([1:ndims(data)],[idxEdim idxTdim]),idxEdim,idxTdim]);
    for jj = 1:length(evntList)
        idx = ( outData.eventID == evntList(jj) );        
        % Calculate average response to the current event:
        avgResp = squeeze(mean(data(:,idx,:),2,'omitnan'));
        % Calculate threshold as a multiple of the standard deviation of the 
        %   baseline period from the average response:        
        avgBsln = mean(avgResp(1:evntFr),'omitnan'); % Average baseline amplitude        
        if strcmpi(opts.ResponsePolarity, 'negative')
            % Find the minimun('trough') response value that crosses the
            % threshold:
            PeakValue = min(avgResp(frOn:frOff),[],1);
        else
            % Find the maximum ('peak') response value that crosses the
            % threshold:
            PeakValue = max(avgResp(frOn:frOff),[],1);
        end
        % Calculate Peak amplitude:
        PeakAmp_arr(jj) = abs(PeakValue) - avgBsln;
    end
    % Put data inside "ROIdata" structure:
    ROIdata{ii} = PeakAmp_arr';    
end
% Update meta data:
outData.eventID = evntList;
outData.label = outData.eventNameList;
outData.data = ROIdata;
outData.dim_names = {'O','E'};

end