function [eventID, state, timestamps] = getEventsFromTTL(TTLsignal, sample_rate, varargin)
% GETEVENTSFROMTTL gets the channel, state and time stamps of a
% TTL signal(TTLSIGNAL) passing a given threshold(THRESHOLD).
% Inputs:
% TTLsignal: a nChannel x Signal matrix containing the analog TTL signal.
% sample_rate : The sample rate of the TTL signal recording.
% threshold (Optional): threshold value in Volts. If not provided, 2.5V is
% used.
% Outputs:
%   eventID: list of channels (index of rows from TTLsignal matrix).
%   state: state of the channel : 1 = rising; 0 = falling.
%   timestamps: time stamps of events in seconds.

default_threshold = 2.5;

p = inputParser;
validAnalogData = @(x) isnumeric(x) && ~isscalar(x);
validNumScal = @(x) isnumeric(x) && isscalar(x);
addRequired(p, 'TTLsignal', validAnalogData);
addRequired(p, 'sample_rate', validNumScal);
addOptional(p,'threshold', default_threshold, validNumScal);
addOptional(p, 'trigType', 'EdgeSet', @(x) ismember(x, {'EdgeSet','EdgeToggle'}));
parse(p,TTLsignal, sample_rate, varargin{:});
% Initialize variables:
data = p.Results.TTLsignal;
sr = p.Results.sample_rate;
thr = p.Results.threshold;
trigType = p.Results.trigType;
clear p
% Flips the matrix to have nChannels x nSamples: assuming that
% there are more samples than channels.
if size(data,1) > size(data,2)
    data = data';
end
% Find samples that cross the threshold (rising and falling):
szdat = size(data);
idx_rise = data < thr & [data(:,2:end) nan(szdat(1),1)] > thr;
idx_fall = data > thr & [data(:,2:end) nan(szdat(1),1)] < thr;
[chanRise,tmRise] = find(idx_rise);
tmRise = tmRise./sr; % transform sample into seconds;
[chanFall,tmFall] = find(idx_fall);

tmFall = tmFall./sr; % transform sample into seconds;
eventID = uint16([chanRise chanFall]);
timestamps = single([tmRise tmFall]);
state = [ones(1,numel(tmRise), 'uint8') zeros(1,numel(tmFall), 'uint8')];
% Sort arrays by time and flip:
[timestamps,idx] = sort(timestamps);
timestamps = timestamps';
eventID = eventID(idx)';
state = state(idx)';
% For Toggle type triggers:
if strcmpi(trigType, 'edgetoggle')
    disp('Setting toggle...');
    % Use 2nd signal onset as the trial "OFF" state:
    eventID = eventID(state==1);
    timestamps = timestamps(state==1);
    % Overwite state:
    state = ones(sum(state),1,'uint8');
    state(2:2:end) = 0;
end
end