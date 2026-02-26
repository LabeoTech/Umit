function genEventsFromStimParameters(saveFolder, varargin)
% GENEVENTSFROMSTIMPARAMETERS is a helper function for run_ImagesClassification
% for data obtained with LabeoTech Optical Imaging Systems.
% This function creates an "events.mat" file based on the STIM info stored in the
% "StimParameters.mat" file.

% Inputs:
%   saveFolder (char): Fullpath of the folder where to save the events.mat
%   file.
%   stimStruct (struct, optional): Structure containing the "Stim"
%   variables with names starting with "Stim_". If provided, this structure
%   will be used to generate the "events.mat" file instead of the
%   "StimParameters.mat" file.

if nargin == 2
    stimFile = varargin{:};
else
    % Create events.mat from StimParameters.mat file:
    stimFile = load(fullfile(saveFolder,'StimParameters.mat'));
end
if isempty(stimFile)
    warning('Aborted event file creation. StimParameters file not found in save folder!')
    return
elseif ~any(startsWith(fieldnames(stimFile), 'stim', 'IgnoreCase',true))
    warning('Aborted event file creation. Triggers not found in StimParameters file!')
    return
end

fn = fieldnames(stimFile);
stimVar = fn(startsWith(fn, 'stim_', 'IgnoreCase', true));
disp('Creating events file...');
% Find stim source:
eventID = [];
state = [];
timestamps = [];
% eventNameList = {};
try
    freq = stimFile.Freq;
catch 
    freq = stimFile.FrameRateHz;
end
for ii = 1:length(stimVar)
    [ID,chanState,tmstmp] = getEventFromChannel(stimFile.(stimVar{ii}), freq);
    eventID = [eventID;ID];
    state = [state;chanState];
    timestamps = [timestamps;tmstmp];    
end
% Reorganize arrays chronologically:
[timestamps, order] = sort(timestamps);
eventID = eventID(order);
state = state(order);
% Save event File:
saveEventsFile(saveFolder, eventID, timestamps, state, num2cell(unique(eventID)))
end

function [ID,state,timestamps] = getEventFromChannel(data, FrameRateHz)

ID = [];
state = [];
timestamps = [];

data = data(:);  % ensure column

% -------------------------------------------------------------
% Detect all transitions in encoded eventID signal
% -------------------------------------------------------------

change_idx = find(diff(data) ~= 0) + 1;

if isempty(change_idx)
    return
end

prev_vals = data(change_idx - 1);
new_vals  = data(change_idx);

% -------------------------------------------------------------
% Build chronological transition list
% -------------------------------------------------------------

for k = 1:length(change_idx)

    idx = change_idx(k);

    % Falling edge (previous ID ends)
    if prev_vals(k) ~= 0
        ID(end+1,1)         = prev_vals(k);
        state(end+1,1)      = false;
        timestamps(end+1,1) = idx / FrameRateHz;
    end

    % Rising edge (new ID begins)
    if new_vals(k) ~= 0
        ID(end+1,1)         = new_vals(k);
        state(end+1,1)      = true;
        timestamps(end+1,1) = idx / FrameRateHz;
    end
end

% -------------------------------------------------------------
% If signal starts already active, prepend rising edge
% -------------------------------------------------------------

if data(1) ~= 0
    ID         = [data(1); ID];
    state      = [true; state];
    timestamps = [1/FrameRateHz; timestamps];
end

% -------------------------------------------------------------
% If signal ends active, append falling edge
% -------------------------------------------------------------

if data(end) ~= 0
    ID(end+1,1)         = data(end);
    state(end+1,1)      = false;
    timestamps(end+1,1) = length(data)/FrameRateHz;
end

% -------------------------------------------------------------
% Ensure strict alternation (true/false/true/false…)
% -------------------------------------------------------------

% Sort chronologically (safety)
[timestamps, order] = sort(timestamps);
ID    = ID(order);
state = state(order);

% Remove accidental duplicates if any (same timestamp & same state)
dup = [false; diff(timestamps)==0 & diff(state)==0];
timestamps(dup) = [];
ID(dup)         = [];
state(dup)      = [];

end