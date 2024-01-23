function saveEventsFile(saveFolder, eventID, timestamps, varargin)
% SAVEEVENTSFILE creates a .MAT file containing a list of events
% identifiers (EVENTID) and their respective time stamps (TIMESTAMPS)
% in seconds. This function mostly performs some validations to the inputs
% before saving them in a .MAT file.

% Inputs:
%   saveFolder: folder where to save the .MAT file.
%   eventID: array containing events indices. Transforms to UINT16.
%   timestamps: time stamps of events in seconds. Transforms to SINGLE.
%   state (optional) : state of the event (ON = 1 or OFF = 0). This is
%   mostly from TTL-type of events. Transforms to LOGICAL.
%   eventNameList (optional) : cell array containing the Names associated
%   with the indices listed in eventID. If not provided, outputs the list
%   of unique event IDs.

p = inputParser;
validState = @(x) (all(ismember(x, [0 1]))) | islogical(x);
addRequired(p, 'saveFolder', @isfolder);
addRequired(p, 'eventID', @isnumeric);
addRequired(p, 'timestamps', @isnumeric);
addOptional(p, 'state', [], validState);
addOptional(p, 'eventNameList', [], @iscell);
parse(p,saveFolder, eventID, timestamps, varargin{:});

saveFolder = p.Results.saveFolder;
eventID = uint16(p.Results.eventID);
timestamps = single(p.Results.timestamps);
state = logical(p.Results.state);
eventNameList = p.Results.eventNameList;
clear p

% Flip arrays:
if size(eventID,1) < size(eventID,2)
    eventID = eventID';
end
if size(timestamps,1) < size(timestamps,2)
    timestamps = timestamps';
end
if isempty(state)
    state = true(size(timestamps));
end
if isempty(eventNameList)
    eventNameList = num2cell(unique(eventID));
end
% Transform numerical eventNameList to string:
eventNameList = cellfun(@num2str,eventNameList, 'UniformOutput', false);
% Checks for equality of lengths of eventID and timestamps
errID = 'Umitoolbox:saveEventsFile:IncompatibleArraySizes';
msg = 'IncompatibleArraySizes: eventID and timestamps must have the same length.';
assert(isequal(length(eventID), length(state), length(timestamps)), errID, msg)
% Check if the number of elements in eventNameList is equal to the unique
if ~isempty(eventNameList)
    msg = 'The unique values of eventID do not match the number of elements in eventNameList.';
    assert(isequal(numel(unique(eventID)), numel(eventNameList)), errID, msg);
end
% Ensure that the event name list is a single column:
if size(eventNameList,1) < size(eventNameList,2)
    eventNameList = eventNameList';
end
if size(state,1) < size(state,2)
    state = state';
end
% Save:
filename = fullfile(saveFolder, 'events.mat');
save(filename, 'eventID', 'state', 'timestamps', 'eventNameList');
disp(['Events MAT file saved in  ' saveFolder]);

end