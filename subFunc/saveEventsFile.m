function saveEventsFile(saveFolder, eventID, timestamps, varargin)
% SAVEEVENTSFILE creates a .MAT file containing a list of events
% identifiers (EVENTID) and their respective time stamps (TIMESTAMPS)
% in seconds. This function simply performs some validations to the inputs
% before saving them in a .MAT file.
% Inputs:
%   saveFolder: folder where to save the .MAT file.
%   eventID: array containing events indices. Transforms to UINT16.
%   timestamps: time stamps of events in seconds. Transforms to SINGLE.
%   state (optional) : state of the event (ON = 1 or OFF = 0). This is
%   mostly from TTL-type of events. Transforms to LOGICAL.
%   eventNameList (optional) : cell array containing the Names associated
%   with the indices listed in eventID. If not provided, outputs an empty
%   cell.
% Outputs:
%   "events.mat" containing eventID, state, timestamps and eventNameList.


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

% Flip arrays:
if size(eventID,1) < size(eventID,2)
    eventID = eventID';
end
if size(timestamps,1) < size(timestamps,2)
    timestamps = timestamps';
end
if isempty(state)
    state = ones(size(timestamps), 'single');
end
if isempty(eventNameList)
    list = unique(eventID);
    eventNameList = arrayfun(@num2str,list, 'UniformOutput', false);
end
% Checks for equality of lengths of eventID and timestamps
errID = 'IsaToolbox:IncompatibleArraySizes';
msg = 'IncompatibleArraySizes: eventID and timestamps must have the same length.';
assert(isequal(length(eventID), length(state), length(timestamps)), errID, msg)
% Check if the number of elements in eventNameList is equal to the unique
if ~isempty(eventNameList)
    msg = 'The unique values of eventID do not match the number of elements in eventNameList.';
    assert(isequal(numel(unique(eventID)), numel(eventNameList)), errID, msg);
end

% Save:
filename = fullfile(saveFolder, 'events.mat');
save(filename, 'eventID', 'state', 'timestamps', 'eventNameList');
disp(['Events MAT file saved in  ' saveFolder]);

end