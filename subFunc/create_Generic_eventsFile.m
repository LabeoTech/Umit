function create_Generic_eventsFile(saveFolder, eventID, timestamps)
% CREATE_GENERIC_EVENTSFILE creates a .MAT file containing a list of events 
% identifiers (EVENTID) and their respective time stamps (TIMESTAMPS) 
% in seconds. This function simply performs some validations to the inputs
% before saving them in a .MAT file.
% Inputs:
%   saveFolder: folder where to save the .MAT file.
%   eventID: cell array containing events IDs. Accepted data types are
%   Text or Numeric.
%   timestamps: time stamps of events in seconds.
% Outputs:
%   "Generic_events.mat" containing the inputs
   

p = inputParser;
validEventList = @(x) iscell(x) && (isnumeric([x{:}]) ||  ischar([x{:}]));
addRequired(p, 'saveFolder', @isfolder);
addRequired(p, 'eventID_list', validEventList);
addRequired(p, 'timestamps', @isnumeric);
parse(p,saveFolder, eventID, timestamps);

saveFolder = p.Results.saveFolder;
eventID = p.Results.eventID_list; % Change name for confirmity with the other createXXXeventsFile functions.
timestamps = p.Results.timestamps;

% Transform timestamps in single:
timestamps = single(timestamps);
% Flip arrays:
if size(eventID,1) < size(eventID,2)
    eventID = eventID';
end
if size(timestamps,1) < size(timestamps,2)
    timestamps = timestamps';
end
% Checks for equality of lengths of eventID_list and timestamps 
errID = 'IsaToolbox:IncompatibleArraySizes';
msg = 'IncompatibleArraySizes: eventID_list and timestamps must have the same length.';
assert(isequal(length(eventID), length(timestamps)), errID, msg)

filename = fullfile(saveFolder, 'Generic_events.mat');
save(filename, 'eventID', 'timestamps');
disp(['Generic_events MAT file saved in  ' saveFolder]);

end