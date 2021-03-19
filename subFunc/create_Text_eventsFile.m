function create_Text_eventsFile(saveFolder, Text, dateAndTime, varargin)
% CREATE_TEXT_EVENTSFILE creates a .MAT file containing a Text and a
% timestamp.
% Inputs:
% saveFolder: folder where to save the .MAT file.
% Text: string.
% dateAndTime: datetime object with date and time associated with TEXT.
% flag (optional):
%   "-w" (default): writes a new file or overwrites existing file.
%   "-a" : appends a new entry to existing file.
% Output:
% "Text_events.mat" containing the Text events:
%   eventID: cell array containing events strings.
%   timestamps: array of datetime objects with date and time of each
%   event.

default_flag = '-w';

p = inputParser;
validText = @(x) ischar(x) || (iscell(x) && (ischar([x{:}])));
validFlag = @(x) ischar(x) && ismember(x, {'-w', '-a'});
addRequired(p, 'saveFolder', @isfolder);
addRequired(p, 'Text', validText);
addRequired(p, 'dateAndTime', @isdatetime);
addOptional(p,'flag', default_flag, validFlag);
parse(p,saveFolder, Text, dateAndTime, varargin{:});

saveFolder = p.Results.saveFolder;
eventID = p.Results.Text;
timestamps = p.Results.dateAndTime;
flag = p.Results.flag;

if ischar(eventID)
    eventID = {eventID};
end
% Flip arrays to have one column.
if size(eventID,1)<size(eventID,2)
    eventID = eventID';
end
if size(timestamps,1)<size(timestamps,2)
    timestamps = timestamps';
end
% Check for equality of lengths of eventID and timestamps:
errID = 'IsaToolbox:IncompatibleArraySizes';
msg = 'Text and DateAndTime arrays must have the same length';
assert(isequal(length(eventID),length(timestamps)), errID, msg)

filename = fullfile(saveFolder, 'Text_events.mat');
switch flag
    case '-w'
        save(filename, 'eventID', 'timestamps');
        disp(['Text events MAT file saved in  ' saveFolder]);
    case '-a'
        try
            a = load(filename);
        catch
            errID = 'IsaToolbox:MissingFile';
            msg = ['Text_events.mat not found in ' saveFolder];
            throw(MException(errID,msg))
        end
        a.eventID = [a.eventID ; eventID];
        a.timestamps = [a.timestamps ; timestamps];
        % Re-organize array in ascending datetime:
        [~,idx] = sort(a.timestamps);
        a.eventID = a.eventID(idx);
        a.timestamps = a.timestamps(idx);
        % Save data:
        save(filename,'-struct','a');
        disp(['Text events MAT file saved in  ' saveFolder]);
end
end