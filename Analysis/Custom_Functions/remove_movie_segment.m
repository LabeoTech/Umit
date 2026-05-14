function [outData, metaData] = remove_movie_segment(data, metaData, SaveFolder, varargin)
%REMOVE_MOVIE_SEGMENT Remove one temporal segment from an image time series.
%
%   [outData, metaData] = remove_movie_segment(data, metaData, SaveFolder)
%   [outData, metaData] = remove_movie_segment(data, metaData, SaveFolder, opts)
%
%   Removes one continuous segment from the time dimension of a 3D image
%   time series. DATA can be either a numeric array or a filename pointing
%   to a disk-backed .dat file. The time dimension is identified from
%   metaData.dim_names.
%
%   Inputs:
%       data       - 3D numeric array or .dat filename.
%       metaData   - Metadata structure or matlab.io.MatFile object.
%       SaveFolder - Folder containing events.mat when event timestamp
%                    updating is requested.
%       opts       - Optional structure with fields:
%                    segment_start - Start of removed segment.
%                                    Default: 0
%                    segment_end   - End of removed segment.
%                                    Default: 0
%                    TimeUnit      - 'Frames' or 'Seconds'.
%                                    Default: 'Frames'
%                    update_events - If true, update events.mat in
%                                    SaveFolder. If events.mat is not found,
%                                    event updating is skipped with a warning.
%                                    Default: false
%                    backup_events - If true, create a timestamped backup
%                                    before overwriting events.mat.
%                                    Default: true
%
%   Outputs:
%       outData    - Segment-removed numeric array, or output .dat filename
%                    in disk-backed mode.
%       metaData   - Updated metadata.
%
%   Notes:
%       - The removed interval is inclusive: segment_start:segment_end.
%       - Disk-backed data are assumed to be stored as single precision,
%         therefore each sample is assumed to use 4 bytes.
%       - When update_events is true, events inside the removed interval are
%         deleted and events after the removed interval are shifted backward.
%       - Event timestamps are assumed to be stored in seconds in the field
%         "timestamps", consistent with downstream functions that convert
%         event timestamps to frames using metaData.Freq.
%       - events.mat is overwritten only when opts.update_events is true.
%         A backup is created first when opts.backup_events is true.

% Defaults
default_Output = 'segment_removed_mov.dat'; %#ok<NASGU>
default_opts = struct('segment_start', 0, 'segment_end', 0, 'TimeUnit', 'Frames','update_events', false, 'backup_events', true);
opts_values = struct('segment_start', [0 Inf],'segment_end', [0 Inf], 'TimeUnit', {{'Seconds','Frames'}}, 'update_events', [true,false], 'backup_events', [true,false]); %#ok<NASGU>

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) || ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'metaData', @(x) isa(x, 'matlab.io.MatFile') || isstruct(x));
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, data, metaData, SaveFolder, varargin{:});

outData    = p.Results.data;
metaData   = p.Results.metaData;
SaveFolder = char(string(p.Results.SaveFolder));
opts       = p.Results.opts;
clear p

% Allow callers to provide partial opts structs.
optNames = fieldnames(default_opts);
for iOpt = 1:numel(optNames)
    if ~isfield(opts, optNames{iOpt}) || isempty(opts.(optNames{iOpt}))
        opts.(optNames{iOpt}) = default_opts.(optNames{iOpt});
    end
end

dim_names = metaData.dim_names;

%%% Validation %%%
errID  = 'umIToolbox:remove_movie_segment:InvalidInput';
errMsg = 'Input must be a 3D matrix with dimension "T".';
assert(all(ismember(dim_names, {'Y','X','T'})), errID, errMsg);

assert(opts.segment_start >= 0, errID, 'segment_start must be >= 0.');
assert(opts.segment_end   >= 0, errID, 'segment_end must be >= 0.');
assert(islogical(opts.update_events) && isscalar(opts.update_events), errID, ...
    'update_events must be a logical scalar.');
assert(islogical(opts.backup_events) && isscalar(opts.backup_events), errID, ...
    'backup_events must be a logical scalar.');

% events.mat has a fixed required filename and must live in SaveFolder.
eventsFile = fullfile(SaveFolder, 'events.mat');

if opts.update_events && ~isfile(eventsFile)
    warning('umIToolbox:remove_movie_segment:MissingEventsFile', ...
        ['Event update will be skipped. The file "events.mat" does not ' ...
        'exist in SaveFolder: %s'], SaveFolder);
    opts.update_events = false;
end

%%% Identify T dimension %%%
idxT = find(strcmp('T', dim_names));

%%% Determine total frames %%%
if isnumeric(outData)
    nT = size(outData, idxT);
else
    inFile = char(string(outData));
    sz = [metaData.datSize, metaData.datLength];
    nT = sz(idxT);
end

%%% Convert segment units to frames %%%
if startsWith(lower(opts.TimeUnit), 'sec')
    fr_start = round(metaData.Freq * opts.segment_start);
    fr_stop  = round(metaData.Freq * opts.segment_end);
else
    fr_start = round(opts.segment_start);
    fr_stop  = round(opts.segment_end);
end

% Treat [0 0] as no-op.
if fr_start == 0 && fr_stop == 0
    disp('No movie segment selected for removal.')
    return
end

fr_start = max(fr_start, 1);
fr_stop  = min(fr_stop, nT);

assert(fr_start <= fr_stop, errID, ...
    sprintf('Invalid segment range [%d %d].', fr_start, fr_stop));

assert(~(fr_start == 1 && fr_stop == nT), errID, ...
    'The selected segment removes the entire movie.');

nRemoved = fr_stop - fr_start + 1;
keepIdx = [1:fr_start-1, fr_stop+1:nT];

%% ==========================================================
% STANDARD MODE (numeric input)
% ==========================================================
if isnumeric(outData)

    orig_dim_indx = 1:numel(dim_names);
    new_dim_indx  = [idxT setdiff(orig_dim_indx, idxT)];

    outData = permute(outData, new_dim_indx);

    data_sz = size(outData);
    outData = reshape(outData, data_sz(1), []);

    disp('Removing movie segment...')
    new_sz = data_sz;
    new_sz(1) = numel(keepIdx);

    outData = outData(keepIdx, :);
    outData = reshape(outData, new_sz);
    outData = permute(outData, [2:numel(dim_names) 1]);

    metaData = genMetaData(outData, dim_names, metaData);

    if opts.update_events
        updateEventsFileForRemovedSegment( ...
            eventsFile, ...
            fr_start, ...
            fr_stop, ...
            nRemoved, ...
            metaData, ...
            opts.backup_events);
    end

    disp('Done')
    return
end

%% ==========================================================
% RAM-SAFE MODE (disk-backed input)
% ==========================================================
disp('Removing movie segment (RAM-Safe mode)...')

% Open input file.
fidIn = fopen(inFile, 'r');
cIn = onCleanup(@() safeFclose(fidIn));
assert(fidIn > 0, errID, 'Failed to open input file.');

% Output file.
outFile = fullfile(fileparts(inFile), 'SEGMENT_REMOVED_DATA.dat');
fidOut = fopen(outFile, 'w');
cOut = onCleanup(@() safeFclose(fidOut));
assert(fidOut > 0, errID, 'Failed to create output file.');

% Data geometry. Disk-backed data are always single precision.
nY = sz(1);
nX = sz(2);
frameStride = nY * nX * 4;

% Stream kept frames.
for t = keepIdx
    fseek(fidIn, (t-1) * frameStride, 'bof');

    frame = fread(fidIn, nY*nX, ['*' metaData.Datatype]);
    assert(numel(frame) == nY*nX, errID, ...
        'Failed to read frame %d from disk-backed input.', t);

    fwrite(fidOut, frame, metaData.Datatype);
end

fclose(fidIn);
fclose(fidOut);

% Update metadata.
metaData.datLength = nT - nRemoved;
outData = outFile;
save(strrep(outFile, '.dat', '.mat'), '-struct', 'metaData')

% Optional events.mat update.
if opts.update_events
    updateEventsFileForRemovedSegment( ...
        eventsFile, ...
        fr_start, ...
        fr_stop, ...
        nRemoved, ...
        metaData, ...
        opts.backup_events);
end

disp('Done')
end

% =========================================================================
% Local helpers
% =========================================================================
function updateEventsFileForRemovedSegment(eventsFile, fr_start, fr_stop, nRemoved, metaData, backupEvents)
%UPDATEEVENTSFILEFORREMOVEDSEGMENT Update events.mat after segment removal.
%
%   updateEventsFileForRemovedSegment(eventsFile, fr_start, fr_stop, ...
%       nRemoved, metaData, backupEvents)
%
%   Updates the unique events.mat file after removing one movie segment.
%   The events.mat file is expected to contain:
%       eventID       - N-by-1 uint16 vector. Event index.
%       timestamps    - N-by-1 single vector. Event timestamps in seconds.
%       state         - N-by-1 logical vector. Trigger state, where
%                       true = HIGH/ON and false = LOW/OFF.
%       eventNameList - M-by-1 cell array of event names. Position in this
%                       array corresponds to eventID.
%
%   Events are removed using paired ON/OFF intervals for each eventID:
%       - Any ON/OFF pair with one or both edges inside the removed movie
%         segment is deleted.
%       - Any ON/OFF pair active across the removed movie segment is also
%         deleted, even when both edges are outside the removed segment.
%       - Remaining events after the removed segment are shifted backward
%         by the removed duration.
%       - eventNameList is compacted when all repetitions of an eventID are
%         removed, and remaining eventID values are remapped accordingly.
%
%   A timestamped backup of events.mat is created before overwriting when
%   backupEvents is true.

errID = 'umIToolbox:remove_movie_segment:InvalidEventsFile';

evntInfo = load(eventsFile);

requiredFields = {'eventID', 'timestamps', 'state', 'eventNameList'};
missingFields = requiredFields(~isfield(evntInfo, requiredFields));

assert(isempty(missingFields), errID, ...
    'events.mat is missing required field(s): %s.', strjoin(missingFields, ', '));

eventID = evntInfo.eventID;
timestamps = evntInfo.timestamps;
state = evntInfo.state;
eventNameList = evntInfo.eventNameList;

assert(isa(eventID, 'uint16') && isvector(eventID), errID, ...
    '"eventID" must be a uint16 vector.');
assert(isa(timestamps, 'single') && isvector(timestamps), errID, ...
    '"timestamps" must be a single vector.');
assert(islogical(state) && isvector(state), errID, ...
    '"state" must be a logical vector.');
assert(iscell(eventNameList) && isvector(eventNameList), errID, ...
    '"eventNameList" must be a cell vector.');

nEvents = numel(timestamps);
assert(numel(eventID) == nEvents && numel(state) == nEvents, errID, ...
    '"eventID", "timestamps", and "state" must have the same number of elements.');

eventID = eventID(:);
timestamps = timestamps(:);
state = state(:);
eventNameList = eventNameList(:);

assert(all(eventID >= 1) && all(double(eventID) <= numel(eventNameList)), errID, ...
    '"eventID" values must be valid indices into "eventNameList".');

% Convert removed frame interval to seconds. The upper bound is exclusive.
% This matches frame-level removal of fr_start:fr_stop.
tStart_s = single((fr_start - 1) / metaData.Freq);
tStop_s  = single(fr_stop / metaData.Freq);
dt_s     = single(nRemoved / metaData.Freq);

removeEvent = false(nEvents, 1);

% Remove complete ON/OFF pairs that overlap the deleted time interval.
uniqueEventIDs = unique(eventID, 'stable');

for iID = 1:numel(uniqueEventIDs)

    thisID = uniqueEventIDs(iID);
    idxThis = find(eventID == thisID);

    % Pairing assumes events are already ordered chronologically within each
    % eventID. If a malformed edge is detected, remove that edge to avoid
    % leaving an invalid events.mat file.
    iEdge = 1;

    while iEdge <= numel(idxThis)

        idxOn = idxThis(iEdge);

        if ~state(idxOn)
            removeEvent(idxOn) = true;
            warning('umIToolbox:remove_movie_segment:MalformedEventPair', ...
                ['Removed an unpaired OFF edge for eventID %d while ' ...
                'updating events.mat.'], thisID);
            iEdge = iEdge + 1;
            continue
        end

        if iEdge == numel(idxThis)
            removeEvent(idxOn) = true;
            warning('umIToolbox:remove_movie_segment:MalformedEventPair', ...
                ['Removed an unpaired ON edge for eventID %d while ' ...
                'updating events.mat.'], thisID);
            iEdge = iEdge + 1;
            continue
        end

        idxOff = idxThis(iEdge + 1);

        if state(idxOff)
            removeEvent(idxOn) = true;
            warning('umIToolbox:remove_movie_segment:MalformedEventPair', ...
                ['Removed an ON edge for eventID %d because the next edge ' ...
                'was not an OFF edge.'], thisID);
            iEdge = iEdge + 1;
            continue
        end

        onTime = timestamps(idxOn);
        offTime = timestamps(idxOff);

        pairOverlapsRemovedSegment = onTime < tStop_s && offTime > tStart_s;

        if pairOverlapsRemovedSegment
            removeEvent([idxOn, idxOff]) = true;
        end

        iEdge = iEdge + 2;
    end
end

% Also remove any event edge that falls inside the removed segment. This is
% redundant for well-formed ON/OFF pairs but protects against malformed data.
insideRemoved = timestamps >= tStart_s & timestamps < tStop_s;
removeEvent = removeEvent | insideRemoved;

% Shift remaining events that occur after the removed segment.
afterRemoved = timestamps >= tStop_s;
timestamps(afterRemoved) = timestamps(afterRemoved) - dt_s;

keepEvent = ~removeEvent;

eventID = eventID(keepEvent);
timestamps = timestamps(keepEvent);
state = state(keepEvent);

% Remove event names for IDs that no longer have any remaining events, then
% remap the remaining eventID values so they continue to index eventNameList.
if isempty(eventID)

    eventNameList = eventNameList([]);
    eventID = uint16.empty(0, 1);

else

    usedEventIDs = unique(eventID, 'stable');

    keepName = false(numel(eventNameList), 1);
    keepName(double(usedEventIDs)) = true;

    oldToNewID = zeros(numel(eventNameList), 1, 'uint16');
    oldToNewID(keepName) = uint16(1:nnz(keepName));

    eventID = oldToNewID(double(eventID));
    eventNameList = eventNameList(keepName);
end

evntInfo.eventID = eventID;
evntInfo.timestamps = timestamps;
evntInfo.state = state;
evntInfo.eventNameList = eventNameList;

if backupEvents
    [eventFolder, eventName, eventExt] = fileparts(eventsFile);
    backupFile = fullfile(eventFolder, ...
        sprintf('%s_before_remove_movie_segment_%s%s', ...
        eventName, datestr(now, 'yyyymmdd_HHMMSS'), eventExt));

    copyfile(eventsFile, backupFile);
end

if isempty(evntInfo.eventID)
    delete(eventsFile);

    warning('umIToolbox:remove_movie_segment:NoEventsRemaining', ...
        ['All events were removed after deleting the movie segment. ' ...
         'The original events.mat file was deleted.']);
else
    save(eventsFile, '-struct', 'evntInfo');
end
end

