function outFile = remove_movie_segment(SaveFolder, varargin)
%REMOVE_MOVIE_SEGMENT Remove one temporal segment from imported .dat movies.
%
%   outFile = remove_movie_segment(SaveFolder)
%   outFile = remove_movie_segment(SaveFolder, opts)
%
%   Removes one continuous temporal segment from all imported disk-backed
%   .dat movies listed in AcqInfos.mat. Channel names are read from
%   AcqInfoStream.Illumination<N>.Color and translated to imported .dat
%   filenames. Each movie is updated in RAM-safe mode by writing a temporary
%   copy, deleting the original file, and renaming the temporary copy to the
%   original filename.
%
%   Inputs:
%       SaveFolder - Folder containing AcqInfos.mat, imported .dat files,
%                    same-name channel metadata .mat files, and optionally
%                    events.mat.
%       opts       - Optional structure with fields:
%                    segment_start_sec - Start time of the removed segment,
%                                        in seconds. Default: 0
%                    segment_end_sec   - End time of the removed segment,
%                                        in seconds. Default: 0
%
%   Output:
%       outFile - Cell array containing the full paths of all imported .dat
%                 files discovered from AcqInfos.mat. This manifest is
%                 returned after validation even when segment_start_sec and
%                 segment_end_sec are both zero and no files are rewritten.
%
%   Notes:
%       - This function does not modify AcqInfos.mat. AcqInfos.mat is used
%         only to discover imported channels.
%       - The same-name metadata file associated with each .dat file is
%         updated after removing the segment, for example red.dat -> red.mat.
%       - All .dat files are assumed to store single-precision samples.
%       - The removed time interval is [segment_start_sec, segment_end_sec).
%       - Frame conversion is done separately for each channel using that
%         channel's metaData.Freq.
%       - The function validates all channel files, metadata, and segment
%         ranges before replacing any original .dat file.
%       - If segment_start_sec and segment_end_sec are both zero, the function
%         performs no file rewrite and returns the full file manifest.
%       - If events.mat exists, it is updated automatically once using
%         seconds because it is shared by the session, not by individual
%         channel files.
%       - If events.mat does not exist, event handling is skipped silently.
%       - Events overlapping the removed interval are deleted. Events after
%         the removed interval are shifted backward by the removed duration.
%       - eventNameList is compacted when all repetitions for an eventID are
%         removed, and remaining eventID values are remapped accordingly.
%       - If all events are removed, events.mat is deleted and a warning is
%         raised.
%
%   Channel filename mapping:
%       Red                         -> red.dat
%       Amber or Yellow             -> yellow.dat
%       Green                       -> green.dat
%       Fluo                        -> fluo.dat
%       Fluo #<N> <wavelength> nm   -> fluo_<wavelength>.dat

% Defaults
default_Output = {'fluo_475.dat', 'fluo_567.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; %#ok<NASGU> % Reference for PipelineManager. Actual outputs are stored in outFile.
default_opts = struct('segment_start_sec', 0, 'segment_end_sec', 0);
opts_values = struct('segment_start_sec', [0 Inf], 'segment_end_sec', [0 Inf]);%#ok<NASGU>

%%% Arguments parsing and validation %%%
p = inputParser;
p.FunctionName = 'remove_movie_segment';
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
opts = p.Results.opts;
clear p

% Allow callers to provide partial opts structs.
optNames = fieldnames(default_opts);
for iOpt = 1:numel(optNames)
    if ~isfield(opts, optNames{iOpt}) || isempty(opts.(optNames{iOpt}))
        opts.(optNames{iOpt}) = default_opts.(optNames{iOpt});
    end
end

%%% Validation %%%
errID = 'umIToolbox:remove_movie_segment:InvalidInput';

assert(isnumeric(opts.segment_start_sec) && isscalar(opts.segment_start_sec) && ...
    isfinite(opts.segment_start_sec) && opts.segment_start_sec >= 0, errID, ...
    'segment_start_sec must be a finite numeric scalar >= 0.');
assert(isnumeric(opts.segment_end_sec) && isscalar(opts.segment_end_sec) && ...
    isfinite(opts.segment_end_sec) && opts.segment_end_sec >= 0, errID, ...
    'segment_end_sec must be a finite numeric scalar >= 0.');
isNoOp = opts.segment_start_sec == 0 && opts.segment_end_sec == 0;

assert(isNoOp || opts.segment_end_sec > opts.segment_start_sec, errID, ...
    'segment_end_sec must be greater than segment_start_sec.');

acqInfoFile = fullfile(SaveFolder, 'AcqInfos.mat');
assert(isfile(acqInfoFile), errID, ...
    'AcqInfos.mat was not found in SaveFolder: %s', SaveFolder);

acqInfo = load(acqInfoFile, 'AcqInfoStream', '-mat');
assert(isfield(acqInfo, 'AcqInfoStream') && isstruct(acqInfo.AcqInfoStream), errID, ...
    'AcqInfos.mat must contain the structure AcqInfoStream.');
AcqInfoStream = acqInfo.AcqInfoStream;

illumFields = fieldnames(AcqInfoStream);
isIllumField = ~cellfun(@isempty, regexp(illumFields, '^Illumination\d+$', 'once'));
illumFields = illumFields(isIllumField);
assert(~isempty(illumFields), errID, ...
    'AcqInfoStream does not contain any Illumination<N> fields.');

% Sort Illumination<N> fields by numeric index for predictable behavior.
illumIdx = zeros(numel(illumFields), 1);
for iField = 1:numel(illumFields)
    token = regexp(illumFields{iField}, '^Illumination(\d+)$', 'tokens', 'once');
    illumIdx(iField) = str2double(token{1});
end
[~, sortIdx] = sort(illumIdx);
illumFields = illumFields(sortIdx);

nChannels = numel(illumFields);
channelInfo = repmat(struct( ...
    'illumField', '', ...
    'colorName', '', ...
    'datFileName', '', ...
    'datFile', '', ...
    'matFile', '', ...
    'metaData', struct(), ...
    'nT', [], ...
    'fr_start', [], ...
    'fr_stop', [], ...
    'nRemoved', [], ...
    'newLength', [], ...
    'tmpDatFile', ''), nChannels, 1);

for iChan = 1:nChannels

    illumField = illumFields{iChan};
    illumInfo = AcqInfoStream.(illumField);

    assert(isstruct(illumInfo) && isfield(illumInfo, 'Color'), errID, ...
        'AcqInfoStream.%s must contain the field Color.', illumField);

    colorName = char(string(illumInfo.Color));
    datFileName = colorNameToDatFileName(colorName);

    channelInfo(iChan).illumField = illumField;
    channelInfo(iChan).colorName = colorName;
    channelInfo(iChan).datFileName = datFileName;
    channelInfo(iChan).datFile = fullfile(SaveFolder, datFileName);

    [datFolder, datBaseName] = fileparts(channelInfo(iChan).datFile);
    channelInfo(iChan).matFile = fullfile(datFolder, [datBaseName '.mat']);
end

% Avoid ambiguous mappings such as Amber and Yellow both targeting yellow.dat.
datFileNames = {channelInfo.datFileName};
for iChan = 1:nChannels
    isDuplicate = strcmp(datFileNames{iChan}, datFileNames);
    assert(nnz(isDuplicate) == 1, errID, ...
        'Multiple illumination channels map to the same .dat file: %s.', ...
        datFileNames{iChan});
end

% Validate all files and metadata before creating any temporary outputs.
for iChan = 1:nChannels

    assert(isfile(channelInfo(iChan).datFile), errID, ...
        'Expected imported .dat file was not found: %s', channelInfo(iChan).datFile);
    assert(isfile(channelInfo(iChan).matFile), errID, ...
        'Expected same-name metadata .mat file was not found: %s', channelInfo(iChan).matFile);

    metaData = load(channelInfo(iChan).matFile, '-mat');

    requiredFields = {'dim_names', 'datSize', 'datLength', 'Datatype', 'Freq'};
    missingFields = requiredFields(~isfield(metaData, requiredFields));
    assert(isempty(missingFields), errID, ...
        'Metadata file %s is missing required field(s): %s.', ...
        channelInfo(iChan).matFile, strjoin(missingFields, ', '));

    dim_names = metaData.dim_names;
    assert(iscell(dim_names) && isequal(dim_names(:).', {'Y','X','T'}), errID, ...
        'Metadata file %s must describe dimensions as {''Y'', ''X'', ''T''}.', ...
        channelInfo(iChan).matFile);

    assert(isnumeric(metaData.datSize) && numel(metaData.datSize) == 2 && ...
        all(isfinite(metaData.datSize(:))) && all(metaData.datSize(:) > 0), errID, ...
        'Invalid datSize in metadata file: %s', channelInfo(iChan).matFile);
    assert(isnumeric(metaData.datLength) && isscalar(metaData.datLength) && ...
        isfinite(metaData.datLength) && metaData.datLength > 0, errID, ...
        'Invalid datLength in metadata file: %s', channelInfo(iChan).matFile);
    assert(isnumeric(metaData.Freq) && isscalar(metaData.Freq) && ...
        isfinite(metaData.Freq) && metaData.Freq > 0, errID, ...
        'Invalid Freq in metadata file: %s', channelInfo(iChan).matFile);

    datatype = lower(strrep(char(string(metaData.Datatype)), '*', ''));
    assert(ismember(datatype, {'single', 'float32'}), errID, ...
        'Imported .dat files must be single precision. Invalid Datatype in %s.', ...
        channelInfo(iChan).matFile);

    nY = double(metaData.datSize(1));
    nX = double(metaData.datSize(2));
    nT = double(metaData.datLength);
    bytesPerSample = 4;
    frameStride = nY * nX * bytesPerSample;

    fileInfo = dir(channelInfo(iChan).datFile);
    expectedBytes = nT * frameStride;
    assert(~isempty(fileInfo) && double(fileInfo.bytes) == expectedBytes, errID, ...
        ['File size does not match metadata for %s. Expected %.0f bytes ' ...
         'but found %.0f bytes.'], ...
        channelInfo(iChan).datFile, expectedBytes, double(fileInfo.bytes));

    if isNoOp
        fr_start = 1;
        fr_stop = 0;
        nRemoved = 0;
        newLength = nT;
    else
        % The removed interval is [segment_start_sec, segment_end_sec). Frames
        % before segment_start_sec are retained. Frames up to segment_end_sec are
        % removed using ceil to avoid retaining a partial frame at the end.
        nFramesBeforeSegment = ceil(opts.segment_start_sec * metaData.Freq);
        fr_start = nFramesBeforeSegment + 1;
        fr_stop = ceil(opts.segment_end_sec * metaData.Freq);

        assert(fr_start >= 1 && fr_start <= nT, errID, ...
            ['segment_start_sec does not select an existing start frame in %s. ' ...
             'Computed start frame: %d. Movie length: %d frame(s).'], ...
            channelInfo(iChan).datFileName, fr_start, nT);
        assert(fr_stop >= 1 && fr_stop <= nT, errID, ...
            ['segment_end_sec does not select an existing end frame in %s. ' ...
             'Computed end frame: %d. Movie length: %d frame(s).'], ...
            channelInfo(iChan).datFileName, fr_stop, nT);
        assert(fr_start <= fr_stop, errID, ...
            ['The selected segment does not remove any frame from %s. ' ...
             'Computed frame range: [%d %d].'], ...
            channelInfo(iChan).datFileName, fr_start, fr_stop);

        nRemoved = fr_stop - fr_start + 1;
        newLength = nT - nRemoved;

        assert(newLength > 0, errID, ...
            ['Invalid segment range for %s. The resulting movie would have %d ' ...
             'frame(s).'], ...
            channelInfo(iChan).datFileName, newLength);
    end

    channelInfo(iChan).metaData = metaData;
    channelInfo(iChan).nT = nT;
    channelInfo(iChan).fr_start = fr_start;
    channelInfo(iChan).fr_stop = fr_stop;
    channelInfo(iChan).nRemoved = nRemoved;
    channelInfo(iChan).newLength = newLength;
end

% Return the dynamic file manifest even when this call is a no-op.
outFile = {channelInfo.datFile};

if isNoOp
    disp('No movie segment selected for removal.')
    return
end

% events.mat has a fixed required filename and must live in SaveFolder.
% When present, it is updated automatically to keep event timestamps aligned
% with the shortened movies. If it is absent, event handling is skipped.
eventsFile = fullfile(SaveFolder, 'events.mat');
hasEventsFile = isfile(eventsFile);

%% ==========================================================
% RAM-SAFE SESSION SEGMENT REMOVAL
% ==========================================================
disp('Removing movie segment from imported movies (RAM-Safe mode)...')

% Phase 1: write all temporary files before replacing any original.
for iChan = 1:nChannels

    metaData = channelInfo(iChan).metaData;
    inFile = channelInfo(iChan).datFile;
    [inFolder, inBaseName] = fileparts(inFile);
    tmpDatFile = fullfile(inFolder, sprintf('%s_remove_segment_tmp_%s.dat', ...
        inBaseName, datestr(now, 'yyyymmdd_HHMMSSFFF')));

    channelInfo(iChan).tmpDatFile = tmpDatFile;

    assert(~strcmpi(inFile, tmpDatFile), errID, ...
        ['Temporary output file resolves to the same path as the input file. ' ...
         'Input: %s Output: %s'], inFile, tmpDatFile);

    fidIn = fopen(inFile, 'rb');
    assert(fidIn > 0, errID, 'Failed to open input file: %s', inFile);
    cIn = onCleanup(@() safeFclose(fidIn));

    fidOut = fopen(tmpDatFile, 'wb');
    assert(fidOut > 0, errID, 'Failed to create temporary output file: %s', tmpDatFile);
    cOut = onCleanup(@() safeFclose(fidOut));

    nY = double(metaData.datSize(1));
    nX = double(metaData.datSize(2));
    bytesPerSample = 4;
    frameStride = nY * nX * bytesPerSample;

    fprintf('Removing segment from %s (%s)...\n', ...
        channelInfo(iChan).datFileName, channelInfo(iChan).colorName);

    % Copy frames before the removed segment and after the removed segment.
    copyRanges = [1, channelInfo(iChan).fr_start - 1; ...
                  channelInfo(iChan).fr_stop + 1, channelInfo(iChan).nT];

    for iRange = 1:size(copyRanges, 1)

        rangeStart = copyRanges(iRange, 1);
        rangeStop = copyRanges(iRange, 2);

        if rangeStart > rangeStop
            continue
        end

        for t = rangeStart:rangeStop

            seekStatus = fseek(fidIn, (t - 1) * frameStride, 'bof');
            assert(seekStatus == 0, errID, ...
                'Failed to seek to frame %d in %s.', t, inFile);

            frame = fread(fidIn, nY*nX, '*single');
            assert(numel(frame) == nY*nX, errID, ...
                'Failed to read frame %d from %s.', t, inFile);

            nWritten = fwrite(fidOut, frame, 'single');
            assert(nWritten == nY*nX, errID, ...
                'Failed to write frame %d to %s.', t, tmpDatFile);
        end
    end

    safeFclose(fidIn);
    clear cIn

    safeFclose(fidOut);
    clear cOut

    tmpInfo = dir(tmpDatFile);
    expectedTmpBytes = channelInfo(iChan).newLength * frameStride;
    assert(~isempty(tmpInfo) && double(tmpInfo.bytes) == expectedTmpBytes, errID, ...
        ['Temporary file size does not match expected output size for %s. ' ...
         'Expected %.0f bytes but found %.0f bytes.'], ...
        channelInfo(iChan).datFileName, expectedTmpBytes, double(tmpInfo.bytes));
end

% Phase 2: replace originals only after all temporary files were created.
for iChan = 1:nChannels
    delete(channelInfo(iChan).datFile);
    movefile(channelInfo(iChan).tmpDatFile, channelInfo(iChan).datFile);
end

% Phase 3: update same-name metadata files. AcqInfos.mat is intentionally
% not modified in this branch.
for iChan = 1:nChannels
    metaData = channelInfo(iChan).metaData;
    metaData.datLength = channelInfo(iChan).newLength;
    save(channelInfo(iChan).matFile, '-struct', 'metaData');
end

% events.mat is updated once because it is shared by the session, not by
% individual channel files. If it does not exist, event handling is skipped.
if hasEventsFile
    updateEventsFileForRemovedSegment( ...
        eventsFile, ...
        opts.segment_start_sec, ...
        opts.segment_end_sec);
end

disp('Done')
end

% =========================================================================
% Local helpers
% =========================================================================

function datFileName = colorNameToDatFileName(colorName)
%COLORNAMETODATFILENAME Convert an illumination color name to a .dat name.

errID = 'umIToolbox:remove_movie_segment:InvalidChannelName';

colorName = strtrim(char(string(colorName)));
colorNameLower = lower(colorName);

switch colorNameLower
    case 'red'
        datFileName = 'red.dat';
    case {'amber', 'yellow'}
        datFileName = 'yellow.dat';
    case 'green'
        datFileName = 'green.dat';
    case 'fluo'
        datFileName = 'fluo.dat';
    otherwise
        tokens = regexp(colorName, '^Fluo\s*#\d+\s+(\d+)\s*nm$', 'tokens', 'once');
        if isempty(tokens)
            error(errID, ...
                'Unsupported illumination Color name: "%s".', colorName);
        end
        datFileName = sprintf('fluo_%s.dat', tokens{1});
end

end

function updateEventsFileForRemovedSegment(eventsFile, segmentStartSec, segmentEndSec)
%UPDATEEVENTSFILEFORREMOVEDSEGMENT Update events.mat after segment removal.
%
%   updateEventsFileForRemovedSegment(eventsFile, segmentStartSec,
%       segmentEndSec)
%
%   Updates the unique events.mat file after removing one movie segment from
%   all imported movies. The events.mat file is expected to contain:
%       eventID       - N-by-1 uint16 vector. Event index.
%       timestamps    - N-by-1 single vector. Event timestamps in seconds.
%       state         - N-by-1 logical vector. Trigger state, where
%                       true = HIGH/ON and false = LOW/OFF.
%       eventNameList - M-by-1 cell array of event names. Position in this
%                       array corresponds to eventID.
%
%   Events are removed using paired ON/OFF intervals for each eventID:
%       - Any ON/OFF pair overlapping the removed segment is deleted.
%       - Remaining events after the removed segment are shifted backward by
%         the removed duration.
%       - eventNameList is compacted when all repetitions of an eventID are
%         removed, and remaining eventID values are remapped accordingly.

errID = 'umIToolbox:remove_movie_segment:InvalidEventsFile';

evntInfo = load(eventsFile, '-mat');

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

tStart_s = single(segmentStartSec);
tStop_s = single(segmentEndSec);
dt_s = single(segmentEndSec - segmentStartSec);

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

    eventNameList = cell(0, 1);
    eventID = uint16.empty(0, 1);
    timestamps = single.empty(0, 1);
    state = logical.empty(0, 1);

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

if isempty(evntInfo.eventID)
    delete(eventsFile);

    warning('umIToolbox:remove_movie_segment:NoEventsRemaining', ...
        ['All events were removed after deleting the movie segment. ' ...
         'The original events.mat file was deleted.']);
else
    save(eventsFile, '-struct', 'evntInfo');
end

end
