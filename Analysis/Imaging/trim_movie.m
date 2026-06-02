function outFile = trim_movie(SaveFolder, varargin)
%TRIM_MOVIE Trim all imported disk-backed .dat movies in a session folder.
%
%   outFile = trim_movie(SaveFolder)
%   outFile = trim_movie(SaveFolder, opts)
%
%   Trims all imported .dat movies listed in AcqInfos.mat. Channel names are
%   read from AcqInfoStream.Illumination<N>.Color and translated to imported
%   .dat filenames. Each movie is trimmed in RAM-safe mode by writing a
%   temporary copy, deleting the original file, and renaming the temporary
%   copy to the original filename.
%
%   Inputs:
%       SaveFolder - Folder containing AcqInfos.mat, imported .dat files,
%                    same-name channel metadata .mat files, and optionally
%                    events.mat.
%       opts       - Optional structure with fields:
%                    crop_start_sec - Seconds to remove from the beginning.
%                                     Default: 0
%                    crop_end_sec   - Seconds to remove from the end.
%                                     Default: 0
%
%   Output:
%       outFile - Cell array containing the full paths of all imported .dat
%                 files discovered from AcqInfos.mat. This manifest is
%                 returned after validation even when crop_start_sec and
%                 crop_end_sec are both zero and no files are rewritten.
%
%   Notes:
%       - This function does not modify AcqInfos.mat. In this metadata model,
%         AcqInfos.mat is used only to discover imported channels.
%       - The same-name metadata file associated with each .dat file is
%         updated after trimming, for example red.dat -> red.mat.
%       - crop_start_sec and crop_end_sec are converted to frames separately
%         for each channel using ceil(seconds * metaData.Freq).
%       - The function validates all channel files, metadata, and crop ranges
%         before replacing any original .dat file.
%       - If crop_start_sec and crop_end_sec are both zero, the function
%         performs no file rewrite and returns the full file manifest.
%       - If events.mat exists, it is updated automatically using seconds.
%         If imported channel durations are not exactly equal, the shortest
%         duration is used as the event recording duration.
%       - Only complete ON/OFF event pairs fully contained inside the kept
%         interval are retained.
%       - Retained event timestamps are shifted so the trimmed movie starts
%         at time zero.
%       - eventNameList is compacted when all repetitions for an eventID are
%         removed, and remaining eventID values are remapped accordingly.
%
%   Channel filename mapping:
%       Red                         -> red.dat
%       Amber or Yellow             -> yellow.dat
%       Green                       -> green.dat
%       Fluo                        -> fluo.dat
%       Fluo #<N> <wavelength> nm   -> fluo_<wavelength>.dat

% Defaults
default_Output = {'fluo_475.dat', 'fluo_567.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; %#ok<NASGU> % Reference for PipelineManager. Actual outputs are stored in outFile.
default_opts = struct( ...
    'crop_start_sec', 0, ...
    'crop_end_sec', 0);
opts_values = struct( ... %#ok<NASGU>
    'crop_start_sec', [0 Inf], ...
    'crop_end_sec', [0 Inf]);

%%% Arguments parsing and validation %%%
p = inputParser;
p.FunctionName = 'trim_movie';
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
errID = 'umIToolbox:trim_movie:InvalidInput';

assert(isnumeric(opts.crop_start_sec) && isscalar(opts.crop_start_sec) && ...
    isfinite(opts.crop_start_sec) && opts.crop_start_sec >= 0, errID, ...
    'crop_start_sec must be a finite numeric scalar >= 0.');
assert(isnumeric(opts.crop_end_sec) && isscalar(opts.crop_end_sec) && ...
    isfinite(opts.crop_end_sec) && opts.crop_end_sec >= 0, errID, ...
    'crop_end_sec must be a finite numeric scalar >= 0.');
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
    'duration_s', [], ...
    'fr_start', [], ...
    'fr_stop', [], ...
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
    assert(all(ismember(dim_names, {'Y','X','T'})), errID, ...
        'Metadata file %s must describe a 3D matrix with dimensions Y, X, and T.', ...
        channelInfo(iChan).matFile);

    idxT = find(strcmp('T', dim_names));
    sz = [metaData.datSize, metaData.datLength];
    nT = sz(idxT);

    assert(isnumeric(nT) && isscalar(nT) && isfinite(nT) && nT > 0, errID, ...
        'Invalid datLength in metadata file: %s', channelInfo(iChan).matFile);
    assert(isnumeric(metaData.Freq) && isscalar(metaData.Freq) && ...
        isfinite(metaData.Freq) && metaData.Freq > 0, errID, ...
        'Invalid Freq in metadata file: %s', channelInfo(iChan).matFile);

    nRemoveStart = ceil(opts.crop_start_sec * metaData.Freq);
    nRemoveEnd = ceil(opts.crop_end_sec * metaData.Freq);

    % Assert that requested removed frames exist. Removing zero frames is valid.
    assert(nRemoveStart <= nT, errID, ...
        ['crop_start_sec removes %d frame(s) from %s, but the movie only ' ...
         'has %d frame(s).'], ...
        nRemoveStart, channelInfo(iChan).datFileName, nT);
    assert(nRemoveEnd <= nT, errID, ...
        ['crop_end_sec removes %d frame(s) from %s, but the movie only ' ...
         'has %d frame(s).'], ...
        nRemoveEnd, channelInfo(iChan).datFileName, nT);

    fr_start = nRemoveStart + 1;
    fr_stop = nT - nRemoveEnd;
    newLength = fr_stop - fr_start + 1;

    assert(newLength > 0, errID, ...
        ['Invalid crop range for %s. The resulting movie would have %d ' ...
         'frame(s).'], ...
        channelInfo(iChan).datFileName, newLength);

    channelInfo(iChan).metaData = metaData;
    channelInfo(iChan).nT = nT;
    channelInfo(iChan).duration_s = nT / metaData.Freq;
    channelInfo(iChan).fr_start = fr_start;
    channelInfo(iChan).fr_stop = fr_stop;
    channelInfo(iChan).newLength = newLength;
end

% Return the dynamic file manifest even when this call is a no-op.
outFile = {channelInfo.datFile};

if opts.crop_start_sec == 0 && opts.crop_end_sec == 0
    disp('No movie trimming requested.')
    return
end

% events.mat has a fixed filename and lives in SaveFolder when present.
eventsFile = fullfile(SaveFolder, 'events.mat');
bUpdateEvents = isfile(eventsFile);

% Use exact equality as requested. If durations differ, use the shortest one
% when updating events.mat.
recordingDurations_s = [channelInfo.duration_s];
recordingDuration_s = min(recordingDurations_s);
if bUpdateEvents && numel(unique(recordingDurations_s)) > 1
    warning('umIToolbox:trim_movie:DifferentChannelDurations', ...
        ['Imported channel durations are not exactly equal. events.mat ' ...
         'will be cropped using the shortest duration: %.9g seconds.'], ...
        recordingDuration_s);
end

assert(opts.crop_start_sec <= recordingDuration_s, errID, ...
    ['crop_start_sec is %.9g seconds, but the shortest imported recording ' ...
     'duration is %.9g seconds.'], ...
    opts.crop_start_sec, recordingDuration_s);
assert(opts.crop_end_sec <= recordingDuration_s, errID, ...
    ['crop_end_sec is %.9g seconds, but the shortest imported recording ' ...
     'duration is %.9g seconds.'], ...
    opts.crop_end_sec, recordingDuration_s);
assert((recordingDuration_s - opts.crop_start_sec - opts.crop_end_sec) > 0, errID, ...
    ['Invalid crop interval. The resulting event time span would be <= 0 ' ...
     'seconds.']);

%% ==========================================================
% RAM-SAFE SESSION TRIMMING
% ==========================================================
disp('Cropping imported movies (RAM-Safe mode)...')

% Phase 1: write all temporary files before replacing any original.
for iChan = 1:nChannels

    metaData = channelInfo(iChan).metaData;
    inFile = channelInfo(iChan).datFile;
    [inFolder, inBaseName] = fileparts(inFile);
    tmpDatFile = fullfile(inFolder, sprintf('%s_trim_movie_tmp_%s.dat', ...
        inBaseName, datestr(now, 'yyyymmdd_HHMMSSFFF')));

    channelInfo(iChan).tmpDatFile = tmpDatFile;

    fidIn = fopen(inFile, 'rb');
    assert(fidIn > 0, errID, 'Failed to open input file: %s', inFile);
    cIn = onCleanup(@() safeFclose(fidIn));

    fidOut = fopen(tmpDatFile, 'wb');
    assert(fidOut > 0, errID, 'Failed to create temporary output file: %s', tmpDatFile);
    cOut = onCleanup(@() safeFclose(fidOut));

    nY = metaData.datSize(1);
    nX = metaData.datSize(2);
    bytesPerSample = 4;
    frameStride = nY * nX * bytesPerSample;

    seekStatus = fseek(fidIn, (channelInfo(iChan).fr_start - 1) * frameStride, 'bof');
    assert(seekStatus == 0, errID, ...
        'Failed to seek to frame %d in %s.', channelInfo(iChan).fr_start, inFile);

    fprintf('Cropping %s (%s)...\n', ...
        channelInfo(iChan).datFileName, channelInfo(iChan).colorName);

    for t = channelInfo(iChan).fr_start:channelInfo(iChan).fr_stop
        frame = fread(fidIn, nY*nX, ['*' metaData.Datatype]);
        assert(numel(frame) == nY*nX, errID, ...
            'Failed to read frame %d from %s.', t, inFile);

        nWritten = fwrite(fidOut, frame, metaData.Datatype);
        assert(nWritten == nY*nX, errID, ...
            'Failed to write frame %d to %s.', t, tmpDatFile);
    end

    safeFclose(fidIn);
    clear cIn

    safeFclose(fidOut);
    clear cOut
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

% Update events.mat once when present because it is shared by the session,
% not by individual channel files.
if bUpdateEvents
    updateEventsFileForTrimmedMovie( ...
        eventsFile, ...
        opts.crop_start_sec, ...
        opts.crop_end_sec, ...
        recordingDuration_s);
end

disp('Done')
end

% =========================================================================
% Local helpers
% =========================================================================

function datFileName = colorNameToDatFileName(colorName)
%COLORNAMETODATFILENAME Convert an illumination color name to a .dat name.

errID = 'umIToolbox:trim_movie:InvalidChannelName';

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

function updateEventsFileForTrimmedMovie(eventsFile, cropStartSec, cropEndSec, recordingDurationSec)
%UPDATEEVENTSFILEFORTRIMMEDMOVIE Update events.mat after session trimming.
%
%   updateEventsFileForTrimmedMovie(eventsFile, cropStartSec, cropEndSec,
%       recordingDurationSec)
%
%   Updates the unique events.mat file after trimming the beginning and/or
%   end of all imported movies. The events.mat file is expected to contain:
%       eventID       - N-by-1 uint16 vector. Event index.
%       timestamps    - N-by-1 single vector. Event timestamps in seconds.
%       state         - N-by-1 logical vector. Trigger state, where
%                       true = HIGH/ON and false = LOW/OFF.
%       eventNameList - M-by-1 cell array of event names. Position in this
%                       array corresponds to eventID.
%
%   Only complete ON/OFF pairs fully contained inside the retained time
%   interval are kept. Remaining timestamps are shifted so that the trimmed
%   movie starts at time zero. eventNameList is compacted when all
%   repetitions of an eventID are removed.

errID = 'umIToolbox:trim_movie:InvalidEventsFile';

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

% Kept interval. The upper bound is treated as inclusive for event edges at
% the end of the retained recording interval.
tStart_s = single(cropStartSec);
tStop_s = single(recordingDurationSec - cropEndSec);

removeEvent = false(nEvents, 1);

% Keep only complete ON/OFF pairs fully contained in the retained interval.
uniqueEventIDs = unique(eventID, 'stable');

for iID = 1:numel(uniqueEventIDs)

    thisID = uniqueEventIDs(iID);
    idxThis = find(eventID == thisID);

    iEdge = 1;

    while iEdge <= numel(idxThis)

        idxOn = idxThis(iEdge);

        if ~state(idxOn)
            removeEvent(idxOn) = true;
            warning('umIToolbox:trim_movie:MalformedEventPair', ...
                ['Removed an unpaired OFF edge for eventID %d while ' ...
                 'updating events.mat.'], thisID);
            iEdge = iEdge + 1;
            continue
        end

        if iEdge == numel(idxThis)
            removeEvent(idxOn) = true;
            warning('umIToolbox:trim_movie:MalformedEventPair', ...
                ['Removed an unpaired ON edge for eventID %d while ' ...
                 'updating events.mat.'], thisID);
            iEdge = iEdge + 1;
            continue
        end

        idxOff = idxThis(iEdge + 1);

        if state(idxOff)
            removeEvent(idxOn) = true;
            warning('umIToolbox:trim_movie:MalformedEventPair', ...
                ['Removed an ON edge for eventID %d because the next edge ' ...
                 'was not an OFF edge.'], thisID);
            iEdge = iEdge + 1;
            continue
        end

        offTime = timestamps(idxOff);
        onTime = timestamps(idxOn);

        pairFullyInsideKeptInterval = onTime >= tStart_s && offTime <= tStop_s;

        if ~pairFullyInsideKeptInterval
            removeEvent([idxOn, idxOff]) = true;
        end

        iEdge = iEdge + 2;
    end
end

% Also remove any single edge outside the kept interval. This is redundant
% for well-formed pairs but protects against malformed data.
outsideKept = timestamps < tStart_s | timestamps > tStop_s;
removeEvent = removeEvent | outsideKept;

keepEvent = ~removeEvent;

eventID = eventID(keepEvent);
timestamps = timestamps(keepEvent) - tStart_s;
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

    warning('umIToolbox:trim_movie:NoEventsRemaining', ...
        ['All events were removed after trimming the movie. ' ...
         'The original events.mat file was deleted.']);
else
    save(eventsFile, '-struct', 'evntInfo');
end

end
