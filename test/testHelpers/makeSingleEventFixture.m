function fx = makeSingleEventFixture(folder, varargin)
%MAKESINGLEEVENTFIXTURE Synthetic recording with exactly one event instance (E = 1).
%
%   fx = makeSingleEventFixture(folder)
%   fx = makeSingleEventFixture(folder, 'Variant', 'shortBaseline', ...)
%
%   Builds, in FOLDER, a small non-square recording (Y ~= X, distinct T and
%   E) whose value encodes its frame, with an events.mat holding one event
%   instance, and returns the continuous data, its expected event split and
%   the same single-event data as .dat files and UMT structs. Nothing is
%   tracked in the repository: everything is created on the fly.
%
%   Name-Value options:
%       'Variant'  - 'plain' (default): the event comes after the baseline;
%                    'shortBaseline': the event comes sooner than the baseline
%                      after the recording start, so the first trial's missing
%                      baseline frames are filled (split_data_by_event);
%                    'nearEnd': the event is near the recording end;
%                    'multiCondition': three conditions A, B, C with one
%                      repetition each; the single event instance is the one
%                      of condition B (fx.condition), reached through the
%                      condition selection of EventsManager.
%       'Rate'     - frame rate in Hz (default 10).
%       'NT'       - number of frames (default 100).
%       'Baseline' - baselinePeriod in seconds (default 2).
%
%   Output struct FX (E = 1 means the split data have size [Y X T 1], held in
%   MATLAB as [Y X T]; size(x, 1:4) recovers the E axis):
%       folder, rate, nT, baseline, variant, condition ('' unless
%       multiCondition), onsetsSec
%       data              - continuous Y-X-T single array
%       datFile           - headered Y-X-T .dat of DATA (in.dat)
%       expectedSplit     - Y-X-T-E expected output of the split with E = 1
%                           (3-D in memory); for multiCondition, the split of
%                           the selected condition
%       expectedSplitAll  - same for every instance of events.mat (equals
%                           expectedSplit unless multiCondition)
%       splitDatFile      - headered Y-X-T-E .dat (E = 1) of expectedSplit
%                           (byEv.dat)
%       imageUMT          - image UMT, Y-X-T, with meta.FrameRateHz
%       imageUMTSplit     - image UMT, Y-X-T-E with E = 1, with eventInfo
%       roiUMT            - ROI UMT, ROI-T (3 ROIs)
%       roiUMTSplit       - ROI UMT, ROI-T-E with E = 1, with eventInfo
%
%   See also writeTestDat, split_data_by_event.

p = inputParser;
addRequired(p, 'folder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Variant', 'plain', @(x) any(strcmp(x, ...
    {'plain', 'shortBaseline', 'nearEnd', 'multiCondition'})));
addParameter(p, 'Rate', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NT', 100, @(x) isnumeric(x) && isscalar(x) && x >= 20);
addParameter(p, 'Baseline', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, folder, varargin{:});

folder = char(folder);
rate = p.Results.Rate;
nT = p.Results.NT;
baseline = p.Results.Baseline;
variant = p.Results.Variant;
nY = 6;
nX = 5;
evOnset = round(baseline * rate);          % baseline columns of a trial

% ---- Recording: the value of a pixel encodes frame, row and column ----------
data = repmat(reshape(single(1:nT), 1, 1, nT), nY, nX, 1) + ...
    repmat(single(0.01 * (1:nY)'), 1, nX, nT) + ...
    repmat(single(0.0001 * (1:nX)), nY, 1, nT);

% ---- Events and the frames of the selected trial ---------------------------
condition = '';
eventNameList = {'A'};
switch variant
    case 'plain'
        onsets = 3;
    case 'shortBaseline'
        onsets = 1;
    case 'nearEnd'
        onsets = (nT - 5) / rate;
    case 'multiCondition'
        onsets = [3 6 9];
        eventNameList = {'A', 'B', 'C'};
        condition = 'B';
end
eventIDs = (1:numel(onsets)).';
if ~strcmp(variant, 'multiCondition')
    eventIDs = ones(numel(onsets), 1);
end
iWriteEvents(folder, onsets, eventIDs, eventNameList, baseline);

frOn = ceil(onsets(:) * rate);
frStart = frOn - evOnset + 1;
if numel(frOn) >= 2
    trialLen = [diff(frOn); max(diff(frOn))];
else
    trialLen = nT - frStart + 1;
end
cropLen = min(trialLen);                   % all-instance split crops to this
selected = 1;
if ~isempty(condition)
    selected = find(strcmp(eventNameList, condition), 1);
end
expectedSplit = iTrial(data, frStart(selected), trialLen(selected), evOnset);
allTrials = arrayfun(@(k) iTrial(data, frStart(k), cropLen, evOnset), ...
    1:numel(frOn), 'UniformOutput', false);
expectedSplitAll = cat(4, allTrials{:});       % 3-D in memory when there is one

% ---- Files and UMT structs -------------------------------------------------
datFile = fullfile(folder, 'in.dat');
writeTestDat(datFile, data, rate);
splitDatFile = fullfile(folder, 'byEv.dat');
writeTestDat(splitDatFile, expectedSplit, rate, 'DimNames', {'Y', 'X', 'T', 'E'});

meta = struct('FrameRateHz', rate);
imageUMT = genUMTStruct(data, 'kind', 'image', 'entryName', 'main', ...
    'dimNames', {'Y', 'X', 'T'}, 'meta', meta);
imageUMTSplit = genUMTStruct(expectedSplit, 'kind', 'image', 'entryName', 'main', ...
    'dimNames', {'Y', 'X', 'T', 'E'}, 'meta', meta);

roiLabels = struct('ROI', {{'r1', 'r2', 'r3'}});
roiTraces = reshape(data(1:3, 1, :), 3, nT);
roiUMT = genUMTStruct(roiTraces, 'kind', 'roi', 'entryName', 'traces', ...
    'dimNames', {'ROI', 'T'}, 'labels', roiLabels, 'meta', meta);
roiSplit = reshape(expectedSplit(1:3, 1, :), 3, size(expectedSplit, 3));
roiUMTSplit = genUMTStruct(roiSplit, 'kind', 'roi', 'entryName', 'traces', ...
    'dimNames', {'ROI', 'T', 'E'}, 'labels', roiLabels, 'meta', meta);

% eventInfo of the single selected instance (the instance of the condition).
ev = EventsManager(folder);
allInfo = ev.exportEventInfo('FrameRateHz', rate, 'IncludeIgnored', true);
oneInfo = iSelectRow(allInfo, selected);
imageUMTSplit = appendUMTEventInfo(imageUMTSplit, 'eventInfo', ...
    iAsUMTEventInfo(oneInfo, eventNameList), 'overwrite', true);
roiUMTSplit = appendUMTEventInfo(roiUMTSplit, 'eventInfo', ...
    iAsUMTEventInfo(oneInfo, eventNameList), 'overwrite', true);

fx = struct('folder', folder, 'rate', rate, 'nT', nT, 'baseline', baseline, ...
    'variant', variant, 'condition', condition, 'onsetsSec', onsets(:), ...
    'data', data, 'datFile', datFile, 'expectedSplit', expectedSplit, ...
    'expectedSplitAll', expectedSplitAll, 'splitDatFile', splitDatFile, ...
    'imageUMT', imageUMT, 'imageUMTSplit', imageUMTSplit, ...
    'roiUMT', roiUMT, 'roiUMTSplit', roiUMTSplit);
end

function trial = iTrial(data, start, len, evOnset)
%ITRIAL One trial (Y-X-T) from frame START, LEN columns; missing leading frames
% (start < 1) take the mean of the baseline frames the trial has.
cols = start:start + len - 1;
valid = cols >= 1 & cols <= size(data, 3);
trial = nan(size(data, 1), size(data, 2), numel(cols), 'single');
trial(:, :, valid) = data(:, :, cols(valid));
nLead = find(valid, 1, 'first') - 1;
if nLead > 0
    trial(:, :, 1:nLead) = repmat(mean(trial(:, :, nLead + 1:evOnset), 3), 1, 1, nLead);
end
end

function row = iSelectRow(info, k)
%ISELECTROW The per-instance fields of struct INFO restricted to instance K.
row = info;
for f = {'eventID', 'repetitionIndex', 'selected', 'selectedEvents', 'durationSec', 'onsetSec'}
    row.(f{1}) = info.(f{1})(k);
end
end

function info = iAsUMTEventInfo(row, eventNameList)
%IASUMTEVENTINFO UMT eventInfo (instances) from an exportEventInfo row.
info = struct('eventID', double(row.eventID(:)), ...
    'repetitionIndex', double(row.repetitionIndex(:)), ...
    'eventName', string(eventNameList(double(row.eventID(:)))).', ...
    'eventAxisMode', 'instances', 'selected', logical(row.selected(:)), ...
    'durationSec', row.durationSec(:), 'onsetSec', row.onsetSec(:), ...
    'baselinePeriod', double(row.baselinePeriod));
info.eventName = info.eventName(:);
end

function iWriteEvents(folder, onsets, eventIDs, eventNameList, baseline)
%IWRITEEVENTS Minimal events.mat (external-signal events, 0.5 s each).
n = numel(onsets);
timestamps = reshape([onsets(:), onsets(:) + 0.5]', [], 1);
state = repmat([true; false], n, 1);
eventID = reshape(repelem(uint16(eventIDs(:)), 2), [], 1);
selectedEvents = true(2 * n, 1);
baselinePeriod = baseline;
b_hasExternalSignal = true;
save(fullfile(folder, 'events.mat'), 'timestamps', 'state', 'eventID', ...
    'eventNameList', 'selectedEvents', 'baselinePeriod', 'b_hasExternalSignal', '-mat');
end
