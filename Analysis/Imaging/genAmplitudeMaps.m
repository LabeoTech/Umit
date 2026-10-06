function [outData, metaData] = genAmplitudeMaps(data, SaveFolder, varargin)
%GENAMPLITUDEMAPS Compute event-wise response amplitude maps from imaging data.
%
%   outData = genAmplitudeMaps(data, SaveFolder)
%   [outData, metaData] = genAmplitudeMaps(data, SaveFolder, 'BaselineMeasure', value, ...)
%
%   This function computes response amplitude maps by subtracting an
%   aggregate baseline value from an aggregate response value along the time
%   dimension for each event condition.
%
%   Supported inputs:
%       1) Numeric image time series with dimensions YXT
%       2) Raw .dat filename containing continuous YXT data
%       3) UMT structure
%       4) .umt filename
%
%   Event handling:
%       - For continuous non-UMT inputs, an "events.mat" file must be
%         available in SaveFolder.
%       - For event-split UMT image inputs (YXTE), top-level UMT eventInfo
%         is used directly.
%
%   Name-Value parameters:
%       'BaselineMeasure' - Aggregate function applied to baseline frames:
%                           'mean' | 'median' | 'min' | 'max'
%                           Default: 'median'
%       'ResponseMeasure' - Aggregate function applied to response frames:
%                           'mean' | 'median' | 'min' | 'max'
%                           Default: 'max'
%       'TimeWindow_sec'  - Response window in seconds relative to stimulus
%                           onset, where 0 is the first frame after the
%                           baseline period. Use:
%                               'all'
%                               [startSec endSec]
%                           Default: 'all'
%       'FrameRateHz'     - Frame rate of DATA (Hz). PipelineManager injects
%                           it from the data; a .dat input's header provides
%                           it otherwise, and a UMT entry's meta.FrameRateHz.
%                           In-RAM arrays need it explicitly; AcqInfos.mat is
%                           not used.
%
%   Output:
%       - Continuous inputs (numeric YXT, .dat, UMT YXT entry): numeric
%         Y x X x E array, one amplitude map per event condition, saved as
%         .dat by PipelineManager; metaData holds dimNames {'Y','X','E'}.
%         The .dat stores no labels: resolveDatEventMapping matches its E
%         axis to events.mat (one slice per condition).
%       - Event-split UMT inputs: UMT struct with the aggregated eventInfo
%         (selected, durationSec, nInstances); metaData is an empty struct.
%
%   Notes:
%       - Conditions follow EventsManager.conditionAggregationPlan (.dat
%         header Phase 8c): ignored instances are excluded, conditions are
%         ordered by first appearance, and a condition whose instances are
%         all ignored gives a NaN map.
%       - BaselineMeasure and ResponseMeasure are applied jointly across
%         frames AND trials within each condition, not per-trial-then-
%         averaged-across-trials. With the default ResponseMeasure='max',
%         this means the single largest response-window value across all
%         trials of a condition, not the mean of each trial's own peak.
%       - Raw .dat input is processed in spatial X slabs. Slab width is
%         derived from the baseline/response frame counts and the largest
%         condition repetition count, so both trial buffers stay within the
%         calculated scratch-memory budget. The final Y x X x E amplitude
%         map remains resident in RAM.

% Legacy pipeline placeholder
default_Output = 'amplitudeMap.umt';
metaData = struct();

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'genAmplitudeMaps';
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || isstring(x));
addParameter(p, 'BaselineMeasure', 'median');
addParameter(p, 'ResponseMeasure', 'max');
addParameter(p, 'TimeWindow_sec', 'all');
addParameter(p, 'FrameRateHz', []);
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
assert(isfolder(SaveFolder), ...
    'Umitoolbox:genAmplitudeMaps:invalidSaveFolder', ...
    'SaveFolder does not exist: "%s".', SaveFolder);

baselineMeasure = char(string(p.Results.BaselineMeasure));
responseMeasure = char(string(p.Results.ResponseMeasure));
timeWindowSec   = p.Results.TimeWindow_sec;

assert(iValidateAggName(baselineMeasure), ...
    'Umitoolbox:genAmplitudeMaps:invalidBaselineMeasure', ...
    'BaselineMeasure must be one of: mean, median, min, max.');
assert(iValidateAggName(responseMeasure), ...
    'Umitoolbox:genAmplitudeMaps:invalidResponseMeasure', ...
    'ResponseMeasure must be one of: mean, median, min, max.');
assert(iValidateTimeWindowInput(timeWindowSec), ...
    'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
    'TimeWindow_sec must be "all" or a numeric [start end] vector.');

src = iResolveInput(data, SaveFolder);

% Frame rate of the data itself: the explicit FrameRateHz (injected by
% PipelineManager), else the .dat header or the UMT entry's
% meta.FrameRateHz. AcqInfos.mat is not used (resolveDataInfoValue).
rateData = [];
if src.isRawDat
    rateData = src.fileName;
end
ownRate = [];
if isstruct(src.entry) && isfield(src.entry, 'meta') && isstruct(src.entry.meta) && ...
        isfield(src.entry.meta, 'FrameRateHz')
    ownRate = src.entry.meta.FrameRateHz;
end
frameRateHz = resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, rateData, ...
    mfilename, 'OwnValue', ownRate, 'OwnSource', 'the UMT entry meta.FrameRateHz');

if src.isRawDat
    outData = iRunChunkedDat(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder, frameRateHz);
else
    outData = iRunStandard(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder, frameRateHz);
end
if isnumeric(outData)
    % One map per condition: saved as .dat (Phase 8c).
    metaData = struct('dimNames', {{'Y','X','E'}});
end

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            mfilename, ...
            ['Compute response amplitude maps by subtracting an ' ...
             'aggregate baseline from an aggregate response period.']);

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'Input imaging data. Supports numeric YXT, raw .dat, UMT struct, and .umt.', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either', ...
            'position', 1, ...
            'callType', 'positional');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing events.mat.', ...
            'isData', false, ...
            'position', 2, ...
            'callType', 'positional');

        info = PipelineManager.addInput( ...
            info, ...
            'BaselineMeasure', ...
            'parameter', ...
            'Aggregate function for baseline frames.', ...
            'kind', 'parameter', ...
            'default', 'median', ...
            'allowed', {'mean','median','min','max'}, ...
            'position', 3, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'ResponseMeasure', ...
            'parameter', ...
            'Aggregate function for response frames.', ...
            'kind', 'parameter', ...
            'default', 'max', ...
            'allowed', {'mean','median','min','max'}, ...
            'position', 4, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'TimeWindow_sec', ...
            'parameter', ...
            'Response time window in seconds relative to stimulus onset.', ...
            'kind', 'parameter', ...
            'default', 'all', ...
            'allowed', {'all',[0 Inf]}, ...
            'position', 5, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            ['YXE amplitude maps, one per event condition: .dat for continuous ' ...
             'inputs, UMT for event-split UMT inputs.'], ...
            default_Output, ...
            1, ...
            'isData', true);

        info = PipelineManager.addOutput( ...
            info, ...
            'metaData', ...
            'metaData', ...
            'data', ...
            'Axes of the .dat output.', ...
            '', ...
            2, ...
            'isData', false);
    end
end

function outData = iRunStandard(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder, frameRateHz)

switch src.representation
    case 'continuous'
        if isfield(src, 'data') && isnumeric(src.data) && ~isempty(src.data)
            dataYXT = single(src.data);
        elseif isfield(src, 'entry') && isstruct(src.entry) && ...
                isfield(src.entry, 'value') && isnumeric(src.entry.value)
            dataYXT = single(src.entry.value);
        elseif isfield(src, 'isRawDat') && src.isRawDat
            dataYXT = loadData(src.fileName);
        else
            error('Umitoolbox:genAmplitudeMaps:missingContinuousData', ...
                'Continuous input data could not be resolved.');
        end

        assert(isfile(fullfile(SaveFolder, 'events.mat')), ...
            'Umitoolbox:genAmplitudeMaps:missingEventsFile', ...
            'The file "events.mat" was not found in SaveFolder.');

        ev = EventsManager(SaveFolder);
        % Every instance is split; the plan excludes ignored ones (8c).
        dataYXTE = single(ev.splitDataByEvents(dataYXT, 'FrameRateHz', frameRateHz, ...
            'IncludeIgnored', true));
        plan = EventsManager.conditionAggregationPlan( ...
            ev.exportEventInfo('FrameRateHz', frameRateHz, 'IncludeIgnored', true));

        baselineFrames = 1:round(double(ev.baselinePeriod) * frameRateHz);
        [baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
            size(dataYXTE, 3), baselineFrames, frameRateHz, timeWindowSec);

        outData = iConditionsToE(EventsManager.reduceByCondition(dataYXTE, plan, ...
            iAmplitudeReducer(baselineFrames, responseFrames, baselineMeasure, responseMeasure), 4));
        return

    case 'eventsplit'
        entry = src.entry;
        dataYXTE = single(entry.value);

        assert(isfield(src.UMT, 'eventInfo') && ~isempty(fieldnames(src.UMT.eventInfo)), ...
            'Umitoolbox:genAmplitudeMaps:missingUMTEventInfo', ...
            'Event-split UMT input must contain top-level eventInfo.');

        eventInfo = src.UMT.eventInfo;
        requiredFields = {'eventID','repetitionIndex','eventName','eventAxisMode'};
        assert(all(isfield(eventInfo, requiredFields)), ...
            'Umitoolbox:genAmplitudeMaps:invalidUMTEventInfo', ...
            'UMT eventInfo is missing required fields.');

        if isfield(eventInfo, 'baselinePeriod') && ~isempty(eventInfo.baselinePeriod)
            baselineFrames = 1:round(double(eventInfo.baselinePeriod) * frameRateHz);
        else
            baselineFrames = 1:min(7, size(dataYXTE, 3));
        end

        [baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
            size(dataYXTE, 3), baselineFrames, frameRateHz, timeWindowSec);

        axisMode = char(string(eventInfo.eventAxisMode));

        if strcmpi(axisMode, 'aggregated_repetitions')
            % Already one slice per condition: one map per row.
            eventIDs = double(eventInfo.eventID(:));
            ampMap = iComputeAmplitudeFromYXTE( ...
                dataYXTE, baselineFrames, responseFrames, ...
                baselineMeasure, responseMeasure, (1:numel(eventIDs)).');
            eventInfoOut = eventInfo;
        elseif strcmpi(axisMode, 'instances')
            plan = EventsManager.conditionAggregationPlan(eventInfo);
            ampMap = iConditionsToE(EventsManager.reduceByCondition(dataYXTE, plan, ...
                iAmplitudeReducer(baselineFrames, responseFrames, baselineMeasure, responseMeasure), 4));
            eventInfoOut = plan.eventInfoOut;
        else
            error('Umitoolbox:genAmplitudeMaps:unsupportedEventAxisMode', ...
                'Unsupported eventAxisMode "%s".', axisMode);
        end

    otherwise
        error('Umitoolbox:genAmplitudeMaps:unknownRepresentation', ...
            'Unknown input representation "%s".', src.representation);
end

outData = iBuildOutputUMT( ...
    ampMap, eventInfoOut, baselineMeasure, responseMeasure, timeWindowSec);

end

function outData = iRunChunkedDat(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder, frameRateHz)
%IRUNCHUNKEDDAT Amplitude maps of a YXT .dat, read in X slabs (Phase 8c plan).
Info = src.Info;
Ny = datAxisSize(Info, 'Y');
Nx = datAxisSize(Info, 'X');
Nt = datAxisSize(Info, 'T');

ev = EventsManager(SaveFolder);
[frMat, conditionIDlist] = ev.getFrameMatrix(Nt, 'FrameRateHz', frameRateHz, 'IncludeIgnored', true);
if isempty(frMat)
    error('Umitoolbox:genAmplitudeMaps:noFrames', ...
        'No event frames were returned by EventsManager.getFrameMatrix.');
end
evInfo = ev.exportEventInfo('FrameRateHz', frameRateHz, 'IncludeIgnored', true);
assert(isequal(double(evInfo.eventID(:)), double(conditionIDlist(:))), ...
    'Umitoolbox:genAmplitudeMaps:eventAxisMismatch', ...
    'The event list does not match the split trials.');
plan = EventsManager.conditionAggregationPlan(evInfo);

baselineFrames = 1:round(double(ev.baselinePeriod) * frameRateHz);
[baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
    size(frMat, 2), baselineFrames, frameRateHz, timeWindowSec);

% Only the baseline and response columns of each trial are read.
usedCols = [baselineFrames, responseFrames];
frUsed = frMat(:, usedCols);
frUsed(frUsed < 1 | frUsed > Nt) = NaN;
needed = unique(frUsed(isfinite(frUsed)));
nB = numel(baselineFrames);
reducer = iAmplitudeReducer(1:nB, nB + (1:numel(responseFrames)), ...
    baselineMeasure, responseMeasure);

nInst = size(frMat, 1);
nCond = numel(plan.conditionID);
outData = zeros(Ny, Nx, nCond, 'single');

slabIn = spatialSlabIO('open', src.fileName, 'Info', Info);
cleanupFid = onCleanup(@() spatialSlabIO('close', slabIn));

scratchBytes = double(Ny) * double(Nx) * ...
    (double(numel(usedCols)) * double(nInst) + double(numel(needed))) * ...
    double(getByteSize('single'));
nChunks = max(1, calculateMaxChunkSize(scratchBytes, 1, 0.2));
chunkX = ceil(Nx / nChunks);

for xStart = 1:chunkX:Nx
    xIdx = xStart:min(xStart + chunkX - 1, Nx);
    frames = single(spatialSlabIO('read', slabIn, xIdx, needed(:).'));
    slabE = nan(Ny, numel(xIdx), numel(usedCols), nInst, 'single');
    for iInst = 1:nInst
        valid = isfinite(frUsed(iInst, :));
        [~, loc] = ismember(frUsed(iInst, valid), needed);
        slabE(:, :, valid, iInst) = frames(:, :, loc);
    end
    outData(:, xIdx, :) = iConditionsToE(EventsManager.reduceByCondition(slabE, plan, reducer, 4));
end
end

function ampMap = iConditionsToE(perCondition)
%ICONDITIONSTOE Y x X x 1 x E (one map per condition along dim 4) to Y x X x E.
ampMap = reshape(perCondition, size(perCondition, 1), size(perCondition, 2), []);
end

function fcn = iAmplitudeReducer(baselineFrames, responseFrames, baselineMeasure, responseMeasure)
%IAMPLITUDEREDUCER Per-condition amplitude: response minus baseline, pooled
%over frames and the condition's selected trials (Y x X x T x E -> Y x X).
fcn = @(x) iAmplitudeOfTrials(x, baselineFrames, responseFrames, baselineMeasure, responseMeasure);
end

function ampMap = iAmplitudeOfTrials(x, baselineFrames, responseFrames, baselineMeasure, responseMeasure)
baselineVals = reshape(x(:, :, baselineFrames, :), size(x, 1), size(x, 2), []);
responseVals = reshape(x(:, :, responseFrames, :), size(x, 1), size(x, 2), []);
ampMap = iApplyAggFcnND(responseVals, responseMeasure) - iApplyAggFcnND(baselineVals, baselineMeasure);
end

function src = iResolveInput(dataIn, SaveFolder)
src = struct();
src.UMT = [];
src.entry = [];
src.fileName = '';
src.Info = [];
src.isRawDat = false;

if isnumeric(dataIn)
    validateattributes(dataIn, {'numeric'}, {'nonempty'}, 'genAmplitudeMaps', 'data');
    assert(ndims(dataIn) == 3, ...
        'Umitoolbox:genAmplitudeMaps:invalidNumericInput', ...
        'Numeric input must be a YXT array.');
    src.representation = 'continuous';
    src.data = single(dataIn);
    return
end

if ischar(dataIn) || (isstring(dataIn) && isscalar(dataIn))
    fileName = char(string(dataIn));
    if ~isfile(fileName)
        altPath = fullfile(SaveFolder, fileName);
        if isfile(altPath)
            fileName = altPath;
        end
    end
    [~,~,ext] = fileparts(fileName);
    ext = lower(ext);
    if strcmp(ext, '.dat')
        src.representation = 'continuous';
        src.fileName = fileName;
        src.Info = loadMetaData(fileName);
        assertDatLayout(src.Info, {{'Y','X','T'}}, 'genAmplitudeMaps');
        src.isRawDat = true;
        src.data = [];
        return
    elseif strcmp(ext, '.umt')
        loadedData = loadData(fileName);
        src.UMT = loadedData;
        [entry, rep] = iSelectUMTImageEntry(loadedData, SaveFolder);
        src.entry = entry;
        src.representation = rep;
        return
    else
        error('Umitoolbox:genAmplitudeMaps:unsupportedExtension', ...
            'Unsupported file extension "%s".', ext);
    end
end

if isstruct(dataIn) && iLooksLikeUMT(dataIn)
    validateUMTStruct(dataIn, 'requireEventInfo', false);
    src.UMT = dataIn;
    [entry, rep] = iSelectUMTImageEntry(dataIn, SaveFolder);
    src.entry = entry;
    src.representation = rep;
    return
end

error('Umitoolbox:genAmplitudeMaps:unsupportedInput', ...
    ['Unsupported input type. Use numeric YXT, raw ".dat", UMT struct, ' ...
     'or ".umt" file.']);
end

function [entry, representation] = iSelectUMTImageEntry(umt, SaveFolder)
entryNames = fieldnames(umt.data);
validIdx = false(size(entryNames));
reps = strings(size(entryNames));
for i = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{i});
    if ~isstruct(thisEntry) || ~isscalar(thisEntry) || ...
            ~isfield(thisEntry, 'value') || ~isfield(thisEntry, 'dimNames')
        continue
    end
    dimNames = cellstr(string(thisEntry.dimNames));
    if isequal(dimNames, {'Y','X','T','E'})
        validIdx(i) = true;
        reps(i) = "eventsplit";
    elseif isequal(dimNames, {'Y','X','T'})
        validIdx(i) = true;
        reps(i) = "continuous";
    end
end
matchNames = entryNames(validIdx);
assert(~isempty(matchNames), ...
    'Umitoolbox:genAmplitudeMaps:noCompatibleUMTEntry', ...
    ['No compatible image entry was found in the UMT input. Supported ' ...
     'dimNames are YXT and YXTE.']);
assert(isscalar(matchNames), ...
    'Umitoolbox:genAmplitudeMaps:multipleCompatibleUMTEntries', ...
    ['Multiple compatible UMT entries were found. The current version ' ...
     'supports exactly one compatible image entry.']);
entry = umt.data.(matchNames{1});
representation = char(reps(find(validIdx, 1, 'first')));
if strcmp(representation, 'continuous')
    assert(isfile(fullfile(SaveFolder, 'events.mat')), ...
        'Umitoolbox:genAmplitudeMaps:missingEventsFile', ...
        ['UMT continuous YXT input requires an "events.mat" file in ' ...
         'SaveFolder.']);
end
end

function ampMap = iComputeAmplitudeFromYXTE(dataYXTE, baselineFrames, responseFrames, baselineMeasure, responseMeasure, eventIDs)
validateattributes(dataYXTE, {'numeric'}, {'nonempty'}, 'iComputeAmplitudeFromYXTE', 'dataYXTE');
assert(ndims(dataYXTE) == 4, ...
    'Umitoolbox:genAmplitudeMaps:invalidEventSplitData', ...
    'Event-split input must have dimensions YXTE.');

eventIDs = double(eventIDs(:));
assert(numel(eventIDs) == size(dataYXTE, 4), ...
    'Umitoolbox:genAmplitudeMaps:eventAxisMismatch', ...
    'Number of event IDs must match the E dimension.');

uniqueIDs = unique(eventIDs, 'stable');
ampMap = zeros(size(dataYXTE,1), size(dataYXTE,2), numel(uniqueIDs), 'single');
for iEv = 1:numel(uniqueIDs)
    idxE = eventIDs == uniqueIDs(iEv);
    thisData = dataYXTE(:,:,:,idxE);
    baselineVals = thisData(:,:,baselineFrames,:);
    responseVals = thisData(:,:,responseFrames,:);
    baselineVals = reshape(baselineVals, size(thisData,1), size(thisData,2), []);
    responseVals = reshape(responseVals, size(thisData,1), size(thisData,2), []);
    baselineMap = iApplyAggFcnND(baselineVals, baselineMeasure);
    responseMap = iApplyAggFcnND(responseVals, responseMeasure);
    ampMap(:,:,iEv) = responseMap - baselineMap;
end
end

function out = iApplyAggFcnND(vals, aggfcn)
switch lower(char(string(aggfcn)))
    case 'mean'
        out = mean(vals, 3, 'omitnan');
    case 'median'
        out = median(vals, 3, 'omitnan');
    case 'max'
        out = max(vals, [], 3, 'omitnan');
    case 'min'
        out = min(vals, [], 3, 'omitnan');
    otherwise
        error('Umitoolbox:genAmplitudeMaps:invalidAggFcn', ...
            'Unknown aggregate function "%s".', char(string(aggfcn)));
end
out = single(out);
end

function outData = iBuildOutputUMT(ampMap, eventInfoOut, baselineMeasure, responseMeasure, timeWindowSec)
meta = struct();
meta.BaselineMeasure = char(string(baselineMeasure));
meta.ResponseMeasure = char(string(responseMeasure));
if ischar(timeWindowSec) || (isstring(timeWindowSec) && isscalar(timeWindowSec))
    meta.TimeWindow_sec = char(string(timeWindowSec));
else
    meta.TimeWindow_sec = double(timeWindowSec(:).');
end

outData = genUMTStruct( ...
    single(ampMap), ...
    'kind', 'image', ...
    'entryName', 'AmplitudeMap', ...
    'dimNames', {'Y','X','E'}, ...
    'meta', meta);

outData = appendUMTEventInfo(outData, ...
    'eventInfo', eventInfoOut, ...
    'overwrite', true);
end

function tf = iValidateAggName(x)
tf = ischar(x) || (isstring(x) && isscalar(x));
if tf
    tf = ismember(lower(char(string(x))), {'mean','median','min','max'});
end
end

function tf = iValidateTimeWindowInput(x)
if ischar(x) || (isstring(x) && isscalar(x))
    tf = strcmpi(char(string(x)), 'all');
    return
end
tf = isnumeric(x) && isvector(x) && numel(x) == 2 && ...
     all(isfinite(x(:))) && all(x(:) >= 0);
end

function [baselineFrames, responseFrames] = iResolveAnalysisFrames(nFramesPerTrial, baselineFrames, frameRateHz, timeWindowSec)
baselineFrames = baselineFrames(:).';
baselineFrames = baselineFrames(baselineFrames >= 1 & baselineFrames <= nFramesPerTrial);
assert(~isempty(baselineFrames), ...
    'Umitoolbox:genAmplitudeMaps:invalidBaselineFrames', ...
    'No valid baseline frames were found.');

if ischar(timeWindowSec) || (isstring(timeWindowSec) && isscalar(timeWindowSec))
    assert(strcmpi(char(string(timeWindowSec)), 'all'), ...
        'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
        'TimeWindow_sec must be "all" or a numeric [start end] vector.');
    responseFrames = (baselineFrames(end) + 1):nFramesPerTrial;
else
    assert(isnumeric(timeWindowSec) && isvector(timeWindowSec) && numel(timeWindowSec) == 2, ...
        'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
        'TimeWindow_sec must be "all" or a numeric [start end] vector.');
    timeWindowSec = double(timeWindowSec(:).');
    assert(all(isfinite(timeWindowSec)) && all(timeWindowSec >= 0), ...
        'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
        'TimeWindow_sec values must be finite and non-negative.');
    assert(timeWindowSec(1) <= timeWindowSec(2), ...
        'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
        'TimeWindow_sec must satisfy start <= end.');
    % Frame baselineFrames(end)+1 is the first response frame, i.e. t = 0
    % after stimulus onset. Anchoring on baselineFrames(end) instead would
    % reject the natural "0 to N seconds" request that this parameter
    % advertises as allowed.
    frOn = baselineFrames(end) + 1 + round(timeWindowSec(1) * frameRateHz);
    frOff = baselineFrames(end) + 1 + round(timeWindowSec(2) * frameRateHz);
    assert(frOn >= baselineFrames(end) + 1, ...
        'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
        'TimeWindow_sec starts before the post-baseline response period.');
    assert(frOff <= nFramesPerTrial, ...
        'Umitoolbox:genAmplitudeMaps:invalidTimeWindow', ...
        'TimeWindow_sec extends beyond the available trial duration.');
    responseFrames = frOn:frOff;
end
assert(~isempty(responseFrames), ...
    'Umitoolbox:genAmplitudeMaps:emptyResponseFrames', ...
    'No response frames were selected for amplitude-map calculation.');
end

function tf = iLooksLikeUMT(x)
tf = isstruct(x) && isscalar(x) && ...
    isfield(x, 'version') && isfield(x, 'kind') && isfield(x, 'data');
end
