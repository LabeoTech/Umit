function outData = genAmplitudeMaps(data, SaveFolder, varargin)
%GENAMPLITUDEMAPS Compute event-wise response amplitude maps from imaging data.
%
%   outData = genAmplitudeMaps(data, SaveFolder)
%   outData = genAmplitudeMaps(data, SaveFolder, 'BaselineMeasure', value, ...)
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

% Legacy pipeline placeholder
default_Output = 'amplitudeMap.umt'; %#ok<NASGU>

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

if src.isRawDat
    outData = iRunLowRAM(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder);
else
    outData = iRunStandard(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder);
end

if nargout == 0 %#ok<UNRCH>
    clear outData
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
            'Folder containing acquisition metadata and/or events.mat.', ...
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

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'UMT structure containing YXE amplitude maps.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

function outData = iRunStandard(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder)

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
        dataYXTE = single(ev.splitDataByEvents(dataYXT));

        frameRateHz = double(ev.AcqInfo.FrameRateHz);
        baselineFrames = 1:round(double(ev.baselinePeriod) * frameRateHz);
        [baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
            size(dataYXTE, 3), baselineFrames, frameRateHz, timeWindowSec);

        eventIDsPerInstance = double(ev.eventID(ev.state));
        uniqueIDs = unique(eventIDsPerInstance, 'stable');

        eventInfoOut = struct();
        eventInfoOut.eventID = reshape(uniqueIDs, [], 1);
        eventInfoOut.repetitionIndex = zeros(numel(uniqueIDs), 1, 'uint16');
        eventInfoOut.eventName = reshape(string(ev.eventNameList(uniqueIDs)), [], 1);
        eventInfoOut.eventAxisMode = 'aggregated_repetitions';
        eventInfoOut.baselinePeriod = double(ev.baselinePeriod);

        ampMap = iComputeAmplitudeFromYXTE( ...
            dataYXTE, baselineFrames, responseFrames, ...
            baselineMeasure, responseMeasure, eventIDsPerInstance);

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

        assert(isfield(entry, 'meta') && isstruct(entry.meta) && ...
            isfield(entry.meta, 'FrameRateHz') && ~isempty(entry.meta.FrameRateHz), ...
            'Umitoolbox:genAmplitudeMaps:missingFrameRate', ...
            'UMT event-split input requires entry.meta.FrameRateHz.');
        frameRateHz = double(entry.meta.FrameRateHz);

        if isfield(eventInfo, 'baselinePeriod') && ~isempty(eventInfo.baselinePeriod)
            baselineFrames = 1:round(double(eventInfo.baselinePeriod) * frameRateHz);
        else
            baselineFrames = 1:min(7, size(dataYXTE, 3));
        end

        [baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
            size(dataYXTE, 3), baselineFrames, frameRateHz, timeWindowSec);

        eventIDsRaw = double(eventInfo.eventID(:));
        eventNamesRaw = string(eventInfo.eventName(:));
        axisMode = char(string(eventInfo.eventAxisMode));

        if strcmpi(axisMode, 'aggregated_repetitions')
            eventIDsForComputation = eventIDsRaw;
            eventIDsOut = eventIDsRaw;
            eventNamesOut = eventNamesRaw;
        elseif strcmpi(axisMode, 'instances')
            eventIDsForComputation = eventIDsRaw;
            [eventIDsOut, ia] = unique(eventIDsRaw, 'stable');
            eventNamesOut = eventNamesRaw(ia);
        else
            error('Umitoolbox:genAmplitudeMaps:unsupportedEventAxisMode', ...
                'Unsupported eventAxisMode "%s".', axisMode);
        end

        eventInfoOut = struct();
        eventInfoOut.eventID = reshape(eventIDsOut, [], 1);
        eventInfoOut.repetitionIndex = zeros(numel(eventIDsOut), 1, 'uint16');
        eventInfoOut.eventName = reshape(eventNamesOut, [], 1);
        eventInfoOut.eventAxisMode = 'aggregated_repetitions';
        if isfield(eventInfo, 'baselinePeriod')
            eventInfoOut.baselinePeriod = double(eventInfo.baselinePeriod);
        end

        ampMap = iComputeAmplitudeFromYXTE( ...
            dataYXTE, baselineFrames, responseFrames, ...
            baselineMeasure, responseMeasure, eventIDsForComputation);

    otherwise
        error('Umitoolbox:genAmplitudeMaps:unknownRepresentation', ...
            'Unknown input representation "%s".', src.representation);
end

outData = iBuildOutputUMT( ...
    ampMap, eventInfoOut, baselineMeasure, responseMeasure, timeWindowSec);

end

function outData = iRunLowRAM(src, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder)
Info = src.Info;
Ny = double(Info.Height);
Nx = double(Info.Width);
Nt = double(Info.datLength);

ev = EventsManager(SaveFolder);
[frMat, conditionIDlist] = ev.getFrameMatrix(Nt);
if isempty(frMat)
    error('Umitoolbox:genAmplitudeMaps:noFrames', ...
        'No event frames were returned by EventsManager.getFrameMatrix.');
end

frameRateHz = double(ev.AcqInfo.FrameRateHz);
baselineFrames = 1:round(double(ev.baselinePeriod) * frameRateHz);
[baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
    size(frMat, 2), baselineFrames, frameRateHz, timeWindowSec);

eventIDs = unique(conditionIDlist(:), 'stable');
nEvents = numel(eventIDs);
ampMap = zeros(Ny, Nx, nEvents, 'single');

fid = fopen(src.fileName, 'r');
if fid == -1
    error('Umitoolbox:genAmplitudeMaps:fileOpenFailed', ...
        'Could not open file "%s".', src.fileName);
end
cleanupFid = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

bytesPerElement = getByteSize('single');
nChunks = calculateMaxChunkSize(Ny * Nx * Nt * bytesPerElement, 2, 0.2);
nChunks = max(1, nChunks);
chunkX = ceil(Nx / nChunks);

for iEv = 1:nEvents
    idxTrials = find(conditionIDlist == eventIDs(iEv));
    trialMat = frMat(idxTrials, :);

    for iChunk = 1:nChunks
        xStart = (iChunk - 1) * chunkX + 1;
        xEnd = min(xStart + chunkX - 1, Nx);
        xIdx = xStart:xEnd;

        baselineVals = zeros(Ny, numel(xIdx), numel(baselineFrames), size(trialMat,1), 'single');
        responseVals = zeros(Ny, numel(xIdx), numel(responseFrames), size(trialMat,1), 'single');

        for iTrial = 1:size(trialMat,1)
            thisFrames = trialMat(iTrial, :);
            for iB = 1:numel(baselineFrames)
                fr = thisFrames(baselineFrames(iB));
                if ~isnan(fr) && fr >= 1 && fr <= Nt
                    baselineVals(:,:,iB,iTrial) = iReadFrameSlab(fid, Ny, Nx, double(fr), xIdx);
                else
                    baselineVals(:,:,iB,iTrial) = NaN;
                end
            end
            for iR = 1:numel(responseFrames)
                fr = thisFrames(responseFrames(iR));
                if ~isnan(fr) && fr >= 1 && fr <= Nt
                    responseVals(:,:,iR,iTrial) = iReadFrameSlab(fid, Ny, Nx, double(fr), xIdx);
                else
                    responseVals(:,:,iR,iTrial) = NaN;
                end
            end
        end

        baselineVals = reshape(baselineVals, Ny, numel(xIdx), []);
        responseVals = reshape(responseVals, Ny, numel(xIdx), []);

        baselineMap = iApplyAggFcnND(baselineVals, baselineMeasure);
        responseMap = iApplyAggFcnND(responseVals, responseMeasure);
        ampMap(:, xIdx, iEv) = responseMap - baselineMap;
    end
end

eventNames = ev.eventNameList(:);
if isempty(eventNames)
    eventNames = arrayfun(@(x) sprintf('Event%d', x), eventIDs, 'UniformOutput', false);
end

eventInfoOut = struct();
eventInfoOut.eventID = reshape(eventIDs, [], 1);
eventInfoOut.repetitionIndex = zeros(nEvents, 1, 'uint16');
eventInfoOut.eventName = reshape(string(eventNames(eventIDs)), [], 1);
eventInfoOut.eventAxisMode = 'aggregated_repetitions';
eventInfoOut.baselinePeriod = double(ev.baselinePeriod);

outData = iBuildOutputUMT( ...
    ampMap, eventInfoOut, baselineMeasure, responseMeasure, timeWindowSec);
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
    [loadedData, Info] = loadData(fileName);
    [~,~,ext] = fileparts(fileName);
    ext = lower(ext);
    if strcmp(ext, '.dat')
        src.representation = 'continuous';
        src.fileName = fileName;
        src.Info = Info;
        src.isRawDat = true;
        src.data = [];
        return
    elseif strcmp(ext, '.umt')
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
assert(numel(matchNames) == 1, ...
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
    idxE = find(eventIDs == uniqueIDs(iEv));
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

function slab = iReadFrameSlab(fid, Ny, Nx, frameIdx, xIdx)
bytesPerElement = getByteSize('single');
slab = zeros(Ny, numel(xIdx), 'single');
for iX = 1:numel(xIdx)
    x = xIdx(iX);
    offset = ((frameIdx - 1) * Ny * Nx + (x - 1) * Ny) * bytesPerElement;
    fseek(fid, offset, 'bof');
    tmp = fread(fid, Ny, '*single');
    if numel(tmp) ~= Ny
        error('Umitoolbox:genAmplitudeMaps:fileReadFailed', ...
            'Unexpected end of file while reading frame %d.', frameIdx);
    end
    slab(:, iX) = tmp;
end
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
