function outData = genAmplitudeMaps(data, SaveFolder, varargin)
%GENAMPLITUDEMAPS Compute event-wise response amplitude maps from event-split data.
%
%   outData = genAmplitudeMaps(data, SaveFolder)
%   outData = genAmplitudeMaps(data, SaveFolder, 'BaselineMeasure', value, ...)
%
%   This function computes response amplitude maps by subtracting an
%   aggregate baseline value from an aggregate response value along the time
%   dimension for each event condition.
%
%   Supported input:
%       Event-split .dat file with axes Y-X-T-E (for example the output of
%       split_data_by_event). Continuous Y-X-T data must be split first;
%       arrays, UMT structs, and .umt files are not supported.
%
%   Event handling:
%       The E axis of the file is matched to the events.mat in SaveFolder
%       (resolveDatEventMapping):
%         - one slice per event instance: the instances of each condition
%           are pooled into one map (ignored instances are excluded);
%         - one slice per condition (an aggregated file): one map per slice.
%       An E axis that cannot be matched to events.mat is rejected. Without
%       an events.mat every slice is treated as a repetition of one
%       condition, with a default baseline of 7 frames.
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
%                           it from the data; the .dat header provides it
%                           otherwise. AcqInfos.mat is not used.
%
%   Output:
%       outData - Full path of the .dat output ("amplitudeMap.dat" in
%                 SaveFolder): axes Y-X-E, one amplitude map per event
%                 condition. The .dat stores no labels: its E axis is
%                 matched to events.mat by resolveDatEventMapping (one slice
%                 per condition).
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
%       - The file is streamed in X slabs and only the baseline and response
%         frames of each trial are read; the output is written slab by slab
%         (Low-RAM mode is always on; the slab width follows the available
%         RAM).

default_Output = 'amplitudeMap.dat';

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

dataFile = iResolveDatFile(data, SaveFolder);
Info = loadMetaData(dataFile);
assertDatLayout(Info, {{'Y','X','T','E'}}, 'genAmplitudeMaps');

% Frame rate of the data itself: the explicit FrameRateHz (injected by
% PipelineManager), else the .dat header. AcqInfos.mat is not used.
frameRateHz = resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, dataFile, mfilename);

outData = iRunChunkedDat(dataFile, Info, baselineMeasure, responseMeasure, ...
    timeWindowSec, SaveFolder, frameRateHz, default_Output);

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            mfilename, ...
            ['Compute response amplitude maps by subtracting an ' ...
             'aggregate baseline from an aggregate response period.']);
        info.version = '2.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'Event-split .dat file with axes Y-X-T-E.', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'file', ...
            'position', 1, ...
            'callType', 'positional');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing events.mat and receiving the .dat output.', ...
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
            'Y-X-E amplitude maps (.dat), one per event condition.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

function dataFile = iResolveDatFile(dataIn, SaveFolder)
%IRESOLVEDATFILE Resolve the input to an existing .dat path.
assert(ischar(dataIn) || (isstring(dataIn) && isscalar(dataIn)), ...
    'Umitoolbox:genAmplitudeMaps:unsupportedInput', ...
    'Unsupported input type. Use the name or path of an event-split ".dat" file.');

dataFile = char(string(dataIn));
if ~isfile(dataFile)
    altPath = fullfile(SaveFolder, dataFile);
    assert(isfile(altPath), ...
        'Umitoolbox:genAmplitudeMaps:inputFileNotFound', ...
        'Input file "%s" was not found.', char(string(dataIn)));
    dataFile = altPath;
end

[~, ~, ext] = fileparts(dataFile);
assert(strcmpi(ext, '.dat'), ...
    'Umitoolbox:genAmplitudeMaps:unsupportedExtension', ...
    'Unsupported file extension "%s". Only .dat files are supported.', ext);
end

function outFile = iRunChunkedDat(dataFile, Info, baselineMeasure, responseMeasure, timeWindowSec, SaveFolder, frameRateHz, defaultOutput)
%IRUNCHUNKEDDAT Amplitude maps of a Y-X-T-E .dat, read and written in X slabs.

Ny = datAxisSize(Info, 'Y');
Nx = datAxisSize(Info, 'X');
Nt = datAxisSize(Info, 'T');
Ne = datAxisSize(Info, 'E');

mapping = resolveDatEventMapping(Info, SaveFolder);
assert(~strcmpi(mapping.status, 'mismatch'), ...
    'Umitoolbox:genAmplitudeMaps:eventMappingMismatch', ...
    'The E axis of "%s" cannot be matched to events.mat: %s', dataFile, mapping.message);
if strcmpi(mapping.status, 'noEvents')
    warning('Umitoolbox:genAmplitudeMaps:eventsNotMatched', '%s', mapping.message);
end

if isfield(mapping.eventInfo, 'baselinePeriod') && ~isempty(mapping.eventInfo.baselinePeriod)
    baselineFrames = 1:round(double(mapping.eventInfo.baselinePeriod) * frameRateHz);
else
    baselineFrames = 1:min(7, Nt);
end
[baselineFrames, responseFrames] = iResolveAnalysisFrames( ...
    Nt, baselineFrames, frameRateHz, timeWindowSec);

% Only the baseline and response frames of each trial are read: the frames
% of slice e are (e-1)*Nt + t, as E is the last axis.
usedCols = [baselineFrames, responseFrames];
nUsed = numel(usedCols);
frameIdx = usedCols(:) + (0:Ne-1) * Nt;   % [nUsed, Ne]
nB = numel(baselineFrames);
bIdx = 1:nB;
rIdx = nB + (1:numel(responseFrames));

if strcmpi(mapping.status, 'aggregated')
    % One slice per condition already: one map per slice.
    nOut = Ne;
    reduceSlab = @(slabE) iPerSliceAmplitude(slabE, bIdx, rIdx, baselineMeasure, responseMeasure);
else
    plan = EventsManager.conditionAggregationPlan(mapping.eventInfo);
    nOut = numel(plan.conditionID);
    reducer = iAmplitudeReducer(bIdx, rIdx, baselineMeasure, responseMeasure);
    reduceSlab = @(slabE) iConditionsToE(EventsManager.reduceByCondition(slabE, plan, reducer, 4));
end

slabIn = spatialSlabIO('open', dataFile, 'Info', Info);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

% Write through a scratch file so the output only appears once the run has
% completed, and so the input can be the very file the output replaces.
outFile = fullfile(SaveFolder, defaultOutput);
[~, outStem, outExt] = fileparts(outFile);
tmpFile = fullfile(SaveFolder, [outStem, '_writing', outExt]);
cTmp = onCleanup(@() iDeleteIfExists(tmpFile));
slabOut = spatialSlabIO('create', tmpFile, datHeaderFromInfo(Info, outStem, ...
    'dataClass', 'single', 'dimNames', {'Y','X','E'}, ...
    'dimSizes', [Ny, Nx, nOut], 'frameRateHz', NaN));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

% Scratch per X column: the frames read and their per-instance copy.
scratchBytes = double(Ny) * double(Nx) * double(nUsed) * double(Ne) * 2 * ...
    double(getByteSize('single'));
nChunks = calculateMaxChunkSize(scratchBytes, 1, 0.2);
chunkX = max(1, ceil(Nx / nChunks));
nChunks = ceil(Nx / chunkX);

for c = 1:nChunks
    xIdx = ((c-1) * chunkX + 1):min(c * chunkX, Nx);

    fprintf('Chunk %i/%i [Reading file ...]\n', c, nChunks)
    frames = single(spatialSlabIO('read', slabIn, xIdx, frameIdx(:).'));
    slabE = reshape(frames, Ny, numel(xIdx), nUsed, Ne);

    fprintf('Chunk %i/%i [Computing amplitude maps ...]\n', c, nChunks)
    amp = reduceSlab(slabE);

    fprintf('Chunk %i/%i [Writing to file ...]\n', c, nChunks)
    spatialSlabIO('write', slabOut, xIdx, reshape(single(amp), Ny, numel(xIdx), nOut));
    fprintf('Chunk %i/%i [Completed]\n', c, nChunks)
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close both files before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:genAmplitudeMaps:outputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end

function iDeleteIfExists(filePath)
%IDELETEIFEXISTS Remove a scratch file left by a failed run.
if isfile(filePath)
    delete(filePath);
end
end

function ampMap = iPerSliceAmplitude(slabE, baselineIdx, responseIdx, baselineMeasure, responseMeasure)
%IPERSLICEAMPLITUDE One amplitude map per E slice (an already aggregated file).
ampMap = zeros(size(slabE, 1), size(slabE, 2), size(slabE, 4), 'single');
for iSlice = 1:size(slabE, 4)
    ampMap(:, :, iSlice) = iAmplitudeOfTrials(slabE(:, :, :, iSlice), ...
        baselineIdx, responseIdx, baselineMeasure, responseMeasure);
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
