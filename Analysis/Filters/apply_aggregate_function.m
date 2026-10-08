function outData = apply_aggregate_function(data, SaveFolder, varargin)
%APPLY_AGGREGATE_FUNCTION Aggregate image data along T or E.
%
%   outData = apply_aggregate_function(data, SaveFolder)
%   outData = apply_aggregate_function(data, SaveFolder, ...
%       'aggregateFcn', aggFcn, 'dimensionName', dimName)
%
% Inputs:
%   data       : One of:
%                1) Filename of a .dat file with axes Y-X-T, Y-X-T-E, or
%                   Y-X-E
%                2) Single-entry image UMT struct
%                3) Filename of a .umt file containing one single-entry
%                   image UMT struct
%                Arrays are not supported.
%
%   SaveFolder : Folder containing events.mat; also the folder of the .dat
%                output.
%
% Name-Value parameters:
%   aggregateFcn  : 'mean','median','std','max','min','sum'
%                   Default: 'mean'
%   dimensionName : 'T' or 'E'. The data must contain that axis.
%                   Default: 'T'
%   FrameRateHz   : Frame rate of DATA (Hz), needed to turn the event times
%                   of events.mat into frames when a continuous Y-X-T .dat
%                   is aggregated along E. PipelineManager injects it from
%                   the data; the .dat header provides it otherwise.
%                   AcqInfos.mat is not used.
%
% Output (same representation as the input):
%   .dat input  -> .dat file ("aggFcn_applied.dat" in SaveFolder); outData
%                  is its full path.
%   UMT input   -> UMT struct (from a .umt file as well).
%
%   Axes of the result:
%       T aggregation   Y-X-T   -> Y-X
%                       Y-X-T-E -> Y-X-E (each trial reduced on its own)
%       E aggregation   Y-X-T-E, Y-X-E -> same axes, one E slice per
%                       condition
%                       Y-X-T   -> Y-X-T-E: the instances of events.mat are
%                       split out of the continuous recording and reduced
%                       per condition (trial length as split_data_by_event)
%
% Notes:
%   - A .dat input is streamed in X slabs and its output written slab by
%     slab, so neither the recording nor the result is ever resident whole
%     (Low-RAM mode is always on; the slab size follows the available RAM).
%   - E aggregation follows EventsManager.conditionAggregationPlan and
%     reduceByCondition (.dat header Phase 8c): ignored instances are
%     excluded, conditions are ordered by first appearance, and a condition
%     whose instances are all ignored gives a NaN slice. The .dat output
%     stores no labels: its E axis is matched to events.mat by
%     resolveDatEventMapping (one slice per condition). An event-split .dat
%     that is already aggregated is refused.
%   - UMT outputs carry the aggregated eventInfo (eventAxisMode
%     'aggregated_repetitions', repetitionIndex 0, selected, durationSec,
%     nInstances). UMT E aggregation uses the frozen eventInfo stored in the
%     UMT. A .umt file is loaded fully into RAM.
%
% See also: spatialSlabIO, genUMTStruct, appendUMTEventInfo

default_Output = 'aggFcn_applied.umt';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'aggregateFcn', 'mean', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'dimensionName', 'T', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'FrameRateHz', []);

parse(p, data, SaveFolder, varargin{:});
explicitRate = p.Results.FrameRateHz;

aggFcn = lower(char(string(p.Results.aggregateFcn)));
dimName = upper(char(string(p.Results.dimensionName)));
SaveFolder = char(string(p.Results.SaveFolder));

validAggFcns = {'mean','median','std','max','min','sum'};
if ~ismember(aggFcn, validAggFcns)
    error('apply_aggregate_function:InvalidAggregateFcn', ...
        'Unsupported aggregateFcn "%s".', aggFcn);
end

if ~ismember(dimName, {'T','E'})
    error('apply_aggregate_function:InvalidDimensionName', ...
        'dimensionName must be either ''T'' or ''E''.');
end

if ~isfolder(SaveFolder)
    error('apply_aggregate_function:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Case 1: file input
% -------------------------------------------------------------------------
if ischar(data) || (isstring(data) && isscalar(data))

    dataFile = char(string(data));

    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('apply_aggregate_function:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);

    switch ext
        case '.dat'
            outData = iAggregateDatFile(dataFile, SaveFolder, default_Output, ...
                aggFcn, dimName, explicitRate);
            return

        case '.umt'
            warning('apply_aggregate_function:UMTFileLoadsInRAM', ...
                ['RAM-Safe mode is not available for data stored in this format. ' ...
                 'Loading the UMT content into RAM.']);
            data = loadData(dataFile);

        otherwise
            error('apply_aggregate_function:UnsupportedInputFile', ...
                'Unsupported input file extension "%s". Only .dat and .umt files are supported.', ext);
    end
end

% -------------------------------------------------------------------------
% Case 2: UMT struct in RAM
% -------------------------------------------------------------------------
if ~(isstruct(data) && isscalar(data))
    error('apply_aggregate_function:UnsupportedInputType', ...
        ['Input "data" must be a .dat filename, a UMT struct, ' ...
         'or a .umt filename containing a UMT struct.']);
end

[entryNames, entryData, entryDims, sourceLabels, sourceEventInfo, entryMetas] = ...
    iExtractValidUMTData(data);

bInputUsesE = false;
for iEntry = 1:numel(entryDims)
    if any(strcmp(entryDims{iEntry}, 'E'))
        bInputUsesE = true;
        break
    end
end

if bInputUsesE && isempty(fieldnames(sourceEventInfo))
    error('apply_aggregate_function:MissingEventInfo', ...
        ['Operation aborted. The input UMT contains entries with an E dimension ' ...
         'but has no shared top-level eventInfo.']);
end

if strcmpi(dimName, 'E')
    if isempty(fieldnames(sourceEventInfo))
        error('apply_aggregate_function:MissingEventInfo', ...
            'UMT E aggregation requires a shared top-level eventInfo field.');
    end

    if ~strcmpi(sourceEventInfo.eventAxisMode, 'instances')
        error('apply_aggregate_function:InvalidEventAxisMode', ...
            ['UMT E aggregation requires eventAxisMode = "instances". ' ...
             'The current eventInfo already represents aggregated repetitions.']);
    end

    plan = EventsManager.conditionAggregationPlan(sourceEventInfo);
end

outEntryData = entryData;
outEntryDims = entryDims;
bAnyAggregated = false;

for iEntry = 1:numel(entryNames)

    thisData = single(entryData{iEntry});
    thisDims = entryDims{iEntry};
    idxTarget = find(strcmp(thisDims, dimName), 1, 'first');

    if isempty(idxTarget)
        continue
    end

    bAnyAggregated = true;

    switch dimName

        case 'T'
            permOrder = [idxTarget setdiff(1:ndims(thisData), idxTarget)];
            dataP = permute(thisData, permOrder);
            szP = size(dataP);
            dataP = reshape(dataP, szP(1), []);

            aggFlat = iCalcAgg(dataP, aggFcn);
            outEntryData{iEntry} = reshape(single(aggFlat), szP(2:end));

            if isvector(outEntryData{iEntry}) && ~isscalar(outEntryData{iEntry})
                outEntryData{iEntry} = outEntryData{iEntry}(:);
            end

            newDims = thisDims;
            newDims(idxTarget) = [];
            outEntryDims{iEntry} = newDims;

        case 'E'
            assert(numel(sourceEventInfo.eventID) == size(thisData, idxTarget), ...
                'apply_aggregate_function:EventAxisMismatch', ...
                ['Entry "%s" has %d elements along dimension "E", but ' ...
                 'sourceEventInfo.eventID has %d elements.'], ...
                entryNames{iEntry}, size(thisData, idxTarget), numel(sourceEventInfo.eventID));
            outEntryData{iEntry} = EventsManager.reduceByCondition(thisData, plan, ...
                @(x) iCalcAgg(x, aggFcn, idxTarget), idxTarget);
            outEntryDims{iEntry} = thisDims;
    end
end

if ~bAnyAggregated
    error('apply_aggregate_function:NoValidEntries', ...
        ['No valid image-backed UMT entry contains the requested ' ...
         'dimension "%s".'], dimName);
end

outLabels = sourceLabels;

outEventInfo = struct();
bOutputUsesE = false;

for iEntry = 1:numel(outEntryDims)
    if any(strcmp(outEntryDims{iEntry}, 'E'))
        bOutputUsesE = true;
        break
    end
end

if bOutputUsesE
    if strcmpi(dimName, 'T')
        outEventInfo = sourceEventInfo;
    else
        outEventInfo = plan.eventInfoOut;
    end
end

outData = iPackageOutputUMT( ...
    entryNames, ...
    outEntryData, ...
    outEntryDims, ...
    outLabels, ...
    outEventInfo, ...
    entryMetas);

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            ['Aggregate image data along T or E. A .dat input gives a .dat ' ...
             'output, a UMT input gives a UMT output.']);

        info.version = '2.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: .dat file (Y-X-T, Y-X-T-E or Y-X-E), ' ...
             'or a .umt file / UMT struct with one image entry.'], ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'file');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing events.mat and receiving the .dat output.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'aggregateFcn', ...
            'parameter', ...
            'Aggregation function applied along the requested dimension.', ...
            'kind', 'parameter', ...
            'default', 'mean', ...
            'allowed', {'mean','median','std','max','min','sum'}, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'dimensionName', ...
            'parameter', ...
            'Dimension to aggregate: T or E.', ...
            'kind', 'parameter', ...
            'default', 'T', ...
            'allowed', {'T','E'}, ...
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
            ['Aggregated output, in the representation of the input: a .dat ' ...
             'file for .dat input, else a UMT struct (saved as .umt).'], ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: streamed aggregation of a .dat file into a .dat file
% =========================================================================
function outFile = iAggregateDatFile(dataFile, SaveFolder, defaultOutput, aggFcn, dimName, explicitRate)
%IAGGREGATEDATFILE Aggregate a .dat along T or E, writing a .dat in X slabs.

meta = loadMetaData(dataFile);
assertDatLayout(meta, {{'Y','X','T'}, {'Y','X','T','E'}, {'Y','X','E'}}, ...
    'apply_aggregate_function');
inDims = cellstr(string(meta.dimNames(:).'));
hasT = any(strcmp(inDims, 'T'));
hasE = any(strcmp(inDims, 'E'));

if ~any(strcmp(inDims, dimName))
    % Only a continuous Y-X-T recording can be aggregated along E without
    % having an E axis (its instances come from events.mat).
    if ~(strcmp(dimName, 'E') && hasT)
        error('apply_aggregate_function:MissingDimension', ...
            'The file "%s" has no %s axis (axes: %s).', dataFile, dimName, strjoin(inDims, '-'));
    end
end

nY = datAxisSize(meta, 'Y');
nX = datAxisSize(meta, 'X');
nT = datAxisSize(meta, 'T');
nE = datAxisSize(meta, 'E');
frameRateHz = [];   % rate of the output header; [] keeps the input's

% Plan: output axes/sizes, and the per-slab reduction.
if strcmp(dimName, 'T')
    outDims = inDims(~strcmp(inDims, 'T'));
    outSizes = [nY, nX, nE(nE > 0)];
    frameRateHz = NaN;   % no T axis left
    reduceSlab = @(slab) iReduceT(slab, aggFcn, nY, nE);
    bytesPerX = nY * max(nT, 1) * max(nE, 1) * 4 * 2;

elseif hasE
    % Event-split file: reduce the E axis per condition.
    mapping = resolveDatEventMapping(meta, SaveFolder);
    if strcmpi(mapping.status, 'aggregated')
        % Raises Umitoolbox:EventsManager:alreadyAggregated.
        EventsManager.conditionAggregationPlan(mapping.eventInfo);
    end
    if ~strcmpi(mapping.status, 'matched')
        warning('Umitoolbox:apply_aggregate_function:eventsNotMatched', '%s', mapping.message);
    end
    plan = EventsManager.conditionAggregationPlan(mapping.eventInfo);
    idxE = find(strcmp(inDims, 'E'), 1);
    outDims = inDims;
    outSizes = meta.dimSizes(:).';
    outSizes(idxE) = numel(plan.conditionID);
    reduceSlab = @(slab) EventsManager.reduceByCondition(slab, plan, ...
        @(x) iCalcAgg(x, aggFcn, idxE), idxE);
    bytesPerX = nY * max(nT, 1) * (nE + numel(plan.conditionID)) * 4 * 2;

else
    % Continuous Y-X-T: split the instances of events.mat and reduce them.
    evObj = EventsManager(SaveFolder);
    frameRateHz = resolveDataInfoValue('frameRateHz', explicitRate, dataFile, 'apply_aggregate_function');
    [frMat, plan] = iEventFramePlan(evObj, nT, frameRateHz);
    nInst = size(frMat, 1);
    nCond = numel(plan.conditionID);
    trialLen = size(frMat, 2);
    outDims = {'Y','X','T','E'};
    outSizes = [nY, nX, trialLen, nCond];
    reduceSlab = @(slab) EventsManager.reduceByCondition( ...
        iInstancesFromFrames(slab, frMat), plan, @(x) iCalcAgg(x, aggFcn, 4), 4);
    % Input slab, all instances (and a permuted copy), and the aggregate.
    bytesPerX = nY * (nT + 2 * trialLen * nInst + trialLen * nCond) * 4;
end

slabIn = spatialSlabIO('open', dataFile, 'Info', meta);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

% Write through a scratch file so the output only appears once the run has
% completed, and so the input can be the very file the output replaces (a
% pipeline re-run).
outFile = fullfile(SaveFolder, strrep(defaultOutput, '.umt', '.dat'));
[~, outStem, outExt] = fileparts(outFile);
tmpFile = fullfile(SaveFolder, [outStem, '_writing', outExt]);
cTmp = onCleanup(@() iDeleteIfExists(tmpFile));
hdrOpts = {'dataClass', 'single', 'dimNames', outDims, 'dimSizes', outSizes};
if ~isempty(frameRateHz)
    hdrOpts = [hdrOpts, {'frameRateHz', frameRateHz}];
end
slabOut = spatialSlabIO('create', tmpFile, datHeaderFromInfo(meta, outStem, hdrOpts{:}));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

nChunks = calculateMaxChunkSize(double(nX) * bytesPerX, 1.5, 0.2);
chunkX = max(1, ceil(nX / nChunks));
nChunks = ceil(nX / chunkX);

for c = 1:nChunks
    xIdx = ((c-1) * chunkX + 1):min(c * chunkX, nX);

    fprintf('Chunk %i/%i [Reading file ...]\n', c, nChunks)
    slab = single(spatialSlabIO('read', slabIn, xIdx));

    fprintf('Chunk %i/%i [Aggregating %s ...]\n', c, nChunks, dimName)
    slab = reduceSlab(slab);

    fprintf('Chunk %i/%i [Writing to file ...]\n', c, nChunks)
    spatialSlabIO('write', slabOut, xIdx, reshape(slab, [nY, numel(xIdx), outSizes(3:end), 1]));
    fprintf('Chunk %i/%i [Completed]\n', c, nChunks)
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close both files before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'apply_aggregate_function:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end

function out = iReduceT(slab, aggFcn, nY, nE)
%IREDUCET Aggregate a [Y, nx, T(, E)] slab along T: [Y, nx(, E)].
% Reduce along the first dimension of a [T, Y*nx*E] view: the same
% accumulation order as the in-memory UMT path, so both give identical
% single-precision sums.
nTloc = size(slab, 3);
out = iCalcAgg(reshape(permute(slab, [3 1 2 4]), nTloc, []), aggFcn);
out = reshape(single(out), nY, size(slab, 2), max(nE, 1));
end

function iDeleteIfExists(filePath)
%IDELETEIFEXISTS Remove a scratch file left by a failed run.
if isfile(filePath)
    delete(filePath);
end
end

% =========================================================================
% Helper: event frames and aggregation plan of a continuous recording
% =========================================================================
function [frMat, plan] = iEventFramePlan(evObj, nT, frameRateHz)
%IEVENTFRAMEPLAN Frame matrix of every instance and its aggregation plan.
%
% Every instance is split, ignored ones included; the plan (EventsManager)
% excludes the ignored ones from each condition's aggregate (Phase 8c).

[frMat, conditionIDlist] = evObj.getFrameMatrix(nT, 'FrameRateHz', frameRateHz, ...
    'IncludeIgnored', true);
if isempty(frMat)
    error('apply_aggregate_function:NoEventsFound', ...
        'No valid events were found in events.mat for E aggregation.');
end
% Same trial length as EventsManager.splitDataByEvents (split_data_by_event):
% crop every trial from the first frame that any instance lacks.
firstNaNCol = find(any(isnan(frMat), 1), 1, 'first');
if ~isempty(firstNaNCol)
    frMat(:, firstNaNCol:end) = [];
end
evInfo = evObj.exportEventInfo('FrameRateHz', frameRateHz, 'IncludeIgnored', true);
assert(isequal(double(evInfo.eventID(:)), double(conditionIDlist(:))), ...
    'apply_aggregate_function:EventAxisMismatch', ...
    'The event list does not match the split trials.');
plan = EventsManager.conditionAggregationPlan(evInfo);
end

% =========================================================================
% Helper: Y x X x T x E trials of every instance from a frame matrix
% =========================================================================
function dataYXTE = iInstancesFromFrames(dataYXT, frMat)
%IINSTANCESFROMFRAMES Trials of DATAYXT (rows of FRMAT); NaN outside the data.

nInst = size(frMat, 1);
dataYXTE = nan(size(dataYXT, 1), size(dataYXT, 2), size(frMat, 2), nInst, 'single');
for iInst = 1:nInst
    validMask = ~isnan(frMat(iInst, :));
    if any(validMask)
        dataYXTE(:, :, validMask, iInst) = dataYXT(:, :, frMat(iInst, validMask));
    end
end
end

% =========================================================================
% Helper: Extract and validate image-backed data from a UMT structure
% =========================================================================
function [entryNames, entryData, entryDims, labels, eventInfo, entryMetas] = iExtractValidUMTData(umt)

validateUMTStruct(umt, 'requireEventInfo', false);

if ~strcmpi(umt.kind, 'image')
    error('apply_aggregate_function:InvalidUMTKind', ...
        ['Operation aborted. UMT input must have kind = "image". ' ...
         'This function does not support non-image UMT structures.']);
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error('apply_aggregate_function:EmptyUMTData', ...
        'Operation aborted. UMT data is empty.');
end

if ~isscalar(entryNames)
    error('apply_aggregate_function:multipleCompatibleUMTEntries', ...
        ['Operation aborted. The input UMT has %d entries; ' ...
         'only a single image entry is supported.'], numel(entryNames));
end

entryData = cell(size(entryNames));
entryDims = cell(size(entryNames));
entryMetas = cell(size(entryNames));

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    if ~all(ismember({'Y','X'}, thisDims))
        error('apply_aggregate_function:NonImageUMTEntry', ...
            ['Operation aborted. All entries in the input UMT must be ' ...
             'image-backed and contain dimensions Y and X.\n' ...
             'Invalid entry: "%s".'], ...
            entryNames{iEntry});
    end

    entryData{iEntry} = single(thisEntry.value);
    entryDims{iEntry} = thisDims;

    if isfield(thisEntry, 'meta') && isstruct(thisEntry.meta) && isscalar(thisEntry.meta)
        entryMetas{iEntry} = thisEntry.meta;
    else
        entryMetas{iEntry} = struct();
    end
end

if isfield(umt, 'labels')
    labels = umt.labels;
else
    labels = struct();
end

if isfield(umt, 'eventInfo')
    eventInfo = umt.eventInfo;
else
    eventInfo = struct();
end
end

% =========================================================================
% Helper: Package aggregated entries into an output UMT structure
% =========================================================================
function outUMT = iPackageOutputUMT(entryNames, entryData, entryDims, labelsIn, eventInfoIn, entryMetasIn)

outUMT = [];

labelsOut = struct();
if ~isempty(labelsIn) && isstruct(labelsIn)
    usedDims = {};
    for iEntry = 1:numel(entryDims)
        usedDims = [usedDims, entryDims{iEntry}]; %#ok<AGROW>
    end
    usedDims = unique(usedDims, 'stable');

    labelFields = fieldnames(labelsIn);
    for iField = 1:numel(labelFields)
        if ismember(labelFields{iField}, usedDims)
            labelsOut.(labelFields{iField}) = labelsIn.(labelFields{iField});
        end
    end
end

for iEntry = 1:numel(entryNames)

    if iEntry == 1
        if isempty(fieldnames(labelsOut))
            outUMT = genUMTStruct( ...
                entryData{iEntry}, ...
                'kind', 'image', ...
                'entryName', entryNames{iEntry}, ...
                'dimNames', entryDims{iEntry}, ...
                'meta', entryMetasIn{iEntry});
        else
            outUMT = genUMTStruct( ...
                entryData{iEntry}, ...
                'kind', 'image', ...
                'entryName', entryNames{iEntry}, ...
                'dimNames', entryDims{iEntry}, ...
                'labels', labelsOut, ...
                'meta', entryMetasIn{iEntry});
        end
    else
        outUMT = genUMTStruct( ...
            outUMT, ...
            'value', entryData{iEntry}, ...
            'entryName', entryNames{iEntry}, ...
            'dimNames', entryDims{iEntry}, ...
            'meta', entryMetasIn{iEntry});
    end
end

if ~isempty(eventInfoIn) && isstruct(eventInfoIn) && ~isempty(fieldnames(eventInfoIn))
    % Struct form: selected, durationSec, nInstances, and baselinePeriod
    % survive (Phase 8c).
    outUMT = appendUMTEventInfo(outUMT, ...
        'eventInfo', eventInfoIn, ...
        'overwrite', true);
else
    validateUMTStruct(outUMT, 'requireEventInfo', true);
end
end

% =========================================================================
% Helper: Core aggregation function
% =========================================================================
function out = iCalcAgg(vals, aggFcn, dim)
%ICALCAGG Reduce VALS along dimension DIM (default 1) using AGGFCN, omitting NaN.
%
% For a pixel/column that is entirely NaN, 'sum' returns 0 (sum's own
% 'omitnan' convention: an empty/all-NaN input sums to 0), while 'mean',
% 'median', 'std', 'max', and 'min' all return NaN. A fully-masked pixel
% therefore reads as a valid zero under 'sum' aggregation, not as missing
% data the way it does under every other aggregateFcn.

if nargin < 3
    dim = 1;
end

switch lower(aggFcn)
    case 'mean'
        out = mean(vals, dim, 'omitnan');

    case 'median'
        out = median(vals, dim, 'omitnan');

    case 'std'
        out = std(vals, 0, dim, 'omitnan');

    case 'max'
        out = max(vals, [], dim, 'omitnan');

    case 'min'
        out = min(vals, [], dim, 'omitnan');

    case 'sum'
        out = sum(vals, dim, 'omitnan');

    otherwise
        error('apply_aggregate_function:InvalidAggregateFcnInternal', ...
            'Unsupported aggregateFcn "%s".', aggFcn);
end
end
