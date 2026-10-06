function [outData, metaData] = apply_aggregate_function(data, SaveFolder, varargin)
%APPLY_AGGREGATE_FUNCTION Aggregate image-backed data along T or E.
%
%   outData = apply_aggregate_function(data, SaveFolder)
%   [outData, metaData] = apply_aggregate_function(data, SaveFolder, ...
%       'aggregateFcn', aggFcn, 'dimensionName', dimName)
%
% Inputs:
%   data       : One of:
%                1) Numeric 3-D array with dimensions Y x X x T
%                2) Filename to a .dat file storing Y x X x T data, or, for
%                   E aggregation, event-split Y-X-T-E / Y-X-E data whose
%                   E axis matches events.mat (resolveDatEventMapping;
%                   loaded fully into RAM)
%                3) UMT struct
%                4) Filename to a .umt file containing a UMT struct
%
%   SaveFolder : Folder containing AcqInfos.mat and events.mat.
%
% Name-Value parameters:
%   aggregateFcn  : 'mean','median','std','max','min','sum'
%                   Default: 'mean'
%   dimensionName : 'T' or 'E'
%                   Default: 'T'
%   FrameRateHz   : Frame rate of DATA (Hz), needed for E aggregation of
%                   raw arrays and .dat files (event times to frames).
%                   PipelineManager injects it from the data; a .dat
%                   input's header provides it otherwise. In-RAM arrays
%                   need it explicitly; AcqInfos.mat is not used.
%
% Output:
%   outData    : - E aggregation of a raw YXT array or a .dat file: numeric
%                  Y x X x T x E (Y x X x E for a Y-X-E .dat) array, one
%                  slice per event condition
%                  (saved as .dat by PipelineManager, .dat header Phase 8c).
%                - Otherwise: output UMT struct.
%   metaData   : struct with dimNames {'Y','X','T','E'} for
%                the numeric output (PipelineManager uses it to save the
%                .dat); empty struct otherwise.
%
% Notes:
%   - Raw YXT arrays and raw .dat files use live event information from
%     events.mat through EventsManager.
%   - UMT inputs use the frozen shared top-level eventInfo stored in the UMT.
%   - Raw .dat input uses spatially chunked reads. The complete aggregate
%     output remains resident in RAM: 4*Y*X bytes for T aggregation, or
%     4*Y*X*trialLen*nConditions bytes for E aggregation (single precision).
%     The E path sizes each slab for the simultaneous input, condition,
%     permutation, and aggregate workspaces.
%   - If a .umt file is provided, its content is loaded fully into RAM.
%   - E aggregation follows EventsManager.conditionAggregationPlan and
%     reduceByCondition (.dat header Phase 8c): ignored instances are
%     excluded, conditions are ordered by first appearance, and a condition
%     whose instances are all ignored gives a NaN slice. UMT outputs carry
%     the aggregated eventInfo (eventAxisMode 'aggregated_repetitions',
%     repetitionIndex 0, selected, durationSec, nInstances). The .dat output
%     stores no labels: its E axis is matched to events.mat by
%     resolveDatEventMapping (one slice per condition).
%   - UMT inputs keep UMT outputs, with their own eventInfo carried through.
%
% See also: spatialSlabIO, genUMTStruct, appendUMTEventInfo

default_Output = 'aggFcn_applied.umt';
metaData = struct();

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
% Case 1: raw YXT array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)

    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, ...
        mfilename, 'data');

    rawData = single(data);

    switch dimName

        case 'T'
            dataP = permute(rawData, [3 1 2]);
            dataP = reshape(dataP, size(rawData,3), []);
            aggFlat = iCalcAgg(dataP, aggFcn);
            aggData = reshape(single(aggFlat), size(rawData,1), size(rawData,2));

            outData = iPackageOutputUMT( ...
                {'main'}, ...
                {aggData}, ...
                {{'Y','X'}}, ...
                struct(), ...
                struct(), ...
                {struct()});

        case 'E'
            evObj = EventsManager(SaveFolder);
            frameRateHz = resolveDataInfoValue('frameRateHz', explicitRate, data, mfilename);
            [frMat, plan] = iEventFramePlan(evObj, size(rawData, 3), frameRateHz);

            dataYXTE = iInstancesFromFrames(rawData, frMat);
            outData = EventsManager.reduceByCondition(dataYXTE, plan, ...
                @(x) iCalcAgg(x, aggFcn, 4), 4);
            metaData = struct('dimNames', {{'Y','X','T','E'}});
    end

    return
end

% -------------------------------------------------------------------------
% Case 2: file input
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
            datMeta = loadMetaData(dataFile);
            if strcmp(dimName, 'E') && any(strcmp(cellstr(string(datMeta.dimNames)), 'E'))
                % Event-split .dat (e.g. split_data_by_event output).
                [outData, metaData] = iAggregateEventSplitDat(dataFile, datMeta, ...
                    SaveFolder, aggFcn);
                return
            end

            [aggData, outDimNames, labels, eventInfo] = ...
                iExecuteChunkedDat(dataFile, SaveFolder, aggFcn, dimName, explicitRate);

            if strcmp(dimName, 'E')
                % One slice per condition: saved as .dat (Phase 8c).
                outData = aggData;
                metaData = struct('dimNames', {outDimNames});
                return
            end

            outData = iPackageOutputUMT( ...
                {'main'}, ...
                {aggData}, ...
                {outDimNames}, ...
                labels, ...
                eventInfo, ...
                {struct()});
            return

        case {'.umt','.mat'}
            warning('apply_aggregate_function:UMTFileLoadsInRAM', ...
                ['RAM-Safe mode is not available for data stored in this format. ' ...
                 'Loading the UMT content into RAM.']);

            try
                loadedUMT = loadData(dataFile);
                if ~(isstruct(loadedUMT) && isscalar(loadedUMT) && ...
                        all(ismember({'version','kind','data'}, fieldnames(loadedUMT))))
                    error('Invalid UMT payload loaded.');
                end
            catch
                S = load(dataFile, '-mat');
                fn = fieldnames(S);
                loadedUMT = [];
                for iField = 1:numel(fn)
                    candidate = S.(fn{iField});
                    if isstruct(candidate) && isscalar(candidate) && ...
                            all(ismember({'version','kind','data'}, fieldnames(candidate)))
                        loadedUMT = candidate;
                        break
                    end
                end
                if isempty(loadedUMT)
                    error('apply_aggregate_function:NoUMTFoundInFile', ...
                        'No scalar UMT struct was found in "%s".', dataFile);
                end
            end

            data = loadedUMT;

        otherwise
            error('apply_aggregate_function:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

% -------------------------------------------------------------------------
% Case 3: UMT struct in RAM
% -------------------------------------------------------------------------
if ~isstruct(data)
    error('apply_aggregate_function:UnsupportedInputType', ...
        ['Input "data" must be a YXT array, a .dat filename, ' ...
         'a UMT struct, or a .umt filename containing a UMT struct.']);
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
            ['Aggregate raw image time-series or UMT image-backed data ' ...
             'along T or E and return a UMT struct.']);

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT array, .dat filename, ' ...
             'UMT struct, or .umt file containing one UMT struct.'], ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing AcqInfos.mat and events.mat.', ...
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
            ['Aggregated output: Y x X x T x E per-condition image data ' ...
             '(.dat) for E aggregation of raw data, else a UMT struct.'], ...
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

% =========================================================================
% Helper: Chunked raw-DAT input execution with an in-memory output
% =========================================================================
function [aggData, outDimNames, labels, eventInfo, frameRateHz] = iExecuteChunkedDat(dataFile, SaveFolder, aggFcn, dimName, explicitRate)

labels = struct();
eventInfo = struct();
frameRateHz = [];

meta = loadMetaData(dataFile);
assertDatLayout(meta, {{'Y','X','T'}}, 'apply_aggregate_function');
if ~isfield(meta, 'dimNames') || ~isfield(meta, 'dimSizes')
    error('apply_aggregate_function:InvalidMetaData', ...
        'loadMetaData did not return dimNames and dimSizes for "%s".', dataFile);
end

nY = datAxisSize(meta, 'Y');
nX = datAxisSize(meta, 'X');
nT = datAxisSize(meta, 'T');

% Conservative fixed chunk-size budget (not derived from calculateMaxChunkSize's
% dynamic available-RAM estimate, to keep this path's chunk sizing predictable).
targetBytes = 128 * 1024 * 1024; % 128 MB
slabIn = spatialSlabIO('open', dataFile, 'Info', meta);
cleanObj = onCleanup(@() spatialSlabIO('close', slabIn));

switch dimName

    case 'T'
        % slabP/reshape can coexist with slab during aggregation.
        bytesPerX = 2 * nY * nT * getByteSize('single');
        xPerSlab = max(1, floor(targetBytes / max(bytesPerX, 1)));

        aggData = zeros(nY, nX, 'single');
        xStart = 1;

        while xStart <= nX
            xEnd = min(nX, xStart + xPerSlab - 1);
            xIdx = xStart:xEnd;

            slab = single(spatialSlabIO('read', slabIn, xIdx));
            slabP = permute(slab, [3 1 2]);
            slabP = reshape(slabP, nT, []);
            aggFlat = iCalcAgg(slabP, aggFcn);
            aggData(:, xIdx) = reshape(single(aggFlat), nY, numel(xIdx));

            xStart = xEnd + 1;
        end

        outDimNames = {'Y','X'};

    case 'E'
        evObj = EventsManager(SaveFolder);
        frameRateHz = resolveDataInfoValue('frameRateHz', explicitRate, dataFile, 'apply_aggregate_function');
        [frMat, plan] = iEventFramePlan(evObj, nT, frameRateHz);

        nInst = size(frMat, 1);
        nCond = numel(plan.conditionID);
        trialLen = size(frMat, 2);

        % Live scratch: the input slab, all instances of the slab, and the
        % per-condition aggregate.
        bytesPerX = nY * (nT + 2 * trialLen * nInst + trialLen * nCond) * ...
            getByteSize('single');
        xPerSlab = max(1, floor(targetBytes / max(bytesPerX, 1)));

        aggData = zeros(nY, nX, trialLen, nCond, 'single');

        xStart = 1;
        while xStart <= nX
            xEnd = min(nX, xStart + xPerSlab - 1);
            xIdx = xStart:xEnd;

            slabData = single(spatialSlabIO('read', slabIn, xIdx));
            slabE = iInstancesFromFrames(slabData, frMat);
            aggData(:, xIdx, :, :) = EventsManager.reduceByCondition(slabE, plan, ...
                @(x) iCalcAgg(x, aggFcn, 4), 4);

            xStart = xEnd + 1;
        end

        eventInfo = plan.eventInfoOut;
        outDimNames = {'Y','X','T','E'};
end
end

% =========================================================================
% Helper: E aggregation of an event-split .dat (Y-X-T-E or Y-X-E)
% =========================================================================
function [outData, metaData] = iAggregateEventSplitDat(dataFile, datMeta, SaveFolder, aggFcn)
%IAGGREGATEEVENTSPLITDAT Per-condition aggregate of an event-split .dat.
%
% The file stores no event labels: its E axis is matched to events.mat by
% resolveDatEventMapping, then reduced through the UMT path in RAM. The
% output is numeric (saved as .dat), one E slice per condition.

assertDatLayout(datMeta, {{'Y','X','T','E'}, {'Y','X','E'}}, 'apply_aggregate_function');
dims = cellstr(string(datMeta.dimNames(:).'));
mapping = resolveDatEventMapping(datMeta, SaveFolder);
if strcmpi(mapping.status, 'aggregated')
    % Raises Umitoolbox:EventsManager:alreadyAggregated.
    EventsManager.conditionAggregationPlan(mapping.eventInfo);
end
if ~strcmpi(mapping.status, 'matched')
    warning('Umitoolbox:apply_aggregate_function:eventsNotMatched', '%s', mapping.message);
end

umt = genUMTStruct(single(loadData(dataFile)), 'kind', 'image', ...
    'entryName', 'main', 'dimNames', dims);
umt = appendUMTEventInfo(umt, 'eventInfo', mapping.eventInfo);
outUMT = apply_aggregate_function(umt, SaveFolder, 'aggregateFcn', aggFcn, ...
    'dimensionName', 'E');

outData = outUMT.data.main.value;
metaData = struct('dimNames', {dims});
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
