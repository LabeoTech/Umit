function outData = apply_aggregate_function(data, SaveFolder, varargin)
%APPLY_AGGREGATE_FUNCTION Aggregate image-backed data along T or E.
%
%   outData = apply_aggregate_function(data, SaveFolder)
%   outData = apply_aggregate_function(data, SaveFolder, ...
%       'aggregateFcn', aggFcn, 'dimensionName', dimName)
%
% Inputs:
%   data       : One of:
%                1) Numeric 3-D array with dimensions Y x X x T
%                2) Filename to a .dat file storing Y x X x T data
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
%
% Output:
%   outData    : Output UMT struct.
%
% Notes:
%   - Raw YXT arrays and raw .dat files use live event information from
%     events.mat through EventsManager.
%   - UMT inputs use the frozen shared top-level eventInfo stored in the UMT.
%   - RAM-safe mode is only implemented for raw .dat input.
%   - If a .umt file is provided, RAM-safe mode is not available
%     and the UMT content is loaded into RAM.
%   - For E-dimension aggregation, RAM-safe mode chunks the input read but
%     still allocates the full Y x X x trialLen x nCond output array before
%     writing it out, so peak RAM is not bounded by condition count
%     (DFR-20260819-012).
%   - Output eventInfo.eventAxisMode is always 'aggregated_repetitions' for
%     E-dimension aggregation: repetitions have been collapsed by aggFcn, so
%     eventInfo.repetitionIndex is a fixed all-zero sentinel, not a real
%     per-entry repetition count.
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

parse(p, data, SaveFolder, varargin{:});

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
            [frMat, conditionIDlist, ~] = evObj.getFrameMatrix(size(rawData, 3));

            if isempty(frMat)
                error('apply_aggregate_function:NoEventsFound', ...
                    'No valid events were found in events.mat for E aggregation.');
            end

            condIDs = unique(conditionIDlist(:), 'stable');
            nCond = numel(condIDs);
            nY = size(rawData, 1);
            nX = size(rawData, 2);
            trialLen = size(frMat, 2);

            aggData = zeros(nY, nX, trialLen, nCond, 'single');
            eventNames = cell(nCond, 1);

            for iCond = 1:nCond
                rowIdx = find(conditionIDlist(:) == condIDs(iCond));
                nRep = numel(rowIdx);

                condBlock = nan(nY, nX, trialLen, nRep, 'single');

                for iRep = 1:nRep
                    validMask = ~isnan(frMat(rowIdx(iRep), :));
                    if any(validMask)
                        frameIdx = frMat(rowIdx(iRep), validMask);
                        condBlock(:, :, validMask, iRep) = rawData(:, :, frameIdx);
                    end
                end

                condP = permute(condBlock, [4 1 2 3]);
                condP = reshape(condP, nRep, []);
                aggFlat = iCalcAgg(condP, aggFcn);
                aggData(:, :, :, iCond) = reshape(single(aggFlat), nY, nX, trialLen);

                eventNames{iCond} = evObj.eventNameList{condIDs(iCond)};
            end

            labels = struct();

            eventInfo = struct();
            eventInfo.eventID = condIDs(:);
            eventInfo.repetitionIndex = zeros(nCond, 1);
            eventInfo.eventName = eventNames;
            eventInfo.eventAxisMode = 'aggregated_repetitions';

            outData = iPackageOutputUMT( ...
                {'main'}, ...
                {aggData}, ...
                {{'Y','X','T','E'}}, ...
                labels, ...
                eventInfo, ...
                {struct()});
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
            [aggData, outDimNames, labels, eventInfo] = ...
                iExecuteLowRAMDat(dataFile, SaveFolder, aggFcn, dimName);

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

    condIDs = unique(sourceEventInfo.eventID(:), 'stable');
    nCond = numel(condIDs);
    eventNames = cell(nCond, 1);

    for iCond = 1:nCond
        idxFirst = find(sourceEventInfo.eventID(:) == condIDs(iCond), 1, 'first');
        eventNames{iCond} = sourceEventInfo.eventName{idxFirst};
    end
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
            permOrder = [idxTarget setdiff(1:ndims(thisData), idxTarget)];
            dataP = permute(thisData, permOrder);
            szP = size(dataP);
            assert(numel(sourceEventInfo.eventID) == szP(1), ...
                'apply_aggregate_function:EventAxisMismatch', ...
                ['Entry "%s" has %d elements along dimension "E", but ' ...
                 'sourceEventInfo.eventID has %d elements.'], ...
                entryNames{iEntry}, szP(1), numel(sourceEventInfo.eventID));
            dataP = reshape(dataP, szP(1), []);

            outFlat = zeros(nCond, size(dataP, 2), 'single');
            for iCond = 1:nCond
                idxCond = sourceEventInfo.eventID(:) == condIDs(iCond);
                outFlat(iCond, :) = single(iCalcAgg(dataP(idxCond, :), aggFcn));
            end

            outP = reshape(outFlat, [nCond szP(2:end)]);
            outEntryData{iEntry} = ipermute(outP, permOrder);
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
        outEventInfo = struct();
        outEventInfo.eventID = condIDs(:);
        outEventInfo.repetitionIndex = zeros(nCond, 1);
        outEventInfo.eventName = eventNames;
        outEventInfo.eventAxisMode = 'aggregated_repetitions';
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

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'Aggregated UMT output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: Low-RAM execution for raw .dat input
% =========================================================================
function [aggData, outDimNames, labels, eventInfo] = iExecuteLowRAMDat(dataFile, SaveFolder, aggFcn, dimName)

labels = struct();
eventInfo = struct();

meta = loadMetaData(dataFile);
if ~isfield(meta, 'Height') || ~isfield(meta, 'Width') || ~isfield(meta, 'datLength')
    error('apply_aggregate_function:InvalidMetaData', ...
        'loadMetaData did not return Height, Width, and datLength for "%s".', dataFile);
end

nY = double(meta.Height);
nX = double(meta.Width);
nT = double(meta.datLength);

% Conservative fixed chunk-size budget (not derived from calculateMaxChunkSize's
% dynamic available-RAM estimate, to keep this path's chunk sizing predictable).
targetBytes = 128 * 1024 * 1024; % 128 MB
fid = fopen(dataFile, 'r');
if fid == -1
    error('apply_aggregate_function:FileOpenFailed', ...
        'Failed to open "%s".', dataFile);
end
cleanObj = onCleanup(@() safeFclose(fid));

switch dimName

    case 'T'
        bytesPerX = nY * nT * getByteSize('single');
        xPerSlab = max(1, floor(targetBytes / max(bytesPerX, 1)));

        aggData = zeros(nY, nX, 'single');
        xStart = 1;

        while xStart <= nX
            xEnd = min(nX, xStart + xPerSlab - 1);
            xIdx = xStart:xEnd;

            slab = spatialSlabIO('read', fid, nY, nX, nT, xIdx, 'single');
            slabP = permute(slab, [3 1 2]);
            slabP = reshape(slabP, nT, []);
            aggFlat = iCalcAgg(slabP, aggFcn);
            aggData(:, xIdx) = reshape(single(aggFlat), nY, numel(xIdx));

            xStart = xEnd + 1;
        end

        outDimNames = {'Y','X'};

    case 'E'
        evObj = EventsManager(SaveFolder);
        [frMat, conditionIDlist, ~] = evObj.getFrameMatrix(nT);

        if isempty(frMat)
            error('apply_aggregate_function:NoEventsFound', ...
                'No valid events were found in events.mat for E aggregation.');
        end

        condIDs = unique(conditionIDlist(:), 'stable');
        nCond = numel(condIDs);
        trialLen = size(frMat, 2);

        maxReps = 0;
        for iCond = 1:nCond
            maxReps = max(maxReps, sum(conditionIDlist(:) == condIDs(iCond)));
        end

        bytesPerX = nY * max(nT, trialLen * max(maxReps,1)) * getByteSize('single');
        xPerSlab = max(1, floor(targetBytes / max(bytesPerX, 1)));

        aggData = zeros(nY, nX, trialLen, nCond, 'single');
        eventNames = cell(nCond, 1);

        xStart = 1;
        while xStart <= nX
            xEnd = min(nX, xStart + xPerSlab - 1);
            xIdx = xStart:xEnd;

            slabData = spatialSlabIO('read', fid, nY, nX, nT, xIdx, 'single');

            for iCond = 1:nCond
                rowIdx = find(conditionIDlist(:) == condIDs(iCond));
                nRep = numel(rowIdx);

                condBlock = nan(nY, numel(xIdx), trialLen, nRep, 'single');
                eventNames{iCond} = evObj.eventNameList{condIDs(iCond)};

                for iRep = 1:nRep
                    validMask = ~isnan(frMat(rowIdx(iRep), :));
                    if any(validMask)
                        frameIdx = frMat(rowIdx(iRep), validMask);
                        condBlock(:, :, validMask, iRep) = slabData(:, :, frameIdx);
                    end
                end

                condP = permute(condBlock, [4 1 2 3]);
                condP = reshape(condP, nRep, []);
                aggFlat = iCalcAgg(condP, aggFcn);
                aggData(:, xIdx, :, iCond) = reshape(single(aggFlat), nY, numel(xIdx), trialLen);
            end

            xStart = xEnd + 1;
        end

        eventInfo.eventID = condIDs(:);
        eventInfo.repetitionIndex = zeros(nCond, 1);
        eventInfo.eventName = eventNames;
        eventInfo.eventAxisMode = 'aggregated_repetitions';

        outDimNames = {'Y','X','T','E'};
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
    outUMT = appendUMTEventInfo(outUMT, ...
        'eventID', eventInfoIn.eventID, ...
        'repetitionIndex', eventInfoIn.repetitionIndex, ...
        'eventName', eventInfoIn.eventName, ...
        'eventAxisMode', eventInfoIn.eventAxisMode, ...
        'overwrite', true);
else
    validateUMTStruct(outUMT, 'requireEventInfo', true);
end
end

% =========================================================================
% Helper: Core aggregation function
% =========================================================================
function out = iCalcAgg(vals, aggFcn)
%ICALCAGG Reduce VALS along dimension 1 using AGGFCN, omitting NaN.
%
% For a pixel/column that is entirely NaN, 'sum' returns 0 (sum's own
% 'omitnan' convention: an empty/all-NaN input sums to 0), while 'mean',
% 'median', 'std', 'max', and 'min' all return NaN. A fully-masked pixel
% therefore reads as a valid zero under 'sum' aggregation, not as missing
% data the way it does under every other aggregateFcn.

switch lower(aggFcn)
    case 'mean'
        out = mean(vals, 1, 'omitnan');

    case 'median'
        out = median(vals, 1, 'omitnan');

    case 'std'
        out = std(vals, 0, 1, 'omitnan');

    case 'max'
        out = max(vals, [], 1, 'omitnan');

    case 'min'
        out = min(vals, [], 1, 'omitnan');

    case 'sum'
        out = sum(vals, 1, 'omitnan');

    otherwise
        error('apply_aggregate_function:InvalidAggregateFcnInternal', ...
            'Unsupported aggregateFcn "%s".', aggFcn);
end
end