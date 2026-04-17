function outData = normalizeBSLN(data, SaveFolder, varargin)
%NORMALIZEBSLN Normalize image data by baseline (DeltaR/R0).
%
%   outData = normalizeBSLN(data, SaveFolder)
%   outData = normalizeBSLN(data, SaveFolder, ...
%       'normalizationMode', mode, ...
%       'baselineMode', baselineMode, ...
%       'b_centerAtOne', tf)
%
% Inputs:
%   data       : One of:
%                1) Numeric 3-D array with dimensions Y x X x T
%                2) Filename to a .dat file storing Y x X x T data
%                3) UMT struct
%                4) Filename to a .umt or .mat file containing a UMT struct
%
%   SaveFolder : Folder containing AcqInfos.mat and, for trial mode,
%                events.mat.
%
% Name-Value parameters:
%   normalizationMode : 'recording' or 'trial'
%                       Default: 'recording'
%
%   baselineMode      : 'auto' or positive numeric scalar (seconds)
%                       Default: 'auto'
%                       Notes:
%                         - recording mode:
%                             * 'auto'  => first 20%% of T
%                             * numeric => first baselineMode seconds
%                         - trial mode:
%                             * must be 'auto'
%                             * baseline uses EventsManager.baselinePeriod
%
%   b_centerAtOne     : Logical scalar. If true, add 1 after DeltaR/R0.
%                       Default: false
%
% Output:
%   outData           : Output UMT struct.
%
% Notes:
%   - Raw YXT arrays and raw .dat files use EventsManager only when
%     normalizationMode = 'trial'.
%   - UMT input must have kind = 'image'.
%   - All UMT entries must contain Y, X, and T.
%   - In recording mode, UMT entries must not contain E.
%   - In trial mode, all UMT entries must contain E, and the UMT must
%     provide shared top-level eventInfo with eventAxisMode = 'instances'.
%   - RAM-safe mode is only implemented for raw .dat files.
%   - If a .umt or .mat file is provided, the UMT content is loaded
%     into RAM and processed there.

default_Output = 'normBSLN.umt'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'normalizationMode', 'recording', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'baselineMode', 'auto', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) || ...
         (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
addParameter(p, 'b_centerAtOne', false, ...
    @(x) islogical(x) && isscalar(x));

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
normalizationMode = lower(char(string(p.Results.normalizationMode)));
baselineMode = p.Results.baselineMode;
b_centerAtOne = p.Results.b_centerAtOne;

if ~ismember(normalizationMode, {'recording','trial'})
    error('normalizeBSLN:InvalidNormalizationMode', ...
        'normalizationMode must be ''recording'' or ''trial''.');
end

if strcmpi(normalizationMode, 'trial')
    if ~(ischar(baselineMode) || (isstring(baselineMode) && isscalar(baselineMode))) || ...
            ~strcmpi(char(string(baselineMode)), 'auto')
        error('normalizeBSLN:InvalidBaselineModeForTrial', ...
            'baselineMode must be ''auto'' when normalizationMode = ''trial''.');
    end
end

if ~isfolder(SaveFolder)
    error('normalizeBSLN:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Case 1: Raw YXT array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)

    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, ...
        mfilename, 'data');

    rawData = single(data);
    freqHz = iGetFrameRateHz(SaveFolder);

    switch normalizationMode
        case 'recording'
            nT = size(rawData, 3);
            nBaseFrames = iResolveRecordingBaselineFrames(nT, freqHz, baselineMode);

            bsln = median(rawData(:,:,1:nBaseFrames), 3, 'omitnan');
            outVal = (rawData - bsln) ./ bsln;
            if b_centerAtOne
                outVal = outVal + 1;
            end

            outData = iPackageOutputUMT( ...
                {'main'}, ...
                {single(outVal)}, ...
                {{'Y','X','T'}}, ...
                struct(), ...
                struct());

        case 'trial'
            evObj = EventsManager(SaveFolder, '', 'csv');
            [frMat, conditionIDlist, repetitionList] = evObj.getFrameMatrix(size(rawData, 3));

            if isempty(frMat)
                error('normalizeBSLN:NoEventsFound', ...
                    'No valid events were found for trial normalization.');
            end

            trialLen = size(frMat, 2);
            nBaseFrames = iResolveTrialBaselineFrames(trialLen, freqHz, double(evObj.baselinePeriod));

            nTrials = size(frMat, 1);
            nY = size(rawData, 1);
            nX = size(rawData, 2);

            outVal = nan(nY, nX, trialLen, nTrials, 'single');
            eventNames = cell(nTrials, 1);

            for iTrial = 1:nTrials
                validMask = ~isnan(frMat(iTrial, :));
                frameIdx = frMat(iTrial, validMask);

                if isempty(frameIdx)
                    continue
                end

                trialData = nan(nY, nX, trialLen, 'single');
                trialData(:,:,validMask) = rawData(:,:,frameIdx);

                bsln = median(trialData(:,:,1:nBaseFrames), 3, 'omitnan');
                trialData = (trialData - bsln) ./ bsln;

                if b_centerAtOne
                    trialData = trialData + 1;
                end

                outVal(:,:,:,iTrial) = trialData;
                eventNames{iTrial} = evObj.eventNameList{conditionIDlist(iTrial)};
            end

            labels = struct();
            labels.E = eventNames;

            eventInfo = struct();
            eventInfo.eventID = conditionIDlist(:);
            eventInfo.repetitionIndex = repetitionList(:);
            eventInfo.eventName = eventNames;
            eventInfo.eventAxisMode = 'instances';

            outData = iPackageOutputUMT( ...
                {'main'}, ...
                {outVal}, ...
                {{'Y','X','T','E'}}, ...
                labels, ...
                eventInfo);
    end

    return
end

% -------------------------------------------------------------------------
% Case 2: File input
% -------------------------------------------------------------------------
if ischar(data) || (isstring(data) && isscalar(data))

    dataFile = char(string(data));

    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('normalizeBSLN:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);

    switch ext
        case '.dat'
            [outVal, outDimNames, labels, eventInfo] = ...
                normalizeBSLN_lowRAMmode( ...
                    dataFile, SaveFolder, normalizationMode, baselineMode, b_centerAtOne);

            outData = iPackageOutputUMT( ...
                {'main'}, ...
                {outVal}, ...
                {outDimNames}, ...
                labels, ...
                eventInfo);
            return

        case {'.umt','.mat'}
            warning('normalizeBSLN:UMTFileLoadsInRAM', ...
                ['RAM-safe mode is not available for data stored in this format. ' ...
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
                    error('normalizeBSLN:NoUMTFoundInFile', ...
                        'No scalar UMT struct was found in "%s".', dataFile);
                end
            end

            data = loadedUMT;

        otherwise
            error('normalizeBSLN:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

% -------------------------------------------------------------------------
% Case 3: UMT struct in RAM
% -------------------------------------------------------------------------
if ~isstruct(data)
    error('normalizeBSLN:UnsupportedInputType', ...
        ['Input "data" must be a YXT array, a .dat filename, ' ...
         'a UMT struct, or a .umt/.mat filename containing a UMT struct.']);
end

[entryNames, entryData, entryDims, sourceLabels, sourceEventInfo, hasE] = ...
    iExtractValidUMTData(data);

switch normalizationMode
    case 'recording'
        if any(hasE)
            error('normalizeBSLN:RecordingModeRequiresContinuousUMT', ...
                ['Recording-mode normalization only supports continuous UMT entries ' ...
                 'without an E dimension.']);
        end

        freqHz = iGetFrameRateHz(SaveFolder);
        baselineSec = [];
        if ~(ischar(baselineMode) || (isstring(baselineMode) && isscalar(baselineMode)))
            baselineSec = double(baselineMode);
        end

        outEntryData = entryData;
        outEntryDims = entryDims;

        for iEntry = 1:numel(entryNames)
            thisData = single(entryData{iEntry});
            thisDims = entryDims{iEntry};

            idxT = find(strcmp(thisDims, 'T'), 1, 'first');

            permOrder = [setdiff(1:ndims(thisData), idxT, 'stable') idxT];
            dataP = permute(thisData, permOrder);
            szP = size(dataP);
            nT = szP(end);

            if isempty(baselineSec)
                nBaseFrames = iResolveRecordingBaselineFrames(nT, freqHz, 'auto');
            else
                nBaseFrames = iResolveRecordingBaselineFrames(nT, freqHz, baselineSec);
            end

            spatialShape = szP(1:end-1);
            data2D = reshape(dataP, prod(spatialShape), nT);

            bsln = median(data2D(:, 1:nBaseFrames), 2, 'omitnan');
            data2D = (data2D - bsln) ./ bsln;

            if b_centerAtOne
                data2D = data2D + 1;
            end

            dataP = reshape(data2D, szP);
            outEntryData{iEntry} = ipermute(dataP, permOrder);
            outEntryDims{iEntry} = thisDims;
        end

        outData = iPackageOutputUMT( ...
            entryNames, outEntryData, outEntryDims, sourceLabels, struct());

    case 'trial'
        if ~all(hasE)
            error('normalizeBSLN:TrialModeRequiresEventSplitUMT', ...
                ['Trial-mode normalization on UMT input requires all entries to ' ...
                 'contain an E dimension.']);
        end

        if isempty(fieldnames(sourceEventInfo))
            error('normalizeBSLN:MissingEventInfo', ...
                'Trial-mode normalization on UMT input requires shared top-level eventInfo.');
        end

        if ~strcmpi(sourceEventInfo.eventAxisMode, 'instances')
            error('normalizeBSLN:InvalidEventAxisMode', ...
                ['Trial-mode normalization on UMT input requires eventAxisMode = ' ...
                 '"instances".']);
        end

        evObj = EventsManager(SaveFolder, '', 'csv');
        freqHz = iGetFrameRateHz(SaveFolder);

        outEntryData = entryData;
        outEntryDims = entryDims;

        for iEntry = 1:numel(entryNames)
            thisData = single(entryData{iEntry});
            thisDims = entryDims{iEntry};

            idxT = find(strcmp(thisDims, 'T'), 1, 'first');
            idxE = find(strcmp(thisDims, 'E'), 1, 'first');

            permOrder = [setdiff(1:ndims(thisData), [idxT idxE], 'stable') idxT idxE];
            dataP = permute(thisData, permOrder);
            szP = size(dataP);

            nT = szP(end-1);
            nTrials = szP(end);
            nBaseFrames = iResolveTrialBaselineFrames(nT, freqHz, double(evObj.baselinePeriod));

            spatialShape = szP(1:end-2);
            data2D = reshape(dataP, prod(spatialShape), nT, nTrials);

            for iTrial = 1:nTrials
                bsln = median(data2D(:, 1:nBaseFrames, iTrial), 2, 'omitnan');
                data2D(:, :, iTrial) = (data2D(:, :, iTrial) - bsln) ./ bsln;

                if b_centerAtOne
                    data2D(:, :, iTrial) = data2D(:, :, iTrial) + 1;
                end
            end

            dataP = reshape(data2D, szP);
            outEntryData{iEntry} = ipermute(dataP, permOrder);
            outEntryDims{iEntry} = thisDims;
        end

        outLabels = sourceLabels;
        if isempty(fieldnames(outLabels))
            outLabels = struct();
        end
        if ~isfield(outLabels, 'E')
            outLabels.E = sourceEventInfo.eventName(:).';
        end

        outData = iPackageOutputUMT( ...
            entryNames, outEntryData, outEntryDims, outLabels, sourceEventInfo);
end

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            ['Normalize raw image time-series or processed event-split data ' ...
             'by baseline and return a UMT struct.']);
        

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT array, .dat filename, ' ...
             'UMT struct, or .umt/.mat file containing one UMT struct.'], ...
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
            'Folder containing AcqInfos.mat and, for trial mode, events.mat.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'normalizationMode', ...
            'parameter', ...
            'Normalization mode: recording or trial.', ...
            'kind', 'parameter', ...
            'default', 'recording', ...
            'allowed', {'recording','trial'}, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'baselineMode', ...
            'parameter', ...
            'Baseline period mode: auto or numeric seconds.', ...
            'kind', 'parameter', ...
            'default', 'auto', ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'b_centerAtOne', ...
            'parameter', ...
            'If true, center normalized data at one.', ...
            'kind', 'parameter', ...
            'default', false, ...
            'allowed', {true,false}, ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ProcessedData','UnknownDataType'}, ...
            'data', ...
            'Baseline-normalized output UMT.', ...
            'normBSLN.umt', ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: Low-RAM execution for raw .dat input
% =========================================================================
function [outVal, outDimNames, labels, eventInfo] = normalizeBSLN_lowRAMmode( ...
    inFile, SaveFolder, normalizationMode, baselineMode, b_centerAtOne)

labels = struct();
eventInfo = struct();

[Ny, Nx, Nt, freqHz] = iGetRawDatInfo(SaveFolder, inFile);

fidIn = fopen(inFile, 'r');
if fidIn == -1
    error('normalizeBSLN:FileOpenFailed', ...
        'Failed to open "%s".', inFile);
end
cleanObj = onCleanup(@() fclose(fidIn)); %#ok<NASGU>

targetBytes = 128 * 1024 * 1024; % 128 MB
bytesPerX = Ny * Nt * 4;
xPerSlab = max(1, floor(targetBytes / max(bytesPerX, 1)));

switch lower(normalizationMode)

    case 'recording'
        nBaseFrames = iResolveRecordingBaselineFrames(Nt, freqHz, baselineMode);

        outVal = zeros(Ny, Nx, Nt, 'single');
        xStart = 1;

        while xStart <= Nx
            xEnd = min(Nx, xStart + xPerSlab - 1);
            xIdx = xStart:xEnd;

            slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, 'single');

            bsln = median(slab(:,:,1:nBaseFrames), 3, 'omitnan');
            slab = (slab - bsln) ./ bsln;

            if b_centerAtOne
                slab = slab + 1;
            end

            outVal(:, xIdx, :) = slab;
            xStart = xEnd + 1;
        end

        outDimNames = {'Y','X','T'};

    case 'trial'
        evObj = EventsManager(SaveFolder, '', 'csv');
        [frMat, conditionIDlist, repetitionList] = evObj.getFrameMatrix(Nt);

        if isempty(frMat)
            error('normalizeBSLN:NoEventsFound', ...
                'No valid events were found for trial normalization.');
        end

        nTrials = size(frMat, 1);
        trialLen = size(frMat, 2);
        nBaseFrames = iResolveTrialBaselineFrames(trialLen, freqHz, double(evObj.baselinePeriod));

        outVal = nan(Ny, Nx, trialLen, nTrials, 'single');
        eventNames = cell(nTrials, 1);

        xStart = 1;
        while xStart <= Nx
            xEnd = min(Nx, xStart + xPerSlab - 1);
            xIdx = xStart:xEnd;

            slabData = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, 'single');

            for iTrial = 1:nTrials
                validMask = ~isnan(frMat(iTrial, :));
                frameIdx = frMat(iTrial, validMask);

                if isempty(frameIdx)
                    continue
                end

                trialData = nan(Ny, numel(xIdx), trialLen, 'single');
                trialData(:,:,validMask) = slabData(:,:,frameIdx);

                bsln = median(trialData(:,:,1:nBaseFrames), 3, 'omitnan');
                trialData = (trialData - bsln) ./ bsln;

                if b_centerAtOne
                    trialData = trialData + 1;
                end

                outVal(:, xIdx, :, iTrial) = trialData;
                eventNames{iTrial} = evObj.eventNameList{conditionIDlist(iTrial)};
            end

            xStart = xEnd + 1;
        end

        labels.E = eventNames;

        eventInfo.eventID = conditionIDlist(:);
        eventInfo.repetitionIndex = repetitionList(:);
        eventInfo.eventName = eventNames;
        eventInfo.eventAxisMode = 'instances';

        outDimNames = {'Y','X','T','E'};
end
end

% =========================================================================
% Helper: Extract and validate image-backed data from a UMT structure
% =========================================================================
function [entryNames, entryData, entryDims, labels, eventInfo, hasE] = iExtractValidUMTData(umt)

validateUMTStruct(umt, 'requireEventInfo', false);

if ~strcmpi(umt.kind, 'image')
    error('normalizeBSLN:InvalidUMTKind', ...
        ['Operation aborted. UMT input must have kind = "image". ' ...
         'This function does not support non-image UMT structures.']);
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error('normalizeBSLN:EmptyUMTData', ...
        'Operation aborted. UMT data is empty.');
end

entryData = cell(size(entryNames));
entryDims = cell(size(entryNames));
hasE = false(size(entryNames));

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    if ~all(ismember({'Y','X','T'}, thisDims))
        error('normalizeBSLN:InvalidUMTEntry', ...
            ['Operation aborted. All entries in the input UMT must be image-backed ' ...
             'and contain dimensions Y, X, and T.\nInvalid entry: "%s".'], ...
            entryNames{iEntry});
    end

    entryData{iEntry} = single(thisEntry.value);
    entryDims{iEntry} = thisDims;
    hasE(iEntry) = any(strcmp(thisDims, 'E'));
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
% Helper: Package processed entries into a UMT output
% =========================================================================
function outUMT = iPackageOutputUMT(entryNames, entryData, entryDims, labelsIn, eventInfoIn)

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
                'dimNames', entryDims{iEntry});
        else
            outUMT = genUMTStruct( ...
                entryData{iEntry}, ...
                'kind', 'image', ...
                'entryName', entryNames{iEntry}, ...
                'dimNames', entryDims{iEntry}, ...
                'labels', labelsOut);
        end
    else
        outUMT = genUMTStruct( ...
            outUMT, ...
            'value', entryData{iEntry}, ...
            'entryName', entryNames{iEntry}, ...
            'dimNames', entryDims{iEntry});
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
% Helper: Acquisition info for raw .dat
% =========================================================================
function [Ny, Nx, Nt, freqHz] = iGetRawDatInfo(SaveFolder, inFile)

acqFile = fullfile(SaveFolder, 'AcqInfos.mat');
if ~isfile(acqFile)
    error('normalizeBSLN:MissingAcqInfos', ...
        'AcqInfos.mat was not found in SaveFolder "%s".', SaveFolder);
end

S = load(acqFile);
if isfield(S, 'AcqInfoStream')
    acqInfo = S.AcqInfoStream;
else
    fn = fieldnames(S);
    acqInfo = S.(fn{1});
end

if ~isfield(acqInfo, 'Height') || ~isfield(acqInfo, 'Width') || ~isfield(acqInfo, 'FrameRateHz')
    error('normalizeBSLN:InvalidAcqInfoStream', ...
        'AcqInfoStream must contain Height, Width, and FrameRateHz.');
end

Ny = double(acqInfo.Height);
Nx = double(acqInfo.Width);
freqHz = double(acqInfo.FrameRateHz);

if isfield(acqInfo, 'Length') && ~isempty(acqInfo.Length)
    Nt = double(acqInfo.Length);
else
    if isempty(inFile)
        error('normalizeBSLN:MissingLength', ...
            'Failed to determine data length from AcqInfos.mat.');
    end
    info = dir(inFile);
    nElem = info.bytes / 4; % single precision
    if mod(nElem, Ny*Nx) ~= 0
        error('normalizeBSLN:IncompatibleRawDatSize', ...
            'Raw .dat file size is incompatible with the inferred image size.');
    end
    Nt = nElem / (Ny*Nx);
end
end

% =========================================================================
% Helper: Frame rate
% =========================================================================
function freqHz = iGetFrameRateHz(SaveFolder)

acqFile = fullfile(SaveFolder, 'AcqInfos.mat');
if ~isfile(acqFile)
    error('normalizeBSLN:MissingAcqInfos', ...
        'AcqInfos.mat was not found in SaveFolder "%s".', SaveFolder);
end

S = load(acqFile);
if isfield(S, 'AcqInfoStream')
    acqInfo = S.AcqInfoStream;
else
    fn = fieldnames(S);
    acqInfo = S.(fn{1});
end

if ~isfield(acqInfo, 'FrameRateHz') || isempty(acqInfo.FrameRateHz)
    error('normalizeBSLN:MissingFrameRate', ...
        'AcqInfoStream must contain FrameRateHz.');
end

freqHz = double(acqInfo.FrameRateHz);
end

% =========================================================================
% Helper: Resolve baseline frames for recording mode
% =========================================================================
function nBaseFrames = iResolveRecordingBaselineFrames(nT, freqHz, baselineMode)

if ischar(baselineMode) || (isstring(baselineMode) && isscalar(baselineMode))
    nBaseFrames = round(0.2 * nT);
else
    nBaseFrames = round(double(baselineMode) * freqHz);
end

nBaseFrames = max(1, nBaseFrames);
nBaseFrames = min(nBaseFrames, nT);
end

% =========================================================================
% Helper: Resolve baseline frames for trial mode
% =========================================================================
function nBaseFrames = iResolveTrialBaselineFrames(trialLen, freqHz, baselinePeriodSec)

nBaseFrames = round(double(baselinePeriodSec) * freqHz);
nBaseFrames = max(1, nBaseFrames);
nBaseFrames = min(nBaseFrames, trialLen);
end