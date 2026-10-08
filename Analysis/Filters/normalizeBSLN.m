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
%   data       : .dat filename (or path) with axes Y-X-T (continuous) or
%                Y-X-T-E (event-split, e.g. split_data_by_event output).
%                Arrays, UMT structs, and .umt files are not supported.
%
%   SaveFolder : Folder of the .dat output; for trial mode it also holds
%                events.mat.
%
% Name-Value parameters:
%   normalizationMode : 'recording' or 'trial'
%                       Default: 'recording'
%                         - 'recording': the baseline is taken from the
%                           start of the recording (Y-X-T), or from the
%                           start of every trial (Y-X-T-E).
%                         - 'trial': Y-X-T-E only. The baseline of every
%                           trial is the baseline period of events.mat.
%
%   baselineMode      : 'auto' or positive numeric scalar (seconds)
%                       Default: 'auto'
%                       Notes:
%                         - recording mode:
%                             * 'auto'  => first 20%% of T
%                             * numeric => first baselineMode seconds
%                         - trial mode:
%                             * must be 'auto'
%                             * baseline uses the baseline period of
%                               events.mat
%
%   b_centerAtOne     : Logical scalar. If true, add 1 after DeltaR/R0.
%                       Default: false
%
%   FrameRateHz       : Frame rate of DATA (Hz), needed for a numeric
%                       baselineMode and for trial mode. PipelineManager
%                       injects it from the data flowing into the step; the
%                       .dat header provides it otherwise. AcqInfos.mat is
%                       not used (resolveDataInfoValue).
%
% Output:
%   outData           : Full path of the .dat output ("normBSLN.dat" in
%                       SaveFolder), with the same axes and sizes as the
%                       input (single precision). A Y-X-T-E file keeps its E
%                       axis: every trial is normalized on its own, ignored
%                       instances included, so the E axis still matches
%                       events.mat.
%
% Notes:
%   - The baseline of each pixel (and trial) is the median over the baseline
%     frames, omitting NaN; a zero baseline is replaced by 1.
%   - The file is streamed in X slabs sized from the available RAM and the
%     output is written slab by slab (Low-RAM mode is always on).
%   - Continuous data are not split here: trial mode needs an event-split
%     file, so split it first with split_data_by_event.

default_Output = 'normBSLN.dat';

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
    @(x) (((ischar(x) || (isstring(x) && isscalar(x))) && strcmpi(char(string(x)), 'auto')) || ...
         (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0)));
addParameter(p, 'b_centerAtOne', false, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'FrameRateHz', []);

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
normalizationMode = lower(char(string(p.Results.normalizationMode)));
baselineMode = p.Results.baselineMode;
b_centerAtOne = p.Results.b_centerAtOne;
explicitRate = p.Results.FrameRateHz;

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

if ~(ischar(data) || (isstring(data) && isscalar(data)))
    error('normalizeBSLN:UnsupportedInputType', ...
        'Input "data" must be the name or path of a .dat file (arrays and UMT inputs are not supported).');
end

dataFile = char(string(data));
if ~isfile(dataFile)
    altPath = fullfile(SaveFolder, dataFile);
    if isfile(altPath)
        dataFile = altPath;
    else
        error('normalizeBSLN:InputFileNotFound', ...
            'Input file "%s" was not found.', char(string(data)));
    end
end

[~, ~, ext] = fileparts(dataFile);
if ~strcmpi(ext, '.dat')
    error('normalizeBSLN:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', ext);
end

outData = iNormalizeDatFile(dataFile, SaveFolder, default_Output, normalizationMode, ...
    baselineMode, b_centerAtOne, explicitRate);

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            ['Normalize image data by baseline (DeltaR/R0): Y-X-T and event-split ' ...
             'Y-X-T-E .dat files in, a .dat file with the same axes out.']);
        info.version = '2.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            'Input .dat file with axes Y-X-T or Y-X-T-E.', ...
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
            'Folder receiving the .dat output and, for trial mode, containing events.mat.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'normalizationMode', ...
            'parameter', ...
            'Normalization mode: recording or trial (event-split data only).', ...
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
            'allowed', {'auto', [0 Inf]}, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'b_centerAtOne', ...
            'parameter', ...
            'If true, center normalized data at one.', ...
            'kind', 'parameter', ...
            'default', false, ...
            'allowed', [true false], ...
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
            {'ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Baseline-normalized .dat output with the same axes as the input.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: streamed normalization of a Y-X-T or Y-X-T-E .dat file
% =========================================================================
function outFile = iNormalizeDatFile(inFile, SaveFolder, defaultOutput, normalizationMode, baselineMode, b_centerAtOne, explicitRate)
%INORMALIZEDATFILE Read in X slabs, normalize along T per trial, write a .dat.

Info = loadMetaData(inFile);
assertDatLayout(Info, {{'Y','X','T'}, {'Y','X','T','E'}}, 'normalizeBSLN');
Ny = datAxisSize(Info, 'Y');
Nx = datAxisSize(Info, 'X');
Nt = datAxisSize(Info, 'T');
Ne = max(1, datAxisSize(Info, 'E'));
hasE = datAxisSize(Info, 'E') > 0;

isTrialMode = strcmp(normalizationMode, 'trial');
if isTrialMode && ~hasE
    error('normalizeBSLN:TrialModeRequiresEventSplit', ...
        ['Trial-mode normalization needs an event-split Y-X-T-E file ("%s" has ' ...
         'axes %s). Split it first with split_data_by_event.'], inFile, ...
        strjoin(cellstr(string(Info.dimNames(:).')), '-'));
end

% The frame rate matters only for a numeric baseline or for trial mode.
needsRate = isTrialMode || ~(ischar(baselineMode) || (isstring(baselineMode) && isscalar(baselineMode)));
freqHz = NaN;
if needsRate
    freqHz = resolveDataInfoValue('frameRateHz', explicitRate, inFile, 'normalizeBSLN');
end

if isTrialMode
    mapping = resolveDatEventMapping(Info, SaveFolder);
    if ~(isfield(mapping.eventInfo, 'baselinePeriod') && ~isempty(mapping.eventInfo.baselinePeriod))
        error('normalizeBSLN:MissingBaselinePeriod', ...
            'Trial-mode normalization needs a baseline period in the events.mat of "%s".', SaveFolder);
    end
    nBaseFrames = iResolveTrialBaselineFrames(Nt, freqHz, double(mapping.eventInfo.baselinePeriod));
else
    nBaseFrames = iResolveRecordingBaselineFrames(Nt, freqHz, baselineMode);
end

slabIn = spatialSlabIO('open', inFile, 'Info', Info);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

% Write through a scratch file so the output only appears once the run has
% completed, and so the input can be the very file the output replaces.
outFile = fullfile(SaveFolder, defaultOutput);
[~, outStem, outExt] = fileparts(outFile);
tmpFile = fullfile(SaveFolder, [outStem, '_writing', outExt]);
cTmp = onCleanup(@() iDeleteIfExists(tmpFile));
slabOut = spatialSlabIO('create', tmpFile, ...
    datHeaderFromInfo(Info, outStem, 'dataClass', 'single'));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

% Peak memory of a slab: the slab, the normalized copy, and a temporary.
nChunks = calculateMaxChunkSize(double(Ny) * Nx * Nt * Ne * 4, 3, 0.2);
chunkX = max(1, ceil(Nx / nChunks));
nChunks = ceil(Nx / chunkX);

for c = 1:nChunks
    xIdx = ((c-1) * chunkX + 1):min(c * chunkX, Nx);

    fprintf('Chunk %i/%i [Reading file ...]\n', c, nChunks)
    slab = single(spatialSlabIO('read', slabIn, xIdx));

    fprintf('Chunk %i/%i [Normalizing data ...]\n', c, nChunks)
    % The baseline of every pixel and trial: median over the first frames.
    bsln = median(slab(:, :, 1:nBaseFrames, :), 3, 'omitnan');
    bsln(bsln == 0) = 1;
    slab = (slab - bsln) ./ bsln;
    if b_centerAtOne
        slab = slab + 1;
    end

    fprintf('Chunk %i/%i [Writing to file ...]\n', c, nChunks)
    spatialSlabIO('write', slabOut, xIdx, slab);
    fprintf('Chunk %i/%i [Completed]\n', c, nChunks)
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close both files before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'normalizeBSLN:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end

function iDeleteIfExists(filePath)
%IDELETEIFEXISTS Remove a scratch file left by a failed run.
if isfile(filePath)
    delete(filePath);
end
end

% =========================================================================
% Helper: Resolve baseline frames for recording mode
% =========================================================================
function nBaseFrames = iResolveRecordingBaselineFrames(nT, freqHz, baselineMode)

if ischar(baselineMode) || (isstring(baselineMode) && isscalar(baselineMode))
    assert(strcmpi(char(string(baselineMode)), 'auto'), ...
        'normalizeBSLN:InvalidBaselineMode', ...
        'baselineMode must be ''auto'' or a positive numeric scalar, got "%s".', ...
        char(string(baselineMode)));
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
