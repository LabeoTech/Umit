function outData = normalizeLPF(data, SaveFolder, varargin)
%NORMALIZELPF Normalize image data by low-pass filtering.
%
%   outData = normalizeLPF(data, SaveFolder)
%   outData = normalizeLPF(data, SaveFolder, ...
%       'BaselineCutoffHz', baselineCutoffHz, ...
%       'SignalCutoffHz', signalCutoffHz, ...
%       'Normalize', tfNormalize, ...
%       'bApplyExpFit', tfApplyExpFit)
%
%   The filtering algorithm consists in creating two low-passed versions of
%   the signal with cut-off frequencies "BaselineCutoffHz" and
%   "SignalCutoffHz", then subtracting them. Optionally, the result can be
%   normalized by the baseline component to express the signal as DeltaR/R.
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X x T
%       2) Numeric array with dimensions Y x X x T x E (event-split data)
%       3) Filename to a .dat file with axes Y-X-T or Y-X-T-E
%       UMT structs and .umt/.mat files are not supported.
%
%   Input/output behavior:
%       - If the input is a numeric array, the output is a numeric array
%         with the same size.
%       - If the input is a .dat filename, the output is a .dat filename
%         ("normLPF.dat" in SaveFolder) with the same axes and sizes.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Output folder (for .dat inputs).
%
%   Name-Value parameters:
%       BaselineCutoffHz - Low cut-off frequency used to estimate the slow
%                          baseline component. 0 disables the baseline
%                          filter. When the recording is shorter than
%                          1/BaselineCutoffHz, the cutoff is raised to
%                          1/(2 x duration), with a warning.
%                          Default: 0.0083
%
%       SignalCutoffHz   - Higher cut-off frequency used to preserve the
%                          signal component. 0 disables the signal
%                          low-pass (the signal is kept unfiltered).
%                          Default: 1
%
%       Normalize        - Logical scalar. If true, express the filtered
%                          signal as DeltaR/R. Default: true
%
%       bApplyExpFit     - Logical scalar. If true, apply a double
%                          exponential decay correction before filtering.
%                          Default: false
%
%       FrameRateHz      - Frame rate of DATA (Hz), used for the filter.
%                          PipelineManager injects it from the data; a
%                          .dat input's header provides it otherwise.
%                          In-RAM arrays need it explicitly; AcqInfos.mat
%                          is not used.
%
%   Output:
%       outData     - Filtered data with the same representation type as
%                     the input.
%
%   Notes:
%       - .dat input must be single precision; other classes are refused
%         (Umitoolbox:normalizeLPF:UnsupportedDatatype).
%       - .dat input is streamed in X slabs (Y-X-T and Y-X-T-E), so the
%         recording is never loaded whole.
%       - The frame rate resolved above (explicit, else the .dat header) is
%         the one used to design the filters, for both input forms.
%       - Event-split data (an E axis) are filtered trial-by-trial along T
%         and keep the E dimension unchanged; with bApplyExpFit each trial
%         gets its own exponential fit.
%       - The exponential fit is one fit on the spatially averaged trace
%         (per trial for event-split data), ignoring NaN pixels.
%       - NaN pixels are replaced by 0 before filtering and restored
%         afterward. This biases the low-pass baseline near mask borders,
%         and with Normalize=true can feed a near-zero divisor there.
%         Pixels that are NaN or 0 over the whole recording (per trial
%         for event-split data) are masked pixels and are not filtered:
%         NaN pixels stay NaN; 0 pixels give NaN with Normalize=true
%         (0/0) and 0 otherwise.
%       - A .dat input may be the normLPF.dat this function writes: the
%         output goes through normLPF_writing.dat and is then moved.

default_Output = 'normLPF.dat'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'BaselineCutoffHz', 0.0083, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'SignalCutoffHz', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'Normalize', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'bApplyExpFit', false, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'FrameRateHz', []);

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
BaselineCutoffHz = double(p.Results.BaselineCutoffHz);
SignalCutoffHz = double(p.Results.SignalCutoffHz);
bNormalize = p.Results.Normalize;
bApplyExpFit = p.Results.bApplyExpFit;
explicitRate = p.Results.FrameRateHz;

if ~isfolder(SaveFolder)
    error('normalizeLPF:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% Resolve a file input's path/extension up front (if any) so the frame
% rate below is read from the file actually being processed, not an
% arbitrary file in SaveFolder.
isFileInput = ischar(data) || (isstring(data) && isscalar(data));
dataFile = '';
ext = '';
if isFileInput
    dataFile = char(string(data));
    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('normalizeLPF:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end
    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);
end

% Frame rate of the data itself: the explicit FrameRateHz (injected by
% PipelineManager), else the .dat header. AcqInfos.mat is not used
% (resolveDataInfoValue).
Fs = [];
nT = [];
fileInfo = [];
if isFileInput && strcmp(ext, '.dat')
    % The layout is checked first: an axis layout without time (for example
    % Y-X-E) has no frame rate to resolve.
    fileInfo = loadMetaData(dataFile);
    assertDatLayout(fileInfo, {{'Y','X','T'}, {'Y','X','T','E'}}, 'normalizeLPF');
    assert(strcmpi(char(string(fileInfo.dataClass)), 'single'), ...
        'Umitoolbox:normalizeLPF:UnsupportedDatatype', ...
        ['normalizeLPF currently expects single-precision .dat files; "%s" is %s. ' ...
         'Filtering writes the input class back, which would round the normalized ' ...
         'ratio of an integer file.'], dataFile, char(string(fileInfo.dataClass)));
    Fs = resolveDataInfoValue('frameRateHz', explicitRate, dataFile, mfilename);
    nT = datAxisSize(fileInfo, 'T');
elseif isnumeric(data) || islogical(data)
    Fs = resolveDataInfoValue('frameRateHz', explicitRate, data, mfilename);
    nT = size(data, 3);
end
if ~isempty(Fs)
    iCheckCutoffs(Fs, BaselineCutoffHz, SignalCutoffHz);
    BaselineCutoffHz = iEffectiveBaselineCutoff(BaselineCutoffHz, ...
        SignalCutoffHz, nT, Fs);
end

% -------------------------------------------------------------------------
% Case 1: YXT or YXTE array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)

    validateattributes(data, {'numeric','logical'}, {'nonempty'}, ...
        mfilename, 'data');
    if ~(ndims(data) == 3 || ndims(data) == 4)
        error('normalizeLPF:InvalidArrayInput', ...
            'Numeric input must be YXT or YXTE.');
    end

    if ndims(data) == 3
        outData = iFilterArray(single(data), BaselineCutoffHz, SignalCutoffHz, ...
            bNormalize, bApplyExpFit, Fs);
    else
        % Event-split data: every trial is filtered along T on its own.
        outData = zeros(size(data), 'single');
        for iTrial = 1:size(data, 4)
            outData(:,:,:,iTrial) = iFilterArray(single(data(:,:,:,iTrial)), ...
                BaselineCutoffHz, SignalCutoffHz, bNormalize, bApplyExpFit, Fs);
        end
    end
    return
end

% -------------------------------------------------------------------------
% Case 2: File input
% -------------------------------------------------------------------------
if isFileInput

    switch ext
        case '.dat'
            outData = fullfile(SaveFolder, 'normLPF.dat');
            iFilterFile(dataFile, fileInfo, outData, ...
                BaselineCutoffHz, SignalCutoffHz, ...
                bNormalize, bApplyExpFit, Fs);
            return

        otherwise
            error('normalizeLPF:UnsupportedInputFile', ...
                'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    end
end

error('normalizeLPF:UnsupportedInputType', ...
    ['Input "data" must be a YXT or YXTE array or a .dat filename. ' ...
     'UMT structs and .umt/.mat files are not supported.']);

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Normalize image data by low-pass filtering.');

        info.version = '2.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT or YXTE array, or a .dat ' ...
             'filename with axes Y-X-T or Y-X-T-E.'], ...
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
            'Output folder.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'BaselineCutoffHz', ...
            'parameter', ...
            ['Low cut-off frequency used to estimate the slow baseline ' ...
             'component. 0 disables the baseline filter. Raised to ' ...
             '1/(2 x duration), with a warning, for recordings shorter ' ...
             'than 1/BaselineCutoffHz.'], ...
            'kind', 'parameter', ...
            'default', 0.0083, ...
            'allowed', [0 Inf], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'SignalCutoffHz', ...
            'parameter', ...
            ['Higher cut-off frequency used to preserve the signal ' ...
             'component. 0 disables the signal low-pass.'], ...
            'kind', 'parameter', ...
            'default', 1, ...
            'allowed', [0 Inf], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'Normalize', ...
            'parameter', ...
            'If true, express the filtered signal as DeltaR/R.', ...
            'kind', 'parameter', ...
            'default', true, ...
            'allowed', [true false], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'bApplyExpFit', ...
            'parameter', ...
            'If true, apply exponential decay correction.', ...
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
            'Low-pass filtered output.', ...
            'normLPF.dat', ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: Filter one in-memory Y-X-T array while preserving NaNs
% =========================================================================
function outArray = iFilterArray(inArray, BaselineCutoffHz, SignalCutoffHz, ...
    bNormalize, bApplyExpFit, Fs)
%IFILTERARRAY Filter one numeric YXT block along T.

if ~isa(inArray, 'single')
    workArray = single(inArray);
else
    workArray = inArray;
end
dims = size(workArray);
Nt = dims(3);

idxNaN = isnan(workArray);

Pred = [];
if bApplyExpFit
    % NaN-safe spatial mean trace, as in the file path.
    flat = reshape(double(workArray), [], Nt);
    isValid = ~isnan(flat);
    flat(~isValid) = 0;
    Pred = iExpFitPredictors(sum(flat, 1) ./ sum(isValid, 1));
    clear flat isValid
end

if any(idxNaN(:))
    workArray(idxNaN) = 0;
end

lpass = iDesignLowpass(BaselineCutoffHz, Fs);
hpass = iDesignLowpass(SignalCutoffHz, Fs);

outArray = workArray;
PrcLims = round(linspace(1, dims(1), 11));
fprintf('Progress: ');
for ind = 1:dims(1)
    Signal = double(reshape(workArray(ind,:,:), dims(2), Nt));
    outArray(ind,:,:) = single(iFilterSignal(Signal, Pred, lpass, hpass, bNormalize));

    if any(ind == PrcLims)
        idx = find(ind == PrcLims);
        fprintf('%d%%..', 10*(idx-1));
    end
end
fprintf('\n');

if any(idxNaN(:))
    outArray(idxNaN) = NaN;
end
end

% =========================================================================
% Helper: Stream a Y-X-T or Y-X-T-E .dat file through the filter
% =========================================================================
function iFilterFile(inFile, fMetaData, outDat, lowFreq, highFreq, bDivide, bExpFit, Freq)
%IFILTERFILE Filter a .dat file in X slabs and write OUTDAT.
%
% FMETADATA is the input's .dat Info (axes Y-X-T or Y-X-T-E, checked by the
% caller). FREQ is the resolved frame rate of the data. Every pixel (of
% every trial) is filtered independently along T, so slabs are exact. The
% exponential fit is one whole-image fit (per trial), computed in a
% pre-pass over the slabs, so the output does not depend on the number of
% slabs (which varies with the available RAM).

storedSize = double(fMetaData.dimSizes);
dataType = char(string(fMetaData.dataClass));
hasEvents = numel(storedSize) == 4;

Ny = storedSize(1);
Nx = storedSize(2);
Nt = storedSize(3);
if hasEvents
    Ne = storedSize(4);
else
    Ne = 1;
end

% Create the headered output file (the input's axes, class, and rate)
% through a fixed-name scratch file, moved onto OUTDAT at the end: the
% input may be the file OUTDAT overwrites (a re-run on its own output).
[outFolder, outName, outExt] = fileparts(outDat);
tmpDat = fullfile(outFolder, [outName '_writing' outExt]);
slabOut = spatialSlabIO('create', tmpDat, datHeaderFromInfo(fMetaData, outName));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

lpass = iDesignLowpass(lowFreq, Freq);
hpass = iDesignLowpass(highFreq, Freq);

fprintf('normalizeLPF (file mode)\n');

slabIn = spatialSlabIO('open', inFile, 'Info', fMetaData);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

if bExpFit
    nChunks = calculateMaxChunkSize(prod(storedSize) * getByteSize(dataType), 2, .3);
else
    nChunks = calculateMaxChunkSize(prod(storedSize) * getByteSize(dataType), 1, .3);
end
chunkX = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkX);
if hasEvents
    fprintf('Filtering event-split data (%i chunk(s))...\n', nChunks)
else
    fprintf('Filtering data (%i chunk(s))...\n', nChunks)
end

% Global exponential-fit pre-pass: stream the whole-image mean trace (per
% trial) across X slabs via a running NaN-safe sum plus a running count of
% non-NaN pixels, so masked pixels do not drag the fit.
PredByTrial = cell(1, Ne);
if bExpFit
    sumS = zeros(Ne, Nt);
    countS = zeros(Ne, Nt);
    for c = 1:nChunks
        xIdx = ((c-1) * chunkX + 1):min(c * chunkX, Nx);
        slab = spatialSlabIO('read', slabIn, xIdx);
        for e = 1:Ne
            trialSlab = reshape(double(slab(:,:,:,e)), [], Nt);
            isValidS = ~isnan(trialSlab);
            trialSlab(~isValidS) = 0;
            sumS(e,:) = sumS(e,:) + sum(trialSlab, 1);
            countS(e,:) = countS(e,:) + sum(isValidS, 1);
        end
    end
    clear slab trialSlab isValidS

    for e = 1:Ne
        PredByTrial{e} = iExpFitPredictors(sumS(e,:) ./ countS(e,:));
    end
end

for c = 1:nChunks
    xIdx = ((c-1) * chunkX + 1):min(c * chunkX, Nx);
    nX = numel(xIdx);

    fprintf('Chunk #%i [Reading data from file...]\n', c)
    slab = spatialSlabIO('read', slabIn, xIdx);

    % Mask NaNs before filtering and restore them afterward, the same
    % policy as the in-RAM path. Per-chunk masking is exact here since the
    % mask is per-element and filtering runs per pixel along T.
    idxNaN = isnan(slab);
    if any(idxNaN(:))
        slab(idxNaN) = 0;
    end

    fprintf('Chunk #%i [Temporal filtering...]\n', c)
    if ~hasEvents
        fprintf('Progress:')
        lastPct = -1;
    end

    for e = 1:Ne
        for y = 1:Ny
            Signal = double(reshape(slab(y,:,:,e), nX, Nt));
            filtered = iFilterSignal(Signal, PredByTrial{e}, lpass, hpass, bDivide);
            slab(y,:,:,e) = reshape(cast(filtered, dataType), 1, nX, Nt);

            if ~hasEvents
                pct = floor(100 * y / Ny);
                if pct ~= lastPct && mod(pct, 10) == 0
                    fprintf('%d%% ', pct);
                    lastPct = pct;
                end
            end
        end
    end

    if any(idxNaN(:))
        slab(idxNaN) = NaN;
    end

    if ~hasEvents
        fprintf('\n');
    end
    fprintf('Chunk #%i [Writing to file...]\n', c)
    spatialSlabIO('write', slabOut, xIdx, slab);
    fprintf('Chunk #%i [Completed]\n', c)
    clear slab idxNaN
end
% Close every handle before the move: on Windows an open handle blocks it.
spatialSlabIO('finalize', slabOut);
spatialSlabIO('close', slabIn);

[moveOk, moveMsg] = movefile(tmpDat, outDat, 'f');
assert(moveOk, 'Umitoolbox:normalizeLPF:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpDat, outDat, moveMsg);
end

% =========================================================================
% Helpers: filtering primitives
% =========================================================================
function out = iFilterSignal(Signal, Pred, lpass, hpass, bDivide)
%IFILTERSIGNAL Filter the traces (rows) of SIGNAL (nPixels x T, double).
%
% PRED is the exponential-fit design matrix (empty: no fit). LPASS and
% HPASS are the baseline and signal low-pass filters (empty: disabled).
% Traces that are 0 over the whole recording (masked pixels, including
% NaN pixels, zeroed by the callers) carry no signal and are not
% filtered: 0/0 in the exponential fit would stop FILTFILT. They get what
% the filter gives a zero trace without the fit: NaN (0/0) with BDIVIDE,
% else 0. NaN pixels are restored to NaN by the callers. Every trace is
% filtered on its own, so leaving rows out does not change the others.

isZero = all(Signal == 0, 2);
out = zeros(size(Signal));
if bDivide
    out(isZero, :) = NaN;
end
keep = ~isZero;
if ~any(keep)
    return
end
Signal = Signal(keep, :);

if ~isempty(Pred)
    B = Pred \ Signal';
    Signal = Signal ./ (Pred * B)';
end

if ~isempty(lpass)
    LP_low = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
else
    LP_low = ones(size(Signal));
end

if ~isempty(hpass)
    LP_high = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
else
    LP_high = Signal;
end

if bDivide
    out(keep, :) = LP_high ./ LP_low;
else
    out(keep, :) = LP_high - LP_low;
end
end

function filt = iDesignLowpass(cutoffHz, Fs)
%IDESIGNLOWPASS 4th-order Butterworth low-pass; empty when CUTOFFHZ is 0.
filt = [];
if cutoffHz > 0
    f = fdesign.lowpass('N,F3dB', 4, cutoffHz, Fs);
    filt = design(f, 'butter');
end
end

function Pred = iExpFitPredictors(S)
%IEXPFITPREDICTORS Design matrix of the double exponential decay correction.
%
% S is the spatially averaged trace (1 x T, double). The decay is fitted
% with FMINSEARCH from a random start (rng('shuffle'), DFR-20260819-001);
% the design matrix holds a constant, a linear drift and the fitted decay.

ExpFun = @(P,x) abs(P(1)).*exp(-abs(P(2)).*x) + ...
                abs(P(3)).*exp(-abs(P(4)).*x) ...
                - abs(P(5))*x + P(6);
Opt = optimset('Display', 'off');
Nt = numel(S);

rng('shuffle');
initB = double(rand(1,6) .* [30 1 20 1 1 mean(S)]);
B = fminsearch(@(P) sum((S - ExpFun(P, 1:Nt)).^2), initB, Opt);
Approx = ExpFun([B(1:4) 0 0], 1:Nt);
Pred = [ones(1,Nt); linspace(0,1,Nt); Approx]';
end

% =========================================================================
% Helpers: cutoff checks
% =========================================================================
function iCheckCutoffs(Fs, BaselineCutoffHz, SignalCutoffHz)
if BaselineCutoffHz < 0 || BaselineCutoffHz > Fs/2
    error('normalizeLPF:InvalidCutoff', ...
        'BaselineCutoffHz must be between 0 and the Nyquist frequency.');
end

if SignalCutoffHz < 0 || SignalCutoffHz > Fs/2
    error('normalizeLPF:InvalidCutoff', ...
        ['SignalCutoffHz must be >= 0 (0 disables the signal low-pass) ' ...
         'and <= the Nyquist frequency.']);
end

if SignalCutoffHz > 0 && SignalCutoffHz < BaselineCutoffHz
    error('normalizeLPF:InvalidCutoff', ...
        'SignalCutoffHz must be >= BaselineCutoffHz.');
end
end

function BaselineCutoffHz = iEffectiveBaselineCutoff(BaselineCutoffHz, SignalCutoffHz, nT, Fs)
%IEFFECTIVEBASELINECUTOFF Raise a baseline cutoff whose period exceeds the recording.
%
% A baseline period longer than the recording (1/BaselineCutoffHz >
% duration) cannot be estimated, so the cutoff is raised to
% 1/(2 x duration), with a warning.

durationSec = double(nT) / double(Fs);
if BaselineCutoffHz > 0 && (1 / BaselineCutoffHz) > durationSec
    raisedHz = 1 / (2 * durationSec);
    msg = sprintf(['BaselineCutoffHz (%g Hz, a %g s period) is raised to ' ...
        '%g Hz: its period is longer than the recording (%d frames at ' ...
        '%g Hz = %g s), so the cutoff is set to 1/(2 x duration).'], ...
        BaselineCutoffHz, 1 / BaselineCutoffHz, raisedHz, nT, Fs, durationSec);
    if SignalCutoffHz > 0 && raisedHz > SignalCutoffHz
        msg = sprintf(['%s The raised cutoff is above SignalCutoffHz ' ...
            '(%g Hz).'], msg, SignalCutoffHz);
    end
    warning('Umitoolbox:normalizeLPF:BaselineCutoffRaised', '%s', msg);
    BaselineCutoffHz = raisedHz;
end
end
