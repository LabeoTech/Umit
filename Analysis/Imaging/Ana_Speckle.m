function outData = Ana_Speckle(SaveFolder, data, varargin)
%ANA_SPECKLE Calculate blood-flow maps from speckle data.
%
%   outData = Ana_Speckle(SaveFolder, data)
%   outData = Ana_Speckle(SaveFolder, data, 'bNormalize', true)
%   outData = Ana_Speckle(SaveFolder, array, 'FrameRateHz', 10, 'ExposureMsec', 5)
%
%   The nature of "data" selects the execution mode:
%
%       Filename (char or string)  -> LOW-RAM mode. The .dat file is read,
%                                     processed, and written in chunks;
%                                     "Flow.dat" is written in SaveFolder.
%       Numeric Y-X-T array        -> STANDARD mode. The array is processed
%                                     in RAM and the flow array is returned.
%
%   Algorithm:
%       1) always-on removal of static structure: each frame is divided by
%          MeanMap (the per-pixel temporal mean of the raw input), which
%          otherwise contaminates the local spatial std/mean ratio below
%       2) local contrast estimation from each frame
%       3) conversion from contrast to flow (iFlowFromContrast)
%       4) temporal median filtering over 0.5 s
%       5) optional output-level normalization of the finished flow map by
%          its own per-pixel temporal mean
%
%   Inputs:
%       SaveFolder - Existing folder. Resolves a bare .dat file name and
%                    receives the Flow.dat output.
%       data       - The speckle data, in one of two forms:
%                      * Name or path of a Y-X-T .dat file (a bare name,
%                        with or without ".dat", is looked up in SaveFolder).
%                        The file is read, processed, and written in chunks
%                        (Low-RAM mode).
%                      * Numeric Y-X-T array with at least 2 frames. It is
%                        processed in RAM (Standard mode).
%                    Other extensions, other layouts, UMT structs, and .umt
%                    files are not supported.
%
%   Name-Value parameters:
%       'bNormalize'   - Logical scalar. Normalize the finished flow map by
%                        its per-pixel temporal mean (step 5). Does not
%                        affect the always-on MeanMap correction (step 1).
%                        Default: false
%       'FrameRateHz'  - Frame rate of the data (Hz). PipelineManager injects
%                        it from the data; a .dat file's header provides it
%                        otherwise. An array needs it explicitly. An
%                        explicit value wins over the header
%                        (resolveDataInfoValue).
%       'ExposureMsec' - Speckle exposure time (ms), with the same rules as
%                        'FrameRateHz'. AcqInfos.mat is not used.
%
%   Output:
%       outData    - File input: full path of the Y-X-T Flow.dat output in
%                    SaveFolder (single precision, same axes, rate, and
%                    exposure as the input).
%                    Array input: single-precision Y-X-T flow array.
%
%   Notes:
%       - The layout of a file is checked first, so a refused input writes
%         nothing.
%       - One flow frame is calculated from each input frame, so the output
%         retains the input temporal length T.
%       - Low-RAM mode writes through Flow_compute.dat and then moves it
%         onto Flow.dat, so the input may be a previous Flow.dat.

% Default output for pipeline management.
default_Output = 'Flow.dat';

if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) && ...
        strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'Ana_Speckle';
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addRequired(p, 'data');
addParameter(p, 'bNormalize', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FrameRateHz', []);
addParameter(p, 'ExposureMsec', []);
parse(p, SaveFolder, data, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
bNormalize = p.Results.bNormalize;

% -------------------------------------------------------------------------
% Resolve and validate the input. The type of "data" selects the mode and
% where the metadata come from.
% -------------------------------------------------------------------------
bLowRAM = ischar(data) || (isstring(data) && isscalar(data));

if bLowRAM
    % Filename: Low-RAM mode.
    [dataFolder, dataStem, dataExt] = fileparts(char(string(data)));
    assert(isempty(dataExt) || strcmpi(dataExt, '.dat'), ...
        'Umitoolbox:Ana_Speckle:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', dataExt);

    % A bare name is looked up in SaveFolder; a path is used as given.
    if isempty(dataFolder)
        datFile = fullfile(SaveFolder, [dataStem '.dat']);
    else
        datFile = fullfile(dataFolder, [dataStem '.dat']);
    end
    assert(isfile(datFile), ...
        'Umitoolbox:Ana_Speckle:FileNotFound', ...
        'Speckle .dat file not found: "%s".', datFile);
    Iptr = loadMetaData(datFile);
    assertDatLayout(Iptr, {{'Y','X','T'}}, 'Ana_Speckle');

    metaSource = datFile;
    mode = 'Low-RAM';
elseif isnumeric(data) || islogical(data)
    % Array: Standard mode.
    assert(ndims(data) == 3 && ~isempty(data), ...
        'Umitoolbox:Ana_Speckle:unsupportedLayout', ...
        'Array input must be a non-empty Y x X x T array.');

    Iptr = struct('dataClass', 'single', 'dimNames', {{'Y','X','T'}}, ...
        'dimSizes', size(data), 'frameRateHz', NaN, 'exposureMsec', NaN);
    metaSource = data;
    mode = 'Standard';
else
    error('Umitoolbox:Ana_Speckle:UnsupportedInputType', ...
        ['Input "data" must be the name of a Y-X-T .dat file or a numeric ' ...
         'Y-X-T array.']);
end

fprintf('Calculating blood flow (%s mode)...\n', mode);

% Explicit value > .dat header > error. AcqInfos.mat is not used.
tFreq = resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, ...
    metaSource, 'Ana_Speckle');
exposureMsec = resolveDataInfoValue('exposureMsec', p.Results.ExposureMsec, ...
    metaSource, 'Ana_Speckle');
assert(exposureMsec > 0, 'Umitoolbox:Ana_Speckle:InvalidExposure', ...
    'The speckle exposure must be positive (got %g ms).', exposureMsec);

assert(datAxisSize(Iptr, 'T') >= 2, 'Umitoolbox:Ana_Speckle:InvalidInputLength', ...
    'Speckle input must contain at least 2 frames.');

if bLowRAM
    % Header of the Flow.dat output: the input's axes and sizes, with the
    % resolved rate and exposure, stored as single.
    [~, outBaseName] = fileparts(default_Output);
    outHeader = datHeaderFromInfo(Iptr, outBaseName, 'dataClass', 'single', ...
        'frameRateHz', tFreq, 'exposureMsec', exposureMsec);
    outData = iFlowLowRAM(datFile, Iptr, fullfile(SaveFolder, default_Output), ...
        outHeader, tFreq, exposureMsec, bNormalize);
else
    outData = iFlowStandard(data, tFreq, exposureMsec, bNormalize);
end

fprintf('Finished Ana_Speckle.\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'Ana_Speckle', ...
            'Calculate blood-flow maps from speckle data.');
        info.version = '3.0.0';

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder resolving a bare file name and receiving Flow.dat.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'data', ...
            'ImageTimeSeries', ...
            'Y-X-T speckle data: a .dat file (Low-RAM mode) or a numeric array (Standard mode).', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'bNormalize', ...
            'parameter', ...
            'If true, normalize the finished flow map by its own temporal mean.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', false, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

        info = PipelineManager.addInput(info, ...
            'ExposureMsec', ...
            'sourceInfo', ...
            'Speckle exposure of the input data (ms), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'exposureMsec', ...
            'required', false);

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            {'ImageTimeSeries', 'ProcessedData'}, ...
            'data', ...
            'Blood-flow output: a single Y-X-T .dat image time series.', ...
            default_Output, ...
            1, ...
            'isData', true, ...
            'saveFileName', default_Output);
    end
end

% =========================================================================
% Low-RAM mode
% =========================================================================
function outFile = iFlowLowRAM(datFile, Iptr, outFile, outHeader, tFreq, exposureMsec, bNormalize)
%IFLOWLOWRAM Stream DATFILE through the flow algorithm into OUTFILE.
%
% Computes through a fixed-name scratch file (bounded RAM via slab I/O),
% then moves the completed file onto OUTFILE. Renaming the declared output
% when it already exists would make every pipeline re-run write to a
% different file and leave the stale original in place.

ny = datAxisSize(Iptr, 'Y');
nx = datAxisSize(Iptr, 'X');
nt = datAxisSize(Iptr, 'T');
speckle_int_time = exposureMsec / 1000;
OPTIONS = iConvOptions();

[outFolder, outBaseName] = fileparts(outFile);
computeScratchFile = fullfile(outFolder, [outBaseName '_compute.dat']);

slabIn = spatialSlabIO('open', datFile, 'Info', Iptr);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

slabOut = spatialSlabIO('create', computeScratchFile, outHeader);
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

% Pass 1: temporal mean, accumulated in double so long recordings do not
% lose precision; the map itself is single, as in Standard mode.
fprintf('PASS 1/3: Calculating temporal mean\n');
MeanMap = zeros(ny, nx);
lastPct = -1;

for t = 1:nt
    frame = single(spatialSlabIO('read', slabIn, 1:nx, t));
    frame = reshape(frame, ny, nx);
    MeanMap = MeanMap + double(frame);

    pct = floor(100 * t / nt);
    if pct ~= lastPct && mod(pct, 10) == 0
        fprintf('%d%% ', pct);
        lastPct = pct;
    end
end
MeanMap = single(MeanMap / nt);

% Pass 2: compute flow frame-by-frame
fprintf('\nPASS 2/3: Computing flow...\n');
lastPct = -1;
speckle_window = fspecial('disk', 2) > 0;

for t = 1:nt
    frameNext = single(spatialSlabIO('read', slabIn, 1:nx, t));
    frameNext = reshape(frameNext, ny, nx);
    frameNext = frameNext ./ MeanMap;

    std_laser  = imgaussfilt(stdfilt(frameNext, speckle_window), 1);
    mean_laser = imgaussfilt(convnfft(frameNext, speckle_window, 'same', 1:2, OPTIONS) / sum(speckle_window(:)), 1);
    contrast   = std_laser ./ mean_laser;
    flow       = single(iFlowFromContrast(contrast, speckle_int_time));

    spatialSlabIO('write', slabOut, 1:nx, flow, t);

    pct = floor(100 * t / nt);
    if pct ~= lastPct && mod(pct, 10) == 0
        fprintf('%d%% ', pct);
        lastPct = pct;
    end
end

% Pass 3: temporal median filter in X chunks
fW = ceil(0.5 * tFreq);
nChunks = calculateMaxChunkSize(ny * nx * nt * getByteSize('single'), 2);
chunkX = ceil(nx / nChunks);
nChunks = ceil(nx / chunkX);

fprintf('\nPASS 3/3: Applying temporal median filter...\n');
lastPct = -1;

for c = 1:nChunks
    xStart = (c-1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, nx);
    xIdx   = xStart:xEnd;

    slab = spatialSlabIO('read', slabOut, xIdx);
    slab = medfilt1(slab, fW, [], 3, 'truncate');
    spatialSlabIO('write', slabOut, xIdx, slab);

    pct = floor(100 * c / nChunks);
    if pct ~= lastPct
        fprintf('%d%% ', pct);
        lastPct = pct;
    end
end

% Pass 4: output-level normalization (streaming), only when requested.
% Reuses the same X-chunking as Pass 3 rather than a new chunk scheme.
if bNormalize
    fprintf('\nNormalizing output by its own temporal mean...\n');
    lastPct = -1;

    for c = 1:nChunks
        xStart = (c-1) * chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, nx);
        xIdx   = xStart:xEnd;

        slab = spatialSlabIO('read', slabOut, xIdx);
        slab = slab ./ mean(slab, 3);
        spatialSlabIO('write', slabOut, xIdx, slab);

        pct = floor(100 * c / nChunks);
        if pct ~= lastPct
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close both handles before replacing the declared output

[moveOk, moveMsg] = movefile(computeScratchFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:Ana_Speckle:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', computeScratchFile, outFile, moveMsg);
fprintf('\nDone!\n');
end

% =========================================================================
% Standard mode
% =========================================================================
function datOut = iFlowStandard(data, tFreq, exposureMsec, bNormalize)
%IFLOWSTANDARD Flow of an in-RAM Y-X-T array.

[ny, nx, nt] = size(data);
speckle_int_time = exposureMsec / 1000;
OPTIONS = iConvOptions();

dat = single(data);
MeanMap = mean(dat, 3);

datOut = zeros(nt, ny, nx, 'single');
speckle_window = fspecial('disk', 2) > 0;

for t = 1:nt
    tmp_laser = dat(:,:,t) ./ MeanMap;

    std_laser = imgaussfilt(stdfilt(tmp_laser, speckle_window), 1);
    mean_laser = imgaussfilt(convnfft(tmp_laser, speckle_window, 'same', 1:2, OPTIONS) / sum(speckle_window(:)), 1);
    contrast = std_laser ./ mean_laser;
    datOut(t,:,:) = single(iFlowFromContrast(contrast, speckle_int_time));
end

% Temporal median filter
fW = ceil(0.5 * tFreq);
datOut = medfilt1(datOut, fW, [], 1, 'truncate');

% Finalize to Y X T
datOut = permute(datOut, [2 3 1]);

% Output-level normalization by the flow map's own per-pixel temporal
% mean. Independent of the always-on MeanMap correction applied above.
if bNormalize
    datOut = datOut ./ mean(datOut, 3);
end

fprintf('Done!\n');
end

% =========================================================================
% Shared helpers
% =========================================================================
function OPTIONS = iConvOptions()
%ICONVOPTIONS CONVNFFT options: CPU FFT, exact lengths, no reuse.
OPTIONS.GPU = 0;
OPTIONS.Power2Flag = 0;
OPTIONS.Brep = 0;
end

function speed = iFlowFromContrast(contrast,T)
%IFLOWFROMCONTRAST Convert speckle contrast to flow (1/tau) for exposure T (s).
contrast(isnan(contrast)|contrast<0)=0;
contrast2 = contrast(3:end-2,3:end-2);
mmean = mean(contrast2(:));
sstd  = std(contrast2(:));
tau=(logspace(-15,0,60).^.5);
K  = ((tau/(2*T)).*(1-exp(-2*T*ones(size(tau))./tau))).^(1/2);
[~, index1] = find(K>(mmean-3*sstd),1);
[~, index2] = find(K>(mmean+3*sstd),1);
if isempty(index1), index1=1; end
if isempty(index2)||index2==index1, index2=60; end
Tau2=(logspace(log10(tau(index1)),log10(tau(index2)),40));
K  = ((Tau2/(2*T)).*(1-exp(-2*T*ones(size(Tau2))./Tau2))).^(1/2);
Tau2=[Tau2(1) Tau2 Tau2(end)];
K=[0 K 1e30];
speed=1./interp1(K,Tau2,contrast);
end
