function varargout = Ana_Speckle(data, SaveFolder, bNormalize, varargin)
%ANA_SPECKLE Calculate blood-flow maps from speckle data.
%
%   Ana_Speckle(data, SaveFolder, bNormalize)
%   out = Ana_Speckle(data, SaveFolder, bNormalize)
%   [out, metaData] = Ana_Speckle(data, SaveFolder, bNormalize, ...)
%
%   The nature of DATA selects the execution mode:
%
%       Filename (char or string)  -> LOW-RAM mode. The .dat file is read,
%                                     processed, and written in chunks;
%                                     "Flow.dat" is written in SaveFolder.
%       Numeric Y-X-T array        -> STANDARD mode. The array is processed
%                                     in RAM and the flow array is returned
%                                     (or saved as "Flow.dat" when no output
%                                     is requested).
%
%   metaData describes the Flow.dat output: in Low-RAM mode it is
%   loadMetaData of the written file; in standard mode (no file written)
%   it holds the .dat Info schema fields the file would have (filePath,
%   format, dataOffset, dataClass, dimNames, dimSizes, frameRateHz,
%   exposureMsec, channelName, writeComplete).
%
%   This function computes blood-flow maps from laser speckle data using
%   the existing algorithm:
%       1) always-on removal of static structure: each frame is divided by
%          MeanMap (the per-pixel temporal mean of the raw input), which
%          otherwise contaminates the local spatial std/mean ratio below
%       2) local contrast estimation from each frame
%       3) conversion from contrast to flow using private_flow_from_contrast
%       4) temporal median filtering
%       5) optional output-level normalization of the finished flow map by
%          its own per-pixel temporal mean
%
%   Inputs:
%       data         - Either the name or path of a Y-X-T speckle .dat file
%                      (a bare name, with or without ".dat", is looked up
%                      in SaveFolder), or a numeric Y-X-T array with at
%                      least 2 frames.
%       SaveFolder   - Existing folder. Resolves a bare file name and
%                      receives "Flow.dat".
%       bNormalize   - Logical scalar. If true, normalize the finished flow
%                      map (after temporal median filtering) by its own
%                      per-pixel temporal mean. Does not affect the
%                      always-on MeanMap correction in step 1, which runs
%                      regardless of this flag.
%
%   Name-Value parameters:
%       'FrameRateHz'  - Frame rate of DATA (Hz). The .dat header provides
%                        it for a file; an array needs it explicitly. An
%                        explicit value wins over the header (resolveDataInfoValue).
%       'ExposureMsec' - Speckle exposure time (ms), with the same rules as
%                        'FrameRateHz'.
%
%   Outputs:
%       out      - Array input returns a single Y x X x T blood-flow array.
%                  File input returns the full path to the Flow.dat output.
%       metaData - Flat compatibility metadata describing the output.
%
%   Notes:
%       - The MeanMap (or per-frame equivalent in Low-RAM mode) correction
%         is always applied, independent of bNormalize; bNormalize affects
%         only the output-level normalization in step 5.
%       - One flow frame is calculated from each input frame, so the output
%         retains the input temporal length T and follows the normal raw
%         .dat timeline contract.
%       - AcqInfos.mat is not used: the frame rate and the speckle exposure
%         come from the .dat header or from the Name-Value parameters.

% Default output for pipeline management.
default_Output = 'Flow.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    varargout{1} = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'Ana_Speckle';
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addRequired(p, 'bNormalize', @(x) islogical(x) && isscalar(x));
addParameter(p, 'FrameRateHz', []);
addParameter(p, 'ExposureMsec', []);
parse(p, data, SaveFolder, bNormalize, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
bNormalize = p.Results.bNormalize;

fprintf('Running Ana_Speckle...\n');

% -------------------------------------------------------------------------
% Control point: the type of DATA selects the mode and where the metadata
% come from.
% -------------------------------------------------------------------------
bLowRAM = ischar(data) || (isstring(data) && isscalar(data));

if bLowRAM
    datFile = localResolveDatFile(char(string(data)), SaveFolder);
    Iptr = loadMetaData(datFile);
    assertDatLayout(Iptr, {{'Y','X','T'}}, 'Ana_Speckle');
    metaSource = datFile;
else
    assert(isnumeric(data) || islogical(data), ...
        'Ana_Speckle:UnsupportedInputType', ...
        'Input "data" must be a .dat filename or a numeric Y-X-T array.');
    assert(ndims(data) == 3 && ~isempty(data), ...
        'Ana_Speckle:UnsupportedLayout', ...
        'Array input must be a non-empty Y x X x T array.');
    Iptr = struct('dataClass', 'single', 'dimNames', {{'Y','X','T'}}, ...
        'dimSizes', size(data), 'frameRateHz', NaN, 'exposureMsec', NaN);
    metaSource = data;
end

ny = datAxisSize(Iptr, 'Y');
nx = datAxisSize(Iptr, 'X');
nt = datAxisSize(Iptr, 'T');

% Explicit value > .dat header > error. AcqInfos.mat is not used.
tFreq = resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, ...
    metaSource, 'Ana_Speckle');
exposureMsec = resolveDataInfoValue('exposureMsec', p.Results.ExposureMsec, ...
    metaSource, 'Ana_Speckle');
assert(exposureMsec > 0, 'Ana_Speckle:InvalidExposure', ...
    'The speckle exposure must be positive (got %g ms).', exposureMsec);
speckle_int_time = exposureMsec / 1000;

OPTIONS.GPU = 0;
OPTIONS.Power2Flag = 0;
OPTIONS.Brep = 0;

% Header of the Flow.dat output: the input's axes and sizes, with the
% resolved rate and exposure, stored as single.
[~, outBaseName] = fileparts(default_Output);
outHeader = datHeaderFromInfo(Iptr, outBaseName, 'dataClass', 'single', ...
    'frameRateHz', tFreq, 'exposureMsec', exposureMsec);

assert(nt >= 2, 'Ana_Speckle:InvalidInputLength', ...
    'Speckle input must contain at least 2 frames.');

%% ------------------------------------------------------------------------
% Low-RAM mode
% -------------------------------------------------------------------------
if bLowRAM
    % Compute through a fixed-name raw scratch file (bounded RAM via slab
    % I/O), then move the completed file onto the declared Flow.dat output.
    % Renaming the declared output when it already exists would make every
    % pipeline re-run write to a different file and leave the stale original
    % in place.
    outFile = fullfile(SaveFolder, default_Output);
    computeScratchFile = fullfile(SaveFolder, [outBaseName '_compute.dat']);

    slabIn = spatialSlabIO('open', datFile, 'Info', Iptr);
    cIn = onCleanup(@() spatialSlabIO('close', slabIn));

    slabOut = spatialSlabIO('create', computeScratchFile, outHeader);
    cOut = onCleanup(@() spatialSlabIO('close', slabOut));

    % Pass 1: temporal mean
    fprintf('PASS 1/3: Calculating temporal mean\n');
    MeanMap = zeros(ny, nx, 'single');
    lastPct = -1;

    for t = 1:nt
        frame = single(spatialSlabIO('read', slabIn, 1:nx, t));
        frame = reshape(frame, ny, nx);
        MeanMap = MeanMap + frame;

        pct = floor(100 * t / nt);
        if pct ~= lastPct && mod(pct, 10) == 0
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end
    MeanMap = MeanMap / nt;

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
        flow       = single(private_flow_from_contrast(contrast, speckle_int_time));

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
    assert(moveOk, 'Ana_Speckle:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', computeScratchFile, outFile, moveMsg);

    if nargout > 0
        varargout{1} = outFile;
        if nargout > 1
            varargout{2} = loadMetaData(outFile);
        end
    end
    fprintf('\nDone!\n');
    return
end

%% ------------------------------------------------------------------------
% Standard mode
% -------------------------------------------------------------------------
dat = single(data);
MeanMap = mean(dat, 3);

datOut = zeros(nt, ny, nx, 'single');
speckle_window = fspecial('disk', 2) > 0;

for t = 1:nt
    tmp_laser = dat(:,:,t) ./ MeanMap;

    std_laser = imgaussfilt(stdfilt(tmp_laser, speckle_window), 1);
    mean_laser = imgaussfilt(convnfft(tmp_laser, speckle_window, 'same', 1:2, OPTIONS) / sum(speckle_window(:)), 1);
    contrast = std_laser ./ mean_laser;
    datOut(t,:,:) = single(private_flow_from_contrast(contrast, speckle_int_time));
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

outFile = fullfile(SaveFolder, default_Output);

if nargout > 0
    varargout{1} = datOut;
    if nargout > 1
        varargout{2} = iPlannedInfo(outFile, outHeader);
    end
else
    fprintf('Saving data to file: "%s"...\n', default_Output);
    saveData(outFile, datOut, 'DimNames', {'Y', 'X', 'T'}, ...
        'Info', struct('frameRateHz', tFreq, 'exposureMsec', exposureMsec));
end

fprintf('Done!\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'Ana_Speckle', ...
            'Calculate blood-flow maps from speckle data.');

        info = PipelineManager.addInput(info, ...
            'data', ...
            'ImageTimeSeries', ...
            'Y-X-T speckle data: a .dat file (Low-RAM mode) or a numeric array (standard mode).', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder resolving a bare file name and receiving Flow.dat.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'bNormalize', ...
            'parameter', ...
            'If true, normalize the finished flow map by its own temporal mean.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'positional', ...
            'default', false, ...
            'allowed', {false,true}, ...
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
            'Blood-flow output: a single Y x X x T image time series.', ...
            default_Output, ...
            1, ...
            'isData', true, ...
            'saveFileName', default_Output);
    end
end

%% Local function
function speed = private_flow_from_contrast(contrast,T)
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

function datFile = localResolveDatFile(data, SaveFolder)
%LOCALRESOLVEDATFILE Resolve a .dat name or path to an existing file.
%   A bare name (with or without ".dat") is looked up in SaveFolder, a path
%   is used as given.

[folder, stem, ext] = fileparts(data);
assert(isempty(ext) || strcmpi(ext, '.dat'), ...
    'Ana_Speckle:UnsupportedInputFile', ...
    'Unsupported input file extension "%s". Only .dat files are supported.', ext);

if isempty(folder)
    datFile = fullfile(SaveFolder, [stem '.dat']);
else
    datFile = fullfile(folder, [stem '.dat']);
end

if ~isfile(datFile)
    error('Ana_Speckle:InputFileNotFound', ...
        'No speckle data file was found: "%s".', datFile);
end
end

function Info = iPlannedInfo(filePath, hdr)
%IPLANNEDINFO .dat Info schema of a headered file that was not written.
[~, codes] = datHeaderSchema(1);
Info = struct('filePath', filePath, 'format', 'header', ...
    'dataOffset', double(codes.constants.headerLength), ...
    'dataClass', hdr.dataClass, 'dimNames', {hdr.dimNames}, ...
    'dimSizes', hdr.dimSizes, 'frameRateHz', hdr.frameRateHz, ...
    'exposureMsec', hdr.exposureMsec, 'channelName', hdr.channelName, ...
    'writeComplete', true);
end
