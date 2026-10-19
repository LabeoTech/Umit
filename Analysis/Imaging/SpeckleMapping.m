function outData = SpeckleMapping(SaveFolder, data, varargin)
%SPECKLEMAPPING Compute speckle contrast maps from a Y-X-T acquisition.
%
%   outData = SpeckleMapping(SaveFolder, data)
%   outData = SpeckleMapping(SaveFolder, data, 'Name', Value, ...)
%   info    = SpeckleMapping('pipelineInfo')
%
%   This function computes the mean spatial or temporal speckle contrast
%   (local standard deviation of the intensity normalized by its temporal
%   mean) of a Laser Speckle Contrast Imaging dataset.
%
%   Inputs:
%       SaveFolder - Existing folder. Resolves a bare .dat file name and
%                    receives the TIFF map when bSaveMap is true.
%       data       - The speckle data, in one of two forms that also select
%                    the execution mode:
%                      * Name or path of a Y-X-T single-precision .dat file
%                        (a bare name is looked up in SaveFolder). The file
%                        is processed in chunks (Low-RAM mode).
%                      * Numeric Y-X-T array. It is processed in RAM
%                        (Standard mode): whole, with one STDFILT call, when
%                        the available RAM allows it, and in chunks
%                        otherwise.
%                    Other extensions, other layouts, UMT structs, and .umt
%                    files are not supported.
%
%   Name-Value parameters:
%       sType     - Standard-deviation mode.
%                   Allowed: 'Spatial', 'Temporal'
%                   Default: 'Temporal'
%
%       bSaveMap  - If true, save the speckle map as a TIFF file.
%                   Default: false
%
%       bLogScale - If true, apply -log10 transform.
%                   Default: false
%
%   Output:
%       outData   - UMT struct of kind "image" containing the speckle map
%                   (one Y-X entry named "SpeckleMap").
%
%   Notes:
%       - A .dat file is always streamed in chunks and never loaded whole.
%       - An array is processed whole, with a single STDFILT call, when
%         calculateMaxChunkSize (which queries the available RAM) reports
%         that the whole computation fits in one chunk. Otherwise it is
%         processed with the same chunked algorithm as a file, reading the
%         chunks from the array. Chunked and whole-array results agree up to
%         floating-point summation order.
%       - The peak-memory factors passed to calculateMaxChunkSize are
%         hard-coded estimates relative to the size of the data in single
%         precision: 1 for the mean pass, 10 for the spatial pass and 12 for
%         the temporal pass. The pass-2 factor also decides whether an array
%         fits whole.
%       - Array input is converted to single precision chunk by chunk when
%         chunked, so a double-precision array is only duplicated whole when
%         it is processed whole.
%       - Metadata of a .dat file are resolved through loadMetaData(...), not
%         legacy per-channel .mat sidecars.
%       - TIFF export, when requested, is saved as:
%             std_speckle.tiff

% Defaults:
default_Output = 'SpeckleMap';

allowedSType = {'Spatial', 'Temporal'};

% -------------------------------------------------------------------------
% pipelineInfo query
% -------------------------------------------------------------------------
if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) ...
        && strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'SaveFolder', @isfolder);
addRequired(p, 'data');

% The allowed values are checked after parsing, so that an unknown mode
% raises Umitoolbox:SpeckleMapping:InvalidSType.
addParameter(p, 'sType', 'Temporal', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'bSaveMap', false, ...
    @(x) islogical(x) && isscalar(x));

addParameter(p, 'bLogScale', false, ...
    @(x) islogical(x) && isscalar(x));

parse(p, SaveFolder, data, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
sTypeName = char(string(p.Results.sType));
bSaveMap = p.Results.bSaveMap;
bLogScale = p.Results.bLogScale;

assert(any(strcmpi(sTypeName, allowedSType)), ...
    'Umitoolbox:SpeckleMapping:InvalidSType', ...
    'Invalid sType. Use "spatial" or "temporal".');
sType = lower(sTypeName);

% -------------------------------------------------------------------------
% Resolve and validate the input. The type of "data" selects the mode.
% -------------------------------------------------------------------------
isFileInput = ischar(data) || (isstring(data) && isscalar(data));

if isFileInput
    % Filename: Low-RAM mode.
    data = char(string(data));

    [dataFolder, dataStem, dataExt] = fileparts(data);
    assert(isempty(dataExt) || strcmpi(dataExt, '.dat'), ...
        'Umitoolbox:SpeckleMapping:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', dataExt);

    % A bare name is looked up in SaveFolder; a path is used as given.
    if isempty(dataFolder)
        datFile = fullfile(SaveFolder, [dataStem '.dat']);
    else
        datFile = fullfile(dataFolder, [dataStem '.dat']);
    end
    assert(isfile(datFile), ...
        'Umitoolbox:SpeckleMapping:FileNotFound', ...
        'Input .dat file not found: "%s".', data);

    mdIn = loadMetaData(datFile);
    assertDatLayout(mdIn, {{'Y','X','T'}}, 'SpeckleMapping');

    mode = 'Low-RAM';
elseif isnumeric(data) || islogical(data)
    % Array: Standard mode.
    assert(ndims(data) == 3 && ~isempty(data), ...
        'Umitoolbox:SpeckleMapping:unsupportedLayout', ...
        'Array input must be a non-empty Y x X x T array.');

    mode = 'Standard';
else
    error('Umitoolbox:SpeckleMapping:UnsupportedInputType', ...
        ['Input "data" must be the name of a Y-X-T .dat file or a numeric ' ...
         'Y-X-T array.']);
end

disp(['Calculating ' sTypeName ' standard deviation (' mode ' mode)...'])

% -------------------------------------------------------------------------
% Readers: the type of DATA selects the reader (and so the mode).
% -------------------------------------------------------------------------
if isFileInput
    % LOW-RAM MODE: stream the .dat file.
    Ny = datAxisSize(mdIn, 'Y');
    Nx = datAxisSize(mdIn, 'X');
    Nt = datAxisSize(mdIn, 'T');
    assert(Nt > 0, ...
        'Umitoolbox:SpeckleMapping:InvalidMetadata', ...
        'Could not resolve the Y, X, and T sizes of "%s".', datFile);
    assert(strcmpi(char(string(mdIn.dataClass)), 'single'), ...
        'Umitoolbox:SpeckleMapping:UnsupportedDatatype', ...
        'SpeckleMapping currently expects single-precision .dat files.');

    slabIn = spatialSlabIO('open', datFile, 'Info', mdIn);
    cIn = onCleanup(@() spatialSlabIO('close', slabIn));

    readFrames = @(tIdx) reshape( ...
        spatialSlabIO('read', slabIn, 1:Nx, tIdx), Ny, Nx, numel(tIdx));
    readColumns = @(xIdx) spatialSlabIO('read', slabIn, xIdx);
else
    % STANDARD MODE: the data are already in RAM.
    [Ny, Nx, Nt] = size(data);

    readFrames = @(tIdx) single(data(:, :, tIdx));
    readColumns = @(xIdx) single(data(:, xIdx, :));
end

totalBytes = double(Ny) * double(Nx) * double(Nt) * getByteSize('single');


% -------------------------------------------------------------------------
% Speckle contrast map
% -------------------------------------------------------------------------
% Peak-memory model of the contrast pass, relative to the data in single
% precision (hard-coded estimates), and the RAM fraction left to the system.
switch sType
    case 'spatial'
        pass2Factor = 10;
        pass2Overhead = .1;
    case 'temporal'
        pass2Factor = 12;
        pass2Overhead = .15;
end
nChunksPass2 = calculateMaxChunkSize(totalBytes, pass2Factor, pass2Overhead);

if ~isFileInput && nChunksPass2 == 1
    % The whole array fits in RAM: one STDFILT call on the whole array.
    disp('Mapping computation (whole array)...');
    frameOut = localWholeArray(data, sType);
else
    frameOut = localChunked(readFrames, readColumns, Ny, Nx, Nt, ...
        totalBytes, sType, nChunksPass2);
end

if bLogScale
    frameOut = -log10(frameOut);
end

speckleMap = single(frameOut);

% -------------------------------------------------------------------------
% Save TIFF map if requested
% -------------------------------------------------------------------------
if bSaveMap
    obj = Tiff(fullfile(SaveFolder, 'std_speckle.tiff'), 'w');
    setTag(obj, 'ImageWidth', Nx);
    setTag(obj, 'ImageLength', Ny);
    setTag(obj, 'Photometric', Tiff.Photometric.MinIsBlack);
    setTag(obj, 'SampleFormat', Tiff.SampleFormat.IEEEFP);
    setTag(obj, 'BitsPerSample', 32);
    setTag(obj, 'SamplesPerPixel', 1);
    setTag(obj, 'Compression', Tiff.Compression.None);
    setTag(obj, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    write(obj, speckleMap);
    close(obj);
end

% -------------------------------------------------------------------------
% Package output as UMT
% -------------------------------------------------------------------------
outData = genUMTStruct( ...
    speckleMap, ...
    'kind', 'image', ...
    'entryName', 'SpeckleMap', ...
    'dimNames', {'Y','X'});

disp('Finished Speckle Mapping.')

    % =====================================================================
    % Local pipelineInfo factory (nested, shares allowedSType/default_Output
    % with the parent scope instead of redefining them)
    % =====================================================================
    function info = localPipelineInfo()

        info = PipelineManager.createPipelineInfo( ...
            mfilename, ...
            'Compute speckle standard-deviation maps.');
        info.version = '3.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder used to resolve the input file and to save the optional TIFF map.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            'ImageTimeSeries', ...
            ['Y-X-T speckle data: a .dat file (Low-RAM mode) or a numeric ' ...
             'array (Standard mode).'], ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput( ...
            info, ...
            'sType', ...
            'parameter', ...
            'Standard-deviation mode.', ...
            'kind', 'parameter', ...
            'default', 'Temporal', ...
            'allowed', allowedSType, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'bSaveMap', ...
            'parameter', ...
            'If true, save the speckle map as a TIFF file.', ...
            'kind', 'parameter', ...
            'default', false, ...
            'allowed', [true false], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'bLogScale', ...
            'parameter', ...
            'If true, apply -log10 transform.', ...
            'kind', 'parameter', ...
            'default', false, ...
            'allowed', [true false], ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'UMT image output containing the computed speckle map.', ...
            default_Output, ...
            1, ...
            'isData', true);

    end

end

function speckleMap = localWholeArray(data, sType)
%LOCALWHOLEARRAY Speckle contrast of an array that fits in RAM.
%   Normalizes by the temporal mean, filters the whole array with STDFILT
%   and averages the contrast over time.

dat = single(data);
dat = dat ./ mean(dat, 3, 'omitnan');

switch sType
    case 'spatial'
        Kernel = single(fspecial('disk', 2) > 0);
    case 'temporal'
        Kernel = ones(1, 1, 5, 'single');
end

speckleMap = stdfilt(dat, Kernel);
speckleMap = mean(speckleMap, 3, 'omitnan');
end

function frameOut = localChunked(readFrames, readColumns, Ny, Nx, Nt, ...
        totalBytes, sType, nChunksPass2)
%LOCALCHUNKED Speckle contrast read in chunks through the given readers.
%   READFRAMES(tIdx) returns the Y-X frames tIdx; READCOLUMNS(xIdx) returns
%   the X columns xIdx for all frames. Both return single precision.
%
%   The 'spatial' mode first needs the temporal mean of every pixel (one
%   full read). The 'temporal' mode normalizes each column slab by its own
%   temporal mean and skips that pass.

frameOut = zeros(Ny, Nx);

switch sType
    case 'spatial'
        % -----------------------------------------------------------------
        % Pass 1: temporal mean
        % -----------------------------------------------------------------
        mData = localTemporalMean(readFrames, Ny, Nx, Nt, totalBytes);

        % -----------------------------------------------------------------
        % Pass 2: speckle contrast
        % -----------------------------------------------------------------
        fprintf('\nPass 2/2 - Calculating Speckle Contrast (%s algorithm)...\n', sType)

        Kernel = single(fspecial('disk', 2) > 0);
        countOut = zeros(Ny, Nx);

        nChunks = nChunksPass2;
        chunkT  = ceil(Nt / nChunks);
        nChunks = ceil(Nt / chunkT);

        lastPct = -1;
        fprintf('0%% ');

        for c = 1:nChunks
            tStart = (c-1) * chunkT + 1;
            tEnd   = min(tStart + chunkT - 1, Nt);

            frameBlock = readFrames(tStart:tEnd);

            frameBlock = frameBlock ./ mData;
            frameBlock = stdfilt(frameBlock, Kernel);
            countOut = countOut + sum(~isnan(frameBlock), 3);
            frameBlock = sum(frameBlock, 3, 'omitnan');
            frameOut = frameOut + frameBlock;

            pct = floor(100 * c / nChunks);
            if pct ~= lastPct
                fprintf('%d%% ', pct);
                lastPct = pct;
            end

            clear frameBlock
        end

        frameOut = frameOut ./ countOut;

    case 'temporal'
        fprintf('Pass 1/1 - Calculating Speckle Contrast (%s algorithm)...\n', sType)

        Kernel = ones(1, 1, 5, 'single');

        nChunks = nChunksPass2;
        chunkX  = ceil(Nx / nChunks);
        nChunks = ceil(Nx / chunkX);

        lastPct = -1;
        fprintf('0%% ');

        for c = 1:nChunks
            xStart = (c-1) * chunkX + 1;
            xEnd   = min(xStart + chunkX - 1, Nx);
            xIdx   = xStart:xEnd;

            slab = readColumns(xIdx);
            slab = slab ./ mean(slab, 3, 'omitnan');
            slab = stdfilt(slab, Kernel);
            frameOut(:, xIdx) = mean(slab, 3, 'omitnan');

            pct = floor(100 * c / nChunks);
            if pct ~= lastPct
                fprintf('%d%% ', pct);
                lastPct = pct;
            end

            clear slab
        end
end
end

function mData = localTemporalMean(readFrames, Ny, Nx, Nt, totalBytes)
%LOCALTEMPORALMEAN Temporal mean of every pixel, streamed in time chunks.

mData = zeros(Ny, Nx, 'single');
countData = zeros(Ny, Nx, 'single');

disp('Pass 1/2 - Calculating temporal mean...')

nChunks = calculateMaxChunkSize(totalBytes, 1, .1);
chunkT  = ceil(Nt / nChunks);
nChunks = ceil(Nt / chunkT);

lastPct = -1;
fprintf('0%% ');

for c = 1:nChunks
    tStart = (c-1) * chunkT + 1;
    tEnd   = min(tStart + chunkT - 1, Nt);

    slab = readFrames(tStart:tEnd);

    % Use omitnan during accumulation.
    mData = mData + sum(slab, 3, 'omitnan');
    countData = countData + sum(~isnan(slab), 3);

    pct = floor(100 * c / nChunks);
    if pct ~= lastPct
        fprintf('%d%% ', pct);
        lastPct = pct;
    end

    clear slab
end

mData = mData ./ countData;
end
