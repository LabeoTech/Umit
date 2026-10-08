function outData = SpeckleMapping(data, SaveFolder, sType, bSaveMap, bLogScale)
%SPECKLEMAPPING Compute speckle contrast maps from a Y-X-T acquisition.
%
%   outData = SpeckleMapping(data, SaveFolder, sType)
%   outData = SpeckleMapping(data, SaveFolder, sType, bSaveMap, bLogScale)
%
%   This function computes the mean spatial or temporal speckle contrast
%   (local standard deviation of the intensity normalized by its temporal
%   mean) of a Y-X-T recording. The nature of DATA selects the execution
%   mode:
%
%       Filename (char or string)  -> LOW-RAM mode. The .dat file is streamed
%                                     in chunks and never loaded whole.
%       Numeric Y-X-T array        -> STANDARD mode. The array is processed
%                                     whole when the available RAM allows
%                                     it, and in chunks otherwise (see
%                                     Notes).
%
%   Regardless of the mode, the output is returned as a UMT struct of kind
%   "image".
%
%   Inputs:
%       data          - Either the name or path of a Y-X-T single-precision
%                       .dat file (a bare name is looked up in SaveFolder), or
%                       a numeric/logical Y-by-X-by-T array.
%       SaveFolder    - Existing folder. Receives the TIFF map when bSaveMap
%                       is true, and resolves a bare .dat file name.
%       sType         - Speckle mapping mode:
%                       'spatial' or 'temporal'
%       bSaveMap      - Logical scalar. If true, save TIFF map.
%                       Default: true
%       bLogScale     - Logical scalar. If true, apply -log10 transform.
%                       Default: true
%
%   Output:
%       outData       - UMT struct of kind "image" containing one entry
%                       named "SpeckleMap" with dimNames = {'Y','X'}.
%
%   Notes:
%       - A .dat file is always streamed in chunks (Low-RAM mode).
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

if nargin < 3
    error('Umitoolbox:SpeckleMapping:MissingInput', ...
        'SpeckleMapping requires DATA, SaveFolder, and sType.');
end
if nargin < 4 || isempty(bSaveMap)
    bSaveMap = true;
end
if nargin < 5 || isempty(bLogScale)
    bLogScale = true;
end

validateattributes(SaveFolder, {'char','string'}, {'nonempty'});
SaveFolder = char(string(SaveFolder));
assert(isfolder(SaveFolder), ...
    'Umitoolbox:SpeckleMapping:InvalidFolder', ...
    'The folder "%s" does not exist.', SaveFolder);

sType = lower(char(string(sType)));
assert(ismember(sType, {'spatial', 'temporal'}), ...
    'Umitoolbox:SpeckleMapping:InvalidSType', ...
    'Invalid sType. Use "spatial" or "temporal".');

% -------------------------------------------------------------------------
% Control point: the type of DATA selects the reader (and so the mode).
% -------------------------------------------------------------------------
isFileInput = ischar(data) || (isstring(data) && isscalar(data));

if isFileInput
    % LOW-RAM MODE: stream the .dat file.
    datFile = localResolveDatFile(char(string(data)), SaveFolder);
    mdIn = loadMetaData(datFile);
    assertDatLayout(mdIn, {{'Y','X','T'}}, 'SpeckleMapping');

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
    assert(isnumeric(data) || islogical(data), ...
        'Umitoolbox:SpeckleMapping:UnsupportedInputType', ...
        'Input "data" must be a .dat filename or a numeric Y-X-T array.');
    assert(ndims(data) == 3 && ~isempty(data), ...
        'Umitoolbox:SpeckleMapping:UnsupportedLayout', ...
        'Array input must be a non-empty Y x X x T array.');

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

disp('Done');
end

function datFile = localResolveDatFile(data, SaveFolder)
%LOCALRESOLVEDATFILE Resolve a .dat name or path to an existing file.
%   A bare name is looked up in SaveFolder first, a path is used as given.

if isempty(fileparts(data))
    datFile = fullfile(SaveFolder, data);
else
    datFile = data;
end

if ~isfile(datFile)
    error('Umitoolbox:SpeckleMapping:FileNotFound', ...
        'The file "%s" was not found.', data);
end

[~, ~, ext] = fileparts(datFile);
assert(strcmpi(ext, '.dat'), ...
    'Umitoolbox:SpeckleMapping:UnsupportedInputFile', ...
    'Unsupported input file extension "%s". Only .dat files are supported.', ext);
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

% -------------------------------------------------------------------------
% Pass 1: temporal mean
% -------------------------------------------------------------------------
frameOut = zeros(Ny, Nx);
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

% -------------------------------------------------------------------------
% Pass 2: speckle contrast
% -------------------------------------------------------------------------
fprintf('\nPass 2/2 - Calculating Speckle Contrast (%s algorithm)...\n', sType)

switch sType
    case 'spatial'
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
