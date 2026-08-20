function outData = SpeckleMapping(folderPath, sType, channel, bSaveMap, bLogScale, bRAMSafeMode)
%SPECKLEMAPPING Compute speckle contrast maps from a 3D acquisition.
%
%   outData = SpeckleMapping(folderPath, sType, channel, bSaveMap, bLogScale, bRAMSafeMode)
%
%   This function computes speckle contrast maps from a channel .dat file
%   stored in "folderPath". The computation can be performed in standard
%   mode or in RAM-safe mode. Regardless of the internal computation mode,
%   the output is returned as a UMT struct of kind "image".
%
%   Inputs:
%       folderPath    - Folder containing the channel .dat file.
%       sType         - Speckle mapping mode:
%                       'spatial' or 'temporal'
%       channel       - Channel name to process.
%                       Default: 'speckle'
%       bSaveMap      - Logical scalar. If true, save TIFF map.
%                       Default: true
%       bLogScale     - Logical scalar. If true, apply -log10 transform.
%                       Default: true
%       bRAMSafeMode  - Logical scalar. If true, use low-RAM processing.
%                       Default: false
%
%   Output:
%       outData       - UMT struct of kind "image" containing one entry
%                       named "SpeckleMap" with dimNames = {'Y','X'}.
%
%   Notes:
%       - The numerical algorithm is preserved from the legacy version.
%       - Metadata are resolved through loadMetaData(...), not legacy
%         per-channel .mat sidecars.
%       - TIFF export, when requested, is saved as:
%             std_speckle.tiff

if nargin < 3 || isempty(channel)
    channel = 'speckle';
end
if nargin < 4 || isempty(bSaveMap)
    bSaveMap = true;
end
if nargin < 5 || isempty(bLogScale)
    bLogScale = true;
end
if nargin < 6 || isempty(bRAMSafeMode)
    bRAMSafeMode = false;
end

validateattributes(folderPath, {'char','string'}, {'nonempty'});
folderPath = char(string(folderPath));
assert(isfolder(folderPath), ...
    'Umitoolbox:SpeckleMapping:InvalidFolder', ...
    'The folder "%s" does not exist.', folderPath);

if ~strcmp(folderPath(end), filesep)
    folderPath = [folderPath filesep]; %#ok<AGROW>
end

channel = lower(char(string(channel)));
datFile = [folderPath channel '.dat'];

assert(isfile(datFile), ...
    'Umitoolbox:SpeckleMapping:MissingChannelFile', ...
    'The file "%s" is missing.', datFile);

mdIn = loadMetaData(datFile);

% Resolve modern metadata fields.
assert(isfield(mdIn, 'Height') && isfield(mdIn, 'Width') && isfield(mdIn, 'Length'), ...
    'Umitoolbox:SpeckleMapping:InvalidMetadata', ...
    'Could not resolve Height/Width/Length from "%s".', datFile);

Ny = double(mdIn.Height);
Nx = double(mdIn.Width);
Nt = double(mdIn.Length);

if isfield(mdIn, 'Datatype') && ~isempty(mdIn.Datatype)
    dataType = char(string(mdIn.Datatype));
else
    dataType = 'single';
end

assert(strcmpi(dataType, 'single'), ...
    'Umitoolbox:SpeckleMapping:UnsupportedDatatype', ...
    'SpeckleMapping currently expects single-precision .dat files.');

if bRAMSafeMode
    % ---------------------------------------------------------------------
    % LOW-RAM MODE
    % ---------------------------------------------------------------------
    frameOut = zeros(Ny, Nx);
    mData = zeros(Ny, Nx, 'single');
    countData = zeros(Ny, Nx, 'single');

    fidIn = fopen(datFile, 'r');
    assert(fidIn ~= -1, ...
        'Umitoolbox:SpeckleMapping:FileOpenFailed', ...
        'Could not open "%s" for reading.', datFile);
    cIn = onCleanup(@() safeFclose(fidIn)); %#ok<NASGU>

    % Pass 1: temporal mean
    disp('Pass 1/2 - Calculating temporal mean...')

    bytesPerFrame = Ny * Nx * getByteSize(dataType);
    totalBytes = Nx * Ny * Nt * getByteSize(dataType);

    nChunks = calculateMaxChunkSize(totalBytes, 1, .1);
    chunkT  = ceil(Nt / nChunks);
    nChunks = ceil(Nt / chunkT);

    lastPct = -1;
    fprintf('0%% ');

    for c = 1:nChunks
        tStart  = (c-1) * chunkT + 1;
        tEnd    = min(tStart + chunkT - 1, Nt);
        nFrames = tEnd - tStart + 1;

        fseek(fidIn, (tStart-1) * bytesPerFrame, 'bof');
        slab = fread(fidIn, Ny * Nx * nFrames, ['*' dataType]);
        slab = reshape(slab, Ny, Nx, nFrames);

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

    fprintf('\nPass 2/2 - Calculating Speckle Contrast (%s algorithm)...\n', sType)

    switch lower(sType)
        case 'spatial'
            Kernel = single(fspecial('disk', 2) > 0);
            countOut = zeros(Ny, Nx);

            nChunks = calculateMaxChunkSize(totalBytes, 10, .1);
            chunkT  = ceil(Nt / nChunks);
            nChunks = ceil(Nt / chunkT);

            lastPct = -1;
            fprintf('0%% ');

            for c = 1:nChunks
                tStart  = (c-1) * chunkT + 1;
                tEnd    = min(tStart + chunkT - 1, Nt);
                nFrames = tEnd - tStart + 1;

                fseek(fidIn, (tStart-1) * bytesPerFrame, 'bof');
                frameBlock = fread(fidIn, Ny * Nx * nFrames, ['*' dataType]);
                frameBlock = reshape(frameBlock, Ny, Nx, nFrames);

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

            nChunks = calculateMaxChunkSize(Nx * Ny * Nt * getByteSize(dataType), 12, .15);
            chunkX  = ceil(Nx / nChunks);
            nChunks = ceil(Nx / chunkX);

            lastPct = -1;
            fprintf('0%% ');

            for c = 1:nChunks
                xStart = (c-1) * chunkX + 1;
                xEnd   = min(xStart + chunkX - 1, Nx);
                xIdx   = xStart:xEnd;

                slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, dataType);
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

        otherwise
            error('Umitoolbox:SpeckleMapping:InvalidSType', ...
                'Invalid sType. Use "spatial" or "temporal".');
    end

    if bLogScale
        frameOut = -log10(frameOut);
    end

    speckleMap = single(frameOut);

else
    % ---------------------------------------------------------------------
    % STANDARD MODE
    % ---------------------------------------------------------------------
    fid = fopen(datFile, 'r');
    assert(fid ~= -1, ...
        'Umitoolbox:SpeckleMapping:FileOpenFailed', ...
        'Could not open "%s" for reading.', datFile);
    c = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

    dat = fread(fid, inf, '*single');
    dat = reshape(dat, Ny, Nx, Nt);
    dat = dat ./ mean(dat, 3, 'omitnan');

    disp('Mapping computation...');

    switch lower(sType)
        case 'spatial'
            Kernel = single(fspecial('disk', 2) > 0);
        case 'temporal'
            Kernel = ones(1, 1, 5, 'single');
        otherwise
            error('Umitoolbox:SpeckleMapping:InvalidSType', ...
                'Invalid sType. Use "spatial" or "temporal".');
    end

    speckleMap = stdfilt(dat, Kernel);
    speckleMap = mean(speckleMap, 3, 'omitnan');

    if bLogScale
        speckleMap = -log10(speckleMap);
    end

    speckleMap = single(speckleMap);
end

% -------------------------------------------------------------------------
% Save TIFF map if requested
% -------------------------------------------------------------------------
if bSaveMap
    obj = Tiff(fullfile(folderPath, 'std_speckle.tiff'), 'w');
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
