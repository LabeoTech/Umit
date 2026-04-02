function DatOut = SpeckleMapping(folderPath, sType, channel, bSaveMap, bLogScale, bRAMSafeMode)
% SPECKLEMAPPING computes speckle contrast maps from a 3D acquisition.
%
% Supports low-RAM processing with bRAMSafeMode flag.
%
% INPUTS:
%   folderPath    : folder containing channel .dat and .mat files
%   sType         : 'spatial' or 'temporal'
%   channel       : channel to process ('speckle' by default)
%   bSaveMap      : boolean, save a .tiff map
%   bLogScale     : boolean, apply -log10 scaling
%   bRAMSafeMode  : boolean, enable low-RAM hybrid processing
%
% OUTPUT:
%   DatOut : speckle map (numeric) or filename if bRAMSafeMode is true

if nargin < 3, channel = 'speckle'; end
if nargin < 4, bSaveMap = true; end
if nargin < 5, bLogScale = true; end
if nargin < 6, bRAMSafeMode = false; end

if ~strcmp(folderPath(end), filesep)
    folderPath = [folderPath filesep];
end

channel = lower(channel);
datFile = [folderPath channel '.dat'];
matFile = [folderPath channel '.mat'];

if ~exist(datFile, 'file')
    error('%s file is missing.', datFile);
end

Infos = load(matFile);
Ny = Infos.datSize(1);
Nx = Infos.datSize(2);
Nt = Infos.datLength;

% Metadata template for map output
md = Infos;
md.datLength = 1;
md.Freq = 0;
md.dim_names = {'Y','X'};

if bRAMSafeMode
    % ---------------------------------------------------------------------
    % LOW-RAM MODE
    % ---------------------------------------------------------------------
    outFilename = [folderPath 'SPECKLEMAP.dat'];
    frameOut = zeros(Ny, Nx, Infos.Datatype);
    mData = zeros(Ny, Nx, Infos.Datatype);

    fidIn = fopen(datFile, 'r');
    cIn = onCleanup(@() safeFclose(fidIn)); %#ok<NASGU>

    % Pass 1: temporal mean
    disp('Pass 1/2 - Calculating temporal mean...')

    bytesPerFrame = Ny * Nx * getByteSize(Infos.Datatype);
    totalBytes = Nx * Ny * Nt * getByteSize(Infos.Datatype);

    nChunks = calculateMaxChunkSize(totalBytes, 1, .1);
    chunkT  = ceil(Nt / nChunks);

    lastPct = -1;
    fprintf('0%% ');

    for c = 1:nChunks
        tStart  = (c-1) * chunkT + 1;
        tEnd    = min(tStart + chunkT - 1, Nt);
        nFrames = tEnd - tStart + 1;

        fseek(fidIn, (tStart-1) * bytesPerFrame, 'bof');
        slab = fread(fidIn, Ny * Nx * nFrames, ['*' Infos.Datatype]);
        slab = reshape(slab, Ny, Nx, nFrames);

        % Use omitnan during accumulation
        mData = mData + sum(slab, 3, 'omitnan');

        pct = floor(100 * c / nChunks);
        if pct ~= lastPct
            fprintf('%d%% ', pct);
            lastPct = pct;
        end

        clear slab
    end

    mData = mData ./ Nt;

    fprintf('\nPass 2/2 - Calculating Speckle Contrast (%s algorithm)...\n', sType)

    switch lower(sType)
        case 'spatial'
            Kernel = single(fspecial('disk', 2) > 0);

            nChunks = calculateMaxChunkSize(totalBytes, 10, .1);
            chunkT  = ceil(Nt / nChunks);

            lastPct = -1;
            fprintf('0%% ');

            for c = 1:nChunks
                tStart  = (c-1) * chunkT + 1;
                tEnd    = min(tStart + chunkT - 1, Nt);
                nFrames = tEnd - tStart + 1;

                fseek(fidIn, (tStart-1) * bytesPerFrame, 'bof');
                frameBlock = fread(fidIn, Ny * Nx * nFrames, ['*' Infos.Datatype]);
                frameBlock = reshape(frameBlock, Ny, Nx, nFrames);

                % Use omitnan normalization semantics
                frameBlock = frameBlock ./ mData;

                frameBlock = stdfilt(frameBlock, Kernel);
                % frameBlock = remOutlier(frameBlock);
                frameBlock = sum(frameBlock, 3, 'omitnan');
                frameOut = frameOut + frameBlock;

                pct = floor(100 * c / nChunks);
                if pct ~= lastPct
                    fprintf('%d%% ', pct);
                    lastPct = pct;
                end

                clear frameBlock
            end

            frameOut = frameOut ./ Nt;

        case 'temporal'
            Kernel = ones(1, 1, 5, 'single');

            nChunks = calculateMaxChunkSize(Nx * Ny * Nt * getByteSize(Infos.Datatype), 12, .15);
            chunkX  = ceil(Nx / nChunks);

            lastPct = -1;
            fprintf('0%% ');

            for c = 1:nChunks
                xStart = (c-1) * chunkX + 1;
                xEnd   = min(xStart + chunkX - 1, Nx);
                xIdx   = xStart:xEnd;

                slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, Infos.Datatype);
                slab = slab ./ mean(slab, 3, 'omitnan');
                slab = stdfilt(slab, Kernel);
                % slab = remOutlier(slab);
                frameOut(:, xIdx) = mean(slab, 3, 'omitnan');

                pct = floor(100 * c / nChunks);
                if pct ~= lastPct
                    fprintf('%d%% ', pct);
                    lastPct = pct;
                end

                clear slab
            end

        otherwise
            error('Invalid sType. Use "spatial" or "temporal".');
    end

    if bLogScale
        frameOut = -log10(frameOut);
    end

    fclose(fidIn);

    md = genMetaData(frameOut, {'Y','X'}, Infos);
    md.datFile = outFilename;
    save2Dat(outFilename, single(frameOut), md);

    DatOut = outFilename;

else
    % ---------------------------------------------------------------------
    % STANDARD MODE
    % ---------------------------------------------------------------------
    fid = fopen(datFile, 'r');
    dat = fread(fid, inf, '*single');
    fclose(fid);

    dat = reshape(dat, Ny, Nx, Nt);
    dat = dat ./ mean(dat, 3, 'omitnan');

    disp('Mapping computation...');

    switch lower(sType)
        case 'spatial'
            Kernel = single(fspecial('disk', 2) > 0);
        case 'temporal'
            Kernel = ones(1, 1, 5, 'single');
        otherwise
            error('Invalid sType. Use "spatial" or "temporal".');
    end

    DatOut = stdfilt(dat, Kernel);
    % DatOut = remOutlier(DatOut);
    DatOut = mean(DatOut, 3, 'omitnan');

    if bLogScale
        DatOut = -log10(DatOut);
    end

    DatOut = single(DatOut);
end

% -------------------------------------------------------------------------
% Save TIFF map if requested
% -------------------------------------------------------------------------
if bSaveMap
    if bRAMSafeMode
        mapHeight = Infos.datSize(1);
        mapWidth  = Infos.datSize(2);

        tmpFid = fopen(DatOut, 'r');
        cTmp = onCleanup(@() safeFclose(tmpFid)); %#ok<NASGU>
        mapData = fread(tmpFid, mapHeight * mapWidth, '*single');
        mapData = reshape(mapData, mapHeight, mapWidth);
    else
        mapHeight = size(DatOut, 1);
        mapWidth  = size(DatOut, 2);
        mapData   = DatOut;
    end

    obj = Tiff(fullfile(folderPath, 'std_speckle.tiff'), 'w');
    setTag(obj, 'ImageWidth', mapWidth);
    setTag(obj, 'ImageLength', mapHeight);
    setTag(obj, 'Photometric', Tiff.Photometric.MinIsBlack);
    setTag(obj, 'SampleFormat', Tiff.SampleFormat.IEEEFP);
    setTag(obj, 'BitsPerSample', 32);
    setTag(obj, 'SamplesPerPixel', 1);
    setTag(obj, 'Compression', Tiff.Compression.None);
    setTag(obj, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    write(obj, mapData);
    close(obj);
end

disp('Done');
end

function data = remOutlier(data)
% REMOUTLIER Replace values above the 99th percentile with the 99th percentile.
dataSort = data(~isnan(data));
dataSort = sort(dataSort(:));
idx99 = floor(0.99 * length(dataSort));
val99th = dataSort(idx99);
data(data > val99th) = val99th;
end