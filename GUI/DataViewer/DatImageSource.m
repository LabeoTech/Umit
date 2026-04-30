classdef DatImageSource < handle
%DATIMAGESOURCE RAM-safe direct-read backend for continuous [Y,X,T] .dat files.
%
%   src = DatImageSource(filePath)
%   src = DatImageSource(filePath, Name, Value)
%
%   This class provides a DataViewer-oriented access layer for legacy/raw
%   .dat files that are stored as single-precision MATLAB-order [Y,X,T]
%   binary arrays. It intentionally avoids MEMMAPFILE. Full frames and
%   temporal-cache blocks are read with explicit FSEEK/FREAD calls.
%
%   Supported layout:
%       dim_names = {'Y','X','T'}
%
%   Unsupported layout:
%       Any .dat file with an event dimension, for example
%       dim_names = {'Y','X','T','E'}. Event-split image data should be
%       stored and loaded as a .umt image structure.
%
%   Name-Value options:
%       cacheRAMFraction - Fraction of conservative usable RAM assigned to
%                          the temporal XY cache. Default: 0.25.
%       RAMoverhead      - Fraction of total RAM reserved for OS/MATLAB and
%                          temporary allocations. Default: 0.30.
%       maxUsableRAMFrac - Maximum fraction of total physical RAM considered
%                          usable by the app. Default: 0.70.
%       maxCacheBytes    - Optional hard cap for cache bytes. Default: [].
%       cacheMode        - 'auto' or 'locked'. Default: 'auto'.
%
%   Main methods:
%       getSize              - Return [Ny Nx Nt 1].
%       getFrame             - Read one full [Y,X] frame.
%       getFrameBlock        - Read one spatial block from one frame.
%       getPixelTrace        - Return full-T pixel trace from the cache.
%       getROIMeanTraceMatrix - Return ROI mean traces for normal/event frames.
%       updateCacheAround    - Rebuild cache centered on one pixel.
%       isInsideCache        - Test if one pixel is inside the cache.
%       getCacheRectangle    - Return rectangle position for overlay.
%       hasPartialTemporalCache - Return true when cache covers only part
%                              of the image.
%
%   Notes:
%       - loadMetaData must be available on the MATLAB path.
%       - queryRAM must be available on the MATLAB path.
%       - Data are assumed to be written in MATLAB column-major order.

    properties (SetAccess = private)
        FilePath char = ''
        Info struct = struct()

        Ny double = 0
        Nx double = 0
        Nt double = 0

        Precision char = 'single'
        BytesPerSample double = 4
        FrameRateHz double = NaN
    end

    properties
        CacheRAMFraction double = 0.25
        RAMoverhead double = 0.30
        MaxUsableRAMFraction double = 0.70
        MaxCacheBytes double = []

        CacheMode char = 'auto'  % 'auto' or 'locked'
    end

    properties (SetAccess = private)
        CacheData = []
        CacheYRange double = []
        CacheXRange double = []
        LastCacheBytes double = 0
    end

    methods
        function obj = DatImageSource(filePath, varargin)
            %DATIMAGESOURCE Construct a .dat image source.
            %
            %   obj = DatImageSource(filePath)
            %   obj = DatImageSource(filePath, 'cacheRAMFraction', 0.25)
            %
            %   The constructor validates that the file is a continuous
            %   [Y,X,T] .dat file. Legacy event-split .dat files are rejected
            %   before any 3D reshape or cache access is attempted.

            p = inputParser;
            p.FunctionName = 'DatImageSource';

            addRequired(p, 'filePath', @(x) ischar(x) || isstring(x));
            addParameter(p, 'cacheRAMFraction', 0.25, ...
                @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
            addParameter(p, 'RAMoverhead', 0.30, ...
                @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 1);
            addParameter(p, 'maxUsableRAMFrac', 0.70, ...
                @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
            addParameter(p, 'maxCacheBytes', [], ...
                @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
            addParameter(p, 'cacheMode', 'auto', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));

            parse(p, filePath, varargin{:});

            filePath = char(string(p.Results.filePath));
            if isempty(fileparts(filePath))
                filePath = fullfile(pwd, filePath);
            end

            if ~isfile(filePath)
                error('DatImageSource:FileNotFound', ...
                    'File not found: "%s".', filePath);
            end

            obj.FilePath = filePath;
            obj.CacheRAMFraction = p.Results.cacheRAMFraction;
            obj.RAMoverhead = p.Results.RAMoverhead;
            obj.MaxUsableRAMFraction = p.Results.maxUsableRAMFrac;
            obj.MaxCacheBytes = p.Results.maxCacheBytes;
            obj.setCacheMode(p.Results.cacheMode);

            obj.Info = loadMetaData(obj.FilePath);
            obj.validateContinuousDatLayout(obj.Info);

            obj.Ny = double(obj.Info.Height);
            obj.Nx = double(obj.Info.Width);
            obj.Nt = double(obj.Info.Length);

            obj.Precision = char(string(obj.Info.Datatype));
            obj.BytesPerSample = obj.getByteSize(obj.Precision);

            if isfield(obj.Info, 'FrameRateHz') && ~isempty(obj.Info.FrameRateHz)
                obj.FrameRateHz = double(obj.Info.FrameRateHz);
            elseif isfield(obj.Info, 'Freq') && ~isempty(obj.Info.Freq)
                obj.FrameRateHz = double(obj.Info.Freq);
            end

            obj.validateFileSize();

            % Initial cache: largest safe cache centered in the image.
            obj.updateCacheAround(round(obj.Ny / 2), round(obj.Nx / 2));
        end

        function sz = getSize(obj)
            %GETSIZE Return canonical viewer size [Y X T E].
            %
            %   sz = obj.getSize()
            %
            %   For .dat files, E is always singleton because only
            %   continuous [Y,X,T] .dat files are supported.

            sz = [obj.Ny, obj.Nx, obj.Nt, 1];
        end

        function info = getInfo(obj)
            %GETINFO Return flat metadata structure.

            info = obj.Info;
        end

        function labels = getLabels(obj)
            %GETLABELS Return empty labels struct for .dat backend.

            labels = struct();
        end

        function eventInfo = getEventInfo(obj)
            %GETEVENTINFO Return empty eventInfo struct for .dat backend.

            eventInfo = struct();
        end

        function frame = getFrame(obj, tIdx, varargin) %#ok<INUSD>
            %GETFRAME Read one full frame from disk.
            %
            %   frame = obj.getFrame(tIdx)
            %
            %   Output:
            %       frame - Numeric matrix of size [Y,X].

            obj.validateFrameIndex(tIdx);

            if obj.isFullImageCached()
                frame = obj.CacheData(:, :, tIdx);
                return
            end

            fid = fopen(obj.FilePath, 'r');
            if fid < 0
                error('DatImageSource:FileOpenFailed', ...
                    'Could not open file: "%s".', obj.FilePath);
            end
            cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

            offset = obj.frameOffsetBytes(tIdx);
            obj.safeFseek(fid, offset, ...
                sprintf('Failed to seek to frame %d.', tIdx));

            frame = fread(fid, [obj.Ny, obj.Nx], ['*' obj.Precision]);

            if numel(frame) ~= obj.Ny * obj.Nx
                error('DatImageSource:UnexpectedEOF', ...
                    'Unexpected end of file while reading frame %d.', tIdx);
            end
        end

        function block = getFrameBlock(obj, tIdx, yRange, xRange)
            %GETFRAMEBLOCK Read one spatial block from one frame.
            %
            %   block = obj.getFrameBlock(tIdx, yRange, xRange)
            %
            %   Inputs:
            %       tIdx   - Frame index.
            %       yRange - Contiguous Y indices.
            %       xRange - Contiguous X indices.
            %
            %   Output:
            %       block  - Numeric matrix [numel(yRange), numel(xRange)].

            obj.validateFrameIndex(tIdx);
            yRange = obj.validateContiguousIndexRange(yRange, obj.Ny, 'Y');
            xRange = obj.validateContiguousIndexRange(xRange, obj.Nx, 'X');

            if obj.isInsideCachedBlock(yRange, xRange)
                yLocal = yRange - obj.CacheYRange(1) + 1;
                xLocal = xRange - obj.CacheXRange(1) + 1;
                block = obj.CacheData(yLocal, xLocal, tIdx);
                return
            end

            fid = fopen(obj.FilePath, 'r');
            if fid < 0
                error('DatImageSource:FileOpenFailed', ...
                    'Could not open file: "%s".', obj.FilePath);
            end
            cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

            nX = numel(xRange);
            x0 = xRange(1);

            offset = obj.frameOffsetBytes(tIdx) + ...
                double(x0 - 1) * obj.Ny * obj.BytesPerSample;

            obj.safeFseek(fid, offset, ...
                sprintf('Failed to seek to frame %d, x=%d.', tIdx, x0));

            slab = fread(fid, [obj.Ny, nX], ['*' obj.Precision]);

            if numel(slab) ~= obj.Ny * nX
                error('DatImageSource:UnexpectedEOF', ...
                    'Unexpected end of file while reading frame block.');
            end

            block = slab(yRange, :);
        end

        function [trace, status] = getPixelTrace(obj, y, x)
            %GETPIXELTRACE Return full temporal trace for one pixel.
            %
            %   trace = obj.getPixelTrace(y, x)
            %   [trace, status] = obj.getPixelTrace(y, x)
            %
            %   If the pixel is inside the temporal cache, the trace is read
            %   directly from memory.
            %
            %   If the pixel is outside the cache:
            %       - cacheMode='auto'   : cache is rebuilt around the pixel.
            %       - cacheMode='locked' : trace is [] and status is
            %                              'outside_locked_cache'.
            %
            %   Status values:
            %       'ok'
            %       'cache_rebuilt'
            %       'outside_locked_cache'

            obj.validatePixelIndex(y, x);

            status = 'ok';

            if ~obj.isInsideCache(y, x)
                if obj.isCacheLocked()
                    trace = [];
                    status = 'outside_locked_cache';
                    return
                end

                obj.updateCacheAround(y, x);
                status = 'cache_rebuilt';
            end

            yLocal = y - obj.CacheYRange(1) + 1;
            xLocal = x - obj.CacheXRange(1) + 1;

            trace = squeeze(obj.CacheData(yLocal, xLocal, :));
            trace = trace(:);
        end

        function tf = isInsideCache(obj, y, x)
            %ISINSIDECACHE Return true if pixel is inside temporal cache.

            tf = ~isempty(obj.CacheData) && ...
                y >= obj.CacheYRange(1) && y <= obj.CacheYRange(end) && ...
                x >= obj.CacheXRange(1) && x <= obj.CacheXRange(end);
        end

        function updateCacheAround(obj, yCenter, xCenter)
            %UPDATECACHEAROUND Rebuild temporal XY cache around one pixel.
            %
            %   obj.updateCacheAround(yCenter, xCenter)
            %
            %   The cache region is the largest [Ycache,Xcache,T] block that
            %   fits the configured RAM budget while approximately preserving
            %   the image aspect ratio. The cache is replaced, not expanded.

            obj.validatePixelIndex(yCenter, xCenter);

            [nYCache, nXCache] = obj.estimateCacheSize();
            [yRange, xRange] = obj.centerRanges(yCenter, xCenter, nYCache, nXCache);

            newCache = obj.readTemporalBlock(yRange, xRange);

            obj.CacheData = newCache;
            obj.CacheYRange = yRange;
            obj.CacheXRange = xRange;
            obj.LastCacheBytes = numel(newCache) * obj.BytesPerSample;
        end

        function rect = getCacheRectangle(obj)
            %GETCACHERECTANGLE Return cache overlay rectangle position.
            %
            %   rect = obj.getCacheRectangle()
            %
            %   Output follows MATLAB rectangle Position convention:
            %       [x y width height]
            %
            %   The position is shifted by -0.5 so the rectangle outlines
            %   pixel boundaries instead of pixel centers.

            if isempty(obj.CacheYRange) || isempty(obj.CacheXRange)
                rect = [];
                return
            end

            rect = [ ...
                obj.CacheXRange(1) - 0.5, ...
                obj.CacheYRange(1) - 0.5, ...
                numel(obj.CacheXRange), ...
                numel(obj.CacheYRange)];
        end


        function [traceMatrix, readMode] = getROIMeanTraceMatrix(obj, roiMasks, frameIdx)
            %GETROIMEANTRACEMATRIX Return mean intensity traces for ROI masks.
            %
            %   [traceMatrix, readMode] = obj.getROIMeanTraceMatrix(roiMasks, frameIdx)
            %
            %   Inputs:
            %       roiMasks - ROI masks as either:
            %                  - cell array of logical [Y,X] masks
            %                  - logical/numeric [Y,X,nROI] stack
            %                  - logical/numeric [Y,X] single mask
            %       frameIdx - Source frame indices. The shape controls the output:
            %                  - vector [1,nFrames] or [nFrames,1]
            %                    returns traceMatrix [nROI,nFrames]
            %                  - matrix [nTrials,nFrames]
            %                    returns traceMatrix [nROI,nTrials,nFrames]
            %
            %   Output:
            %       traceMatrix - ROI spatial means for requested frames.
            %       readMode    - 'cache', 'direct_full_frame', or 'none'.
            %
            %   Notes:
            %       - If all ROI masks are fully inside the active temporal cache,
            %         traces are computed from CacheData.
            %       - If any ROI mask is outside the active temporal cache, traces
            %         are computed by RAM-bounded full-frame direct reads.
            %       - Direct reads do not rebuild or move the temporal cache.

            if nargin < 3 || isempty(frameIdx)
                frameIdx = 1:obj.Nt;
            end

            [maskStack, emptyROI] = obj.normalizeROIMasks(roiMasks);
            frameIdx = double(frameIdx);

            if isempty(maskStack)
                traceMatrix = [];
                readMode = 'none';
                return
            end

            obj.validateROIFrameIndex(frameIdx);

            nROI = size(maskStack, 3);
            bEventMatrix = ~isvector(frameIdx) && ismatrix(frameIdx);

            if bEventMatrix
                traceMatrix = nan(nROI, size(frameIdx, 1), size(frameIdx, 2));
            else
                traceMatrix = nan(nROI, numel(frameIdx));
                frameIdx = frameIdx(:).';
            end

            if ~any(isfinite(frameIdx(:)))
                readMode = 'none';
                return
            end

            if obj.areROIMasksInsideCache(maskStack)
                traceMatrix = obj.computeROIMeanTraceMatrixFromCache( ...
                    maskStack, emptyROI, frameIdx, bEventMatrix);
                readMode = 'cache';
            else
                traceMatrix = obj.computeROIMeanTraceMatrixDirect( ...
                    maskStack, emptyROI, frameIdx, bEventMatrix);
                readMode = 'direct_full_frame';
            end
        end

        function setCacheMode(obj, mode)
            %SETCACHEMODE Set temporal cache mode.
            %
            %   obj.setCacheMode('auto')
            %   obj.setCacheMode('locked')

            mode = lower(char(string(mode)));

            if ~ismember(mode, {'auto', 'locked'})
                error('DatImageSource:InvalidCacheMode', ...
                    'Cache mode must be "auto" or "locked".');
            end

            obj.CacheMode = mode;
        end

        function setCacheLocked(obj, tf)
            %SETCACHELOCKED Convenience setter for cache lock state.

            validateattributes(tf, {'logical'}, {'scalar'}, ...
                'setCacheLocked', 'tf');

            if tf
                obj.CacheMode = 'locked';
            else
                obj.CacheMode = 'auto';
            end
        end

        function tf = isCacheLocked(obj)
            %ISCACHELOCKED Return true when cache mode is locked.

            tf = strcmpi(obj.CacheMode, 'locked');
        end

        function tf = hasPartialTemporalCache(obj)
            %HASPARTIALTEMPORALCACHE Return true for partial XYT cache mode.
            %
            %   tf = obj.hasPartialTemporalCache()
            %
            %   Returns true when the active temporal cache covers only part
            %   of the image. Returns false when the cache is empty or when
            %   the whole image is cached in memory.

            if isempty(obj.CacheData)
                tf = false;
                return
            end

            tf = ~obj.isFullImageCached();
        end

        function txt = getCacheStatusText(obj)
            %GETCACHESTATUSTEXT Return compact user-facing cache status.

            if isempty(obj.CacheData)
                txt = 'Cache: empty';
                return
            end

            txt = sprintf('Cache: %s | %d x %d x %d', ...
                obj.CacheMode, ...
                numel(obj.CacheYRange), ...
                numel(obj.CacheXRange), ...
                obj.Nt);
        end
    end

    methods (Access = private)

        function [maskStack, emptyROI] = normalizeROIMasks(obj, roiMasks)
            %NORMALIZEROIMASKS Validate ROI masks and return [Y,X,nROI] stack.

            if isempty(roiMasks)
                maskStack = false(obj.Ny, obj.Nx, 0);
                emptyROI = false(0, 1);
                return
            end

            if iscell(roiMasks)
                nROI = numel(roiMasks);
                maskStack = false(obj.Ny, obj.Nx, nROI);

                for iROI = 1:nROI
                    thisMask = roiMasks{iROI};
                    obj.validateOneROIMask(thisMask, iROI);
                    maskStack(:, :, iROI) = logical(thisMask);
                end
            else
                if ~isnumeric(roiMasks) && ~islogical(roiMasks)
                    error('DatImageSource:InvalidROIMasks', ...
                        'roiMasks must be a cell array or a numeric/logical mask array.');
                end

                if ismatrix(roiMasks)
                    obj.validateOneROIMask(roiMasks, 1);
                    maskStack = logical(roiMasks);
                    maskStack = reshape(maskStack, obj.Ny, obj.Nx, 1);
                elseif ndims(roiMasks) == 3
                    if size(roiMasks, 1) ~= obj.Ny || size(roiMasks, 2) ~= obj.Nx
                        error('DatImageSource:InvalidROIMaskSize', ...
                            'ROI mask stack must have size [Y,X,nROI]=[%d,%d,nROI].', ...
                            obj.Ny, obj.Nx);
                    end
                    maskStack = logical(roiMasks);
                else
                    error('DatImageSource:InvalidROIMasks', ...
                        'roiMasks must be [Y,X], [Y,X,nROI], or a cell array of [Y,X] masks.');
                end
            end

            nROI = size(maskStack, 3);
            emptyROI = false(nROI, 1);

            for iROI = 1:nROI
                emptyROI(iROI) = ~any(maskStack(:, :, iROI), 'all');
            end
        end

        function validateOneROIMask(obj, mask, roiIdx)
            %VALIDATEONEROIMASK Validate one ROI mask.

            if (~isnumeric(mask) && ~islogical(mask)) || ~ismatrix(mask)
                error('DatImageSource:InvalidROIMask', ...
                    'ROI mask %d must be a numeric or logical 2D array.', roiIdx);
            end

            if ~isequal(size(mask), [obj.Ny, obj.Nx])
                error('DatImageSource:InvalidROIMaskSize', ...
                    'ROI mask %d must have size [Y,X]=[%d,%d].', ...
                    roiIdx, obj.Ny, obj.Nx);
            end
        end

        function validateROIFrameIndex(obj, frameIdx)
            %VALIDATEROIFRAMEINDEX Validate ROI trace frame indices.

            if isempty(frameIdx)
                return
            end

            if ~isnumeric(frameIdx) || ~ismatrix(frameIdx)
                error('DatImageSource:InvalidROIFrameIndex', ...
                    'frameIdx must be a numeric vector or numeric matrix.');
            end

            finiteFrameIdx = frameIdx(isfinite(frameIdx));

            if isempty(finiteFrameIdx)
                return
            end

            if any(mod(finiteFrameIdx, 1) ~= 0) || ...
                    any(finiteFrameIdx < 1) || any(finiteFrameIdx > obj.Nt)
                error('DatImageSource:InvalidROIFrameIndex', ...
                    'frameIdx contains indices outside the valid range [1,%d].', obj.Nt);
            end
        end

        function tf = areROIMasksInsideCache(obj, maskStack)
            %AREROIMASKSINSIDECACHE True if all non-empty ROIs fit in cache.

            tf = false;

            if isempty(obj.CacheData) || isempty(obj.CacheYRange) || isempty(obj.CacheXRange)
                return
            end

            nROI = size(maskStack, 3);

            for iROI = 1:nROI
                thisMask = maskStack(:, :, iROI);

                if ~any(thisMask(:))
                    continue
                end

                yRows = find(any(thisMask, 2));
                xCols = find(any(thisMask, 1));

                if isempty(yRows) || isempty(xCols)
                    continue
                end

                if yRows(1) < obj.CacheYRange(1) || yRows(end) > obj.CacheYRange(end) || ...
                        xCols(1) < obj.CacheXRange(1) || xCols(end) > obj.CacheXRange(end)
                    return
                end
            end

            tf = true;
        end

        function traceMatrix = computeROIMeanTraceMatrixFromCache(obj, maskStack, emptyROI, frameIdx, bEventMatrix)
            %COMPUTEROIMEANTRACEMATRIXFROMCACHE Compute ROI traces from CacheData.

            yRange = obj.CacheYRange;
            xRange = obj.CacheXRange;
            [roiWeights, emptyLocalROI] = obj.makeROIWeightMatrix(maskStack, yRange, xRange);
            emptyRows = emptyROI | emptyLocalROI;

            nROI = size(maskStack, 3);
            nPix = numel(yRange) * numel(xRange);

            if bEventMatrix
                nTrials = size(frameIdx, 1);
                nFrames = size(frameIdx, 2);
                traceMatrix = nan(nROI, nTrials, nFrames);

                for iTrial = 1:nTrials
                    trialFrames = frameIdx(iTrial, :);
                    validPos = find(isfinite(trialFrames));

                    if isempty(validPos)
                        continue
                    end

                    cacheBlock = obj.CacheData(:, :, trialFrames(validPos));
                    cacheBlock2D = reshape(cacheBlock, nPix, numel(validPos));
                    traceChunk = roiWeights * double(cacheBlock2D);
                    traceChunk(emptyRows, :) = NaN;
                    traceMatrix(:, iTrial, validPos) = traceChunk;
                end
            else
                frameIdx = frameIdx(:).';
                nFrames = numel(frameIdx);
                traceMatrix = nan(nROI, nFrames);
                validPos = find(isfinite(frameIdx));

                if ~isempty(validPos)
                    cacheBlock = obj.CacheData(:, :, frameIdx(validPos));
                    cacheBlock2D = reshape(cacheBlock, nPix, numel(validPos));
                    traceChunk = roiWeights * double(cacheBlock2D);
                    traceChunk(emptyRows, :) = NaN;
                    traceMatrix(:, validPos) = traceChunk;
                end
            end
        end

        function traceMatrix = computeROIMeanTraceMatrixDirect(obj, maskStack, emptyROI, frameIdx, bEventMatrix)
            %COMPUTEROIMEANTRACEMATRIXDIRECT Compute ROI traces from full-frame chunks.

            [roiWeights, emptyLocalROI] = obj.makeROIWeightMatrix(maskStack, 1:obj.Ny, 1:obj.Nx);
            emptyRows = emptyROI | emptyLocalROI;

            nROI = size(maskStack, 3);
            nPix = obj.Ny * obj.Nx;

            fid = fopen(obj.FilePath, 'r');
            if fid < 0
                error('DatImageSource:FileOpenFailed', ...
                    'Could not open file: "%s".', obj.FilePath);
            end
            cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

            if bEventMatrix
                nTrials = size(frameIdx, 1);
                nFrames = size(frameIdx, 2);
                traceMatrix = nan(nROI, nTrials, nFrames);

                for iTrial = 1:nTrials
                    trialFrames = frameIdx(iTrial, :);
                    validPos = find(isfinite(trialFrames));

                    if isempty(validPos)
                        continue
                    end

                    framesPerChunk = obj.estimateROIFramesPerChunk(numel(validPos));

                    for firstPos = 1:framesPerChunk:numel(validPos)
                        lastPos = min(firstPos + framesPerChunk - 1, numel(validPos));
                        chunkPos = validPos(firstPos:lastPos);
                        chunkFrames = trialFrames(chunkPos);

                        frameBlock = obj.readFullFrameListFromOpenFile(fid, chunkFrames);
                        frameBlock2D = reshape(frameBlock, nPix, numel(chunkFrames));
                        traceChunk = roiWeights * double(frameBlock2D);
                        traceChunk(emptyRows, :) = NaN;
                        traceMatrix(:, iTrial, chunkPos) = traceChunk;
                    end
                end
            else
                frameIdx = frameIdx(:).';
                nFrames = numel(frameIdx);
                traceMatrix = nan(nROI, nFrames);
                validPos = find(isfinite(frameIdx));

                if isempty(validPos)
                    return
                end

                framesPerChunk = obj.estimateROIFramesPerChunk(numel(validPos));

                for firstPos = 1:framesPerChunk:numel(validPos)
                    lastPos = min(firstPos + framesPerChunk - 1, numel(validPos));
                    chunkPos = validPos(firstPos:lastPos);
                    chunkFrames = frameIdx(chunkPos);

                    frameBlock = obj.readFullFrameListFromOpenFile(fid, chunkFrames);
                    frameBlock2D = reshape(frameBlock, nPix, numel(chunkFrames));
                    traceChunk = roiWeights * double(frameBlock2D);
                    traceChunk(emptyRows, :) = NaN;
                    traceMatrix(:, chunkPos) = traceChunk;
                end
            end
        end

        function [roiWeights, emptyROI] = makeROIWeightMatrix(obj, maskStack, yRange, xRange)
            %MAKEROIWEIGHTMATRIX Build sparse ROI averaging matrix.

            yRange = double(yRange(:).');
            xRange = double(xRange(:).');

            nROI = size(maskStack, 3);
            nPix = numel(yRange) * numel(xRange);

            rowIdx = [];
            colIdx = [];
            values = [];
            emptyROI = false(nROI, 1);

            for iROI = 1:nROI
                thisMask = maskStack(yRange, xRange, iROI);
                pixelIdx = find(thisMask(:));

                if isempty(pixelIdx)
                    emptyROI(iROI) = true;
                    continue
                end

                nMaskPixels = numel(pixelIdx);
                rowIdx = [rowIdx; repmat(iROI, nMaskPixels, 1)]; %#ok<AGROW>
                colIdx = [colIdx; pixelIdx(:)]; %#ok<AGROW>
                values = [values; repmat(1 ./ nMaskPixels, nMaskPixels, 1)]; %#ok<AGROW>
            end

            roiWeights = sparse(rowIdx, colIdx, values, nROI, nPix);
        end

        function framesPerChunk = estimateROIFramesPerChunk(obj, nFramesRequested)
            %ESTIMATEROIFRAMESPERCHUNK Estimate RAM-safe full-frame chunk size.

            nFramesRequested = max(0, round(double(nFramesRequested)));

            if nFramesRequested <= 0
                framesPerChunk = 0;
                return
            end

            requiredBytes = double(obj.Ny) * double(obj.Nx) * ...
                double(nFramesRequested) * double(obj.BytesPerSample);

            sizeFactor = 2.5;

            try
                nChunks = calculateMaxChunkSize(requiredBytes, sizeFactor, obj.RAMoverhead);
            catch
                nChunks = nFramesRequested;
            end

            if isempty(nChunks) || ~isfinite(nChunks) || nChunks < 1
                nChunks = nFramesRequested;
            end

            nChunks = min(nFramesRequested, max(1, ceil(nChunks)));
            framesPerChunk = max(1, ceil(nFramesRequested ./ nChunks));
        end

        function frameBlock = readFullFrameListFromOpenFile(obj, fid, frameList)
            %READFULLFRAMELISTFROMOPENFILE Read full frames from an open file.

            frameList = double(frameList(:).');

            if isempty(frameList)
                frameBlock = zeros(obj.Ny, obj.Nx, 0, obj.Precision);
                return
            end

            if any(~isfinite(frameList)) || any(mod(frameList, 1) ~= 0)
                error('DatImageSource:InvalidFrameList', ...
                    'Frame list contains invalid frame indices.');
            end

            if numel(frameList) == 1 || all(diff(frameList) == 1)
                frameBlock = obj.readFullFrameChunkFromOpenFile(fid, frameList(1), numel(frameList));
                return
            end

            frameBlock = zeros(obj.Ny, obj.Nx, numel(frameList), obj.Precision);

            for iFrame = 1:numel(frameList)
                frameBlock(:, :, iFrame) = obj.readFullFrameChunkFromOpenFile(fid, frameList(iFrame), 1);
            end
        end

        function frameBlock = readFullFrameChunkFromOpenFile(obj, fid, startFrame, nFrames)
            %READFULLFRAMECHUNKFROMOPENFILE Read consecutive full frames from an open file.

            startFrame = double(startFrame);
            nFrames = double(nFrames);

            obj.validateFrameIndex(startFrame);
            obj.validateFrameIndex(startFrame + nFrames - 1);

            offset = obj.frameOffsetBytes(startFrame);
            obj.safeFseek(fid, offset, ...
                sprintf('Failed to seek to frame %d.', startFrame));

            nElements = obj.Ny * obj.Nx * nFrames;
            [raw, count] = fread(fid, nElements, ['*' obj.Precision]);

            if count ~= nElements
                error('DatImageSource:UnexpectedEOF', ...
                    'Unexpected end of file while reading full-frame chunk.');
            end

            frameBlock = reshape(raw, obj.Ny, obj.Nx, nFrames);
        end

        function validateContinuousDatLayout(obj, Info) %#ok<INUSL>
            %VALIDATECONTINUOUSDATLAYOUT Reject non-YXT .dat files.
            %
            %   This method deliberately fails fast for legacy event-split
            %   .dat files. Those files should be converted to .umt image
            %   structures before being opened in the viewer.

            if ~isfield(Info, 'FileType') || ~strcmpi(char(string(Info.FileType)), '.dat')
                error('DatImageSource:InvalidFileType', ...
                    'DatImageSource can only open .dat files.');
            end

            if ~isfield(Info, 'dim_names') || isempty(Info.dim_names)
                error('DatImageSource:MissingDimNames', ...
                    '.dat metadata must contain dim_names.');
            end

            dimNames = upper(cellstr(string(Info.dim_names(:).')));
            expected = {'Y', 'X', 'T'};

            if any(strcmp(dimNames, 'E'))
                error('DatImageSource:EventSplitDatUnsupported', ...
                    ['This .dat file declares an event dimension and cannot be ' ...
                     'opened by DatImageSource. Legacy event-split .dat files ' ...
                     'are no longer supported by the DataViewer backend. ' ...
                     'Convert this dataset to a .umt image structure with ' ...
                     'dimNames {''Y'',''X'',''T'',''E''} and shared eventInfo.']);
            end

            if numel(dimNames) ~= 3 || ~all(strcmp(dimNames, expected))
                error('DatImageSource:UnsupportedDatLayout', ...
                    ['Unsupported .dat layout: {%s}. DatImageSource only ' ...
                     'supports continuous files with dim_names={''Y'',''X'',''T''}.'], ...
                    strjoin(dimNames, ', '));
            end

            if ~isfield(Info, 'Datatype') || isempty(Info.Datatype)
                error('DatImageSource:MissingDatatype', ...
                    '.dat metadata must contain Datatype.');
            end

            if ~strcmpi(char(string(Info.Datatype)), 'single')
                error('DatImageSource:UnsupportedDatatype', ...
                    'Only single-precision .dat files are currently supported.');
            end

            if ~isfield(Info, 'Height') || ~isfield(Info, 'Width') || ...
                    ~isfield(Info, 'Length')
                error('DatImageSource:MissingCoreMetadata', ...
                    '.dat metadata must contain Height, Width, and Length.');
            end

            validateattributes(double(Info.Height), {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
                'DatImageSource', 'Height');
            validateattributes(double(Info.Width), {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
                'DatImageSource', 'Width');
            validateattributes(double(Info.Length), {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
                'DatImageSource', 'Length');

            if isfield(Info, 'datSize') && ~isempty(Info.datSize)
                datSize = double(Info.datSize(:).');
                expectedYX = [double(Info.Height), double(Info.Width)];
                expectedYXT = [double(Info.Height), double(Info.Width), double(Info.Length)];

                if numel(datSize) == 2
                    if ~isequal(datSize, expectedYX)
                        error('DatImageSource:InvalidDatSize', ...
                            'datSize is incompatible with Height and Width.');
                    end
                elseif numel(datSize) == 3
                    if ~isequal(datSize, expectedYXT)
                        error('DatImageSource:InvalidDatSize', ...
                            'datSize is incompatible with Height, Width, and Length.');
                    end
                else
                    error('DatImageSource:InvalidDatSize', ...
                        ['Continuous .dat files must have datSize=[Y X] or ' ...
                         'datSize=[Y X T].']);
                end
            end
        end

        function validateFileSize(obj)
            %VALIDATEFILESIZE Ensure file bytes match continuous YXT layout.

            fileInfo = dir(obj.FilePath);

            expectedBytes = obj.Ny * obj.Nx * obj.Nt * obj.BytesPerSample;

            if fileInfo.bytes ~= expectedBytes
                error('DatImageSource:FileSizeMismatch', ...
                    ['File size mismatch for "%s". Expected %.0f bytes for ' ...
                     '[Y,X,T]=[%d,%d,%d] with precision "%s", found %.0f bytes.'], ...
                    obj.FilePath, ...
                    expectedBytes, ...
                    obj.Ny, ...
                    obj.Nx, ...
                    obj.Nt, ...
                    obj.Precision, ...
                    fileInfo.bytes);
            end
        end

        function validateFrameIndex(obj, tIdx)
            %VALIDATEFRAMEINDEX Validate one frame index.

            validateattributes(tIdx, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Nt}, ...
                'DatImageSource', 'tIdx');
        end

        function validatePixelIndex(obj, y, x)
            %VALIDATEPIXELINDEX Validate one image pixel coordinate.

            validateattributes(y, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Ny}, ...
                'DatImageSource', 'y');

            validateattributes(x, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Nx}, ...
                'DatImageSource', 'x');
        end

        function idx = validateContiguousIndexRange(obj, idx, maxIdx, dimName) %#ok<INUSL>
            %VALIDATECONTIGUOUSINDEXRANGE Validate contiguous positive indices.

            if isempty(idx)
                error('DatImageSource:EmptyIndexRange', ...
                    '%s range cannot be empty.', dimName);
            end

            idx = double(idx(:).');

            if any(~isfinite(idx)) || any(mod(idx, 1) ~= 0) || ...
                    any(idx < 1) || any(idx > maxIdx)
                error('DatImageSource:InvalidIndexRange', ...
                    '%s range contains invalid indices.', dimName);
            end

            if any(diff(idx) ~= 1)
                error('DatImageSource:NonContiguousIndexRange', ...
                    '%s range must be contiguous.', dimName);
            end
        end

        function [nYCache, nXCache] = estimateCacheSize(obj)
            %ESTIMATECACHESIZE Estimate largest safe XY cache block.

            [availBytes, totalBytes] = queryRAM();

            if isempty(availBytes) || isempty(totalBytes)
                % Conservative fallback when OS RAM query fails.
                budgetBytes = 512 * 1024^2;
            else
                usableBytes = min(double(availBytes), ...
                    double(totalBytes) * obj.MaxUsableRAMFraction);

                usableBytes = usableBytes - double(totalBytes) * obj.RAMoverhead;
                usableBytes = max(0, usableBytes);

                budgetBytes = usableBytes * obj.CacheRAMFraction;
            end

            if ~isempty(obj.MaxCacheBytes)
                budgetBytes = min(budgetBytes, obj.MaxCacheBytes);
            end

            bytesPerPixelTrace = obj.Nt * obj.BytesPerSample;
            maxCachePixels = floor(budgetBytes / bytesPerPixelTrace);

            if maxCachePixels < 1
                error('DatImageSource:InsufficientRAM', ...
                    ['Not enough available RAM to cache one full temporal pixel ' ...
                     'trace for this dataset.']);
            end

            maxCachePixels = min(maxCachePixels, obj.Ny * obj.Nx);

            imageAspect = obj.Nx / obj.Ny;

            nYCache = floor(sqrt(maxCachePixels / imageAspect));
            nXCache = floor(nYCache * imageAspect);

            nYCache = max(1, min(obj.Ny, nYCache));
            nXCache = max(1, min(obj.Nx, nXCache));

            while nYCache * nXCache > maxCachePixels
                if nXCache >= nYCache && nXCache > 1
                    nXCache = nXCache - 1;
                elseif nYCache > 1
                    nYCache = nYCache - 1;
                else
                    break
                end
            end
        end

        function [yRange, xRange] = centerRanges(obj, yCenter, xCenter, nY, nX)
            %CENTERRANGES Return clipped contiguous ranges around a point.

            yStart = round(yCenter - (nY - 1) / 2);
            xStart = round(xCenter - (nX - 1) / 2);

            yStart = max(1, min(yStart, obj.Ny - nY + 1));
            xStart = max(1, min(xStart, obj.Nx - nX + 1));

            yRange = yStart:(yStart + nY - 1);
            xRange = xStart:(xStart + nX - 1);
        end

        function cache = readTemporalBlock(obj, yRange, xRange)
            %READTEMPORALBLOCK Read [Ycache,Xcache,T] from the .dat file.
            %
            %   This reads one contiguous X slab per frame and immediately
            %   crops Y, avoiding an intermediate [fullY,Xcache,T] array.

            yRange = obj.validateContiguousIndexRange(yRange, obj.Ny, 'Y');
            xRange = obj.validateContiguousIndexRange(xRange, obj.Nx, 'X');

            nY = numel(yRange);
            nX = numel(xRange);

            cache = zeros(nY, nX, obj.Nt, obj.Precision);

            fid = fopen(obj.FilePath, 'r');
            if fid < 0
                error('DatImageSource:FileOpenFailed', ...
                    'Could not open file: "%s".', obj.FilePath);
            end
            cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

            x0 = xRange(1);

            for tIdx = 1:obj.Nt
                offset = obj.frameOffsetBytes(tIdx) + ...
                    double(x0 - 1) * obj.Ny * obj.BytesPerSample;

                obj.safeFseek(fid, offset, ...
                    sprintf('Failed to seek to frame %d, x=%d.', tIdx, x0));

                slab = fread(fid, [obj.Ny, nX], ['*' obj.Precision]);

                if numel(slab) ~= obj.Ny * nX
                    error('DatImageSource:UnexpectedEOF', ...
                        ['Unexpected end of file while reading temporal cache ' ...
                         'at frame %d.'], ...
                        tIdx);
                end

                cache(:, :, tIdx) = slab(yRange, :);
            end
        end

        function tf = isFullImageCached(obj)
            %ISFULLIMAGECACHED Return true when cache covers full [Y,X,T].

            tf = ~isempty(obj.CacheData) && ...
                numel(obj.CacheYRange) == obj.Ny && ...
                numel(obj.CacheXRange) == obj.Nx && ...
                obj.CacheYRange(1) == 1 && obj.CacheYRange(end) == obj.Ny && ...
                obj.CacheXRange(1) == 1 && obj.CacheXRange(end) == obj.Nx;
        end

        function tf = isInsideCachedBlock(obj, yRange, xRange)
            %ISINSIDECACHEDBLOCK Return true if a spatial block is cached.

            tf = ~isempty(obj.CacheData) && ...
                yRange(1) >= obj.CacheYRange(1) && ...
                yRange(end) <= obj.CacheYRange(end) && ...
                xRange(1) >= obj.CacheXRange(1) && ...
                xRange(end) <= obj.CacheXRange(end);
        end

        function offset = frameOffsetBytes(obj, tIdx)
            %FRAMEOFFSETBYTES Return byte offset for one frame.

            offset = double(tIdx - 1) * ...
                double(obj.Ny) * ...
                double(obj.Nx) * ...
                double(obj.BytesPerSample);
        end

        function safeFseek(obj, fid, offset, failMessage) %#ok<INUSL>
            %SAFEFSEEK Execute FSEEK and error on failure.

            status = fseek(fid, offset, 'bof');

            if status ~= 0
                error('DatImageSource:FseekFailed', '%s', failMessage);
            end
        end

        function bytes = getByteSize(obj, precision) %#ok<INUSL>
            %GETBYTESIZE Return bytes per element for one precision string.

            switch lower(char(string(precision)))
                case {'single', 'float32'}
                    bytes = 4;
                case {'double', 'float64'}
                    bytes = 8;
                case {'uint8', 'int8', 'char'}
                    bytes = 1;
                case {'uint16', 'int16'}
                    bytes = 2;
                case {'uint32', 'int32'}
                    bytes = 4;
                case {'uint64', 'int64'}
                    bytes = 8;
                otherwise
                    error('DatImageSource:UnsupportedPrecision', ...
                        'Unsupported precision: "%s".', precision);
            end
        end
    end
end
