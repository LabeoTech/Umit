classdef DatImageSource < handle
%DATIMAGESOURCE RAM-safe direct-read backend for .dat image files.
%
%   src = DatImageSource(filePath)
%   src = DatImageSource(filePath, Name, Value)
%
%   This class provides a DataViewer-oriented access layer for .dat files
%   that hold single-precision MATLAB-order image arrays: headered files
%   and legacy sidecar files. Metadata come from
%   loadMetaData once, in the constructor. Full frames and temporal-cache
%   blocks are read through spatialSlabIO, reopening the file for each read
%   with the already resolved Info, so the file is never held open between
%   reads.
%
%   Supported layouts:
%       {'Y','X','T'}       continuous time series
%       {'Y','X'}           single frame (shown as one frame)
%       {'Y','X','T','E'}   event-split time series (.dat header Phase 8b)
%       {'Y','X','E'}       one map per event (T = 1)
%
%   Event-split files: the E axis is mapped onto the SaveFolder's
%   events.mat by resolveDatEventMapping (EventMapping property). getSize
%   returns [Y X T E]; getFrame, getFrameBlock, and getPixelTrace take an
%   optional event index, like UMTImageSource. Internally the trailing axes
%   are read as flat frames f = t + (e-1)*T (spatialSlabIO order), so the
%   cache holds all T*E frames of the cached pixels.
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
%       getSize              - Return [Ny Nx Nt Ne].
%       getFrame             - Read one full [Y,X] frame (optionally of event e).
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
        Ne double = 1
        NFrames double = 0          % flat frames read from disk: Nt * Ne

        HasEventAxis logical = false
        EventMapping struct = struct()  % resolveDatEventMapping output (E files)

        Precision char = 'single'
        BytesPerSample double = 4
        FrameRateHz double = NaN

        % Non-empty when the dataset's stored Rig UUID/ID could not be
        % resolved against UMITRigStore at load time (e.g. deleted Rig,
        % dataset moved from another install). Holds the failure reason;
        % callers use this to surface a non-blocking warning to the user.
        RigAssociationIssue char = ''
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
            %   The constructor validates the layout (Y-X-T, Y-X, Y-X-T-E, or
            %   Y-X-E, single class) before any cache access, and maps the E
            %   axis of event-split files onto events.mat.

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

            dataFolder = fileparts(obj.FilePath);
            if isfile(fullfile(dataFolder, 'AcqInfos.mat')) && ...
                    exist('UMITRigStore', 'class') == 8
                % Legacy .dat folders acquire the current Active Rig once;
                % existing UUID/ID associations, including archived history,
                % are preserved by the backend. This is a best-effort
                % convenience step, not a precondition for reading pixel
                % data, so an unresolvable/stale Rig reference (deleted
                % Rig, dataset moved from another install, etc.) must not
                % block the load.
                try
                    UMITRigStore.ensureDatasetRigAssociation(dataFolder);
                catch ME
                    obj.RigAssociationIssue = ME.message;
                    warning('Umitoolbox:DatImageSource:rigAssociationSkipped', ...
                        'Could not resolve/assign a Rig for "%s": %s', ...
                        dataFolder, ME.message);
                end
            end
            obj.Info = loadMetaData(obj.FilePath);
            obj.validateDatLayout(obj.Info);

            obj.Ny = datAxisSize(obj.Info, 'Y');
            obj.Nx = datAxisSize(obj.Info, 'X');
            % A single-frame Y-X file is shown as one frame (.dat header
            % Phase 6b-1); its frame rate is NaN.
            obj.Nt = max(1, datAxisSize(obj.Info, 'T'));
            obj.HasEventAxis = any(strcmpi(cellstr(string(obj.Info.dimNames)), 'E'));
            if obj.HasEventAxis
                obj.Ne = datAxisSize(obj.Info, 'E');
            end
            obj.NFrames = obj.Nt * obj.Ne;

            obj.Precision = char(string(obj.Info.dataClass));
            obj.BytesPerSample = obj.getByteSize(obj.Precision);

            if ~isempty(obj.Info.frameRateHz)
                obj.FrameRateHz = double(obj.Info.frameRateHz);
            end

            obj.validateFileSize();

            if obj.HasEventAxis
                obj.EventMapping = resolveDatEventMapping(obj.Info, dataFolder);
            end

            % Initial cache: largest safe cache centered in the image.
            obj.updateCacheAround(round(obj.Ny / 2), round(obj.Nx / 2));
        end

        function sz = getSize(obj)
            %GETSIZE Return canonical viewer size [Y X T E].
            %
            %   sz = obj.getSize()
            %
            %   E is 1 for continuous files; T is 1 for Y-X and Y-X-E files.

            sz = [obj.Ny, obj.Nx, obj.Nt, obj.Ne];
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
            %GETEVENTINFO Return the UMT-style eventInfo of the E axis.
            %
            %   For event-split files: the mapping onto events.mat (see
            %   resolveDatEventMapping and EventMapping.status). For continuous
            %   files: an empty struct (events come from EventsManager).

            eventInfo = struct();
            if obj.HasEventAxis && isfield(obj.EventMapping, 'eventInfo')
                eventInfo = obj.EventMapping.eventInfo;
            end
        end

        function mapping = refreshEventMapping(obj)
            %REFRESHEVENTMAPPING Re-map the E axis onto the current events.mat.
            %
            %   Call after events.mat changed (for example after editing
            %   events in Events Manager), so ignore flags and names follow it.
            %   Continuous files return an empty struct.

            mapping = struct();
            if ~obj.HasEventAxis
                return
            end
            obj.EventMapping = resolveDatEventMapping(obj.Info, fileparts(obj.FilePath));
            mapping = obj.EventMapping;
        end

        function f = flatFrameIndex(obj, tIdx, eIdx)
            %FLATFRAMEINDEX Flat frame index of (t, e): t + (e-1)*Nt.
            if nargin < 3 || isempty(eIdx)
                eIdx = 1;
            end
            obj.validateFrameIndex(tIdx);
            obj.validateEventIndex(eIdx);
            f = double(tIdx) + (double(eIdx) - 1) * obj.Nt;
        end

        function frame = getFrame(obj, tIdx, eIdx)
            %GETFRAME Read one full frame from disk.
            %
            %   frame = obj.getFrame(tIdx)
            %   frame = obj.getFrame(tIdx, eIdx)   % event-split files
            %
            %   Output:
            %       frame - Numeric matrix of size [Y,X].

            if nargin < 3
                eIdx = [];
            end
            f = obj.flatFrameIndex(tIdx, eIdx);

            if obj.isFullImageCached()
                frame = obj.CacheData(:, :, f);
                return
            end

            reader = obj.openReader();
            cleanupObj = onCleanup(@() spatialSlabIO('close', reader)); %#ok<NASGU>

            frame = spatialSlabIO('read', reader, 1:obj.Nx, f);
        end

        function block = getFrameBlock(obj, tIdx, yRange, xRange, eIdx)
            %GETFRAMEBLOCK Read one spatial block from one frame.
            %
            %   block = obj.getFrameBlock(tIdx, yRange, xRange)
            %   block = obj.getFrameBlock(tIdx, yRange, xRange, eIdx)
            %
            %   Inputs:
            %       tIdx   - Frame index.
            %       yRange - Contiguous Y indices.
            %       xRange - Contiguous X indices.
            %       eIdx   - Event index (event-split files; default 1).
            %
            %   Output:
            %       block  - Numeric matrix [numel(yRange), numel(xRange)].

            if nargin < 5
                eIdx = [];
            end
            tIdx = obj.flatFrameIndex(tIdx, eIdx);
            yRange = obj.validateContiguousIndexRange(yRange, obj.Ny, 'Y');
            xRange = obj.validateContiguousIndexRange(xRange, obj.Nx, 'X');

            if obj.isInsideCachedBlock(yRange, xRange)
                yLocal = yRange - obj.CacheYRange(1) + 1;
                xLocal = xRange - obj.CacheXRange(1) + 1;
                block = obj.CacheData(yLocal, xLocal, tIdx);
                return
            end

            reader = obj.openReader();
            cleanupObj = onCleanup(@() spatialSlabIO('close', reader)); %#ok<NASGU>

            slab = spatialSlabIO('read', reader, xRange, tIdx);
            block = slab(yRange, :);
        end

        function [trace, status] = getPixelTrace(obj, y, x, eIdx)
            %GETPIXELTRACE Return full temporal trace for one pixel.
            %
            %   trace = obj.getPixelTrace(y, x)
            %   [trace, status] = obj.getPixelTrace(y, x)
            %   [trace, status] = obj.getPixelTrace(y, x, eIdx)
            %
            %   For event-split files the trace covers T of event eIdx
            %   (default 1).
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
            if nargin < 4 || isempty(eIdx)
                eIdx = 1;
            end
            obj.validateEventIndex(eIdx);

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

            frames = (double(eIdx) - 1) * obj.Nt + (1:obj.Nt);
            trace = squeeze(obj.CacheData(yLocal, xLocal, frames));
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
            %       frameIdx - Source frame indices (flat indices for event-split
            %                  files, see flatFrameIndex). The shape controls the output:
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
                obj.NFrames);
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
                    any(finiteFrameIdx < 1) || any(finiteFrameIdx > obj.NFrames)
                error('DatImageSource:InvalidROIFrameIndex', ...
                    'frameIdx contains indices outside the valid range [1,%d].', obj.NFrames);
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

            reader = obj.openReader();
            cleanupObj = onCleanup(@() spatialSlabIO('close', reader)); %#ok<NASGU>

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

                        frameBlock = obj.readFullFrameList(reader, chunkFrames);
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

                    frameBlock = obj.readFullFrameList(reader, chunkFrames);
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

        function frameBlock = readFullFrameList(obj, reader, frameList)
            %READFULLFRAMELIST Read full frames through an open spatialSlabIO reader.
            %
            %   Frames are returned in the order given; repeated frame indices
            %   are read once and copied.

            frameList = double(frameList(:).');

            if isempty(frameList)
                frameBlock = zeros(obj.Ny, obj.Nx, 0, obj.Precision);
                return
            end

            if any(~isfinite(frameList)) || any(mod(frameList, 1) ~= 0)
                error('DatImageSource:InvalidFrameList', ...
                    'Frame list contains invalid frame indices.');
            end

            obj.validateFlatFrameIndex(min(frameList));
            obj.validateFlatFrameIndex(max(frameList));

            [uniqueFrames, ~, whichFrame] = unique(frameList, 'stable');
            frameBlock = spatialSlabIO('read', reader, 1:obj.Nx, uniqueFrames);
            if numel(uniqueFrames) ~= numel(frameList)
                frameBlock = frameBlock(:, :, whichFrame);
            end
        end

        function validateDatLayout(obj, Info) %#ok<INUSL>
            %VALIDATEDATLAYOUT Accept Y-X-T, Y-X, Y-X-T-E, and Y-X-E single files.
            %
            %   Uses the .dat Info schema from loadMetaData. Single-frame Y-X
            %   files are shown as one frame; event-split files (E last) are
            %   supported since .dat header Phase 8b. Other layouts and data
            %   classes other than single fail fast.

            if ~isfield(Info, 'format') || ~isfield(Info, 'filePath')
                error('DatImageSource:InvalidFileType', ...
                    'DatImageSource can only open .dat files.');
            end

            if ~isfield(Info, 'dimNames') || isempty(Info.dimNames)
                error('DatImageSource:MissingDimNames', ...
                    '.dat metadata must contain dimNames.');
            end

            dimNames = upper(cellstr(string(Info.dimNames(:).')));

            supported = {{'Y', 'X', 'T'}, {'Y', 'X'}, {'Y', 'X', 'T', 'E'}, {'Y', 'X', 'E'}};
            if ~any(cellfun(@(c) isequal(dimNames, c), supported))
                error('DatImageSource:UnsupportedDatLayout', ...
                    ['Unsupported .dat layout: {%s}. DatImageSource supports dimNames ' ...
                     '{Y,X,T}, {Y,X}, {Y,X,T,E}, and {Y,X,E}.'], ...
                    strjoin(dimNames, ', '));
            end

            if ~isfield(Info, 'dataClass') || isempty(Info.dataClass)
                error('DatImageSource:MissingDatatype', ...
                    '.dat metadata must contain dataClass.');
            end

            if ~strcmpi(char(string(Info.dataClass)), 'single')
                error('DatImageSource:UnsupportedDatatype', ...
                    'Only single-precision .dat files are currently supported.');
            end

            if ~isfield(Info, 'dimSizes') || numel(Info.dimSizes) ~= numel(dimNames)
                error('DatImageSource:MissingCoreMetadata', ...
                    '.dat metadata must contain the size of every axis.');
            end

            for k = 1:numel(dimNames)
                validateattributes(double(Info.dimSizes(k)), {'numeric'}, ...
                    {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
                    'DatImageSource', [dimNames{k} ' size']);
            end
        end

        function validateFileSize(obj)
            %VALIDATEFILESIZE Ensure the file holds the whole array.
            %
            %   The header of a headered file (Info.dataOffset bytes) is not
            %   counted as data.

            fileInfo = dir(obj.FilePath);

            dataOffset = double(obj.Info.dataOffset);
            expectedBytes = dataOffset + obj.Ny * obj.Nx * obj.NFrames * obj.BytesPerSample;

            if fileInfo.bytes < expectedBytes
                error('DatImageSource:FileSizeMismatch', ...
                    ['File size mismatch for "%s". Expected %.0f bytes for ' ...
                     '[Y,X,frames]=[%d,%d,%d] with precision "%s" after a %d-byte ' ...
                     'header, found %.0f bytes.'], ...
                    obj.FilePath, ...
                    expectedBytes, ...
                    obj.Ny, ...
                    obj.Nx, ...
                    obj.NFrames, ...
                    obj.Precision, ...
                    dataOffset, ...
                    fileInfo.bytes);
            end
        end

        function validateFrameIndex(obj, tIdx)
            %VALIDATEFRAMEINDEX Validate one frame index (T axis).

            validateattributes(tIdx, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Nt}, ...
                'DatImageSource', 'tIdx');
        end

        function validateEventIndex(obj, eIdx)
            %VALIDATEEVENTINDEX Validate one event index (E axis; 1 without E).

            validateattributes(eIdx, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Ne}, ...
                'DatImageSource', 'eIdx');
        end

        function validateFlatFrameIndex(obj, f)
            %VALIDATEFLATFRAMEINDEX Validate one flat frame index (1..Nt*Ne).

            validateattributes(f, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.NFrames}, ...
                'DatImageSource', 'frame');
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

            bytesPerPixelTrace = obj.NFrames * obj.BytesPerSample;
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
            %READTEMPORALBLOCK Read [Ycache,Xcache,T*E] (flat frames) from the .dat file.
            %
            %   Reads the cached X columns in chunks of consecutive frames
            %   through spatialSlabIO and crops Y per chunk, so peak memory
            %   stays close to the cache size instead of holding a
            %   full-height [Ny,Xcache,T] copy.

            yRange = obj.validateContiguousIndexRange(yRange, obj.Ny, 'Y');
            xRange = obj.validateContiguousIndexRange(xRange, obj.Nx, 'X');

            nY = numel(yRange);
            nX = numel(xRange);

            cache = zeros(nY, nX, obj.NFrames, obj.Precision);

            reader = obj.openReader();
            cleanupObj = onCleanup(@() spatialSlabIO('close', reader)); %#ok<NASGU>

            chunkBytes = 64 * 1024^2;
            framesPerChunk = max(1, floor(chunkBytes / (obj.Ny * nX * obj.BytesPerSample)));

            for t1 = 1:framesPerChunk:obj.NFrames
                t2 = min(t1 + framesPerChunk - 1, obj.NFrames);
                slab = spatialSlabIO('read', reader, xRange, t1:t2);
                cache(:, :, t1:t2) = slab(yRange, :, :);
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

        function reader = openReader(obj)
            %OPENREADER Open a spatialSlabIO reader with the resolved Info.
            %
            %   The caller closes it (typically with onCleanup). Reusing
            %   obj.Info avoids calling loadMetaData on every read.

            reader = spatialSlabIO('open', obj.FilePath, 'Info', obj.Info);
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
