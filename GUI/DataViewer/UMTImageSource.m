classdef UMTImageSource < handle
    %UMTIMAGESOURCE In-memory image backend for .umt image structures.
    %
    %   src = UMTImageSource(filePath)
    %   src = UMTImageSource(filePath, entryName)
    %
    %   This class provides a DataViewer-oriented access layer for .umt files
    %   with kind='image'. It loads one image entry into memory and exposes a
    %   common backend interface compatible with DatImageSource.
    %
    %   Supported image layouts are defined by getUMTSchema:
    %       {'Y','X'}
    %       {'Y','X','T'}
    %       {'Y','X','E'}
    %       {'Y','X','T','E'}
    %
    %   Event metadata, when present, is read from top-level umt.eventInfo.
    %
    %   Notes:
    %       - loadMetaData must be available on the MATLAB path.
    %       - validateUMTStruct must be available on the MATLAB path.
    %       - This class intentionally treats UMT image entries as fully cached.
    %         A lazy MAT-file backend can be added later if needed.

    properties (SetAccess = private)
        FilePath char = ''
        Info struct = struct()

        UMT struct = struct()
        EntryName char = ''
        EntryNames cell = {}

        Data
        DimNames cell = {}

        Labels struct = struct()
        EventInfo struct = struct()

        Ny double = 0
        Nx double = 0
        Nt double = 1
        Ne double = 1
        FrameRateHz double = NaN
    end

    methods
        function obj = UMTImageSource(filePath, entryName)
            %UMTIMAGESOURCE Construct a .umt image source.
            %
            %   obj = UMTImageSource(filePath)
            %   obj = UMTImageSource(filePath, entryName)
            %
            %   If entryName is omitted and the UMT contains exactly one
            %   image entry, that entry is selected automatically. If the
            %   file has multiple entries, an error is raised so the GUI can
            %   request explicit selection from the user.

            if nargin < 2
                entryName = '';
            end

            validateattributes(filePath, {'char', 'string'}, {'scalartext'}, ...
                'UMTImageSource', 'filePath');

            filePath = char(string(filePath));
            if isempty(fileparts(filePath))
                filePath = fullfile(pwd, filePath);
            end

            if ~isfile(filePath)
                error('UMTImageSource:FileNotFound', ...
                    'File not found: "%s".', filePath);
            end

            [~, ~, ext] = fileparts(filePath);
            if ~strcmpi(ext, '.umt')
                error('UMTImageSource:InvalidExtension', ...
                    'UMTImageSource can only open .umt files.');
            end

            obj.FilePath = filePath;
            obj.Info = loadMetaData(obj.FilePath);
            obj.UMT = obj.loadUMTStruct(obj.FilePath);

            validateUMTStruct(obj.UMT, 'requireEventInfo', true);

            if ~strcmpi(obj.UMT.kind, 'image')
                error('UMTImageSource:UnsupportedKind', ...
                    'Only UMT files with kind="image" can be loaded by UMTImageSource.');
            end

            obj.EntryNames = fieldnames(obj.UMT.data);
            obj.EntryName = obj.resolveEntryName(entryName);
            obj.loadSelectedEntry();

            if isfield(obj.UMT, 'labels')
                obj.Labels = obj.UMT.labels;
            end

            if isfield(obj.UMT, 'eventInfo')
                obj.EventInfo = obj.UMT.eventInfo;
            end

            if isfield(obj.Info, 'FrameRateHz') && ~isempty(obj.Info.FrameRateHz)
                obj.FrameRateHz = double(obj.Info.FrameRateHz);
            elseif isfield(obj.Info, 'Freq') && ~isempty(obj.Info.Freq)
                obj.FrameRateHz = double(obj.Info.Freq);
            end
        end

        function sz = getSize(obj)
            %GETSIZE Return canonical viewer size [Y X T E].

            sz = [obj.Ny, obj.Nx, obj.Nt, obj.Ne];
        end

        function info = getInfo(obj)
            %GETINFO Return flat metadata structure.

            info = obj.Info;
        end

        function labels = getLabels(obj)
            %GETLABELS Return top-level UMT labels.

            labels = obj.Labels;
        end

        function eventInfo = getEventInfo(obj)
            %GETEVENTINFO Return top-level UMT eventInfo.

            eventInfo = obj.EventInfo;
        end

        function entryName = getEntryName(obj)
            %GETENTRYNAME Return the currently loaded UMT entry name.

            entryName = obj.EntryName;

        end

        function entryNames = getEntryNames(obj)
            %GETENTRYNAMES Return all available UMT entry names.

            entryNames = obj.EntryNames;

        end

        function summaryTable = getEntrySummary(obj)
            %GETENTRYSUMMARY Return summary table for all entries in the loaded UMT.

            summaryTable = UMTImageSource.buildEntrySummaryTable(obj.UMT);

        end


        function frame = getFrame(obj, tIdx, eIdx)
            %GETFRAME Return one image frame.
            %
            %   frame = obj.getFrame(tIdx)
            %   frame = obj.getFrame(tIdx, eIdx)
            %
            %   For entries with no T dimension, tIdx is ignored but must be
            %   1 if supplied. For entries with no E dimension, eIdx is
            %   ignored but must be 1 if supplied.

            if nargin < 2 || isempty(tIdx)
                tIdx = 1;
            end
            if nargin < 3 || isempty(eIdx)
                eIdx = 1;
            end

            obj.validateFrameIndex(tIdx);
            obj.validateEventIndex(eIdx);

            dimMap = obj.getDimMap();

            idx = repmat({':'}, 1, ndims(obj.Data));

            if isfield(dimMap, 'T')
                idx{dimMap.T} = tIdx;
            end

            if isfield(dimMap, 'E')
                idx{dimMap.E} = eIdx;
            end

            frame = obj.Data(idx{:});
            frame = squeeze(frame);

            if ~isequal(size(frame), [obj.Ny, obj.Nx])
                frame = reshape(frame, obj.Ny, obj.Nx);
            end
        end

        function block = getFrameBlock(obj, tIdx, yRange, xRange, eIdx)
            %GETFRAMEBLOCK Return one spatial block from one frame.
            %
            %   block = obj.getFrameBlock(tIdx, yRange, xRange)
            %   block = obj.getFrameBlock(tIdx, yRange, xRange, eIdx)

            if nargin < 5 || isempty(eIdx)
                eIdx = 1;
            end

            yRange = obj.validateIndexRange(yRange, obj.Ny, 'Y');
            xRange = obj.validateIndexRange(xRange, obj.Nx, 'X');

            frame = obj.getFrame(tIdx, eIdx);
            block = frame(yRange, xRange);
        end

        function [trace, status] = getPixelTrace(obj, y, x, eIdx)
            %GETPIXELTRACE Return temporal trace for one pixel.
            %
            %   trace = obj.getPixelTrace(y, x)
            %   trace = obj.getPixelTrace(y, x, eIdx)
            %   [trace, status] = obj.getPixelTrace(...)
            %
            %   UMT entries are fully in memory, so status is always 'ok'.

            if nargin < 4 || isempty(eIdx)
                eIdx = 1;
            end

            obj.validatePixelIndex(y, x);
            obj.validateEventIndex(eIdx);

            dimMap = obj.getDimMap();
            idx = repmat({':'}, 1, ndims(obj.Data));

            idx{dimMap.Y} = y;
            idx{dimMap.X} = x;

            if isfield(dimMap, 'E')
                idx{dimMap.E} = eIdx;
            end

            if ~isfield(dimMap, 'T')
                trace = obj.Data(idx{:});
                trace = trace(:);
            else
                trace = squeeze(obj.Data(idx{:}));
                trace = trace(:);
            end

            status = 'ok';
        end

        function tf = isInsideCache(obj, varargin) %#ok<INUSD>
            %ISINSIDECACHE Return true because UMT data are fully cached.

            tf = true;
        end

        function tf = hasPartialTemporalCache(obj)
            %HASPARTIALTEMPORALCACHE Return false for in-memory UMT source.

            tf = false;
        end

        function updateCacheAround(obj, varargin) %#ok<INUSD>
            %UPDATECACHEAROUND No-op because UMT data are fully cached.

        end

        function rect = getCacheRectangle(obj)
            %GETCACHERECTANGLE Return [] because there is no partial cache.

            rect = [];
        end

        function setCacheMode(obj, varargin) %#ok<INUSD>
            %SETCACHEMODE No-op for UMT source.

        end

        function setCacheLocked(obj, varargin) %#ok<INUSD>
            %SETCACHELOCKED No-op for UMT source.

        end

        function tf = isCacheLocked(obj)
            %ISCACHELOCKED Return false for UMT source.

            tf = false;
        end

        function txt = getCacheStatusText(obj)
            %GETCACHESTATUSTEXT Return compact cache status text.

            txt = 'Cache: full UMT entry';
        end
    end

    methods (Static)
        function summaryTable = inspectEntries(filePath)
            %INSPECTENTRIES Return a summary table of image entries in a .umt file.
            %
            %   summaryTable = UMTImageSource.inspectEntries(filePath)
            %
            %   The returned table contains:
            %       Name
            %       DimNames
            %       Size
            %       Class
            %       HasT
            %       HasE

            validateattributes(filePath, {'char', 'string'}, {'scalartext'}, ...
                'UMTImageSource.inspectEntries', 'filePath');

            filePath = char(string(filePath));

            if isempty(fileparts(filePath))
                filePath = fullfile(pwd, filePath);
            end

            if ~isfile(filePath)
                error('UMTImageSource:FileNotFound', ...
                    'File not found: "%s".', filePath);
            end

            [~, ~, ext] = fileparts(filePath);
            if ~strcmpi(ext, '.umt')
                error('UMTImageSource:InvalidExtension', ...
                    'UMTImageSource.inspectEntries can only inspect .umt files.');
            end

            umt = UMTImageSource.loadUMTStruct(filePath);

            validateUMTStruct(umt, 'requireEventInfo', true);

            if ~strcmpi(umt.kind, 'image')
                error('UMTImageSource:UnsupportedKind', ...
                    'Only UMT files with kind="image" can be loaded by UMTImageSource.');
            end

            summaryTable = UMTImageSource.buildEntrySummaryTable(umt);
        end
    end

    methods (Access = private)
        function entryName = resolveEntryName(obj, entryName)
            %RESOLVEENTRYNAME Select an entry or raise ambiguity error.

            if isempty(entryName)
                if numel(obj.EntryNames) == 1
                    entryName = obj.EntryNames{1};
                    return
                end

                error('UMTImageSource:AmbiguousEntry', ...
                    ['The UMT file contains multiple image entries: %s. ' ...
                    'Please select which entry to load.'], ...
                    strjoin(obj.EntryNames(:).', ', '));
            end

            entryName = char(string(entryName));

            if ~isfield(obj.UMT.data, entryName)
                error('UMTImageSource:InvalidEntryName', ...
                    'Entry "%s" was not found in the UMT file.', entryName);
            end
        end

        function loadSelectedEntry(obj)
            %LOADSELECTEDENTRY Load selected UMT data entry into properties.

            entry = obj.UMT.data.(obj.EntryName);

            obj.Data = entry.value;
            obj.DimNames = cellstr(string(entry.dimNames(:).'));

            dimMap = obj.getDimMap();
            dataSize = obj.getDeclaredSize();

            obj.Ny = dataSize(dimMap.Y);
            obj.Nx = dataSize(dimMap.X);

            if isfield(dimMap, 'T')
                obj.Nt = dataSize(dimMap.T);
            else
                obj.Nt = 1;
            end

            if isfield(dimMap, 'E')
                obj.Ne = dataSize(dimMap.E);
            else
                obj.Ne = 1;
            end
        end

        function dimMap = getDimMap(obj)
            %GETDIMMAP Return dimension-name to index mapping.

            dimMap = struct();

            for iDim = 1:numel(obj.DimNames)
                thisDim = obj.DimNames{iDim};
                dimMap.(thisDim) = iDim;
            end

            if ~isfield(dimMap, 'Y') || ~isfield(dimMap, 'X')
                error('UMTImageSource:InvalidImageEntry', ...
                    'Image entries must contain Y and X dimensions.');
            end
        end

        function sz = getDeclaredSize(obj)
            %GETDECLAREDSIZE Return array size compatible with DimNames.

            nDimsExpected = numel(obj.DimNames);
            sz = size(obj.Data);

            if numel(sz) < nDimsExpected
                sz(end+1:nDimsExpected) = 1;
            elseif numel(sz) > nDimsExpected
                sz = sz(1:nDimsExpected);
            end
        end

        function validateFrameIndex(obj, tIdx)
            %VALIDATEFRAMEINDEX Validate one T index.

            validateattributes(tIdx, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Nt}, ...
                'UMTImageSource', 'tIdx');
        end

        function validateEventIndex(obj, eIdx)
            %VALIDATEEVENTINDEX Validate one E index.

            validateattributes(eIdx, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Ne}, ...
                'UMTImageSource', 'eIdx');
        end

        function validatePixelIndex(obj, y, x)
            %VALIDATEPIXELINDEX Validate one image pixel coordinate.

            validateattributes(y, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Ny}, ...
                'UMTImageSource', 'y');

            validateattributes(x, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', '>=', 1, '<=', obj.Nx}, ...
                'UMTImageSource', 'x');
        end

        function idx = validateIndexRange(obj, idx, maxIdx, dimName) %#ok<INUSL>
            %VALIDATEINDEXRANGE Validate positive integer vector.

            if isempty(idx)
                error('UMTImageSource:EmptyIndexRange', ...
                    '%s range cannot be empty.', dimName);
            end

            idx = double(idx(:).');

            if any(~isfinite(idx)) || any(mod(idx, 1) ~= 0) || ...
                    any(idx < 1) || any(idx > maxIdx)
                error('UMTImageSource:InvalidIndexRange', ...
                    '%s range contains invalid indices.', dimName);
            end
        end
    end

    methods (Access = private, Static)
        function umt = loadUMTStruct(fileName)
            %LOADUMTSTRUCT Load the first scalar UMT struct from a .umt file.
            %
            %   This intentionally mirrors the robust lookup behavior used
            %   by loadMetaData rather than assuming a fixed variable name.

            S = load(fileName, '-mat');

            if isstruct(S) && isscalar(S) && ...
                    all(ismember({'version', 'kind', 'data'}, fieldnames(S)))
                umt = S;
                return
            end

            fields = fieldnames(S);
            for iField = 1:numel(fields)
                candidate = S.(fields{iField});

                if isstruct(candidate) && isscalar(candidate) && ...
                        all(ismember({'version', 'kind', 'data'}, fieldnames(candidate)))
                    umt = candidate;
                    return
                end
            end

            error('UMTImageSource:InvalidUMTFile', ...
                'No scalar UMT struct was found in "%s".', fileName);
        end

        function summaryTable = buildEntrySummaryTable(umt)
            %BUILDENTRYSUMMARYTABLE Build display summary for UMT entries.

            entryNames = fieldnames(umt.data);
            nEntries = numel(entryNames);

            Name = strings(nEntries, 1);
            DimNames = strings(nEntries, 1);
            SizeText = strings(nEntries, 1);
            DataClass = strings(nEntries, 1);
            HasT = strings(nEntries, 1);
            HasE = strings(nEntries, 1);

            for iEntry = 1:nEntries
                entryName = entryNames{iEntry};
                entry = umt.data.(entryName);

                dimNames = cellstr(string(entry.dimNames(:).'));
                value = entry.value;

                nDimsExpected = numel(dimNames);
                dataSize = size(value);

                if numel(dataSize) < nDimsExpected
                    dataSize(end+1:nDimsExpected) = 1;
                elseif numel(dataSize) > nDimsExpected
                    dataSize = dataSize(1:nDimsExpected);
                end

                bHasT = any(strcmp(dimNames, 'T'));
                bHasE = any(strcmp(dimNames, 'E'));

                Name(iEntry) = string(entryName);
                DimNames(iEntry) = strjoin(string(dimNames), ' x ');
                SizeText(iEntry) = "[" + strjoin(string(dataSize), " x ") + "]";
                DataClass(iEntry) = string(class(value));

                if bHasT
                    HasT(iEntry) = "Yes";
                else
                    HasT(iEntry) = "No";
                end

                if bHasE
                    HasE(iEntry) = "Yes";
                else
                    HasE(iEntry) = "No";
                end
            end

            summaryTable = table( ...
                Name, ...
                DimNames, ...
                SizeText, ...
                DataClass, ...
                HasT, ...
                HasE, ...
                'VariableNames', {'Name', 'DimNames', 'Size', 'Class', 'HasT', 'HasE'});

        end
    end
end
