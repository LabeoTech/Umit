classdef TestIOHelpers
    %TESTIOHELPERS Helper utilities for UMT I/O unit tests.
    %
    %   This helper class provides:
    %       - temporary-folder management
    %       - loading or synthesizing base 3D image data
    %       - generation of valid test payloads for all allowed .umt
    %         dimension combinations for kind='image' and kind='roi'
    %       - generation of minimal AcqInfoStream/ImportedChannels fixtures

    methods (Static)
        function folderPath = makeTempFolder(prefix)
            %MAKETEMPFOLDER Create a unique temporary folder.

            if nargin < 1 || isempty(prefix)
                prefix = 'tUMTIO_';
            end

            folderPath = fullfile(tempdir, [prefix char(java.util.UUID.randomUUID)]);
            mkdir(folderPath);
        end

        function removeFolderIfExists(folderPath)
            %REMOVEFOLDERIFEXISTS Remove folder recursively if it exists.

            if ~isempty(folderPath) && isfolder(folderPath)
                rmdir(folderPath, 's');
            end
        end

        function [data, sourceName] = getTest3DArray(sourceFolder)
            %GETTEST3DARRAY Get a 3D single array from MAT files or fallback.
            %
            %   Searches sourceFolder recursively for MAT files containing a
            %   numeric 3D array. If none is found, returns a synthetic 3D
            %   single array.

            data = [];
            sourceName = 'syntheticFallback';

            if nargin >= 1 && strlength(string(sourceFolder)) > 0 && isfolder(sourceFolder)
                [data, sourceName] = TestIOHelpers.tryLoadReal3DArray(char(string(sourceFolder)));
            end

            if isempty(data)
                data = TestIOHelpers.makeSynthetic3DArray();
                sourceName = 'syntheticFallback';
            end

            data = single(data);
        end

        function [data, sourceName] = tryLoadReal3DArray(sourceFolder)
            %TRYLOADREAL3DARRAY Search recursively for a MAT file with a 3D array.

            data = [];
            sourceName = '';

            files = dir(fullfile(sourceFolder, '**', '*.mat'));
            for k = 1:numel(files)
                filePath = fullfile(files(k).folder, files(k).name);

                try
                    S = load(filePath);
                catch
                    continue
                end

                vars = fieldnames(S);
                for iVar = 1:numel(vars)
                    thisVar = S.(vars{iVar});

                    if isnumeric(thisVar) && ndims(thisVar) == 3 && ~isempty(thisVar)
                        data = single(thisVar);
                        sourceName = filePath;
                        return
                    end
                end
            end
        end

        function data = makeSynthetic3DArray()
            %MAKESYNTHETIC3DARRAY Create synthetic single 3D test data.

            Ny = 16;
            Nx = 12;
            Nt = 9;

            data = reshape(single(1:(Ny * Nx * Nt)), [Ny, Nx, Nt]);
        end

        function AcqInfoStream = makeAcqInfoStreamForData(data, frameRateHz)
            %MAKEACQINFOSTREAMFORDATA Create minimal AcqInfoStream for .dat tests.

            if nargin < 2 || isempty(frameRateHz)
                frameRateHz = 10;
            end

            [Ny, Nx, Nt] = size(data);

            AcqInfoStream = struct();
            AcqInfoStream.Height = Ny;
            AcqInfoStream.Width = Nx;
            AcqInfoStream.Length = Nt;
            AcqInfoStream.FrameRateHz = frameRateHz;
            AcqInfoStream.Datatype = 'single';
            AcqInfoStream.dim_names = {'Y','X','T'};
            AcqInfoStream.ExposureMsec = 5;
            AcqInfoStream.MultiCam = false;
        end

        function channelInfo = makeImportedChannelInfo(datFile, lengthFrames, frameRateHz, exposureMsec, camIdx)
            %MAKEIMPORTEDCHANNELINFO Create one ImportedChannels test entry.

            if nargin < 4 || isempty(exposureMsec)
                exposureMsec = 5;
            end

            channelInfo = struct();
            channelInfo.DatFile = char(string(datFile));
            channelInfo.Tag = erase(lower(channelInfo.DatFile), '.dat');
            channelInfo.Color = channelInfo.Tag;
            channelInfo.Length = double(lengthFrames);
            channelInfo.FrameRateHz = double(frameRateHz);
            channelInfo.ExposureMsec = double(exposureMsec);

            if nargin >= 5 && ~isempty(camIdx)
                channelInfo.CamIdx = double(camIdx);
            end
        end

        function [value, labels] = makeUMTPayload(kind, dimNames, baseData)
            %MAKEUMTPAYLOAD Create a valid payload/labels pair for one schema layout.
            %
            %   [value, labels] = makeUMTPayload(kind, dimNames, baseData)
            %
            %   Inputs:
            %       kind     - 'image' or 'roi'
            %       dimNames - allowed .umt dimension-name cell array
            %       baseData - 3D single image array used as numeric seed
            %
            %   Outputs:
            %       value    - numeric/logical payload with the requested layout
            %       labels   - optional shared top-level labels struct

            kind = lower(char(string(kind)));
            dimNames = cellstr(string(dimNames));
            labels = struct();

            switch kind
                case 'image'
                    [value, labels] = TestIOHelpers.makeImagePayload(dimNames, baseData);

                case 'roi'
                    [value, labels] = TestIOHelpers.makeROIPayload(dimNames, baseData);

                otherwise
                    error('TestIOHelpers:makeUMTPayload:invalidKind', ...
                        'Unsupported kind "%s".', kind);
            end
        end

        function [value, labels] = makeImagePayload(dimNames, baseData)
            %MAKEIMAGEPAYLOAD Create valid image-kind payloads.

            [Ny, Nx, ~] = size(baseData);
            dimKey = strjoin(cellstr(string(dimNames)), '|');
            labels = struct();

            nROI = 4;
            nF = 2;
            baseMap = mean(baseData, 3);

            switch dimKey
                case ''
                    value = mean(baseData, 'all');

                case 'Y|X'
                    value = baseMap;

                case 'Y|X|T'
                    value = baseData;

                case 'Y|X|ROI'
                    value = zeros(Ny, Nx, nROI, 'single');
                    for iROI = 1:nROI
                        value(:,:,iROI) = baseMap + iROI;
                    end
                    labels.ROI = arrayfun(@(k) sprintf('ROI_%02d', k), 1:nROI, 'UniformOutput', false);

                case 'Y|X|E'
                    value = cat(3, ...
                        baseMap, ...
                        baseMap + 1, ...
                        baseMap + 2);
                    labels.E = {'Event_1', 'Event_2', 'Event_3'};

                case 'Y|X|F'
                    value = cat(3, ...
                        baseMap, ...
                        baseMap + max(baseMap(:)) + 1);
                    labels.F = {'Feature_1', 'Feature_2'};

                case 'Y|X|T|E'
                    value = cat(4, ...
                        baseData, ...
                        baseData + max(baseData(:)) + 1);
                    labels.E = {'Event_1', 'Event_2'};

                otherwise
                    error('TestIOHelpers:makeImagePayload:invalidDimNames', ...
                        'Unsupported image dimNames pattern: %s.', dimKey);
            end

            % Ensure output uses expected image XY size when applicable.
            if ~isempty(dimNames) && numel(dimNames) >= 2
                assert(size(value, 1) == Ny && size(value, 2) == Nx);
            end
        end

        function [value, labels] = makeROIPayload(dimNames, baseData)
            %MAKEROIPAYLOAD Create valid roi-kind payloads.

            Nt = size(baseData, 3);
            Nroi = 4;
            Ne = 3;
            Nmeasure = 2;
            Npixel = 5;

            roiLabels = arrayfun(@(k) sprintf('ROI_%02d', k), 1:Nroi, 'UniformOutput', false);
            eventLabels = {'Event_1', 'Event_2', 'Event_3'};
            measureLabels = {'Mean', 'Peak'};
            pixelLabels = arrayfun(@(k) sprintf('Pixel_%02d', k), 1:Npixel, 'UniformOutput', false);

            dimKey = strjoin(cellstr(string(dimNames)), '|');
            labels = struct();

            switch dimKey
                case ''
                    value = mean(baseData, 'all');

                case 'ROI'
                    value = reshape(single(1:Nroi), [Nroi, 1]);
                    labels.ROI = roiLabels;

                case 'ROI|T'
                    value = reshape(single(1:(Nroi * Nt)), [Nroi, Nt]);
                    labels.ROI = roiLabels;

                case 'ROI|E'
                    value = reshape(single(1:(Nroi * Ne)), [Nroi, Ne]);
                    labels.ROI = roiLabels;
                    labels.E = eventLabels;

                case 'ROI|F'
                    value = reshape(single(1:(Nroi * Nmeasure)), [Nroi, Nmeasure]);
                    labels.ROI = roiLabels;
                    labels.F = {'F1', 'F2'};

                case 'ROI|Measure'
                    value = reshape(single(1:(Nroi * Nmeasure)), [Nroi, Nmeasure]);
                    labels.ROI = roiLabels;
                    labels.Measure = measureLabels;

                case 'ROI|Measure|E'
                    value = reshape(single(1:(Nroi * Nmeasure * Ne)), [Nroi, Nmeasure, Ne]);
                    labels.ROI = roiLabels;
                    labels.Measure = measureLabels;
                    labels.E = eventLabels;

                case 'ROI|T|E'
                    value = reshape(single(1:(Nroi * Nt * Ne)), [Nroi, Nt, Ne]);
                    labels.ROI = roiLabels;
                    labels.E = eventLabels;

                case 'ROI|Pixel'
                    value = reshape(single(1:(Nroi * Npixel)), [Nroi, Npixel]);
                    labels.ROI = roiLabels;
                    labels.Pixel = pixelLabels;

                case 'ROI|Pixel|T'
                    value = reshape(single(1:(Nroi * Npixel * Nt)), [Nroi, Npixel, Nt]);
                    labels.ROI = roiLabels;
                    labels.Pixel = pixelLabels;

                case 'ROI|Pixel|E'
                    value = reshape(single(1:(Nroi * Npixel * Ne)), [Nroi, Npixel, Ne]);
                    labels.ROI = roiLabels;
                    labels.Pixel = pixelLabels;
                    labels.E = eventLabels;

                case 'ROI|Pixel|T|E'
                    value = reshape(single(1:(Nroi * Npixel * Nt * Ne)), [Nroi, Npixel, Nt, Ne]);
                    labels.ROI = roiLabels;
                    labels.Pixel = pixelLabels;
                    labels.E = eventLabels;

                case 'ROI|ROI'
                    value = reshape(single(1:(Nroi * Nroi)), [Nroi, Nroi]);
                    labels.ROI = roiLabels;

                otherwise
                    error('TestIOHelpers:makeROIPayload:invalidDimNames', ...
                        'Unsupported roi dimNames pattern: %s.', dimKey);
            end
        end

        function umt = appendDefaultEventInfoIfNeeded(umt)
            %APPENDDEFAULTEVENTINFOIFNEEDED Finalize E-based UMT structures.

            if ~TestIOHelpers.usesEventDimension(umt)
                return
            end

            eLen = TestIOHelpers.getSharedEventLength(umt);
            eventID = 1:eLen;
            repetitionIndex = ones(1, eLen);
            eventName = arrayfun(@(k) sprintf('Event_%d', k), 1:eLen, 'UniformOutput', false);

            umt = appendUMTEventInfo(umt, ...
                'eventID', eventID, ...
                'repetitionIndex', repetitionIndex, ...
                'eventName', eventName, ...
                'eventAxisMode', 'instances', ...
                'overwrite', true);
        end

        function tf = usesEventDimension(umt)
            %USESEVENTDIMENSION Return true when any entry uses E.

            tf = false;
            if ~isstruct(umt) || ~isfield(umt, 'data') || ~isstruct(umt.data)
                return
            end

            entryNames = fieldnames(umt.data);
            for iEntry = 1:numel(entryNames)
                dimNames = cellstr(string(umt.data.(entryNames{iEntry}).dimNames));
                if any(strcmp(dimNames, 'E'))
                    tf = true;
                    return
                end
            end
        end

        function eLen = getSharedEventLength(umt)
            %GETSHAREDEVENTLENGTH Return shared E length from UMT entries.

            entryNames = fieldnames(umt.data);
            eLengths = [];

            for iEntry = 1:numel(entryNames)
                entry = umt.data.(entryNames{iEntry});
                dimNames = cellstr(string(entry.dimNames));
                idxE = find(strcmp(dimNames, 'E'), 1, 'first');
                if isempty(idxE)
                    continue
                end

                sz = size(entry.value);
                if numel(sz) < idxE
                    eLengths(end+1) = 1; %#ok<AGROW>
                else
                    eLengths(end+1) = sz(idxE); %#ok<AGROW>
                end
            end

            assert(~isempty(eLengths), 'No E dimension was found.');
            assert(isscalar(unique(eLengths)), 'E dimension length is not shared across entries.');
            eLen = eLengths(1);
        end
    end
end
