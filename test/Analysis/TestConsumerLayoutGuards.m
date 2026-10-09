classdef TestConsumerLayoutGuards < matlab.unittest.TestCase
    %TESTCONSUMERLAYOUTGUARDS Analysis consumers refuse unsupported .dat layouts (Phase 6b-1).
    %
    %   Since writers emit any supported layout, every consumer that assumes
    %   Y-X-T checks the axes of its .dat input and raises
    %   Umitoolbox:<fn>:unsupportedLayout before reading data or writing
    %   outputs. Each consumer is given a Y-X-T-E file and a Y-X-E file; the
    %   SaveFolder must be left unchanged. run_ConvertToTiff accepts Y-X-T-E
    %   (one TIFF per event, see TestRunConvertToTiff) and exports a
    %   single-frame Y-X file as one TIFF page.

    properties
        Folder
    end

    properties (TestParameter)
        consumer = struct( ...
            'genCorrelationMatrix', {{'genCorrelationMatrix', @(f, sv) genCorrelationMatrix(f, sv, 'ROImasks_filename', 'guard.roi')}}, ...
            'split_data_by_event', {{'split_data_by_event', @(f, sv) split_data_by_event(f, sv)}}, ...
            'genRetinotopyMaps', {{'genRetinotopyMaps', @(f, sv) genRetinotopyMaps(f, sv)}})
        % Consumers that take Y-X-T and Y-X-T-E: only Y-X-E is refused.
        yxteConsumer = struct( ...
            'GSR', {{'GSR', @(f, sv) GSR(f, sv)}}, ...
            'normalizeZScore', {{'normalizeZScore', @(f, sv) normalizeZScore(f, sv)}}, ...
            'apply_detrend', {{'apply_detrend', @(f, sv) apply_detrend(f, sv)}}, ...
            'normalizeBSLN', {{'normalizeBSLN', @(f, sv) normalizeBSLN(f, sv)}}, ...
            'run_ConvertToTiff', {{'run_ConvertToTiff', @(f, sv) run_ConvertToTiff(f, sv)}})
        layout = struct( ...
            'YXTE', struct('names', {{'Y', 'X', 'T', 'E'}}, 'size', [6 5 4 2]), ...
            'YXE', struct('names', {{'Y', 'X', 'E'}}, 'size', [6 5 3]))
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            AcqInfoStream = struct('Height', 6, 'Width', 5, 'Length', 4, 'FrameRateHz', 10, ...
                'ExposureMsec', 1, 'Datatype', 'single', 'MultiCam', false);
            save(fullfile(testCase.Folder, 'AcqInfos.mat'), 'AcqInfoStream');
            % split_data_by_event checks that events.mat exists before reading
            % its input; the layout check runs before the file is loaded.
            placeholder = 1;
            save(fullfile(testCase.Folder, 'events.mat'), 'placeholder');
            roi = iMakeRectROI('ROI1', [6 5], 1:3, 1:3, [1 0 0]);
            saveROIFile(fullfile(testCase.Folder, 'guard.roi'), ...
                createROIFile([6 5], 'ROIs', roi));
        end
    end

    methods (Test)
        function unsupportedLayoutIsRefused(testCase, consumer, layout)
            name = consumer{1};
            call = consumer{2};
            f = fullfile(testCase.Folder, 'input.dat');
            data = reshape(single(1:prod(layout.size)), layout.size);
            writeTestDat(f, data, 10, 'DimNames', layout.names);
            before = iListing(testCase.Folder);

            testCase.verifyError(@() call(f, testCase.Folder), ...
                sprintf('Umitoolbox:%s:unsupportedLayout', name));
            testCase.verifyEqual(iListing(testCase.Folder), before, ...
                'no output may be written for a refused input');
        end

        function yxeIsRefusedByYxteConsumers(testCase, yxteConsumer)
            name = yxteConsumer{1};
            call = yxteConsumer{2};
            f = fullfile(testCase.Folder, 'input.dat');
            writeTestDat(f, reshape(single(1:90), [6 5 3]), 10, 'DimNames', {'Y', 'X', 'E'});
            before = iListing(testCase.Folder);

            testCase.verifyError(@() call(f, testCase.Folder), ...
                sprintf('Umitoolbox:%s:unsupportedLayout', name));
            testCase.verifyEqual(iListing(testCase.Folder), before, ...
                'no output may be written for a refused input');
        end

        function convertToTiffExportsASingleFrame(testCase)
            f = fullfile(testCase.Folder, 'frame.dat');
            img = reshape(single(1:30), 6, 5);
            writeTestDat(f, img, 10, 'DimNames', {'Y', 'X'});

            outFile = run_ConvertToTiff(f, testCase.Folder);

            testCase.verifyEqual(outFile, {'img_frame.tif'});
            tifPath = fullfile(testCase.Folder, 'img_frame.tif');
            testCase.assertTrue(isfile(tifPath));
            testCase.verifyNumElements(imfinfo(tifPath), 1, 'one TIFF page');
            testCase.verifyEqual(sort(single(reshape(imread(tifPath), [], 1))), ...
                sort(img(:)), 'the page holds the frame values');
        end
    end
end

function names = iListing(folder)
l = dir(folder);
names = sort({l(~[l.isdir]).name});
end

function roi = iMakeRectROI(name, imageSizeYX, rowRange, colRange, color)
%IMAKERECTROI Deterministic schema-valid rectangle ROI (as in TestGetDataFromROI).
mask = false(imageSizeYX(1), imageSizeYX(2));
mask(rowRange, colRange) = true;
x0 = min(colRange) - 0.5; x1 = max(colRange) + 0.5;
y0 = min(rowRange) - 0.5; y1 = max(rowRange) + 0.5;
vertices = [x0 y0; x1 y0; x1 y1; x0 y1];
pgon = polyshape(vertices(:,1), vertices(:,2));
stats = struct('computedOn', datetime('now'), 'NPixels', nnz(mask), 'areaPx2', nnz(mask), ...
    'areaMM2', [], 'centroidXY_px', [], 'centroidXY_mm', [], 'distanceFromOrigin_px', [], ...
    'distanceFromOrigin_mm', [], 'spatialMean', [], 'spatialStd', [], 'spatialMedian', [], ...
    'spatialMin', [], 'spatialMax', []);
roi = struct('name', name, 'type', 'polygon', 'DOC', datetime('now'), ...
    'modifiedOn', datetime('now'), 'color', color, 'notes', '', ...
    'geometry', struct('polyshape', pgon, 'verticesXY_px', vertices, 'ROIType', 'polygon', ...
        'ROIParameters', struct('ROIType', 'polygon', 'Position', vertices, 'Vertices', vertices)), ...
    'mask', logical(mask), 'stats', stats);
end
