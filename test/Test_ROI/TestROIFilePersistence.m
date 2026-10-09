classdef TestROIFilePersistence < matlab.unittest.TestCase
%TESTROIFILEPERSISTENCE Persistence, schema, and migration tests for ROI files.

    properties
        TempFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            testCase.TempFolder = tempname;
            mkdir(testCase.TempFolder);
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolder(testCase)
            if isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testSaveLoadRoundTripPreservesGeometryAndStripsRuntime(testCase)
            ROIFile = createROIFile([8 10], ...
                'dataFile', 'green.dat', ...
                'statsImage', zeros(8, 10, 'single'), ...
                'ROIs', iMakeRectangleROI([8 10]));
            ROIFile.ROIs.ID = 'app-only-id';
            ROIFile.ROIs.runtime = struct('selected', true);

            filePath = fullfile(testCase.TempFolder, 'example.roi');
            saveROIFile(filePath, ROIFile);
            loaded = loadROIFile(filePath);

            testCase.verifyEqual(loaded.schemaName, 'UMIT_ROI');
            testCase.verifyEqual(loaded.version, 1);
            testCase.verifyEqual(loaded.imageInfo.imageSizeYX, [8 10]);
            testCase.verifyEqual(loaded.imageInfo.dataFile, 'green.dat');
            testCase.verifyClass(loaded.ROIs.geometry.polyshape, 'polyshape');
            testCase.verifyEqual(loaded.ROIs.mask, ROIFile.ROIs.mask);
            testCase.verifyEqual(loaded.ROIs.stats.NPixels, ...
                ROIFile.ROIs.stats.NPixels);
            testCase.verifyFalse(isfield(loaded.ROIs, 'ID'));
            testCase.verifyFalse(isfield(loaded.ROIs, 'runtime'));
            testCase.verifyWarningFree(@() validateROIFile(loaded));
        end

        function testInvalidMaskPayloadIsRejected(testCase)
            ROIFile = createROIFile([8 10], ...
                'ROIs', iMakeRectangleROI([8 10]));
            ROIFile.ROIs.mask = false(7, 10);

            testCase.verifyError(@() validateROIFile(ROIFile), ...
                'validateROIFile:MaskSizeMismatch');
        end

        function testDuplicateNamesAndInvalidSchemaAreRejected(testCase)
            roi = iMakeRectangleROI([8 10]);
            roi(2) = roi(1);
            roi(2).name = roi(1).name;
            ROIFile = createROIFile([8 10]);
            ROIFile.ROIs = roi;

            testCase.verifyError(@() validateROIFile(ROIFile), ...
                'validateROIFile:DuplicateROINames');

            ROIFile = createROIFile([8 10]);
            ROIFile.schemaName = 'wrong-schema';
            testCase.verifyError(@() validateROIFile(ROIFile), ...
                'validateROIFile:InvalidSchemaName');
        end

        function testCurrentMigrationIsStableAndUnsupportedVersionsFail(testCase)
            ROIFile = createROIFile([8 10]);
            migrated = migrateROIFile(ROIFile);
            testCase.verifyEqual(migrated, ROIFile);

            legacy = rmfield(ROIFile, 'version');
            testCase.verifyError(@() migrateROIFile(legacy), ...
                'migrateROIFile:MissingVersion');

            unsupported = ROIFile;
            unsupported.version = 2;
            testCase.verifyError(@() migrateROIFile(unsupported), ...
                'migrateROIFile:UnsupportedVersion');
        end

        function testMalformedROIFileVariableFailsOnLoad(testCase)
            filePath = fullfile(testCase.TempFolder, 'malformed.roi');
            NotROIFile = struct('version', 1);
            save(filePath, '-mat', 'NotROIFile');

            testCase.verifyError(@() loadROIFile(filePath), ...
                'loadROIFile:MissingROIFileVariable');
        end
    end
end

function roi = iMakeRectangleROI(imageSizeYX)
%IMAKERECTANGLEROI Build one deterministic schema-valid rectangle ROI.

height = imageSizeYX(1);
width = imageSizeYX(2);
vertices = [2 2; width - 1 2; width - 1 height - 2; 2 height - 2];
mask = poly2mask(vertices(:, 1), vertices(:, 2), height, width);
pgon = polyshape(vertices(:, 1), vertices(:, 2));

stats = struct( ...
    'computedOn', datetime('now'), ...
    'NPixels', nnz(mask), ...
    'areaPx2', nnz(mask), ...
    'areaMM2', [], ...
    'centroidXY_px', [], ...
    'centroidXY_mm', [], ...
    'distanceFromOrigin_px', [], ...
    'distanceFromOrigin_mm', [], ...
    'spatialMean', [], ...
    'spatialStd', [], ...
    'spatialMedian', [], ...
    'spatialMin', [], ...
    'spatialMax', []);

roi = struct( ...
    'name', 'RectangleROI', ...
    'type', 'polygon', ...
    'DOC', datetime('now'), ...
    'modifiedOn', datetime('now'), ...
    'color', [1 0 0], ...
    'notes', '', ...
    'geometry', struct( ...
        'polyshape', pgon, ...
        'verticesXY_px', vertices, ...
        'ROIType', 'polygon', ...
        'ROIParameters', struct('ROIType', 'polygon', ...
            'Position', vertices, 'Vertices', vertices)), ...
    'mask', logical(mask), ...
    'stats', stats);
end
