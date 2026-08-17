classdef TestDataParams < matlab.unittest.TestCase
%TESTDATAPARAMS Focused persistence and invariant tests for DataParams.

    properties
        TempFolder
        RawFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolders(testCase)
            testCase.TempFolder = tempname;
            mkdir(testCase.TempFolder);
            testCase.RawFolder = fullfile(testCase.TempFolder, 'Raw');
            mkdir(testCase.RawFolder);
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolders(testCase)
            if isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testCreatesSpatialDefaultsFromAcquisitionMetadata(testCase)
            testCase.writeAcquisitionInfo(5, 7);

            DataParams = createDataParams(testCase.TempFolder);

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, ...
                'DataParams.mat')));
            testCase.verifyEqual(DataParams.view.imageSizeYX, [5 7]);
            testCase.verifySize(DataParams.mask.logical, [5 7]);
            testCase.verifyTrue(all(DataParams.mask.logical, 'all'));
            testCase.verifyEqual(DataParams.folders.RawFolder, 'Missing');
            testCase.verifyEqual(DataParams.folders.RawFolderStatus, 'missing');
            testCase.verifyWarningFree(@() validateDataParams(DataParams));
        end

        function testCreatesGracefullyWithoutOptionalAcquisitionMetadata(testCase)
            DataParams = createDataParams(testCase.TempFolder);

            testCase.verifyEmpty(DataParams.view.imageSizeYX);
            testCase.verifyEmpty(DataParams.mask.logical);
            testCase.verifyEqual(DataParams.folders.RawFolder, 'Missing');
            testCase.verifyEqual(DataParams.folders.RawFolderStatus, 'missing');
            testCase.verifyWarningFree(@() validateDataParams(DataParams));
        end

        function testLoadsLegacyRawFolderAndNormalizesFolderMetadata(testCase)
            DataParams = createDataParams(testCase.TempFolder);
            DataParams = rmfield(DataParams, 'folders');
            DataParams.RawFolder = testCase.RawFolder;
            save(fullfile(testCase.TempFolder, 'DataParams.mat'), ...
                'DataParams', '-mat');

            loaded = loadDataParams(testCase.TempFolder);

            testCase.verifyEqual(loaded.folders.RawFolder, testCase.RawFolder);
            testCase.verifyEqual(loaded.folders.RawFolderStatus, 'valid');
            testCase.verifyEmpty(loaded.folders.RawFolderSetOn);
            testCase.verifyEqual(loaded.folders.RawFolderSetBy, '');
            testCase.verifyWarningFree(@() validateDataParams(loaded));
        end

        function testRawFolderUpdateTracksStatusAndRoundTrips(testCase)
            testCase.writeAcquisitionInfo(4, 6);
            createDataParams(testCase.TempFolder);

            mask = false(4, 6);
            mask(2:3, 2:5) = true;
            updateDataParam(testCase.TempFolder, 'mask.logical', mask);
            updateDataParam(testCase.TempFolder, 'folders.RawFolder', ...
                testCase.RawFolder);

            loaded = loadDataParams(testCase.TempFolder);
            testCase.verifyEqual(loaded.mask.logical, mask);
            testCase.verifyEqual(loaded.folders.RawFolder, testCase.RawFolder);
            testCase.verifyEqual(loaded.folders.RawFolderStatus, 'valid');
            testCase.verifyEqual(loaded.folders.RawFolderSetBy, ...
                'updateDataParam');
            testCase.verifyClass(loaded.folders.RawFolderSetOn, 'datetime');

            updateDataParam(testCase.TempFolder, 'folders.RawFolder', ...
                fullfile(testCase.TempFolder, 'does-not-exist'));
            loaded = loadDataParams(testCase.TempFolder);
            testCase.verifyEqual(loaded.folders.RawFolderStatus, 'invalid');
        end

        function testRegistrationProvenanceRequiresCompletePair(testCase)
            DataParams = createDataParams(testCase.TempFolder);
            DataParams.registration.imageReferenceUUID = ...
                '550e8400-e29b-41d4-a716-446655440000';

            testCase.verifyError(@() saveDataParams( ...
                testCase.TempFolder, DataParams), ...
                'validateDataParams:IncompleteImageReferenceProvenance');

            DataParams.registration.imageReferenceChecksum = repmat('a', 1, 64);
            testCase.verifyWarningFree(@() saveDataParams( ...
                testCase.TempFolder, DataParams));

            loaded = loadDataParams(testCase.TempFolder);
            testCase.verifyEqual(loaded.registration.imageReferenceUUID, ...
                DataParams.registration.imageReferenceUUID);
            testCase.verifyEqual(loaded.registration.imageReferenceChecksum, ...
                DataParams.registration.imageReferenceChecksum);
        end

        function testRejectsMaskSizeMismatch(testCase)
            testCase.writeAcquisitionInfo(5, 7);
            DataParams = createDataParams(testCase.TempFolder);
            DataParams.mask.logical = false(4, 7);

            testCase.verifyError(@() validateDataParams(DataParams), ...
                'validateDataParams:MaskSizeMismatch');
        end

        function testMissingDataParamsFileFailsClearly(testCase)
            testCase.verifyError(@() loadDataParams(testCase.TempFolder), ...
                'loadDataParams:MissingFile');
        end
    end

    methods (Access = private)
        function writeAcquisitionInfo(testCase, height, width)
            AcqInfoStream = struct( ...
                'Height', height, ...
                'Width', width);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
        end
    end
end
