classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestDataParams < matlab.unittest.TestCase
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

        function testImageSizeComesFromTheDataNotAcqInfos(testCase)
            % .dat header Phase 7a: AcqInfos.mat Height/Width describe the
            % raw acquisition; the default size comes from the .dat headers.
            testCase.writeAcquisitionInfo(50, 70);
            writeTestDat(fullfile(testCase.TempFolder, 'green.dat'), ...
                zeros(5, 7, 3, 'single'), 10);

            DataParams = testCase.verifyWarningFree(@() createDataParams(testCase.TempFolder));

            testCase.verifyEqual(DataParams.view.imageSizeYX, [5 7]);
            testCase.verifySize(DataParams.mask.logical, [5 7]);
        end

        function testAcqInfosSizeIsAWarnedLegacyFallback(testCase)
            testCase.writeAcquisitionInfo(5, 7);
            DataParams = testCase.verifyWarning(@() createDataParams(testCase.TempFolder), ...
                'createDataParams:SizeFromAcqInfos');
            testCase.verifyEqual(DataParams.view.imageSizeYX, [5 7]);
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

        function testRegistrationRecordVersion2Layout(testCase)
            % Project redesign Phase 2: DataParams.registration is the
            % version-2 record {version, reference, estimate, apply}; the
            % version-1 provenance fields are gone.
            DataParams = createDataParams(testCase.TempFolder);
            registration = DataParams.registration;

            testCase.verifyEqual(registration.version, 2);
            testCase.verifyEmpty(registration.reference);
            testCase.verifyEmpty(registration.estimate);
            testCase.verifyFalse(registration.apply.isRegistered);
            testCase.verifyEqual(sort(fieldnames(registration))', ...
                sort({'version', 'reference', 'estimate', 'apply'}));

            DataParams.registration.imageReferenceUUID = '550e8400-e29b-41d4-a716-446655440000';
            testCase.verifyError(@() saveDataParams(testCase.TempFolder, DataParams), ...
                'validateDataParams:InvalidRegistration');
        end

        function testRegistrationRecordRules(testCase)
            DataParams = createDataParams(testCase.TempFolder);
            [~, referenceTemplate, estimateTemplate] = newRegistrationRecord();
            image = single(magic(6));
            reference = referenceTemplate;
            [reference.contentHash, reference.image] = computeImageContentHash(image);
            reference.imageSizeYX = [6 6];
            reference.name = 'test';
            reference.source = referenceSourceDescriptor('file', 'C:\ref.tif');
            estimate = estimateTemplate;
            estimate.tform = affine2d(eye(3));
            estimate.referenceHash = reference.contentHash;
            DataParams.registration.reference = reference;
            DataParams.registration.estimate = estimate;
            testCase.verifyWarningFree(@() saveDataParams(testCase.TempFolder, DataParams));

            % Apply fields only once aligned.
            bad = DataParams;
            bad.registration.apply.appliedOn = '2026-10-15 10:00:00';
            testCase.verifyError(@() validateDataParams(bad), ...
                'validateDataParams:AppliedFieldWithoutRegistration');

            % Aligned: the commitment must equal the reference hash.
            aligned = DataParams;
            aligned.view.imageSizeYX = [6 6];
            aligned.mask.logical = true(6);
            aligned.registration.apply.isRegistered = true;
            aligned.registration.apply.appliedReferenceHash = reference.contentHash;
            testCase.verifyWarningFree(@() validateDataParams(aligned));
            wrongHash = aligned;
            wrongHash.registration.apply.appliedReferenceHash = ['sha256:' repmat('0', 1, 64)];
            testCase.verifyError(@() validateDataParams(wrongHash), ...
                'validateDataParams:AppliedReferenceMismatch');
            noEstimate = aligned;
            noEstimate.registration.estimate = [];
            testCase.verifyError(@() validateDataParams(noEstimate), ...
                'validateDataParams:MissingTform');

            badHash = DataParams;
            badHash.registration.reference.contentHash = 'abc';
            testCase.verifyError(@() validateDataParams(badHash), ...
                'validateDataParams:InvalidContentHash');
            badSource = DataParams;
            badSource.registration.reference.source = rmfield(reference.source, 'warning');
            testCase.verifyError(@() validateDataParams(badSource), ...
                'validateDataParams:InvalidReferenceSource');
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
