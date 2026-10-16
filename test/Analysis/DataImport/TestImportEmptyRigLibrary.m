classdef (SharedTestFixtures = {UMITTemporaryRootFixture('emptyRigLibrary')}) TestImportEmptyRigLibrary < matlab.unittest.TestCase
%TESTIMPORTEMPTYRIGLIBRARY Imports on a machine whose rig library is empty.
%
%   Rig-ledger redesign, Phase 3: no import path may resolve, create or
%   activate a rig on a machine without rigs. Single-camera and TIFF
%   imports record an explicit "no rig"; dual-camera data are imported
%   without coregistration, with a warning (decision Q1).
%
%   The class has its own UMITTemporaryRootFixture (tag 'emptyRigLibrary'),
%   so the rig library is a fresh, empty temporary folder whatever the
%   machine holds. The assumeEqual guard remains as a safety net.

    properties (Constant)
        DualCameraRawFolder = 'D:\UMIT-DEV\Test2Cam';
        DualCameraFiles = {'info.txt', 'ai_00000.bin', 'img_00000.bin', 'imgCam2_00000.bin'};
    end

    properties
        TempRoot
    end

    methods (TestMethodSetup)
        function requireEmptyLibrary(testCase)
            testCase.assumeEqual(height(UMITRigStore.listRigs()), 0, ...
                'Skipped: the rig library is not empty on this machine.');
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolders(testCase)
            if ~isempty(testCase.TempRoot) && isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function testSingleCameraImportWritesNoRig(testCase)
            fixtureFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'TestingData_with_events');
            saveFolder = fullfile(testCase.TempRoot, 'SingleCam');
            mkdir(saveFolder);

            [~, rigResolution] = ImagesClassification(fixtureFolder, saveFolder, ...
                'backupOpts', 'ERASE');

            testCase.verifyEqual(UMITRigLedger.read(saveFolder).state, 'NoRig');
            testCase.verifyEqual(rigResolution.ledgerState, 'NoRig');
            testCase.verifyEqual(height(UMITRigStore.listRigs()), 0);
        end

        function testDualCameraImportIsUncoregisteredWithNoRig(testCase)
            testCase.assumeTrue(isfolder(testCase.DualCameraRawFolder), ...
                'Skipped: dual-camera test data not found.');
            rawFolder = fullfile(testCase.TempRoot, 'Raw2Cam');
            mkdir(rawFolder);
            for iFile = 1:numel(testCase.DualCameraFiles)
                copyfile(fullfile(testCase.DualCameraRawFolder, testCase.DualCameraFiles{iFile}), ...
                    rawFolder);
            end
            saveFolder = fullfile(testCase.TempRoot, 'DualCam');
            mkdir(saveFolder);

            [~, rigResolution] = testCase.verifyWarning( ...
                @() ImagesClassification(rawFolder, saveFolder, 'backupOpts', 'ERASE'), ...
                'Umitoolbox:ImagesClassification:ImportedWithoutCoregistration');

            testCase.verifyEqual(UMITRigLedger.read(saveFolder).state, 'NoRig');
            testCase.verifyEqual(rigResolution.cameraCoregistration.skipReason, 'noActiveRig');
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'tformDualCam.mat')));
            testCase.verifyEqual(height(UMITRigStore.listRigs()), 0);
        end

        function testTifImportWritesNoRig(testCase)
            rawFolder = fullfile(testCase.TempRoot, 'RawTif');
            saveFolder = fullfile(testCase.TempRoot, 'Tif');
            mkdir(rawFolder);
            mkdir(saveFolder);
            data = reshape(uint16(1:4*5*3), [4, 5, 3]);
            for iFrame = 1:size(data, 3)
                writeMode = 'append';
                if iFrame == 1
                    writeMode = 'overwrite';
                end
                imwrite(data(:, :, iFrame), fullfile(rawFolder, 'img_red.tif'), 'tif', ...
                    'WriteMode', writeMode, 'Compression', 'none');
            end
            meta = struct('DateTime', '20261014_120000', 'Tiffiles', struct( ...
                'filename', 'img_red.tif', 'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, 'IlluminationColor', 'Red'));
            fid = fopen(fullfile(rawFolder, 'info.json'), 'w');
            fwrite(fid, jsonencode(meta), 'char');
            fclose(fid);

            [outFile, rigResolution] = importFromTif(rawFolder, saveFolder);

            testCase.verifyEqual(outFile, {'red.dat'});
            testCase.verifyEqual(rigResolution.ledgerState, 'NoRig');
            testCase.verifyEqual(height(UMITRigStore.listRigs()), 0);
        end
    end
end
