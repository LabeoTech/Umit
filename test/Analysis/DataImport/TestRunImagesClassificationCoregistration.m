classdef TestRunImagesClassificationCoregistration < matlab.unittest.TestCase
%TESTRUNIMAGESCLASSIFICATIONCOREGISTRATION Functional tests for
%run_ImagesClassification.m's dual-camera coregistration path added this
%session: preferring the rig's active cameraCoregistration managed resource
%over the legacy flat tform file, and correctly recording resourceUUID/rigID.
%
%   Uses a copy of the real dual-camera raw acquisition at
%   D:\UMIT-DEV\Test2Cam (never reads/writes the original -- only the raw
%   .bin/info.txt files are copied into an isolated temp RawFolder per test).
%
%   UMITRigStore.getRigsRoot() is a single, non-redirectable real folder (see
%   test/RigManagement/TestUMITRigStore.m for the same constraint). Each test
%   snapshots UMITRigStore.listRigs() before/after and removes exactly the
%   rig folder(s) that appeared, so real pre-existing rigs are never touched.

    properties (Constant)
        SourceRawFolder = 'D:\UMIT-DEV\Test2Cam';
        RawFileNames = {'info.txt', 'ai_00000.bin', 'img_00000.bin', 'imgCam2_00000.bin'};
    end

    properties
        TempRoot
        RawFolder
        RigRootsToClean = {};
        RigFixture
    end

    methods (TestClassSetup)
        function checkSourceDataAvailable(testCase)
        %CHECKSOURCEDATAAVAILABLE Skip the whole class if the real test data
        %isn't present on this machine (e.g. a different dev machine/CI).

            testCase.assumeTrue(isfolder(testCase.SourceRawFolder), ...
                sprintf('Skipped: test data folder not found: %s', ...
                testCase.SourceRawFolder));
        end
    end

    methods (TestMethodSetup)
        function createIsolatedRawFolder(testCase)
        %CREATEISOLATEDRAWFOLDER Copy only the raw acquisition files needed
        %by ImagesClassification into a fresh temp folder.

            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.RawFolder = fullfile(testCase.TempRoot, 'Raw');
            mkdir(testCase.RawFolder);

            for iFile = 1:numel(testCase.RawFileNames)
                sourceFile = fullfile(testCase.SourceRawFolder, testCase.RawFileNames{iFile});
                testCase.assumeTrue(isfile(sourceFile), ...
                    sprintf('Skipped: expected raw file missing: %s', sourceFile));
                copyfile(sourceFile, fullfile(testCase.RawFolder, testCase.RawFileNames{iFile}));
            end

            testCase.RigRootsToClean = {};

            % DFR-20260819-010: run_ImagesClassification resolves its Rig via
            % UMITRigStore.getOrCreateDefaultRig(), which requires exactly
            % one Active Rig. Guarantee that regardless of ambient state,
            % and restore whatever was Active before in teardown.
            testCase.RigFixture = setupIsolatedActiveRigFixture();
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolders(testCase)
        %REMOVETEMPORARYFOLDERS Remove the temp raw/save folders and exactly
        %the rig folder(s) this test caused to be created.

            teardownIsolatedActiveRigFixture(testCase.RigFixture);
            testCase.removeOwnedDefaultPointer();
            for iRig = 1:numel(testCase.RigRootsToClean)
                rigRoot = testCase.RigRootsToClean{iRig};
                if isfolder(rigRoot)
                    rmdir(rigRoot, 's');
                end
            end

            if isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function testFallsBackGracefullyWithNoActiveRigResource(testCase)
        %TESTFALLSBACKGRACEFULLYWITHNOACTIVERIGRESOURCE Verify import
        %succeeds and records no coregistration when the resolved
        %(auto-provisioned) rig has no active cameraCoregistration resource
        %and no legacy flat file exists either.

            rigsBefore = UMITRigStore.listRigs();
            saveFolder = fullfile(testCase.TempRoot, 'Save1');
            mkdir(saveFolder);

            classifiedFiles = run_ImagesClassification( ...
                testCase.RawFolder, saveFolder, 'backupOpts', 'ERASE');
            testCase.trackNewRigs(rigsBefore);

            testCase.verifyNotEmpty(classifiedFiles);
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'DataParams.mat')));

            S = load(fullfile(saveFolder, 'DataParams.mat'));
            testCase.verifyEqual(S.DataParams.folders.RawFolder, testCase.RawFolder);
            testCase.verifyFalse(S.DataParams.cameraCoregistration.isCoregistered);
            testCase.verifyEqual(S.DataParams.cameraCoregistration.resourceUUID, '');
        end

        function testPrefersActiveRigResourceOverLegacyFallback(testCase)
        %TESTPREFERSACTIVERIGRESOURCEOVERLEGACYFALLBACK Core regression test
        %for this session's run_ImagesClassification.m fix: when the rig has
        %an active cameraCoregistration resource, import must use it (and
        %record its real resourceUUID/rigID), not the legacy flat file.

            % DFR-20260819-010: TestMethodSetup's RigFixture already
            % guarantees a restorable Active Rig, so it is safe to switch
            % the Default Rig again here regardless of ambient state.
            rigInfo = struct();
            rigInfo.rigID = ['UnitTestRig_' char(java.util.UUID.randomUUID())];
            rigStore = UMITRigStore.create(rigInfo);
            testCase.RigRootsToClean{end+1} = rigStore.RigRoot;
            UMITRigStore.setDefaultRig(rigStore.getRigInfo().uuid);

            [tformSourceFile, expectedTform] = testCase.createSyntheticTformSource();
            resourceUUID = rigStore.addCameraCoregistration(tformSourceFile, ...
                struct('displayName', 'Unit test coregistration'));
            rigStore.setActiveCameraCoregistration(resourceUUID);

            saveFolder = fullfile(testCase.TempRoot, 'Save2');
            mkdir(saveFolder);

            classifiedFiles = run_ImagesClassification( ...
                testCase.RawFolder, saveFolder, 'backupOpts', 'ERASE'); %#ok<NASGU>

            S = load(fullfile(saveFolder, 'DataParams.mat'));
            cc = S.DataParams.cameraCoregistration;

            testCase.verifyTrue(cc.isCoregistered);
            testCase.verifyEqual(cc.resourceUUID, resourceUUID);
            testCase.verifyEqual(cc.rigID, rigInfo.rigID);
            testCase.verifyEqual(cc.method, 'run_ImagesClassification');
            testCase.verifyEqual(cc.confirmationMode, 'automatic-import-application');

            % Confirm the transform was actually applied to Camera 2 data:
            % applyTform2Cams writes tformDualCam.mat only on successful
            % application.
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'tformDualCam.mat')));

            appliedT = load(fullfile(saveFolder, 'tformDualCam.mat'));
            testCase.verifyClass(appliedT.tform, 'affine2d');
            testCase.verifyEqual(appliedT.tform.T, expectedTform.T, 'AbsTol', 1e-9);
        end
    end

    methods (Access = private)
        function trackNewRigs(testCase, rigsBefore)
        %TRACKNEWRIGS Diff listRigs() against a prior snapshot and track any
        %newly-appeared rig folder(s) for teardown cleanup (covers
        %getOrCreateDefaultRig auto-provisioning a real rig on this machine).

            rigsAfter = UMITRigStore.listRigs();
            if isempty(rigsBefore)
                newRoots = cellstr(rigsAfter.RigRoot);
            else
                newRoots = setdiff(cellstr(rigsAfter.RigRoot), cellstr(rigsBefore.RigRoot));
            end
            testCase.RigRootsToClean = [testCase.RigRootsToClean, newRoots(:)'];
        end

        function removeOwnedDefaultPointer(testCase)
            schema = getUMITRigSchema();
            defaultFile = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.store.internalFolder, schema.store.defaultFile);
            if ~isfile(defaultFile)
                return
            end
            loaded = load(defaultFile, schema.store.defaultVariable, '-mat');
            if ~isfield(loaded, schema.store.defaultVariable)
                return
            end
            defaultInfo = loaded.(schema.store.defaultVariable);
            for iRig = 1:numel(testCase.RigRootsToClean)
                rigFile = fullfile(testCase.RigRootsToClean{iRig}, schema.files.rigMetadata);
                if ~isfile(rigFile)
                    continue
                end
                rigLoaded = load(rigFile, schema.metadataVariables.rig, '-mat');
                if strcmpi(defaultInfo.rigUUID, rigLoaded.RigInfo.uuid)
                    delete(defaultFile);
                    return
                end
            end
        end

        function [filePath, tform] = createSyntheticTformSource(testCase)
        %CREATESYNTHETICTFORMSOURCE Build a small, real geometric transform
        %(matching genTform2Cams's output shape) and save it as a
        %tform/tformInfo .mat file suitable for addCameraCoregistration.

            tform = affine2d([1 0 0; 0 1 0; 3 -2 1]); % small pure translation
            tformInfo = struct( ...
                'RegisteredImages', {{single(zeros(4,5)), single(zeros(4,5))}}, ...
                'OriginalImages', {{single(zeros(4,5)), single(zeros(4,5))}}, ...
                'Binning', 1, 'BinningSpatial', 1, ...
                'Rotation', 0, 'X_Offset', 0, 'Y_Offset', 0);

            filePath = fullfile(testCase.TempRoot, 'synthetic_tform.mat');
            save(filePath, 'tform', 'tformInfo', '-mat');
        end
    end
end
