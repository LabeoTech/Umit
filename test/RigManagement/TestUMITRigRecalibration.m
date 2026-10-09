classdef TestUMITRigRecalibration < matlab.unittest.TestCase
    %TESTUMITRIGRECALIBRATION Focused active-resource transition tests.

    properties
        RigRoot = ''
        SourceFolder = ''
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            testCase.SourceFolder = tempname;
            mkdir(testCase.SourceFolder);

            suffix = char(java.util.UUID.randomUUID());
            store = UMITRigStore.create(struct( ...
                'rigID', ['RecalibrationTest_' suffix(1:8)]));
            testCase.RigRoot = store.RigRoot;
        end
    end

    methods (TestMethodTeardown)
        function removeFixture(testCase)
            if isfolder(testCase.RigRoot)
                rmdir(testCase.RigRoot, 's');
            end
            if isfolder(testCase.SourceFolder)
                rmdir(testCase.SourceFolder, 's');
            end
        end
    end

    methods (Test)
        function testDeactivateAndRestore(testCase)
            store = testCase.openStore();
            oldUUID = store.addCameraCoregistration( ...
                testCase.createTform('old.mat', 1), struct());

            store.clearActiveCameraCoregistration();
            testCase.verifyEmpty(store.getActiveCameraCoregistration());
            testCase.verifyEqual(store.getResource(oldUUID).status, 'available');

            store.setActiveCameraCoregistration(oldUUID);
            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, oldUUID);
        end

        function testSuccessfulReplacementPreservesHistory(testCase)
            store = testCase.openStore();
            oldUUID = store.addCameraCoregistration( ...
                testCase.createTform('old.mat', 1), struct());
            store.clearActiveCameraCoregistration();
            newUUID = store.addCameraCoregistration( ...
                testCase.createTform('new.mat', 2), struct());
            store.setActiveCameraCoregistration(newUUID);

            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, newUUID);
            testCase.verifyEqual(store.getResource(oldUUID).status, 'available');
            testCase.verifyTrue(isfile(store.resolveResourcePath(oldUUID)));
            testCase.verifyEqual(numel(store.listResources( ...
                'Type', 'cameraCoregistration')), 2);
        end

        function testManagerResourceLifecycle(testCase)
            % The GUI manager deliberately composes the existing store APIs;
            % exercise the exact state transitions it exposes.
            store = testCase.openStore();
            activeUUID = store.addCameraCoregistration( ...
                testCase.createTform('active.mat', 1), struct( ...
                'displayName', 'Active calibration'));
            availableUUID = store.addCameraCoregistration( ...
                testCase.createTform('available.mat', 2), struct( ...
                'displayName', 'Available calibration'));
            store.setActiveCameraCoregistration(activeUUID);

            testCase.verifyError(@() store.archiveResource(activeUUID), ...
                'Umitoolbox:UMITRigStore:archiveFailed');

            store.archiveResource(availableUUID);
            archivedPath = store.resolveResourcePath(availableUUID);
            testCase.verifyEqual(store.getResource(availableUUID).status, ...
                'archived');
            testCase.verifyTrue(contains(archivedPath, ...
                [filesep 'archive' filesep]));
            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, activeUUID);

            store.restoreResource(availableUUID);
            restoredPath = store.resolveResourcePath(availableUUID);
            testCase.verifyEqual(store.getResource(availableUUID).status, ...
                'available');
            testCase.verifyFalse(contains(restoredPath, ...
                [filesep 'archive' filesep]));
            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, activeUUID);

            store.setActiveCameraCoregistration(availableUUID);
            testCase.verifyEqual(store.getResource(activeUUID).status, ...
                'available');
            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, availableUUID);

            store.archiveResource(activeUUID);
            store.purgeResource(activeUUID);
            remaining = store.listResources( ...
                'Type', 'cameraCoregistration');
            testCase.verifyFalse(any(strcmp({remaining.uuid}, activeUUID)));
        end
    end

    methods (Access = private)
        function store = openStore(testCase)
            info = load(fullfile(testCase.RigRoot, 'rig.mat'), 'RigInfo');
            store = UMITRigStore.open(info.RigInfo.uuid);
        end

        function path = createTform(testCase, fileName, seed)
            path = fullfile(testCase.SourceFolder, fileName);
            tform = affine2d([1 0 0; 0 1 0; seed 0 1]);
            tformInfo = struct('OriginalImages', ...
                {{single(ones(3)), single(ones(3))}});
            save(path, 'tform', 'tformInfo', '-mat');
        end

    end
end
