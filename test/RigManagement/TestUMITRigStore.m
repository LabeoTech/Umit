classdef TestUMITRigStore < matlab.unittest.TestCase
%TESTUMITRIGSTORE Unit tests for centralized rig management (RigManagement/UMITRigStore.m).
%
%   Covers rig creation/discovery, the cameraCoregistration managed-resource
%   lifecycle (add/setActive/archive/restore), and the new purgeResource method
%   added for the DataViewer_Coreg2Cams UMITRigStore migration.
%
%   NOTE: UMITRigStore.getRigsRoot() resolves to a single, non-configurable real
%   folder (getUmitFolder('rigs')), same constraint as UMITProjectStore. Like
%   test/ProjectManagement/TestUMITProjectStore.m, these tests create real rigs
%   under that real root and track the exact RigRoot/RigUUID created so teardown
%   can remove precisely that folder -- never anything else under the root.
%   getOrCreateDefaultRig tests additionally guard with assumeTrue against the
%   real environment already having rigs, so they self-skip rather than
%   interfere with real user rigs if run on a machine that already has some.

    properties
        TempRoot
        RigRoot
        RigRoot2
        OriginalActiveRigUUID
        SourceFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolders(testCase)
        %CREATETEMPORARYFOLDERS Create isolated folders for each test.

            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.SourceFolder = fullfile(testCase.TempRoot, 'Sources');
            mkdir(testCase.SourceFolder);
            testCase.RigRoot = '';
            testCase.RigRoot2 = '';
            testCase.OriginalActiveRigUUID = '';
            lifecycle = UMITRigStore.validateLifecycle();
            if lifecycle.activeCount == 1
                testCase.OriginalActiveRigUUID = ...
                    UMITRigStore.getActiveRig().getRigInfo().uuid;
            end
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolders(testCase)
        %REMOVETEMPORARYFOLDERS Remove exactly the rig folder(s) this test
        %created, plus the temp source folder. Never touches the rest of the
        %real rigs root.

            if ~isempty(testCase.OriginalActiveRigUUID) && ...
                    UMITRigStore.rigExists(testCase.OriginalActiveRigUUID)
                try
                    originalActive = UMITRigStore.open( ...
                        testCase.OriginalActiveRigUUID);
                    if strcmp(originalActive.getRigInfo().status, 'archived')
                        originalActive.restoreRig();
                    end
                    originalActive.activateRig();
                catch
                    % Preserve the test failure; cleanup is best effort.
                end
            end

            testCase.removeOwnedDefaultPointer();
            if ~isempty(testCase.RigRoot) && isfolder(testCase.RigRoot)
                rmdir(testCase.RigRoot, 's');
            end
            if ~isempty(testCase.RigRoot2) && isfolder(testCase.RigRoot2)
                rmdir(testCase.RigRoot2, 's');
            end
            if isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function testCreateRig(testCase)
        %TESTCREATERIG Verify canonical rig creation.

            store = testCase.createRig();
            RigInfo = store.getRigInfo();
            schema = getUMITRigSchema();

            testCase.verifyEqual(fileparts(testCase.RigRoot), ...
                UMITRigStore.getRigsRoot());
            testCase.verifyTrue(isfile(fullfile( ...
                testCase.RigRoot, schema.files.rigMetadata)));
            testCase.verifyFalse(store.IsReadOnly);
            testCase.verifyNotEmpty(RigInfo.uuid);
            testCase.verifyEqual(RigInfo.activeCoregistrationUUID, '');
            testCase.verifyFalse(isfield(RigInfo, 'activeCalibrationFileUUID'));
            testCase.verifyFalse(isfield(RigInfo, 'isDefault'));

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testOpenAndOpenByRigID(testCase)
        %TESTOPENANDOPENBYRIGID Verify both static lookup paths resolve the
        %same rig.

            store = testCase.createRig();
            RigInfo = store.getRigInfo();
            clear store

            reopenedByUUID = UMITRigStore.open(RigInfo.uuid);
            testCase.verifyEqual(reopenedByUUID.RigRoot, testCase.RigRoot);

            reopenedByID = UMITRigStore.openByRigID(RigInfo.rigID);
            testCase.verifyEqual(reopenedByID.RigRoot, testCase.RigRoot);
        end

        function testTemporaryActivationRestoresActiveRigStatus(testCase)
        %TESTTEMPORARYACTIVATIONRESTORESACTIVERIGSTATUS Regression for test
        %fixtures that used to restore only defaultRig.mat after changing
        %the Active Rig, leaving the former Active Rig persisted as Available.

            original = testCase.createRig();
            original.activateRig();
            originalInfo = original.getRigInfo();

            temporary = testCase.createRig(2);
            temporaryInfo = temporary.getRigInfo();
            fixture = activateRigTemporarily(temporaryInfo.uuid);
            testCase.verifyEqual( ...
                UMITRigStore.getActiveRig().getRigInfo().uuid, ...
                temporaryInfo.uuid);

            deactivateRigTemporarily(fixture);

            restored = UMITRigStore.open(originalInfo.uuid).getRigInfo();
            demoted = UMITRigStore.open(temporaryInfo.uuid).getRigInfo();
            testCase.verifyEqual(restored.status, 'active');
            testCase.verifyEqual(demoted.status, 'available');
            testCase.verifyEqual( ...
                UMITRigStore.getDefaultRig().getRigInfo().uuid, ...
                originalInfo.uuid);
        end

        function testDuplicateSingleCameraRigCopiesConfigurationOnly(testCase)
        %TESTDUPLICATESINGLECAMERARIGCOPIESCONFIGURATIONONLY Verify that a
        %duplicate receives copied editable/configuration data but fresh
        %identity and no default or managed-resource state.

            source = testCase.createConfiguredRigForDuplication(1);
            sourceInfo = source.getRigInfo();
            UMITRigStore.setDefaultRig(sourceInfo.uuid);
            resourceFile = testCase.createTformSource('source-coreg.mat', 1);
            source.addCameraCoregistration(resourceFile, struct());
            sourceWithResource = source.getRigInfo();

            duplicate = source.duplicate( ...
                ['Duplicated_' sourceInfo.rigID], 'DisplayName', 'Duplicated Rig');
            testCase.RigRoot2 = duplicate.RigRoot;
            duplicateInfo = duplicate.getRigInfo();

            testCase.verifyNotEqual(duplicateInfo.uuid, sourceInfo.uuid);
            testCase.verifyEqual(duplicateInfo.displayName, 'Duplicated Rig');
            testCase.verifyEqual(duplicateInfo.description, sourceInfo.description);
            testCase.verifyEqual(duplicateInfo.metadata, sourceInfo.metadata);
            testCase.verifyEqual(duplicateInfo.cameras, sourceInfo.cameras);
            testCase.verifyEqual(duplicateInfo.illuminations, sourceInfo.illuminations);
            testCase.verifyEqual(duplicateInfo.status, 'available');
            testCase.verifyTrue(isnat(duplicateInfo.archivedOn));
            testCase.verifyEqual(duplicateInfo.activeCoregistrationUUID, '');
            testCase.verifyEmpty(duplicateInfo.resourceRegistry);
            copiedResourcePath = fullfile(duplicate.RigRoot, strrep( ...
                sourceWithResource.resourceRegistry.relativePath, '/', filesep));
            testCase.verifyFalse(isfile(copiedResourcePath));

            rigs = UMITRigStore.listRigs();
            duplicateRow = strcmp(rigs.RigUUID, duplicateInfo.uuid);
            testCase.verifyFalse(rigs.IsDefault(duplicateRow));
            testCase.verifyEqual(UMITRigStore.getDefaultRig().getRigInfo().uuid, ...
                sourceInfo.uuid);
            report = duplicate.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testDuplicateDualCameraRigUsesDefaultDisplayName(testCase)
        %TESTDUPLICATEDUALCAMERARIGUSESDEFAULTDISPLAYNAME Verify both camera
        %records survive duplication and the default copy label is sensible.

            source = testCase.createConfiguredRigForDuplication(2);
            sourceInfo = source.getRigInfo();
            duplicate = source.duplicate(['DualCopy_' sourceInfo.rigID]);
            testCase.RigRoot2 = duplicate.RigRoot;
            duplicateInfo = duplicate.getRigInfo();

            testCase.verifyEqual(numel(duplicateInfo.cameras), 2);
            testCase.verifyEqual(duplicateInfo.cameras, sourceInfo.cameras);
            testCase.verifyEqual(duplicateInfo.displayName, ...
                [sourceInfo.displayName ' Copy']);
            testCase.verifyEqual(duplicateInfo.status, 'available');
        end

        function testDuplicateRejectsInvalidOrCollidingRigID(testCase)
        %TESTDUPLICATEREJECTSINVALIDORCOLLIDINGRIGID Verify duplication uses
        %the existing filesystem-backed Rig ID validation and uniqueness.

            source = testCase.createConfiguredRigForDuplication(1);
            sourceInfo = source.getRigInfo();
            testCase.verifyError(@() source.duplicate('invalid rig id'), ...
                'Umitoolbox:UMITRigStore:invalidID');
            testCase.verifyError(@() source.duplicate(sourceInfo.rigID), ...
                'Umitoolbox:UMITRigStore:createFailed');
        end

        function testGetBuiltInRigDefinitionIsHydratedAndReadOnly(testCase)
        %TESTGETBUILTINRIGDEFINITIONISHYDRATEDANDREADONLY Verify the GUI
        %template query returns canonical hardware without creating a Rig.

            before = UMITRigStore.listRigs();
            definition = UMITRigStore.getBuiltInRigDefinition('OiS200');
            after = UMITRigStore.listRigs();

            testCase.verifyEqual(definition.rigID, 'OiS200');
            testCase.verifyEqual(numel(definition.cameras), 2);
            testCase.verifyEqual([definition.cameras.index], [1 2]);
            testCase.verifyEqual( ...
                {definition.illuminations.name}, {'red', 'green', 'yellow'});
            testCase.verifyTrue(all(arrayfun( ...
                @(record) isfield(record, 'manufacturer') && ...
                isfield(record, 'serialNumber'), definition.cameras)));
            testCase.verifyEqual(height(after), height(before));
            testCase.verifyError( ...
                @() UMITRigStore.getBuiltInRigDefinition('Unknown'), ...
                'Umitoolbox:UMITRigStore:unknownBuiltInRig');
        end

        function testSetHardwareConfigurationUpdatesBothCollections(testCase)
        %TESTSETHARDWARECONFIGURATIONUPDATESBOTHCOLLECTIONS Verify one
        %staged hardware apply commits cameras and illuminations together.

            store = testCase.createRig();
            definition = UMITRigStore.getBuiltInRigDefinition('OiS200');
            cameras = definition.cameras(1);
            cameras.displayName = 'Primary Camera';
            illuminations = definition.illuminations;
            illuminations(1).displayName = 'Primary Red';

            store.setHardwareConfiguration(cameras, illuminations);
            actual = store.getRigInfo();

            testCase.verifyEqual(actual.cameras, cameras);
            testCase.verifyEqual(actual.illuminations, illuminations);
            testCase.verifyTrue(store.validate('Mode', 'full').isValid);
        end

        function testSetHardwareConfigurationRejectsAtomically(testCase)
        %TESTSETHARDWARECONFIGURATIONREJECTSATOMICALLY Verify invalid
        %illumination input cannot partially persist valid camera edits.

            store = testCase.createConfiguredRigForDuplication(1);
            before = store.getRigInfo();
            cameras = before.cameras;
            cameras(1).model = 'Should Not Persist';
            illuminations = before.illuminations;
            illuminations(1).spectrumID = 'MissingSpectrum';

            testCase.verifyError(@() store.setHardwareConfiguration( ...
                cameras, illuminations), ...
                'Umitoolbox:UMITRigStore:spectrumNotFound');

            after = store.getRigInfo();
            testCase.verifyEqual(after.cameras, before.cameras);
            testCase.verifyEqual(after.illuminations, before.illuminations);
        end

        function testSetHardwareConfigurationRejectsArchivedRig(testCase)
        %TESTSETHARDWARECONFIGURATIONREJECTSARCHIVEDRIG Verify archived
        %Rigs remain immutable through the combined GUI-facing operation.

            store = testCase.createConfiguredRigForDuplication(1);
            info = store.getRigInfo();
            replacement = testCase.createRig(2);
            store.archiveRig(replacement.getRigInfo().uuid);

            testCase.verifyError(@() store.setHardwareConfiguration( ...
                info.cameras, info.illuminations), ...
                'Umitoolbox:UMITRigStore:archivedRigReadOnly');
        end

        function testListRigsIncludesCreatedRig(testCase)
        %TESTLISTRIGSINCLUDESCREATEDRIG Verify rig enumeration.

            store = testCase.createRig();
            RigInfo = store.getRigInfo();

            rigs = UMITRigStore.listRigs();
            idx = find(strcmp(rigs.RigUUID, RigInfo.uuid), 1, 'first');

            testCase.verifyNotEmpty(idx);
            testCase.verifyTrue(rigs.IsReadable(idx));
            testCase.verifyEqual(char(rigs.RigRoot(idx)), testCase.RigRoot);
        end

        function testGetOrCreateDefaultRigCreatesWhenNoneExist(testCase)
        %TESTGETORCREATEDEFAULTRIGCREATESWHENNONEEXIST Verify auto-provisioning.
        %
        %   Guarded: only runs if the real environment currently has zero
        %   rigs, since getOrCreateDefaultRig operates on the real,
        %   non-redirectable root.

            existingRigs = UMITRigStore.listRigs();
            testCase.assumeTrue(height(existingRigs) == 0, ...
                ['Skipped: real environment already has rig(s); ' ...
                'this scenario is not safely testable without disturbing them.']);

            [store, wasCreated, resolution] = UMITRigStore.getOrCreateDefaultRig();
            testCase.RigRoot = store.RigRoot;

            RigInfo = store.getRigInfo();
            testCase.verifyTrue(wasCreated);
            testCase.verifyEqual(resolution, 'createdDefault');
            testCase.verifyEqual(RigInfo.displayName, 'OiS200 LightTrack Imaging System');

            rigs = UMITRigStore.listRigs();
            testCase.verifyEqual(height(rigs), 1);
            testCase.verifyTrue(rigs.IsDefault(1));
        end

        function testGetOrCreateDefaultRigReusesCreatedDefault(testCase)
        %TESTGETORCREATE... Verify a second imaging-resolution request reuses
        %the default created by the first request.

            existingRigs = UMITRigStore.listRigs();
            testCase.assumeTrue(height(existingRigs) == 0, ...
                ['Skipped: real environment already has rig(s); ' ...
                 'this clean-store scenario is not safely testable without ' ...
                 'disturbing them.']);

            [firstStore, firstWasCreated] = UMITRigStore.getOrCreateDefaultRig();
            testCase.RigRoot = firstStore.RigRoot;
            firstInfo = firstStore.getRigInfo();

            [secondStore, secondWasCreated, secondResolution] = ...
                UMITRigStore.getOrCreateDefaultRig();
            secondInfo = secondStore.getRigInfo();

            testCase.verifyTrue(firstWasCreated);
            testCase.verifyFalse(secondWasCreated);
            testCase.verifyEqual(secondResolution, 'existingDefault');
            testCase.verifyEqual(secondInfo.uuid, firstInfo.uuid);
            testCase.verifyEqual(height(UMITRigStore.listRigs()), 1);
        end

        function testGetOrCreateDefaultRigPromotesSingleExisting(testCase)
        %TESTGETORCREATEDEFAULTRIGPROMOTESSINGLEEXISTING Verify promotion of
        %the sole existing rig to default.
        %
        %   Guarded: only runs if the real environment currently has zero
        %   rigs (so that, after creating exactly one, it is guaranteed to
        %   be the sole rig getOrCreateDefaultRig sees).

            existingRigs = UMITRigStore.listRigs();
            testCase.assumeTrue(height(existingRigs) == 0, ...
                'Skipped: real environment already has rig(s).');

            testCase.createRig();
            [promoted, wasCreated, resolution] = UMITRigStore.getOrCreateDefaultRig();
            testCase.verifyEqual(promoted.RigRoot, testCase.RigRoot);
            testCase.verifyFalse(wasCreated);
            testCase.verifyEqual(resolution, 'promotedExisting');
            testCase.verifyEqual(UMITRigStore.getDefaultRig().RigRoot, testCase.RigRoot);
        end

        function testGetOrCreateDefaultRigErrorsWithMultipleNonDefault(testCase)
        %TESTGETORCREATEDEFAULTRIGERRORSWITHMULTIPLENONDEFAULT Verify ambiguity
        %error when >1 rig exists and none is marked default.
        %
        %   Guarded: only runs if the real environment currently has zero
        %   rigs (so that, after creating exactly two, the "multiple, none
        %   default" precondition is unambiguous).

            existingRigs = UMITRigStore.listRigs();
            testCase.assumeTrue(height(existingRigs) == 0, ...
                'Skipped: real environment already has rig(s).');

            store1 = testCase.createRig();
            store2 = testCase.createRig(2);
            testCase.verifyNotEqual(store1.getRigInfo().uuid, store2.getRigInfo().uuid);

            testCase.verifyError(@() UMITRigStore.getOrCreateDefaultRig(), ...
                'Umitoolbox:UMITRigStore:noDefaultRig');
        end

        function testAddCameraCoregistrationAutoActivatesFirst(testCase)
        %TESTADDCAMERACOREGISTRATIONAUTOACTIVATESFIRST Verify the first
        %imported resource of a type is automatically made active.

            store = testCase.createRig();
            source = testCase.createTformSource('tform1.mat', 1);

            uuid = store.addCameraCoregistration(source, ...
                struct('displayName', 'Coreg 1'));

            active = store.getActiveCameraCoregistration();
            testCase.verifyEqual(active.uuid, uuid);
            testCase.verifyEqual(active.status, 'active');
            testCase.verifyEqual(active.type, 'cameraCoregistration');

            resolvedPath = store.resolveResourcePath(uuid);
            testCase.verifyTrue(isfile(resolvedPath));
            [~, fileName, ext] = fileparts(resolvedPath);
            testCase.verifyTrue(startsWith(fileName, 'CameraCoregistration_'));
            testCase.verifyEqual(ext, '.mat');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testSetActiveCameraCoregistrationDemotesPrevious(testCase)
        %TESTSETACTIVECAMERACOREGISTRATIONDEMOTESPREVIOUS Verify activating a
        %second resource demotes the first (auto-activated on import, being
        %the first ever added) to 'available' (not archived, not deleted).

            store = testCase.createRig();
            source1 = testCase.createTformSource('tform1.mat', 1);
            source2 = testCase.createTformSource('tform2.mat', 2);

            uuid1 = store.addCameraCoregistration(source1, struct());
            uuid2 = store.addCameraCoregistration(source2, struct());

            % uuid1 is auto-activated (first resource of this type ever
            % added); uuid2 lands as 'available' since one is already active.
            testCase.verifyEqual(store.getResource(uuid1).status, 'active');
            testCase.verifyEqual(store.getResource(uuid2).status, 'available');

            store.setActiveCameraCoregistration(uuid2);

            active = store.getActiveCameraCoregistration();
            testCase.verifyEqual(active.uuid, uuid2);
            resource1After = store.getResource(uuid1);
            testCase.verifyEqual(resource1After.status, 'available');
            testCase.verifyTrue(isfile(store.resolveResourcePath(uuid1)));
        end

        function testClearAndRestoreActiveCameraCoregistration(testCase)
        %TESTCLEARANDRESTOREACTIVECAMERACOREGISTRATION Verify recalibration can
        %temporarily deactivate a transform without losing resource history.

            store = testCase.createRig();
            source1 = testCase.createTformSource('tform1.mat', 1);
            source2 = testCase.createTformSource('tform2.mat', 2);

            uuid1 = store.addCameraCoregistration(source1, struct());
            uuid2 = store.addCameraCoregistration(source2, struct());
            store.clearActiveCameraCoregistration();

            testCase.verifyEmpty(store.getActiveCameraCoregistration());
            testCase.verifyEqual(store.getResource(uuid1).status, 'available');
            testCase.verifyEqual(store.getResource(uuid2).status, 'available');
            testCase.verifyTrue(isfile(store.resolveResourcePath(uuid1)));

            store.setActiveCameraCoregistration(uuid1);

            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, uuid1);
            testCase.verifyEqual(store.getResource(uuid1).status, 'active');
            testCase.verifyEqual(store.getResource(uuid2).status, 'available');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testClearActiveCameraCoregistrationIsIdempotent(testCase)
        %TESTCLEARACTIVECAMERACOREGISTRATIONISIDEMPOTENT Clearing an already
        %empty active pointer must not affect available resource history.

            store = testCase.createRig();
            source = testCase.createTformSource('tform1.mat', 1);
            uuid = store.addCameraCoregistration(source, struct());

            store.clearActiveCameraCoregistration();
            store.clearActiveCameraCoregistration();

            testCase.verifyEmpty(store.getActiveCameraCoregistration());
            testCase.verifyEqual(store.getResource(uuid).status, 'available');
        end

        function testSuccessfulRecalibrationPreservesHistory(testCase)
        %TESTSUCCESSFULRECALIBRATIONPRESERVESHISTORY A new validated resource
        %becomes active while the previous transform remains available.

            store = testCase.createRig();
            oldSource = testCase.createTformSource('old_tform.mat', 1);
            newSource = testCase.createTformSource('new_tform.mat', 2);

            oldUUID = store.addCameraCoregistration(oldSource, struct());
            store.clearActiveCameraCoregistration();
            newUUID = store.addCameraCoregistration(newSource, struct());
            store.setActiveCameraCoregistration(newUUID);

            testCase.verifyEqual( ...
                store.getActiveCameraCoregistration().uuid, newUUID);
            testCase.verifyEqual(store.getResource(newUUID).status, 'active');
            testCase.verifyEqual(store.getResource(oldUUID).status, 'available');
            testCase.verifyTrue(isfile(store.resolveResourcePath(oldUUID)));

            history = store.listResources('Type', 'cameraCoregistration');
            testCase.verifyEqual(numel(history), 2);
        end

        function testResourceMovesRecreateMissingStateFolders(testCase)
        %TESTRESOURCEMOVESRECREATEMISSINGSTATEFOLDERS A rig that lost its
        %transforms tree still installs and archives resources (2026-10-05:
        %the real OiS200 rig had no transforms folder).

            store = testCase.createRig();
            schema = getUMITRigSchema();
            transformsRoot = fullfile(store.RigRoot, schema.folders.transforms);
            testCase.assertTrue(isfolder(transformsRoot));
            rmdir(transformsRoot, 's');

            uuid1 = store.addCameraCoregistration( ...
                testCase.createTformSource('tform1.mat', 1), struct());
            uuid2 = store.addCameraCoregistration( ...
                testCase.createTformSource('tform2.mat', 2), struct());
            testCase.verifyTrue(isfile(store.resolveResourcePath(uuid1)));
            testCase.verifyTrue(contains(store.getResource(uuid1).relativePath, 'active'));

            archiveFolder = fullfile(transformsRoot, ...
                schema.folders.cameraCoregistration, schema.folders.archive);
            testCase.assertFalse(isfolder(archiveFolder));
            store.archiveResource(uuid1, 'ReplacementUUID', uuid2);
            testCase.verifyTrue(isfile(store.resolveResourcePath(uuid1)));
            testCase.verifyTrue(contains(store.getResource(uuid1).relativePath, 'archive'));

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testArchiveResourceRequiresReplacementForActive(testCase)
        %TESTARCHIVERESOURCEREQUIRESREPLACEMENTFORACTIVE Verify the active
        %resource cannot be archived without a valid replacement, and that a
        %valid replacement succeeds.

            store = testCase.createRig();
            source1 = testCase.createTformSource('tform1.mat', 1);
            source2 = testCase.createTformSource('tform2.mat', 2);

            uuid1 = store.addCameraCoregistration(source1, struct());
            uuid2 = store.addCameraCoregistration(source2, struct());

            testCase.verifyError( ...
                @() store.archiveResource(uuid1), ...
                'Umitoolbox:UMITRigStore:archiveFailed');

            store.archiveResource(uuid1, 'ReplacementUUID', uuid2);

            resource1 = store.getResource(uuid1);
            resource2 = store.getResource(uuid2);
            testCase.verifyEqual(resource1.status, 'archived');
            testCase.verifyEqual(resource2.status, 'active');
            testCase.verifyEqual(store.getActiveCameraCoregistration().uuid, uuid2);
            testCase.verifyTrue(contains(resource1.relativePath, 'archive'));
        end

        function testRestoreResource(testCase)
        %TESTRESTORERESOURCE Verify an archived resource can be restored to
        %'available' (not automatically active).

            store = testCase.createRig();
            source1 = testCase.createTformSource('tform1.mat', 1);
            source2 = testCase.createTformSource('tform2.mat', 2);

            uuid1 = store.addCameraCoregistration(source1, struct());
            uuid2 = store.addCameraCoregistration(source2, struct());
            store.archiveResource(uuid1, 'ReplacementUUID', uuid2);

            store.restoreResource(uuid1);

            resource1 = store.getResource(uuid1);
            testCase.verifyEqual(resource1.status, 'available');
            testCase.verifyTrue(contains(resource1.relativePath, 'active'));
            testCase.verifyEqual(store.getActiveCameraCoregistration().uuid, uuid2);
        end

        function testPurgeResourceRejectsNonArchived(testCase)
        %TESTPURGERESOURCEREJECTSNONARCHIVED Verify purge refuses 'active' and
        %'available' resources.

            store = testCase.createRig();
            source1 = testCase.createTformSource('tform1.mat', 1);
            source2 = testCase.createTformSource('tform2.mat', 2);

            uuid1 = store.addCameraCoregistration(source1, struct());
            uuid2 = store.addCameraCoregistration(source2, struct());

            % uuid1 is auto-activated (first added); uuid2 is 'available'.
            testCase.verifyEqual(store.getResource(uuid1).status, 'active');
            testCase.verifyError(@() store.purgeResource(uuid1), ...
                'Umitoolbox:UMITRigStore:purgeFailed');

            testCase.verifyEqual(store.getResource(uuid2).status, 'available');
            testCase.verifyError(@() store.purgeResource(uuid2), ...
                'Umitoolbox:UMITRigStore:purgeFailed');

            % Neither rejected purge should have touched the registry or files.
            testCase.verifyTrue(isfile(store.resolveResourcePath(uuid1)));
            testCase.verifyTrue(isfile(store.resolveResourcePath(uuid2)));
        end

        function testPurgeResourceRejectsUnowned(testCase)
        %TESTPURGERESOURCEREJECTSUNOWNED Verify purge errors for a UUID the
        %rig doesn't own.

            store = testCase.createRig();
            randomUUID = lower(char(java.util.UUID.randomUUID()));

            testCase.verifyError(@() store.purgeResource(randomUUID), ...
                'Umitoolbox:UMITRigStore:purgeFailed');
        end

        function testPurgeResourceDeletesFileAndRegistryEntry(testCase)
        %TESTPURGERESOURCEDELETESFILEANDREGISTRYENTRY Core new-feature test:
        %purging an archived resource removes both its file and its registry
        %row, and is the only sanctioned deletion path.

            store = testCase.createRig();
            source1 = testCase.createTformSource('tform1.mat', 1);
            source2 = testCase.createTformSource('tform2.mat', 2);

            uuid1 = store.addCameraCoregistration(source1, struct());
            uuid2 = store.addCameraCoregistration(source2, struct());
            store.archiveResource(uuid1, 'ReplacementUUID', uuid2);

            resourcePath = store.resolveResourcePath(uuid1);
            testCase.verifyTrue(isfile(resourcePath));

            store.purgeResource(uuid1);

            testCase.verifyFalse(isfile(resourcePath));
            testCase.verifyError(@() store.getResource(uuid1), ...
                'Umitoolbox:UMITRigStore:getResourceFailed');

            remaining = store.listResources('Type', 'cameraCoregistration');
            testCase.verifyEqual(numel(remaining), 1);
            testCase.verifyEqual(remaining(1).uuid, uuid2);

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testResolveResourcePathMustExist(testCase)
        %TESTRESOLVERESOURCEPATHMUSTEXIST Verify resolveResourcePath raises
        %when the pointed-to file is missing on disk, and that
        %'MustExist', false suppresses that check.

            store = testCase.createRig();
            source = testCase.createTformSource('tform1.mat', 1);
            uuid = store.addCameraCoregistration(source, struct());

            resourcePath = store.resolveResourcePath(uuid);
            delete(resourcePath);

            testCase.verifyError( ...
                @() store.resolveResourcePath(uuid), ...
                'Umitoolbox:UMITRigStore:resourceFileMissing');

            testCase.verifyEqual( ...
                store.resolveResourcePath(uuid, 'MustExist', false), ...
                resourcePath);
        end

        function testGetActiveCameraCoregistrationEmptyWhenNoneActive(testCase)
        %TESTGETACTIVECAMERACOREGISTRATIONEMPTYWHENNONEACTIVE Verify a fresh
        %rig with no resources reports no active coregistration.

            store = testCase.createRig();
            testCase.verifyEmpty(store.getActiveCameraCoregistration());
        end

        function testUpdateRigMetadataUpdatesBasicFields(testCase)
        %TESTUPDATERIGMETADATAUPDATESBASICFIELDS Verify displayName/description
        %are updated in place and modifiedOn advances.

            store = testCase.createRig();
            before = store.getRigInfo();

            store.updateRigMetadata(struct( ...
                'displayName', 'Renamed Rig', ...
                'description', 'Updated description.'));

            after = store.getRigInfo();
            testCase.verifyEqual(after.displayName, 'Renamed Rig');
            testCase.verifyEqual(after.description, 'Updated description.');
            testCase.verifyGreaterThanOrEqual(after.modifiedOn, before.modifiedOn);
        end

        function testUpdateRigMetadataMergesExtensibleMetadataField(testCase)
        %TESTUPDATERIGMETADATAMERGESEXTENSIBLEMETADATAFIELD Verify the open
        %'metadata' bucket (reserved for future camera names/info, filter
        %spectral profiles, etc.) is shallow-merged rather than replaced, so
        %adding one key never clobbers keys set by an earlier update.

            store = testCase.createRig();
            RigInfo = store.getRigInfo();
            testCase.verifyEqual(RigInfo.metadata, struct());

            store.updateRigMetadata(struct( ...
                'metadata', struct('cameraNames', {{'cam1', 'cam2'}})));
            afterFirst = store.getRigInfo();
            testCase.verifyEqual(afterFirst.metadata.cameraNames, {'cam1', 'cam2'});

            store.updateRigMetadata(struct( ...
                'metadata', struct('filterProfile', 'GFP')));
            afterSecond = store.getRigInfo();
            testCase.verifyEqual(afterSecond.metadata.cameraNames, {'cam1', 'cam2'});
            testCase.verifyEqual(afterSecond.metadata.filterProfile, 'GFP');
        end

        function testUpdateRigMetadataRejectsNonStructMetadata(testCase)
        %TESTUPDATERIGMETADATAREJECTSNONSTRUCTMETADATA Verify a non-struct
        %'metadata' value is rejected rather than silently coerced.

            store = testCase.createRig();
            testCase.verifyError( ...
                @() store.updateRigMetadata(struct('metadata', 'not-a-struct')), ...
                'Umitoolbox:UMITRigStore:updateRigFailed');
        end

        function testArchiveRigMarksArchivedAndPreservesData(testCase)
        %TESTARCHIVERIGMARKSARCHIVEDANDPRESERVESDATA Verify archiveRig is a
        %soft-delete: status flips to 'archived', archivedOn is set, and the
        %rig folder/metadata are left in place (nothing deleted or moved).

            store = testCase.createRig();
            replacement = testCase.createRig(2);
            RigInfo = store.getRigInfo();
            rigID = RigInfo.rigID;
            testCase.verifyEqual(RigInfo.status, 'active');
            testCase.verifyTrue(isnat(RigInfo.archivedOn));

            store.archiveRig(replacement.getRigInfo().uuid);

            after = store.getRigInfo();
            testCase.verifyEqual(after.status, 'archived');
            testCase.verifyFalse(isnat(after.archivedOn));
            testCase.verifyTrue(isfolder(testCase.RigRoot));
            testCase.verifyTrue(isfile(fullfile(testCase.RigRoot, 'rig.mat')));

            reopened = UMITRigStore.openByRigID(rigID);
            testCase.verifyEqual(reopened.getRigInfo().status, 'archived');

            rigs = UMITRigStore.listRigs();
            idx = find(strcmp(rigs.RigUUID, RigInfo.uuid), 1, 'first');
            testCase.verifyEqual(char(rigs.Status(idx)), 'archived');
        end

        function testArchiveRigRejectsAlreadyArchived(testCase)
        %TESTARCHIVERIGREJECTSALREADYARCHIVED Verify archiving twice errors
        %instead of silently no-op'ing, matching archiveResource's style.

            store = testCase.createRig();
            replacement = testCase.createRig(2);
            store.archiveRig(replacement.getRigInfo().uuid);

            testCase.verifyError(@() store.archiveRig(), ...
                'Umitoolbox:UMITRigStore:archivedRigReadOnly');
        end

        function testArchivedRigCannotBecomeDefault(testCase)
        %TESTARCHIVEDRIGCANNOTBECOMEDEFAULT Verify archived rigs cannot be
        %selected for new imaging operations.

            store = testCase.createRig();
            replacement = testCase.createRig(2);
            rigUUID = store.getRigInfo().uuid;
            store.archiveRig(replacement.getRigInfo().uuid);
            testCase.verifyError(@() UMITRigStore.setDefaultRig(rigUUID), ...
                'Umitoolbox:UMITRigStore:setDefaultRigFailed');
        end

        function testRigLifecycleActivationAndAssignmentStates(testCase)
        %TESTRIGLIFECYCLEACTIVATIONANDASSIGNMENTSTATES Verify one Active
        %Rig, Available alternatives, and activation demotion.

            first = testCase.createRig();
            second = testCase.createRig(2);
            testCase.verifyEqual(first.getRigInfo().status, 'active');
            testCase.verifyEqual(second.getRigInfo().status, 'available');
            testCase.verifyTrue(UMITRigStore.validateLifecycle().isValid);

            second.activateRig();
            testCase.verifyEqual(first.getRigInfo().status, 'available');
            testCase.verifyEqual(second.getRigInfo().status, 'active');
            lifecycle = UMITRigStore.validateLifecycle();
            testCase.verifyEqual(lifecycle.activeCount, 1);
            testCase.verifyEqual(lifecycle.activeRigUUID, ...
                second.getRigInfo().uuid);
        end

        function testActiveRigArchiveRequiresReplacement(testCase)
        %TESTACTIVERIGARCHIVEREQUIRESREPLACEMENT Verify archival transition.

            active = testCase.createRig();
            replacement = testCase.createRig(2);
            testCase.verifyError(@() active.archiveRig(), ...
                'Umitoolbox:UMITRigStore:archiveRigFailed');

            active.archiveRig(replacement.getRigInfo().uuid);
            testCase.verifyEqual(active.getRigInfo().status, 'archived');
            testCase.verifyEqual(replacement.getRigInfo().status, 'active');
            testCase.verifyTrue(UMITRigStore.validateLifecycle().isValid);

            active.restoreRig();
            testCase.verifyEqual(active.getRigInfo().status, 'available');
            active.activateRig();
            testCase.verifyEqual(active.getRigInfo().status, 'active');
            testCase.verifyEqual(replacement.getRigInfo().status, 'available');
        end

        function testLegacyDatasetRigMigrationIsOneTime(testCase)
        %TESTLEGACYDATASETRIGMIGRATIONIS ONETIME Verify UUID persistence.

            active = testCase.createRig();
            dataFolder = fullfile(testCase.TempRoot, 'LegacyDat');
            mkdir(dataFolder);
            AcqInfoStream = struct();
            save(fullfile(dataFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');

            [firstInfo, migrated] = ...
                UMITRigStore.ensureDatasetRigAssociation(dataFolder);
            testCase.verifyTrue(migrated);
            testCase.verifyEqual(firstInfo.uuid, active.getRigInfo().uuid);

            second = testCase.createRig(2);
            second.activateRig();
            [secondInfo, migratedAgain] = ...
                UMITRigStore.ensureDatasetRigAssociation(dataFolder);
            testCase.verifyFalse(migratedAgain);
            testCase.verifyEqual(secondInfo.uuid, firstInfo.uuid);

            persisted = load(fullfile(dataFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
            testCase.verifyEqual(persisted.AcqInfoStream.rigUUID, firstInfo.uuid);
        end

        function testAssignDatasetRigOverwritesExistingAssociation(testCase)
        %TESTASSIGNDATASETRIGOVERWRITESEXISTINGASSOCIATION Verify the explicit
        %(re)assignment primitive, unlike ensureDatasetRigAssociation, always
        %overwrites -- and that readDatasetRigAssociation/
        %ensureDatasetRigAssociation both observe the pinned Rig afterward.

            first = testCase.createRig();
            second = testCase.createRig(2);
            dataFolder = fullfile(testCase.TempRoot, 'AssignedDataset');
            mkdir(dataFolder);
            AcqInfoStream = struct();
            save(fullfile(dataFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');

            firstInfo = UMITRigStore.ensureDatasetRigAssociation(dataFolder);
            testCase.verifyEqual(firstInfo.uuid, first.getRigInfo().uuid);

            secondInfo = UMITRigStore.assignDatasetRig( ...
                dataFolder, second.getRigInfo().uuid);
            testCase.verifyEqual(secondInfo.uuid, second.getRigInfo().uuid);

            [rigUUID, rigID] = UMITRigStore.readDatasetRigAssociation(dataFolder);
            testCase.verifyEqual(rigUUID, second.getRigInfo().uuid);
            testCase.verifyEqual(rigID, second.getRigInfo().rigID);

            % A subsequent ensure-call must preserve the explicit assignment,
            % not silently fall back to whichever Rig is Active.
            [preservedInfo, wasMigrated] = ...
                UMITRigStore.ensureDatasetRigAssociation(dataFolder);
            testCase.verifyFalse(wasMigrated);
            testCase.verifyEqual(preservedInfo.uuid, second.getRigInfo().uuid);
        end

        function testAssignDatasetRigRejectsArchivedRig(testCase)
        %TESTASSIGNDATASETRIGREJECTSARCHIVEDRIG Verify an Archived Rig cannot
        %be (re)assigned to a dataset. (Relocated from the equivalent
        %addSession-level scenario in TestUMITProjectStore, now that
        %UMITProjectStore no longer accepts caller-supplied Rig fields at all
        %-- UMITRigStore.assignDatasetRig is the only place left where a
        %caller supplies an explicit target Rig.)

            active = testCase.createRig();
            archived = testCase.createRig(2);
            archived.archiveRig(active.getRigInfo().uuid);
            archivedUUID = archived.getRigInfo().uuid;

            dataFolder = fullfile(testCase.TempRoot, 'ArchivedAssignmentDataset');
            mkdir(dataFolder);
            AcqInfoStream = struct();
            save(fullfile(dataFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');

            testCase.verifyError(@() UMITRigStore.assignDatasetRig( ...
                dataFolder, archivedUUID), ...
                'Umitoolbox:UMITRigStore:archivedRigAssignment');
        end

        function testReadDatasetRigAssociationHandlesMissingOrEmptyDataset(testCase)
        %TESTREADDATASETRIGASSOCIATIONHANDLESMISSINGOREMPTYDATASET Verify the
        %read-only accessor never errors: it returns empty for a folder with
        %no AcqInfos.mat, and empty for one with no Rig association yet.

            missingFolder = fullfile(testCase.TempRoot, 'NoAcqInfosHere');
            mkdir(missingFolder);
            [rigUUID, rigID] = ...
                UMITRigStore.readDatasetRigAssociation(missingFolder);
            testCase.verifyEqual(rigUUID, '');
            testCase.verifyEqual(rigID, '');

            unassociatedFolder = fullfile(testCase.TempRoot, 'UnassociatedDataset');
            mkdir(unassociatedFolder);
            AcqInfoStream = struct();
            save(fullfile(unassociatedFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
            [rigUUID, rigID] = ...
                UMITRigStore.readDatasetRigAssociation(unassociatedFolder);
            testCase.verifyEqual(rigUUID, '');
            testCase.verifyEqual(rigID, '');
        end

        function testRigExistsForLiveRig(testCase)
        %TESTRIGEXISTSFORLIVERIG Verify rigExists resolves both by rigID
        %(what DataParams.cameraCoregistration.rigID stores) and by UUID.

            store = testCase.createRig();
            RigInfo = store.getRigInfo();

            testCase.verifyTrue(UMITRigStore.rigExists(RigInfo.rigID));
            testCase.verifyTrue(UMITRigStore.rigExists(RigInfo.uuid));
        end

        function testRigExistsForArchivedRig(testCase)
        %TESTRIGEXISTSFORARCHIVEDRIG Verify an archived rig still resolves --
        %old sessions that reference it by rigID must not appear broken.

            store = testCase.createRig();
            RigInfo = store.getRigInfo();
            replacement = testCase.createRig(2);
            store.archiveRig(replacement.getRigInfo().uuid);

            testCase.verifyTrue(UMITRigStore.rigExists(RigInfo.rigID));
            testCase.verifyTrue(UMITRigStore.rigExists(RigInfo.uuid));
        end

        function testRigExistsReturnsFalseForUnknownOrInvalidInput(testCase)
        %TESTRIGEXISTSRETURNSFALSEFORUNKNOWNORINVALIDINPUT Verify rigExists
        %never throws -- unresolvable input just returns false.

            randomUUID = lower(char(java.util.UUID.randomUUID()));
            testCase.verifyFalse(UMITRigStore.rigExists('NoSuchRigID'));
            testCase.verifyFalse(UMITRigStore.rigExists(randomUUID));
            testCase.verifyFalse(UMITRigStore.rigExists(''));
            testCase.verifyFalse(UMITRigStore.rigExists(123));
        end

        function testGetOrCreateDefaultRigSkipsArchivedRig(testCase)
        %TESTGETORCREATEDEFAULTRIGSKIPSARCHIVEDRIG Verify an archived rig is
        %never auto-selected as the default. With only an archived rig on
        %disk, getOrCreateDefaultRig treats that as "no rig exists" and
        %auto-provisions a fresh default rig rather than promoting the
        %archived one.
        %
        %   Guarded: only runs if the real environment currently has zero
        %   rigs, like the other getOrCreateDefaultRig tests.

            existingRigs = UMITRigStore.listRigs();
            testCase.assumeTrue(height(existingRigs) == 0, ...
                'Skipped: real environment already has rig(s).');

            store = testCase.createRig();
            replacement = testCase.createRig(2);
            store.archiveRig(replacement.getRigInfo().uuid);

            [newDefault, wasCreated] = UMITRigStore.getOrCreateDefaultRig();
            testCase.RigRoot2 = newDefault.RigRoot;

            newInfo = newDefault.getRigInfo();
            testCase.verifyFalse(wasCreated);
            testCase.verifyEqual(newInfo.status, 'active');
            testCase.verifyEqual(newInfo.uuid, replacement.getRigInfo().uuid);
        end

        function testValidationRejectsEscapingManagedResourcePath(testCase)
            store = testCase.createRig();
            source = testCase.createTformSource('valid.mat', 1);
            resourceUUID = store.addCameraCoregistration(source, struct());
            schema = getUMITRigSchema();
            rigFile = fullfile(store.RigRoot, schema.files.rigMetadata);
            loaded = load(rigFile, schema.metadataVariables.rig, '-mat');
            RigInfo = loaded.(schema.metadataVariables.rig);
            index = find(strcmp({RigInfo.resourceRegistry.uuid}, resourceUUID), 1);
            RigInfo.resourceRegistry(index).relativePath = '../outside.mat';
            save(rigFile, 'RigInfo', '-mat');

            report = store.validate('Mode', 'full');
            testCase.verifyFalse(report.isValid);
            testCase.verifyTrue(any(strcmp({report.errors.code}, 'invalidResourcePath')));
        end

        function testMutationRemovesRetiredCalibrationFolder(testCase)
            store = testCase.createRig();
            retiredFolder = fullfile(store.RigRoot, 'calibration-files', 'active');
            mkdir(retiredFolder);
            fid = fopen(fullfile(retiredFolder, 'obsolete.mat'), 'w');
            fclose(fid);

            store.updateRigMetadata(struct('description', 'Trigger v3 persistence.'));
            testCase.verifyFalse(isfolder(fullfile(store.RigRoot, 'calibration-files')));
        end

        function testV2UpgradePreservesCoregistration(testCase)
        %TESTV2UPGRADE... Verify lazy schema migration survives normal writes.

            store = testCase.createRig();
            coregSource = testCase.createTformSource('coreg.mat', 1);
            coregUUID = store.addCameraCoregistration(coregSource, struct());
            original = store.getRigInfo();

            schema = getUMITRigSchema();
            legacyFolder = fullfile(store.RigRoot, 'calibration-files', 'active');
            mkdir(legacyFolder);
            legacyFile = fullfile(legacyFolder, 'LegacyCalibration.mat');
            legacyPayload = 1;
            save(legacyFile, 'legacyPayload', '-mat');

            legacyUUID = lower(char(java.util.UUID.randomUUID()));
            calibrationRecord = original.resourceRegistry(1);
            calibrationRecord.uuid = legacyUUID;
            calibrationRecord.type = 'calibrationFile';
            calibrationRecord.fileName = 'LegacyCalibration.mat';
            calibrationRecord.relativePath = 'calibration-files/active/LegacyCalibration.mat';
            calibrationRecord.status = 'available';
            calibrationRecord.checksum = computeFileChecksum(legacyFile);
            calibrationRecord.sourceFile = legacyFile;

            RigInfo = original;
            RigInfo.schemaVersion = 2;
            RigInfo.isDefault = false;
            RigInfo.activeCalibrationFileUUID = legacyUUID;
            RigInfo.resourceRegistry = [original.resourceRegistry; calibrationRecord];
            rigFile = fullfile(store.RigRoot, schema.files.rigMetadata);
            save(rigFile, 'RigInfo', '-mat');

            reopened = UMITRigStore.open(original.uuid);
            testCase.verifyFalse(reopened.IsReadOnly);
            reopened.updateRigMetadata(struct('description', 'Migrated v2 Rig.'));

            persisted = load(rigFile, schema.metadataVariables.rig, '-mat');
            migrated = persisted.(schema.metadataVariables.rig);
            testCase.verifyEqual(migrated.schemaVersion, 3);
            testCase.verifyEqual(migrated.uuid, original.uuid);
            testCase.verifyEqual(migrated.activeCoregistrationUUID, coregUUID);
            testCase.verifyFalse(isfield(migrated, 'isDefault'));
            testCase.verifyFalse(isfield(migrated, 'activeCalibrationFileUUID'));
            testCase.verifyEqual(numel(migrated.resourceRegistry), 1);
            testCase.verifyEqual(migrated.resourceRegistry.uuid, coregUUID);
            testCase.verifyFalse(isfolder(fullfile(store.RigRoot, 'calibration-files')));

            report = reopened.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, testCase.issueMessages(report.errors));
        end

        function testV1UpgradePreservesCoregistration(testCase)
        %TESTV1UPGRADE... Verify the oldest supported schema upgrades lazily.

            store = testCase.createRig();
            coregSource = testCase.createTformSource('coreg-v1.mat', 1);
            coregUUID = store.addCameraCoregistration(coregSource, struct());
            original = store.getRigInfo();
            schema = getUMITRigSchema();

            RigInfo = original;
            RigInfo.schemaVersion = 1;
            RigInfo = rmfield(RigInfo, ...
                {'status', 'archivedOn', 'metadata', 'cameras', 'illuminations'});
            rigFile = fullfile(store.RigRoot, schema.files.rigMetadata);
            save(rigFile, 'RigInfo', '-mat');

            reopened = UMITRigStore.open(original.uuid);
            testCase.verifyFalse(reopened.IsReadOnly);
            reopened.updateRigMetadata(struct('description', 'Migrated v1 Rig.'));

            persisted = load(rigFile, schema.metadataVariables.rig, '-mat');
            migrated = persisted.(schema.metadataVariables.rig);
            testCase.verifyEqual(migrated.schemaVersion, 3);
            testCase.verifyEqual(migrated.uuid, original.uuid);
            testCase.verifyEqual(migrated.activeCoregistrationUUID, coregUUID);
            testCase.verifyEqual(migrated.status, 'active');
            testCase.verifyTrue(isfield(migrated, 'metadata'));
            testCase.verifyTrue(isfield(migrated, 'cameras'));
            testCase.verifyTrue(isfield(migrated, 'illuminations'));

            report = reopened.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, testCase.issueMessages(report.errors));
        end
    end

    methods (Access = private)
        function store = createConfiguredRigForDuplication(testCase, cameraCount)
        %CREATECONFIGUREDRIGFORDUPLICATION Create one optical Rig fixture.

            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            cameras = struct( ...
                'index', {1, 2}, ...
                'displayName', {'Camera 1', 'Camera 2'}, ...
                'manufacturer', {'Maker A', 'Maker B'}, ...
                'model', {'Model A', 'Model B'}, ...
                'serialNumber', {'Serial A', 'Serial B'}, ...
                'spectrumID', {'PF1024', 'PF1312'});
            cameras = cameras(1:cameraCount);
            illuminations = struct( ...
                'name', {'red', 'green', 'yellow'}, ...
                'displayName', {'Red', 'Green', 'Yellow'}, ...
                'manufacturer', {'', '', ''}, ...
                'model', {'', '', ''}, ...
                'spectrumID', {'LED_632nm', 'LED_521nm', 'LED_593nm'});
            store = UMITRigStore.create(struct( ...
                'rigID', ['DupSource_' suffix(1:8)], ...
                'displayName', 'Duplication Source', ...
                'description', 'Rig duplication fixture.', ...
                'metadata', struct('facility', 'Unit Test', 'room', 7), ...
                'cameras', cameras, ...
                'illuminations', illuminations));
            testCase.RigRoot = store.RigRoot;
        end

        function store = createRig(testCase, variant)
        %CREATERIG Create one non-default rig for a test.

            if nargin < 2
                variant = 1;
            end

            uniqueSuffix = char(java.util.UUID.randomUUID());
            rigID = sprintf('UnitTestRig_%d_%s', variant, uniqueSuffix(1:8));

            rigInfo = struct();
            rigInfo.rigID = rigID;
            rigInfo.displayName = ['Unit Test Rig ' uniqueSuffix(1:8)];
            rigInfo.description = 'Temporary rig created by tests.';

            store = UMITRigStore.create(rigInfo);

            if variant == 1 && strcmp(store.getRigInfo().status, 'available')
                store.activateRig();
            end

            if variant == 2
                testCase.RigRoot2 = store.RigRoot;
            else
                testCase.RigRoot = store.RigRoot;
            end
        end

        function filePath = createTformSource(testCase, fileName, seed)
        %CREATETFORMSOURCE Create a minimal tform/tformInfo .mat payload,
        %mirroring what DataViewer_Coreg2Cams.saveTform actually saves.

            filePath = fullfile(testCase.SourceFolder, fileName);

            tform = affine2d([1 0 0; 0 1 0; seed 0 1]);
            tformInfo = struct( ...
                'OriginalImages', {{single(seed .* ones(4, 5)), single(seed .* ones(4, 5))}}, ...
                'CandidateSource', 'unit-test');

            save(filePath, 'tform', 'tformInfo', '-mat');
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
            ownedUUIDs = strings(0, 1);
            for root = string({testCase.RigRoot, testCase.RigRoot2})
                rigFile = fullfile(root, schema.files.rigMetadata);
                if isfile(rigFile)
                    rigLoaded = load(rigFile, schema.metadataVariables.rig, '-mat');
                    ownedUUIDs(end+1, 1) = string(rigLoaded.RigInfo.uuid); %#ok<AGROW>
                end
            end
            if isfield(defaultInfo, 'rigUUID') && ...
                    any(strcmpi(string(defaultInfo.rigUUID), ownedUUIDs))
                delete(defaultFile);
            end
        end

        function message = issueMessages(~, issues)
        %ISSUEMESSAGES Combine validation errors for test diagnostics.

            if isempty(issues)
                message = '';
            else
                message = strjoin({issues.message}, ' | ');
            end
        end
    end
end
