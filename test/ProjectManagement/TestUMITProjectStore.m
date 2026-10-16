classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestUMITProjectStore < matlab.unittest.TestCase
%TESTUMITPROJECTSTORE Unit tests for centralized UMIT project management.

    properties
        TempRoot
        ProjectRoot
        ProjectUUID
        ProjectRoot2
        ProjectUUID2
        SourceFolder
        RigRoot
        RigRoot2
        OriginalActiveRigUUID
    end

    methods (TestMethodSetup)
        function createTemporaryFolders(testCase)
        %CREATETEMPORARYFOLDERS Create isolated folders for each test.

            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.ProjectRoot = '';
            testCase.ProjectUUID = '';
            testCase.ProjectRoot2 = '';
            testCase.ProjectUUID2 = '';
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
            if lifecycle.activeCount == 0
                rigInfo = struct('rigID', ['ProjectStoreTestRig_' ...
                    strrep(char(java.util.UUID.randomUUID()), '-', '')]);
                rigStore = UMITRigStore.create(rigInfo);
                testCase.RigRoot = rigStore.RigRoot;
            end
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolders(testCase)
        %REMOVETEMPORARYFOLDERS Remove all files created by a test.
        %
        %   Every step below is best-effort and independently guarded: this
        %   the temporary UMIT root (UMITTemporaryRootFixture) is shared by
        %   consecutive test classes, so one locked/failed removal must never
        %   abort the remaining cleanup steps and leave other fixtures behind.

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

            testCase.removeRigFixture(testCase.RigRoot);
            testCase.removeRigFixture(testCase.RigRoot2);

            if ~isempty(testCase.ProjectRoot) && ...
                    isfolder(testCase.ProjectRoot)
                try
                    rmdir(testCase.ProjectRoot, 's');
                catch
                end
            end

            if ~isempty(testCase.ProjectRoot2) && ...
                    isfolder(testCase.ProjectRoot2)
                try
                    rmdir(testCase.ProjectRoot2, 's');
                catch
                end
            end

            if isfolder(testCase.TempRoot)
                try
                    rmdir(testCase.TempRoot, 's');
                catch
                end
            end
        end
    end

    methods (Access = private)
        function removeRigFixture(~, rigRoot)
        %REMOVERIGFIXTURE Best-effort removal of one fixture Rig this test
        %created in the (temporary) Rig store.
        %
        %   If the fixture Rig is still the store's default, repoint the
        %   default to another existing Rig first. Deleting a default Rig's
        %   folder without repointing leaves UMITRigStore.getDefaultRig
        %   permanently broken for every later caller,
        %   not just this test. If no other Rig exists to repoint to, leave
        %   the folder in place rather than corrupt the default pointer.

            if isempty(rigRoot) || ~isfolder(rigRoot)
                return
            end

            try
                rigs = UMITRigStore.listRigs();
                idx = find(strcmp(rigs.RigRoot, rigRoot), 1, 'first');
                if ~isempty(idx) && rigs.IsDefault(idx)
                    otherUUID = rigs.RigUUID(rigs.RigUUID ~= rigs.RigUUID(idx));
                    if isempty(otherUUID)
                        return
                    end
                    UMITRigStore.setDefaultRig(otherUUID(1));
                end
            catch
                % Best effort; fall through to folder removal regardless.
            end

            try
                rmdir(rigRoot, 's');
            catch
            end
        end
    end

    methods (Test)
        function testCreateProject(testCase)
        %TESTCREATEPROJECT Verify canonical top-level project creation.

            store = testCase.createProject();
            schema = getUMITProjectSchema();
            projectsRoot = UMITProjectStore.getProjectsRoot();

            testCase.verifyEqual( ...
                fileparts(testCase.ProjectRoot), projectsRoot);
            testCase.verifyTrue(endsWith( ...
                testCase.ProjectRoot, testCase.ProjectUUID));

            testCase.verifyTrue(isfile(fullfile( ...
                testCase.ProjectRoot, schema.files.projectMetadata)));
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, schema.folders.subjects)));
            testCase.verifyFalse(store.IsReadOnly);

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testResolveAndOpenProjectByUUID(testCase)
        %TESTRESOLVEANDOPENPROJECTBYUUID Verify static-root discovery.

            store = testCase.createProject();
            ProjectInfo = store.getProjectInfo();
            clear store

            resolvedRoot = UMITProjectStore.resolveProjectRoot( ...
                ProjectInfo.projectUUID);
            testCase.verifyEqual(resolvedRoot, testCase.ProjectRoot);

            reopened = UMITProjectStore.open(ProjectInfo.projectUUID);
            reopenedInfo = reopened.getProjectInfo();
            testCase.verifyEqual( ...
                reopenedInfo.projectUUID, ProjectInfo.projectUUID);
            testCase.verifyEqual( ...
                reopened.ProjectRoot, testCase.ProjectRoot);
        end

        function testListProjectsIncludesCreatedProject(testCase)
        %TESTLISTPROJECTSINCLUDESCAPTUREDPROJECT Verify project enumeration.

            store = testCase.createProject();
            ProjectInfo = store.getProjectInfo();

            projects = UMITProjectStore.listProjects();
            idx = find(strcmp( ...
                projects.ProjectUUID, ProjectInfo.projectUUID), ...
                1, 'first');

            testCase.verifyNotEmpty(idx);
            testCase.verifyTrue(projects.IsReadable(idx));
            testCase.verifyEqual( ...
                char(projects.ProjectRoot(idx)), ...
                testCase.ProjectRoot);
        end

        function testInvalidIDsAreReportedNotRewritten(testCase)
        %TESTINVALIDIDSAREREPORTEDNOTREWRITTEN Audit item 24: a new or
        %renamed ID that is not a valid name is refused with the valid form
        %in the message, never rewritten silently. Lookups still tolerate
        %unnormalized input.

            store = testCase.createProject();
            message = testCase.errorMessage(@() store.addSubject( ...
                struct('subjectID', 'Mouse 104')), ...
                'Umitoolbox:UMITProjectStore:invalidID');
            testCase.verifySubstring(message, 'Use "Mouse104" instead');
            testCase.verifyEmpty(store.getProjectInfo().subjectRegistry);

            subjectUUID = store.addSubject(struct('subjectID', 'Mouse104'));
            testCase.verifyEqual(store.getSubjectInfo('Mouse 104').uuid, subjectUUID, ...
                'lookups normalize their input');

            message = testCase.errorMessage(@() testCase.addSession(store, 'Mouse104', ...
                struct('sessionID', '2026-07-10 Session-01')), ...
                'Umitoolbox:UMITProjectStore:invalidID');
            testCase.verifySubstring(message, 'x2026_07_10Session_01');
            testCase.verifyError(@() store.renameSubjectID('Mouse104', 'Mouse 105'), ...
                'Umitoolbox:UMITProjectStore:invalidID');
            testCase.verifyError(@() store.addSubject( ...
                struct('subjectID', 'CON')), ...
                'Umitoolbox:UMITProjectStore:invalidID');
            testCase.verifyError(@() store.addSubject( ...
                struct('subjectID', 'Mouse105', 'colour', 'brown')), ...
                'Umitoolbox:UMITProjectStore:unsupportedField');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testCreateAndRenameHierarchy(testCase)
        %TESTCREATEANDRENAMEHIERARCHY Verify subject/session folder renames.

            store = testCase.createProject();
            store.addSubject(struct( ...
                'subjectID', 'Mouse_104', ...
                'displayName', 'Mouse 104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', '2026-07-10_S01'));

            store.renameSessionID( ...
                'Mouse_104', '2026-07-10_S01', '2026-07-10_S02');
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', '2026-07-10_S02')));
            testCase.verifyFalse(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', '2026-07-10_S01')));

            store.renameSubjectID('Mouse_104', 'Mouse_104A');
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104A')));
            testCase.verifyFalse(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104')));

            SessionInfo = store.getSessionInfo( ...
                'Mouse_104A', '2026-07-10_S02');
            testCase.verifyEqual(SessionInfo.subjectID, 'Mouse_104A');
            testCase.verifyEqual(SessionInfo.sessionID, '2026-07-10_S02');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testArchiveRestoreAndReplacement(testCase)
        %TESTARCHIVERESTOREANDREPLACEMENT Verify managed resource states.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            source1 = testCase.createMATSource('reference1.mat', 1);
            source2 = testCase.createMATSource('reference2.mat', 2);

            uuid1 = store.addImageReference('Mouse_104', source1, ...
                struct('displayName', 'Reference 1'));
            uuid2 = store.addImageReference('Mouse_104', source2, ...
                struct('displayName', 'Reference 2'));

            testCase.verifyError(@() store.archiveResource(uuid1), ...
                'Umitoolbox:UMITProjectStore:archiveFailed');

            store.archiveResource(uuid1, 'ReplacementUUID', uuid2);
            SubjectInfo = store.getSubjectInfo('Mouse_104');

            idx1 = find(strcmp({SubjectInfo.resourceRegistry.uuid}, uuid1), 1);
            idx2 = find(strcmp({SubjectInfo.resourceRegistry.uuid}, uuid2), 1);
            testCase.verifyEqual( ...
                SubjectInfo.resourceRegistry(idx1).status, 'archived');
            testCase.verifyEqual( ...
                SubjectInfo.resourceRegistry(idx2).status, 'active');
            testCase.verifyEqual(SubjectInfo.activeImageReferenceUUID, uuid2);
            testCase.verifyTrue(contains( ...
                SubjectInfo.resourceRegistry(idx1).relativePath, '/archive/'));

            store.restoreResource(uuid1);
            SubjectInfo = store.getSubjectInfo('Mouse_104');
            idx1 = find(strcmp({SubjectInfo.resourceRegistry.uuid}, uuid1), 1);
            testCase.verifyEqual( ...
                SubjectInfo.resourceRegistry(idx1).status, 'available');
            testCase.verifyTrue(contains( ...
                SubjectInfo.resourceRegistry(idx1).relativePath, '/active/'));
            testCase.verifyEqual(SubjectInfo.activeImageReferenceUUID, uuid2);

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testSubjectRenameUpdatesResourcePaths(testCase)
        %TESTSUBJECTRENAMEUPDATESRESOURCEPATHS Verify path propagation.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            source = testCase.createMATSource('reference.mat', 1);
            resourceUUID = store.addImageReference( ...
                'Mouse_104', source, struct());

            store.renameSubjectID('Mouse_104', 'Mouse_105');
            SubjectInfo = store.getSubjectInfo('Mouse_105');
            idx = find(strcmp( ...
                {SubjectInfo.resourceRegistry.uuid}, resourceUUID), 1);

            testCase.verifyTrue(startsWith( ...
                SubjectInfo.resourceRegistry(idx).relativePath, ...
                'subjects/Mouse_105/'));
            testCase.verifyTrue(isfile(testCase.resolveRelative( ...
                SubjectInfo.resourceRegistry(idx).relativePath)));

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testSessionsKeepNoRig(testCase)
        %TESTSESSIONSKEEPNORIG Project redesign 1.12 (T1-T5): a session has
        %no rig fields; the folder's rig lives only in its rig ledger, which
        %addSession neither reads into the session nor writes. Rig fields
        %are refused like any unknown field.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            boundRig = UMITRigStore.create(struct('rigID', ...
                ['NoMirrorRig_' strrep(char(java.util.UUID.randomUUID()), '-', '')]));
            testCase.RigRoot2 = boundRig.RigRoot;
            saveFolder = fullfile(testCase.TempRoot, 'RigBoundSaveFolder');
            testCase.writeValidAcqInfos(saveFolder);
            testCase.bindFolderRig(saveFolder, boundRig);
            ledgerBefore = fileread(UMITRigLedger.filePath(saveFolder));
            acqBefore = load(fullfile(saveFolder, 'AcqInfos.mat'));

            store.addSession('Mouse_104', struct( ...
                'sessionID', 'Session_01', 'processedDataFolder', saveFolder));
            SessionInfo = store.getSessionInfo('Mouse_104', 'Session_01');
            testCase.verifyFalse(any(isfield(SessionInfo, {'rigUUID', 'rigID', ...
                'activeRegistrationTransformUUID'})));
            testCase.verifyEqual(fileread(UMITRigLedger.filePath(saveFolder)), ledgerBefore);
            testCase.verifyEqual(load(fullfile(saveFolder, 'AcqInfos.mat')), acqBefore);
            context = store.findSessionByDataFolder(saveFolder);
            testCase.verifyFalse(any(isfield(context, {'rigUUID', 'rigID'})));
            testCase.verifyFalse(ismethod(store, 'migrateSessionRigAssignment'));

            otherFolder = fullfile(testCase.TempRoot, 'OtherSaveFolder');
            testCase.writeValidAcqInfos(otherFolder);
            testCase.verifyError(@() store.addSession('Mouse_104', struct( ...
                'sessionID', 'Session_02', 'rigID', boundRig.getRigInfo().rigID, ...
                'processedDataFolder', otherFolder)), ...
                'Umitoolbox:UMITProjectStore:unsupportedField');
            testCase.verifyError(@() store.updateSessionMetadata( ...
                'Mouse_104', 'Session_01', struct('notes', 'x', ...
                'rigUUID', boundRig.getRigInfo().uuid)), ...
                'Umitoolbox:UMITProjectStore:updateSessionFailed');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, testCase.issueMessages(report.errors));
        end

        function testChecksumCorruptionIsDetected(testCase)
        %TESTCHECKSUMCORRUPTIONISDETECTED Detect manual resource modification.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            source = testCase.createMATSource('reference.mat', 1);
            resourceUUID = store.addImageReference( ...
                'Mouse_104', source, struct());

            SubjectInfo = store.getSubjectInfo('Mouse_104');
            idx = find(strcmp( ...
                {SubjectInfo.resourceRegistry.uuid}, resourceUUID), 1);
            resourcePath = testCase.resolveRelative( ...
                SubjectInfo.resourceRegistry(idx).relativePath);

            tampered = 42;
            save(resourcePath, 'tampered', '-append');

            quickReport = store.validate('Mode', 'quick');
            fullReport = store.validate('Mode', 'full');

            testCase.verifyTrue(quickReport.isValid);
            testCase.verifyFalse(fullReport.isValid);
            testCase.verifyTrue(any(strcmp( ...
                {fullReport.errors.code}, 'resource_checksum_mismatch')));
        end

        function testBrokenSessionBlocksOnlyItself(testCase)
        %TESTBROKENSESSIONBLOCKSONLYITSELF Brief 1.11: errors in one session
        %block mutations of that session only; the project stays writable.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_105'));
            brokenFolder = fullfile(testCase.TempRoot, 'BrokenFolder');
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Broken', 'processedDataFolder', brokenFolder));
            testCase.addSession(store, 'Mouse_104', struct('sessionID', 'Healthy'));
            delete(fullfile(brokenFolder, 'AcqInfos.mat'));

            reopened = UMITProjectStore.open(testCase.ProjectUUID);
            report = reopened.LastValidationReport;
            testCase.verifyFalse(report.isValid);
            testCase.verifyFalse(reopened.IsReadOnly);
            broken = report.errors(strcmp({report.errors.code}, 'not_a_dataset'));
            testCase.verifyNotEmpty(broken);
            testCase.verifyFalse(any([broken.storeLevel]));

            reopened.updateSessionMetadata('Mouse_104', 'Healthy', struct('notes', 'ok'));
            reopened.updateSubjectMetadata('Mouse_104', struct('notes', 'ok'));
            reopened.addSubject(struct('subjectID', 'Mouse_106'));
            reopened.updateProjectMetadata(struct('notes', 'ok'));
            testCase.verifyError(@() reopened.updateSessionMetadata( ...
                'Mouse_104', 'Broken', struct('notes', 'blocked')), ...
                'Umitoolbox:UMITProjectStore:invalidTarget');
            testCase.verifyFalse(reopened.IsReadOnly);

            % Narrow repair: removing the broken session is allowed.
            reopened.removeSessionFromProject('Mouse_104', 'Broken');
            testCase.verifyTrue(reopened.validate('Mode', 'full').isValid);
        end

        function testMutationsValidateOnlyWhatTheyTouch(testCase)
        %TESTMUTATIONSVALIDATEONLYWHATTHEYTOUCH Brief 1.11: mutations run a
        %quick validation scoped to their subject/session, so a damaged
        %resource elsewhere (found only by a full checksum pass) does not
        %block them.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_105'));
            resourceUUID = store.addImageReference('Mouse_104', ...
                testCase.createMATSource('reference.mat', 1), struct());
            resource = store.getResource(resourceUUID);
            tampered = 42;
            save(resource.absolutePath, 'tampered', '-append');
            testCase.verifyFalse(store.validate('Mode', 'full').isValid);

            store.updateSubjectMetadata('Mouse_105', struct('notes', 'ok'));
            store.updateSubjectMetadata('Mouse_104', struct('notes', 'ok'));
            testCase.verifyEqual(store.getSubjectInfo('Mouse_105').notes, 'ok');
            testCase.verifyFalse(store.IsReadOnly);
        end

        function testHeaderedDatIsADatasetWithoutAcqInfos(testCase)
        %TESTHEADEREDDATISADATASETWITHOUTACQINFOS Audit 12.1 item 5: the
        %"is a dataset" test is format-neutral (isUmitDataset).

            saveFolder = fullfile(testCase.TempRoot, 'HeaderOnly');
            mkdir(saveFolder);
            writeTestDat(fullfile(saveFolder, 'green.dat'), ones(4, 5, 3, 'single'), 10, 5);
            testCase.verifyTrue(isUmitDataset(saveFolder));
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSession('Mouse_104', struct( ...
                'sessionID', 'Session_01', 'processedDataFolder', saveFolder));
            report = store.validate('Mode', 'quick');
            testCase.verifyTrue(report.isValid, testCase.issueMessages(report.errors));
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
        end

        function testUnregisteredFolderOpensReadOnly(testCase)
        %TESTUNREGISTEREDFOLDEROPENSREADONLY Detect direct folder edits.

            store = testCase.createProject(); %#ok<NASGU>
            clear store

            mkdir(fullfile(testCase.ProjectRoot, ...
                'subjects', 'ManuallyCreatedMouse'));

            reopened = UMITProjectStore.open(testCase.ProjectUUID);
            testCase.verifyTrue(reopened.IsReadOnly);
            testCase.verifyFalse(reopened.LastValidationReport.isValid);
            errors = reopened.LastValidationReport.errors;
            unregistered = errors(strcmp({errors.code}, 'unregistered_folder'));
            testCase.verifyNotEmpty(unregistered);
            testCase.verifyTrue(all([unregistered.storeLevel]), ...
                'a folder that belongs to no registered subject is store-level');
        end

        function testDeleteInvalidReadOnlyProject(testCase)
        %TESTDELETEINVALIDREADONLYPROJECT Allow explicit cleanup of bad metadata.

            store = testCase.createProject(); %#ok<NASGU>
            clear store
            mkdir(fullfile(testCase.ProjectRoot, ...
                'subjects', 'ManuallyCreatedMouse'));

            reopened = UMITProjectStore.open(testCase.ProjectUUID);
            testCase.verifyTrue(reopened.IsReadOnly);

            result = reopened.deleteProject(testCase.ProjectUUID);

            testCase.verifyEqual(result.projectUUID, testCase.ProjectUUID);
            testCase.verifyFalse(isfolder(testCase.ProjectRoot));
        end

        function testCaseInsensitiveDuplicateIDIsRejected(testCase)
        %TESTCASEINSENSITIVEDUPLICATEIDISREJECTED Enforce Windows-safe IDs.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            testCase.verifyError(@() store.addSubject( ...
                struct('subjectID', 'mouse_104')), ...
                'Umitoolbox:UMITProjectStore:duplicateID');
        end

        function testCaseOnlyRenameUsesStaging(testCase)
        %TESTCASEONLYRENAMEUSESSTAGING Verify case-only ID changes.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.renameSubjectID('Mouse_104', 'mouse_104');

            SubjectInfo = store.getSubjectInfo('mouse_104');
            testCase.verifyEqual(SubjectInfo.subjectID, 'mouse_104');
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'mouse_104')));

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testRejectDataFolderInsideProject(testCase)
        %TESTREJECTDATAFOLDERINSIDEPROJECT Keep imaging paths external.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', '2026-07-10_S01'));

            internalDataPath = fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104');

            testCase.verifyError(@() store.bindRawDataFolder( ...
                'Mouse_104', '2026-07-10_S01', internalDataPath), ...
                'Umitoolbox:UMITProjectStore:internalDataPath');
        end

        function testClearStaleLockArchivesLockFile(testCase)
        %TESTCLEARSTALELOCKARCHIVESLOCKFILE Verify managed lock recovery.

            store = testCase.createProject();
            schema = getUMITProjectSchema();
            lockFolder = fullfile(testCase.ProjectRoot, ...
                schema.folders.internal, schema.folders.lock);
            mkdir(lockFolder);
            lockPath = fullfile(lockFolder, schema.files.lockMetadata);

            LockInfo = struct();
            LockInfo.operation = 'interruptedTest';
            LockInfo.createdOn = datetime('now') - hours(1);
            LockInfo.processID = -1;
            LockInfo.userName = 'test';
            LockInfo.hostName = 'test';
            saveMatAtomic(lockPath, 'LockInfo', LockInfo);

            store.clearStaleLock('MinimumAgeMinutes', 10);
            testCase.verifyFalse(isfolder(lockFolder));

            archivedLocks = dir(fullfile(testCase.ProjectRoot, ...
                schema.folders.internal, schema.folders.logs, ...
                'stale_lock_*'));
            testCase.verifyEqual(numel(archivedLocks), 1);
        end

        function testListSubjectResourcesAndArchiveFilters(testCase)
        %TESTLISTSUBJECTRESOURCESANDARCHIVEFILTERS Verify query filtering.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            source1 = testCase.createMATSource('reference1.mat', 1);
            source2 = testCase.createMATSource('reference2.mat', 2);
            uuid1 = store.addImageReference( ...
                'Mouse_104', source1, struct());
            uuid2 = store.addImageReference( ...
                'Mouse_104', source2, struct());
            store.archiveResource(uuid2);

            resources = store.listSubjectResources( ...
                'Mouse_104', 'Type', 'imageReference');
            testCase.verifyEqual(numel(resources), 1);
            testCase.verifyEqual(resources.uuid, uuid1);
            testCase.verifyEqual(resources.ownerType, 'subject');
            testCase.verifyTrue(resources.fileExists);
            testCase.verifyTrue(isfile(resources.absolutePath));

            archived = store.listSubjectResources( ...
                'Mouse_104', ...
                'Type', 'imageReference', ...
                'Status', 'archived', ...
                'VerifyFiles', true);
            testCase.verifyEqual(numel(archived), 1);
            testCase.verifyEqual(archived.uuid, uuid2);
            testCase.verifyEqual(archived.status, 'archived');

            allResources = store.listSubjectResources( ...
                'Mouse_104', 'Status', {});
            testCase.verifyEqual(sort({allResources.uuid}), ...
                sort({uuid1, uuid2}));
        end

        function testRetiredTransformResourcesAreTolerated(testCase)
        %TESTRETIREDTRANSFORMRESOURCESARETOLERATED Section 4 Q5/Q9: the
        %registration-transform API is gone; records written before stay
        %valid (a warning), are never listed, and new sessions have no
        %transforms folder.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', ...
                struct('sessionID', '2026-07-10_S01'));
            sessionRel = 'subjects/Mouse_104/sessions/2026-07-10_S01';
            sessionPath = testCase.resolveRelative(sessionRel);
            testCase.verifyFalse(isfolder(fullfile(sessionPath, 'transforms')));
            for name = {'addRegistrationTransform', 'setActiveRegistrationTransform', ...
                    'getActiveRegistrationTransform'}
                testCase.verifyFalse(ismethod(store, name{1}), name{1});
            end

            % A legacy session written before the retirement.
            legacyRel = [sessionRel '/transforms/registration/active/SessionToReference_1.mat'];
            mkdir(fileparts(testCase.resolveRelative(legacyRel)));
            payload = 1;
            save(testCase.resolveRelative(legacyRel), 'payload');
            sessionMat = fullfile(sessionPath, 'session.mat');
            loaded = load(sessionMat, 'SessionInfo');
            SessionInfo = loaded.SessionInfo;
            legacyUUID = lower(char(java.util.UUID.randomUUID()));
            SessionInfo.activeRegistrationTransformUUID = legacyUUID;
            SessionInfo.resourceRegistry = struct('uuid', legacyUUID, ...
                'type', 'registrationTransform', 'fileName', 'SessionToReference_1.mat', ...
                'relativePath', legacyRel, 'displayName', 'legacy', 'description', '', ...
                'createdOn', datetime('now'), 'modifiedOn', datetime('now'), ...
                'status', 'active', 'archivedOn', NaT, 'checksum', '', 'sourceFile', '');
            save(sessionMat, 'SessionInfo');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, testCase.issueMessages(report.errors));
            testCase.verifyTrue(any(strcmp({report.warnings.code}, 'retired_resource_type')));
            testCase.verifyEmpty(store.listSessionResources('Mouse_104', '2026-07-10_S01'));
            store.updateSessionMetadata('Mouse_104', '2026-07-10_S01', struct('notes', 'still writable'));
            testCase.verifyEqual(store.getSessionInfoByUUID(sessionUUID).notes, 'still writable');
        end

        function testRigResourcesLiveInUMITRigStoreNotProject(testCase)
        %TESTRIGRESOURCESLIVEINUMITRIGSTORENOTPROJECT Verify rig resources
        %(cameraCoregistration) are owned by the
        %independent UMITRigStore, not by UMITProjectStore, and that a
        %session keeps no pointer to it (the folder's rig ledger does).

            rigStore = UMITRigStore.create(struct('rigID', 'WidefieldRig01'));
            testCase.RigRoot = rigStore.RigRoot;

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', '2026-07-10_S01', ...
                'rigID', 'WidefieldRig01'));

            coregSource = testCase.createMATSource('coregistration.mat', 2);
            coregUUID = rigStore.addCameraCoregistration(coregSource, struct());

            rigResources = rigStore.listResources();
            testCase.verifyEqual(numel(rigResources), 1);
            testCase.verifyEqual(rigResources.uuid, coregUUID);

            coregResources = rigStore.listResources('Type', 'cameraCoregistration');
            testCase.verifyEqual(numel(coregResources), 1);
            testCase.verifyEqual(coregResources.uuid, coregUUID);

            SessionInfo = store.getSessionInfo('Mouse_104', '2026-07-10_S01');
            testCase.verifyFalse(any(isfield(SessionInfo, {'rigID', 'rigUUID'})), ...
                'Sessions keep no rig (project redesign 1.12).');
        end

        function testGetResourceAndResolvePath(testCase)
        %TESTGETRESOURCEANDRESOLVEPATH Verify UUID-based resource lookup.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            source = testCase.createMATSource('reference.mat', 1);
            resourceUUID = store.addImageReference( ...
                'Mouse_104', source, struct());

            SubjectInfo = store.getSubjectInfo('Mouse_104');
            resource = store.getResource(resourceUUID);

            testCase.verifyEqual(resource.uuid, resourceUUID);
            testCase.verifyEqual(resource.ownerType, 'subject');
            testCase.verifyEqual(resource.ownerUUID, SubjectInfo.uuid);
            testCase.verifyTrue(resource.fileExists);
            testCase.verifyEqual( ...
                store.resolveResourcePath(resourceUUID), ...
                resource.absolutePath);

            testCase.verifyError(@() store.getResource( ...
                '00000000-0000-0000-0000-000000000000'), ...
                'Umitoolbox:UMITProjectStore:resourceNotFound');
        end


        function testGetResourceRejectsDuplicateUUID(testCase)
        %TESTGETRESOURCEREJECTSDUPLICATEUUID Reject ambiguous UUID lookup.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            uuid1 = store.addImageReference( ...
                'Mouse_104', ...
                testCase.createMATSource('reference1.mat', 1), struct());
            store.addImageReference( ...
                'Mouse_104', ...
                testCase.createMATSource('reference2.mat', 2), struct());

            schema = getUMITProjectSchema();
            subjectPath = fullfile(testCase.ProjectRoot, ...
                schema.folders.subjects, 'Mouse_104');
            SubjectInfo = store.getSubjectInfo('Mouse_104');
            SubjectInfo.resourceRegistry(2).uuid = uuid1;
            saveMatAtomic( ...
                fullfile(subjectPath, schema.files.subjectMetadata), ...
                schema.metadataVariables.subject, SubjectInfo);

            testCase.verifyError(@() store.getResource(uuid1), ...
                'Umitoolbox:UMITProjectStore:duplicateResourceUUID');
        end

        function testMissingResourceFileQueryBehavior(testCase)
        %TESTMISSINGRESOURCEFILEQUERYBEHAVIOR Verify optional existence checks.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            source = testCase.createMATSource('reference.mat', 1);
            resourceUUID = store.addImageReference( ...
                'Mouse_104', source, struct());

            resource = store.getResource(resourceUUID);
            delete(resource.absolutePath);

            staleResource = store.getResource(resourceUUID);
            testCase.verifyFalse(staleResource.fileExists);
            testCase.verifyEqual( ...
                store.resolveResourcePath( ...
                resourceUUID, 'RequireFile', false), ...
                staleResource.absolutePath);

            testCase.verifyError(@() ...
                store.resolveResourcePath(resourceUUID), ...
                'Umitoolbox:UMITProjectStore:missingResourceFile');
            testCase.verifyError(@() ...
                store.listSubjectResources( ...
                'Mouse_104', 'VerifyFiles', true), ...
                'Umitoolbox:UMITProjectStore:missingResourceFile');
        end

        function testActiveResourceGetters(testCase)
        %TESTACTIVERESOURCEGETTERS Verify typed active-resource queries
        %(project-owned resource types only; rig-owned types are covered by
        %testRigActiveResourceGetters against UMITRigStore directly).

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct('sessionID', '2026-07-10_S01'));

            testCase.verifyEmpty( ...
                store.getActiveImageReference('Mouse_104'));

            referenceUUID = store.addImageReference( ...
                'Mouse_104', ...
                testCase.createMATSource('reference.mat', 1), struct());

            activeReference = ...
                store.getActiveImageReference('Mouse_104');

            testCase.verifyEqual(activeReference.uuid, referenceUUID);
        end

        function testRigActiveResourceGetters(testCase)
        %TESTRIGACTIVERESOURCEGETTERS Verify typed active-resource queries
        %against the independent UMITRigStore.

            rigStore = UMITRigStore.create(struct('rigID', 'WidefieldRig01'));
            testCase.RigRoot = rigStore.RigRoot;

            testCase.verifyEmpty(rigStore.getActiveCameraCoregistration());

            coregUUID = rigStore.addCameraCoregistration( ...
                testCase.createMATSource('coregistration.mat', 1), struct());

            activeCoregistration = rigStore.getActiveCameraCoregistration();

            testCase.verifyEqual(activeCoregistration.uuid, coregUUID);
        end

        function testRigActiveGetterRejectsWrongResourceType(testCase)
        %TESTRIGACTIVEGETTERREJECTSWRONGRESOURCETYPE Detect pointer corruption
        %on the independent UMITRigStore.

            rigStore = UMITRigStore.create(struct('rigID', 'WidefieldRig01'));
            testCase.RigRoot = rigStore.RigRoot;

            rigStore.addCameraCoregistration( ...
                testCase.createMATSource('coregistration.mat', 1), struct());
            schema = getUMITRigSchema();
            RigInfo = rigStore.getRigInfo();
            RigInfo.activeCoregistrationUUID = lower(char(java.util.UUID.randomUUID()));
            saveMatAtomic( ...
                fullfile(rigStore.RigRoot, schema.files.rigMetadata), ...
                schema.metadataVariables.rig, RigInfo);

            testCase.verifyError(@() ...
                rigStore.getActiveCameraCoregistration(), ...
                'Umitoolbox:UMITRigStore:invalidActiveResource');
        end

        function testFindSessionByRawAndProcessedFolder(testCase)
        %TESTFINDSESSIONBYRAWANDPROCESSEDFOLDER Resolve project-ledger context.

            rawFolder = fullfile(testCase.TempRoot, 'RawData');
            processedFolder = fullfile(testCase.TempRoot, 'ProcessedData');
            unknownFolder = fullfile(testCase.TempRoot, 'UnknownData');
            mkdir(rawFolder);
            mkdir(processedFolder);
            mkdir(unknownFolder);
            testCase.writeRawDataFile(rawFolder);
            testCase.writeValidAcqInfos(processedFolder);

            rigStore = UMITRigStore.create(struct('rigID', 'WidefieldRig01'));
            testCase.RigRoot = rigStore.RigRoot;

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', '2026-07-10_S01', ...
                'rigID', 'WidefieldRig01'));

            store.bindRawDataFolder( ...
                'Mouse_104', '2026-07-10_S01', rawFolder);
            store.bindProcessedDataFolder( ...
                'Mouse_104', '2026-07-10_S01', processedFolder);

            rawContext = store.findSessionByDataFolder(rawFolder);
            testCase.verifyEqual(rawContext.subjectID, 'Mouse_104');
            testCase.verifyEqual(rawContext.sessionID, '2026-07-10_S01');
            testCase.verifyFalse(isfield(rawContext, 'rigID'));
            testCase.verifyEqual(rawContext.matchedField, 'rawDataFolder');

            processedContext = ...
                store.findSessionByDataFolder(processedFolder);
            testCase.verifyEqual( ...
                processedContext.matchedField, 'processedDataFolder');
            testCase.verifyEqual( ...
                processedContext.sessionUUID, rawContext.sessionUUID);

            testCase.verifyEmpty( ...
                store.findSessionByDataFolder(unknownFolder));
        end

        function testProjectLedgersOfRawAndSaveFolders(testCase)
        %TESTPROJECTLEDGERSOFRAWANDSAVEFOLDERS Brief section 1.9 / Phase 3.
        %   A RawFolder distinct from the SaveFolder gets its own ledger;
        %   resolving never writes; an identity that does not match the
        %   store is detected.

            rawFolder = fullfile(testCase.TempRoot, 'RawData');
            saveFolder = fullfile(testCase.TempRoot, 'SaveData');
            mkdir(rawFolder);
            mkdir(saveFolder);
            testCase.writeRawDataFile(rawFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct('sessionID', 'S01'));
            store.bindProcessedDataFolder('Mouse_104', 'S01', saveFolder);
            store.bindRawDataFolder('Mouse_104', 'S01', rawFolder);

            rawLedger = UMITProjectLedger.read(rawFolder);
            saveLedger = UMITProjectLedger.read(saveFolder);
            testCase.verifyEqual(rawLedger.binding.folderRole, 'rawDataFolder');
            testCase.verifyEqual(saveLedger.binding.folderRole, 'processedDataFolder');
            testCase.verifyEqual(rawLedger.binding.sessionUUID, saveLedger.binding.sessionUUID);
            testCase.verifyNotEqual(rawLedger.binding.bindingUUID, saveLedger.binding.bindingUUID);

            path = UMITProjectLedger.filePath(saveFolder);
            before = fileread(path);
            UMITProjectStore.resolveProjectBinding(saveFolder);
            store.findSessionByDataFolder(saveFolder);
            store.validate('Mode', 'full');
            testCase.verifyEqual(fileread(path), before, 'reading never writes');

            % A structurally valid ledger whose session is not this one.
            otherSession = lower(char(java.util.UUID.randomUUID()));
            testCase.writeText(path, strrep(before, saveLedger.binding.sessionUUID, otherSession));
            testCase.verifyTrue(UMITProjectLedger.validate(saveFolder).isValid);
            testCase.verifyError(@() UMITProjectStore.resolveProjectBinding(saveFolder), ...
                ?MException);
            report = store.validate('Mode', 'quick');
            testCase.verifySubstring(strjoin({report.errors.code}, ' '), 'invalid_project_binding');
            testCase.writeText(path, before);

            % Unbinding the raw folder keeps its ledger and history.
            store.unbindRawDataFolder('Mouse_104', 'S01');
            rawLedger = UMITProjectLedger.read(rawFolder);
            testCase.verifyEqual(rawLedger.state, 'Unbound');
            testCase.verifyEqual(rawLedger.events{end}.type, 'unbind');
            testCase.verifyEqual(rawLedger.events{end}.reason, 'folderUnbound');
        end

        function testRejectBindingFolderAlreadyAssigned(testCase)
        %TESTREJECTBINDINGFOLDERALREADYASSIGNED Prevent duplicate ownership.

            sharedFolder = fullfile(testCase.TempRoot, 'SharedData');
            mkdir(sharedFolder);
            testCase.writeRawDataFile(sharedFolder);
            testCase.writeValidAcqInfos(sharedFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_02'));

            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', sharedFolder);

            testCase.verifyError(@() ...
                store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_02', sharedFolder), ...
                'Umitoolbox:UMITProjectStore:dataFolderAlreadyBound');
        end

        function testFindSessionPathComparisonIsCaseInsensitiveOnWindows(testCase)
        %TESTFINDSESSIONPATHCOMPARISONISCASEINSENSITIVEONWINDOWS Verify Windows paths.

            testCase.assumeTrue(ispc, ...
                'Case-insensitive path test applies to Windows.');

            rawFolder = fullfile(testCase.TempRoot, 'MixedCaseRawData');
            mkdir(rawFolder);
            testCase.writeRawDataFile(rawFolder, '.TIF');

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', rawFolder);

            context = store.findSessionByDataFolder(lower(rawFolder));
            testCase.verifyEqual(context.sessionID, 'Session_01');
        end

        function testProjectBindingFileAndSessionReciprocity(testCase)
        %TESTPROJECTBINDINGFILEANDSESSIONRECIPROCITY Verify project-ledger binding.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            ProjectBinding = store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            % The binding is the folder's project ledger (brief section 1.9).
            ledger = UMITProjectLedger.read(saveFolder);
            testCase.verifyEqual(ledger.state, 'Bound');
            testCase.verifyEqual(ledger.binding.bindingUUID, ProjectBinding.bindingUUID);
            testCase.verifyEqual(ledger.binding.folderRole, 'processedDataFolder');
            ProjectInfo = store.getProjectInfo();
            testCase.verifyEqual(ledger.binding.projectUUID, ProjectInfo.projectUUID);
            testCase.verifyEqual(ledger.events{end}.type, 'bind');
            testCase.verifyEqual(ledger.snapshot.subjectID, 'Mouse_104');
            testCase.verifyEqual(ledger.snapshot.sessionID, 'Session_01');
            testCase.verifyEqual(ledger.snapshot.project.projectName, ProjectInfo.projectName);

            SessionInfo = store.getSessionInfo( ...
                'Mouse_104', 'Session_01');
            testCase.verifyEqual(SessionInfo.processedDataFolder, ...
                saveFolder);
            testCase.verifyEqual( ...
                SessionInfo.processedDataBindingUUID, ...
                ProjectBinding.bindingUUID);

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testResolveBindingAndUUIDLookups(testCase)
        %TESTRESOLVEBINDINGANDUUIDLOOKUPS Resolve all hierarchy levels by UUID.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct( ...
                'subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            ProjectInfo = store.getProjectInfo();
            [ProjectInfoByUUID, projectRoot] = ...
                UMITProjectStore.getProjectInfoByUUID( ...
                ProjectInfo.projectUUID);
            testCase.verifyEqual(ProjectInfoByUUID.projectUUID, ...
                ProjectInfo.projectUUID);
            testCase.verifyEqual(projectRoot, testCase.ProjectRoot);

            SubjectInfo = store.getSubjectInfoByUUID(subjectUUID);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SubjectInfo.subjectID, 'Mouse_104');
            testCase.verifyEqual(SessionInfo.sessionID, 'Session_01');

            uuidContext = store.getSessionContextByUUID(sessionUUID);
            testCase.verifyEqual(uuidContext.subjectUUID, subjectUUID);
            testCase.verifyEqual(uuidContext.sessionUUID, sessionUUID);

            [bindingContext, resolvedStore] = ...
                UMITProjectStore.resolveProjectBinding(saveFolder);
            testCase.verifyEqual(bindingContext.projectUUID, ...
                ProjectInfo.projectUUID);
            testCase.verifyEqual(bindingContext.subjectUUID, subjectUUID);
            testCase.verifyEqual(bindingContext.sessionUUID, sessionUUID);
            testCase.verifyEqual(resolvedStore.ProjectRoot, ...
                testCase.ProjectRoot);
        end

        function testUnbindProcessedDataFolderIsDeprecated(testCase)
        %TESTUNBINDPROCESSEDDATAFOLDERISDEPRECATED Preserve the invariant.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            testCase.verifyError(@() ...
                store.unbindProcessedDataFolder( ...
                'Mouse_104', 'Session_01'), ...
                ['Umitoolbox:UMITProjectStore:' ...
                 'unbindSaveFolderDeprecated']);

            SessionInfo = store.getSessionInfo( ...
                'Mouse_104', 'Session_01');
            testCase.verifyEqual( ...
                SessionInfo.processedDataFolder, saveFolder);
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(saveFolder), 'bound'));

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testMissingSaveFolderPreservesSessionAndStoredPath(testCase)
        %TESTMISSINGSAVEFOLDERPRESERVESSESSIONANDSTOREDPATH Runtime-only state.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            movedFolder = fullfile(testCase.TempRoot, 'TemporarilyOffline');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));
            movefile(saveFolder, movedFolder);

            status = store.getSaveFolderBindingStatus( ...
                'Mouse_104', 'Session_01');
            testCase.verifyEqual(status.state, 'missing');
            testCase.verifyFalse(status.isAvailable);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.processedDataFolder, saveFolder);
            SubjectInfo = store.getSubjectInfo('Mouse_104');
            testCase.verifyEqual(numel(SubjectInfo.sessionRegistry), 1);
        end

        function testAddSessionRequiresValidAcqInfos(testCase)
        %TESTADDSESSIONREQUIRESVALIDACQINFOS Reject non-dataset folders.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            testCase.verifyError(@() store.addSession( ...
                'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder)), ...
                'Umitoolbox:UMITProjectStore:notADataset');

            payload = 1;
            save(fullfile(saveFolder, 'AcqInfos.mat'), ...
                'payload', '-mat');
            testCase.verifyError(@() store.addSession( ...
                'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder)), ...
                'Umitoolbox:UMITProjectStore:invalidAcqInfos');

            AcqInfoStream = repmat(struct(), 2, 1);
            save(fullfile(saveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
            testCase.verifyError(@() store.addSession( ...
                'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder)), ...
                'Umitoolbox:UMITProjectStore:invalidAcqInfos');

            testCase.verifyEmpty( ...
                store.getSubjectInfo('Mouse_104').sessionRegistry);
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder), 'bound'));
        end

        function testBlankDisplayNamesDefaultToIDs(testCase)
        %TESTBLANKDISPLAYNAMESDEFAULTTOIDS Apply backend creation defaults.

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct( ...
                'subjectID', 'Mouse_104', ...
                'displayName', '   '));
            SubjectInfo = store.getSubjectInfoByUUID(subjectUUID);
            testCase.verifyEqual( ...
                SubjectInfo.displayName, 'Mouse_104');

            sessionUUID = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'displayName', ''));
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.displayName, 'Session_01');
        end

        function testExpandedMetadataDefaultsUpdatesAndLegacyNormalization(testCase)
        %TESTEXPANDEDMETADATAPersist expanded fields through supported APIs.

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            projectFields = {'institution', 'lab', 'principalInvestigator', ...
                'experimenters', 'projectStartDate', 'projectEndDate', 'notes'};
            subjectFields = {'species', 'strain', 'genotype', 'sex', ...
                'dateOfBirth', 'ageAtProjectEntry', 'ageReference', 'weight_g', ...
                'housing', 'diet', 'healthStatus', 'earTag', 'rfid', 'notes'};
            sessionFields = {'behavioralTask', 'stimulusNotes', 'pharmacology', ...
                'surgery', 'virus', 'anesthesia', 'notes'};
            testCase.verifyTrue(all(isfield(store.getProjectInfo(), projectFields)));
            testCase.verifyTrue(all(isfield(store.getSubjectInfoByUUID(subjectUUID), subjectFields)));
            testCase.verifyTrue(all(isfield(store.getSessionInfoByUUID(sessionUUID), sessionFields)));

            store.updateProjectMetadata(struct( ...
                'institution', 'Example Institute', ...
                'experimenters', {{'A. Researcher', 'B. Operator'}}, ...
                'projectStartDate', datetime(2026, 1, 2), ...
                'projectEndDate', datetime(2026, 12, 31), ...
                'notes', 'Project notes'));
            store.updateSubjectMetadata('Mouse_104', struct( ...
                'species', 'Mus musculus', 'dateOfBirth', datetime(2025, 1, 2), ...
                'weight_g', 24.5, 'housing', 'Pair housed', 'rfid', 'RFID-104'));
            store.updateSessionMetadata('Mouse_104', 'Session_01', struct( ...
                'behavioralTask', 'Whisker task', 'stimulusNotes', 'LED pulses', ...
                'anesthesia', 'None', 'notes', 'Session notes'));
            testCase.verifyError(@() store.updateSubjectMetadata( ...
                'Mouse_104', struct('weight_g', 'heavy')), ...
                'Umitoolbox:UMITProjectStore:updateSubjectFailed');

            projectUUID = store.getProjectInfo().projectUUID;
            reopened = UMITProjectStore.open(projectUUID);
            testCase.verifyEqual(reopened.getProjectInfo().experimenters, ...
                {'A. Researcher'; 'B. Operator'});
            testCase.verifyEqual(reopened.getSubjectInfoByUUID(subjectUUID).weight_g, 24.5);
            testCase.verifyEqual(reopened.getSessionInfoByUUID(sessionUUID).behavioralTask, ...
                'Whisker task');

            projectDescriptions = reopened.getMetadataFieldDescriptions('project');
            subjectDescriptions = reopened.getMetadataFieldDescriptions('subject');
            sessionDescriptions = reopened.getMetadataFieldDescriptions('session');
            testCase.verifyTrue(all(isfield(projectDescriptions, projectFields)));
            testCase.verifyTrue(all(isfield(subjectDescriptions, subjectFields)));
            testCase.verifyTrue(all(isfield(sessionDescriptions, sessionFields)));
            testCase.verifyEmpty(fieldnames(reopened.getMetadataFieldDescriptions('unknown')));

            ProjectInfo = reopened.getProjectInfo();
            ProjectInfo = rmfield(ProjectInfo, projectFields);
            save(fullfile(testCase.ProjectRoot, 'project.mat'), 'ProjectInfo', '-mat');
            SubjectInfo = reopened.getSubjectInfoByUUID(subjectUUID);
            SubjectInfo = rmfield(SubjectInfo, subjectFields);
            save(fullfile(testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'subject.mat'), 'SubjectInfo', '-mat');
            SessionInfo = reopened.getSessionInfoByUUID(sessionUUID);
            SessionInfo = rmfield(SessionInfo, sessionFields);
            save(fullfile(testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01', 'session.mat'), 'SessionInfo', '-mat');
            reopened = UMITProjectStore.open(projectUUID);
            testCase.verifyTrue(all(isfield(reopened.getProjectInfo(), projectFields)));
            testCase.verifyTrue(all(isfield( ...
                reopened.getSubjectInfoByUUID(subjectUUID), subjectFields)));
            testCase.verifyTrue(all(isfield( ...
                reopened.getSessionInfoByUUID(sessionUUID), sessionFields)));
        end

        function testRelocateSaveFolderPreservesSessionIdentity(testCase)
        %TESTRELOCATESAVEFOLDERPRESERVESSESSIONIDENTITY Rebind atomically.

            oldFolder = fullfile(testCase.TempRoot, 'OldSaveFolder');
            newFolder = fullfile(testCase.TempRoot, 'NewSaveFolder');
            mkdir(oldFolder);
            mkdir(newFolder);
            testCase.writeValidAcqInfos(oldFolder);
            testCase.writeValidAcqInfos(newFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'displayName', 'Original metadata', ...
                'processedDataFolder', oldFolder));

            newBinding = store.relocateSaveFolder( ...
                'Mouse_104', 'Session_01', newFolder);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.uuid, sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.displayName, 'Original metadata');
            testCase.verifyEqual( ...
                SessionInfo.processedDataFolder, newFolder);
            testCase.verifyEqual( ...
                SessionInfo.processedDataBindingUUID, ...
                newBinding.bindingUUID);
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(oldFolder), 'bound'));
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(newFolder), 'bound'));
        end

        function testRelocateFromUnavailableOldFolder(testCase)
        %TESTRELOCATEFROMUNAVAILABLEOLDFOLDER Preserve ledger until success.

            oldFolder = fullfile(testCase.TempRoot, 'OldSaveFolder');
            parkedFolder = fullfile(testCase.TempRoot, 'OfflineSaveFolder');
            newFolder = fullfile(testCase.TempRoot, 'NewSaveFolder');
            mkdir(oldFolder);
            mkdir(newFolder);
            testCase.writeValidAcqInfos(oldFolder);
            testCase.writeValidAcqInfos(newFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', oldFolder));
            movefile(oldFolder, parkedFolder);

            store.relocateSaveFolder( ...
                'Mouse_104', 'Session_01', newFolder);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.processedDataFolder, newFolder);
            testCase.verifyTrue(isfolder(parkedFolder));
        end

        function testRelocationConflictRollsBackLedgerAndLink(testCase)
        %TESTRELOCATIONCONFLICTROLLSBACKLEDGERANDLINK Reject another owner.

            folderOne = fullfile(testCase.TempRoot, 'SaveFolderOne');
            folderTwo = fullfile(testCase.TempRoot, 'SaveFolderTwo');
            mkdir(folderOne);
            mkdir(folderTwo);
            testCase.writeValidAcqInfos(folderOne);
            testCase.writeValidAcqInfos(folderTwo);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', folderOne));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_02', ...
                'processedDataFolder', folderTwo));
            before = store.getSessionInfo( ...
                'Mouse_104', 'Session_02');

            testCase.verifyError(@() store.relocateSaveFolder( ...
                'Mouse_104', 'Session_02', folderOne), ...
                'Umitoolbox:UMITProjectStore:dataFolderAlreadyBound');

            after = store.getSessionInfo('Mouse_104', 'Session_02');
            testCase.verifyEqual(after.processedDataFolder, ...
                before.processedDataFolder);
            testCase.verifyEqual(after.processedDataBindingUUID, ...
                before.processedDataBindingUUID);
            binding = store.getProcessedDataFolderBinding( ...
                'Mouse_104', 'Session_02');
            testCase.verifyEqual(binding.bindingUUID, ...
                before.processedDataBindingUUID);
        end

        function testInvalidSaveFolderRelocationRollsBack(testCase)
        %TESTINVALIDSAVEFOLDERRELOCATIONROLLSBACK Preserve old binding.

            oldFolder = fullfile(testCase.TempRoot, 'OldSaveFolder');
            invalidFolder = fullfile(testCase.TempRoot, ...
                'InvalidSaveFolder');
            mkdir(oldFolder);
            mkdir(invalidFolder);
            testCase.writeValidAcqInfos(oldFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', oldFolder));
            before = store.getSessionInfoByUUID(sessionUUID);

            testCase.verifyError(@() store.relocateSaveFolder( ...
                'Mouse_104', 'Session_01', invalidFolder), ...
                'Umitoolbox:UMITProjectStore:notADataset');

            after = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(after.processedDataFolder, ...
                before.processedDataFolder);
            testCase.verifyEqual(after.processedDataBindingUUID, ...
                before.processedDataBindingUUID);
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(oldFolder), 'bound'));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(invalidFolder), 'bound'));
        end

        function testValidationReportsMissingAcqInfos(testCase)
        %TESTVALIDATIONREPORTSMISSINGACQINFOS Flag accessible invalid datasets.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));

            delete(fullfile(saveFolder, 'AcqInfos.mat'));
            report = store.validate('Mode', 'quick');

            testCase.verifyFalse(report.isValid);
            testCase.verifyTrue(any(strcmp( ...
                {report.errors.code}, 'not_a_dataset')));
        end

        function testRelocateRepairsMissingAcqInfosState(testCase)
        %TESTRELOCATEREPAIRSMISSINGACQINFOSSTATE Escape invalid old dataset.

            oldFolder = fullfile(testCase.TempRoot, 'OldSaveFolder');
            newFolder = fullfile(testCase.TempRoot, 'NewSaveFolder');
            mkdir(oldFolder);
            mkdir(newFolder);
            testCase.writeValidAcqInfos(oldFolder);
            testCase.writeValidAcqInfos(newFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', oldFolder));
            delete(fullfile(oldFolder, 'AcqInfos.mat'));

            store.relocateSaveFolder( ...
                'Mouse_104', 'Session_01', newFolder);

            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.processedDataFolder, newFolder);
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(newFolder), 'bound'));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(oldFolder), 'bound'));
        end

        function testRelocateRawDataFolderIsAtomic(testCase)
        %TESTRELOCATERAWDATAFOLDERISATOMIC Preserve one reciprocal binding.

            oldRawFolder = fullfile(testCase.TempRoot, 'OldRawFolder');
            newRawFolder = fullfile(testCase.TempRoot, 'NewRawFolder');
            mkdir(oldRawFolder);
            mkdir(newRawFolder);
            testCase.writeRawDataFile(oldRawFolder);
            testCase.writeRawDataFile(newRawFolder, '.tiff');

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', oldRawFolder);

            newBinding = store.relocateRawDataFolder( ...
                'Mouse_104', 'Session_01', newRawFolder);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);

            testCase.verifyEqual( ...
                SessionInfo.rawDataFolder, newRawFolder);
            testCase.verifyEqual( ...
                SessionInfo.rawDataBindingUUID, ...
                newBinding.bindingUUID);
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(oldRawFolder), 'bound'));
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(newRawFolder), 'bound'));
        end

        function testRelocateRawFolderFromUnavailableOldFolder(testCase)
        %TESTRELOCATERAWFOLDERFROMUNAVAILABLEOLDFOLDER Keep raw data intact.

            oldRawFolder = fullfile(testCase.TempRoot, 'OldRawFolder');
            parkedRawFolder = fullfile(testCase.TempRoot, ...
                'ParkedRawFolder');
            newRawFolder = fullfile(testCase.TempRoot, 'NewRawFolder');
            mkdir(oldRawFolder);
            mkdir(newRawFolder);
            testCase.writeRawDataFile(oldRawFolder);
            testCase.writeRawDataFile(newRawFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', oldRawFolder);
            movefile(oldRawFolder, parkedRawFolder);

            store.relocateRawDataFolder( ...
                'Mouse_104', 'Session_01', newRawFolder);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.rawDataFolder, newRawFolder);
            testCase.verifyTrue(isfolder(parkedRawFolder));
        end

        function testRawFolderRelocationConflictRollsBack(testCase)
        %TESTRAWFOLDERRELOCATIONCONFLICTROLLSBACK Keep original raw binding.

            rawFolderOne = fullfile(testCase.TempRoot, 'RawFolderOne');
            rawFolderTwo = fullfile(testCase.TempRoot, 'RawFolderTwo');
            mkdir(rawFolderOne);
            mkdir(rawFolderTwo);
            testCase.writeRawDataFile(rawFolderOne);
            testCase.writeRawDataFile(rawFolderTwo);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionOne = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_02'));
            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', rawFolderOne);
            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_02', rawFolderTwo);
            before = store.getSessionInfoByUUID(sessionOne);

            testCase.verifyError(@() ...
                store.relocateRawDataFolder( ...
                'Mouse_104', 'Session_01', rawFolderTwo), ...
                'Umitoolbox:UMITProjectStore:dataFolderAlreadyBound');

            after = store.getSessionInfoByUUID(sessionOne);
            testCase.verifyEqual( ...
                after.rawDataFolder, before.rawDataFolder);
            testCase.verifyEqual( ...
                after.rawDataBindingUUID, ...
                before.rawDataBindingUUID);
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(rawFolderOne), 'bound'));
        end

        function testRawFolderRequiresRecognizableImagingFile(testCase)
        %TESTRAWFOLDERREQUIRESRECOGNIZABLEIMAGINGFILE Reject empty folders.

            rawFolder = fullfile(testCase.TempRoot, 'EmptyRawFolder');
            mkdir(rawFolder);
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            testCase.verifyError(@() store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', rawFolder), ...
                'Umitoolbox:UMITProjectStore:missingRawDataFiles');
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEmpty(SessionInfo.rawDataFolder);
            testCase.verifyEmpty(SessionInfo.rawDataBindingUUID);
        end

        function testRawFolderMayEqualSaveFolder(testCase)
        %TESTRAWFOLDERMAYEQUALSAVEFOLDER Share one reciprocal link.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            before = store.getSessionInfoByUUID(sessionUUID);
            testCase.writeRawDataFile(before.processedDataFolder, '.tif');

            binding = store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', ...
                before.processedDataFolder);
            after = store.getSessionInfoByUUID(sessionUUID);

            testCase.verifyEqual( ...
                after.rawDataFolder, after.processedDataFolder);
            testCase.verifyEqual( ...
                after.rawDataBindingUUID, ...
                after.processedDataBindingUUID);
            testCase.verifyEqual( ...
                binding.bindingUUID, after.rawDataBindingUUID);
            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testRelocateSaveFolderAwayFromSharedRawFolder(testCase)
        %TESTRELOCATESAVEFOLDERAWAYFROMSHAREDRAWFOLDER Preserve both links.

            newSaveFolder = fullfile( ...
                testCase.TempRoot, 'RelocatedSaveFolder');
            mkdir(newSaveFolder);
            testCase.writeValidAcqInfos(newSaveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession( ...
                store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            before = store.getSessionInfoByUUID(sessionUUID);
            oldSharedFolder = before.processedDataFolder;
            testCase.writeRawDataFile(oldSharedFolder);
            store.bindRawDataFolder( ...
                'Mouse_104', 'Session_01', oldSharedFolder);

            store.relocateSaveFolder( ...
                'Mouse_104', 'Session_01', newSaveFolder);
            after = store.getSessionInfoByUUID(sessionUUID);
            rawBinding = store.getRawDataFolderBinding( ...
                'Mouse_104', 'Session_01');

            testCase.verifyEqual( ...
                after.rawDataFolder, oldSharedFolder);
            testCase.verifyEqual( ...
                after.processedDataFolder, newSaveFolder);
            testCase.verifyEqual( ...
                rawBinding.folderRole, 'rawDataFolder');
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(oldSharedFolder), 'bound'));
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(newSaveFolder), 'bound'));
            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testRepairMissingSaveFolderBinding(testCase)
        %TESTREPAIRMISSINGSAVEFOLDERBINDING Restore store-owned link.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));
            linkPath = fullfile(saveFolder, ...
                UMITProjectStore.getProjectBindingFileName());
            delete(linkPath);

            status = store.getSaveFolderBindingStatus( ...
                'Mouse_104', 'Session_01');
            testCase.verifyEqual(status.state, 'invalid');
            repaired = store.repairSaveFolderBinding( ...
                'Mouse_104', 'Session_01');
            testCase.verifyTrue(isfile(linkPath));
            testCase.verifyEqual(repaired.bindingUUID, ...
                store.getSessionInfo( ...
                'Mouse_104', 'Session_01').processedDataBindingUUID);
        end

        function testRemoveSessionDoesNotDeleteScientificData(testCase)
        %TESTREMOVESESSIONDOESNOTDELETESCIENTIFICDATA Explicit ledger removal.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            scientificFile = fullfile(saveFolder, 'green.dat');
            fid = fopen(scientificFile, 'w');
            fwrite(fid, uint8([1, 2, 3]));
            fclose(fid);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));

            removed = store.removeSessionFromProject( ...
                'Mouse_104', 'Session_01');
            testCase.verifyEqual(removed.uuid, sessionUUID);
            testCase.verifyTrue(isfile(scientificFile));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder), 'bound'));
            SubjectInfo = store.getSubjectInfo('Mouse_104');
            testCase.verifyEmpty(SubjectInfo.sessionRegistry);
            testCase.verifyError(@() ...
                store.getSessionInfoByUUID(sessionUUID), ...
                'Umitoolbox:UMITProjectStore:sessionNotFound');
        end

        function testRemoveSessionWithUnavailableFolder(testCase)
        %TESTREMOVESESSIONWITHUNAVAILABLEFOLDER Remove ledger without data access.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            parkedFolder = fullfile(testCase.TempRoot, 'OfflineSaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));
            movefile(saveFolder, parkedFolder);

            store.removeSessionFromProject('Mouse_104', 'Session_01');
            testCase.verifyTrue(isfolder(parkedFolder));
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(parkedFolder), 'bound'));
            testCase.verifyEmpty( ...
                store.getSubjectInfo('Mouse_104').sessionRegistry);
        end

        function testRemoveEmptySubjectFromProject(testCase)
        %TESTREMOVEEMPTYSUBJECTFROMPROJECT Remove an empty subject by UUID.

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse_104'));

            result = store.removeSubjectFromProject(subjectUUID);

            testCase.verifyEqual(result.subjectUUID, subjectUUID);
            testCase.verifyEqual(result.sessionsRemoved, 0);
            testCase.verifyEmpty(store.getProjectInfo().subjectRegistry);
            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testRemoveSubjectUnbindsMultipleSessions(testCase)
        %TESTREMOVESUBJECTUNBINDSMULTIPLESESSIONS Preserve external data.

            saveFolder1 = fullfile(testCase.TempRoot, 'SaveFolder_1');
            saveFolder2 = fullfile(testCase.TempRoot, 'SaveFolder_2');
            mkdir(saveFolder1);
            mkdir(saveFolder2);
            testCase.writeValidAcqInfos(saveFolder1);
            testCase.writeValidAcqInfos(saveFolder2);
            scientificFile = fullfile(saveFolder1, 'green.dat');
            fclose(fopen(scientificFile, 'w'));

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder1));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_02', ...
                'processedDataFolder', saveFolder2));
            otherUUID = store.addSubject(struct('subjectID', 'Mouse_105'));

            result = store.removeSubjectFromProject(subjectUUID);

            testCase.verifyEqual(result.sessionsRemoved, 2);
            testCase.verifyTrue(isfile(scientificFile));
            testCase.verifyTrue(isfolder(saveFolder2));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder1), 'bound'));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder2), 'bound'));
            testCase.verifyEqual(store.getSubjectInfoByUUID(otherUUID).uuid, ...
                otherUUID);
            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testRemoveSubjectFailureLeavesLedgerUnchanged(testCase)
        %TESTREMOVESUBJECTFAILURELEAVESLEDGERUNCHANGED Reject unsafe cascade.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            offlineFolder = fullfile(testCase.TempRoot, 'OfflineSaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));
            movefile(saveFolder, offlineFolder);

            testCase.verifyError(@() store.removeSubjectFromProject(subjectUUID), ...
                'Umitoolbox:UMITProjectStore:bindingFolderUnavailable');
            testCase.verifyEqual(store.getSubjectInfoByUUID(subjectUUID).uuid, ...
                subjectUUID);
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(offlineFolder), 'bound'));
        end

        function testDeleteEmptyProject(testCase)
        %TESTDELETEEMPTYPROJECT Remove only centralized project metadata.

            store = testCase.createProject();
            result = store.deleteProject(testCase.ProjectUUID);

            testCase.verifyEqual(result.subjectsRemoved, 0);
            testCase.verifyFalse(isfolder(testCase.ProjectRoot));
            projects = UMITProjectStore.listProjects();
            testCase.verifyFalse(any(strcmpi( ...
                projects.ProjectUUID, testCase.ProjectUUID)));
        end

        function testDeleteProjectWithEmptySubjectsAndManagedResources(testCase)
        %TESTDELETEPROJECTWITHEMPTYSUBJECTSANDMANAGEDRESOURCES Cascade safely.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_105'));
            source = testCase.createMATSource('reference.mat', 1);
            store.addImageReference('Mouse_104', source, struct());

            result = store.deleteProject(testCase.ProjectUUID);

            testCase.verifyEqual(result.subjectsRemoved, 2);
            testCase.verifyEqual(result.sessionsRemoved, 0);
            testCase.verifyTrue(isfile(source));
            testCase.verifyFalse(isfolder(testCase.ProjectRoot));
        end

        function testDeleteProjectUnbindsMultipleSubjectsAndSessions(testCase)
        %TESTDELETEPROJECTUNBINDSMULTIPLESUBJECTSANDSESSIONS Preserve data.

            saveFolder1 = fullfile(testCase.TempRoot, 'SaveFolder_1');
            saveFolder2 = fullfile(testCase.TempRoot, 'SaveFolder_2');
            mkdir(saveFolder1);
            mkdir(saveFolder2);
            testCase.writeValidAcqInfos(saveFolder1);
            testCase.writeValidAcqInfos(saveFolder2);
            scientificFile = fullfile(saveFolder2, 'science.dat');
            fclose(fopen(scientificFile, 'w'));
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_105'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder1));
            testCase.addSession(store, 'Mouse_105', struct( ...
                'sessionID', 'Session_02', ...
                'processedDataFolder', saveFolder2));
            otherStore = testCase.createSecondProject();

            result = store.deleteProject(testCase.ProjectUUID);

            testCase.verifyEqual(result.subjectsRemoved, 2);
            testCase.verifyEqual(result.sessionsRemoved, 2);
            testCase.verifyTrue(isfile(scientificFile));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder1), 'bound'));
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder2), 'bound'));
            testCase.verifyEqual(otherStore.getProjectInfo().projectUUID, ...
                testCase.ProjectUUID2);
        end

        function testDeleteProjectSkipsOfflineAndBrokenSessions(testCase)
        %TESTDELETEPROJECTSKIPSOFFLINEANDBROKENSESSIONS Audit item 23:
        %deleteProject works on corrupted or offline projects. Offline
        %folders and unreadable sessions are skipped and reported; a skipped
        %folder keeps its ledger, which then reads as an orphan binding.
        %Staged-delete leftovers never appear in listings.

            onlineFolder = fullfile(testCase.TempRoot, 'OnlineSaveFolder');
            offlineFolder = fullfile(testCase.TempRoot, 'OfflineSaveFolder');
            movedFolder = fullfile(testCase.TempRoot, 'MovedSaveFolder');
            brokenFolder = fullfile(testCase.TempRoot, 'BrokenSessionFolder');
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Online', 'processedDataFolder', onlineFolder));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Offline', 'processedDataFolder', offlineFolder));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Broken', 'processedDataFolder', brokenFolder));
            movefile(offlineFolder, movedFolder);
            testCase.writeText(testCase.resolveRelative( ...
                'subjects/Mouse_104/sessions/Broken/session.mat'), 'not a MAT file');

            reopened = UMITProjectStore.open(testCase.ProjectUUID);
            testCase.verifyFalse(reopened.IsReadOnly, ...
                'session errors do not make the project read-only');
            result = reopened.deleteProject(testCase.ProjectUUID);

            testCase.verifyFalse(isfolder(testCase.ProjectRoot));
            testCase.verifyEqual(numel(result.skipped), 2);
            testCase.verifyTrue(any(strcmp({result.skipped.folder}, offlineFolder)));
            testCase.verifyEqual(UMITProjectStore.getProjectBindingState(onlineFolder), 'unbound');
            testCase.verifyEqual(UMITProjectStore.getProjectBindingState(brokenFolder), 'bound', ...
                'unreadable session: its folder was not touched');
            testCase.verifyTrue(UMITProjectStore.isOrphanProjectBinding(movedFolder));
            testCase.verifyTrue(UMITProjectStore.isOrphanProjectBinding(brokenFolder));

            % A crashed delete leaves a dot-folder: never listed or resolved.
            leftover = fullfile(UMITProjectStore.getProjectsRoot(), '.deleteProjectStage_test');
            mkdir(leftover);
            testCase.addTeardown(@() rmdir(leftover, 's'));
            ProjectInfo = struct('projectUUID', testCase.ProjectUUID, 'schemaVersion', 1);
            save(fullfile(leftover, 'project.mat'), 'ProjectInfo');
            projects = UMITProjectStore.listProjects();
            testCase.verifyFalse(any(strcmpi(projects.ProjectUUID, testCase.ProjectUUID)));
            testCase.verifyError(@() UMITProjectStore.resolveProjectRoot(testCase.ProjectUUID), ...
                'Umitoolbox:UMITProjectStore:resolveProjectFailed');
            testCase.ProjectRoot = '';
            testCase.ProjectUUID = '';
        end

        function testDetectAndRemoveOrphanProjectBinding(testCase)
        %TESTDETECTANDREMOVEORPHANPROJECTBINDING Remove link to missing project.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            originalBinding = store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            clear store
            rmdir(testCase.ProjectRoot, 's');
            testCase.ProjectRoot = '';
            testCase.ProjectUUID = '';

            [isOrphan, detectedBinding] = ...
                UMITProjectStore.isOrphanProjectBinding(saveFolder);
            testCase.verifyTrue(isOrphan);
            testCase.verifyEqual(detectedBinding.bindingUUID, ...
                originalBinding.bindingUUID);

            removedBinding = ...
                UMITProjectStore.removeOrphanProjectBinding(saveFolder);
            testCase.verifyEqual(removedBinding.bindingUUID, ...
                originalBinding.bindingUUID);
            testCase.verifyFalse(strcmp(UMITProjectStore.getProjectBindingState(saveFolder), 'bound'));
            % The ledger and its history stay; the removal is an event.
            ledger = UMITProjectLedger.read(saveFolder);
            testCase.verifyEqual(ledger.state, 'Unbound');
            testCase.verifyEqual(ledger.events{end}.type, 'orphanRemoval');
            testCase.verifyEqual(ledger.events{end}.previousBinding.bindingUUID, ...
                originalBinding.bindingUUID);
        end

        function testRefuseOrphanRemovalWhenProjectIsAvailable(testCase)
        %TESTREFUSEORPHANREMOVALWHENPROJECTISAVAILABLE Preserve reciprocal state.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            isOrphan = ...
                UMITProjectStore.isOrphanProjectBinding(saveFolder);
            testCase.verifyFalse(isOrphan);

            testCase.verifyError(@() ...
                UMITProjectStore.removeOrphanProjectBinding(saveFolder), ...
                'Umitoolbox:UMITProjectStore:bindingProjectAvailable');
        end

        function testReplaceOrphanBindingTransactionally(testCase)
        %TESTREPLACEORPHANBINDINGTRANSACTIONALLY Reassign to another project.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            oldStore = testCase.createProject();
            oldStore.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(oldStore, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            oldBinding = oldStore.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            clear oldStore
            rmdir(testCase.ProjectRoot, 's');
            testCase.ProjectRoot = '';
            testCase.ProjectUUID = '';

            newStore = testCase.createProject();
            newStore.addSubject(struct('subjectID', 'Mouse_205'));
            testCase.addSession(newStore, 'Mouse_205', struct( ...
                'sessionID', 'Session_02'));

            testCase.verifyError(@() ...
                newStore.bindProcessedDataFolder( ...
                'Mouse_205', 'Session_02', saveFolder), ...
                'Umitoolbox:UMITProjectStore:bindingConflict');

            % Audit item 13: Session_02 is still healthily bound to its own
            % SaveFolder, so it cannot take over the orphan folder...
            testCase.verifyError(@() newStore.bindProcessedDataFolder( ...
                'Mouse_205', 'Session_02', saveFolder, ...
                'ReplaceOrphanBinding', true), ...
                'Umitoolbox:UMITProjectStore:sessionAlreadyBound');
            % ...until that folder is gone (the data were moved here).
            ownFolder = newStore.getSessionInfo('Mouse_205', 'Session_02').processedDataFolder;
            movefile(ownFolder, [ownFolder '_gone']);

            newBinding = newStore.bindProcessedDataFolder( ...
                'Mouse_205', 'Session_02', saveFolder, ...
                'ReplaceOrphanBinding', true);

            testCase.verifyNotEqual(newBinding.bindingUUID, ...
                oldBinding.bindingUUID);

            [context, resolvedStore] = ...
                UMITProjectStore.resolveProjectBinding(saveFolder);
            NewProjectInfo = newStore.getProjectInfo();
            testCase.verifyEqual(context.projectUUID, ...
                NewProjectInfo.projectUUID);
            testCase.verifyEqual(context.subjectID, 'Mouse_205');
            testCase.verifyEqual(context.sessionID, 'Session_02');
            testCase.verifyEqual(resolvedStore.ProjectRoot, ...
                newStore.ProjectRoot);

            SessionInfo = newStore.getSessionInfo( ...
                'Mouse_205', 'Session_02');
            testCase.verifyEqual( ...
                SessionInfo.processedDataBindingUUID, ...
                newBinding.bindingUUID);
        end

        function testAddSessionRequiresSaveFolderAndRejectsDirectUpdate(testCase)
        %TESTADDSESSIONREQUIRESSAVEFOLDERANDREJECTSDIRECTUPDATE Enforce owner.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            testCase.verifyError(@() store.addSession( ...
                'Mouse_104', struct( ...
                'sessionID', 'Session_01')), ...
                'Umitoolbox:UMITProjectStore:saveFolderRequired');

            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual( ...
                SessionInfo.processedDataFolder, saveFolder);
            testCase.verifyTrue(strcmp(UMITProjectStore.getProjectBindingState(saveFolder), 'bound'));

            testCase.verifyError(@() store.updateSessionMetadata( ...
                'Mouse_104', 'Session_01', struct( ...
                'processedDataFolder', fullfile(saveFolder, 'other'))), ...
                'Umitoolbox:UMITProjectStore:folderBindingRequired');
        end

        function testCopiedBindingPathMismatchIsDetected(testCase)
        %TESTCOPIEDBINDINGPATHMISMATCHISDETECTED Detect copied bound folders.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            copiedFolder = fullfile(testCase.TempRoot, 'CopiedFolder');
            mkdir(saveFolder);
            mkdir(copiedFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            bindingName = ...
                UMITProjectStore.getProjectBindingFileName();
            copyfile(fullfile(saveFolder, bindingName), ...
                fullfile(copiedFolder, bindingName), 'f');

            testCase.verifyError(@() ...
                UMITProjectStore.resolveProjectBinding(copiedFolder), ...
                'Umitoolbox:UMITProjectStore:bindingPathMismatch');
        end

        function testRejectInvalidResourceFilters(testCase)
        %TESTREJECTINVALIDRESOURCEFILTERS Reject incompatible query filters.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            testCase.verifyError(@() ...
                store.listSubjectResources( ...
                'Mouse_104', 'Type', 'cameraCoregistration'), ...
                'Umitoolbox:UMITProjectStore:invalidResourceFilter');
            testCase.verifyError(@() ...
                store.listSubjectResources( ...
                'Mouse_104', 'Status', 'deleted'), ...
                'Umitoolbox:UMITProjectStore:invalidResourceFilter');
        end


        function testRejectMalformedImageReferencePayload(testCase)
        %TESTREJECTMALFORMEDIMAGEREFERENCEPAYLOAD Enforce resource contract.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            badFile = fullfile( ...
                testCase.SourceFolder, 'malformed_image.mat');
            payload = 42;
            save(badFile, 'payload', '-mat');

            testCase.verifyError(@() store.addImageReference( ...
                'Mouse_104', badFile, struct()), ...
                'Umitoolbox:UMITProjectStore:invalidResourcePayload');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testAcceptCanonicalImageReferencePayload(testCase)
        %TESTACCEPTCANONICALIMAGEREFERENCEPAYLOAD Import validated reference.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));

            sourceFile = testCase.createMATSource( ...
                'reference_valid.mat', 7);

            resourceUUID = store.addImageReference( ...
                'Mouse_104', sourceFile, ...
                struct('displayName', 'Validated Reference'));

            resource = store.getResource(resourceUUID);
            testCase.verifyTrue(resource.fileExists);

            loaded = load( ...
                resource.absolutePath, 'ImageReference', '-mat');
            testCase.verifyTrue(isfield(loaded, 'ImageReference'));
            testCase.verifyWarningFree(@() ...
                validateImageReferenceStruct(loaded.ImageReference));
            testCase.verifyEqual(loaded.ImageReference.version, 2);
        end

        function testImageReferenceImportChecksIdentityAndVersion(testCase)
        %TESTIMAGEREFERENCEIMPORTCHECKSIDENTITYANDVERSION Audit item 21.
        %   Identity fields present in the payload must name this project
        %   and the owner subject; libraries store version 2 only.

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_105'));
            ProjectInfo = store.getProjectInfo();
            image = single(magic(5));

            matching = genImageReferenceStruct(image, 'Name', 'own', ...
                'ProjectUUID', ProjectInfo.projectUUID, 'ProjectName', ProjectInfo.projectName, ...
                'SubjectUUID', subjectUUID, 'SubjectID', 'Mouse_104');
            ImageReference = matching;
            ownFile = fullfile(testCase.SourceFolder, 'own_ref.mat');
            save(ownFile, 'ImageReference', '-mat');
            testCase.verifyNotEmpty(store.addImageReference('Mouse_104', ownFile, struct()));
            % The same payload is refused by another subject of the project.
            testCase.verifyError(@() store.addImageReference('Mouse_105', ownFile, struct()), ...
                'Umitoolbox:UMITProjectStore:imageReferenceIdentityMismatch');

            ImageReference = matching;
            ImageReference.projectUUID = lower(char(java.util.UUID.randomUUID()));
            ImageReference.subjectUUID = '';
            ImageReference.subjectID = '';
            foreignFile = fullfile(testCase.SourceFolder, 'foreign_ref.mat');
            save(foreignFile, 'ImageReference', '-mat');
            testCase.verifyError(@() store.addImageReference('Mouse_104', foreignFile, struct()), ...
                'Umitoolbox:UMITProjectStore:imageReferenceIdentityMismatch');
            % Cleared identity: accepted (adoption of a foreign reference).
            testCase.verifyNotEmpty(importImageReferenceFile(store, 'Mouse_104', foreignFile, ...
                'ClearIdentity', true));

            ImageReference = rmfield(matching, 'contentHash');
            ImageReference.version = 1;
            v1File = fullfile(testCase.SourceFolder, 'v1_ref.mat');
            save(v1File, 'ImageReference', '-mat');
            testCase.verifyError(@() store.addImageReference('Mouse_104', v1File, struct()), ...
                'Umitoolbox:UMITProjectStore:imageReferenceVersion');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, testCase.issueMessages(report.errors));
        end

        function testRenameEntityByUUID(testCase)
        %TESTRENAMEENTITYBYUUID Verify UUID-keyed rename dispatch.

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', '2026-07-10_S01'));

            store.renameEntity('subject', subjectUUID, 'Mouse_104A');
            SubjectInfo = store.getSubjectInfoByUUID(subjectUUID);
            testCase.verifyEqual(SubjectInfo.subjectID, 'Mouse_104A');
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104A')));

            store.renameEntity('session', sessionUUID, '2026-07-10_S02');
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.sessionID, '2026-07-10_S02');
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104A', ...
                'sessions', '2026-07-10_S02')));

            testCase.verifyError(@() ...
                store.renameEntity('rig', subjectUUID, 'x'), ...
                'Umitoolbox:UMITProjectStore:renameEntityFailed');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testMoveSessionToSubjectHappyPath(testCase)
        %TESTMOVESESSIONTOSUBJECTHAPPYPATH Verify within-project re-parenting.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_205'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            dataFolder = fullfile(testCase.TempRoot, 'RawData');
            mkdir(dataFolder);
            testCase.writeValidAcqInfos(dataFolder);
            store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', dataFolder);

            targetSubjectUUID = store.getSubjectInfo('Mouse_205').uuid;

            store.moveSessionToSubject(sessionUUID, targetSubjectUUID);

            testCase.verifyFalse(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01')));
            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_205', ...
                'sessions', 'Session_01')));

            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.subjectID, 'Mouse_205');
            testCase.verifyEqual(SessionInfo.subjectUUID, targetSubjectUUID);
            testCase.verifyEqual(SessionInfo.sessionID, 'Session_01');

            OldSubjectInfo = store.getSubjectInfo('Mouse_104');
            testCase.verifyEmpty(OldSubjectInfo.sessionRegistry);
            NewSubjectInfo = store.getSubjectInfo('Mouse_205');
            testCase.verifyEqual( ...
                numel(NewSubjectInfo.sessionRegistry), 1);

            ProjectBinding = store.getProcessedDataFolderBinding( ...
                'Mouse_205', 'Session_01');
            testCase.verifyEqual( ...
                ProjectBinding.subjectUUID, targetSubjectUUID);

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testMoveSessionToSubjectIsNoOpForSameSubject(testCase)
        %TESTMOVESESSIONTOSUBJECTISNOOPFORSAMESUBJECT No-op re-parent guard.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            subjectUUID = store.getSubjectInfo('Mouse_104').uuid;

            testCase.verifyWarningFree(@() ...
                store.moveSessionToSubject(sessionUUID, subjectUUID));

            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.subjectID, 'Mouse_104');
        end

        function testMoveSessionToSubjectRejectsInvalidTarget(testCase)
        %TESTMOVESESSIONTOSUBJECTREJECTSINVALIDTARGET Reject bad UUIDs.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            bogusUUID = lower(char(java.util.UUID.randomUUID()));
            testCase.verifyError(@() ...
                store.moveSessionToSubject(sessionUUID, bogusUUID), ...
                'Umitoolbox:UMITProjectStore:moveSessionFailed');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testMoveSessionToSubjectRejectsSessionIDCollision(testCase)
        %TESTMOVESESSIONTOSUBJECTREJECTSSESSIONIDCOLLISION Reject ID clash.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_205'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            testCase.addSession(store, 'Mouse_205', struct( ...
                'sessionID', 'Session_01'));

            targetSubjectUUID = store.getSubjectInfo('Mouse_205').uuid;
            testCase.verifyError(@() ...
                store.moveSessionToSubject(sessionUUID, targetSubjectUUID), ...
                'Umitoolbox:UMITProjectStore:duplicateID');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testMoveSessionToSubjectRollbackOnFailure(testCase)
        %TESTMOVESESSIONTOSUBJECTROLLBACKONFAILURE Verify full rollback.
        %
        %   Forces the final metadata write (the destination subject's
        %   subject.mat) to fail by holding it open, after the session
        %   folder has already been physically relocated and its own
        %   SessionInfo/bindings/source-subject registry have already been
        %   rewritten. This exercises the catch-block rollback of every
        %   step, not just an upfront precondition rejection.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            store.addSubject(struct('subjectID', 'Mouse_205'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            targetSubjectUUID = store.getSubjectInfo('Mouse_205').uuid;

            schema = getUMITProjectSchema();
            targetMetadataPath = fullfile(testCase.ProjectRoot, ...
                'subjects', 'Mouse_205', schema.files.subjectMetadata);

            fid = fopen(targetMetadataPath, 'r');
            testCase.assertNotEqual(fid, -1);
            fidCleanup = onCleanup(@() fclose(fid));

            testCase.verifyError(@() ...
                store.moveSessionToSubject(sessionUUID, targetSubjectUUID), ...
                'Umitoolbox:saveMatAtomic:backupFailed');

            clear fidCleanup

            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01')));
            testCase.verifyFalse(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_205', ...
                'sessions', 'Session_01')));

            OldSubjectInfo = store.getSubjectInfo('Mouse_104');
            testCase.verifyEqual( ...
                numel(OldSubjectInfo.sessionRegistry), 1);
            NewSubjectInfo = store.getSubjectInfo('Mouse_205');
            testCase.verifyEmpty(NewSubjectInfo.sessionRegistry);

            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.subjectID, 'Mouse_104');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testRebindSessionToProjectHappyPath(testCase)
        %TESTREBINDSESSIONTOPROJECTHAPPYPATH Verify cross-project move.

            sourceStore = testCase.createProject();
            sourceStore.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(sourceStore, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            dataFolder = fullfile(testCase.TempRoot, 'RawData');
            mkdir(dataFolder);
            testCase.writeValidAcqInfos(dataFolder);
            sourceStore.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', dataFolder);

            targetStore = testCase.createSecondProject();

            sourceStore.rebindSessionToProject( ...
                sessionUUID, testCase.ProjectUUID2);

            testCase.verifyFalse(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01')));
            OldSubjectInfo = sourceStore.getSubjectInfo('Mouse_104');
            testCase.verifyEmpty(OldSubjectInfo.sessionRegistry);

            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot2, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01')));

            [SessionInfo, NewSubjectInfo] = ...
                targetStore.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.sessionID, 'Session_01');
            testCase.verifyEqual(NewSubjectInfo.subjectID, 'Mouse_104');
            testCase.verifyEqual( ...
                SessionInfo.subjectUUID, NewSubjectInfo.uuid);

            ProjectBinding = targetStore.getProcessedDataFolderBinding( ...
                'Mouse_104', 'Session_01');
            NewProjectInfo = targetStore.getProjectInfo();
            testCase.verifyEqual( ...
                ProjectBinding.projectUUID, NewProjectInfo.projectUUID);
            testCase.verifyEqual( ...
                ProjectBinding.subjectUUID, NewSubjectInfo.uuid);

            sourceReport = sourceStore.validate('Mode', 'full');
            testCase.verifyTrue(sourceReport.isValid, ...
                testCase.issueMessages(sourceReport.errors));
            targetReport = targetStore.validate('Mode', 'full');
            testCase.verifyTrue(targetReport.isValid, ...
                testCase.issueMessages(targetReport.errors));
        end

        function testRebindSessionToProjectIsNoOpForSameProject(testCase)
        %TESTREBINDSESSIONTOPROJECTISNOOPFORSAMEPROJECT No-op guard.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            testCase.verifyWarningFree(@() ...
                store.rebindSessionToProject(sessionUUID, testCase.ProjectUUID));

            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.subjectID, 'Mouse_104');
        end

        function testRebindSessionToProjectRejectsInvalidTarget(testCase)
        %TESTREBINDSESSIONTOPROJECTREJECTSINVALIDTARGET Reject bad target.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            bogusProjectUUID = lower(char(java.util.UUID.randomUUID()));
            testCase.verifyError(@() ...
                store.rebindSessionToProject(sessionUUID, bogusProjectUUID), ...
                'Umitoolbox:UMITProjectStore:resolveProjectFailed');

            OldSubjectInfo = store.getSubjectInfo('Mouse_104');
            testCase.verifyEqual( ...
                numel(OldSubjectInfo.sessionRegistry), 1);
        end

        function testRebindSessionToProjectRollbackOnAttachFailure(testCase)
        %TESTREBINDSESSIONTOPROJECTROLLBACKONATTACHFAILURE Verify phase-2
        %failure restores the session into the source project.
        %
        %   The target subject is pre-created so its subject.mat already
        %   exists and can be locked open, forcing the attach phase's final
        %   metadata write to fail AFTER the session has already been
        %   detached from the source project (phase 1 committed) and moved
        %   into the target project's staging location. Confirms the
        %   documented compensating behavior: the target is rolled back to
        %   not referencing the session, and the session is automatically
        %   restored into the source project rather than left stranded.

            sourceStore = testCase.createProject();
            sourceStore.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(sourceStore, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            targetStore = testCase.createSecondProject();
            targetStore.addSubject(struct('subjectID', 'Mouse_104'));

            schema = getUMITProjectSchema();
            targetMetadataPath = fullfile(testCase.ProjectRoot2, ...
                'subjects', 'Mouse_104', schema.files.subjectMetadata);

            fid = fopen(targetMetadataPath, 'r');
            testCase.assertNotEqual(fid, -1);
            fidCleanup = onCleanup(@() fclose(fid));

            testCase.verifyError(@() sourceStore.rebindSessionToProject( ...
                sessionUUID, testCase.ProjectUUID2), ...
                'Umitoolbox:UMITProjectStore:rebindSessionFailed');

            clear fidCleanup

            testCase.verifyTrue(isfolder(fullfile( ...
                testCase.ProjectRoot, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01')));
            OldSubjectInfo = sourceStore.getSubjectInfo('Mouse_104');
            testCase.verifyEqual( ...
                numel(OldSubjectInfo.sessionRegistry), 1);
            SessionInfo = sourceStore.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.subjectID, 'Mouse_104');

            testCase.verifyFalse(isfolder(fullfile( ...
                testCase.ProjectRoot2, 'subjects', 'Mouse_104', ...
                'sessions', 'Session_01')));
            NewSubjectInfo = targetStore.getSubjectInfo('Mouse_104');
            testCase.verifyEmpty(NewSubjectInfo.sessionRegistry);

            sourceReport = sourceStore.validate('Mode', 'full');
            testCase.verifyTrue(sourceReport.isValid, ...
                testCase.issueMessages(sourceReport.errors));
            targetReport = targetStore.validate('Mode', 'full');
            testCase.verifyTrue(targetReport.isValid, ...
                testCase.issueMessages(targetReport.errors));
        end

    end

    methods (Access = private)
        function sessionUUID = addSession( ...
                testCase, store, subjectID, sessionInfo)
        %ADDSESSION Create a dataset-backed session fixture.
        %
        %   Sessions keep no rig. A rigID/rigUUID field given to this
        %   fixture binds that rig in the SaveFolder's rig ledger
        %   (UMITRigLedger) and is removed before addSession is called.

            if ~isfield(sessionInfo, 'processedDataFolder')
                saveFolder = fullfile(testCase.TempRoot, ...
                    ['SaveFolder_' char(java.util.UUID.randomUUID())]);
                mkdir(saveFolder);
                sessionInfo.processedDataFolder = saveFolder;
            end
            testCase.writeValidAcqInfos( ...
                sessionInfo.processedDataFolder);

            requestedRigFields = intersect(fieldnames(sessionInfo), ...
                {'rigID', 'rigUUID'});
            if ~isempty(requestedRigFields)
                if isfield(sessionInfo, 'rigUUID')
                    rigStore = UMITRigStore.open(sessionInfo.rigUUID);
                else
                    rigStore = UMITRigStore.openByRigID(sessionInfo.rigID);
                end
                testCase.bindFolderRig(sessionInfo.processedDataFolder, rigStore);
                sessionInfo = rmfield(sessionInfo, requestedRigFields);
            end

            sessionUUID = store.addSession(subjectID, sessionInfo);
        end

        function message = errorMessage(testCase, fcn, expectedID)
        %ERRORMESSAGE Run FCN, verify it throws EXPECTEDID, return the message.

            message = '';
            try
                fcn();
                testCase.verifyFail(sprintf('Expected error %s.', expectedID));
            catch ME
                testCase.verifyEqual(ME.identifier, expectedID);
                message = ME.message;
            end
        end

        function bindFolderRig(~, saveFolder, rigStore)
        %BINDFOLDERRIG Bind (or rebind) a SaveFolder's rig ledger to a library rig.

            binding = UMITRigLedger.bindingFromRig(rigStore, 'IncludeCoregistration', false);
            if UMITRigLedger.exists(saveFolder)
                UMITRigLedger.rebind(saveFolder, binding, 'Origin', 'TestUMITProjectStore');
            else
                UMITRigLedger.bind(saveFolder, binding, 'Origin', 'TestUMITProjectStore');
            end
        end

        function writeValidAcqInfos(~, saveFolder)
        %WRITEVALIDACQINFOS Create minimum canonical dataset metadata.

            if ~isfolder(saveFolder)
                mkdir(saveFolder);
            end
            AcqInfoStream = struct();
            save(fullfile(saveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
        end

        function writeText(~, path, text)
        %WRITETEXT Overwrite a text file (UTF-8).

            fid = fopen(path, 'w', 'n', 'UTF-8');
            fprintf(fid, '%s', text);
            fclose(fid);
        end

        function writeRawDataFile(~, rawFolder, extension)
        %WRITERAWDATAFILE Create a minimal recognizable raw-data fixture.

            if nargin < 3
                extension = '.bin';
            end
            if ~isfolder(rawFolder)
                mkdir(rawFolder);
            end
            filePath = fullfile(rawFolder, ['raw_data', extension]);
            [fid, message] = fopen(filePath, 'w');
            if fid < 0
                error('TestUMITProjectStore:fixtureWriteFailed', ...
                    'Could not create raw-data fixture: %s', message);
            end
            cleanupFile = onCleanup(@() fclose(fid));
            fwrite(fid, uint8(0), 'uint8');
            clear cleanupFile
        end

        function store = createProject(testCase)
        %CREATEPROJECT Create one valid project for a test.

            uniqueSuffix = char(java.util.UUID.randomUUID());

            projectInfo = struct();
            projectInfo.projectName = [ ...
                'Unit Test Project ', uniqueSuffix(1:8)];
            projectInfo.description = ...
                'Temporary project created by tests.';

            store = UMITProjectStore.create(projectInfo);

            ProjectInfo = store.getProjectInfo();
            testCase.ProjectRoot = store.ProjectRoot;
            testCase.ProjectUUID = ProjectInfo.projectUUID;
        end

        function store = createSecondProject(testCase)
        %CREATESECONDPROJECT Create a second, independent project.

            uniqueSuffix = char(java.util.UUID.randomUUID());

            projectInfo = struct();
            projectInfo.projectName = [ ...
                'Unit Test Project B ', uniqueSuffix(1:8)];
            projectInfo.description = ...
                'Second temporary project created by tests.';

            store = UMITProjectStore.create(projectInfo);

            ProjectInfo = store.getProjectInfo();
            testCase.ProjectRoot2 = store.ProjectRoot;
            testCase.ProjectUUID2 = ProjectInfo.projectUUID;
        end

        function filePath = createMATSource(testCase, fileName, value)
        %CREATEMATSOURCE Create a small valid managed-resource source file.
        %
        %   Filenames containing "reference" produce a canonical
        %   ImageReference payload. Other filenames retain the generic MAT
        %   payload used by transform and calibration tests.

            filePath = fullfile(testCase.SourceFolder, fileName);

            if contains(lower(fileName), 'reference')
                % Identity-less (ImageReference v2): any subject library
                % accepts it (identity fields, when present, must match).
                ImageReference = genImageReferenceStruct( ...
                    single(value .* ones(4, 5)), ...
                    'Name', sprintf('Reference %g', value), ...
                    'SourceType', 'unit-test', ...
                    'CreatedBy', 'TestUMITProjectStore');
                save(filePath, 'ImageReference', '-mat');
            elseif contains(lower(fileName), 'coregistration')
                tform = affine2d([1 0 0; 0 1 0; value 0 1]);
                tformInfo = struct('CandidateSource', 'unit-test');
                save(filePath, 'tform', 'tformInfo', '-mat');
            else
                payload = value;
                save(filePath, 'payload', '-mat');
            end
        end

        function path = resolveRelative(testCase, relativePath)
        %RESOLVERELATIVE Resolve one canonical project-relative path.

            parts = strsplit(relativePath, '/');
            path = fullfile(testCase.ProjectRoot, parts{:});
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
