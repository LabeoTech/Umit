classdef TestUMITProjectStore < matlab.unittest.TestCase
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
        %   teardown runs against the real, non-redirectable Rig store (see
        %   createTemporaryFolders), so one locked/failed removal must never
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
        %created on the real Rig store.
        %
        %   If the fixture Rig is still the store's default, repoint the
        %   default to another existing Rig first. Deleting a default Rig's
        %   folder without repointing leaves UMITRigStore.getDefaultRig /
        %   getOrCreateDefaultRig permanently broken for every later caller,
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

        function testNormalizeSubjectAndSessionIDsToValidMATLABNames(testCase)
        %TESTNORMALIZESUBJECTANDSESSIONIDSTOVALIDMATLABNAMES Sanitize IDs.

            store = testCase.createProject();
            subjectUUID = store.addSubject(struct('subjectID', 'Mouse 104'));
            SubjectInfo = store.getSubjectInfo('Mouse 104');
            testCase.verifyEqual(SubjectInfo.uuid, subjectUUID);
            testCase.verifyEqual(SubjectInfo.subjectID, 'Mouse104');

            sessionUUID = testCase.addSession(store, 'Mouse 104', struct( ...
                'sessionID', '2026-07-10 Session-01'));
            SessionInfo = store.getSessionInfo( ...
                'Mouse104', '2026-07-10 Session-01');
            testCase.verifyEqual(SessionInfo.uuid, sessionUUID);
            testCase.verifyEqual(SessionInfo.sessionID, ...
                'x2026_07_10Session_01');

            testCase.verifyError(@() store.addSubject( ...
                struct('subjectID', 'CON')), ...
                'Umitoolbox:UMITProjectStore:invalidID');

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

        function testSessionRigPointerIsUUIDBackedAcrossExternalRigRename(testCase)
        %TESTSESSIONRIGPOINTERISUUIDBACKEDACROSSEXTERNALRIGRENAME Verify a
        %session's rig pointer survives an external rig rename.
        %
        %   Rigs are independent of any project (UMITRigStore), so
        %   UMITProjectStore has no way to cascade-update a session's cached
        %   rigID when the rig is renamed externally -- only the immutable
        %   rigUUID is guaranteed current. This is an accepted limitation:
        %   the UUID remains authoritative and can always re-resolve the
        %   rig's current human-readable ID via UMITRigStore.

            rigStore = UMITRigStore.create(struct('rigID', 'WidefieldRig01'));
            testCase.RigRoot = rigStore.RigRoot;
            rigUUID = rigStore.getRigInfo().uuid;

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', '2026-07-10_S01', ...
                'rigID', 'WidefieldRig01'));

            SessionInfoBefore = store.getSessionInfo( ...
                'Mouse_104', '2026-07-10_S01');
            testCase.verifyEqual(SessionInfoBefore.rigUUID, rigUUID);
            testCase.verifyEqual(SessionInfoBefore.rigID, 'WidefieldRig01');

            rigStore.renameRigID('WidefieldRig02');
            % renameRigID updates rigStore.RigRoot in place; refresh the
            % tracked cleanup path so teardown removes the folder at its
            % actual (post-rename) location, not the stale pre-rename path.
            testCase.RigRoot = rigStore.RigRoot;

            SessionInfoAfter = store.getSessionInfo( ...
                'Mouse_104', '2026-07-10_S01');
            testCase.verifyEqual(SessionInfoAfter.rigUUID, rigUUID);
            testCase.verifyEqual(SessionInfoAfter.rigID, 'WidefieldRig01', ...
                'Cached rigID is not cascade-updated by an external rig rename; only rigUUID is guaranteed current.');

            reresolvedRig = UMITRigStore.open(SessionInfoAfter.rigUUID);
            testCase.verifyEqual(reresolvedRig.getRigInfo().rigID, 'WidefieldRig02');

            report = store.validate('Mode', 'full');
            testCase.verifyTrue(report.isValid, ...
                testCase.issueMessages(report.errors));
        end

        function testNewSessionRequiresAssignableRig(testCase)
        %TESTNEWSESSIONREQUIRESASSIGNABLERIG Auto-resolve the Active Rig when
        %the caller supplies none, and reject caller-supplied Rig fields
        %outright -- UMITRigStore is the sole owner of Rig assignment, so
        %addSession never accepts or validates a caller-asserted Rig anymore.
        %(The Archived-Rig-rejection scenario previously tested here now
        %belongs to UMITRigStore.assignDatasetRig -- see
        %test/RigManagement/TestUMITRigStore.m,
        %testAssignDatasetRigRejectsArchivedRig.)

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            saveFolderA = fullfile(testCase.TempRoot, 'RigRequiredSaveFolderA');
            mkdir(saveFolderA);
            testCase.writeValidAcqInfos(saveFolderA);

            activeRigInfo = UMITRigStore.getActiveRig().getRigInfo();

            store.addSession('Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolderA));
            SessionInfo = store.getSessionInfo('Mouse_104', 'Session_01');
            testCase.verifyEqual(SessionInfo.rigUUID, activeRigInfo.uuid);
            testCase.verifyEqual(SessionInfo.rigID, activeRigInfo.rigID);

            saveFolderB = fullfile(testCase.TempRoot, 'RigRequiredSaveFolderB');
            mkdir(saveFolderB);
            testCase.writeValidAcqInfos(saveFolderB);

            testCase.verifyError(@() store.addSession('Mouse_104', struct( ...
                'sessionID', 'Session_02', ...
                'rigID', activeRigInfo.rigID, ...
                'processedDataFolder', saveFolderB)), ...
                'Umitoolbox:UMITProjectStore:rigAssignmentNotSupported');
        end

        function testArchivedHistoricalSessionReferenceRemainsValid(testCase)
        %TESTARCHIVEDHISTORICALSESSIONREFERENCEREMAINSVALID Preserve history.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            historical = UMITRigStore.getActiveRig();
            replacement = UMITRigStore.create(struct('rigID', ...
                ['ReplacementHistory_' strrep(char(java.util.UUID.randomUUID()), '-', '')]));
            testCase.RigRoot2 = replacement.RigRoot;
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'rigID', historical.getRigInfo().rigID));

            historical.archiveRig(replacement.getRigInfo().uuid);
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(SessionInfo.rigUUID, ...
                historical.getRigInfo().uuid);
            testCase.verifyTrue(UMITRigStore.rigExists(SessionInfo.rigUUID));
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

        function testUnregisteredFolderOpensReadOnly(testCase)
        %TESTUNREGISTEREDFOLDEROPENSREADONLY Detect direct folder edits.

            store = testCase.createProject(); %#ok<NASGU>
            clear store

            mkdir(fullfile(testCase.ProjectRoot, ...
                'subjects', 'ManuallyCreatedMouse'));

            reopened = UMITProjectStore.open(testCase.ProjectUUID);
            testCase.verifyTrue(reopened.IsReadOnly);
            testCase.verifyFalse(reopened.LastValidationReport.isValid);
            testCase.verifyTrue(any(strcmp( ...
                {reopened.LastValidationReport.errors.code}, ...
                'unregistered_folder')));
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

        function testListSessionResources(testCase)
        %TESTLISTSESSIONRESOURCES Verify session-owned resource queries.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct('sessionID', '2026-07-10_S01'));

            transformSource = testCase.createMATSource( ...
                'registration.mat', 1);
            transformUUID = store.addRegistrationTransform( ...
                'Mouse_104', '2026-07-10_S01', ...
                transformSource, struct());

            sessionResources = store.listSessionResources( ...
                'Mouse_104', '2026-07-10_S01');
            testCase.verifyEqual(numel(sessionResources), 1);
            testCase.verifyEqual(sessionResources.uuid, transformUUID);
            testCase.verifyEqual(sessionResources.ownerType, 'session');
        end

        function testRigResourcesLiveInUMITRigStoreNotProject(testCase)
        %TESTRIGRESOURCESLIVEINUMITRIGSTORENOTPROJECT Verify rig resources
        %(cameraCoregistration) are owned by the
        %independent UMITRigStore, not by UMITProjectStore, and that a
        %session merely keeps a rigID/rigUUID pointer to that external rig.

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
            testCase.verifyEqual(SessionInfo.rigID, 'WidefieldRig01');
            testCase.verifyEqual(SessionInfo.rigUUID, rigStore.getRigInfo().uuid);
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
            testCase.verifyEmpty( ...
                store.getActiveRegistrationTransform( ...
                'Mouse_104', '2026-07-10_S01'));

            referenceUUID = store.addImageReference( ...
                'Mouse_104', ...
                testCase.createMATSource('reference.mat', 1), struct());
            transformUUID = store.addRegistrationTransform( ...
                'Mouse_104', '2026-07-10_S01', ...
                testCase.createMATSource('registration.mat', 2), struct());

            activeReference = ...
                store.getActiveImageReference('Mouse_104');
            activeTransform = store.getActiveRegistrationTransform( ...
                'Mouse_104', '2026-07-10_S01');

            testCase.verifyEqual(activeReference.uuid, referenceUUID);
            testCase.verifyEqual(activeTransform.uuid, transformUUID);
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
        %TESTFINDSESSIONBYRAWANDPROCESSEDFOLDER Resolve .umitlink context.

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
            testCase.verifyEqual(rawContext.rigID, 'WidefieldRig01');
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
        %TESTPROJECTBINDINGFILEANDSESSIONRECIPROCITY Verify .umitlink creation.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            ProjectBinding = store.bindProcessedDataFolder( ...
                'Mouse_104', 'Session_01', saveFolder);

            schema = getUMITProjectSchema();
            bindingPath = fullfile(saveFolder, ...
                schema.files.projectBinding);
            testCase.verifyTrue(isfile(bindingPath));

            loaded = load(bindingPath, ...
                schema.metadataVariables.projectBinding, '-mat');
            testCase.verifyTrue(isfield(loaded, ...
                schema.metadataVariables.projectBinding));
            testCase.verifyEqual(loaded.ProjectBinding.bindingUUID, ...
                ProjectBinding.bindingUUID);
            testCase.verifyEqual(loaded.ProjectBinding.folderRole, ...
                'processedDataFolder');

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
            testCase.verifyTrue(isfile(fullfile(saveFolder, ...
                UMITProjectStore.getProjectBindingFileName())));

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
                'Umitoolbox:UMITProjectStore:missingAcqInfos');

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
            testCase.verifyFalse(isfile(fullfile(saveFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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

        function testMigrateSessionRigAssignmentRepairsLegacySession(testCase)
        %TESTMIGRATESESSIONRIGASSIGNMENTREPAIRSLEGACYSESSION DFR-20260818-002:
        %a session with blank rigUUID/rigID (pre-dating Rig ownership) is
        %exactly the case migrateSessionRigAssignment exists to repair, and
        %must actually succeed against it -- not be rejected by the same
        %full-validation gate it is trying to fix.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            sessionMatPath = fullfile(testCase.ProjectRoot, 'subjects', ...
                'Mouse_104', 'sessions', 'Session_01', 'session.mat');
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            SessionInfo.rigUUID = '';
            SessionInfo.rigID = '';
            save(sessionMatPath, 'SessionInfo', '-mat');

            preRepairReport = store.validate('Mode', 'full');
            testCase.verifyFalse(preRepairReport.isValid);
            testCase.verifyTrue(any(strcmp( ...
                {preRepairReport.errors.code}, 'session_rig_required')));

            activeRigUUID = UMITRigStore.getActiveRig().getRigInfo().uuid;
            repaired = store.migrateSessionRigAssignment('Mouse_104', 'Session_01');
            testCase.verifyEqual(repaired.rigUUID, activeRigUUID);
            testCase.verifyNotEmpty(repaired.rigID);

            postRepairReport = store.validate('Mode', 'full');
            testCase.verifyTrue(postRepairReport.isValid, ...
                testCase.issueMessages(postRepairReport.errors));
        end

        function testMigrateSessionRigAssignmentBlockedByUnrelatedError(testCase)
        %TESTMIGRATESESSIONRIGASSIGNMENTBLOCKEDBYUNRELATEDERROR
        %The narrowed mutation-health gate exempts only the specific
        %pre-diagnosed rig issue on the targeted session -- an unrelated
        %error on that same session must still block the repair.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));

            sessionUUID = store.getSessionInfo('Mouse_104', 'Session_01').uuid;
            sessionMatPath = fullfile(testCase.ProjectRoot, 'subjects', ...
                'Mouse_104', 'sessions', 'Session_01', 'session.mat');
            SessionInfo = store.getSessionInfoByUUID(sessionUUID);
            SessionInfo.rigUUID = '';
            SessionInfo.rigID = '';
            SessionInfo.sessionID = 'Session_01_Renamed';
            save(sessionMatPath, 'SessionInfo', '-mat');

            testCase.verifyError(@() store.migrateSessionRigAssignment( ...
                'Mouse_104', 'Session_01'), ...
                'Umitoolbox:UMITProjectStore:invalidProject');
        end

        function testUpdateSessionMetadataRejectsRigFields(testCase)
        %TESTUPDATESESSIONMETADATAREJECTSRIGFIELDS Rig assignment is owned
        %exclusively by UMITRigStore -- updateSessionMetadata must never
        %accept rigID/rigUUID, even alongside otherwise-valid updates.

            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            sessionUUID = testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01'));
            originalRigUUID = store.getSessionInfoByUUID(sessionUUID).rigUUID;

            replacementRig = UMITRigStore.create(struct('rigID', ...
                ['RejectedRigUpdate_' strrep(char(java.util.UUID.randomUUID()), '-', '')]));
            testCase.RigRoot2 = replacementRig.RigRoot;

            testCase.verifyError(@() store.updateSessionMetadata( ...
                'Mouse_104', 'Session_01', struct( ...
                'rigUUID', replacementRig.getRigInfo().uuid)), ...
                'Umitoolbox:UMITProjectStore:rigAssignmentNotSupported');
            testCase.verifyError(@() store.updateSessionMetadata( ...
                'Mouse_104', 'Session_01', struct( ...
                'notes', 'Unrelated update', ...
                'rigID', replacementRig.getRigInfo().rigID)), ...
                'Umitoolbox:UMITProjectStore:rigAssignmentNotSupported');

            testCase.verifyEqual( ...
                store.getSessionInfoByUUID(sessionUUID).rigUUID, originalRigUUID);
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
            testCase.verifyFalse(isfile(fullfile(oldFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyTrue(isfile(fullfile(newFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
                'Umitoolbox:UMITProjectStore:missingAcqInfos');

            after = store.getSessionInfoByUUID(sessionUUID);
            testCase.verifyEqual(after.processedDataFolder, ...
                before.processedDataFolder);
            testCase.verifyEqual(after.processedDataBindingUUID, ...
                before.processedDataBindingUUID);
            testCase.verifyTrue(isfile(fullfile(oldFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyFalse(isfile(fullfile(invalidFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
                {report.errors.code}, 'missing_acq_infos')));
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
            testCase.verifyTrue(isfile(fullfile(newFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyFalse(isfile(fullfile(oldFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyFalse(isfile(fullfile(oldRawFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyTrue(isfile(fullfile(newRawFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyTrue(isfile(fullfile(rawFolderOne, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyTrue(isfile(fullfile(oldSharedFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyTrue(isfile(fullfile(newSaveFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyFalse(isfile(fullfile(saveFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyTrue(isfile(fullfile(parkedFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyFalse(isfile(fullfile(saveFolder1, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyFalse(isfile(fullfile(saveFolder2, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyTrue(isfile(fullfile(offlineFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyFalse(isfile(fullfile(saveFolder1, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyFalse(isfile(fullfile(saveFolder2, ...
                UMITProjectStore.getProjectBindingFileName())));
            testCase.verifyEqual(otherStore.getProjectInfo().projectUUID, ...
                testCase.ProjectUUID2);
        end

        function testDeleteProjectFailureLeavesProjectRecoverable(testCase)
        %TESTDELETEPROJECTFAILURELEAVESPROJECTRECOVERABLE Reject unavailable link.

            saveFolder = fullfile(testCase.TempRoot, 'SaveFolder');
            offlineFolder = fullfile(testCase.TempRoot, 'OfflineSaveFolder');
            mkdir(saveFolder);
            testCase.writeValidAcqInfos(saveFolder);
            store = testCase.createProject();
            store.addSubject(struct('subjectID', 'Mouse_104'));
            testCase.addSession(store, 'Mouse_104', struct( ...
                'sessionID', 'Session_01', ...
                'processedDataFolder', saveFolder));
            movefile(saveFolder, offlineFolder);

            testCase.verifyError(@() store.deleteProject(testCase.ProjectUUID), ...
                'Umitoolbox:UMITProjectStore:bindingFolderUnavailable');
            testCase.verifyTrue(isfolder(testCase.ProjectRoot));
            testCase.verifyTrue(isfile(fullfile(offlineFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyFalse(isfile(fullfile(saveFolder, ...
                UMITProjectStore.getProjectBindingFileName())));
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
            testCase.verifyTrue(isfile(fullfile(saveFolder, ...
                UMITProjectStore.getProjectBindingFileName())));

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
        %   Rig assignment is owned exclusively by UMITRigStore -- addSession
        %   itself now rejects caller-supplied rigID/rigUUID and always
        %   auto-resolves via UMITRigStore.ensureDatasetRigAssociation. A
        %   caller of this fixture that wants a SPECIFIC (non-Active) Rig
        %   pinned assigns it to the SaveFolder first via
        %   UMITRigStore.assignDatasetRig, then lets auto-resolution pick up
        %   that existing association.

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
                    targetRigUUID = sessionInfo.rigUUID;
                else
                    targetRigUUID = UMITRigStore.openByRigID( ...
                        sessionInfo.rigID).getRigInfo().uuid;
                end
                UMITRigStore.assignDatasetRig( ...
                    sessionInfo.processedDataFolder, targetRigUUID);
                sessionInfo = rmfield(sessionInfo, requestedRigFields);
            end

            sessionUUID = store.addSession(subjectID, sessionInfo);
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
                ImageReference = genImageReferenceStruct( ...
                    single(value .* ones(4, 5)), ...
                    'Name', sprintf('Reference %g', value), ...
                    'ProjectUUID', lower(char(java.util.UUID.randomUUID())), ...
                    'ProjectName', 'Unit Test Project', ...
                    'SubjectUUID', lower(char(java.util.UUID.randomUUID())), ...
                    'SubjectID', 'Mouse_104', ...
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
