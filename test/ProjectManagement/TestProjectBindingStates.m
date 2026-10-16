classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestProjectBindingStates < matlab.unittest.TestCase
    %TESTPROJECTBINDINGSTATES Binding states and recovery (project redesign Phase 6).
    %
    %   describeProjectBinding (classification by identifiers and metadata,
    %   never message text; audit item 14), clearing stale links (item 9),
    %   orphan rebind with a new session (item 12), and refusing to take
    %   over a healthy session (item 13). Runs on the temporary UMIT root.

    properties
        Root = ''
        Store = []
    end

    methods (TestMethodSetup)
        function createProject(testCase)
            testCase.Root = tempname;
            mkdir(testCase.Root);
            testCase.addTeardown(@() rmdir(testCase.Root, 's'));
            suffix = char(java.util.UUID.randomUUID());
            testCase.Store = UMITProjectStore.create(struct( ...
                'projectName', ['Binding States ' suffix(1:8)], ...
                'description', 'Temporary project created by tests.'));
            projectRoot = testCase.Store.ProjectRoot;
            testCase.addTeardown(@() iRemoveFolder(projectRoot));
            testCase.Store.addSubject(struct('subjectID', 'Mouse_1'));
        end
    end

    methods (Test)
        function testUnboundBoundInvalidAndPathMismatch(testCase)
            folder = testCase.dataset('A');
            info = describeProjectBinding(folder);
            testCase.verifyEqual(info.state, 'unbound');
            testCase.verifyEqual(info.actions, {'bindToSession'});

            testCase.addSession('S1', folder);
            info = describeProjectBinding(folder);
            testCase.verifyEqual(info.state, 'bound');
            testCase.verifyEqual(info.context.sessionID, 'S1');

            copy = fullfile(testCase.Root, 'A_copy');
            copyfile(folder, copy);
            info = describeProjectBinding(copy);
            testCase.verifyEqual(info.state, 'pathMismatch');
            testCase.verifyEqual(info.identifier, 'Umitoolbox:UMITProjectStore:bindingPathMismatch');

            fid = fopen(UMITProjectLedger.filePath(copy), 'w');
            fprintf(fid, '{"broken": true}');
            fclose(fid);
            info = describeProjectBinding(copy);
            testCase.verifyEqual(info.state, 'invalid');
            testCase.verifyEqual(info.actions, {'restoreLedger'});
        end

        function testOrphanIsDetectedWithoutMessageText(testCase)
            folder = testCase.dataset('A');
            testCase.addSession('S1', folder);
            rmdir(testCase.Store.ProjectRoot, 's');
            info = describeProjectBinding(folder);
            testCase.verifyEqual(info.state, 'orphan');
            testCase.verifyTrue(ismember('removeOrphanBinding', info.actions));
            testCase.verifyTrue(UMITProjectStore.isOrphanProjectBinding(folder));
        end

        function testStaleLinkAfterOfflineRemovalCanBeCleared(testCase)
            % Item 9: a session removed while its folder was offline leaves
            % a link to an available project that no session matches.
            folder = testCase.dataset('A');
            offline = fullfile(testCase.Root, 'A_offline');
            testCase.addSession('S1', folder);
            movefile(folder, offline);
            testCase.Store.removeSessionFromProject('Mouse_1', 'S1');
            movefile(offline, folder);

            info = describeProjectBinding(folder);
            testCase.verifyEqual(info.state, 'stale');
            testCase.verifyTrue(ismember('removeStaleBinding', info.actions));
            testCase.verifyFalse(UMITProjectStore.isOrphanProjectBinding(folder));
            testCase.verifyError(@() UMITProjectStore.removeOrphanProjectBinding(folder), ...
                'Umitoolbox:UMITProjectStore:bindingProjectAvailable');

            UMITProjectStore.removeStaleProjectBinding(folder);
            testCase.verifyEqual(describeProjectBinding(folder).state, 'unbound');
            ledger = UMITProjectLedger.read(folder);
            testCase.verifyEqual(ledger.events{end}.type, 'unbind');
            testCase.verifyEqual(ledger.events{end}.reason, 'staleBinding');

            % The folder can be added again.
            testCase.addSession('S2', folder);
            testCase.verifyEqual(describeProjectBinding(folder).state, 'bound');
        end

        function testStaleLinkAfterOfflineRelocation(testCase)
            folder = testCase.dataset('A');
            newFolder = testCase.dataset('B');
            offline = fullfile(testCase.Root, 'A_offline');
            testCase.addSession('S1', folder);
            movefile(folder, offline);
            testCase.Store.relocateSaveFolder('Mouse_1', 'S1', newFolder);
            movefile(offline, folder);

            testCase.verifyEqual(describeProjectBinding(newFolder).state, 'bound');
            testCase.verifyEqual(describeProjectBinding(folder).state, 'stale');
            testCase.verifyError(@() UMITProjectStore.removeStaleProjectBinding(newFolder), ...
                'Umitoolbox:UMITProjectStore:bindingNotStale');
            UMITProjectStore.removeStaleProjectBinding(folder);
            testCase.verifyEqual(describeProjectBinding(folder).state, 'unbound');
        end

        function testOrphanFolderCanStartANewSession(testCase)
            % Item 12: an orphan link is replaced (a rebind event) when the
            % folder is added as a new session with ReplaceOrphanBinding.
            folder = testCase.dataset('A');
            other = UMITProjectStore.create(struct('projectName', 'Gone', ...
                'description', 'removed below'));
            other.addSubject(struct('subjectID', 'X'));
            other.addSession('X', struct('sessionID', 'Old', 'processedDataFolder', folder));
            rmdir(other.ProjectRoot, 's');
            testCase.verifyEqual(describeProjectBinding(folder).state, 'orphan');

            testCase.verifyError(@() testCase.addSession('S1', folder), ...
                'Umitoolbox:UMITProjectStore:bindingFileExists');
            testCase.Store.addSession('Mouse_1', struct('sessionID', 'S1', ...
                'processedDataFolder', folder), 'ReplaceOrphanBinding', true);
            info = describeProjectBinding(folder);
            testCase.verifyEqual(info.state, 'bound');
            ledger = UMITProjectLedger.read(folder);
            testCase.verifyEqual(cellfun(@(e) e.type, ledger.events, 'UniformOutput', false), ...
                {'bind', 'rebind'});
        end

        function testRelocationOntoAnOrphanDoesNotTakeOverAHealthySession(testCase)
            % Item 13: a session that is healthily bound elsewhere is not
            % moved onto an orphan folder.
            healthy = testCase.dataset('Healthy');
            testCase.addSession('S1', healthy);
            orphan = testCase.dataset('Orphan');
            other = UMITProjectStore.create(struct('projectName', 'Gone', ...
                'description', 'removed below'));
            other.addSubject(struct('subjectID', 'X'));
            other.addSession('X', struct('sessionID', 'Old', 'processedDataFolder', orphan));
            rmdir(other.ProjectRoot, 's');

            testCase.verifyError(@() testCase.Store.relocateSaveFolder('Mouse_1', 'S1', ...
                orphan, 'ReplaceOrphanBinding', true), ...
                'Umitoolbox:UMITProjectStore:sessionAlreadyBound');
            testCase.verifyEqual(describeProjectBinding(healthy).state, 'bound');
            testCase.verifyEqual(describeProjectBinding(orphan).state, 'orphan');

            % Once the healthy folder is gone (offline), the takeover is a
            % legitimate relocation.
            movefile(healthy, [healthy '_offline']);
            testCase.Store.relocateSaveFolder('Mouse_1', 'S1', orphan, ...
                'ReplaceOrphanBinding', true);
            testCase.verifyEqual(describeProjectBinding(orphan).state, 'bound');
        end
    end

    methods (Test)
        function testBindingWithNewProjectIsAllOrNothing(testCase)
            % Item 11: a failure removes what the call created, and a retry
            % creates exactly one project.
            folder = fullfile(testCase.Root, 'NotYetImported');
            mkdir(folder);
            projectName = ['Bind Retry ' char(java.util.UUID.randomUUID())];
            selection = struct('project', struct('projectName', projectName, ...
                'description', 'test'), 'subject', struct('subjectID', 'M1'), ...
                'session', struct('sessionID', 'S1'));

            testCase.verifyError(@() bindFolderToProjectSession(folder, selection), ...
                'Umitoolbox:UMITProjectStore:bindFolderFailed');
            testCase.verifyEqual(testCase.countProjects(projectName), 0);
            testCase.verifyEqual(describeProjectBinding(folder).state, 'unbound');

            AcqInfoStream = struct('Width', 4, 'Height', 4);
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            [context, store] = bindFolderToProjectSession(folder, selection);
            testCase.addTeardown(@() iRemoveFolder(store.ProjectRoot));
            testCase.verifyEqual(context.sessionID, 'S1');
            testCase.verifyEqual(testCase.countProjects(projectName), 1);
        end

        function testBindingWithNewSubjectRemovesItOnFailure(testCase)
            folder = fullfile(testCase.Root, 'Empty');
            mkdir(folder);
            ProjectInfo = testCase.Store.getProjectInfo();
            selection = struct('project', ProjectInfo.projectUUID, ...
                'subject', struct('subjectID', 'Mouse_2'), ...
                'session', struct('sessionID', 'S1'));
            testCase.verifyError(@() bindFolderToProjectSession(folder, selection), ...
                'Umitoolbox:UMITProjectStore:bindFolderFailed');
            registry = testCase.Store.getProjectInfo().subjectRegistry;
            testCase.verifyEqual({registry.id}, {'Mouse_1'});

            % Existing project, subject and session.
            existing = testCase.dataset('Existing');
            testCase.addSession('S9', existing);
            moved = testCase.dataset('Moved');
            movefile(existing, [existing '_gone']);
            context = bindFolderToProjectSession(moved, struct('project', ...
                ProjectInfo.projectUUID, 'subject', 'Mouse_1', 'session', 'S9'));
            testCase.verifyEqual(context.sessionID, 'S9');
        end
    end

    methods (Test)
        function testRawFolderFollowsDataParams(testCase)
            % Item 7: DataParams.folders.RawFolder is the single RawFolder
            % authority; the session's raw binding follows it.
            saveFolder = testCase.dataset('Save');
            rawA = testCase.rawFolder('RawA');
            rawB = testCase.rawFolder('RawB');
            saveDataParams(saveFolder, createDataParams(saveFolder));
            updateDataParam(saveFolder, 'folders.RawFolder', rawA);
            testCase.addSession('S1', saveFolder);

            testCase.verifyError(@() testCase.Store.bindRawDataFolder('Mouse_1', 'S1', rawB), ...
                'Umitoolbox:UMITProjectStore:rawFolderMismatch');
            testCase.Store.bindRawDataFolder('Mouse_1', 'S1', '');
            SessionInfo = testCase.Store.getSessionInfo('Mouse_1', 'S1');
            testCase.verifyEqual(lower(SessionInfo.rawDataFolder), lower(rawA));
            testCase.verifyEqual(describeProjectBinding(rawA).state, 'bound');

            report = testCase.Store.validate('Mode', 'quick');
            testCase.verifyFalse(any(strcmp({report.warnings.code}, 'raw_folder_mismatch')));
            updateDataParam(saveFolder, 'folders.RawFolder', rawB);
            report = testCase.Store.validate('Mode', 'quick');
            testCase.verifyTrue(report.isValid, 'a mismatch is a warning, not an error');
            testCase.verifyTrue(any(strcmp({report.warnings.code}, 'raw_folder_mismatch')));

            testCase.Store.relocateRawDataFolder('Mouse_1', 'S1', rawB);
            testCase.verifyEqual(describeProjectBinding(rawB).state, 'bound');
        end
    end

    methods (Access = private)
        function folder = rawFolder(testCase, name)
            folder = fullfile(testCase.Root, name);
            mkdir(folder);
            fid = fopen(fullfile(folder, 'img_00000.bin'), 'w');
            fwrite(fid, uint8(1:10));
            fclose(fid);
        end

        function count = countProjects(~, projectName)
            projects = UMITProjectStore.listProjects();
            count = sum(strcmp(cellstr(projects.ProjectName), projectName));
        end

        function folder = dataset(testCase, name)
            folder = fullfile(testCase.Root, name);
            mkdir(folder);
            AcqInfoStream = struct('Width', 4, 'Height', 4);
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
        end

        function addSession(testCase, sessionID, folder)
            testCase.Store.addSession('Mouse_1', struct( ...
                'sessionID', sessionID, 'processedDataFolder', folder));
        end
    end
end

function iRemoveFolder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
