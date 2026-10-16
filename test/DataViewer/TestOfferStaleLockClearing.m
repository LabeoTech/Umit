classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestOfferStaleLockClearing < matlab.unittest.TestCase
    %TESTOFFERSTALELOCKCLEARING Clearing a stale lock from a "locked" error.
    %
    %   offerStaleLockClearing recognises the three lock errors (also as a
    %   cause), clears only stale locks, keeps recent ones, and ignores any
    %   other error. Run non-interactively; the dialogs are checked by hand.

    properties
        Folder = ''
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            testCase.Folder = tempname;
            mkdir(testCase.Folder);
            testCase.addTeardown(@() rmdir(testCase.Folder, 's'));
        end
    end

    methods (Test)
        function testOtherErrorsAreIgnored(testCase)
            [cleared, message, isLockError] = offerStaleLockClearing([], ...
                MException('Some:other:error', 'x'), 'Folder', testCase.Folder, ...
                'Interactive', false);
            testCase.verifyFalse(cleared);
            testCase.verifyEmpty(message);
            testCase.verifyFalse(isLockError, 'the caller shows its own alert');
        end

        function testStaleProjectLedgerLockIsCleared(testCase)
            lockFolder = testCase.makeLock(UMITProjectLedger.LockFolderName, hours(2));
            [cleared, message] = offerStaleLockClearing([], ...
                MException('Umitoolbox:ProjectLedger:ledgerLocked', 'locked'), ...
                'Folder', testCase.Folder, 'Interactive', false);
            testCase.verifyTrue(cleared);
            testCase.verifyFalse(isfolder(lockFolder));
            testCase.verifySubstring(message, 'Try the operation again');
        end

        function testRecentLockIsKept(testCase)
            lockFolder = testCase.makeLock(UMITProjectLedger.LockFolderName, minutes(1));
            [cleared, message, isLockError] = offerStaleLockClearing([], ...
                MException('Umitoolbox:ProjectLedger:ledgerLocked', 'locked'), ...
                'Folder', testCase.Folder, 'Interactive', false);
            testCase.verifyFalse(cleared);
            testCase.verifyTrue(isLockError);
            testCase.verifyTrue(isfolder(lockFolder));
            testCase.verifySubstring(message, 'recent');
        end

        function testLockErrorAsACauseOfARigLedgerWrite(testCase)
            lockFolder = testCase.makeLock(UMITRigLedger.LockFolderName, hours(1));
            wrapper = MException('Umitoolbox:UMITProjectStore:bindFolderFailed', 'failed');
            wrapper = addCause(wrapper, MException('Umitoolbox:RigLedger:ledgerLocked', 'locked'));
            cleared = offerStaleLockClearing([], wrapper, 'Folder', testCase.Folder, ...
                'Interactive', false);
            testCase.verifyTrue(cleared);
            testCase.verifyFalse(isfolder(lockFolder));
        end

        function testFolderIsTakenFromTheLockError(testCase)
            % A real ledger "locked" error names its folder; no 'Folder'
            % option is needed (the Project Manager does not know it).
            lockFolder = testCase.makeLock(UMITProjectLedger.LockFolderName, hours(1));
            fclose(fopen(UMITProjectLedger.filePath(testCase.Folder), 'w'));
            lockedError = [];
            try
                UMITProjectLedger.recordFreshImport(testCase.Folder, 'Origin', 'test');
            catch lockedError
            end
            testCase.assertNotEmpty(lockedError);
            testCase.assertEqual(lockedError.identifier, 'Umitoolbox:ProjectLedger:ledgerLocked');
            wrapper = addCause(MException('Some:write:failed', 'failed'), lockedError);

            [cleared, message] = offerStaleLockClearing([], wrapper, ...
                'Folder', tempdir, 'Interactive', false);
            testCase.verifyTrue(cleared, 'the folder named in the error wins');
            testCase.verifyFalse(isfolder(lockFolder));
            testCase.verifySubstring(message, testCase.Folder);
        end

        function testStaleProjectLockIsCleared(testCase)
            suffix = char(java.util.UUID.randomUUID());
            store = UMITProjectStore.create(struct('projectName', ['Locks ' suffix(1:8)], ...
                'description', 'test'));
            projectRoot = store.ProjectRoot;
            testCase.addTeardown(@() rmdir(projectRoot, 's'));
            lockFolder = fullfile(projectRoot, '.umit', 'project.lock');
            mkdir(lockFolder);
            LockInfo = struct('operation', 'addSession', ...
                'createdOn', datetime('now') - hours(3), 'processID', 1, ...
                'userName', 'someone', 'hostName', 'somewhere');
            save(fullfile(lockFolder, 'lockInfo.mat'), 'LockInfo');

            cleared = offerStaleLockClearing([], ...
                MException('Umitoolbox:UMITProjectStore:projectLocked', 'locked'), ...
                'Store', store, 'Interactive', false);
            testCase.verifyTrue(cleared);
            testCase.verifyFalse(isfolder(lockFolder));
            store.addSubject(struct('subjectID', 'M1')); % writable again
        end
    end

    methods (Access = private)
        function lockFolder = makeLock(testCase, name, age)
            lockFolder = fullfile(testCase.Folder, name);
            mkdir(lockFolder);
            LockInfo = struct('operation', 'test', 'createdOn', datetime('now') - age, ...
                'processID', 1);
            save(fullfile(lockFolder, 'lockInfo.mat'), 'LockInfo');
        end
    end
end
