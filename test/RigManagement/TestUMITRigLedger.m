classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestUMITRigLedger < matlab.unittest.TestCase
    %TESTUMITRIGLEDGER Phase-1 tests of the folder-bound rig ledger.
    %
    %   Covers state transitions, canonical hashing, filter-set rules,
    %   fresh-import reset, conflicts, acknowledgement, hadConflict
    %   permanence, integrity failures and atomic-write behaviour. Rigs are
    %   plain RigInfo structs over built-in spectra, so the tests never
    %   touch the per-machine rig library.

    properties
        Folder = ''
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.Folder = tempname;
            mkdir(testCase.Folder);
        end
    end

    methods (TestMethodTeardown)
        function removeFolder(testCase)
            if isfolder(testCase.Folder)
                rmdir(testCase.Folder, 's');
            end
        end
    end

    %% State transitions
    methods (Test)
        function testMissingLedgerReadsUnknownWithoutWriting(testCase)
            info = UMITRigLedger.read(testCase.Folder);
            report = UMITRigLedger.validate(testCase.Folder);

            testCase.verifyEqual(info.state, 'Unknown');
            testCase.verifyFalse(info.hadConflict);
            testCase.verifyTrue(report.isValid);
            testCase.verifyFalse(UMITRigLedger.exists(testCase.Folder));
            testCase.verifyEmpty(testCase.folderEntries());
        end

        function testNoRigIsIdempotentAndNeverCommits(testCase)
            testCase.verifyTrue(UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test'));
            testCase.verifyFalse(UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test'));

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'NoRig');
            testCase.verifyNumElements(info.events, 1);
            testCase.verifyError(@() UMITRigLedger.commit(testCase.Folder, ...
                'filterSet', UMITRigLedger.buildFilterSetSnapshot('none'), 'Origin', 'test'), ...
                'Umitoolbox:RigLedger:noRigBound');
            testCase.verifyError(@() UMITRigLedger.assertCanUseOptics(testCase.Folder), ...
                'Umitoolbox:RigLedger:noRigBound');
        end

        function testBindFromNoRigThenNoRigIsRefused(testCase)
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'BoundUncommitted');
            testCase.verifyEqual(info.binding.sourceRig.rigID, 'A');
            testCase.verifyError(@() UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test'), ...
                'Umitoolbox:RigLedger:alreadyBound');
            testCase.verifyError(@() UMITRigLedger.bind(testCase.Folder, ...
                testCase.binding('B'), 'Origin', 'test'), 'Umitoolbox:RigLedger:alreadyBound');
        end

        function testRebindWhileUncommittedIsFree(testCase)
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');
            conflict = UMITRigLedger.rebind(testCase.Folder, ...
                testCase.binding('B', 'LED_593nm', 2), 'Origin', 'test');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyFalse(conflict.hasConflict);
            testCase.verifyEqual(info.state, 'BoundUncommitted');
            testCase.verifyEqual(info.binding.sourceRig.rigID, 'B');
            testCase.verifyFalse(info.hadConflict);
        end

        function testRebindToIdenticalBindingWritesNothing(testCase)
            binding = testCase.binding('A');
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');
            before = fileread(UMITRigLedger.filePath(testCase.Folder));

            UMITRigLedger.rebind(testCase.Folder, binding, 'Origin', 'test');

            testCase.verifyEqual(fileread(UMITRigLedger.filePath(testCase.Folder)), before);
        end

        function testCommitMovesToBoundCommittedOnce(testCase)
            binding = testCase.binding('A');
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');

            testCase.verifyTrue(UMITRigLedger.commit(testCase.Folder, 'optics', ...
                binding.optics.hash, 'Origin', 'test', 'FirstUsedBy', 'HemoCompute'));
            testCase.verifyFalse(UMITRigLedger.commit(testCase.Folder, 'optics', ...
                binding.optics, 'Origin', 'test'));

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'BoundCommitted');
            testCase.verifyEqual(info.commitments.optics, binding.optics.hash);
            testCase.verifyEqual(info.events{end}.firstUsedBy, 'HemoCompute');
            other = testCase.binding('B', 'LED_593nm');
            testCase.verifyError(@() UMITRigLedger.commit(testCase.Folder, 'optics', ...
                other.optics.hash, 'Origin', 'test'), 'Umitoolbox:RigLedger:commitMismatch');
        end

        function testEventsCarryNoUserOrHost(testCase)
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'unitTest');
            info = UMITRigLedger.read(testCase.Folder);
            event = info.events{1};

            testCase.verifyEqual(event.origin, 'unitTest');
            testCase.verifyEqual(event.toolboxVersion, UMITRigLedger.ToolboxVersion);
            testCase.verifyFalse(any(isfield(event, {'user', 'host'})));
            testCase.verifyError(@() UMITRigLedger.writeNoRig(testCase.Folder), ...
                'Umitoolbox:RigLedger:invalidInput');
        end
    end

    %% Conflicts, acknowledgement, hadConflict
    methods (Test)
        function testRebindOfCommittedOpticsConflictsAndBlocks(testCase)
            testCase.bindAndCommitOptics('A');
            preview = UMITRigLedger.previewRebind(testCase.Folder, ...
                testCase.binding('B', 'LED_593nm'));
            conflict = UMITRigLedger.rebind(testCase.Folder, ...
                testCase.binding('B', 'LED_593nm'), 'Origin', 'test');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyTrue(preview.hasConflict);
            testCase.verifyEqual(conflict.parts, {'optics'});
            testCase.verifyEqual(info.state, 'Conflicted');
            testCase.verifyTrue(info.hadConflict);
            testCase.verifyEqual(info.events{end}.type, 'conflict');

            filterSet = UMITRigLedger.buildFilterSetSnapshot('none');
            testCase.verifyError(@() UMITRigLedger.assertCanUseOptics(testCase.Folder), ...
                'Umitoolbox:RigLedger:conflicted');
            testCase.verifyError(@() UMITRigLedger.assertCanUseFilterSet( ...
                testCase.Folder, filterSet), 'Umitoolbox:RigLedger:conflicted');
            testCase.verifyError(@() UMITRigLedger.commit(testCase.Folder, ...
                'filterSet', filterSet, 'Origin', 'test'), 'Umitoolbox:RigLedger:conflicted');
        end

        function testConflictMessageNamesRigsAndFixesWithoutHashes(testCase)
            % Bruno (2026-10-15): users are told which rigs differ, never
            % snapshot hashes.
            testCase.bindAndCommitOptics('A');
            requested = testCase.binding('B', 'LED_593nm');
            UMITRigLedger.rebind(testCase.Folder, requested, 'Origin', 'test');

            try
                UMITRigLedger.assertCanUseOptics(testCase.Folder, 'Operation', 'HemoCompute');
                testCase.verifyFail('Expected a conflict error.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:RigLedger:conflicted');
                testCase.verifySubstring(ME.message, 'HemoCompute was refused');
                testCase.verifySubstring(ME.message, testCase.Folder);
                testCase.verifySubstring(ME.message, 'Rig A');
                testCase.verifySubstring(ME.message, 'Rig B');
                testCase.verifySubstring(ME.message, 'whose optics differ');
                testCase.verifyFalse(contains(ME.message, UMITRigLedger.HashPrefix));
                testCase.verifySubstring(ME.message, 'revertBinding');
                testCase.verifySubstring(ME.message, 'fresh save folder');
                testCase.verifySubstring(ME.message, 'acknowledgeConflict');
            end
        end

        function testRevertResolvesConflictButHadConflictStays(testCase)
            committed = testCase.bindAndCommitOptics('A');
            UMITRigLedger.rebind(testCase.Folder, testCase.binding('B', 'LED_593nm'), 'Origin', 'test');

            conflict = UMITRigLedger.revertBinding(testCase.Folder, 'Origin', 'test');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyFalse(conflict.hasConflict);
            testCase.verifyEqual(info.state, 'BoundCommitted');
            testCase.verifyEqual(info.binding.opticsHash, committed.optics.hash);
            testCase.verifyEqual(info.binding.sourceRig.rigID, 'A');
            testCase.verifyTrue(info.hadConflict);
            result = UMITRigLedger.assertCanUseOptics(testCase.Folder);
            testCase.verifyFalse(result.needsCommit);
        end

        function testAcknowledgeAcceptsBindingAndKeepsHadConflict(testCase)
            testCase.bindAndCommitOptics('A');
            testCase.verifyError(@() UMITRigLedger.acknowledgeConflict( ...
                testCase.Folder, 'Origin', 'test'), 'Umitoolbox:RigLedger:noConflict');
            requested = testCase.binding('B', 'LED_593nm');
            UMITRigLedger.rebind(testCase.Folder, requested, 'Origin', 'test');

            UMITRigLedger.acknowledgeConflict(testCase.Folder, 'Origin', 'test', ...
                'Comment', 'LED replaced before acquisition');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'BoundCommitted');
            testCase.verifyEqual(info.commitments.optics, requested.optics.hash);
            testCase.verifyTrue(info.hadConflict);
            testCase.verifyEqual(info.events{end}.comment, 'LED replaced before acquisition');
            result = UMITRigLedger.assertCanUseOptics(testCase.Folder);
            testCase.verifyEqual(result.hash, requested.optics.hash);

            % A later rebind back to the old optics is a new conflict.
            UMITRigLedger.rebind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');
            testCase.verifyEqual(UMITRigLedger.read(testCase.Folder).state, 'Conflicted');
        end

        function testCoregistrationConflictAndRequestedTransform(testCase)
            binding = testCase.binding('A', 'LED_632nm', 1);
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');
            UMITRigLedger.commit(testCase.Folder, 'coregistration', ...
                binding.coregistration.hash, 'Origin', 'test');
            other = testCase.binding('A', 'LED_632nm', 2);

            result = UMITRigLedger.assertCanUseCoregistration(testCase.Folder, ...
                'RequestedHash', binding.coregistration.hash);
            testCase.verifyFalse(result.needsCommit);
            testCase.verifyEqual(result.snapshot.transformMatrix, binding.coregistration.transformMatrix);
            testCase.verifyError(@() UMITRigLedger.assertCanUseCoregistration(testCase.Folder, ...
                'RequestedHash', other.coregistration.hash), ...
                'Umitoolbox:RigLedger:coregistrationConflict');

            conflict = UMITRigLedger.rebind(testCase.Folder, other, 'Origin', 'test');
            testCase.verifyEqual(conflict.parts, {'coregistration'});
            testCase.verifyEqual(UMITRigLedger.read(testCase.Folder).state, 'Conflicted');
        end

        function testCoregistrationGuardWithoutBoundTransform(testCase)
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');
            testCase.verifyError(@() UMITRigLedger.assertCanUseCoregistration(testCase.Folder), ...
                'Umitoolbox:RigLedger:noCoregistrationBound');

            UMITRigLedger.rebind(testCase.Folder, testCase.binding('A', 'LED_632nm', 1), 'Origin', 'test');
            other = testCase.binding('A', 'LED_632nm', 2);
            testCase.verifyError(@() UMITRigLedger.assertCanUseCoregistration(testCase.Folder, ...
                'RequestedHash', other.coregistration.hash), ...
                'Umitoolbox:RigLedger:coregistrationNotBound');
        end

        function testGuardOnFolderWithoutLedger(testCase)
            testCase.verifyError(@() UMITRigLedger.assertCanUseOptics(testCase.Folder), ...
                'Umitoolbox:RigLedger:noLedger');
            testCase.verifyEmpty(testCase.folderEntries());
        end
    end

    %% Fresh import
    methods (Test)
        function testFreshImportResetsAllCommitmentsAndKeepsBinding(testCase)
            binding = testCase.binding('A', 'LED_632nm', 1);
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');
            UMITRigLedger.commit(testCase.Folder, 'optics', binding.optics.hash, 'Origin', 'test');
            UMITRigLedger.commit(testCase.Folder, 'coregistration', ...
                binding.coregistration.hash, 'Origin', 'test');
            filterSet = UMITRigLedger.buildFilterSetSnapshot('GCaMP');
            UMITRigLedger.commit(testCase.Folder, 'filterSet', filterSet, 'Origin', 'test');
            eventsBefore = numel(UMITRigLedger.read(testCase.Folder).events);

            UMITRigLedger.recordFreshImport(testCase.Folder, 'Origin', 'ImagesClassification');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'BoundUncommitted');
            testCase.verifyEmpty(info.commitments.optics);
            testCase.verifyEmpty(info.commitments.coregistration);
            testCase.verifyEmpty(info.commitments.filterSets);
            testCase.verifyEqual(info.binding.opticsHash, binding.optics.hash);
            testCase.verifyEqual(info.binding.coregistrationHash, binding.coregistration.hash);
            % History is retained; the event records what was reset.
            testCase.verifyNumElements(info.events, eventsBefore + 1);
            event = info.events{end};
            testCase.verifyEqual(event.type, 'freshImport');
            testCase.verifyEqual(event.previousOpticsHash, binding.optics.hash);
            testCase.verifyEqual(event.previousCoregistrationHash, binding.coregistration.hash);
            testCase.verifyEqual(event.previousFilterSets.id, 'GCaMP');
            testCase.verifySubstring(regexprep(fileread(UMITRigLedger.filePath(testCase.Folder)), ...
                '\s', ''), '"previousFilterSets":[{');

            % Different optics, a new transform and a modified filter set are
            % all allowed afterwards, without conflict.
            other = testCase.binding('B', 'LED_593nm', 3);
            conflict = UMITRigLedger.rebind(testCase.Folder, other, 'Origin', 'test');
            testCase.verifyFalse(conflict.hasConflict);
            testCase.verifyTrue(UMITRigLedger.commit(testCase.Folder, 'optics', ...
                other.optics.hash, 'Origin', 'test'));
            testCase.verifyTrue(UMITRigLedger.commit(testCase.Folder, 'coregistration', ...
                other.coregistration.hash, 'Origin', 'test'));
            modified = filterSet;
            modified.emissionResponse(100) = 0.5;
            modified.hash = UMITRigLedger.computeHash(modified);
            testCase.verifyTrue(UMITRigLedger.assertCanUseFilterSet(testCase.Folder, modified));
            testCase.verifyEqual(UMITRigLedger.read(testCase.Folder).state, 'BoundCommitted');
            testCase.verifyFalse(UMITRigLedger.read(testCase.Folder).hadConflict);
        end

        function testFreshImportEndsConflictButHadConflictStays(testCase)
            % Decision (Bruno, 2026-10-15): hadConflict stays true by history.
            testCase.bindAndCommitOptics('A');
            UMITRigLedger.rebind(testCase.Folder, testCase.binding('B', 'LED_593nm'), 'Origin', 'test');
            testCase.verifyEqual(UMITRigLedger.read(testCase.Folder).state, 'Conflicted');

            UMITRigLedger.recordFreshImport(testCase.Folder, 'Origin', 'test');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'BoundUncommitted');
            testCase.verifyEqual(info.binding.sourceRig.rigID, 'B');
            testCase.verifyEmpty(info.conflictParts);
            testCase.verifyTrue(info.hadConflict);
            result = UMITRigLedger.assertCanUseOptics(testCase.Folder);
            testCase.verifyTrue(result.needsCommit);
        end

        function testFreshImportOfNoRigFolderStaysNoRig(testCase)
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');

            UMITRigLedger.recordFreshImport(testCase.Folder, 'Origin', 'test');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'NoRig');
            testCase.verifyEqual(cellfun(@(e) e.type, info.events, 'UniformOutput', false), ...
                {'noRig', 'freshImport'});
        end

        function testLegacyReimportEventIsInvalid(testCase)
            % Decision (Bruno, 2026-10-15): pre-release, no migration. A
            % ledger holding the retired 'reimport' event is invalid.
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            UMITRigLedger.recordFreshImport(testCase.Folder, 'Origin', 'test');
            path = UMITRigLedger.filePath(testCase.Folder);
            testCase.writeText(path, strrep(fileread(path), '"freshImport"', '"reimport"'));

            report = UMITRigLedger.validate(testCase.Folder);
            testCase.verifyFalse(report.isValid);
            testCase.verifySubstring(strjoin(report.issues, ' '), 'unknown event type');
            testCase.verifyError(@() UMITRigLedger.read(testCase.Folder), ...
                'Umitoolbox:RigLedger:invalidLedger');
        end

        function testDropCoregistrationBindingKeepsOptics(testCase)
            binding = testCase.binding('A', 'LED_632nm', 1);
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');
            UMITRigLedger.commit(testCase.Folder, 'coregistration', ...
                binding.coregistration.hash, 'Origin', 'test');
            UMITRigLedger.recordFreshImport(testCase.Folder, 'Origin', 'test');

            conflict = UMITRigLedger.dropCoregistrationBinding(testCase.Folder, 'Origin', 'test');
            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyFalse(conflict.hasConflict);
            testCase.verifyEqual(info.binding.opticsHash, binding.optics.hash);
            testCase.verifyEqual(info.binding.sourceRig.rigID, 'A');
            testCase.verifyEmpty(info.binding.coregistrationHash);
            testCase.verifyEqual(info.state, 'BoundUncommitted');

            % Nothing left to drop: no write.
            before = fileread(UMITRigLedger.filePath(testCase.Folder));
            UMITRigLedger.dropCoregistrationBinding(testCase.Folder, 'Origin', 'test');
            testCase.verifyEqual(fileread(UMITRigLedger.filePath(testCase.Folder)), before);
        end

        function testFreshImportRequiresLedger(testCase)
            testCase.verifyError(@() UMITRigLedger.recordFreshImport(testCase.Folder, ...
                'Origin', 'test'), 'Umitoolbox:RigLedger:noLedger');
            testCase.verifyEmpty(testCase.folderEntries());
        end
    end

    %% Filter-set rules
    methods (Test)
        function testFilterSetRules(testCase)
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');
            gcamp = UMITRigLedger.buildFilterSetSnapshot('GCaMP');
            jrgeco = UMITRigLedger.buildFilterSetSnapshot('jrGECO');

            % New ID: allowed, commit written.
            testCase.verifyTrue(UMITRigLedger.assertCanUseFilterSet(testCase.Folder, gcamp));
            testCase.verifyTrue(UMITRigLedger.commit(testCase.Folder, 'filterSet', gcamp, 'Origin', 'test'));
            % Same ID, same content: allowed, nothing written.
            before = fileread(UMITRigLedger.filePath(testCase.Folder));
            testCase.verifyFalse(UMITRigLedger.assertCanUseFilterSet(testCase.Folder, gcamp));
            testCase.verifyFalse(UMITRigLedger.commit(testCase.Folder, 'filterSet', gcamp, 'Origin', 'test'));
            testCase.verifyEqual(fileread(UMITRigLedger.filePath(testCase.Folder)), before);
            % Another new ID: allowed.
            testCase.verifyTrue(UMITRigLedger.commit(testCase.Folder, 'filterSet', jrgeco, 'Origin', 'test'));

            % Same ID, modified content: error.
            modified = gcamp;
            modified.emissionResponse(100) = 0.5;
            modified.hash = UMITRigLedger.computeHash(modified);
            testCase.verifyError(@() UMITRigLedger.assertCanUseFilterSet(testCase.Folder, modified), ...
                'Umitoolbox:RigLedger:filterSetModified');
            testCase.verifyError(@() UMITRigLedger.commit(testCase.Folder, 'filterSet', ...
                modified, 'Origin', 'test'), 'Umitoolbox:RigLedger:filterSetModified');

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(sort({info.commitments.filterSets.id}), {'GCaMP', 'jrGECO'});
        end

        function testFilterSetSnapshotWithWrongHashIsRejected(testCase)
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');
            tampered = UMITRigLedger.buildFilterSetSnapshot('GCaMP');
            tampered.emissionResponse(1) = 0.25;

            testCase.verifyError(@() UMITRigLedger.assertCanUseFilterSet(testCase.Folder, tampered), ...
                'Umitoolbox:RigLedger:hashMismatch');
        end
    end

    %% Canonical hashing
    methods (Test)
        function testCanonicalFormIsExact(testCase)
            value = struct('b', [1 0.5], 'a', 'x"y', 'hash', 'ignored', ...
                'c', struct('z', true, 'y', -1234567), 'd', 0);

            text = UMITRigLedger.canonicalJSON(value);

            testCase.verifyEqual(text, ...
                '{"a":"x\"y","b":[1e0,5e-1],"c":{"y":-1.23457e6,"z":true},"d":0}');
        end

        function testHashIsStableAndKeyOrderIndependent(testCase)
            snapshotA = UMITRigLedger.buildOpticsSnapshot(testCase.rigInfo('A'));
            snapshotB = UMITRigLedger.buildOpticsSnapshot(testCase.rigInfo('A'));
            reordered = orderfields(snapshotA, fliplr(fieldnames(snapshotA)'));
            reordered.cameras = orderfields(reordered.cameras, ...
                fliplr(fieldnames(reordered.cameras)'));

            testCase.verifyEqual(snapshotA.hash, snapshotB.hash);
            testCase.verifyEqual(UMITRigLedger.computeHash(reordered), snapshotA.hash);
            testCase.verifyMatches(snapshotA.hash, '^sha256:[0-9a-f]{64}$');
        end

        function testHashIgnoresRigIdentityButTracksContent(testCase)
            sameOptics = UMITRigLedger.buildOpticsSnapshot(testCase.rigInfo('B'));
            original = UMITRigLedger.buildOpticsSnapshot(testCase.rigInfo('A'));
            changed = UMITRigLedger.buildOpticsSnapshot(testCase.rigInfo('A', 'LED_593nm'));
            libraryOrder = testCase.rigInfo('A');
            libraryOrder.cameras = libraryOrder.cameras([2 1]);
            libraryOrder.illuminations = libraryOrder.illuminations([2 1]);

            testCase.verifyEqual(sameOptics.hash, original.hash);
            testCase.verifyNotEqual(changed.hash, original.hash);
            testCase.verifyEqual(UMITRigLedger.buildOpticsSnapshot(libraryOrder).hash, original.hash);
        end

        function testSnapshotsEmbedSpectraAndSurviveRoundTrip(testCase)
            binding = testCase.binding('A', 'LED_632nm', 1);
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');

            result = UMITRigLedger.assertCanUseOptics(testCase.Folder);
            testCase.verifyEqual(result.hash, binding.optics.hash);
            testCase.verifyEqual(UMITRigLedger.computeHash(result.snapshot), binding.optics.hash);
            testCase.verifySize(result.snapshot.cameras(1).response, [1 301]);
            testCase.verifyEqual(result.snapshot.cameras(1).spectrumSource, 'builtIn');
            red = UMITRigStore.getSpectrum('illumination', 'LED_632nm');
            redRow = strcmp({result.snapshot.illuminations.name}, 'red');
            testCase.verifyEqual(result.snapshot.illuminations(redRow).response(:), ...
                red.response(:), 'RelTol', 1e-5);
        end

        function testSingleCameraSnapshotStaysAList(testCase)
            info = testCase.rigInfo('A');
            info.cameras = info.cameras(1);
            snapshot = UMITRigLedger.buildOpticsSnapshot(info);
            UMITRigLedger.bind(testCase.Folder, UMITRigLedger.makeBinding( ...
                UMITRigLedger.sourceRigFromRigInfo(info), snapshot, []), 'Origin', 'test');

            text = fileread(UMITRigLedger.filePath(testCase.Folder));
            decoded = jsondecode(text);
            testCase.verifyTrue(UMITRigLedger.validate(testCase.Folder).isValid);
            testCase.verifyTrue(isstruct(decoded.snapshots.optics.cameras) || ...
                iscell(decoded.snapshots.optics.cameras));
            testCase.verifySubstring(regexprep(text, '\s', ''), '"cameras":[{');
        end
    end

    %% Integrity
    methods (Test)
        function testTamperedSnapshotInvalidatesLedger(testCase)
            testCase.bindAndCommitOptics('A');
            path = UMITRigLedger.filePath(testCase.Folder);
            text = fileread(path);
            testCase.writeText(path, regexprep(text, '"model": "D1024"', '"model": "D9999"', 'once'));

            report = UMITRigLedger.validate(testCase.Folder);
            testCase.verifyFalse(report.isValid);
            testCase.verifyNotEmpty(report.issues);
            testCase.verifyError(@() UMITRigLedger.read(testCase.Folder), ...
                'Umitoolbox:RigLedger:invalidLedger');
            testCase.verifyError(@() UMITRigLedger.assertCanUseOptics(testCase.Folder), ...
                'Umitoolbox:RigLedger:invalidLedger');
            testCase.verifyError(@() UMITRigLedger.recordFreshImport(testCase.Folder, 'Origin', 'test'), ...
                'Umitoolbox:RigLedger:invalidLedger');
        end

        function testInconsistentEventsInvalidateLedger(testCase)
            testCase.bindAndCommitOptics('A');
            UMITRigLedger.rebind(testCase.Folder, testCase.binding('B', 'LED_593nm'), 'Origin', 'test');
            path = UMITRigLedger.filePath(testCase.Folder);
            raw = jsondecode(fileread(path));
            types = cellfun(@(e) e.type, raw.events, 'UniformOutput', false);
            raw.events(strcmp(types, 'conflict')) = [];
            testCase.writeText(path, jsonencode(raw));

            report = UMITRigLedger.validate(testCase.Folder);
            testCase.verifyFalse(report.isValid);
            testCase.verifySubstring(strjoin(report.issues, ' '), 'conflict');
        end

        function testMalformedAndForeignFilesAreInvalid(testCase)
            path = UMITRigLedger.filePath(testCase.Folder);
            testCase.writeText(path, '{not json');
            testCase.verifyFalse(UMITRigLedger.validate(testCase.Folder).isValid);
            testCase.verifyError(@() UMITRigLedger.read(testCase.Folder), ...
                'Umitoolbox:RigLedger:invalidLedger');
            testCase.verifyError(@() UMITRigLedger.bind(testCase.Folder, ...
                testCase.binding('A'), 'Origin', 'test'), 'Umitoolbox:RigLedger:invalidLedger');
            testCase.verifyEqual(fileread(path), '{not json');

            testCase.writeText(path, '{"format":"other","schemaVersion":1}');
            testCase.verifyFalse(UMITRigLedger.validate(testCase.Folder).isValid);
        end

        function testUnsupportedSchemaVersionIsInvalid(testCase)
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            path = UMITRigLedger.filePath(testCase.Folder);
            raw = jsondecode(fileread(path));
            raw.schemaVersion = 2;
            testCase.writeText(path, jsonencode(raw));

            report = UMITRigLedger.validate(testCase.Folder);
            testCase.verifyFalse(report.isValid);
            testCase.verifySubstring(strjoin(report.issues, ' '), 'schemaVersion');
        end
    end

    %% Atomic writes and locking
    methods (Test)
        function testWritesLeaveNoTemporaryFiles(testCase)
            testCase.bindAndCommitOptics('A');
            UMITRigLedger.rebind(testCase.Folder, testCase.binding('B', 'LED_593nm'), 'Origin', 'test');
            UMITRigLedger.acknowledgeConflict(testCase.Folder, 'Origin', 'test');

            testCase.verifyEqual(testCase.folderEntries(), {UMITRigLedger.FileName});
        end

        function testHeldLockBlocksWritesAndKeepsFile(testCase)
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            path = UMITRigLedger.filePath(testCase.Folder);
            before = fileread(path);
            lockFolder = fullfile(testCase.Folder, UMITRigLedger.LockFolderName);
            mkdir(lockFolder);

            testCase.verifyError(@() UMITRigLedger.bind(testCase.Folder, ...
                testCase.binding('A'), 'Origin', 'test'), 'Umitoolbox:RigLedger:ledgerLocked');
            testCase.verifyEqual(fileread(path), before);
            testCase.verifyEqual(UMITRigLedger.read(testCase.Folder).state, 'NoRig');

            testCase.verifyError(@() UMITRigLedger.clearStaleLock(testCase.Folder), ...
                'Umitoolbox:RigLedger:lockNotStale');
            testCase.verifyTrue(UMITRigLedger.clearStaleLock(testCase.Folder, ...
                'MinimumAge', seconds(0)));
            testCase.verifyFalse(isfolder(lockFolder));
            UMITRigLedger.bind(testCase.Folder, testCase.binding('A'), 'Origin', 'test');
            testCase.verifyEqual(UMITRigLedger.read(testCase.Folder).state, 'BoundUncommitted');
        end

        function testRejectedWriteKeepsPreviousLedger(testCase)
            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            path = UMITRigLedger.filePath(testCase.Folder);
            before = fileread(path);
            binding = testCase.binding('A');
            binding.optics.cameras(1).response(1) = 0.123; % stale hash

            testCase.verifyError(@() UMITRigLedger.bind(testCase.Folder, binding, ...
                'Origin', 'test'), 'Umitoolbox:RigLedger:hashMismatch');
            testCase.verifyEqual(fileread(path), before);
            testCase.verifyEqual(testCase.folderEntries(), {UMITRigLedger.FileName});
        end
    end

    %% Helpers
    methods
        function info = rigInfo(~, rigID, redSpectrum)
            if nargin < 3
                redSpectrum = 'LED_632nm';
            end
            info = struct( ...
                'uuid', ['uuid-' rigID], ...
                'rigID', rigID, ...
                'displayName', ['Rig ' rigID], ...
                'modifiedOn', datetime(2026, 10, 1), ...
                'cameras', struct('index', {1, 2}, 'spectrumID', {'PF1024', 'PF1024'}, ...
                    'model', {'D1024', 'D1024'}), ...
                'illuminations', struct('name', {'red', 'green'}, ...
                    'spectrumID', {redSpectrum, 'LED_521nm'}));
        end

        function binding = binding(testCase, rigID, redSpectrum, coregShift)
            if nargin < 3
                redSpectrum = 'LED_632nm';
            end
            info = testCase.rigInfo(rigID, redSpectrum);
            coregistration = [];
            if nargin >= 4
                matrix = eye(3);
                matrix(3, 1) = coregShift;
                coregistration = UMITRigLedger.coregistrationSnapshotFromTransform( ...
                    sprintf('resource-%d', coregShift), 'abc123', affine2d(matrix));
            end
            binding = UMITRigLedger.makeBinding(UMITRigLedger.sourceRigFromRigInfo(info), ...
                UMITRigLedger.buildOpticsSnapshot(info), coregistration);
        end

        function binding = bindAndCommitOptics(testCase, rigID)
            binding = testCase.binding(rigID);
            UMITRigLedger.bind(testCase.Folder, binding, 'Origin', 'test');
            UMITRigLedger.commit(testCase.Folder, 'optics', binding.optics.hash, 'Origin', 'test');
        end

        function names = folderEntries(testCase)
            listing = dir(testCase.Folder);
            names = setdiff({listing.name}, {'.', '..'});
        end

        function writeText(~, path, text)
            fid = fopen(path, 'w', 'n', 'UTF-8');
            fprintf(fid, '%s', text);
            fclose(fid);
        end
    end
end
