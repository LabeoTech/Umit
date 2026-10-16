classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestUMITProjectLedger < matlab.unittest.TestCase
    %TESTUMITPROJECTLEDGER Folder project ledger (project redesign Phase 3).
    %
    %   Histories of bind, rebind, unbind and orphan removal; frozen,
    %   unvalidated snapshots; strict identity; event-order integrity;
    %   fresh-import events; read-only reads. No project store is used.

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

    methods (Test)
        function testBindRebindUnbindHistory(testCase)
            first = testCase.binding();
            UMITProjectLedger.bind(testCase.Folder, first, testCase.snapshot('S1'), 'Origin', 'test');
            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'Bound');
            testCase.verifyEqual(info.binding.sessionUUID, first.sessionUUID);
            testCase.verifyError(@() UMITProjectLedger.bind(testCase.Folder, testCase.binding(), ...
                testCase.snapshot('S1'), 'Origin', 'test'), 'Umitoolbox:ProjectLedger:alreadyBound');

            second = first;
            second.subjectUUID = testCase.uuid();
            UMITProjectLedger.rebind(testCase.Folder, second, testCase.snapshot('S2'), 'Origin', 'test');
            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.binding.subjectUUID, second.subjectUUID);
            testCase.verifyEqual(info.snapshot.subjectID, 'S2');
            rebind = info.events{end};
            testCase.verifyEqual(rebind.type, 'rebind');
            testCase.verifyEqual(rebind.previousBinding.subjectUUID, first.subjectUUID);
            testCase.verifyNotEmpty(rebind.machine);

            UMITProjectLedger.unbind(testCase.Folder, 'Reason', 'sessionRemoved', 'Origin', 'test');
            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'Unbound');
            testCase.verifyEmpty(info.binding);
            testCase.verifyEqual(cellfun(@(e) e.type, info.events, 'UniformOutput', false), ...
                {'bind', 'rebind', 'unbind'});
            testCase.verifyError(@() UMITProjectLedger.rebind(testCase.Folder, second, ...
                testCase.snapshot('S2'), 'Origin', 'test'), 'Umitoolbox:ProjectLedger:notBound');

            % Binding again after an unbind keeps the whole history.
            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), testCase.snapshot('S3'), ...
                'Origin', 'test');
            testCase.verifyNumElements(UMITProjectLedger.read(testCase.Folder).events, 4);
        end

        function testOrphanRemovalIsAnEvent(testCase)
            original = testCase.binding();
            UMITProjectLedger.bind(testCase.Folder, original, testCase.snapshot('S1'), 'Origin', 'test');

            UMITProjectLedger.recordOrphanRemoval(testCase.Folder, 'Origin', 'test');

            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'Unbound');
            event = info.events{end};
            testCase.verifyEqual(event.type, 'orphanRemoval');
            testCase.verifyEqual(event.previousBinding.bindingUUID, original.bindingUUID);
            testCase.verifyNotEmpty(event.machine);
            testCase.verifyTrue(isfile(UMITProjectLedger.filePath(testCase.Folder)));
        end

        function testSnapshotContentIsNeverValidated(testCase)
            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), testCase.snapshot('S1'), ...
                'Origin', 'test');
            path = UMITProjectLedger.filePath(testCase.Folder);
            text = fileread(path);
            % Edit, retype and add snapshot content: still valid.
            text = strrep(text, '"subjectID": "S1"', '"subjectID": 42, "anything": {"x": [1, null]}');
            text = strrep(text, '"projectName": "Project A"', '"projectName": "Renamed by hand"');
            testCase.writeText(path, text);

            report = UMITProjectLedger.validate(testCase.Folder);
            testCase.verifyTrue(report.isValid, strjoin(report.issues, ' | '));
            testCase.verifyEqual(report.info.snapshot.subjectID, 42);
        end

        function testIdentityIsStrictlyValidated(testCase)
            binding = testCase.binding();
            UMITProjectLedger.bind(testCase.Folder, binding, testCase.snapshot('S1'), 'Origin', 'test');
            path = UMITProjectLedger.filePath(testCase.Folder);
            original = fileread(path);

            testCase.writeText(path, strrep(original, binding.sessionUUID, 'not-a-uuid'));
            testCase.verifyFalse(UMITProjectLedger.validate(testCase.Folder).isValid);
            testCase.verifyError(@() UMITProjectLedger.read(testCase.Folder), ...
                'Umitoolbox:ProjectLedger:invalidLedger');

            testCase.writeText(path, strrep(original, '"folderRole": "processedDataFolder"', ...
                '"folderRole": "processedDataFolder", "extra": 1'));
            testCase.verifyFalse(UMITProjectLedger.validate(testCase.Folder).isValid);

            testCase.writeText(path, strrep(original, '"processedDataFolder"', '"elsewhere"'));
            testCase.verifyFalse(UMITProjectLedger.validate(testCase.Folder).isValid);

            testCase.verifyError(@() UMITProjectLedger.bind(testCase.Folder, ...
                rmfield(testCase.binding(), 'folderRole'), struct(), 'Origin', 'test'), ...
                'Umitoolbox:ProjectLedger:invalidBinding');
        end

        function testEventOrderIsChecked(testCase)
            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), testCase.snapshot('S1'), ...
                'Origin', 'test');
            UMITProjectLedger.unbind(testCase.Folder, 'Reason', 'x', 'Origin', 'test');
            path = UMITProjectLedger.filePath(testCase.Folder);
            raw = jsondecode(fileread(path));
            raw.events = raw.events(2); % an unbind without a binding
            testCase.writeText(path, jsonencode(raw));

            report = UMITProjectLedger.validate(testCase.Folder);
            testCase.verifyFalse(report.isValid);
            testCase.verifySubstring(strjoin(report.issues, ' '), 'not bound');
        end

        function testFreshImportOnlyForBoundFolders(testCase)
            testCase.verifyFalse(UMITProjectLedger.recordFreshImport(testCase.Folder, ...
                'Origin', 'ImagesClassification'));
            testCase.verifyFalse(UMITProjectLedger.exists(testCase.Folder), ...
                'an unbound folder gets no project ledger');

            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), testCase.snapshot('S1'), ...
                'Origin', 'test');
            testCase.verifyTrue(UMITProjectLedger.recordFreshImport(testCase.Folder, ...
                'Origin', 'ImagesClassification'));
            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.events{end}.type, 'freshImport');
            testCase.verifyEqual(info.state, 'Bound');

            UMITProjectLedger.unbind(testCase.Folder, 'Reason', 'x', 'Origin', 'test');
            before = fileread(UMITProjectLedger.filePath(testCase.Folder));
            testCase.verifyFalse(UMITProjectLedger.recordFreshImport(testCase.Folder, ...
                'Origin', 'ImagesClassification'));
            testCase.verifyEqual(fileread(UMITProjectLedger.filePath(testCase.Folder)), before);
        end

        function testReadingNeverWrites(testCase)
            info = UMITProjectLedger.read(testCase.Folder);
            testCase.verifyEqual(info.state, 'Unknown');
            testCase.verifyEqual(UMITProjectStore.getProjectBindingState(testCase.Folder), 'unbound');
            testCase.verifyEmpty(testCase.entries());

            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), testCase.snapshot('S1'), ...
                'Origin', 'test');
            path = UMITProjectLedger.filePath(testCase.Folder);
            before = fileread(path);
            listing = dir(path);
            UMITProjectLedger.read(testCase.Folder);
            UMITProjectLedger.validate(testCase.Folder);
            UMITProjectLedger.isBound(testCase.Folder);
            UMITProjectStore.getProjectBindingState(testCase.Folder);
            UMITProjectStore.readProjectBinding(testCase.Folder);
            testCase.verifyEqual(fileread(path), before);
            testCase.verifyEqual(dir(path).datenum, listing.datenum);
            testCase.verifyEqual(testCase.entries(), {UMITProjectLedger.FileName});
        end

        function testLegacyUmitlinkIsIgnored(testCase)
            % Brief section 4 Q8: a folder with only a .umitlink reads as unbound.
            ProjectBinding = struct('version', 1);
            save(fullfile(testCase.Folder, 'UMITProjectBinding.umitlink'), 'ProjectBinding', '-mat');
            testCase.verifyEqual(UMITProjectStore.getProjectBindingState(testCase.Folder), 'unbound');
            testCase.verifyError(@() UMITProjectStore.readProjectBinding(testCase.Folder), ...
                'Umitoolbox:UMITProjectStore:bindingReadFailed');
            testCase.verifyTrue(isfile(fullfile(testCase.Folder, 'UMITProjectBinding.umitlink')), ...
                'never deleted automatically');
        end

        function testHeldLockBlocksWrites(testCase)
            mkdir(fullfile(testCase.Folder, UMITProjectLedger.LockFolderName));
            testCase.verifyError(@() UMITProjectLedger.bind(testCase.Folder, testCase.binding(), ...
                testCase.snapshot('S1'), 'Origin', 'test'), 'Umitoolbox:ProjectLedger:ledgerLocked');
            testCase.verifyTrue(UMITProjectLedger.clearStaleLock(testCase.Folder, ...
                'MinimumAge', seconds(0)));
            UMITProjectLedger.bind(testCase.Folder, testCase.binding(), testCase.snapshot('S1'), ...
                'Origin', 'test');
            testCase.verifyEqual(testCase.entries(), {UMITProjectLedger.FileName});
        end
    end

    methods
        function value = uuid(~)
            value = lower(char(java.util.UUID.randomUUID()));
        end

        function binding = binding(testCase)
            binding = struct('bindingUUID', testCase.uuid(), 'projectUUID', testCase.uuid(), ...
                'subjectUUID', testCase.uuid(), 'sessionUUID', testCase.uuid(), ...
                'folderRole', 'processedDataFolder');
        end

        function snapshot = snapshot(~, subjectID)
            ProjectInfo = struct('projectName', 'Project A', 'principalInvestigator', 'PI', ...
                'experimenters', {{'E1', 'E2'}}, 'projectStartDate', datetime(2026, 1, 1), ...
                'subjectRegistry', struct('uuid', {}));
            snapshot = UMITProjectLedger.makeSnapshot(ProjectInfo, subjectID, 'Session_01');
        end

        function names = entries(testCase)
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
