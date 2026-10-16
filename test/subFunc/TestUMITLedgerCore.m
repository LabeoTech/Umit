classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestUMITLedgerCore < matlab.unittest.TestCase
    %TESTUMITLEDGERCORE Tests of the shared ledger engine (UMITLedgerCore).
    %
    %   Covers canonical serialization, hash stability, byte hashing, the
    %   ledger header and event integrity checks, atomic writes (success
    %   and rejected candidates) and locking. A private test spec is used,
    %   so these tests are independent of the rig and project ledgers.

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

    %% Canonical serialization and hashing
    methods (Test)
        function testCanonicalFormIsExact(testCase)
            value = struct('b', [1 0.5], 'a', 'x"y', 'hash', 'ignored', ...
                'c', struct('z', true, 'y', -1234567), 'd', 0, ...
                'e', string(missing), 'f', [1 2; 3 4], 'g', sprintf('a\tb'));

            text = UMITLedgerCore.canonicalJSON(value);

            testCase.verifyEqual(text, ['{"a":"x\"y","b":[1e0,5e-1],' ...
                '"c":{"y":-1.23457e6,"z":true},"d":0,"e":null,' ...
                '"f":[[1e0,2e0],[3e0,4e0]],"g":"a\u0009b"}']);
        end

        function testListFieldsStayListsWithOneElement(testCase)
            value = struct('items', struct('k', 1), 'other', struct('k', 1));
            options = struct('listFields', {{'items'}});

            text = UMITLedgerCore.canonicalJSON(value, options);

            testCase.verifyEqual(text, '{"items":[{"k":1e0}],"other":{"k":1e0}}');
        end

        function testHashIsStableAndKeyOrderIndependent(testCase)
            a = struct('x', 1.5, 'y', {{'p', 'q'}}, 'z', struct('m', 'n', 'o', 2));
            b = struct('z', struct('o', 2, 'm', 'n'), 'y', {{'p', 'q'}}, 'x', 1.5);
            a.hash = 'anything';

            hashA = UMITLedgerCore.computeHash(a);

            testCase.verifyEqual(hashA, UMITLedgerCore.computeHash(b));
            testCase.verifyEqual(hashA, UMITLedgerCore.computeHash(a));
            testCase.verifyMatches(hashA, '^sha256:[0-9a-f]{64}$');
            testCase.verifyNotEqual(UMITLedgerCore.computeHash(setfield(b, 'x', 1.6)), hashA); %#ok<SFLD>
            % Values equal at 6 significant digits hash identically.
            testCase.verifyEqual(UMITLedgerCore.computeHash(setfield(b, 'x', 1.5000001)), hashA); %#ok<SFLD>
        end

        function testHashSurvivesJsonRoundTrip(testCase)
            value = struct('w', UMITLedgerCore.quantize(linspace(0.1, 7.3, 11)), ...
                'name', 'abc', 'flag', false);

            decoded = jsondecode(jsonencode(value));

            testCase.verifyEqual(UMITLedgerCore.computeHash(decoded), ...
                UMITLedgerCore.computeHash(value));
        end

        function testByteAndTextHashesMatchKnownVectors(testCase)
            abc = 'sha256:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad';
            empty = 'sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855';

            testCase.verifyEqual(UMITLedgerCore.hashBytes(uint8('abc')), abc);
            testCase.verifyEqual(UMITLedgerCore.hashText("abc"), abc);
            testCase.verifyEqual(UMITLedgerCore.hashBytes(uint8([])), empty);
            testCase.verifyEqual(UMITLedgerCore.shortHash(abc), 'sha256:ba7816bf8f01');
            testCase.verifyEqual(UMITLedgerCore.shortHash(''), '(none)');
            testCase.verifyError(@() UMITLedgerCore.hashBytes(single(1)), ...
                'Umitoolbox:LedgerCore:invalidInput');
        end

        function testUnhashableValuesUseTheCallerPrefix(testCase)
            options = struct('errorPrefix', 'Umitoolbox:TestLedger');

            testCase.verifyError(@() UMITLedgerCore.computeHash(struct('v', NaN), options), ...
                'Umitoolbox:TestLedger:hashFailed');
            testCase.verifyError(@() UMITLedgerCore.computeHash(struct('v', ones(2, 2, 2))), ...
                'Umitoolbox:LedgerCore:hashFailed');
            testCase.verifyError(@() UMITLedgerCore.computeHash(struct('v', @sin), options), ...
                'Umitoolbox:TestLedger:hashFailed');
        end
    end

    %% Header and event integrity
    methods (Test)
        function testNewLedgerAndEventsAreValid(testCase)
            spec = testCase.spec();
            ledger = UMITLedgerCore.newLedger(spec, struct('payload', 'x'));
            ledger = UMITLedgerCore.appendEvent(ledger, UMITLedgerCore.newEvent('first', 'test'));

            testCase.verifyEqual(fieldnames(ledger)', ...
                {'format', 'schemaVersion', 'ledgerUUID', 'createdOn', 'payload', 'events'});
            testCase.verifyEmpty(UMITLedgerCore.headerIssues(ledger, spec));
            testCase.verifyEmpty(UMITLedgerCore.eventIssues(ledger));
            event = ledger.events{1};
            testCase.verifyEqual(fieldnames(event)', UMITLedgerCore.BaseEventFields);
            testCase.verifyEqual(event.toolboxVersion, UMITLedgerCore.ToolboxVersion);
            testCase.verifyFalse(any(isfield(event, {'user', 'host'})));
            testCase.verifyTrue(isdatetime(UMITLedgerCore.localTime(event.timestamp)));
        end

        function testHeaderIssuesDetectForeignFormatAndVersion(testCase)
            spec = testCase.spec();
            ledger = UMITLedgerCore.newLedger(spec, struct('payload', 'x'));
            ledger.format = 'other';
            ledger.schemaVersion = 9;
            ledger.ledgerUUID = '';

            issues = strjoin(UMITLedgerCore.headerIssues(ledger, spec), ' ');

            testCase.verifySubstring(issues, 'format');
            testCase.verifySubstring(issues, 'schemaVersion 9');
            testCase.verifySubstring(issues, 'ledgerUUID');
        end

        function testEventIssuesDetectMissingFieldsAndDuplicates(testCase)
            ledger = UMITLedgerCore.newLedger(testCase.spec(), struct('payload', 'x'));
            testCase.verifySubstring(strjoin(UMITLedgerCore.eventIssues(ledger)), 'no events');

            event = UMITLedgerCore.newEvent('a', 'test');
            ledger = UMITLedgerCore.appendEvent(ledger, event);
            ledger = UMITLedgerCore.appendEvent(ledger, event);
            ledger = UMITLedgerCore.appendEvent(ledger, rmfield(event, 'origin'));

            issues = strjoin(UMITLedgerCore.eventIssues(ledger), ' ');
            testCase.verifySubstring(issues, 'Event 2 (a) reuses eventUUID');
            testCase.verifySubstring(issues, 'Event 3 is missing field(s): origin');
        end

        function testDecodeAndNormalizeRejectMalformedFiles(testCase)
            spec = testCase.spec();
            path = fullfile(testCase.Folder, spec.fileName);
            testCase.writeText(path, '{not json');
            testCase.verifyError(@() UMITLedgerCore.decodeFile(path, spec), ...
                'Umitoolbox:TestLedger:readFailed');
            testCase.verifyError(@() UMITLedgerCore.decodeFile( ...
                fullfile(testCase.Folder, 'absent.json'), spec), 'Umitoolbox:TestLedger:readFailed');

            testCase.writeText(path, '{"format":"x","schemaVersion":1}');
            try
                UMITLedgerCore.normalizeHeader(UMITLedgerCore.decodeFile(path, spec), spec);
                testCase.verifyFail('Expected missing fields.');
            catch ME
                testCase.verifySubstring(ME.message, ...
                    'Test ledger is missing field(s): ledgerUUID, createdOn, payload, events.');
            end
            testCase.writeText(path, '[1,2]');
            testCase.verifyError(@() UMITLedgerCore.normalizeHeader( ...
                UMITLedgerCore.decodeFile(path, spec), spec), ?MException);
        end
    end

    %% Atomic writes and locking
    methods (Test)
        function testWriteInstallsVerifiedLedgerWithoutLeftovers(testCase)
            spec = testCase.spec();
            ledger = testCase.oneEventLedger(spec);

            UMITLedgerCore.writeLedger(testCase.Folder, spec, ledger, @(p) testCase.verifier(p, spec));
            ledger = UMITLedgerCore.appendEvent(ledger, UMITLedgerCore.newEvent('second', 'test'));
            UMITLedgerCore.writeLedger(testCase.Folder, spec, ledger, @(p) testCase.verifier(p, spec));

            testCase.verifyEqual(testCase.folderEntries(), {spec.fileName});
            written = UMITLedgerCore.normalizeHeader(UMITLedgerCore.decodeFile( ...
                UMITLedgerCore.filePath(testCase.Folder, spec), spec), spec);
            testCase.verifyNumElements(written.events, 2);
            testCase.verifyEqual(written.ledgerUUID, ledger.ledgerUUID);
        end

        function testRejectedCandidateKeepsPreviousFile(testCase)
            spec = testCase.spec();
            ledger = testCase.oneEventLedger(spec);
            UMITLedgerCore.writeLedger(testCase.Folder, spec, ledger, @(p) testCase.verifier(p, spec));
            path = UMITLedgerCore.filePath(testCase.Folder, spec);
            before = fileread(path);
            bad = ledger;
            bad.format = 'other';

            testCase.verifyError(@() UMITLedgerCore.writeLedger(testCase.Folder, spec, bad, ...
                @(p) testCase.verifier(p, spec)), 'Umitoolbox:TestLedger:writeFailed');

            testCase.verifyEqual(fileread(path), before);
            testCase.verifyEqual(testCase.folderEntries(), {spec.fileName});
        end

        function testLockBlocksSecondHolderAndIsReleased(testCase)
            spec = testCase.spec();
            lockFolder = fullfile(testCase.Folder, spec.lockFolderName);

            lock = UMITLedgerCore.acquireLock(testCase.Folder, spec, 'first'); %#ok<NASGU> released by clear below
            testCase.verifyTrue(isfile(fullfile(lockFolder, 'lockInfo.mat')));
            try
                UMITLedgerCore.acquireLock(testCase.Folder, spec, 'second');
                testCase.verifyFail('Expected a held lock.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:TestLedger:ledgerLocked');
                testCase.verifySubstring(ME.message, 'TestLedger.clearStaleLock(folder)');
            end
            clear lock
            testCase.verifyFalse(isfolder(lockFolder));
        end

        function testClearStaleLockHonoursMinimumAge(testCase)
            spec = testCase.spec();
            testCase.verifyFalse(UMITLedgerCore.clearStaleLock(testCase.Folder, spec, minutes(10)));
            mkdir(fullfile(testCase.Folder, spec.lockFolderName));

            testCase.verifyError(@() UMITLedgerCore.clearStaleLock(testCase.Folder, spec, ...
                minutes(10)), 'Umitoolbox:TestLedger:lockNotStale');
            testCase.verifyTrue(UMITLedgerCore.clearStaleLock(testCase.Folder, spec, seconds(0)));
            testCase.verifyEmpty(testCase.folderEntries());
        end

        function testInputChecksUseTheSpecPrefix(testCase)
            spec = testCase.spec();

            testCase.verifyError(@() UMITLedgerCore.existingFolder( ...
                fullfile(testCase.Folder, 'missing'), spec), 'Umitoolbox:TestLedger:invalidInput');
            testCase.verifyError(@() UMITLedgerCore.folderText('', spec), ...
                'Umitoolbox:TestLedger:invalidInput');
            testCase.verifyError(@() UMITLedgerCore.requireOrigin('', spec), ...
                'Umitoolbox:TestLedger:invalidInput');
            testCase.verifyEqual(UMITLedgerCore.folderText([testCase.Folder filesep], spec), ...
                testCase.Folder);
            testCase.verifyError(@() UMITLedgerCore.makeSpec(struct('fileName', 'x')), ...
                'Umitoolbox:LedgerCore:invalidSpec');
        end
    end

    %% Helpers
    methods
        function spec = spec(~)
            spec = UMITLedgerCore.makeSpec(struct( ...
                'fileName', 'TestLedger.json', ...
                'lockFolderName', 'TestLedger.lock', ...
                'format', 'umit-test-ledger', ...
                'schemaVersion', 1, ...
                'errorPrefix', 'Umitoolbox:TestLedger', ...
                'label', 'test ledger', ...
                'clearLockCall', 'TestLedger.clearStaleLock(folder)', ...
                'bodyFields', {{'payload'}}));
        end

        function ledger = oneEventLedger(~, spec)
            ledger = UMITLedgerCore.newLedger(spec, struct('payload', 'x'));
            ledger = UMITLedgerCore.appendEvent(ledger, UMITLedgerCore.newEvent('first', 'test'));
        end

        function verifier(~, path, spec)
            ledger = UMITLedgerCore.normalizeHeader(UMITLedgerCore.decodeFile(path, spec), spec);
            issues = [UMITLedgerCore.headerIssues(ledger, spec), UMITLedgerCore.eventIssues(ledger)];
            if ~isempty(issues)
                error('%s', strjoin(issues, ' | '));
            end
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
