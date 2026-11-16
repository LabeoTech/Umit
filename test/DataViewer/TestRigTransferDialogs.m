classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestRigTransferDialogs < matlab.unittest.TestCase
    %TESTRIGTRANSFERDIALOGS DataViewer Export Rig / Import Rig workflows.
    %
    %   promptExportRig and promptImportRig are exercised through their
    %   non-interactive mode (no dialog is opened); the dialogs themselves
    %   (uigetfile, uiputfile, uiconfirm, uialert) are checked manually. The
    %   store behavior behind them is covered by TestUMITRigStore. These tests
    %   guard the helper contract: a file written by Export Rig installs
    %   through Import Rig, a duplicate or ID clash never overwrites a rig,
    %   the "copy" path finds a free rig ID, and every failure comes back as
    %   an error result instead of an exception.

    properties
        TempRoot
        ExtraRigRoots = {}
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.addTeardown(@() iRemoveFolder(testCase.TempRoot));
            testCase.ExtraRigRoots = {};
            testCase.addTeardown(@() cellfun(@iRemoveFolder, testCase.ExtraRigRoots));
        end
    end

    methods (Test)
        function testExportThenImportRestoresTheRig(testCase)
            source = testCase.createLibraryRig('Source');
            sourceInfo = source.getRigInfo();
            exportFile = fullfile(testCase.TempRoot, 'rigExport');

            exported = promptExportRig([], 'Interactive', false, ...
                'RigUUID', sourceInfo.uuid, 'FilePath', exportFile);
            testCase.verifyEqual(exported.severity, 'info', exported.message);
            testCase.verifyFalse(exported.changed);
            testCase.verifyEqual(exported.filePath, [exportFile '.umitrig']);
            testCase.verifyTrue(isfile(exported.filePath), '.umitrig is added to the name.');
            testCase.verifySubstring(exported.message, sourceInfo.rigID);

            % The library no longer has the rig (as on the target machine).
            sourceHash = UMITRigLedger.buildOpticsSnapshot(source).hash;
            rmdir(source.RigRoot, 's');
            testCase.verifyFalse(UMITRigStore.rigExists(sourceInfo.uuid));

            imported = promptImportRig([], 'Interactive', false, ...
                'FilePath', exported.filePath);
            testCase.verifyEqual(imported.severity, 'info', imported.message);
            testCase.verifyTrue(imported.changed);
            testCase.verifyEqual(imported.rigUUID, sourceInfo.uuid);
            testCase.verifyEqual(imported.rigID, sourceInfo.rigID);
            installed = UMITRigStore.open(sourceInfo.uuid);
            testCase.ExtraRigRoots{end+1} = installed.RigRoot;
            testCase.verifyTrue(installed.validate('Mode', 'full').isValid);
            testCase.verifyEqual(UMITRigLedger.buildOpticsSnapshot(installed).hash, sourceHash);
        end

        function testImportOfAnExistingRigIsRefusedByDefault(testCase)
            source = testCase.createLibraryRig('Dup');
            exportFile = testCase.exportRigFile(source);
            before = UMITRigStore.listRigs();

            result = promptImportRig([], 'Interactive', false, 'FilePath', exportFile);
            testCase.verifyFalse(result.changed);
            testCase.verifyEqual(result.severity, 'info');
            testCase.verifySubstring(result.message, 'cancelled');
            testCase.verifyEqual(height(UMITRigStore.listRigs()), height(before), ...
                'Nothing is installed or overwritten.');
        end

        function testImportOfAnExistingRigAsCopyFindsFreeRigIDs(testCase)
            source = testCase.createLibraryRig('Dup');
            sourceInfo = source.getRigInfo();
            exportFile = testCase.exportRigFile(source);

            first = promptImportRig([], 'Interactive', false, ...
                'FilePath', exportFile, 'IfExists', 'copy');
            testCase.verifyEqual(first.severity, 'info', first.message);
            testCase.verifyTrue(first.changed);
            testCase.verifyNotEqual(first.rigUUID, sourceInfo.uuid, 'A copy has a new identity.');
            testCase.verifyEqual(first.rigID, [sourceInfo.rigID '_copy']);
            testCase.trackRig(first.rigUUID);

            second = promptImportRig([], 'Interactive', false, ...
                'FilePath', exportFile, 'IfExists', 'copy');
            testCase.verifyTrue(second.changed, second.message);
            testCase.verifyEqual(second.rigID, [sourceInfo.rigID '_copy2']);
            testCase.verifyNotEqual(second.rigUUID, first.rigUUID);
            testCase.trackRig(second.rigUUID);

            testCase.verifyTrue(UMITRigStore.open(first.rigUUID).validate('Mode', 'full').isValid);
            testCase.verifyEqual(UMITRigStore.open(sourceInfo.uuid).getRigInfo().rigID, ...
                sourceInfo.rigID, 'The original rig is unchanged.');
        end

        function testRigIDTakenByAnotherRigKeepsTheExportedIdentity(testCase)
            source = testCase.createLibraryRig('Clash');
            sourceInfo = source.getRigInfo();
            exportFile = testCase.exportRigFile(source);

            % Another rig (different UUID) now owns the exported rig ID.
            rmdir(source.RigRoot, 's');
            other = testCase.createLibraryRig('Other', sourceInfo.rigID);
            testCase.verifyNotEqual(other.getRigInfo().uuid, sourceInfo.uuid);

            refused = promptImportRig([], 'Interactive', false, 'FilePath', exportFile);
            testCase.verifyFalse(refused.changed);
            testCase.verifySubstring(refused.message, 'cancelled');

            renamed = promptImportRig([], 'Interactive', false, ...
                'FilePath', exportFile, 'IfExists', 'copy');
            testCase.verifyTrue(renamed.changed, renamed.message);
            testCase.verifyEqual(renamed.rigID, [sourceInfo.rigID '_copy']);
            testCase.verifyEqual(renamed.rigUUID, sourceInfo.uuid, ...
                'Only the rig ID changes when just the ID clashes.');
            testCase.trackRig(renamed.rigUUID);
        end

        function testImportFailuresAreReportedNotThrown(testCase)
            before = height(UMITRigStore.listRigs());

            notARig = fullfile(testCase.TempRoot, 'notARig.umitrig');
            somethingElse = 1;
            save(notARig, 'somethingElse', '-mat');
            result = promptImportRig([], 'Interactive', false, 'FilePath', notARig);
            testCase.verifyEqual(result.severity, 'error');
            testCase.verifyFalse(result.changed);
            testCase.verifySubstring(result.message, 'Could not import the rig');

            missing = promptImportRig([], 'Interactive', false, ...
                'FilePath', fullfile(testCase.TempRoot, 'missing.umitrig'));
            testCase.verifyEqual(missing.severity, 'error');
            testCase.verifyFalse(missing.changed);

            noPath = promptImportRig([], 'Interactive', false);
            testCase.verifyEqual(noPath.severity, 'error');

            testCase.verifyEqual(height(UMITRigStore.listRigs()), before, ...
                'A failed import installs nothing.');
        end

        function testExportFailuresAreReportedNotThrown(testCase)
            source = testCase.createLibraryRig('Exp');
            uuid = source.getRigInfo().uuid;

            badFolder = promptExportRig([], 'Interactive', false, 'RigUUID', uuid, ...
                'FilePath', fullfile(testCase.TempRoot, 'noSuchFolder', 'rig.umitrig'));
            testCase.verifyEqual(badFolder.severity, 'error');
            testCase.verifyEqual(badFolder.filePath, '');

            noInput = promptExportRig([], 'Interactive', false);
            testCase.verifyEqual(noInput.severity, 'error');

            unknownRig = promptExportRig([], 'Interactive', false, ...
                'RigUUID', char(java.util.UUID.randomUUID()), ...
                'FilePath', fullfile(testCase.TempRoot, 'unknown.umitrig'));
            testCase.verifyEqual(unknownRig.severity, 'error');
            testCase.verifyFalse(isfile(fullfile(testCase.TempRoot, 'unknown.umitrig')));

            unknownSet = promptExportRig([], 'Interactive', false, 'RigUUID', uuid, ...
                'FilePath', fullfile(testCase.TempRoot, 'unknownSet.umitrig'), ...
                'FilterSetIDs', {'NoSuchSet_dlg'});
            testCase.verifyEqual(unknownSet.severity, 'error');
            testCase.verifyFalse(isfile(fullfile(testCase.TempRoot, 'unknownSet.umitrig')));
        end

        function testChosenFilterSetsTravelWithTheRig(testCase)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            setID = ['DlgSet_' suffix(1:8)];
            UMITRigStore.importFilterSet(struct('id', setID, 'displayName', 'Dialog set', ...
                'excitationSpectrumID', 'FF0148017', 'emissionSpectrumID', 'FF01593LP'));
            testCase.addTeardown(@() iRemoveFilterSet(setID));
            source = testCase.createLibraryRig('FS');
            sourceInfo = source.getRigInfo();

            plain = promptExportRig([], 'Interactive', false, 'RigUUID', sourceInfo.uuid, ...
                'FilePath', fullfile(testCase.TempRoot, 'plain.umitrig'));
            testCase.verifyEmpty(plain.filterSetIDs);
            testCase.verifyFalse(contains(plain.message, 'filter set', 'IgnoreCase', true));

            withSet = promptExportRig([], 'Interactive', false, 'RigUUID', sourceInfo.uuid, ...
                'FilePath', fullfile(testCase.TempRoot, 'withSet.umitrig'), ...
                'FilterSetIDs', setID);
            testCase.verifyEqual(withSet.severity, 'info', withSet.message);
            testCase.verifyEqual(withSet.filterSetIDs, {setID});
            testCase.verifySubstring(withSet.message, setID);

            % Target library: the rig and the set are gone.
            rmdir(source.RigRoot, 's');
            UMITRigStore.removeFilterSet(setID);

            imported = promptImportRig([], 'Interactive', false, 'FilePath', withSet.filePath);
            testCase.verifyEqual(imported.severity, 'info', imported.message);
            testCase.verifyTrue(imported.changed);
            testCase.verifyEqual(imported.filterSetsAdded, {setID});
            testCase.verifySubstring(imported.message, setID);
            testCase.ExtraRigRoots{end+1} = UMITRigStore.open(imported.rigUUID).RigRoot;
            testCase.verifyEqual(UMITRigStore.getFilterSet(setID).emissionSpectrumID, 'FF01593LP');
        end

        function testDeleteRigNeedsAnArchivedRigAndTheTypedRigID(testCase)
        %TESTDELETERIGNEEDSANARCHIVEDRIGANDTHETYPEDRIGID promptDeleteRig deletes
        %only an archived rig and only when the rig ID is typed exactly.

            keep = testCase.createLibraryRig('Keep');
            target = testCase.createLibraryRig('Gone');
            if strcmp(target.getRigInfo().status, 'active')
                keep.activateRig();
            end
            info = target.getRigInfo();

            notArchived = promptDeleteRig([], info.uuid, 'Interactive', false, ...
                'ConfirmRigID', info.rigID);
            testCase.verifyEqual(notArchived.severity, 'warning');
            testCase.verifyFalse(notArchived.changed);
            testCase.verifyTrue(UMITRigStore.rigExists(info.uuid));

            target.archiveRig();
            wrongText = promptDeleteRig([], info.uuid, 'Interactive', false, ...
                'ConfirmRigID', 'wrong');
            testCase.verifyFalse(wrongText.changed);
            testCase.verifySubstring(wrongText.message, 'cancelled');
            testCase.verifyTrue(UMITRigStore.rigExists(info.uuid));
            noText = promptDeleteRig([], info.uuid, 'Interactive', false);
            testCase.verifyFalse(noText.changed);

            deleted = promptDeleteRig([], info.uuid, 'Interactive', false, ...
                'ConfirmRigID', info.rigID);
            testCase.verifyTrue(deleted.changed, deleted.message);
            testCase.verifyEqual(deleted.severity, 'success');
            testCase.verifyEqual(deleted.rigID, info.rigID);
            testCase.verifyFalse(UMITRigStore.rigExists(info.uuid));
            testCase.verifyTrue(UMITRigStore.rigExists(keep.getRigInfo().uuid));
            testCase.verifyTrue(UMITRigStore.validateLifecycle().isValid);

            unknown = promptDeleteRig([], char(java.util.UUID.randomUUID()), ...
                'Interactive', false, 'ConfirmRigID', 'x');
            testCase.verifyEqual(unknown.severity, 'error');
        end

        function testConflictingFilterSetFailsTheImportWithoutInstallingTheRig(testCase)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            setID = ['DlgConflict_' suffix(1:8)];
            UMITRigStore.importFilterSet(struct('id', setID, 'displayName', 'Original', ...
                'excitationSpectrumID', 'FF0148017', 'emissionSpectrumID', 'FF01593LP'));
            testCase.addTeardown(@() iRemoveFilterSet(setID));
            source = testCase.createLibraryRig('FSC');
            exported = promptExportRig([], 'Interactive', false, ...
                'RigUUID', source.getRigInfo().uuid, ...
                'FilePath', fullfile(testCase.TempRoot, 'conflict.umitrig'), ...
                'FilterSetIDs', setID);
            testCase.assertEqual(exported.severity, 'info', exported.message);

            rmdir(source.RigRoot, 's');
            UMITRigStore.removeFilterSet(setID);
            UMITRigStore.importFilterSet(struct('id', setID, 'displayName', 'Other', ...
                'excitationSpectrumID', 'none', 'emissionSpectrumID', 'none'));
            before = height(UMITRigStore.listRigs());

            result = promptImportRig([], 'Interactive', false, 'FilePath', exported.filePath);
            testCase.verifyEqual(result.severity, 'error');
            testCase.verifyFalse(result.changed);
            testCase.verifySubstring(result.message, setID);
            testCase.verifyEqual(height(UMITRigStore.listRigs()), before);
        end
    end

    methods (Access = private)
        function store = createLibraryRig(testCase, label, rigID)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            if nargin < 3
                rigID = ['Transfer' label '_' suffix(1:8)];
            end
            store = UMITRigStore.create(struct( ...
                'rigID', rigID, ...
                'displayName', ['Transfer rig ' label], ...
                'cameras', struct('index', 1, 'displayName', 'Cam', 'manufacturer', '', ...
                    'model', 'D1024', 'serialNumber', '', 'spectrumID', 'PF1024'), ...
                'illuminations', struct('name', 'red', 'displayName', 'Red', ...
                    'manufacturer', '', 'model', '', 'spectrumID', 'LED_632nm')));
            testCase.ExtraRigRoots{end+1} = store.RigRoot;
        end

        function filePath = exportRigFile(testCase, store)
            result = promptExportRig([], 'Interactive', false, ...
                'RigUUID', store.getRigInfo().uuid, ...
                'FilePath', fullfile(testCase.TempRoot, [store.getRigInfo().rigID '.umitrig']));
            testCase.assertEqual(result.severity, 'info', result.message);
            filePath = result.filePath;
        end

        function trackRig(testCase, rigUUID)
            testCase.ExtraRigRoots{end+1} = UMITRigStore.open(rigUUID).RigRoot;
        end
    end
end

function iRemoveFilterSet(filterSetID)
try
    UMITRigStore.removeFilterSet(filterSetID);
catch
    % Already removed by the test.
end
end

function iRemoveFolder(folder)
if ~isempty(folder) && isfolder(folder)
    rmdir(folder, 's');
end
end
