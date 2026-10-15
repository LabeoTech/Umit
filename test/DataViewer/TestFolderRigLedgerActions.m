classdef TestFolderRigLedgerActions < matlab.unittest.TestCase
    %TESTFOLDERRIGLEDGERACTIONS DataViewer folder-rig backend (Phase 5).
    %
    %   Covers describeFolderRigLedger (display), planFolderRigBinding
    %   (read-only preview incl. conflicts), applyFolderRigBinding (writes)
    %   and formatRigLedgerHistory. The dialogs that wrap them
    %   (runFolderRigBindingDialog, promptAcknowledgeRigConflict,
    %   showRigLedgerHistoryDialog) are checked manually.

    properties
        Folder
        RigA
        RigB
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.Folder = tempname;
            mkdir(testCase.Folder);
            testCase.addTeardown(@() rmdir(testCase.Folder, 's'));
            testCase.RigA = testCase.createLibraryRig('A', 'LED_632nm');
            testCase.RigB = testCase.createLibraryRig('B', 'LED_593nm');
        end
    end

    methods (Test)
        function testDescribeEveryState(testCase)
            summary = describeFolderRigLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Unknown');
            testCase.verifySubstring(summary.displayText, 'no rig ledger');
            testCase.verifyTrue(summary.canBind);
            testCase.verifyFalse(summary.canViewHistory);
            testCase.verifyFalse(UMITRigLedger.exists(testCase.Folder), 'Describing never writes.');

            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            testCase.verifyEqual(describeFolderRigLedger(testCase.Folder).displayText, 'none (no rig)');

            applyFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigA));
            summary = describeFolderRigLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'BoundUncommitted');
            testCase.verifyEqual(summary.displayText, testCase.RigA.getRigInfo().displayName);

            testCase.commitBoundOptics();
            summary = describeFolderRigLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'BoundCommitted');
            testCase.verifySubstring(summary.displayText, '(committed)');
            testCase.verifyNumElements(summary.committedLines, 1);

            applyFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigB));
            summary = describeFolderRigLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Conflicted');
            testCase.verifySubstring(summary.displayText, 'CONFLICT');
            testCase.verifyTrue(summary.canAcknowledge);
            testCase.verifyEqual(summary.conflictParts, {'optics'});

            UMITRigLedger.acknowledgeConflict(testCase.Folder, 'Origin', 'test');
            summary = describeFolderRigLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'BoundCommitted');
            testCase.verifyTrue(summary.hadConflict);
            testCase.verifySubstring(summary.displayText, '[had conflict]');
            testCase.verifyFalse(summary.canAcknowledge);
        end

        function testDescribeInvalidLedgerDoesNotThrow(testCase)
            testCase.writeText(UMITRigLedger.filePath(testCase.Folder), '{not json');
            summary = describeFolderRigLedger(testCase.Folder);
            testCase.verifyEqual(summary.state, 'Invalid');
            testCase.verifyFalse(summary.canBind);
            testCase.verifyNotEmpty(summary.issues);
            testCase.verifyEqual(describeFolderRigLedger('').displayText, '<none>');
        end

        function testPlanIsReadOnlyAndPreviewsConflicts(testCase)
            plan = planFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigA));
            testCase.verifyEqual(plan.action, 'bind');
            testCase.verifyFalse(UMITRigLedger.exists(testCase.Folder));

            applyFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigA));
            testCase.verifyEqual(planFolderRigBinding(testCase.Folder, ...
                testCase.uuid(testCase.RigA)).action, 'none');

            plan = planFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigB));
            testCase.verifyEqual(plan.action, 'rebind');
            testCase.verifyFalse(plan.conflict.hasConflict, 'Uncommitted: rebinding is free.');

            testCase.commitBoundOptics();
            before = fileread(UMITRigLedger.filePath(testCase.Folder));
            plan = planFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigB));
            testCase.verifyTrue(plan.conflict.hasConflict);
            testCase.verifyEqual(plan.conflict.parts, {'optics'});
            testCase.verifySubstring(plan.confirmText, 'CONFLICTED');
            testCase.verifySubstring(plan.confirmText, testCase.RigA.getRigInfo().rigID);
            testCase.verifySubstring(plan.confirmText, testCase.RigB.getRigInfo().rigID);
            testCase.verifyEqual(fileread(UMITRigLedger.filePath(testCase.Folder)), before, ...
                'Previewing never writes.');
        end

        function testRebindKeepsTheFolderCoregistration(testCase)
            % The transform applied to the data only changes by re-import.
            infoA = testCase.RigA.getRigInfo();
            coreg = UMITRigLedger.coregistrationSnapshotFromTransform( ...
                'resource-1', 'abc', affine2d([1 0 0; 0 1 0; 4 2 1]));
            UMITRigLedger.bind(testCase.Folder, UMITRigLedger.makeBinding( ...
                UMITRigLedger.sourceRigFromRigInfo(infoA), ...
                UMITRigLedger.buildOpticsSnapshot(infoA), coreg), 'Origin', 'test');
            UMITRigLedger.commit(testCase.Folder, 'coregistration', coreg.hash, 'Origin', 'test');

            result = applyFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigB));

            info = UMITRigLedger.read(testCase.Folder);
            testCase.verifyEqual(result.action, 'rebind');
            testCase.verifyFalse(result.conflict.hasConflict);
            testCase.verifyEqual(info.binding.coregistrationHash, coreg.hash);
            testCase.verifyEqual(info.binding.sourceRig.uuid, testCase.uuid(testCase.RigB));
            testCase.verifyEqual(info.state, 'BoundCommitted');
        end

        function testArchivedRigCannotBeBound(testCase)
            testCase.RigB.archiveRig();
            testCase.verifyError(@() planFolderRigBinding(testCase.Folder, ...
                testCase.uuid(testCase.RigB)), 'Umitoolbox:DataViewer:archivedRig');
            testCase.verifyError(@() applyFolderRigBinding(testCase.Folder, ...
                testCase.uuid(testCase.RigB)), 'Umitoolbox:DataViewer:archivedRig');
            testCase.verifyFalse(UMITRigLedger.exists(testCase.Folder));
        end

        function testHistoryListsEveryEvent(testCase)
            testCase.verifySubstring(strjoin(formatRigLedgerHistory(testCase.Folder)), ...
                'no rig ledger');

            applyFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigA));
            testCase.commitBoundOptics();
            applyFolderRigBinding(testCase.Folder, testCase.uuid(testCase.RigB));
            UMITRigLedger.acknowledgeConflict(testCase.Folder, 'Origin', 'DataViewer', ...
                'Comment', 'LED swapped before acquisition');

            text = strjoin(formatRigLedgerHistory(testCase.Folder), newline);
            for type = ["bind", "commit", "rebind", "conflict", "acknowledgeConflict"]
                testCase.verifySubstring(text, char(type));
            end
            testCase.verifySubstring(text, 'LED swapped before acquisition');
            testCase.verifySubstring(text, 'had a rig conflict');
            testCase.verifySubstring(text, testCase.RigB.getRigInfo().rigID);
        end
    end

    methods (Access = private)
        function store = createLibraryRig(testCase, label, redSpectrumID)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            store = UMITRigStore.create(struct( ...
                'rigID', ['FolderRig' label '_' suffix(1:8)], ...
                'displayName', ['Folder rig ' label], ...
                'cameras', struct('index', 1, 'displayName', 'Cam', 'manufacturer', '', ...
                    'model', 'D1024', 'serialNumber', '', 'spectrumID', 'PF1024'), ...
                'illuminations', struct('name', 'red', 'displayName', 'Red', ...
                    'manufacturer', '', 'model', '', 'spectrumID', redSpectrumID)));
            rigRoot = store.RigRoot;
            testCase.addTeardown(@() iRemoveFolder(rigRoot));
        end

        function value = uuid(~, store)
            value = store.getRigInfo().uuid;
        end

        function commitBoundOptics(testCase)
            info = UMITRigLedger.read(testCase.Folder);
            UMITRigLedger.commit(testCase.Folder, 'optics', info.binding.opticsHash, ...
                'Origin', 'test');
        end

        function writeText(~, path, text)
            fid = fopen(path, 'w');
            fprintf(fid, '%s', text);
            fclose(fid);
        end
    end
end

function iRemoveFolder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
