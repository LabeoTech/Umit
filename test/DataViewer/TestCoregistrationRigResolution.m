classdef TestCoregistrationRigResolution < matlab.unittest.TestCase
    %TESTCOREGISTRATIONRIGRESOLUTION Target rig of DataViewer_Coreg2Cams.
    %
    %   Rig-ledger redesign, Phase 4: the calibration target is the rig bound
    %   in the folder's rig ledger, else a caller-supplied rig, else a user
    %   choice with the Active rig pre-selected -- never the machine default.
    %   Covers resolveCoregistrationRig, buildLibraryRigChoices and that the
    %   promptSelectLibraryRig dialog is closed on return (the dialog's look
    %   is checked manually).

    properties
        Folder
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            testCase.Folder = tempname;
            mkdir(testCase.Folder);
            testCase.addTeardown(@() rmdir(testCase.Folder, 's'));
        end
    end

    methods (Test)
        function testUnboundFolderNeedsSelection(testCase)
            target = resolveCoregistrationRig(testCase.Folder);
            testCase.verifyTrue(target.needsSelection);
            testCase.verifyEmpty(target.rigUUID);
            testCase.verifyNotEmpty(target.message);

            UMITRigLedger.writeNoRig(testCase.Folder, 'Origin', 'test');
            testCase.verifyTrue(resolveCoregistrationRig(testCase.Folder).needsSelection);
            testCase.verifyTrue(resolveCoregistrationRig('').needsSelection);
        end

        function testBoundFolderUsesItsRigOverCallerRig(testCase)
            bound = testCase.createLibraryRig('Bound');
            other = testCase.createLibraryRig('Other');
            UMITRigLedger.bind(testCase.Folder, UMITRigLedger.bindingFromRig(bound, ...
                'IncludeCoregistration', false), 'Origin', 'test');

            target = resolveCoregistrationRig(testCase.Folder, ...
                'RigUUID', other.getRigInfo().uuid);

            testCase.verifyFalse(target.needsSelection);
            testCase.verifyEqual(target.origin, 'folderLedger');
            testCase.verifyEqual(target.rigUUID, bound.getRigInfo().uuid);
        end

        function testCallerRigUsedWhenFolderUnbound(testCase)
            other = testCase.createLibraryRig('Caller');
            target = resolveCoregistrationRig(testCase.Folder, ...
                'RigUUID', other.getRigInfo().uuid);
            testCase.verifyEqual(target.origin, 'caller');
            testCase.verifyEqual(target.rigID, other.getRigInfo().rigID);
        end

        function testArchivedOrMissingBoundRigIsAnError(testCase)
            bound = testCase.createLibraryRig('Archived');
            UMITRigLedger.bind(testCase.Folder, UMITRigLedger.bindingFromRig(bound, ...
                'IncludeCoregistration', false), 'Origin', 'test');
            bound.archiveRig();
            testCase.verifyError(@() resolveCoregistrationRig(testCase.Folder), ...
                'Umitoolbox:Coreg2Cams:boundRigArchived');

            rmdir(bound.RigRoot, 's');
            testCase.verifyError(@() resolveCoregistrationRig(testCase.Folder), ...
                'Umitoolbox:Coreg2Cams:boundRigNotInLibrary');
        end

        function testInvalidLedgerIsAnError(testCase)
            fid = fopen(UMITRigLedger.filePath(testCase.Folder), 'w');
            fprintf(fid, '{not json');
            fclose(fid);
            testCase.verifyError(@() resolveCoregistrationRig(testCase.Folder), ...
                'Umitoolbox:RigLedger:invalidLedger');
        end

        function testLibraryRigChoicesPreselection(testCase)
            rigs = table( ...
                ["u1"; "u2"; "u3"; "u4"], ["R1"; "R2"; "R3"; "R4"], ...
                ["One"; ""; "Three"; "Four"], false(4, 1), ...
                ["available"; "active"; "archived"; "available"], ...
                [true; true; true; false], strings(4, 1), ...
                'VariableNames', {'RigUUID', 'RigID', 'DisplayName', 'IsDefault', ...
                'Status', 'IsReadable', 'RigRoot'});

            [labels, uuids, defaultIndex] = buildLibraryRigChoices(rigs);
            testCase.verifyEqual(uuids, {'u1', 'u2'});
            testCase.verifyEqual(labels, {'One (R1)', 'R2 - Active'});
            testCase.verifyEqual(defaultIndex, 2);

            [~, ~, defaultIndex] = buildLibraryRigChoices(rigs, 'u1');
            testCase.verifyEqual(defaultIndex, 1);
            [~, ~, defaultIndex] = buildLibraryRigChoices(rigs, 'u3');
            testCase.verifyEqual(defaultIndex, 2, 'Archived rigs are never offered.');

            [labels, uuids, defaultIndex] = buildLibraryRigChoices(rigs(3:4, :));
            testCase.verifyEmpty(labels);
            testCase.verifyEmpty(uuids);
            testCase.verifyEqual(defaultIndex, 0);
        end

        function testPickerDialogIsClosedOnReturn(testCase)
            % Regression (2026-10-15): the modal picker stayed open after OK
            % and hid DataViewer's Rebind confirmation.
            testCase.createLibraryRig('Picker');
            title = ['Picker test ' strrep(char(java.util.UUID.randomUUID()), '-', '')];
            clicker = timer('StartDelay', 1, 'TimerFcn', @(~, ~) iPressButton(title, 'OK'));
            testCase.addTeardown(@() delete(clicker));
            start(clicker);

            rigUUID = promptSelectLibraryRig([], 'Title', title);

            testCase.verifyNotEmpty(rigUUID);
            testCase.verifyEmpty(findall(groot, 'Type', 'figure', 'Name', title), ...
                'The picker dialog must be deleted when promptSelectLibraryRig returns.');
        end
    end

    methods (Access = private)
        function store = createLibraryRig(testCase, label)
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            store = UMITRigStore.create(struct( ...
                'rigID', ['CoregTarget' label '_' suffix(1:8)], ...
                'displayName', ['Coreg target ' label], ...
                'cameras', struct('index', 1, 'displayName', 'Cam', 'manufacturer', '', ...
                    'model', 'D1024', 'serialNumber', '', 'spectrumID', 'PF1024'), ...
                'illuminations', struct('name', 'red', 'displayName', 'Red', ...
                    'manufacturer', '', 'model', '', 'spectrumID', 'LED_632nm')));
            rigRoot = store.RigRoot;
            testCase.addTeardown(@() iRemoveFolder(rigRoot));
        end
    end
end

function iPressButton(figureName, buttonText)
dialog = findall(groot, 'Type', 'figure', 'Name', figureName);
if isempty(dialog)
    return
end
button = findall(dialog(1), 'Type', 'uibutton', 'Text', buttonText);
button(1).ButtonPushedFcn(button(1), []);
end

function iRemoveFolder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
