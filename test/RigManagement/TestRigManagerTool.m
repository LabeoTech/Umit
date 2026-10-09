classdef TestRigManagerTool < matlab.unittest.TestCase
%TESTRIGMANAGERTOOL Focused non-interactive Rig Manager GUI tests.

    properties
        TempRoot char
        RigRoots cell
        RigUUIDs cell
        Apps cell
        RigActivationFixture = struct()
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.RigRoots = {};
            testCase.RigUUIDs = {};
            testCase.Apps = {};
            testCase.RigActivationFixture = struct();
        end
    end

    methods (TestMethodTeardown)
        function removeFixture(testCase)
            for iApp = 1:numel(testCase.Apps)
                app = testCase.Apps{iApp};
                if ~isempty(app) && isvalid(app)
                    delete(app);
                end
            end
            deactivateRigTemporarily(testCase.RigActivationFixture);
            rigs = UMITRigStore.listRigs();
            for iRig = 1:numel(testCase.RigUUIDs)
                row = find(strcmpi(rigs.RigUUID, ...
                    string(testCase.RigUUIDs{iRig})), 1);
                if ~isempty(row)
                    rigRoot = char(rigs.RigRoot(row));
                    if isfolder(rigRoot)
                        rmdir(rigRoot, 's');
                    end
                end
            end
            for iRig = 1:numel(testCase.RigRoots)
                if isfolder(testCase.RigRoots{iRig})
                    rmdir(testCase.RigRoots{iRig}, 's');
                end
            end
            if isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function testStartupContractAndRealTabs(testCase)
            store = testCase.createRig(1, 'Startup');
            info = store.getRigInfo();
            app = testCase.openApp(info.uuid);

            testCase.verifyEqual(app.SelectedRigUUID, info.uuid);
            testCase.verifyEqual(app.SelectedRigID, info.rigID);
            testCase.verifyEqual(numel(app.RigTabGroup.Children), 3);
            testCase.verifyEqual(sort(string( ...
                {app.RigTabGroup.Children.Title})), ...
                sort(["Overview", "Hardware", "Coregistration"]));
            testCase.verifyEqual(string(app.UseSelectedRigButton.Visible), "off");
        end

        function testSingleCameraKeepsCamera2Optional(testCase)
            store = testCase.createRig(1, 'Single');
            app = testCase.openApp(store.getRigInfo().uuid);

            testCase.verifyEqual(string(app.Camera2Panel.Visible), "off");
            testCase.verifyEqual(string(app.AddCamera2Button.Visible), "on");
            testCase.verifySubstring(app.CoregistrationMessageLabel.Text, ...
                'only applies to Rigs with Camera 2');
            testCase.verifyEqual(string(app.ImportCoregistrationButton.Enable), "off");
        end

        function testHardwareEditingIsStagedAndCanonical(testCase)
            store = testCase.createRig(1, 'Hardware');
            app = testCase.openApp(store.getRigInfo().uuid);
            originalModel = app.Camera1ModelField.Value;

            testCase.verifyEqual(app.IlluminationTable.Data(:, 1), ...
                {'Red'; 'Green'; 'Yellow'});
            testCase.verifyFalse(any(app.IlluminationTable.ColumnEditable));

            testCase.push(app.EditHardwareButton);
            testCase.verifyEqual(app.IlluminationTable.ColumnEditable, ...
                [false true true true false]);
            testCase.verifyEqual(string(app.Camera1ModelField.Editable), "on");
            testCase.verifyEqual(string(app.Camera1SpectrumField.Editable), "off");
            app.Camera1ModelField.Value = 'Staged Only';

            persisted = store.getRigInfo();
            testCase.verifyEqual(persisted.cameras(1).model, originalModel);
            testCase.push(app.RevertHardwareButton);
            testCase.verifyEqual(app.Camera1ModelField.Value, originalModel);
            testCase.verifyEqual(string(app.Camera1ModelField.Editable), "off");
        end

        function testArchivedRigIsReadOnlyAndRestorable(testCase)
            store = testCase.createRig(1, 'Archived');
            info = store.getRigInfo();
            replacement = testCase.createRig(1, 'ArchiveReplacement');
            store.archiveRig(replacement.getRigInfo().uuid);
            app = testCase.openApp(info.uuid);

            testCase.verifyEqual(string(app.DisplayNameField.Editable), "off");
            testCase.verifyEqual(string(app.EditHardwareButton.Enable), "off");
            testCase.verifyEqual(string(app.ArchiveRigButton.Visible), "off");
            testCase.verifyEqual(string(app.RestoreRigButton.Visible), "on");
            testCase.verifyEqual(string(app.RestoreRigButton.Enable), "on");
            testCase.verifySubstring(app.RigStatusLabel.Text, 'Read-only');
        end

        function testInvalidRigOpensReadOnlyWithDiagnostics(testCase)
            store = testCase.createRig(1, 'Invalid');
            info = store.getRigInfo();
            schema = getUMITRigSchema();
            rigFile = fullfile(store.RigRoot, schema.files.rigMetadata);
            loaded = load(rigFile, schema.metadataVariables.rig, '-mat');
            RigInfo = loaded.RigInfo;
            RigInfo.rigID = ['Mismatched_' info.rigID];
            save(rigFile, 'RigInfo', '-mat');

            app = testCase.openApp(info.uuid);
            testCase.verifyEqual(string(app.DisplayNameField.Editable), "off");
            testCase.verifyEqual(string(app.EditHardwareButton.Enable), "off");
            testCase.verifySubstring(app.RigStatusLabel.Text, 'Invalid');
            testCase.verifyTrue(any(contains(string(app.AdvancedArea.Value), ...
                'invalidRigIdentity', 'IgnoreCase', true)));
        end

        function testDefaultIndicatorAndPublicShowRig(testCase)
            first = testCase.createRig(1, 'First');
            second = testCase.createRig(2, 'Second');
            testCase.RigActivationFixture = ...
                activateRigTemporarily(first.getRigInfo().uuid);
            app = testCase.openApp(second.getRigInfo().uuid);

            testCase.verifyTrue(any(app.RigListTable.Data.Default == "Default"));
            app.showRig(first.getRigInfo().rigID);
            testCase.verifyEqual(app.SelectedRigUUID, first.getRigInfo().uuid);
            testCase.verifyEqual(app.DefaultStatusField.Value, 'Yes');
        end

        function testOpticalLibraryIsSeparateWorkspace(testCase)
            store = testCase.createRig(1, 'Library');
            app = RigManagerTool([], 'manage', ...
                'InitialRig', store.getRigInfo().uuid, ...
                'StartWorkspace', 'opticalLibrary', 'Visible', false);
            testCase.Apps{end+1} = app;

            testCase.verifyEqual(string(app.RigWorkspaceGrid.Visible), "off");
            testCase.verifyEqual(string(app.LibraryWorkspaceGrid.Visible), "on");
            testCase.verifyGreaterThan(height(app.LibraryTable.Data), 0);
            testCase.verifyTrue(ismember('Origin', ...
                app.LibraryTable.Data.Properties.VariableNames));
            testCase.verifyEqual(string(app.LibraryRemoveButton.Enable), "off");
        end

        function testSelectModeExposesSelectionContract(testCase)
            store = testCase.createRig(1, 'Select');
            app = RigManagerTool([], 'select', ...
                'InitialRig', store.getRigInfo().uuid, 'Visible', false);
            testCase.Apps{end+1} = app;

            testCase.verifyEqual(string(app.UseSelectedRigButton.Visible), "on");
            testCase.verifyEqual(app.UseSelectedRigButton.Text, 'Use Selected Rig');
            testCase.verifyFalse(app.WasSelectionConfirmed);
            testCase.verifyEqual(app.OutputRigUUID, '');
        end

        function testAcceptedLayoutAndRequiredCameraControls(testCase)
            store = testCase.createRig(1, 'Layout');
            app = testCase.openApp(store.getRigInfo().uuid);

            testCase.verifyEqual(app.RigWorkspaceGrid.RowHeight{1}, 175);
            testCase.verifyClass(app.RigListTable, 'matlab.ui.control.Table');
            testCase.verifyFalse(isprop(app, 'DuplicateForHardwareButton'));
            testCase.verifyEqual(string(app.OverviewTab.Children.Scrollable), "on");
            testCase.verifyEqual(string(app.HardwareTab.Children.Scrollable), "on");
            testCase.verifyEqual(string(app.CoregistrationTab.Children.Scrollable), "on");
            cameraLabels = findall(app.Camera1Panel, 'Type', 'uilabel');
            testCase.verifyTrue(any(contains(string({cameraLabels.Text}), ...
                'Spectral Profile *')));
            testCase.verifyClass(app.RemoveIlluminationSpectrumButton, ...
                'matlab.ui.control.Button');
        end
    end

    methods (Access = private)
        function store = createRig(testCase, cameraCount, label)
            definition = UMITRigStore.getBuiltInRigDefinition('OiS200');
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            definition.rigID = sprintf('RigManager_%s_%s', ...
                label, suffix(1:8));
            definition.displayName = sprintf('Rig Manager %s', label);
            definition.cameras = definition.cameras(1:cameraCount);
            store = UMITRigStore.create(definition);
            testCase.RigRoots{end+1} = store.RigRoot;
            testCase.RigUUIDs{end+1} = store.getRigInfo().uuid;
        end

        function app = openApp(testCase, rigUUID)
            app = RigManagerTool([], 'manage', ...
                'InitialRig', rigUUID, 'Visible', false);
            testCase.Apps{end+1} = app;
        end

        function push(~, button)
            callback = button.ButtonPushedFcn;
            callback(button, []);
            drawnow;
        end
    end
end
