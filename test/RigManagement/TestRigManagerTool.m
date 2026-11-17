classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestRigManagerTool < matlab.unittest.TestCase
%TESTRIGMANAGERTOOL Focused non-interactive Rig Manager GUI tests.

    properties
        TempRoot char
        RigRoots cell
        RigUUIDs cell
        Apps cell
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.RigRoots = {};
            testCase.RigUUIDs = {};
            testCase.Apps = {};
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
            % A test may have made one of its (now removed) rigs the Default
            % rig; do not leave a dangling pointer in the shared temp root.
            schema = getUMITRigSchema();
            defaultFile = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.store.internalFolder, schema.store.defaultFile);
            if isfile(defaultFile)
                delete(defaultFile);
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
            UMITRigStore.setDefaultRig(first.getRigInfo().uuid);
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

        function testCombineButtonForSpectraNeedsSelectionAndHidesForFilterSets(testCase)
            store = testCase.createRig(1, 'Combine');
            app = RigManagerTool([], 'manage', ...
                'InitialRig', store.getRigInfo().uuid, ...
                'StartWorkspace', 'opticalLibrary', 'Visible', false);
            testCase.Apps{end+1} = app;

            expectedVisible = struct('illumination', "on", 'filter', "on", ...
                'camera', "on", 'filterSets', "off");
            for category = fieldnames(expectedVisible)'
                app.LibraryCategoryDropDown.Value = category{1};
                app.LibraryCategoryDropDown.ValueChangedFcn([], []);
                testCase.verifyEqual(string(app.LibraryCombineButton.Visible), ...
                    expectedVisible.(category{1}), category{1});
                if expectedVisible.(category{1}) == "on"
                    % Combining starts from the selected table item.
                    testCase.verifyEqual(string(app.LibraryCombineButton.Enable), ...
                        "off", category{1});
                    app.LibraryTable.CellSelectionCallback([], struct('Indices', [1 1]));
                    testCase.verifyEqual(string(app.LibraryCombineButton.Enable), ...
                        "on", category{1});
                end
            end
        end

        function testLibrarySelectionListsCompositionOfCombinedSpectrum(testCase)
            store = testCase.createRig(1, 'CombineInfo');
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            compositeID = ['CombineInfo_' suffix(1:8)];
            UMITRigStore.createCompositeSpectrum('filter', compositeID, ...
                struct('category', {'filter'; 'filter'}, ...
                'spectrumID', {'FF0154080'; 'FF0151444'}), struct());
            testCase.addTeardown(@() UMITRigStore.removeSpectrum('filter', compositeID));
            app = RigManagerTool([], 'manage', ...
                'InitialRig', store.getRigInfo().uuid, ...
                'StartWorkspace', 'opticalLibrary', 'Visible', false);
            testCase.Apps{end+1} = app;

            app.LibraryCategoryDropDown.Value = 'filter';
            app.LibraryCategoryDropDown.ValueChangedFcn([], []);
            row = find(app.LibraryTable.Data.ID == string(compositeID));
            testCase.assertNotEmpty(row);
            app.LibraryTable.CellSelectionCallback([], struct('Indices', [row 1]));

            text = strjoin(string(app.LibraryMetadataArea.Value), newline);
            testCase.verifySubstring(text, "Combined from");
            testCase.verifySubstring(text, "FF0154080 x FF0151444");
            testCase.verifyEqual(string(app.LibraryRemoveButton.Enable), "on");
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

        function testNewRigFromTemplateCreatesSingleCameraRig(testCase)
        %TESTNEWRIGFROMTEMPLATECREATESSINGLECAMERARIG Drive the real New Rig
        %button with the OiS200 template; the modal dialogs are replaced by
        %non-blocking stubs. The single-camera template must create a
        %one-camera Rig without any Camera 2 prompt, and Camera 2 must stay
        %addable afterwards.

            host = testCase.createRig(1, 'NewRigHost');
            app = testCase.openApp(host.getRigInfo().uuid);
            newRigID = ['RigManager_FromTemplate_' ...
                char(extractBefore(string(java.util.UUID.randomUUID()), 9))];
            testCase.installDialogStubs(newRigID);

            testCase.push(app.NewRigButton);

            titles = getappdata(0, 'TestRigManagerDialogTitles');
            testCase.verifyEqual(titles, {'Active Rig', 'Review New Rig'});
            rigs = UMITRigStore.listRigs();
            row = find(strcmp(rigs.RigID, newRigID), 1);
            testCase.assertNotEmpty(row, 'The template Rig was not created.');
            testCase.RigUUIDs{end+1} = char(rigs.RigUUID(row));
            info = UMITRigStore.open(char(rigs.RigUUID(row))).getRigInfo();

            testCase.verifyEqual(numel(info.cameras), 1);
            testCase.verifyEqual(info.cameras.manufacturer, 'Thorlabs');
            testCase.verifyEqual(info.cameras.model, 'CS2100M');
            testCase.verifyEqual(info.cameras.spectrumID, 'ThorQLux');
            testCase.verifyEqual({info.illuminations.name}, ...
                {'red', 'green', 'yellow'});
            testCase.verifyEqual(app.SelectedRigID, newRigID);
            testCase.verifyEqual(string(app.Camera2Panel.Visible), "off");
            testCase.verifyEqual(string(app.AddCamera2Button.Visible), "on");
        end
    end

    methods (Access = private)
        function store = createRig(testCase, cameraCount, label)
            definition = UMITRigStore.getBuiltInRigDefinition('OiS200');
            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            definition.rigID = sprintf('RigManager_%s_%s', ...
                label, suffix(1:8));
            definition.displayName = sprintf('Rig Manager %s', label);
            % The built-in template is single-camera; replicate Camera 1
            % as Camera 2 when a dual-camera Rig is requested.
            cameras = definition.cameras(1);
            for iCamera = 2:cameraCount
                cameras(iCamera) = definition.cameras(1);
                cameras(iCamera).index = iCamera;
                cameras(iCamera).displayName = sprintf( ...
                    '%s %d', definition.cameras(1).displayName, iCamera);
            end
            definition.cameras = cameras;
            store = UMITRigStore.create(definition);
            testCase.RigRoots{end+1} = store.RigRoot;
            testCase.RigUUIDs{end+1} = store.getRigInfo().uuid;
        end

        function app = openApp(testCase, rigUUID)
            app = RigManagerTool([], 'manage', ...
                'InitialRig', rigUUID, 'Visible', false);
            testCase.Apps{end+1} = app;
        end

        function installDialogStubs(testCase, newRigID)
            %INSTALLDIALOGSTUBS Shadow the modal dialogs used by onNewRig.
            stubFolder = tempname;
            mkdir(stubFolder);
            testCase.addTeardown(@() rmdir(stubFolder, 's'));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(stubFolder));
            testCase.writeStub(stubFolder, 'listdlg', { ...
                'function [selection, ok] = listdlg(varargin)', ...
                'selection = 1; ok = true;', 'end'});
            testCase.writeStub(stubFolder, 'inputdlg', { ...
                'function answer = inputdlg(~, ~, ~, definput)', ...
                'answer = definput;', ...
                'answer{1} = getappdata(0, ''TestRigManagerNewRigID'');', ...
                'end'});
            testCase.writeStub(stubFolder, 'uiconfirm', { ...
                'function answer = uiconfirm(~, ~, titleText, varargin)', ...
                'titles = getappdata(0, ''TestRigManagerDialogTitles'');', ...
                'setappdata(0, ''TestRigManagerDialogTitles'', [titles, {titleText}]);', ...
                'switch titleText', ...
                '    case ''Active Rig'', answer = ''No'';', ...
                '    case ''Review New Rig'', answer = ''Create Rig'';', ...
                '    otherwise, error(''TestRigManagerTool:unexpectedDialog'', ''Unexpected dialog: %s'', titleText);', ...
                'end', 'end'});
            setappdata(0, 'TestRigManagerNewRigID', newRigID);
            setappdata(0, 'TestRigManagerDialogTitles', {});
            testCase.addTeardown(@() testCase.clearStubAppData());
        end

        function clearStubAppData(~)
            for name = {'TestRigManagerNewRigID', 'TestRigManagerDialogTitles'}
                if isappdata(0, name{1})
                    rmappdata(0, name{1});
                end
            end
        end

        function writeStub(~, folder, name, lines)
            fid = fopen(fullfile(folder, [name '.m']), 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end

        function push(~, button)
            callback = button.ButtonPushedFcn;
            callback(button, []);
            drawnow;
        end
    end
end
