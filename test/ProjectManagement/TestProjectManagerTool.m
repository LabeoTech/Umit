classdef TestProjectManagerTool < matlab.unittest.TestCase
    %TESTPROJECTMANAGERTOOL Integration tests for ProjectManagerTool.mlapp.

    properties
        App
        Store
        ProjectRoot = ''
        SaveFolder = ''
        SubjectUUID = ''
        SessionUUID = ''
    end

    methods (TestMethodSetup)
        function createProjectAndApp(testCase)
            projectName = sprintf( ...
                'ProjectManagerToolTest_%s', ...
                char(java.util.UUID.randomUUID));
            testCase.Store = UMITProjectStore.create(struct( ...
                'projectName', projectName, ...
                'description', 'Temporary GUI integration test.'));
            testCase.ProjectRoot = testCase.Store.ProjectRoot;

            testCase.SubjectUUID = testCase.Store.addSubject(struct( ...
                'subjectID', 'Subject_01', ...
                'displayName', 'Subject One'));
            testCase.SaveFolder = tempname;
            mkdir(testCase.SaveFolder);
            AcqInfoStream = struct();
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
            testCase.SessionUUID = testCase.Store.addSession( ...
                'Subject_01', struct( ...
                'sessionID', 'Session_01', ...
                'displayName', 'Session One', ...
                'processedDataFolder', testCase.SaveFolder));

            testCase.App = ProjectManagerTool(testCase.Store);
            drawnow
        end
    end

    methods (TestMethodTeardown)
        function removeProjectAndApp(testCase)
            if ~isempty(testCase.App) && isvalid(testCase.App)
                delete(testCase.App);
            end

            if ~isempty(testCase.ProjectRoot) && ...
                    isfolder(testCase.ProjectRoot)
                projectsRoot = UMITProjectStore.getProjectsRoot();
                projectParent = fileparts(testCase.ProjectRoot);
                if strcmpi(projectParent, projectsRoot)
                    rmdir(testCase.ProjectRoot, 's');
                end
            end
            if ~isempty(testCase.SaveFolder) && ...
                    isfolder(testCase.SaveFolder)
                rmdir(testCase.SaveFolder, 's');
            end
        end
    end

    methods (Test)
        function testStartupAndHierarchyRendering(testCase)
            testCase.verifyEmpty(testCase.App.UIFigure.UserData);

            roots = testCase.App.ProjectTree.Children;
            testCase.verifyNumElements(roots, 1);
            testCase.verifyEqual(roots(1).NodeData.kind, 'project');

            subjectNode = testCase.findNode( ...
                roots(1), 'subject', testCase.SubjectUUID);
            testCase.verifyNotEmpty(subjectNode);
            sessionNode = testCase.findNode( ...
                subjectNode, 'session', testCase.SessionUUID);
            testCase.verifyNotEmpty(sessionNode);

            testCase.selectNode(subjectNode);
            testCase.verifyTrue(startsWith( ...
                testCase.App.ItemHeaderLabel.Text, 'Subject:'));
            testCase.verifyEqual(char( ...
                testCase.App.AddSessionButton.Visible), 'on');
            testCase.verifyGreaterThanOrEqual( ...
                numel(testCase.App.MetadataTabGroup.Children), 2);

            testCase.selectNode(sessionNode);
            testCase.verifyTrue(startsWith( ...
                testCase.App.ItemHeaderLabel.Text, 'Session:'));
            testCase.verifyEqual(char( ...
                testCase.App.RebindSaveFolderButton.Visible), 'on');
            testCase.verifyEqual(char( ...
                testCase.App.RepairSaveFolderButton.Visible), 'off');
            testCase.verifyEqual(char( ...
                testCase.App.OpenSaveFolderButton.Visible), 'on');
            testCase.verifyEqual(char( ...
                testCase.App.BindRawFolderButton.Visible), 'on');
            testCase.verifyEqual( ...
                testCase.App.BindRawFolderButton.Text, ...
                'Bind Raw Folder');
            testCase.verifyEqual(char( ...
                testCase.App.RemoveSessionFromProjectButton.Visible), ...
                'on');
            testCase.verifyFalse(isprop( ...
                testCase.App, 'RetrySaveFolderButton'));
            testCase.verifyFalse(isprop( ...
                testCase.App, 'AddRegistrationTransformButton'));
        end

        function testFourColumnLayoutAndScrollableMetadata(testCase)
            testCase.verifyEqual(numel( ...
                testCase.App.SelectedItemBottomGrid.ColumnWidth), 4);

            root = testCase.App.ProjectTree.Children(1);
            subjectNode = testCase.findNode( ...
                root, 'subject', testCase.SubjectUUID);
            sessionNode = testCase.findNode( ...
                subjectNode, 'session', testCase.SessionUUID);
            testCase.selectNode(sessionNode);

            visibleButtons = testCase.visibleActionButtons();
            columns = arrayfun(@(button) ...
                button.Layout.Column, visibleButtons);
            testCase.verifyLessThanOrEqual(max(columns), 4);
            testCase.verifyTrue(any(columns == 4));

            tabs = testCase.App.MetadataTabGroup.Children;
            testCase.verifyNotEmpty(tabs);
            for iTab = 1:numel(tabs)
                testCase.verifyEqual(char(tabs(iTab).Scrollable), 'on');
                grids = findobj(tabs(iTab), ...
                    'Type', 'uigridlayout');
                for iGrid = 1:numel(grids)
                    testCase.verifyEqual( ...
                        char(grids(iGrid).Scrollable), 'on');
                end
            end
        end

        function testRegistryColumnsAndResourcesStayOutOfTree(testCase)
            projectInfo = testCase.Store.getProjectInfo();
            sourceFile = fullfile( ...
                testCase.SaveFolder, 'reference_for_gui_test.mat');
            ImageReference = genImageReferenceStruct( ...
                single(ones(4, 5)), ...
                'Name', 'GUI Test Reference', ...
                'ProjectUUID', projectInfo.projectUUID, ...
                'ProjectName', projectInfo.projectName, ...
                'SubjectUUID', testCase.SubjectUUID, ...
                'SubjectID', 'Subject_01', ...
                'SourceType', 'unit-test', ...
                'CreatedBy', 'TestProjectManagerTool');
            save(sourceFile, 'ImageReference', '-mat');
            testCase.Store.addImageReference( ...
                'Subject_01', sourceFile, struct());

            testCase.App.RefreshButton.ButtonPushedFcn( ...
                testCase.App.RefreshButton, []);
            drawnow

            root = testCase.App.ProjectTree.Children(1);
            resourceNode = testCase.findNode(root, 'resource', '');
            testCase.verifyEmpty(resourceNode);

            testCase.selectNode(root);
            registryTab = testCase.findTab('Subject Registry');
            testCase.assertNotEmpty(registryTab);
            registryTable = testCase.findTable(registryTab);
            testCase.assertNotEmpty(registryTable);
            testCase.verifyClass(registryTable.Data, 'table');
            testCase.verifyEqual( ...
                registryTable.Data.Properties.VariableNames(1:3), ...
                {'id', 'displayName', 'uuid'});
            testCase.verifyFalse(ismember( ...
                'Value', registryTable.Data.Properties.VariableNames));

            subjectNode = testCase.findNode( ...
                root, 'subject', testCase.SubjectUUID);
            testCase.selectNode(subjectNode);
            resourceTab = testCase.findTab('Resource Registry');
            testCase.assertNotEmpty(resourceTab);
            resourceTable = testCase.findTable(resourceTab);
            testCase.assertNotEmpty(resourceTable);
            testCase.verifyEqual(height(resourceTable.Data), 1);
        end

        function testValidationHelpWithoutButtonIcon(testCase)
            testCase.verifySubstring( ...
                testCase.App.ValidationModeDropDown.Tooltip, ...
                'without recomputing file checksums');
            testCase.verifySubstring( ...
                testCase.App.ValidationModeDropDown.Tooltip, ...
                'SHA-256');

            testCase.Store.validate('Mode', 'quick');
            testCase.App.RefreshButton.ButtonPushedFcn( ...
                testCase.App.RefreshButton, []);
            drawnow
            testCase.verifyEqual( ...
                string(testCase.App.ValidateButton.Icon), "");
        end

        function testRawFolderPlaceholderAndContextualButton(testCase)
            root = testCase.App.ProjectTree.Children(1);
            subjectNode = testCase.findNode( ...
                root, 'subject', testCase.SubjectUUID);
            sessionNode = testCase.findNode( ...
                subjectNode, 'session', testCase.SessionUUID);
            testCase.selectNode(sessionNode);

            missingControl = testCase.findControlByValue( ...
                testCase.App.MetadataTabGroup, 'MISSING');
            testCase.verifyNotEmpty(missingControl);
            testCase.verifyEqual( ...
                testCase.App.BindRawFolderButton.Text, ...
                'Bind Raw Folder');

            rawFolder = tempname;
            mkdir(rawFolder);
            testCase.writeRawDataFile(rawFolder);
            cleanupRawFolder = onCleanup(@() ...
                testCase.removeFolderIfPresent(rawFolder));
            testCase.Store.bindRawDataFolder( ...
                'Subject_01', 'Session_01', rawFolder);
            testCase.App.RefreshButton.ButtonPushedFcn( ...
                testCase.App.RefreshButton, []);
            drawnow

            root = testCase.App.ProjectTree.Children(1);
            subjectNode = testCase.findNode( ...
                root, 'subject', testCase.SubjectUUID);
            sessionNode = testCase.findNode( ...
                subjectNode, 'session', testCase.SessionUUID);
            testCase.selectNode(sessionNode);
            pathControl = testCase.findControlByValue( ...
                testCase.App.MetadataTabGroup, rawFolder);
            testCase.verifyNotEmpty(pathControl);
            testCase.verifyEqual( ...
                testCase.App.BindRawFolderButton.Text, ...
                'Change Raw Folder');
            clear cleanupRawFolder
        end

        function testValidAcquisitionDatePrefillDropsTime(testCase)
            AcqInfoStream = struct( ...
                'DateTime', '20200630_113506');
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');

            actual = getSessionAcquisitionDatePrefill( ...
                testCase.SaveFolder);

            testCase.verifyEqual(actual, datetime(2020, 6, 30));
            testCase.verifyEqual(timeofday(actual), hours(0));
        end

        function testMissingOrMalformedDatePrefillReturnsNaT(testCase)
            actual = getSessionAcquisitionDatePrefill( ...
                testCase.SaveFolder);
            testCase.verifyTrue(isnat(actual));

            AcqInfoStream = struct('DateTime', 'not-a-date');
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), ...
                'AcqInfoStream', '-mat');
            actual = getSessionAcquisitionDatePrefill( ...
                testCase.SaveFolder);
            testCase.verifyTrue(isnat(actual));
        end

        function testSingularAcquisitionMetadataNameIsIgnored(testCase)
            AcqInfoStream = struct( ...
                'DateTime', '20200630_113506');
            save(fullfile(testCase.SaveFolder, 'AcqInfo.mat'), ...
                'AcqInfoStream', '-mat');

            actual = getSessionAcquisitionDatePrefill( ...
                testCase.SaveFolder);
            testCase.verifyTrue(isnat(actual));
        end

        function testSessionRigChoicesFilterUnavailableRigs(testCase)
            rigs = table( ...
                ["uuid-a"; "uuid-b"; "uuid-c"], ...
                ["Rig_A"; "Rig_B"; "Rig_C"], ...
                ["Primary Rig"; "Archived Rig"; "Unreadable Rig"], ...
                [true; false; false], ...
                ["active"; "archived"; "active"], ...
                [true; true; false], ...
                ["A"; "B"; "C"], ...
                'VariableNames', { ...
                'RigUUID', 'RigID', 'DisplayName', ...
                'IsDefault', 'Status', 'IsReadable', 'RigRoot'});

            [labels, ids, defaultRigID] = ...
                buildSessionRigChoices(rigs);

            testCase.verifyEqual(labels, ...
                {'(No rig)', 'Primary Rig (Rig_A)'});
            testCase.verifyEqual(ids, {'', 'Rig_A'});
            testCase.verifyEqual(defaultRigID, 'Rig_A');
        end

        function testNoAvailableRigsLeavesNoRigChoice(testCase)
            [labels, ids, defaultRigID] = ...
                buildSessionRigChoices(table());

            testCase.verifyEqual(labels, {'(No rig)'});
            testCase.verifyEqual(ids, {''});
            testCase.verifyEqual(defaultRigID, '');
        end

        function testDirtyStateAndRevert(testCase)
            root = testCase.App.ProjectTree.Children(1);
            subjectNode = testCase.findNode( ...
                root, 'subject', testCase.SubjectUUID);
            testCase.selectNode(subjectNode);

            displayNameControl = testCase.findControlByValue( ...
                testCase.App.MetadataTabGroup, 'Subject One');
            testCase.assertNotEmpty(displayNameControl);
            originalValue = displayNameControl.Value;
            displayNameControl.Value = 'Unsaved Name';
            displayNameControl.ValueChangedFcn( ...
                displayNameControl, []);

            testCase.verifyEqual(char( ...
                testCase.App.SaveChangesButton.Enable), 'on');

            testCase.App.RevertButton.ButtonPushedFcn( ...
                testCase.App.RevertButton, []);

            revertedControl = testCase.findControlByValue( ...
                testCase.App.MetadataTabGroup, originalValue);
            testCase.assertNotEmpty(revertedControl);
            testCase.verifyEqual( ...
                revertedControl.Value, ...
                originalValue);
        end

        function testLaunchByProjectUUID(testCase)
            delete(testCase.App);
            testCase.App = ProjectManagerTool( ...
                testCase.Store.getProjectInfo().projectUUID);
            drawnow

            testCase.verifyEmpty(testCase.App.UIFigure.UserData);
            testCase.verifyNumElements( ...
                testCase.App.ProjectTree.Children, 1);
        end
    end

    methods (Access = private)
        function selectNode(testCase, node)
            testCase.App.ProjectTree.SelectedNodes = node;
            callback = ...
                testCase.App.ProjectTree.SelectionChangedFcn;
            callback(testCase.App.ProjectTree, ...
                struct('SelectedNodes', node));
            drawnow
        end

        function node = findNode(testCase, parent, kind, uuid)
            node = [];
            children = parent.Children;
            for iChild = 1:numel(children)
                child = children(iChild);
                matchesKind = isstruct(child.NodeData) && ...
                    isfield(child.NodeData, 'kind') && ...
                    strcmp(child.NodeData.kind, kind);
                matchesUUID = isempty(uuid) || ...
                    (isfield(child.NodeData, 'uuid') && ...
                    strcmp(child.NodeData.uuid, uuid));
                if matchesKind && matchesUUID
                    node = child;
                    return
                end
                node = testCase.findNode(child, kind, uuid);
                if ~isempty(node)
                    return
                end
            end
        end

        function control = findControlByValue( ...
                testCase, parent, expectedValue)
            control = [];
            if ~isprop(parent, 'Children')
                return
            end

            children = parent.Children;
            for iChild = 1:numel(children)
                child = children(iChild);
                if isprop(child, 'Value')
                    value = child.Value;
                    if (ischar(value) || ...
                            (isstring(value) && isscalar(value))) && ...
                            strcmp(char(string(value)), ...
                            char(string(expectedValue)))
                        control = child;
                        return
                    end
                end

                control = testCase.findControlByValue( ...
                    child, expectedValue);
                if ~isempty(control)
                    return
                end
            end
        end

        function buttons = visibleActionButtons(testCase)
            children = testCase.App.SelectedItemBottomGrid.Children;
            isButton = arrayfun(@(child) ...
                isa(child, 'matlab.ui.control.Button'), children);
            buttons = children(isButton);
            isVisible = arrayfun(@(button) ...
                strcmp(char(button.Visible), 'on'), buttons);
            buttons = buttons(isVisible);
        end

        function tab = findTab(testCase, titleText)
            tab = [];
            tabs = testCase.App.MetadataTabGroup.Children;
            for iTab = 1:numel(tabs)
                if strcmp(tabs(iTab).Title, titleText)
                    tab = tabs(iTab);
                    return
                end
            end
        end

        function tableControl = findTable(testCase, parent)
            tableControl = [];
            if isa(parent, 'matlab.ui.control.Table')
                tableControl = parent;
                return
            end
            if ~isprop(parent, 'Children')
                return
            end
            children = parent.Children;
            for iChild = 1:numel(children)
                tableControl = testCase.findTable(children(iChild));
                if ~isempty(tableControl)
                    return
                end
            end
        end

        function removeFolderIfPresent(~, folderPath)
            if isfolder(folderPath)
                rmdir(folderPath, 's');
            end
        end

        function writeRawDataFile(~, rawFolder)
            filePath = fullfile(rawFolder, 'raw_data.bin');
            [fid, message] = fopen(filePath, 'w');
            if fid < 0
                error('TestProjectManagerTool:fixtureWriteFailed', ...
                    'Could not create raw-data fixture: %s', message);
            end
            cleanupFile = onCleanup(@() fclose(fid));
            fwrite(fid, uint8(0), 'uint8');
            clear cleanupFile
        end
    end
end
