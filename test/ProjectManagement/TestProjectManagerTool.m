classdef TestProjectManagerTool < matlab.unittest.TestCase
    %TESTPROJECTMANAGERTOOL Integration tests for ProjectManagerTool.mlapp.

    properties
        App
        Store
        ProjectRoot = ''
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
            testCase.SessionUUID = testCase.Store.addSession( ...
                'Subject_01', struct( ...
                'sessionID', 'Session_01', ...
                'displayName', 'Session One'));

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
                testCase.App.BindRawFolderButton.Visible), 'on');
            testCase.verifyEqual(char( ...
                testCase.App.AddRegistrationTransformButton.Visible), ...
                'on');
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
                if strcmp(child.NodeData.kind, kind) && ...
                        strcmp(child.NodeData.uuid, uuid)
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
    end
end
