classdef TestDataViewerManagementMenu < matlab.unittest.TestCase
    %TESTDATAVIEWERMANAGEMENTMENU Verify the DataViewer management hierarchy.

    methods (Test)
        function testManagementMenuContainsSessionActions(testCase)
            app = DataViewer();
            testCase.addTeardown(@() iDeleteApp(app));

            testCase.verifyEqual(app.ManagementMenu.Text, 'Management');
            testCase.verifyEqual(app.SetRawFolderMenu.Parent, ...
                app.ManagementMenu);
            testCase.verifyFalse(any(strcmp( ...
                {app.FileMenu.Children.Text}, 'Set Raw Folder...')));

            testCase.verifyEqual(app.ProjectMenu.Parent, app.ManagementMenu);
            testCase.verifyEqual(char(app.ProjectMenu.Separator), 'on');
            testCase.verifyEqual(app.BindCurrentSessionMenu.Parent, ...
                app.ProjectMenu);
            testCase.verifyEqual(app.BindCurrentSessionMenu.Text, ...
                'Bind Current Session...');
            testCase.verifyEqual(app.ProjectManagerMenu.Parent, ...
                app.ProjectMenu);

            testCase.verifyEqual(app.RigMenu.Parent, app.ManagementMenu);
            testCase.verifyEqual(app.AssignRigMenu.Parent, app.RigMenu);
            testCase.verifyEqual(app.AssignRigMenu.Text, 'Assign Rig');
            testCase.verifyEqual(app.RigManagerMenu.Parent, app.RigMenu);

            testCase.verifyEqual(app.PreferencesMenu.Parent, app.FileMenu);
            testCase.verifyEqual(app.PreferencesMenu.Text, 'Preferences...');
            testCase.verifyEqual(char(app.PreferencesMenu.Enable), 'on');
            testCase.verifyNotEmpty(app.PreferencesMenu.MenuSelectedFcn);
            testCase.verifyEmpty(findall(app.UIFigure, ...
                'Type', 'uimenu', 'Tag', 'UmitThemeMenu'));
            testCase.verifyEqual(getappdata(app.UIFigure, 'UmitTheme'), ...
                getUmitTheme());
        end

        function testPreferencesDialogOpens(testCase)
            app = DataViewer();
            testCase.addTeardown(@() iDeleteApp(app));

            dialogFigure = openDataViewerPreferencesDialog( ...
                app.UIFigure, ...
                app.ColormapDropDown.Items, ...
                app.ColormapDropDown.ItemsData, ...
                []);
            testCase.addTeardown(@() iDeleteApp(dialogFigure));

            testCase.verifyEqual(dialogFigure.Name, ...
                'DataViewer Preferences');
            testCase.verifyEqual(char(dialogFigure.WindowStyle), 'modal');
            testCase.verifyNumElements(findall(dialogFigure, ...
                'Tag', 'ThemeDropDown'), 1);
            testCase.verifyNumElements(findall(dialogFigure, ...
                'Tag', 'DefaultColormapDropDown'), 1);
            testCase.verifyNumElements(findall(dialogFigure, ...
                'Tag', 'ReopenLastFileCheckBox'), 1);
            testCase.verifyNumElements(findall(dialogFigure, ...
                'Tag', 'RememberLastSaveFolderCheckBox'), 1);
            defaultFolderField = findall(dialogFigure, ...
                'Tag', 'DefaultFolderField');
            testCase.verifyNumElements(defaultFolderField, 1);
            testCase.verifyEqual(char(defaultFolderField.Editable), 'off');
            testCase.verifyNumElements(findall(dialogFigure, ...
                'Tag', 'BrowseDefaultFolderButton'), 1);
            testCase.verifyNumElements(findall(dialogFigure, ...
                'Tag', 'SaveDataViewerPreferencesButton'), 1);
        end
    end
end

function iDeleteApp(app)
if ~isempty(app) && isvalid(app)
    delete(app);
end
end
