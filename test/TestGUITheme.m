classdef TestGUITheme < matlab.unittest.TestCase
%TESTGUITHEME Tests for umIT GUI preferences and R2024b theme handling.

    properties
        RepoRoot
        TemporaryFolders cell = {}
        Figures cell = {}
    end

    methods (TestClassSetup)
        function addRepositoryToPath(testCase)
            testCase.RepoRoot = fileparts(fileparts(mfilename('fullpath')));
            requiredFolders = {
                testCase.RepoRoot
                fullfile(testCase.RepoRoot, 'GUI')
                fullfile(testCase.RepoRoot, 'subFunc')
                };
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                requiredFolders));
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryArtifacts(testCase)
            for idx = 1:numel(testCase.Figures)
                figureHandle = testCase.Figures{idx};
                if isvalid(figureHandle)
                    delete(figureHandle);
                end
            end

            for idx = 1:numel(testCase.TemporaryFolders)
                folder = testCase.TemporaryFolders{idx};
                if isfolder(folder)
                    rmdir(folder, 's');
                end
            end

            testCase.Figures = {};
            testCase.TemporaryFolders = {};
        end
    end

    methods (Test)
        function testDefaultPreferenceFileIsCreated(testCase)
            preferencesFolder = testCase.newTemporaryFolder(false);

            mode = getUmitTheme( ...
                'PreferencesFolder', preferencesFolder);

            testCase.verifyEqual(mode, 'light');
            preferencesFile = fullfile(preferencesFolder, ...
                'appPreferences.json');
            testCase.verifyTrue(isfile(preferencesFile));

            preferences = jsondecode(fileread(preferencesFile));
            testCase.verifyEqual(preferences.schemaVersion, 1);
            testCase.verifyEqual(preferences.theme, 'light');
        end

        function testThemeUpdatePreservesOtherPreferences(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            getUmitTheme('PreferencesFolder', preferencesFolder);

            preferencesFile = fullfile(preferencesFolder, ...
                'appPreferences.json');
            preferences = jsondecode(fileread(preferencesFile));
            preferences.previewFrameCount = 25;
            preferences.viewer = struct('showEvents', true);
            testCase.writeJson(preferencesFile, preferences);

            returnedFile = setUmitTheme('dark', ...
                'PreferencesFolder', preferencesFolder);
            updated = jsondecode(fileread(returnedFile));

            testCase.verifyEqual(updated.theme, 'dark');
            testCase.verifyEqual(updated.previewFrameCount, 25);
            testCase.verifyTrue(updated.viewer.showEvents);
            testCase.verifyError(@() setUmitTheme('system', ...
                'PreferencesFolder', preferencesFolder), ...
                'Umitoolbox:setUmitTheme:InvalidTheme');

            folderFiles = dir(preferencesFolder);
            folderFiles = folderFiles(~[folderFiles.isdir]);
            testCase.verifyEqual({folderFiles.name}, ...
                {'appPreferences.json'});
        end

        function testInvalidPreferencesFallBackToLight(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            preferencesFile = fullfile(preferencesFolder, ...
                'appPreferences.json');

            testCase.writeText(preferencesFile, '{not valid JSON');
            testCase.verifyWarning(@() getUmitTheme( ...
                'PreferencesFolder', preferencesFolder), ...
                'Umitoolbox:Preferences:InvalidFile');

            testCase.writeText(preferencesFile, ...
                '{"schemaVersion":1,"theme":"system"}');
            testCase.verifyWarning(@() getUmitTheme( ...
                'PreferencesFolder', preferencesFolder), ...
                'Umitoolbox:Preferences:InvalidTheme');

            testCase.verifyEqual(testCase.callGetThemeWithoutWarning( ...
                preferencesFolder), 'light');
        end

        function testColorPaletteCoversCommonUIComponents(testCase)
            figureHandle = testCase.newFigure();
            button = uibutton(figureHandle);
            panel = uipanel(figureHandle);
            label = uilabel(panel);
            axesHandle = uiaxes(figureHandle);
            tableHandle = uitable(figureHandle, ...
                'Data', [1 2; 3 4]);
            tableBackground = [0.20 0.30 0.40; 0.30 0.40 0.50];
            tableForeground = [0.60 0.70 0.80];
            tableHandle.BackgroundColor = tableBackground;
            tableHandle.ForegroundColor = tableForeground;

            setGUIcolorScheme(figureHandle, 'dark');

            testCase.verifyEqual(figureHandle.Color, ...
                [0.27 0.27 0.27], 'AbsTol', 1e-12);
            testCase.verifyEqual(button.BackgroundColor, ...
                [0.50 0.50 0.50], 'AbsTol', 1e-12);
            testCase.verifyEqual(button.FontColor, ...
                [0.90 0.90 0.90], 'AbsTol', 1e-12);
            testCase.verifyEqual(panel.BackgroundColor, ...
                [0.27 0.27 0.27], 'AbsTol', 1e-12);
            testCase.verifyEqual(panel.ForegroundColor, ...
                [0.90 0.90 0.90], 'AbsTol', 1e-12);
            testCase.verifyEqual(label.BackgroundColor, 'none');
            testCase.verifyEqual(label.FontColor, ...
                [0.90 0.90 0.90], 'AbsTol', 1e-12);
            testCase.verifyEqual(axesHandle.Color, ...
                [0.43 0.43 0.43], 'AbsTol', 1e-12);
            testCase.verifyEqual(axesHandle.XColor, ...
                [0.90 0.90 0.90], 'AbsTol', 1e-12);
            testCase.verifyEqual(tableHandle.BackgroundColor, ...
                tableBackground, 'AbsTol', 1e-12);
            testCase.verifyEqual(tableHandle.ForegroundColor, ...
                tableForeground, 'AbsTol', 1e-12);

            setGUIcolorScheme(figureHandle, 'light');

            testCase.verifyEqual(figureHandle.Color, ...
                [0.94 0.94 0.94], 'AbsTol', 1e-12);
            testCase.verifyEqual(button.BackgroundColor, ...
                [0.90 0.90 0.90], 'AbsTol', 1e-12);
            testCase.verifyEqual(button.FontColor, ...
                [0 0 0], 'AbsTol', 1e-12);
            testCase.verifyEqual(axesHandle.Color, ...
                [1 1 1], 'AbsTol', 1e-12);
            testCase.verifyEqual(getappdata(figureHandle, 'UmitTheme'), ...
                'light');
        end

        function testTaggedIconsSwitchOnlyWhenPairExists(testCase)
            figureHandle = testCase.newFigure();
            button = uibutton(figureHandle);
            imageHandle = uiimage(figureHandle);
            iconsFolder = fullfile(testCase.RepoRoot, 'GUI', 'icons');

            blackIcon = fullfile(iconsFolder, 'icon_plusBlack.png');
            whiteIcon = fullfile(iconsFolder, 'icon_plusWhite.png');
            untaggedIcon = fullfile(iconsFolder, 'icon_check.png');
            unpairedIcon = fullfile(iconsFolder, ...
                'icon_stopButton_Black.png');

            button.Icon = blackIcon;
            imageHandle.ImageSource = blackIcon;
            setGUIcolorScheme(figureHandle, 'dark');
            testCase.verifyEqual(button.Icon, whiteIcon);
            testCase.verifyEqual(imageHandle.ImageSource, whiteIcon);

            setGUIcolorScheme(figureHandle, 'light');
            testCase.verifyEqual(button.Icon, blackIcon);
            testCase.verifyEqual(imageHandle.ImageSource, blackIcon);

            button.Icon = untaggedIcon;
            setGUIcolorScheme(figureHandle, 'dark');
            testCase.verifyEqual(button.Icon, untaggedIcon);

            button.Icon = unpairedIcon;
            setGUIcolorScheme(figureHandle, 'dark');
            testCase.verifyEqual(button.Icon, unpairedIcon);
        end

        function testMenusExcludedAndTabFontsPreserved(testCase)
            figureHandle = testCase.newFigure();
            menuHandle = uimenu(figureHandle, 'Text', 'Menu');
            tabGroup = uitabgroup(figureHandle);
            tabHandle = uitab(tabGroup, 'Title', 'Tab');
            secondTabHandle = uitab(tabGroup, 'Title', 'Second tab');
            button = uibutton(tabHandle);
            secondButton = uibutton(secondTabHandle);

            menuColor = [0.10 0.20 0.30];
            tabBackground = [0.15 0.25 0.35];
            tabForeground = [0.45 0.55 0.65];
            secondTabBackground = [0.20 0.30 0.40];
            secondTabForeground = [0.70 0.30 0.20];
            menuHandle.ForegroundColor = menuColor;
            tabHandle.BackgroundColor = tabBackground;
            tabHandle.ForegroundColor = tabForeground;
            secondTabHandle.BackgroundColor = secondTabBackground;
            secondTabHandle.ForegroundColor = secondTabForeground;

            setGUIcolorScheme(figureHandle, 'dark');

            testCase.verifyClass(tabGroup, ...
                'matlab.ui.container.TabGroup');
            testCase.verifyEqual(menuHandle.ForegroundColor, ...
                menuColor, 'AbsTol', 1e-12);
            testCase.verifyEqual(tabHandle.BackgroundColor, ...
                [0.27 0.27 0.27], 'AbsTol', 1e-12);
            testCase.verifyEqual(tabHandle.ForegroundColor, ...
                tabForeground, 'AbsTol', 1e-12);
            testCase.verifyEqual(secondTabHandle.BackgroundColor, ...
                [0.27 0.27 0.27], 'AbsTol', 1e-12);
            testCase.verifyEqual(secondTabHandle.ForegroundColor, ...
                secondTabForeground, 'AbsTol', 1e-12);
            testCase.verifyEqual(button.BackgroundColor, ...
                [0.50 0.50 0.50], 'AbsTol', 1e-12);
            testCase.verifyEqual(secondButton.BackgroundColor, ...
                [0.50 0.50 0.50], 'AbsTol', 1e-12);

            setGUIcolorScheme(figureHandle, 'light');
            testCase.verifyEqual(tabHandle.BackgroundColor, ...
                [0.94 0.94 0.94], 'AbsTol', 1e-12);
            testCase.verifyEqual(tabHandle.ForegroundColor, ...
                tabForeground, 'AbsTol', 1e-12);
            testCase.verifyEqual(secondTabHandle.BackgroundColor, ...
                [0.94 0.94 0.94], 'AbsTol', 1e-12);
            testCase.verifyEqual(secondTabHandle.ForegroundColor, ...
                secondTabForeground, 'AbsTol', 1e-12);
        end

        function testLightStartupLeavesDesignerColorsUnchanged(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            setUmitTheme('light', ...
                'PreferencesFolder', preferencesFolder);

            figureHandle = testCase.newFigure();
            button = uibutton(figureHandle);
            figureHandle.Color = [0.20 0.30 0.40];
            button.BackgroundColor = [0.60 0.70 0.80];

            mode = initializeUmitTheme(figureHandle, ...
                'PreferencesFolder', preferencesFolder);

            testCase.verifyEqual(mode, 'light');
            testCase.verifyEqual(figureHandle.Color, ...
                [0.20 0.30 0.40], 'AbsTol', 1e-12);
            testCase.verifyEqual(button.BackgroundColor, ...
                [0.60 0.70 0.80], 'AbsTol', 1e-12);
            testCase.verifyEqual(getappdata(figureHandle, 'UmitTheme'), ...
                'light');
        end

        function testInitializationAndThemeMenu(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            setUmitTheme('dark', ...
                'PreferencesFolder', preferencesFolder);
            figureHandle = testCase.newFigure();

            mode = initializeUmitTheme(figureHandle, ...
                'ShowMenu', true, ...
                'ThemeChangedFcn', @(theme) setappdata( ...
                    figureHandle, 'TestCallbackTheme', theme), ...
                'PreferencesFolder', preferencesFolder);
            initializeUmitTheme(figureHandle, ...
                'ShowMenu', true, ...
                'PreferencesFolder', preferencesFolder);

            testCase.verifyEqual(mode, 'dark');
            testCase.verifyEqual(getappdata(figureHandle, ...
                'TestCallbackTheme'), 'dark');
            testCase.verifyEqual(figureHandle.Color, ...
                [0.27 0.27 0.27], 'AbsTol', 1e-12);

            themeMenus = findall(figureHandle, 'Type', 'uimenu', ...
                'Tag', 'UmitThemeMenu');
            lightMenu = findall(figureHandle, 'Type', 'uimenu', ...
                'Tag', 'UmitLightThemeMenu');
            darkMenu = findall(figureHandle, 'Type', 'uimenu', ...
                'Tag', 'UmitDarkThemeMenu');

            testCase.verifyNumElements(themeMenus, 1);
            testCase.verifyNumElements(lightMenu, 1);
            testCase.verifyNumElements(darkMenu, 1);
            testCase.verifyEqual(char(lightMenu.Checked), 'off');
            testCase.verifyEqual(char(darkMenu.Checked), 'on');

            callback = lightMenu.MenuSelectedFcn;
            callback(lightMenu, []);

            testCase.verifyEqual(getUmitTheme( ...
                'PreferencesFolder', preferencesFolder), 'light');
            testCase.verifyEqual(getappdata(figureHandle, ...
                'TestCallbackTheme'), 'light');
            testCase.verifyEqual(figureHandle.Color, ...
                [0.94 0.94 0.94], 'AbsTol', 1e-12);
            testCase.verifyEqual(char(lightMenu.Checked), 'on');
            testCase.verifyEqual(char(darkMenu.Checked), 'off');
        end

        function testInitializationWithoutMenuAppliesDark(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            setUmitTheme('dark', ...
                'PreferencesFolder', preferencesFolder);
            figureHandle = testCase.newFigure();

            initializeUmitTheme(figureHandle, ...
                'PreferencesFolder', preferencesFolder);

            testCase.verifyEqual(figureHandle.Color, ...
                [0.27 0.27 0.27], 'AbsTol', 1e-12);
            testCase.verifyEmpty(findall(figureHandle, ...
                'Type', 'uimenu', 'Tag', 'UmitThemeMenu'));
        end

        function testStatusMessageUsesInitializedFigureTheme(testCase)
            lightFigure = testCase.newFigure();
            lightLabel = uilabel(lightFigure);
            setStatusMessage(lightLabel, 'Problem', ...
                'Severity', 'error', 'Refresh', false);
            testCase.verifyEqual(lightLabel.FontColor, ...
                [0.80 0 0], 'AbsTol', 1e-12);

            preferencesFolder = testCase.newTemporaryFolder(true);
            setUmitTheme('dark', ...
                'PreferencesFolder', preferencesFolder);
            darkFigure = testCase.newFigure();
            darkLabel = uilabel(darkFigure);
            initializeUmitTheme(darkFigure, ...
                'PreferencesFolder', preferencesFolder);

            setStatusMessage(darkLabel, 'Problem', ...
                'Severity', 'error', 'Refresh', false);
            testCase.verifyEqual(darkLabel.FontColor, ...
                [1.00 0.35 0.35], 'AbsTol', 1e-12);
        end
    end

    methods (Access = private)
        function folder = newTemporaryFolder(testCase, createNow)
            folder = tempname;
            testCase.TemporaryFolders{end+1} = folder;
            if createNow
                mkdir(folder);
            end
        end

        function figureHandle = newFigure(testCase)
            figureHandle = uifigure( ...
                'Visible', 'off', ...
                'Tag', 'UmitThemeTestFigure');
            testCase.Figures{end+1} = figureHandle;
        end

        function writeJson(testCase, filePath, value)
            testCase.writeText(filePath, ...
                jsonencode(value, 'PrettyPrint', true));
        end

        function writeText(~, filePath, value)
            fileId = fopen(filePath, 'w', 'n', 'UTF-8');
            cleanupFile = onCleanup(@() fclose(fileId));
            fwrite(fileId, value, 'char');
            clear cleanupFile
        end

        function mode = callGetThemeWithoutWarning(~, preferencesFolder)
            warningState = warning('off', ...
                'Umitoolbox:Preferences:InvalidTheme');
            cleanupWarning = onCleanup(@() warning(warningState));
            mode = getUmitTheme('PreferencesFolder', preferencesFolder);
            clear cleanupWarning
        end
    end
end
