classdef TestDataViewerPreferences < matlab.unittest.TestCase
%TESTDATAVIEWERPREFERENCES Focused persistence and startup-path tests.

    properties
        RepoRoot
        TemporaryFolders cell = {}
    end

    methods (TestClassSetup)
        function addRepositoryToPath(testCase)
            testCase.RepoRoot = fileparts(fileparts( ...
                fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture({ ...
                testCase.RepoRoot
                fullfile(testCase.RepoRoot, 'GUI')
                fullfile(testCase.RepoRoot, 'subFunc')}));
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolders(testCase)
            for idx = 1:numel(testCase.TemporaryFolders)
                folder = testCase.TemporaryFolders{idx};
                if isfolder(folder)
                    rmdir(folder, 's');
                end
            end
            testCase.TemporaryFolders = {};
        end
    end

    methods (Test)
        function testFirstRunReturnsDefaults(testCase)
            preferencesFolder = testCase.newTemporaryFolder(false);

            preferences = DataViewerPreferences.load( ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'gray'});

            testCase.verifyEqual(preferences, ...
                DataViewerPreferences.defaults());
            saved = jsondecode(fileread(fullfile( ...
                preferencesFolder, 'appPreferences.json')));
            testCase.verifyTrue(isfield(saved, 'dataViewer'));
        end

        function testInvalidAndPartialValuesFallBackSafely(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            raw = struct( ...
                'schemaVersion', 1, ...
                'theme', 'system', ...
                'dataViewer', struct( ...
                    'defaultColormap', 'missingMap', ...
                    'reopenLastFile', 'yes', ...
                    'rememberLastSaveFolder', 1, ...
                    'defaultFolder', 42, ...
                    'lastFile', 42));
            testCase.writePreferences(preferencesFolder, raw);

            preferences = DataViewerPreferences.load( ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'gray', 'hot'});

            testCase.verifyEqual(preferences.theme, 'light');
            testCase.verifyEqual(preferences.defaultColormap, 'gray');
            testCase.verifyFalse(preferences.reopenLastFile);
            testCase.verifyFalse(preferences.rememberLastSaveFolder);
            testCase.verifyEmpty(preferences.defaultFolder);
            testCase.verifyEmpty(preferences.lastFile);
            testCase.verifyEmpty(preferences.lastSaveFolder);
        end

        function testSavePreservesUnrelatedPreferences(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            raw = struct('schemaVersion', 1, 'theme', 'dark', ...
                'futureSetting', 17);
            testCase.writePreferences(preferencesFolder, raw);

            preferences = DataViewerPreferences.defaults();
            preferences.theme = 'dark';
            preferences.defaultColormap = 'hot';
            preferences.reopenLastFile = true;
            DataViewerPreferences.save(preferences, ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'hot'});

            saved = jsondecode(fileread(fullfile( ...
                preferencesFolder, 'appPreferences.json')));
            testCase.verifyEqual(saved.theme, 'dark');
            testCase.verifyEqual(saved.futureSetting, 17);
            testCase.verifyFalse(isfield(saved.dataViewer, 'theme'));
            testCase.verifyEqual(saved.dataViewer.defaultColormap, 'hot');
            testCase.verifyTrue(saved.dataViewer.reopenLastFile);
        end

        function testDialogFolderPrecedenceAndMissingPaths(testCase)
            root = testCase.newTemporaryFolder(true);
            remembered = fullfile(root, 'remembered');
            configured = fullfile(root, 'configured');
            fallback = fullfile(root, 'fallback');
            mkdir(remembered);
            mkdir(configured);
            mkdir(fallback);

            preferences = DataViewerPreferences.defaults();
            preferences.rememberLastSaveFolder = true;
            preferences.lastSaveFolder = remembered;
            preferences.defaultFolder = configured;

            testCase.verifyEqual( ...
                DataViewerPreferences.resolveDialogFolder( ...
                fallback, preferences), remembered);

            rmdir(remembered);
            testCase.verifyEqual( ...
                DataViewerPreferences.resolveDialogFolder( ...
                fallback, preferences), configured);

            rmdir(configured);
            testCase.verifyEqual( ...
                DataViewerPreferences.resolveDialogFolder( ...
                fallback, preferences), fallback);
        end

        function testStartupFileResolutionIsNonblocking(testCase)
            root = testCase.newTemporaryFolder(true);
            lastFile = fullfile(root, 'last.dat');
            testCase.writeText(lastFile, '');

            preferences = DataViewerPreferences.defaults();
            preferences.reopenLastFile = true;
            preferences.lastFile = lastFile;

            testCase.verifyEqual( ...
                DataViewerPreferences.resolveStartupFile('', preferences), ...
                lastFile);

            delete(lastFile);
            testCase.verifyEmpty( ...
                DataViewerPreferences.resolveStartupFile('', preferences));

            explicitFile = fullfile(root, 'explicit.umt');
            testCase.verifyEqual( ...
                DataViewerPreferences.resolveStartupFile( ...
                explicitFile, preferences), explicitFile);
        end

        function testSuccessfulOpenRecordsFileAndSaveFolder(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            saveFolder = testCase.newTemporaryFolder(true);
            dataFile = fullfile(saveFolder, 'recording.umt');
            testCase.writeText(dataFile, '');

            enabled = DataViewerPreferences.defaults();
            enabled.reopenLastFile = true;
            enabled.rememberLastSaveFolder = true;
            DataViewerPreferences.save(enabled, ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula'});

            DataViewerPreferences.recordSuccessfulOpen(dataFile, ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula'});
            preferences = DataViewerPreferences.load( ...
                'PreferencesFolder', preferencesFolder);

            testCase.verifyEqual(preferences.lastFile, dataFile);
            testCase.verifyEqual(preferences.lastSaveFolder, saveFolder);
            testCase.verifyError(@() ...
                DataViewerPreferences.recordSuccessfulOpen( ...
                fullfile(saveFolder, 'missing.dat'), ...
                'PreferencesFolder', preferencesFolder), ...
                'Umitoolbox:DataViewerPreferences:InvalidLastFile');
        end

        function testResetRestoresThemeAndPreservesUnknownRootSettings(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            raw = struct('schemaVersion', 1, 'theme', 'dark', ...
                'futureSetting', 17);
            testCase.writePreferences(preferencesFolder, raw);

            changed = DataViewerPreferences.defaults();
            changed.theme = 'dark';
            changed.defaultColormap = 'hot';
            changed.reopenLastFile = true;
            DataViewerPreferences.save(changed, ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'hot'});
            DataViewerPreferences.reset( ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'hot'});

            saved = jsondecode(fileread(fullfile( ...
                preferencesFolder, 'appPreferences.json')));
            defaults = DataViewerPreferences.defaults();
            testCase.verifyEqual(saved.theme, 'light');
            testCase.verifyEqual(saved.futureSetting, 17);
            testCase.verifyEqual(saved.dataViewer, ...
                rmfield(defaults, 'theme'));
        end

        function testDialogCancelDoesNotPersist(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            DataViewerPreferences.save( ...
                DataViewerPreferences.defaults(), ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'hot'});

            parentFigure = uifigure('Visible', 'off');
            testCase.addTeardown(@() iDeleteFigure(parentFigure));
            dialogFigure = openDataViewerPreferencesDialog( ...
                parentFigure, {'Parula', 'Hot'}, {'parula', 'hot'}, [], ...
                'PreferencesFolder', preferencesFolder);
            testCase.addTeardown(@() iDeleteFigure(dialogFigure));

            colormapDropDown = findall(dialogFigure, ...
                'Tag', 'DefaultColormapDropDown');
            themeDropDown = findall(dialogFigure, ...
                'Tag', 'ThemeDropDown');
            reopenCheckBox = findall(dialogFigure, ...
                'Tag', 'ReopenLastFileCheckBox');
            themeDropDown.Value = 'dark';
            colormapDropDown.Value = 'hot';
            reopenCheckBox.Value = true;
            cancelButton = findall(dialogFigure, ...
                'Tag', 'CancelDataViewerPreferencesButton');
            cancelFunction = cancelButton.ButtonPushedFcn;
            cancelFunction(cancelButton, []);

            saved = DataViewerPreferences.load( ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'hot'});
            testCase.verifyEqual(saved, ...
                DataViewerPreferences.defaults());
        end

        function testDialogOKPersistsAndApplies(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            selectedFolder = testCase.newTemporaryFolder(true);
            parentFigure = uifigure('Visible', 'off');
            testCase.addTeardown(@() iDeleteFigure(parentFigure));
            dialogFigure = openDataViewerPreferencesDialog( ...
                parentFigure, {'Parula', 'Hot'}, {'parula', 'hot'}, ...
                @(saved) setappdata(parentFigure, ...
                'AppliedPreferences', saved), ...
                'PreferencesFolder', preferencesFolder);
            testCase.addTeardown(@() iDeleteFigure(dialogFigure));

            colormapDropDown = findall(dialogFigure, ...
                'Tag', 'DefaultColormapDropDown');
            themeDropDown = findall(dialogFigure, ...
                'Tag', 'ThemeDropDown');
            reopenCheckBox = findall(dialogFigure, ...
                'Tag', 'ReopenLastFileCheckBox');
            rememberCheckBox = findall(dialogFigure, ...
                'Tag', 'RememberLastSaveFolderCheckBox');
            defaultFolderField = findall(dialogFigure, ...
                'Tag', 'DefaultFolderField');
            themeDropDown.Value = 'dark';
            colormapDropDown.Value = 'hot';
            reopenCheckBox.Value = true;
            rememberCheckBox.Value = true;
            defaultFolderField.Value = selectedFolder;
            saveButton = findall(dialogFigure, ...
                'Tag', 'SaveDataViewerPreferencesButton');
            saveFunction = saveButton.ButtonPushedFcn;
            saveFunction(saveButton, []);

            saved = DataViewerPreferences.load( ...
                'PreferencesFolder', preferencesFolder, ...
                'AvailableColormaps', {'parula', 'hot'});
            testCase.verifyEqual(saved.theme, 'dark');
            testCase.verifyEqual(getUmitTheme( ...
                'PreferencesFolder', preferencesFolder), 'dark');
            testCase.verifyEqual(saved.defaultColormap, 'hot');
            testCase.verifyTrue(saved.reopenLastFile);
            testCase.verifyTrue(saved.rememberLastSaveFolder);
            testCase.verifyEqual(saved.defaultFolder, selectedFolder);
            testCase.verifyEqual(getappdata(parentFigure, ...
                'AppliedPreferences'), saved);
        end

        function testPreferenceOperationsDoNotChangeWorkingFolder(testCase)
            preferencesFolder = testCase.newTemporaryFolder(true);
            fallbackFolder = testCase.newTemporaryFolder(true);
            originalFolder = pwd;

            preferences = DataViewerPreferences.defaults();
            preferences.defaultFolder = fullfile( ...
                fallbackFolder, 'missing');
            DataViewerPreferences.save(preferences, ...
                'PreferencesFolder', preferencesFolder);
            loaded = DataViewerPreferences.load( ...
                'PreferencesFolder', preferencesFolder);
            DataViewerPreferences.resolveDialogFolder( ...
                fallbackFolder, loaded);

            testCase.verifyEqual(pwd, originalFolder);
        end
    end

    methods (Access = private)
        function folder = newTemporaryFolder(testCase, createNow)
            folder = tempname;
            testCase.TemporaryFolders{end + 1} = folder;
            if createNow
                mkdir(folder);
            end
        end

        function writePreferences(testCase, folder, preferences)
            testCase.writeText(fullfile(folder, 'appPreferences.json'), ...
                jsonencode(preferences, 'PrettyPrint', true));
        end

        function writeText(~, filePath, contents)
            fileId = fopen(filePath, 'w', 'n', 'UTF-8');
            cleanupFile = onCleanup(@() fclose(fileId));
            fwrite(fileId, contents, 'char');
            clear cleanupFile
        end
    end
end

function iDeleteFigure(figureHandle)
if ~isempty(figureHandle) && isvalid(figureHandle)
    delete(figureHandle);
end
end
