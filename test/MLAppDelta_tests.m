classdef MLAppDelta_tests < matlab.unittest.TestCase
%MLAPPDELTA_TESTS Unit tests for the .mlapp delta-tracking tool chain
% (_mlapp_tools/exportMlappToM, parseMatlabExport, diffMatlabBlocks,
% generateDeltaFile, confirmAndUpdateReference, exportMatlabAppToDelta).
%
% Reference app: GUI/DataViewer/DataViewer_EventsManager.mlapp. Tests
% never modify this file; where a "changed app" is needed, either the
% parsed struct is mutated in memory, or a throwaway temp copy of the
% .mlapp is mutated at the zip level.
%
% Run with: runtests('test/MLAppDelta_tests.m')

    properties (Constant)
        SourceApp = fullfile('GUI', 'DataViewer', 'DataViewer_EventsManager.mlapp');
        AppNamePrefix = 'MLAppDeltaTest_';
    end

    properties
        RepoRoot
        AppsExportedDir
        TempFiles
        TestAppNames
    end

    methods (TestClassSetup)
        function addToolsToPath(testCase)
            here = fileparts(mfilename('fullpath'));
            testCase.RepoRoot = fileparts(here);
            testCase.AppsExportedDir = fullfile(testCase.RepoRoot, '_apps_exported');
            toolsDir = fullfile(testCase.RepoRoot, '_mlapp_tools');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(toolsDir));
        end
    end

    methods (TestMethodSetup)
        function resetTracking(testCase)
            testCase.TempFiles = {};
            testCase.TestAppNames = {};
        end
    end

    methods (TestMethodTeardown)
        function cleanupArtifacts(testCase)
            for k = 1:numel(testCase.TempFiles)
                f = testCase.TempFiles{k};
                if isfile(f)
                    delete(f);
                end
            end
            for k = 1:numel(testCase.TestAppNames)
                testCase.removeAppArtifacts(testCase.TestAppNames{k});
            end
        end
    end

    methods (Test)
        function testParseDataViewerEventsManagerExport(testCase)
        %TESTPARSEDATAVIEWEREVENTSMANAGEREXPORT Test Case 1: parse a real export.

            mlappPath = fullfile(testCase.RepoRoot, testCase.SourceApp);
            testCase.assumeTrue(isfile(mlappPath), 'Reference app not found in this checkout.');

            tempM = testCase.newTempFile('.m');
            exportMlappToM(mlappPath, tempM);

            parsed = parseMatlabExport(tempM);

            hasCreateComponents = any(cellfun( ...
                @(m) strcmpi(m.name, 'createComponents'), parsed.methods));
            testCase.verifyFalse(hasCreateComponents, ...
                'createComponents must be skipped, not returned as a method.');

            structuralFailure = any(contains(parsed.warnings, 'failed structurally'));
            testCase.verifyFalse(structuralFailure, 'Parser must not fail structurally on a real export.');

            testCase.verifyGreaterThanOrEqual(numel(parsed.properties), 5);
            testCase.verifyGreaterThanOrEqual(numel(parsed.methods), 5);

            fprintf('[Test1] properties=%d methods=%d callbacks=%d warnings=%d\n', ...
                numel(parsed.properties), numel(parsed.methods), numel(parsed.callbacks), ...
                numel(parsed.warnings));
        end

        function testDeltaOnUnchangedApp(testCase)
        %TESTDELTAONUNCHANGEDAPP Test Case 2: delta against itself is empty.

            appName = testCase.registerTestApp();
            mlappPath = testCase.copySourceAppAs(appName);

            s1 = exportMatlabAppToDelta(appName, mlappPath); %#ok<NASGU> bootstrap run
            s2 = exportMatlabAppToDelta(appName, mlappPath);

            testCase.verifyFalse(s2.referenceUpdated, ...
                'Second run must diff against the existing reference, not recreate it.');
            testCase.verifyEqual(s2.changesSummary.numNew, 0);
            testCase.verifyEqual(s2.changesSummary.numModified, 0);
            testCase.verifyEqual(s2.changesSummary.numMissing, 0);

            testCase.verifyTrue(isfile(s2.deltaFilePath));
            info = dir(s2.deltaFilePath);
            testCase.verifyLessThan(info.bytes, 1024, ...
                'An unchanged-app delta should only contain headers.');

            fprintf('[Test2] delta size = %d bytes (< 1024 required)\n', info.bytes);
        end

        function testDeltaOnSimulatedChange(testCase)
        %TESTDELTAONSIMULATEDCHANGE Test Case 3: mutate parsed structs directly
        % and confirm diffMatlabBlocks/generateDeltaFile surface exactly the
        % changed items.

            mlappPath = fullfile(testCase.RepoRoot, testCase.SourceApp);
            testCase.assumeTrue(isfile(mlappPath), 'Reference app not found in this checkout.');

            tempM = testCase.newTempFile('.m');
            exportMlappToM(mlappPath, tempM);

            parsedRef = parseMatlabExport(tempM);
            parsedChanged = parsedRef;

            testCase.assumeGreaterThanOrEqual(numel(parsedChanged.methods), 1);
            bodyLengths = cellfun(@(m) numel(m.body), parsedChanged.methods);
            [~, shortestIdx] = min(bodyLengths);
            modifiedName = parsedChanged.methods{shortestIdx}.name;
            parsedChanged.methods{shortestIdx}.body = [parsedChanged.methods{shortestIdx}.body, ...
                sprintf('\n            %% simulated change comment')];

            newCallback = struct('name', 'test_NewCallback_Callback', ...
                'attributes', '(Access = private)', ...
                'signature', 'function test_NewCallback_Callback(app, event)', ...
                'startLine', -1, 'endLine', -1, ...
                'body', '            disp(''simulated'');', 'isCallback', true);
            parsedChanged.methods{end+1} = newCallback;

            delta = diffMatlabBlocks(parsedRef, parsedChanged);

            testCase.verifyTrue(ismember('test_NewCallback_Callback', delta.newItems));
            testCase.verifyTrue(ismember(modifiedName, delta.modifiedItems));

            deltaFile = testCase.newTempFile('.txt');
            generateDeltaFile(delta, deltaFile, ...
                struct('appName', 'DataViewer_EventsManager', 'referenceInfo', 'unit-test'));

            text = fileread(deltaFile);
            testCase.verifyFalse(contains(text, 'createComponents'), ...
                'Delta file must never mention createComponents.');
            testCase.verifyTrue(contains(text, 'test_NewCallback_Callback'));
            testCase.verifyTrue(contains(text, modifiedName));

            info = dir(deltaFile);
            testCase.verifyLessThan(info.bytes, 5 * 1024);

            fprintf('[Test3] delta size = %d bytes (< 5120 required); new=%s modified=%s\n', ...
                info.bytes, 'test_NewCallback_Callback', modifiedName);
        end

        function testEndToEndWorkflow(testCase)
        %TESTENDTOENDWORKFLOW Test Case 4: bootstrap -> real zip-level app
        % mutation -> delta detects it -> confirmAndUpdateReference ->
        % delta is empty again.

            appName = testCase.registerTestApp();
            mlappPath = testCase.copySourceAppAs(appName);

            s1 = exportMatlabAppToDelta(appName, mlappPath);
            testCase.verifyTrue(s1.referenceUpdated);
            testCase.verifyEqual(s1.changesSummary.numNew, 0);
            testCase.verifyEqual(s1.changesSummary.numModified, 0);

            testCase.mutateMlappInPlace(mlappPath);

            s2 = exportMatlabAppToDelta(appName, mlappPath);
            testCase.verifyFalse(s2.referenceUpdated);
            testCase.verifyGreaterThanOrEqual(s2.changesSummary.numNew, 1);
            testCase.verifyGreaterThanOrEqual(s2.changesSummary.numModified, 1);

            freshExport = testCase.newTempFile('.m');
            exportMlappToM(mlappPath, freshExport);
            confirmResult = confirmAndUpdateReference(appName, freshExport);
            testCase.verifyEqual(confirmResult.status, 'success');

            s3 = exportMatlabAppToDelta(appName, mlappPath);
            testCase.verifyFalse(s3.referenceUpdated);
            testCase.verifyEqual(s3.changesSummary.numNew, 0);
            testCase.verifyEqual(s3.changesSummary.numModified, 0);
            testCase.verifyEqual(s3.changesSummary.numMissing, 0);

            fprintf('[Test4] bootstrap->change->confirm->clean cycle verified for "%s"\n', appName);
        end

        function testDataViewerManagementMenuModelMatchesRuntimeSource(testCase)
        %TESTDATAVIEWERMANAGEMENTMENUMODELMATCHESRUNTIMESOURCE Guard MLAPP metadata.

            appPath = fullfile(testCase.RepoRoot, ...
                'GUI', 'DataViewer', 'DataViewer.mlapp');
            testCase.assumeTrue(isfile(appPath), 'DataViewer app not found in this checkout.');

            [codeData, ~, ~, appData] = ...
                appdesigner.internal.comparison.getAppData(appPath);
            figureHandle = appData.components.UIFigure;
            cleanupFigure = onCleanup(@() delete(figureHandle));
            components = findall(figureHandle, '-property', 'DesignTimeProperties');
            codeNames = string(arrayfun(@(component) ...
                component.DesignTimeProperties.CodeName, components, ...
                'UniformOutput', false));

            expectedComponents = ["ManagementMenu", "SetRawFolderMenu", ...
                "ProjectMenu", "BindCurrentSessionMenu", "ProjectManagerMenu", ...
                "RigMenu", "AssignRigMenu", "RigManagerMenu"];
            testCase.verifyTrue(all(ismember(expectedComponents, codeNames)), ...
                'App Designer metadata must contain the complete Management hierarchy.');
            testCase.verifyFalse(any(ismember( ...
                ["ProjectManagerLabel", "ProjectManagerToolButton"], codeNames)), ...
                'Removed Project Manager launcher components must not remain in metadata.');

            management = components(codeNames == "ManagementMenu");
            setRaw = components(codeNames == "SetRawFolderMenu");
            project = components(codeNames == "ProjectMenu");
            rig = components(codeNames == "RigMenu");
            testCase.verifyEqual(setRaw.Parent, management);
            testCase.verifyEqual(project.Parent, management);
            testCase.verifyEqual(rig.Parent, management);

            callbackNames = string({codeData.Callbacks.Name});
            expectedCallbacks = ["ManagementMenuSelected", ...
                "AssignRigMenuSelected", "BindCurrentSessionMenuSelected", ...
                "ProjectManagerMenuSelected", "RigManagerMenuSelected"];
            testCase.verifyTrue(all(ismember(expectedCallbacks, callbackNames)));
            testCase.verifyFalse(ismember( ...
                "ProjectManagerToolButtonPushed", callbackNames));

            editableText = strjoin(string(codeData.EditableSectionCode), newline);
            testCase.verifyTrue(all(contains(editableText, ...
                ["function refreshRigAssignmentMenu", ...
                "function SelectRigMenuSelected", ...
                "function openProjectManagerTool"])));
        end

        function testDataViewerGraphicalCallbacksAreBound(testCase)
        %TESTDATAVIEWERGRAPHICALCALLBACKSAREBOUND Guard component bindings.

            appPath = fullfile(testCase.RepoRoot, ...
                'GUI', 'DataViewer', 'DataViewer.mlapp');
            report = auditMlappCallbackBindings(appPath);
            testCase.verifyEmpty(report.UnboundCallbacks);
            testCase.verifyEmpty(report.UnknownBindings);
            testCase.verifyEmpty(report.DuplicateCallbackOwners);
            testCase.verifyEqual(report.BindingCount, report.CallbackCount);

            bindings = string({report.Bindings.Component}) + "." + ...
                string({report.Bindings.Property}) + "=" + ...
                string({report.Bindings.Callback});
            expected = [ ...
                "Slider.ValueChangedFcn=SliderValueChanged"
                "Slider.ValueChangingFcn=SliderValueChanging"
                "PreviousFrameButton.ButtonPushedFcn=PreviousFrameButtonPushed"
                "NextFrameButton.ButtonPushedFcn=NextFrameButtonPushed"
                "PlayMovieButton.ButtonPushedFcn=PlayMovieButtonPushed"
                "MovieSpeedDropDown.ValueChangedFcn=MovieSpeedDropDownValueChanged"
                "ColormapDropDown.ValueChangedFcn=ColormapDropDownValueChanged"
                "InvertCheckBox.ValueChangedFcn=InvertCheckBoxValueChanged"
                "ClipSliderRange.ValueChangedFcn=ClipSliderRangeValueChanged"
                "ClipSliderRange.ValueChangingFcn=ClipSliderRangeValueChanging"
                "AutoButton.ButtonPushedFcn=AutoButtonPushed"
                "SetClipButton.ButtonPushedFcn=SetClipButtonPushed"
                "HidecrosshairCheckBox.ValueChangedFcn=HidecrosshairCheckBoxValueChanged"
                "Switch.ValueChangedFcn=SwitchValueChanged"
                "ConditionDropDown.ValueChangedFcn=ConditionDropDownValueChanged"
                "RepetitionDropDown.ValueChangedFcn=RepetitionDropDownValueChanged"
                "DeleteConditionButton.ButtonPushedFcn=DeleteConditionButtonPushed"
                "DeleteRepetitionButton.ButtonPushedFcn=DeleteRepetitionButtonPushed"
                "RestoreButton.ButtonPushedFcn=RestoreButtonPushed"];
            testCase.verifyTrue(all(ismember(expected, bindings)));
        end

        function testSafeMlappSaveIsIdempotent(testCase)
        %TESTSAFEMLAPPSAVEISIDEMPOTENT Guard serializer-derived metadata.

            sourcePath = fullfile(testCase.RepoRoot, ...
                'GUI', 'ProjectManagerTool.mlapp');
            testCase.assumeTrue(isfile(sourcePath), ...
                'ProjectManagerTool app not found in this checkout.');
            appPath = testCase.newTempFile('.mlapp');
            copyfile(sourcePath, appPath);
            expectedSource = exportMlappToM(appPath);

            for iPass = 1:2
                [codeData, compatibilityData, metadata, appData] = ...
                    appdesigner.internal.comparison.getAppData(appPath);
                figureHandle = appData.components.UIFigure;
                cleanupFigure = onCleanup(@() delete(figureHandle));
                report = saveMlappSafely(appPath, figureHandle, ...
                    expectedSource, codeData, compatibilityData, ...
                    metadata, appData);
                testCase.verifyTrue(report.IsValid);
                testCase.verifyEqual(exportMlappToM(appPath), expectedSource);
                testCase.verifyEmpty( ...
                    report.CallbackAudit.DuplicateCallbackOwners);
                clear cleanupFigure
            end
        end

        function testDataViewerOverwriteImportForcesSelectedSteps(testCase)
        %TESTDATAVIEWEROVERWRITEIMPORTFORCESSELECTEDSTEPS Guard stale-file handling.

            helperBody = testCase.getDataViewerMethodBody( ...
                'executeDataImportPipeline', 'updateDataImportProgress');
            binCallback = testCase.getDataViewerMethodBody( ...
                'frombinMenuSelected', 'fromtifMenuSelected');

            testCase.verifySubstring(helperBody, 'ppln.b_skipSteps = false;', ...
                'DataViewer imports must ignore stale destination cache state.');
            testCase.verifySubstring(binCallback, 'ppln.addStep(''getEvents'')');
            testCase.verifySubstring(binCallback, ...
                'app.executeDataImportPipeline(ppln, ''Data import'')');

            saveFolder = tempname;
            mkdir(saveFolder);
            cleanupFolder = onCleanup(@() rmdir(saveFolder, 's'));

            pm = PipelineManager(saveFolder, saveFolder, testCase.RepoRoot);
            testCase.verifyTrue(pm.b_skipSteps, ...
                'Ordinary PipelineManager runs must retain cache-based skipping by default.');
        end

        function testDataViewerImportUsesModalProgressAndCleanup(testCase)
        %TESTDATAVIEWERIMPORTUSESMODALPROGRESSANDCLEANUP Guard interaction blocking.

            helperBody = testCase.getDataViewerMethodBody( ...
                'executeDataImportPipeline', 'frombinMenuSelected');

            testCase.verifySubstring(helperBody, 'uiprogressdlg(app.UIFigure');
            testCase.verifySubstring(helperBody, '''Cancelable'', ''on''');
            testCase.verifySubstring(helperBody, '''ProgressFcn''');
            testCase.verifySubstring(helperBody, '''CancelFcn''');
            testCase.verifySubstring(helperBody, 'onCleanup');
            testCase.verifySubstring(helperBody, 'close(progressDlg)');
            testCase.verifySubstring(helperBody, 'app.setInteractionMode(previousMode)');
        end

        function testDataViewerImportCallbacksHandleCancellation(testCase)
        %TESTDATAVIEWERIMPORTCALLBACKSHANDLECANCELLATION Guard terminal status paths.

            binCallback = testCase.getDataViewerMethodBody( ...
                'frombinMenuSelected', 'fromtifMenuSelected');
            tifCallback = testCase.getDataViewerMethodBody( ...
                'fromtifMenuSelected', 'totifMenuSelected');

            testCase.verifySubstring(binCallback, ...
                'strcmpi(string(importResult.status), ''cancelled'')');
            testCase.verifySubstring(binCallback, 'Data import cancelled.');

            testCase.verifySubstring(tifCallback, ...
                'strcmpi(string(importResult.status), ''cancelled'')');
            testCase.verifySubstring(tifCallback, 'TIFF import cancelled.');
        end

        function testDataViewerDataSwitchValidatesRetainedROIs(testCase)
        %TESTDATAVIEWERDATASWITCHVALIDATESRETAINEDROIS Guard ROI state at file switches.

            loadBody = testCase.getDataViewerMethodBody( ...
                'loadDataSource', 'loadDataParamsForCurrentFile');
            switchBody = testCase.getDataViewerMethodBody( ...
                'reconcileROIsForNewDataSize', 'promptROIDataSwitchSizeAction');

            testCase.verifySubstring(loadBody, 'previousSizeYX');
            testCase.verifySubstring(loadBody, ...
                'app.reconcileROIsForNewDataSize(previousSizeYX)');
            testCase.verifySubstring(switchBody, ...
                'app.rebuildLoadedROIListForCurrentImage');
            testCase.verifySubstring(switchBody, 'app.clearAllROIsWithoutPrompt()');
            testCase.verifySubstring(switchBody, 'report.nSkipped');
        end

        function testDataViewerROISelectionDisablesProjectManagerTool(testCase)
        %TESTDATAVIEWERROISELECTIONDISABLESPROJECTMANAGERTOOL Guard edit-mode UI state.

            guiBody = testCase.getDataViewerMethodBody( ...
                'updateGUIEnabledState', 'setContainerEnabled');

            testCase.verifySubstring(guiBody, 'app.ProjectManagerToolButton');
            testCase.verifySubstring(guiBody, 'caps.hasData && isIdle');
        end

        function testDataViewerCropsROIsToEffectiveImageMask(testCase)
        %TESTDATAVIEWERCROPSROISTOEFFECTIVEIMAGEMASK Guard image-boundary clipping.

            createBody = testCase.getDataViewerMethodBody( ...
                'createNewROI', 'selectROIForEditingByID');
            editBody = testCase.getDataViewerMethodBody( ...
                'commitROIEditingByIndex', 'cancelROIEditingByIndex');
            groupBody = testCase.getDataViewerMethodBody( ...
                'commitGroupROIEdit', 'approximateGroupEditedROIPrimitive');
            rebuildBody = testCase.getDataViewerMethodBody( ...
                'rebuildOneLoadedROIForCurrentImage', 'normalizeLoadedROIRecord');
            boundsBody = testCase.getDataViewerMethodBody( ...
                'roiVerticesExtendOutsideImage', 'clipROIMaskToActiveLogicalMask');

            testCase.verifySubstring(createBody, 'app.roiVerticesExtendOutsideImage');
            testCase.verifySubstring(editBody, 'app.roiVerticesExtendOutsideImage');
            testCase.verifySubstring(groupBody, 'app.roiVerticesExtendOutsideImage');
            testCase.verifySubstring(rebuildBody, 'app.roiVerticesExtendOutsideImage');
            testCase.verifySubstring(boundsBody, 'xBounds = [0.5');
            testCase.verifySubstring(boundsBody, 'yBounds = [0.5');
        end
    end

    methods (Access = private)
        function appName = registerTestApp(testCase)
            appName = sprintf('%s%s', testCase.AppNamePrefix, char(java.util.UUID.randomUUID));
            appName = strrep(appName, '-', '_');
            testCase.TestAppNames{end+1} = appName;
        end

        function f = newTempFile(testCase, ext)
            f = fullfile(tempdir, [char(java.util.UUID.randomUUID) ext]);
            testCase.TempFiles{end+1} = f;
        end

        function mlappPath = copySourceAppAs(testCase, appName)
            src = fullfile(testCase.RepoRoot, testCase.SourceApp);
            testCase.assumeTrue(isfile(src), 'Reference app not found in this checkout.');
            mlappPath = fullfile(tempdir, [appName '.mlapp']);
            copyfile(src, mlappPath);
            testCase.TempFiles{end+1} = mlappPath;
        end

        function body = getDataViewerMethodBody(testCase, methodName, nextMethodName)
        %GETDATAVIEWERMETHODBODY Return one focused method span from DataViewer.

            appPath = fullfile(testCase.RepoRoot, ...
                'GUI', 'DataViewer', 'DataViewer.mlapp');
            code = exportMlappToM(appPath);

            startPattern = ['function[^\r\n]*\<' methodName '\s*\('];
            startIndex = regexp(code, startPattern, 'start');
            testCase.assertNumElements(startIndex, 1, ...
                sprintf('Expected one %s method.', methodName));

            endPattern = ['function[^\r\n]*\<' nextMethodName '\s*\('];
            endIndex = regexp(code, endPattern, 'start');
            endIndex = endIndex(endIndex > startIndex);
            testCase.assertNotEmpty(endIndex, ...
                sprintf('Could not find the method following %s.', methodName));

            body = code(startIndex:endIndex(1)-1);
        end

        function removeAppArtifacts(testCase, appName)
            d = testCase.AppsExportedDir;
            refFile = fullfile(d, [appName '_reference.m']);
            if isfile(refFile)
                delete(refFile);
            end
            deltaFiles = dir(fullfile(d, [appName '_delta_*.txt']));
            for k = 1:numel(deltaFiles)
                delete(fullfile(deltaFiles(k).folder, deltaFiles(k).name));
            end
            archiveFiles = dir(fullfile(d, 'archive', [appName '_reference_*.m']));
            for k = 1:numel(archiveFiles)
                delete(fullfile(archiveFiles(k).folder, archiveFiles(k).name));
            end
        end

        function mutateMlappInPlace(~, mlappPath)
        %MUTATEMLAPPINPLACE Rewrite the CDATA code block inside a .mlapp's
        % matlab/document.xml to add a callback and tweak an existing
        % method, then rezip in place. Used only on throwaway temp copies.

            scratchDir = fullfile(tempdir, ['mlapp_mutate_' char(java.util.UUID.randomUUID)]);
            mkdir(scratchDir);
            cleanupObj = onCleanup(@() rmdir(scratchDir, 's'));

            unzip(mlappPath, scratchDir);
            docPath = fullfile(scratchDir, 'matlab', 'document.xml');
            xml = fileread(docPath);
            tok = regexp(xml, '<!\[CDATA\[(.*)\]\]>', 'tokens', 'once');
            if isempty(tok)
                error('MLAppDelta_tests:MutateFailed', 'No CDATA block found to mutate.');
            end
            code = tok{1};

            classNameTok = regexp(code, '^classdef\s+(\w+)', 'tokens', 'once');
            if isempty(classNameTok)
                error('MLAppDelta_tests:MutateFailed', 'Could not determine classdef name.');
            end
            marker = sprintf('function app = %s(varargin)', classNameTok{1});
            markerIdx = strfind(code, marker);
            if isempty(markerIdx)
                error('MLAppDelta_tests:MutateFailed', 'No insertion point found.');
            end

            insertion = sprintf(['        function test_NewCallback_Callback(app, event)\n' ...
                '            disp(''hello from simulated callback'');\n' ...
                '        end\n\n        ']);
            mutatedCode = [code(1:markerIdx(1)-1) insertion code(markerIdx(1):end)];

            % Tweak one existing method body to also exercise "modified" detection.
            knownLine = 'selectedChannels = {app.TriggerChannelCheckBoxes(checked).Text};';
            if contains(mutatedCode, knownLine)
                mutatedCode = strrep(mutatedCode, knownLine, [knownLine ' % simulated tweak']);
            else
                error('MLAppDelta_tests:MutateFailed', ...
                    'Expected line to tweak was not found — app source may have changed.');
            end

            enc = strrep(mutatedCode, '&', '&amp;');
            enc = strrep(enc, '<', '&lt;');
            enc = strrep(enc, '>', '&gt;');
            newXml = strrep(xml, tok{1}, enc);

            fid = fopen(docPath, 'w');
            fprintf(fid, '%s', newXml);
            fclose(fid);

            zipOut = fullfile(tempdir, ['mutated_' char(java.util.UUID.randomUUID) '.zip']);
            origDir = pwd;
            cd(scratchDir);
            cleanupCd = onCleanup(@() cd(origDir));
            zip(zipOut, {'[Content_Types].xml', '_rels', 'appdesigner', 'matlab', 'metadata'});
            cd(origDir);

            delete(mlappPath);
            movefile(zipOut, mlappPath);
        end
    end
end
