classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) MLAppDelta_tests < matlab.unittest.TestCase
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
                "RigMenu", "FolderRigMenu", "RigManagerMenu"];
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
                "FolderRigMenuSelected", "BindCurrentSessionMenuSelected", ...
                "ProjectManagerMenuSelected", "RigManagerMenuSelected"];
            testCase.verifyTrue(all(ismember(expectedCallbacks, callbackNames)));
            testCase.verifyFalse(ismember( ...
                "ProjectManagerToolButtonPushed", callbackNames));

            editableText = strjoin(string(codeData.EditableSectionCode), newline);
            testCase.verifyTrue(all(contains(editableText, ...
                ["function refreshRigAssignmentMenu", ...
                "function SelectRigMenuSelected", ...
                "function openProjectManagerTool"])));
            % The component was renamed to FolderRigMenu; no method body may
            % still reference the old name.
            testCase.verifyFalse(contains(editableText, "AssignRigMenu"));
        end

        function testDataViewerGraphicalCallbacksAreBound(testCase)
        %TESTDATAVIEWERGRAPHICALCALLBACKSAREBOUND Guard component bindings.

            appPath = fullfile(testCase.RepoRoot, ...
                'GUI', 'DataViewer', 'DataViewer.mlapp');
            report = auditMlappCallbackBindings(appPath);
            testCase.verifyEmpty(report.UnboundCallbacks);
            testCase.verifyEmpty(report.UnknownBindings);
            testCase.verifyEmpty(report.DuplicateCallbackOwners);

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
                "RestoreButton.ButtonPushedFcn=RestoreButtonPushed"
                "BatchPipelineMenu.MenuSelectedFcn=BatchPipelineMenuSelected"
                "PipeLauncherBatchMenu.MenuSelectedFcn=BatchPipelineMenuSelected"];
            testCase.verifyTrue(all(ismember(expected, bindings)));

            % BatchPipelineMenuSelected is deliberately shared by the File
            % menu entry and the pipeline launcher context menu, so there is
            % exactly one more binding than there are callbacks.
            sharedCallbacks = "BatchPipelineMenuSelected";
            callbackOfBinding = string({report.Bindings.Callback});
            testCase.verifyEqual(report.BindingCount, ...
                report.CallbackCount + nnz(callbackOfBinding == sharedCallbacks) - 1, ...
                'Only BatchPipelineMenuSelected may be bound by several components.');
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
                'reconcileROIsForNewDataSize', 'roiVerticesExtendOutsideImage');

            testCase.verifySubstring(loadBody, 'previousSizeYX');
            testCase.verifySubstring(loadBody, ...
                'app.reconcileROIsForNewDataSize(previousSizeYX)');
            testCase.verifySubstring(switchBody, ...
                'app.rebuildLoadedROIListForCurrentImage');
            testCase.verifySubstring(switchBody, 'app.clearAllROIsWithoutPrompt()');
            testCase.verifySubstring(switchBody, 'report.nSkipped');
        end

        function testDataViewerROISelectionGatesLaunchersOnIdle(testCase)
        %TESTDATAVIEWERROISELECTIONGATESLAUNCHERSONIDLE Guard edit-mode UI state.

            guiBody = testCase.getDataViewerMethodBody( ...
                'updateGUIEnabledState', 'setContainerEnabled');

            % The Project Manager launcher button was removed; the launchers
            % that remain must still be gated on the idle state.
            testCase.verifyFalse(contains(guiBody, 'ProjectManagerToolButton'));
            testCase.verifySubstring(guiBody, ...
                'app.setComponentEnabled(app.BatchPipelineMenu, isIdle)');
            testCase.verifySubstring(guiBody, 'caps.hasData && isIdle');
        end

        function testPipelineManagerToolIsNonModalAndDataViewerBlocksItself(testCase)
        %TESTPIPELINEMANAGERTOOLISNONMODALANDDATAVIEWERBLOCKSITSELF Guard PM Tool blocking.
        %
        % A modal PM Tool hides its menu bar, so DataViewer blocks itself with a
        % progress dialog instead. Every tool-close and launch-failure path must
        % release that dialog or DataViewer stays locked.

            toolPath = fullfile(testCase.RepoRoot, 'GUI', 'PipelineManagerTool.mlapp');
            startupBody = testCase.getMlappMethodBody(toolPath, 'startupFcn');
            testCase.verifyFalse(contains(startupBody, 'WindowStyle'), ...
                'PipelineManagerTool must not set a modal WindowStyle at startup.');

            dataViewerPath = fullfile(testCase.RepoRoot, ...
                'GUI', 'DataViewer', 'DataViewer.mlapp');
            blockBody = testCase.getMlappMethodBody(dataViewerPath, ...
                'blockForPipelineManager');
            testCase.verifySubstring(blockBody, 'uiprogressdlg(app.UIFigure');
            testCase.verifySubstring(blockBody, '''Indeterminate'', ''on''');

            for launcher = ["PipeLauncherButtonPushed", "BatchPipelineMenuSelected"]
                body = testCase.getMlappMethodBody(dataViewerPath, launcher);
                testCase.verifySubstring(body, 'app.blockForPipelineManager()');
                testCase.verifySubstring(body, 'app.guardPipelineManagerRelease()');
                testCase.verifySubstring(body, 'app.releasePipelineManagerBlock()', ...
                    sprintf('%s must release the block when launching fails.', launcher));
            end

            for closedFcn = ["onPipelineManagerClosed", "onBatchPipelineManagerClosed"]
                body = testCase.getMlappMethodBody(dataViewerPath, closedFcn);
                testCase.verifySubstring(body, 'app.releasePipelineManagerBlock()');
            end

            destroyedBody = testCase.getMlappMethodBody(dataViewerPath, ...
                'onPipelineManagerFigureDestroyed');
            testCase.verifySubstring(destroyedBody, 'app.releasePipelineManagerBlock()');
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
                'roiVerticesExtendOutsideImage', 'alertNoSupportedROIImport');

            % ROI creation/edit/rebuild paths rasterize into image coordinates
            % and clip through the active logical mask; they no longer call
            % roiVerticesExtendOutsideImage directly.
            testCase.verifySubstring(createBody, 'app.clipROIMaskToActiveLogicalMask');
            testCase.verifySubstring(editBody, 'app.clipROIMaskToActiveLogicalMask');
            testCase.verifySubstring(groupBody, 'app.clipROIMaskToActiveLogicalMask');
            testCase.verifySubstring(rebuildBody, 'app.clipROIMaskToActiveLogicalMask');
            testCase.verifySubstring(boundsBody, 'xBounds = [0.5');
            testCase.verifySubstring(boundsBody, 'yBounds = [0.5');
        end

        function testParserBalancesNestingAroundTrickyConstructs(testCase)
        %TESTPARSERBALANCESNESTINGAROUNDTRICKYCONSTRUCTS Regression: one-line
        % blocks (`if c; return; end`) used to leave an unclosed block, so every
        % method after them (most of DataViewer) was silently dropped.

            src = {
                'classdef Tricky'
                '    properties (Access = public)'
                '        Items cell = {1, ''end''}'
                '    end'
                '    methods (Access = private)'
                '        function oneLiner(app, h)'
                '            if strcmp(h.Color,''none'');return;end'
                '            for k = 1:3, app.Items{end+1} = k; end'
                '            try, x = 1; catch, x = 2; end'
                '        end'
                '        function endIndexAndStrings(app, v)'
                '            last = v(end);'
                '            f = @(q) q(end);'
                '            fprintf(''100% if for end while'');'
                '            s = "if end switch";  % end trailing comment'
                '            w = v'';  % transpose then quote'
                '            m = [1 2 ...'
                '                 3 end];'
                '        end'
                '        function withArguments(app, a)'
                '            arguments'
                '                app'
                '                a (1,1) double = 0'
                '            end'
                '            switch a'
                '                case 1'
                '                    disp(1)'
                '                otherwise'
                '                    disp(0)'
                '            end'
                '        end'
                '%{'
                'end end end'
                '%}'
                '        function lastOne(app)'
                '            while true, break, end'
                '        end'
                '    end'
                'end'};
            mPath = testCase.newTempFile('.m');
            fid = fopen(mPath, 'w');
            fwrite(fid, strjoin(src, newline));
            fclose(fid);

            parsed = parseMatlabExport(mPath);

            testCase.verifyEmpty(parsed.warnings, ...
                'Tricky but valid constructs must not trigger unbalanced-nesting warnings.');
            names = cellfun(@(m) m.name, parsed.methods, 'UniformOutput', false);
            testCase.verifyEqual(names, ...
                {'oneLiner', 'endIndexAndStrings', 'withArguments', 'lastOne'});
            testCase.verifyEqual(numel(parsed.properties), 1);
        end

        function testExportKeepsXmlEntitiesInCdataVerbatim(testCase)
        %TESTEXPORTKEEPSXMLENTITIESINCDATAVERBATIM Regression: CDATA is literal,
        % so `strrep(t, '&', '&amp;')` must not be decoded to `'&', '&'`
        % (an unterminated-string corruption that broke parsing downstream).

            workDir = tempname;
            mkdir(fullfile(workDir, 'matlab'));
            cleanup = onCleanup(@() rmdir(workDir, 's'));
            code = sprintf(['classdef E\n    methods\n        function f(app, t)\n' ...
                '            t = strrep(t, ''&'', ''&amp;'');\n        end\n' ...
                '        function g(app)\n        end\n    end\nend\n']);
            fid = fopen(fullfile(workDir, 'matlab', 'document.xml'), 'w');
            fwrite(fid, ['<w:document><![CDATA[' code ']]></w:document>']);
            fclose(fid);
            zipPath = [tempname '.zip'];
            zip(zipPath, 'matlab', workDir);
            mlappPath = [tempname '.mlapp'];
            movefile(zipPath, mlappPath);
            testCase.TempFiles{end+1} = mlappPath;

            exported = exportMlappToM(mlappPath);

            testCase.verifySubstring(exported, '''&amp;''');
            mPath = testCase.newTempFile('.m');
            fid = fopen(mPath, 'w');
            fwrite(fid, exported);
            fclose(fid);
            parsed = parseMatlabExport(mPath);
            testCase.verifyEmpty(parsed.warnings);
            testCase.verifyEqual(numel(parsed.methods), 2);
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

        function body = getMlappMethodBody(testCase, appPath, methodName)
        %GETMLAPPMETHODBODY Return one method span, ending at the next function.
        %
        % Only valid for methods without nested functions.

            code = exportMlappToM(appPath);

            startPattern = ['function[^\r\n]*\<' char(methodName) '\s*\('];
            startIndex = regexp(code, startPattern, 'start');
            testCase.assertNumElements(startIndex, 1, ...
                sprintf('Expected one %s method.', methodName));

            nextIndex = regexp(code, '\n\s*function\s', 'start');
            nextIndex = nextIndex(nextIndex > startIndex);
            if isempty(nextIndex)
                body = code(startIndex:end);
            else
                body = code(startIndex:nextIndex(1));
            end
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
