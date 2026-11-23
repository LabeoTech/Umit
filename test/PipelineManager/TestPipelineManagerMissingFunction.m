classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerMissingFunction < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERMISSINGFUNCTION A step whose function no longer exists is invalid.
    %
    %   A saved .pipe can name a function that was renamed or removed after
    %   the pipeline was saved (for example "run_ImagesClassification" after
    %   the rename to "ImagesClassification"). loadPipe already warns about it;
    %   diagnosePipeline and executePipeline must agree with that verdict
    %   instead of reporting the pipeline as valid and failing mid-run.
    %
    %   The saved pipelines are built from funcTemplate (copied into a
    %   temporary Analysis category) with the function name rewritten, so the
    %   tests need no dataset.

    properties
        ProjectRoot char
        SaveFolder char
        PipeFolder char
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);

            % funcTemplate plus a renamed copy: a branch cannot hold the same
            % function twice, and the chain tests need two different steps.
            templateFile = fullfile(testCase.ProjectRoot, 'Analysis', 'funcTemplate.m');
            categoryFolder = createAnalysisFixtureCategory(testCase.ProjectRoot, ...
                struct('name', {'funcTemplate', 'funcTemplateB'}, ...
                'source', {templateFile, ''}, ...
                'text', {'', strrep(fileread(templateFile), 'funcTemplate', 'funcTemplateB')}));
            testCase.addTeardown(@() removeAnalysisFixtureCategory( ...
                categoryFolder, testCase.ProjectRoot));

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            testCase.PipeFolder = fixture.Folder;
            mkdir(testCase.SaveFolder);

            inputData = single(reshape(-12:11, [2 3 4]));
            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                'FrameRateHz', 20, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                'DatFile', 'input.dat', 'Length', 4, 'FrameRateHz', 20));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            saveData(fullfile(testCase.SaveFolder, 'input.dat'), inputData, ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20);
        end
    end

    methods (Test)
        function availableFunctionKeepsPipelineValid(testCase)
            pm = testCase.newManagerWithStep();
            report = pm.diagnosePipeline();

            testCase.verifyTrue(report.isValid);
            testCase.verifyTrue(report.functionsAreValid);
            testCase.verifyFalse(report.hasInvalidBlocking);

            pipeFile = fullfile(testCase.PipeFolder, 'current.pipe');
            pm.savePipe(pipeFile);
            restored = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            testCase.verifyWarningFree(@() restored.loadPipe(pipeFile));
            testCase.verifyTrue(restored.diagnosePipeline().isValid);
        end

        function missingFunctionMakesLoadedPipelineInvalid(testCase)
            restored = testCase.loadPipeWithRenamedFunction('run_funcTemplate');
            report = restored.diagnosePipeline();

            testCase.verifyFalse(report.isValid);
            testCase.verifyFalse(report.functionsAreValid);
            testCase.verifyTrue(report.graphIsValid, ...
                'The graph itself is fine; only the function is missing.');
            testCase.verifyTrue(report.hasInvalidBlocking);
            testCase.verifySubstring(strjoin(report.errors, newline), 'run_funcTemplate');

            row = report.nodeStatus(report.nodeStatus.functionName == "run_funcTemplate", :);
            testCase.assertEqual(height(row), 1);
            testCase.verifyEqual(row.status, "invalid_blocking");
        end

        function messageSuggestsTheCurrentFunction(testCase)
            restored = testCase.loadPipeWithRenamedFunction('run_funcTemplate');
            [tf, report] = restored.validateNodeFunctions();

            testCase.verifyFalse(tf);
            testCase.assertNumElements(report.errors, 1);
            testCase.verifySubstring(report.errors{1}, 'was not found');
            testCase.verifySubstring(report.errors{1}, 'Did you mean "funcTemplate"?');
            testCase.assertNumElements(report.nodeErrors, 1);
        end

        function unrelatedNameGetsNoSuggestion(testCase)
            restored = testCase.loadPipeWithRenamedFunction('removedStep');
            [tf, report] = restored.validateNodeFunctions();

            testCase.verifyFalse(tf);
            testCase.assertNumElements(report.errors, 1);
            testCase.verifySubstring(report.errors{1}, 'removedStep');
            testCase.verifyEmpty(strfind(report.errors{1}, 'Did you mean'));
        end

        function throwOnErrorRaisesMissingFunctionsError(testCase)
            restored = testCase.loadPipeWithRenamedFunction('run_funcTemplate');

            testCase.verifyError( ...
                @() restored.validateNodeFunctions('throwOnError', true), ...
                'PipelineManager:validateNodeFunctions:MissingFunctions');
        end

        function executePipelineRefusesMissingFunctionBeforeRunning(testCase)
            restored = testCase.loadPipeWithRenamedFunction('run_funcTemplate');
            restored.b_skipSteps = false;

            testCase.verifyError(@() restored.executePipeline('PrintSummary', false), ...
                'PipelineManager:validateNodeFunctions:MissingFunctions');
            testCase.verifyFalse(isfile(fullfile(testCase.SaveFolder, 'pipelineTemplate.dat')));
        end

        function removeInvalidNodesRemovesStepWithMissingFunction(testCase)
            restored = testCase.loadPipeWithRenamedFunction('run_funcTemplate');

            removed = restored.removeInvalidNodes();

            % The step and the file source it left unused.
            testCase.verifyNumElements(removed, 2);
            report = restored.diagnosePipeline();
            testCase.verifyTrue(report.isValid);
            testCase.verifyTrue(report.functionsAreValid);
            testCase.verifyEqual(height(report.nodeStatus), 0);
        end

        function removeInvalidNodesRemovesDownstreamStepThatLosesItsInput(testCase)
            % A is missing; B reads A, so B loses its input and goes too.
            restored = testCase.loadChainWithRenamedFunction('run_funcTemplate', 1);
            removed = restored.removeInvalidNodes();

            testCase.verifyNumElements(removed, 3, 'Both steps and the unused file source.');
            report = restored.diagnosePipeline();
            testCase.verifyTrue(report.isValid);
            testCase.verifyEqual(height(report.nodeStatus), 0);
        end

        function removeInvalidNodesKeepsUpstreamStepWhenDownstreamIsMissing(testCase)
            % B is missing; A (its source) stays and the pipeline is valid.
            restored = testCase.loadChainWithRenamedFunction('run_funcTemplateB', 2);

            removed = restored.removeInvalidNodes();

            testCase.verifyNumElements(removed, 1);
            report = restored.diagnosePipeline();
            testCase.verifyTrue(report.isValid);
            testCase.verifyEqual(nnz(report.nodeStatus.functionName == "funcTemplate"), 1);
            testCase.verifyEqual(nnz(report.nodeStatus.functionName == "run_funcTemplateB"), 0);
        end

        function removeInvalidNodesLeavesAValidPipelineUntouched(testCase)
            pm = testCase.newManagerWithStep();

            removed = pm.removeInvalidNodes();

            testCase.verifyEmpty(removed);
            testCase.verifyTrue(pm.diagnosePipeline().isValid);
        end
    end

    methods (Access = private)
        function pm = newManagerWithStep(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'pipelineTemplate.dat');
        end

        function pm = newChainedManager(testCase)
            %NEWCHAINEDMANAGER Two funcTemplate steps: B reads A.
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            tagA = pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'stepA.dat');
            pm.addStep('funcTemplateB', ...
                'input', sprintf('%s:outData', tagA), ...
                'saveas', 'stepB.dat');
        end

        function restored = loadPipeWithRenamedFunction(testCase, newName)
            %LOADPIPEWITHRENAMEDFUNCTION Save a funcTemplate pipeline, rename its function, reload.
            %
            %   Mimics a .pipe saved before the function was renamed or removed.

            restored = testCase.reloadWithRenamedStep( ...
                testCase.newManagerWithStep(), newName, 1);
        end

        function restored = loadChainWithRenamedFunction(testCase, newName, stepNumber)
            %LOADCHAINWITHRENAMEDFUNCTION Same, for the two-step chain (rename step 1 or 2).

            restored = testCase.reloadWithRenamedStep( ...
                testCase.newChainedManager(), newName, stepNumber);
        end

        function restored = reloadWithRenamedStep(testCase, pm, newName, stepNumber)
            %RELOADWITHRENAMEDSTEP Save PM, rename the function of one stream step, reload.

            pipeFile = fullfile(testCase.PipeFolder, sprintf('%s_%d.pipe', newName, stepNumber));
            pm.savePipe(pipeFile);

            loaded = load(pipeFile, 'pipeStruct', '-mat');
            pipeStruct = loaded.pipeStruct;
            streamIdx = find(strcmpi({pipeStruct.nodes.kind}, 'stream'));
            testCase.assertGreaterThanOrEqual(numel(streamIdx), stepNumber);
            idx = streamIdx(stepNumber);
            oldName = char(string(pipeStruct.nodes(idx).info.name));

            for fieldName = {'name', 'tag'}
                if isfield(pipeStruct.nodes, fieldName{1})
                    pipeStruct.nodes(idx).(fieldName{1}) = strrep( ...
                        char(string(pipeStruct.nodes(idx).(fieldName{1}))), ...
                        oldName, newName);
                end
            end
            pipeStruct.nodes(idx).info.name = newName;
            if isfield(pipeStruct, 'savedFuncName')
                % Swap exactly this step's function name in the saved list.
                saved = string(pipeStruct.savedFuncName(:));
                saved(saved == string(oldName)) = newName;
                pipeStruct.savedFuncName = saved;
            end
            save(pipeFile, 'pipeStruct', '-mat');

            restored = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            testCase.verifyWarning(@() restored.loadPipe(pipeFile), ...
                'PipelineManager:loadPipe:MissingFunctions');
        end
    end
end
