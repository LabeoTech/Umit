classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerDataViewerSourceRule < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERDATAVIEWERSOURCERULE Which file is the pipeline input when DataViewer reopens it.
    %
    %   chooseDataViewerSource decides between keeping the stored source (the pipeline
    %   stays a recipe rooted at the original file) and swapping to the file DataViewer
    %   shows. The tests apply the decision the way PipelineManagerTool does (swap =
    %   ensureFileSourceStep + reconnectStep + rmStep) and check what a run then
    %   produces, so a wrong decision shows as a result that was filtered twice.
    %
    %   Data is 1:24 as 2 x 3 x 4 (Y, X, T): non-square. pmSkipStepA / pmSkipStepB add 1
    %   and pmReuseGainStep multiplies by its Gain; all count their real executions.

    properties
        ProjectRoot char
        Folder char
        Base single
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.Folder);

            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                'FrameRateHz', 10, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                struct('DatFile', 'green.dat', 'Length', 4, 'FrameRateHz', 10));
            save(fullfile(testCase.Folder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.Base = single(reshape(1:24, [2 3 4]));
            saveData(fullfile(testCase.Folder, 'green.dat'), testCase.Base, ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));

            setappdata(0, 'pmSkipLog', struct());
            testCase.addTeardown(@() setappdata(0, 'pmSkipLog', struct()));
        end
    end

    methods (Test)
        function keepsTheOriginalWhenTheShownFileIsAResultOfThePipeline(testCase)
            [pm, shown] = testCase.firstRun({'pmSkipStepA', 'pmSkipStepB'});

            [action, reason] = pm.chooseDataViewerSource(shown, iSourceID(pm));

            testCase.verifyEqual(action, "keep");
            testCase.verifyEqual(reason, "currentIsResultOfStoredSource");
        end

        function keepsTheStoredSourceWhenItIsTheShownFile(testCase)
            [pm, ~] = testCase.firstRun({'pmSkipStepA'});

            [action, reason] = pm.chooseDataViewerSource( ...
                fullfile(testCase.Folder, 'green.dat'), iSourceID(pm));

            testCase.verifyEqual(action, "keep");
            testCase.verifyEqual(reason, "sameFile");
        end

        function swapsToAnUnrelatedShownFile(testCase)
            [pm, ~] = testCase.firstRun({'pmSkipStepA'});
            other = fullfile(testCase.Folder, 'other.dat');
            saveData(other, testCase.Base * 5, 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);

            [action, reason] = pm.chooseDataViewerSource(other, iSourceID(pm));

            testCase.verifyEqual(action, "swap");
            testCase.verifyEqual(reason, "currentNotDerived");
        end

        function swapsWhenThereIsNoStoredSource(testCase)
            [pm, shown] = testCase.firstRun({'pmSkipStepA'});
            streamID = pm.nodes(strcmp({pm.nodes.kind}, 'stream')).id;

            [a1, r1] = pm.chooseDataViewerSource(shown, NaN);
            [a2, r2] = pm.chooseDataViewerSource(shown, []);
            [a3, r3] = pm.chooseDataViewerSource(shown, 9999);
            [a4, r4] = pm.chooseDataViewerSource(shown, streamID);

            testCase.verifyEqual([a1 a2 a3 a4], repmat("swap", 1, 4));
            testCase.verifyEqual([r1 r2 r3 r4], repmat("noStoredSource", 1, 4));
        end

        function swapsWhenTheOriginalFileIsGoneAndTheResultStillRuns(testCase)
            [pm, shown] = testCase.firstRun({'pmSkipStepA', 'pmSkipStepB'});
            delete(fullfile(testCase.Folder, 'green.dat'));

            [action, reason] = testCase.applyChoice(pm, shown);

            testCase.verifyEqual(action, "swap");
            testCase.verifyEqual(reason, "storedSourceMissing");
            result = testCase.runPipeline(pm, true);
            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'The steps are already in the shown file.');
            testCase.verifyEqual(iCalls('pmSkipStepB'), 1);
            testCase.verifyEqual(testCase.outputOf(result), testCase.Base + 2);
        end

        function editedSettingsAreRecomputedFromTheOriginalOnEveryRun(testCase)
            % Edit a step before the shown temporary file on a second, third, ... run.
            for reuse = [true false]
                setappdata(0, 'pmSkipLog', struct());
                [pm, shown] = testCase.firstRun({'pmReuseGainStep'}, 2);
                testCase.verifyEqual(iFactor(testCase.fileData(shown), testCase.Base), 2);

                for gain = [3 4 5]
                    [action, ~] = testCase.applyChoice(pm, shown);
                    testCase.verifyEqual(action, "keep");
                    pm = testCase.setGain(pm, gain);

                    result = testCase.runPipeline(pm, reuse);

                    testCase.verifyEqual(iFactor(testCase.outputOf(result), testCase.Base), gain, ...
                        sprintf('Gain %d, reuse %d: the original must be filtered once.', gain, reuse));
                    shown = testCase.outputPath(result);
                end
            end
        end

        function settingsEditedWithoutRunningAreNotAppliedOnTopOfTheResult(testCase)
            % The tool reopens after the edit but before any run: the shown file no longer
            % matches the pipeline, yet it was made from the same original.
            [pm, shown] = testCase.firstRun({'pmReuseGainStep'}, 2);
            pm = testCase.setGain(pm, 3);

            [action, reason] = testCase.applyChoice(pm, shown);
            result = testCase.runPipeline(pm, true);

            testCase.verifyEqual(action, "keep");
            testCase.verifyEqual(reason, "currentIsResultOfStoredSource");
            testCase.verifyEqual(iFactor(testCase.outputOf(result), testCase.Base), 3);
        end

        function aStepRemovedBeforeReopeningIsRecomputedFromTheOriginal(testCase)
            [~, shown] = testCase.firstRun({'pmSkipStepA', 'pmSkipStepB'});

            % A new pipeline over the same folder holds only step B.
            pmB = PipelineManager(testCase.Folder, '', testCase.ProjectRoot);
            pmB.addStep('pmSkipStepB', 'input', 'green.dat');
            setappdata(0, 'pmSkipLog', struct());

            [action, reason] = testCase.applyChoice(pmB, shown);
            result = testCase.runPipeline(pmB, true);

            testCase.verifyEqual(action, "keep");
            testCase.verifyEqual(reason, "currentIsResultOfStoredSource");
            testCase.verifyEqual(testCase.outputOf(result), testCase.Base + 1, ...
                'Step B alone applies to the original, not on top of A and B.');
        end

        function aStepAddedOnTopRunsOnlyTheNewStepWhenReusing(testCase)
            [pm, shown] = testCase.firstRun({'pmSkipStepA', 'pmSkipStepB'});
            [action, ~] = testCase.applyChoice(pm, shown);
            pm.addStep('pmReuseGainStep', 'input', 'pmSkipStepB_1:outData');
            pm = testCase.setGain(pm, 2);
            setappdata(0, 'pmSkipLog', struct());

            result = testCase.runPipeline(pm, true);

            testCase.verifyEqual(action, "keep");
            testCase.verifyEqual(iCalls('pmSkipStepA'), 0, 'The finished steps come from the shown file.');
            testCase.verifyEqual(iCalls('pmSkipStepB'), 0);
            testCase.verifyEqual(iCalls('pmReuseGainStep'), 1);
            testCase.verifyEqual(testCase.outputOf(result), (testCase.Base + 2) * 2);
        end

        function aStepAddedOnTopRunsEverythingFromTheOriginalWhenNotReusing(testCase)
            [pm, shown] = testCase.firstRun({'pmSkipStepA', 'pmSkipStepB'});
            testCase.applyChoice(pm, shown);
            pm.addStep('pmReuseGainStep', 'input', 'pmSkipStepB_1:outData');
            pm = testCase.setGain(pm, 2);
            setappdata(0, 'pmSkipLog', struct());

            result = testCase.runPipeline(pm, false);

            testCase.verifyEqual(iCalls('pmSkipStepA'), 1);
            testCase.verifyEqual(iCalls('pmSkipStepB'), 1);
            testCase.verifyEqual(testCase.outputOf(result), (testCase.Base + 2) * 2, ...
                'A and B must be applied once, to the original.');
        end

        function aPipelineRootedAtATemporaryResultKeepsItsSource(testCase)
            % The user cleared the pipeline and continued from the shown result T1.
            [~, t1] = testCase.firstRun({'pmSkipStepA'});
            pmChain = PipelineManager(testCase.Folder, '', testCase.ProjectRoot);
            [~, t1Name, t1Ext] = fileparts(t1);
            pmChain.addStep('pmSkipStepB', 'input', [t1Name t1Ext]);
            setappdata(0, 'pmSkipLog', struct());
            result = testCase.runPipeline(pmChain, true);
            t2 = testCase.outputPath(result);
            testCase.verifyEqual(testCase.fileData(t2), testCase.Base + 2);

            [action, reason] = pmChain.chooseDataViewerSource(t2, iSourceID(pmChain));
            result2 = testCase.runPipeline(pmChain, false);

            testCase.verifyEqual(action, "keep");
            testCase.verifyEqual(reason, "currentIsResultOfStoredSource");
            testCase.verifyEqual(testCase.outputOf(result2), testCase.Base + 2, ...
                'Re-running B from T1 must not add a third step.');
        end
    end

    methods
        function [pm, shown] = firstRun(testCase, funcs, gain)
            % Run the recipe once on green.dat under the viewer's temporary-output policy.
            pm = PipelineManager(testCase.Folder, '', testCase.ProjectRoot);
            prev = 'green.dat';
            for k = 1:numel(funcs)
                tag = pm.addStep(funcs{k}, 'input', prev);
                prev = sprintf('%s:outData', tag);
            end
            if nargin > 2
                pm = testCase.setGain(pm, gain);
            end
            setappdata(0, 'pmSkipLog', struct());
            result = pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);
            shown = testCase.outputPath(result);
        end

        function [action, reason] = applyChoice(~, pm, shown)
            % What PipelineManagerTool does with the decision.
            [action, reason] = pm.chooseDataViewerSource(shown, iSourceID(pm));
            if action == "swap"
                iSwapSource(pm, shown);
            end
        end

        function result = runPipeline(~, pm, reuse)
            pm.b_skipSteps = reuse;
            result = pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);
        end

        function p = outputPath(testCase, result)
            out = result.viewerCompatibleFinalOutputs;
            p = fullfile(testCase.Folder, char(out.ActualFileName(end)));
        end

        function data = outputOf(testCase, result)
            data = testCase.fileData(testCase.outputPath(result));
        end

        function data = fileData(~, filePath)
            data = single(loadData(filePath));
        end

        function pm = setGain(testCase, pm, gain)
            % Parameters are edited through a dialog in the app: tests apply them through
            % the supported save/load .pipe representation.
            pipeFile = [tempname(testCase.Folder) '.pipe'];
            pm.savePipe(pipeFile);
            cleanupFile = onCleanup(@() iDeleteIfFile(pipeFile));
            loaded = load(pipeFile, 'pipeStruct', '-mat');
            pipeStruct = loaded.pipeStruct;
            for iNode = 1:numel(pipeStruct.nodes)
                node = pipeStruct.nodes(iNode);
                if strcmpi(node.kind, 'stream') && strcmpi(node.info.name, 'pmReuseGainStep')
                    params = node.info.parameters;
                    idx = find(strcmpi({params.name}, 'Gain'), 1);
                    params(idx).value = gain;
                    pipeStruct.nodes(iNode).info.parameters = params;
                end
            end
            save(pipeFile, 'pipeStruct', '-mat');
            pm.loadPipe(pipeFile);
            clear cleanupFile
        end
    end
end

function id = iSourceID(pm)
id = pm.nodes(find(strcmp({pm.nodes.kind}, 'folder'), 1, 'first')).id;
end

function iSwapSource(pm, shown)
% PipelineManagerTool.ensureDataViewerPresetFileSource with the public APIs.
oldIdx = find(strcmp({pm.nodes.kind}, 'folder'), 1, 'first');
oldID = pm.nodes(oldIdx).id;
oldTag = pm.nodes(oldIdx).name;
[newTag, ~] = pm.ensureFileSourceStep(shown);
conns = pm.connections([pm.connections.sourceNodeID] == oldID);
for c = conns(:)'
    target = pm.nodes([pm.nodes.id] == c.targetNodeID);
    pm.reconnectStep(target.name, 'input', sprintf('%s:data', newTag));
end
pm.rmStep(oldTag);
end

function f = iFactor(data, base)
q = data ./ base;
if all(abs(q(:) - q(1)) < 1e-4)
    f = double(q(1));
else
    f = NaN;
end
end

function n = iCalls(name)
callLog = getappdata(0, 'pmSkipLog');
if isstruct(callLog) && isfield(callLog, name)
    n = callLog.(name);
else
    n = 0;
end
end

function iDeleteIfFile(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function parserFolder = iCreateFixtureCategory(projectRoot)
% PipelineManager discovers functions under Analysis/: copy the fixtures into a
% temporary category there (createAnalysisFixtureCategory).
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
names = {'pmSkipStepA', 'pmSkipStepB', 'pmReuseGainStep'};
files = struct('name', names, 'source', fullfile(fixtureFolder, strcat(names, '.m')));
parserFolder = createAnalysisFixtureCategory(projectRoot, files);
end

function iRemoveFixtureCategory(parserFolder, projectRoot)
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(parserFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerDataViewerSourceRule:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmSkipStepA pmSkipStepB pmReuseGainStep
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
end
