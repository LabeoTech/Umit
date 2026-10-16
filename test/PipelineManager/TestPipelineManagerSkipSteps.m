classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerSkipSteps < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERSKIPSTEPS History-based step skipping (b_skipSteps).
    %
    %   Re-running an unchanged pipeline must not execute a step whose result is
    %   already on disk and recorded in dataHistory.mat, and the run result must
    %   still report that existing file as the leaf output, so DataViewer can show
    %   it again. This holds for permanent outputs and for the PMTMP_* files of the
    %   'viewerTemp' policy.
    %
    %   The fixtures pmSkipStepA / pmSkipStepB add 1 to the data and count their
    %   real executions in getappdata(0, 'pmSkipLog'), so a skipped step shows up
    %   as a counter that does not move.

    properties
        ProjectRoot char
        SaveFolder char
        SourceFile char
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
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);

            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                'FrameRateHz', 10, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                struct('DatFile', 'green.dat', 'Length', 4, 'FrameRateHz', 10));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.SourceFile = fullfile(testCase.SaveFolder, 'green.dat');
            saveData(testCase.SourceFile, single(reshape(1:24, [2 3 4])), ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));

            setappdata(0, 'pmSkipLog', struct());
            testCase.addTeardown(@() setappdata(0, 'pmSkipLog', struct()));
        end
    end

    methods (Test)
        function unchangedRerunSkipsPermanentLeaf(testCase)
            pm = testCase.newPipeline();
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);

            r1 = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            filesAfterFirst = iDatFiles(testCase.SaveFolder);
            r2 = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(r1.status, "completed");
            testCase.verifyEqual(r2.status, "completed");
            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'The second run must skip the step.');
            testCase.verifyEqual(iDatFiles(testCase.SaveFolder), filesAfterFirst, ...
                'A skipped run must not create another output file.');

            testCase.verifyEqual(height(r2.viewerCompatibleFinalOutputs), 1);
            testCase.verifyEqual(r2.viewerCompatibleFinalOutputs.ActualFileName(1), ...
                r1.viewerCompatibleFinalOutputs.ActualFileName(1));
            testCase.verifyTrue(r2.viewerCompatibleFinalOutputs.FileExists(1));
        end

        function unchangedRerunSkipsViewerTempLeafAndReportsExistingFile(testCase)
            pm = testCase.newPipeline();
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);

            r1 = pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);
            tempAfterFirst = iTempFiles(testCase.SaveFolder);
            r2 = pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);

            testCase.verifyEqual(numel(tempAfterFirst), 1);
            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'The second run must skip the step.');
            testCase.verifyEqual(iTempFiles(testCase.SaveFolder), tempAfterFirst, ...
                'A skipped run must not create another PMTMP_* file.');

            % The skipped run re-presents the existing temporary file.
            testCase.verifyEqual(height(r2.viewerCompatibleFinalOutputs), 1);
            row = r2.viewerCompatibleFinalOutputs(1, :);
            testCase.verifyEqual(char(row.ActualFileName), tempAfterFirst{1});
            testCase.verifyTrue(row.FileExists);
            testCase.verifyTrue(row.IsTemporary);
            testCase.verifyEqual(r2.viewerCompatibleFinalOutputs.ActualFileName(1), ...
                r1.viewerCompatibleFinalOutputs.ActualFileName(1));
        end

        function skippingDisabledAlwaysRunsAndWritesAgain(testCase)
            pm = testCase.newPipeline();
            pm.b_skipSteps = false;
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);

            pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);
            pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);

            testCase.verifyEqual(iCalls('pmSkipStepA'), 2);
            testCase.verifyEqual(numel(iTempFiles(testCase.SaveFolder)), 2);
        end

        function addingAStepReusesTheCachedPrefix(testCase)
            pm = testCase.newPipeline();
            tagA = pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            testCase.verifyEqual(iCalls('pmSkipStepA'), 1);

            pm.addStep('pmSkipStepB', 'input', sprintf('%s:outData', tagA));
            r2 = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(r2.status, "completed");
            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'Step A is already cached.');
            testCase.verifyEqual(iCalls('pmSkipStepB'), 1, 'Only the new step runs.');

            outputs = r2.viewerCompatibleFinalOutputs;
            testCase.verifyEqual(height(outputs), 1);
            result = loadData(fullfile(testCase.SaveFolder, char(outputs.ActualFileName(1))));
            expected = single(reshape(1:24, [2 3 4])) + 2;
            testCase.verifyEqual(single(result), expected, ...
                'B must have run on the cached output of A.');
        end
    end

    methods (Test)
        function producedOutputIsRecognisedAsCachedResult(testCase)
            pm = testCase.newPipeline();
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);
            tempNames = iTempFiles(testCase.SaveFolder);
            tempFile = fullfile(testCase.SaveFolder, tempNames{1});

            testCase.verifyTrue(pm.isCachedResultOfPipeline(tempFile), ...
                'The PMTMP_* output of the pipeline is a cached result of it.');
            testCase.verifyFalse(pm.isCachedResultOfPipeline(testCase.SourceFile), ...
                'The raw source has no history entry.');
            testCase.verifyFalse(pm.isCachedResultOfPipeline( ...
                fullfile(testCase.SaveFolder, 'does_not_exist.dat')));
            testCase.verifyFalse(pm.isCachedResultOfPipeline(''));
        end

        function cachedResultStaysRecognisedWhenAStepIsAdded(testCase)
            pm = testCase.newPipeline();
            tagA = pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            r1 = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            outputA = fullfile(testCase.SaveFolder, char(r1.viewerCompatibleFinalOutputs.ActualFileName(1)));
            testCase.verifyTrue(pm.isCachedResultOfPipeline(outputA));

            pm.addStep('pmSkipStepB', 'input', sprintf('%s:outData', tagA));

            testCase.verifyTrue(pm.isCachedResultOfPipeline(outputA), ...
                'The output of step A is still a reusable prefix of the longer pipeline.');
        end

        function outputOfADifferentPipelineIsNotACachedResult(testCase)
            first = testCase.newPipeline();
            first.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            r1 = first.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            outputA = fullfile(testCase.SaveFolder, char(r1.viewerCompatibleFinalOutputs.ActualFileName(1)));

            other = testCase.newPipeline();
            other.addStep('pmSkipStepB', 'input', testCase.SourceFile);

            testCase.verifyFalse(other.isCachedResultOfPipeline(outputA), ...
                'A file produced by a different step is not reusable by this pipeline.');
        end
    end

    methods
        function pm = newPipeline(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
        end
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

function names = iDatFiles(folder)
listing = dir(fullfile(folder, '*.dat'));
names = sort({listing.name});
end

function names = iTempFiles(folder)
listing = dir(fullfile(folder, 'PMTMP_*'));
names = sort({listing.name});
end

function parserFolder = iCreateFixtureCategory(projectRoot)
% PipelineManager discovers functions under Analysis/: copy the fixtures into a
% temporary category there (createAnalysisFixtureCategory).
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
names = {'pmSkipStepA', 'pmSkipStepB'};
files = struct('name', names, 'source', fullfile(fixtureFolder, strcat(names, '.m')));
parserFolder = createAnalysisFixtureCategory(projectRoot, files);
end

function iRemoveFixtureCategory(parserFolder, projectRoot)
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(parserFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerSkipSteps:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmSkipStepA pmSkipStepB
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
end
