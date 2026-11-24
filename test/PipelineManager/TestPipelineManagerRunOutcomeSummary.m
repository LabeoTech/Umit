classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerRunOutcomeSummary < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERRUNOUTCOMESUMMARY One set of folder outcomes for every run report.
    %
    %   Every selected folder ends in exactly one category (completed, completed
    %   with errors, failed, skipped, cancelled, not run), so the counts of
    %   result.folderSummary add up to the number of selected folders, whatever
    %   mix of outcomes the run had. "Created files" counts the outputs that exist
    %   after the run (deleted intermediates are not counted; temporary files are
    %   counted apart). A warning raised by a step is kept with its step and
    %   folder and never changes the run status.
    %
    %   pmSkipStepA counts its calls and pmOutcomeStep can warn or fail on demand.
    %   The data is 2 x 3 x 4 (Y, X, T): non-square.

    properties
        ProjectRoot char
        Root char
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
            testCase.Root = fixture.Folder;

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));

            setappdata(0, 'pmSkipLog', struct());
            testCase.addTeardown(@() setappdata(0, 'pmSkipLog', struct()));
        end
    end

    methods (Test)
        function mixedBatchCountsAddUpToTheSelectedFolders(testCase)
            % Folders 1 and 2: step 1 runs, step 2 fails (completed with errors).
            % Folder 3 has no input file: nothing runs (failed).
            folders = testCase.makeFolders(3, [true true false]);
            pm = PipelineManager(folders, repmat({''}, 1, 3), testCase.ProjectRoot);
            tagA = pm.addStep('pmSkipStepA', 'input', 'green.dat');
            pm.addStep('pmOutcomeStep', 'input', sprintf('%s:outData', tagA));
            pm = testCase.setMode(pm, 'pmOutcomeStep', 'fail', folders{1});

            result = testCase.runQuietly(pm);

            s = result.folderSummary;
            testCase.verifyEqual(result.status, "completed_with_errors");
            testCase.verifyEqual(s.nSelected, 3);
            testCase.verifyEqual(s.nCompleted, 0);
            testCase.verifyEqual(s.nCompletedWithErrors, 2, ...
                'A folder whose first step ran and second failed must be counted.');
            testCase.verifyEqual(s.nFailed, 1);
            testCase.verifyEqual(iTotal(s), s.nSelected, 'The categories must add up to the selected folders.');

            text = PipelineManager.formatFolderSummary(result);
            testCase.verifySubstring(text, '3 selected');
            testCase.verifySubstring(text, '2 completed with errors');
            testCase.verifySubstring(text, '1 failed');
        end

        function stopRecordsTheFoldersThatNeverStarted(testCase)
            folders = testCase.makeFolders(3, [true true true]);
            pm = PipelineManager(folders, repmat({''}, 1, 3), testCase.ProjectRoot);
            pm.addStep('pmSkipStepA', 'input', 'green.dat');

            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false, ...
                'CancelFcn', @() iCalls('pmSkipStepA') >= 1);

            s = result.folderSummary;
            testCase.verifyEqual(result.status, "cancelled");
            testCase.verifyEqual(height(result.globalPipeLog), 3, ...
                'Every selected folder needs a row, also the ones never reached.');
            testCase.verifyEqual(string(result.globalPipeLog.Status), ["Cancelled"; "Not run"; "Not run"]);
            testCase.verifyEqual(s.nSelected, 3);
            testCase.verifyEqual(s.nCancelled, 1);
            testCase.verifyEqual(s.nNotRun, 2);
            testCase.verifyEqual(iTotal(s), 3);
            testCase.verifySubstring(PipelineManager.formatFolderSummary(result), '2 not run');
            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'Folders 2 and 3 must not run.');
        end

        function stoppedFolderReportsTheFilesItCreated(testCase)
            folders = testCase.makeFolders(2, [true true]);
            pm = PipelineManager(folders, repmat({''}, 1, 2), testCase.ProjectRoot);
            pm.addStep('pmSkipStepA', 'input', 'green.dat');

            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false, ...
                'CancelFcn', @() iCalls('pmSkipStepA') >= 1);

            testCase.verifyTrue(isfile(fullfile(folders{1}, 'pmSkipStepA.dat')));
            testCase.verifyEqual(result.folderSummary.nCreatedFiles, 1, ...
                'The output the cancelled folder wrote before the stop is a created file.');
        end

        function createdFilesCountOnlyOutputsThatExist(testCase)
            folders = testCase.makeFolders(1, true);
            pm = PipelineManager(folders, {''}, testCase.ProjectRoot);
            pm.ramMode = 'ramsafe';
            tagA = pm.addStep('pmSkipStepA', 'input', 'green.dat');
            pm.addStep('pmSkipStepB', 'input', sprintf('%s:outData', tagA));

            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyGreaterThan(height(result.createdFiles), 1, ...
                'The RAM-safe run records its intermediate among the created files.');
            testCase.verifyEqual(result.folderSummary.nCreatedFiles, 1, ...
                'Only the output that exists after the run is counted.');
            testCase.verifyEqual(result.folderSummary.nTemporaryFiles, 0);
            testCase.verifySubstring(PipelineManager.formatFolderSummary(result), 'Created files: 1.');
        end

        function stepWarningIsKeptWithItsStepAndFolder(testCase)
            folders = testCase.makeFolders(1, true);
            pm = PipelineManager(folders, {''}, testCase.ProjectRoot);
            pm.addStep('pmOutcomeStep', 'input', 'green.dat');
            pm = testCase.setMode(pm, 'pmOutcomeStep', 'warn', folders{1});

            lastwarn('earlier warning', 'Test:earlier');
            result = testCase.runQuietly(pm);
            [msg, id] = lastwarn;

            testCase.verifyEqual(result.status, "completed", 'A warning must not change the status.');
            testCase.verifyEqual(height(result.stepWarnings), 1);
            testCase.verifyEqual(result.stepWarnings.Identifier(1), "Umitoolbox:pmOutcomeStep:Test");
            testCase.verifyEqual(result.stepWarnings.FunctionName(1), "pmOutcomeStep");
            testCase.verifySubstring(result.stepWarnings.StepTag(1), 'pmOutcomeStep');
            testCase.verifyEqual(result.stepWarnings.SaveFolder(1), string(folders{1}));
            testCase.verifySubstring(result.stepWarnings.Message(1), 'Test warning');
            testCase.verifyEqual(result.folderSummary.nStepWarnings, 1);
            testCase.verifySubstring(PipelineManager.formatFolderSummary(result), 'Step warnings: 1.');
            testCase.verifyEqual(msg, 'earlier warning', 'The caller''s lastwarn must be restored.');
            testCase.verifyEqual(id, 'Test:earlier');
        end

        function runOutcomeSentencePointsAtStepWarnings(testCase)
            folders = testCase.makeFolders(1, true);
            pm = PipelineManager(folders, {''}, testCase.ProjectRoot);
            pm.addStep('pmOutcomeStep', 'input', 'green.dat');
            pm = testCase.setMode(pm, 'pmOutcomeStep', 'warn', folders{1});

            result = testCase.runQuietly(pm);
            info = PipelineManager.describeExecutionResult(result);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyEqual(info.severity, 'warning', ...
                'A clean run with a step warning must not read as a plain success.');
            testCase.verifySubstring(info.message, 'Pipeline completed successfully.');
            testCase.verifySubstring(info.message, '1 step warning(s)');
        end

        function warningOfAFailingStepIsNotLost(testCase)
            folders = testCase.makeFolders(1, true);
            pm = PipelineManager(folders, {''}, testCase.ProjectRoot);
            pm.addStep('pmOutcomeStep', 'input', 'green.dat');
            pm = testCase.setMode(pm, 'pmOutcomeStep', 'warnfail', folders{1});

            result = testCase.runQuietly(pm);

            testCase.verifyEqual(result.status, "failed");
            testCase.verifyEqual(height(result.stepWarnings), 1);
            testCase.verifyEqual(result.stepWarnings.Identifier(1), "Umitoolbox:pmOutcomeStep:Test");
        end

        function aRunWithoutWarningsHasAnEmptyWarningTable(testCase)
            folders = testCase.makeFolders(1, true);
            pm = PipelineManager(folders, {''}, testCase.ProjectRoot);
            pm.addStep('pmSkipStepA', 'input', 'green.dat');

            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(height(result.stepWarnings), 0);
            testCase.verifyEqual(result.folderSummary.nStepWarnings, 0);
            testCase.verifyEmpty(strfind(PipelineManager.formatFolderSummary(result), 'warning'));
        end

        function countFolderOutcomesAssignsEveryFolderOneCategory(testCase)
            s = PipelineManager.countFolderOutcomes( ...
                ["Completed"; "Partially executed"; "Failed"; "Skipped"; "Cancelled"; "Not run"; "Odd"], 9);

            testCase.verifyEqual([s.nCompleted s.nCompletedWithErrors s.nFailed s.nSkipped ...
                s.nCancelled s.nOther], [1 1 1 1 1 1]);
            testCase.verifyEqual(s.nNotRun, 1 + 2, 'Two folders without a log row were never reached.');
            testCase.verifyEqual(iTotal(s), 9);
        end

        function formatFolderSummaryCountsAHandBuiltResult(testCase)
            result = struct();
            result.globalPipeLog = table({'a'; 'b'; 'c'; 'd'}, ...
                {'Completed'; 'Partially executed'; 'Failed'; 'Skipped'}, ...
                'VariableNames', {'SaveFolder', 'Status'});

            text = PipelineManager.formatFolderSummary(result);

            testCase.verifyEqual(text, ['Folders: 4 selected, 1 completed, 1 completed with errors, ' ...
                '1 failed, 1 skipped, 0 cancelled, 0 not run. Created files: 0.']);
            testCase.verifyEqual(PipelineManager.formatFolderSummary(struct()), ...
                'Folders: no run information.');
        end
    end

    methods
        function folders = makeFolders(testCase, n, hasInput)
            % n folders with AcqInfos and, where hasInput is true, a green.dat.
            folders = cell(1, n);
            for k = 1:n
                folders{k} = fullfile(testCase.Root, sprintf('data%d', k));
                mkdir(folders{k});
                AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                    'FrameRateHz', 10, 'ExposureMsec', 1);
                AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                    struct('DatFile', 'green.dat', 'Length', 4, 'FrameRateHz', 10));
                save(fullfile(folders{k}, 'AcqInfos.mat'), 'AcqInfoStream');
                if hasInput(k)
                    saveData(fullfile(folders{k}, 'green.dat'), single(reshape(1:24, [2 3 4])), ...
                        'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);
                end
            end
        end

        function pm = setMode(~, pm, funcName, mode, scratchFolder)
            % Parameters are edited through a dialog in the app: tests apply them
            % through the supported save/load .pipe representation.
            pipeFile = [tempname(scratchFolder) '.pipe'];
            pm.savePipe(pipeFile);
            cleanupFile = onCleanup(@() iDeleteIfFile(pipeFile));
            loaded = load(pipeFile, 'pipeStruct', '-mat');
            pipeStruct = loaded.pipeStruct;
            for iNode = 1:numel(pipeStruct.nodes)
                node = pipeStruct.nodes(iNode);
                if strcmpi(node.kind, 'stream') && strcmpi(node.info.name, funcName)
                    params = node.info.parameters;
                    idx = find(strcmpi({params.name}, 'Mode'), 1);
                    params(idx).value = mode;
                    pipeStruct.nodes(iNode).info.parameters = params;
                end
            end
            save(pipeFile, 'pipeStruct', '-mat');
            pm.loadPipe(pipeFile);
            clear cleanupFile
        end

        function result = runQuietly(~, pm)
            % The step warnings and failures are expected; keep the log clean.
            warnState = warning('off', 'all');
            restore = onCleanup(@() warning(warnState));
            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            clear restore
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

function n = iTotal(s)
n = s.nCompleted + s.nCompletedWithErrors + s.nFailed + s.nSkipped + s.nCancelled + s.nNotRun + s.nOther;
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
names = {'pmSkipStepA', 'pmSkipStepB', 'pmOutcomeStep'};
files = struct('name', names, 'source', fullfile(fixtureFolder, strcat(names, '.m')));
parserFolder = createAnalysisFixtureCategory(projectRoot, files);
end

function iRemoveFixtureCategory(parserFolder, projectRoot)
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(parserFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerRunOutcomeSummary:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmSkipStepA pmSkipStepB pmOutcomeStep
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
end
