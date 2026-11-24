classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerReuseExistingFiles < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERREUSEEXISTINGFILES Reuse of existing outputs and in-place replacement.
    %
    %   Reuse on (pm.b_skipSteps = true): an identical re-run is satisfied by the
    %   outputs on disk and writes nothing; a changed step is not reused.
    %   Reuse off: the output the pipeline itself wrote earlier (listed in
    %   dataHistory.mat) is replaced in place, never copied to name_1, name_2...
    %   A same-named file the pipeline never wrote, or one it reads as an input,
    %   keeps its content and the new output gets a numbered name.
    %
    %   pmSkipStepA / pmSkipStepB / pmReuseGainStep count their real executions in
    %   getappdata(0, 'pmSkipLog'). The data is 2 x 3 x 4 (Y, X, T): non-square.

    properties
        ProjectRoot char
        SaveFolder char
        SourceFile char
        SourceData single
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
            testCase.SourceData = single(reshape(1:24, [2 3 4]));
            saveData(testCase.SourceFile, testCase.SourceData, ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));

            setappdata(0, 'pmSkipLog', struct());
            testCase.addTeardown(@() setappdata(0, 'pmSkipLog', struct()));
        end
    end

    methods (Test)
        function identicalRerunWithReuseOnLeavesOutputsUntouched(testCase)
            pm = testCase.newPipeline();
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);

            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            before = iDatInfo(testCase.SaveFolder);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            after = iDatInfo(testCase.SaveFolder);

            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'The second run must reuse the output.');
            testCase.verifyEqual({after.name}, {before.name}, 'No new file may appear.');
            testCase.verifyEqual([after.datenum], [before.datenum], 'No file may be rewritten.');
            testCase.verifyEqual([after.bytes], [before.bytes]);
        end

        function changedParameterWithReuseOnIsNotReusedAndReplacesInPlace(testCase)
            pm = testCase.newGainPipeline(1);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            testCase.verifyEqual(iCalls('pmReuseGainStep'), 1);

            pm2 = testCase.newGainPipeline(3);
            pm2.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(iCalls('pmReuseGainStep'), 2, ...
                'A changed parameter must not be satisfied by the old output.');
            testCase.verifyEqual(iDatNames(testCase.SaveFolder), ...
                {'green.dat', 'pmReuseGainStep.dat'}, ...
                'The new result replaces the old output in place.');
            result = loadData(fullfile(testCase.SaveFolder, 'pmReuseGainStep.dat'));
            testCase.verifyEqual(single(result), testCase.SourceData * 3, ...
                'The file must hold the result of the new parameter.');
            testCase.verifyEqual(iHistoryFileCount(testCase.SaveFolder, 'pmReuseGainStep.dat'), 1);
        end

        function reuseOffReplacesOutputInPlaceEveryRun(testCase)
            pm = testCase.newPipeline();
            pm.b_skipSteps = false;
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);

            for k = 1:3
                pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            end

            testCase.verifyEqual(iCalls('pmSkipStepA'), 3);
            testCase.verifyEqual(iDatNames(testCase.SaveFolder), {'green.dat', 'pmSkipStepA.dat'}, ...
                'No name_1, name_2 copies may accumulate.');
            testCase.verifyEqual(iHistoryFileCount(testCase.SaveFolder, 'pmSkipStepA.dat'), 1);
            testCase.verifyEmpty(iTempFiles(testCase.SaveFolder), ...
                'No scratch file may remain after replacing an output.');
        end

        function reuseOffInRamSafeModeReplacesOutputInPlace(testCase)
            pm = testCase.newPipeline();
            pm.b_skipSteps = false;
            pm.ramMode = 'ramsafe';
            tagA = pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            pm.addStep('pmSkipStepB', 'input', sprintf('%s:outData', tagA));

            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(iDatNames(testCase.SaveFolder), {'green.dat', 'pmSkipStepB.dat'});
            result = loadData(fullfile(testCase.SaveFolder, 'pmSkipStepB.dat'));
            testCase.verifyEqual(single(result), testCase.SourceData + 2);
            testCase.verifyEmpty(iTempFiles(testCase.SaveFolder));
        end

        function unregisteredSameNamedFileIsNeverOverwritten(testCase)
            foreign = fullfile(testCase.SaveFolder, 'pmSkipStepA.dat');
            foreignData = testCase.SourceData * 100;
            saveData(foreign, foreignData, 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);

            pm = testCase.newPipeline();
            pm.b_skipSteps = false;
            pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(iDatNames(testCase.SaveFolder), ...
                {'green.dat', 'pmSkipStepA.dat', 'pmSkipStepA_1.dat'});
            testCase.verifyEqual(single(loadData(foreign)), foreignData, ...
                'A file the pipeline never wrote must keep its content.');
            testCase.verifyEqual(single(loadData(fullfile(testCase.SaveFolder, 'pmSkipStepA_1.dat'))), ...
                testCase.SourceData + 1);
        end

        function pipelineInputIsNeverReplacedByItsOwnOutput(testCase)
            seed = testCase.newPipeline();
            seed.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            seed.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            inputFile = fullfile(testCase.SaveFolder, 'pmSkipStepA.dat');
            inputBefore = single(loadData(inputFile));

            pm = testCase.newPipeline();
            pm.b_skipSteps = false;
            pm.addStep('pmSkipStepA', 'input', inputFile);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(single(loadData(inputFile)), inputBefore, ...
                'The file the step reads must not be overwritten while it is read.');
            testCase.verifyEqual(single(loadData(fullfile(testCase.SaveFolder, 'pmSkipStepA_1.dat'))), ...
                inputBefore + 1);
        end

        function copiedFolderReusesWhenInputHistoryRecordsAnotherFolder(testCase)
            % The input file carries the history of the folder it was made in. A
            % copy of that folder keeps the history, whose folder paths no longer
            % equal the current folder: they must not stop the reuse.
            seed = testCase.newPipeline();
            seed.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            seed.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            iRelocateHistoryFolders(testCase.SaveFolder, 'pmSkipStepA.dat', ...
                fullfile(tempdir, 'UMITElsewhere'));
            inputFile = fullfile(testCase.SaveFolder, 'pmSkipStepA.dat');

            pm = testCase.newPipeline();
            pm.addStep('pmSkipStepB', 'input', inputFile);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);
            namesFirst = iDatNames(testCase.SaveFolder);
            pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(iCalls('pmSkipStepB'), 1, 'The second run must reuse the output.');
            testCase.verifyEqual(iDatNames(testCase.SaveFolder), namesFirst);
        end

        function cancelledRunLeavesNoIntermediateFile(testCase)
            pm = testCase.newPipeline();
            pm.ramMode = 'ramsafe';
            tagA = pm.addStep('pmSkipStepA', 'input', testCase.SourceFile);
            pm.addStep('pmSkipStepB', 'input', sprintf('%s:outData', tagA));

            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', ...
                'PrintSummary', false, ...
                'CancelFcn', @() iCalls('pmSkipStepA') >= 1);

            testCase.verifyEqual(iCalls('pmSkipStepA'), 1, 'The run must stop after the first step.');
            testCase.verifyEqual(iCalls('pmSkipStepB'), 0);
            testCase.verifyNotEqual(result.status, "completed");
            testCase.verifyEmpty(iTempFiles(testCase.SaveFolder), ...
                'A cancelled run must not leave its intermediate on disk.');
            testCase.verifyEqual(iDatNames(testCase.SaveFolder), {'green.dat'});
        end
    end

    methods
        function pm = newPipeline(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
        end

        function pm = newGainPipeline(testCase, gain)
            pm = testCase.newPipeline();
            pm.addStep('pmReuseGainStep', 'input', testCase.SourceFile);
            pm = testCase.setGain(pm, 'pmReuseGainStep', gain);
        end

        function pm = setGain(testCase, pm, funcName, gain)
            % Parameters are edited through a dialog in the app: tests apply them
            % through the supported save/load .pipe representation.
            pipeFile = [tempname(testCase.SaveFolder) '.pipe'];
            pm.savePipe(pipeFile);
            cleanupFile = onCleanup(@() iDeleteIfFile(pipeFile));
            loaded = load(pipeFile, 'pipeStruct', '-mat');
            pipeStruct = loaded.pipeStruct;
            for iNode = 1:numel(pipeStruct.nodes)
                node = pipeStruct.nodes(iNode);
                if strcmpi(node.kind, 'stream') && strcmpi(node.info.name, funcName)
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

function n = iCalls(name)
callLog = getappdata(0, 'pmSkipLog');
if isstruct(callLog) && isfield(callLog, name)
    n = callLog.(name);
else
    n = 0;
end
end

function info = iDatInfo(folder)
info = dir(fullfile(folder, '*.dat'));
[~, order] = sort({info.name});
info = info(order);
end

function names = iDatNames(folder)
names = {iDatInfo(folder).name};
end

function names = iTempFiles(folder)
listing = dir(fullfile(folder, 'PMTMP_*'));
names = sort({listing.name});
end

function n = iHistoryFileCount(folder, fileName)
stored = load(fullfile(folder, 'dataHistory.mat'), 'dataHistory');
n = nnz(strcmpi(string({stored.dataHistory.file}), fileName));
end

function iRelocateHistoryFolders(folder, fileName, otherFolder)
% Rewrite the SaveFolder/RawFolder of every recorded step of one file's history
% to another folder, as if the file had been made elsewhere and then copied.
historyFile = fullfile(folder, 'dataHistory.mat');
stored = load(historyFile, 'dataHistory');
dataHistory = stored.dataHistory;
idx = find(strcmpi(string({dataHistory.file}), fileName), 1);
for k = 1:numel(dataHistory(idx).info)
    inputs = dataHistory(idx).info(k).inputs;
    for j = 1:numel(inputs)
        if any(strcmpi(inputs(j).name, {'SaveFolder', 'RawFolder'}))
            inputs(j).source = otherFolder;
        end
    end
    dataHistory(idx).info(k).inputs = inputs;
end
save(historyFile, 'dataHistory');
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
    error('Umitoolbox:TestPipelineManagerReuseExistingFiles:UnsafeCleanupPath', ...
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
