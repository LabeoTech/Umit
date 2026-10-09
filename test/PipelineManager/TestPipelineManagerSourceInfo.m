classdef TestPipelineManagerSourceInfo < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERSOURCEINFO PipelineManager propagates the source Info.
    %
    %   .dat header Phase 4c-1: .dat files that PipelineManager saves from
    %   in-memory step outputs (and its temporary spills) carry the frame
    %   rate and exposure of the file the branch was read from, unless a
    %   step updated them through a "metaData" output. The same holds for
    %   scripts produced by generateScript. input.dat is headered at 25 Hz
    %   with a 3 ms exposure, while AcqInfos.mat says 10 Hz, so each source
    %   is recognizable.

    properties
        ProjectRoot char
        SaveFolder char
        InputData single
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createSaveFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);
            testCase.InputData = single(reshape(1:24, [2 3 4]));

            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                'FrameRateHz', 10, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                'DatFile', 'input.dat', 'Length', 4, 'FrameRateHz', 10));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            saveData(fullfile(testCase.SaveFolder, 'input.dat'), testCase.InputData, ...
                'Info', struct('frameRateHz', 25, 'exposureMsec', 3), 'DimNames', {'Y', 'X', 'T'});

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));
        end
    end

    methods (Test)
        function ramChainInheritsFileInfo(testCase)
            pm = testCase.newManager();
            firstTag = pm.addStep('pmSourceInfoScale', 'input', 'input.dat', 'saveas', 'scaled.dat');
            pm.addStep('pmSourceInfoScaleB', 'input', firstTag, 'saveas', 'scaled2.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            testCase.verifyHeader('scaled.dat', 25, 3, 2 .* testCase.InputData);
            testCase.verifyHeader('scaled2.dat', 25, 3, 4 .* testCase.InputData);
        end

        function ramSafeSpillsInheritFileInfo(testCase)
            pm = testCase.newManager();
            pm.ramMode = 'ramsafe';
            firstTag = pm.addStep('pmSourceInfoScale', 'input', 'input.dat');
            pm.addStep('pmSourceInfoScaleB', 'input', firstTag, 'saveas', 'spilled.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            testCase.verifyHeader('spilled.dat', 25, 3, 4 .* testCase.InputData);
        end

        function metaDataOutputUpdatesInfo(testCase)
            pm = testCase.newManager();
            firstTag = pm.addStep('pmSourceInfoRate', 'input', 'input.dat', 'saveas', 'rated.dat');
            pm.addStep('pmSourceInfoScale', 'input', firstTag, 'saveas', 'after.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            testCase.verifyHeader('rated.dat', 5, 9, testCase.InputData + 1);
            testCase.verifyHeader('after.dat', 5, 9, 2 .* (testCase.InputData + 1));
        end

        function nonStructMetaDataIsIgnored(testCase)
            pm = testCase.newManager();
            pm.addStep('pmSourceInfoBadMeta', 'input', 'input.dat', 'saveas', 'bad.dat');

            result = testCase.verifyWarningFree(@() pm.executePipeline('PrintSummary', false));

            testCase.assertEqual(result.status, "completed");
            testCase.verifyHeader('bad.dat', 25, 3, testCase.InputData);
        end

        function noSourceFileRequiresFrameRate(testCase)
            % .dat header Phase 6a: without a source file, Y-X-T data has no
            % frame rate; saveData no longer falls back to AcqInfos.mat
            % (10 Hz here), so the save fails and nothing is written.
            pm = testCase.newManager();
            pm.addStep('pmSourceInfoMake', 'saveas', 'made.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyNotEqual(result.status, "completed");
            testCase.verifyFalse(isfile(fullfile(testCase.SaveFolder, 'made.dat')));
        end

        % ------------------------------------------------ axes (Phase 6a)
        function imageOutputIsSavedAsYX(testCase)
            pm = testCase.newManager();
            pm.ramMode = 'ramsafe';
            pm.addStep('pmSourceInfoMean', 'input', 'input.dat', 'saveas', 'mean.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            testCase.verifyAxes('mean.dat', {'Y', 'X'}, [2 3]);
            testCase.verifyTrue(isnan(readDatHeader(fullfile(testCase.SaveFolder, 'mean.dat')).frameRateHz));
            testCase.verifyEqual(loadData(fullfile(testCase.SaveFolder, 'mean.dat')), ...
                mean(testCase.InputData, 3));
        end

        function eventSplitOutputIsSavedAsYXTE(testCase)
            pm = testCase.newManager();
            pm.ramMode = 'ramsafe';
            pm.addStep('pmSourceInfoEvents', 'input', 'input.dat', 'saveas', 'events.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            expected = cat(4, testCase.InputData(:, :, 1:2), testCase.InputData(:, :, 3:4));
            testCase.verifyAxes('events.dat', {'Y', 'X', 'T', 'E'}, [2 3 2 2]);
            testCase.verifyHeader('events.dat', 25, 3, expected);
        end

        function downstreamStepInheritsSavedAxes(testCase)
            % In RAM-safe mode the event-split value goes through a file; the
            % next step's output inherits its Y-X-T-E axes.
            pm = testCase.newManager();
            pm.ramMode = 'ramsafe';
            firstTag = pm.addStep('pmSourceInfoEvents', 'input', 'input.dat');
            pm.addStep('pmSourceInfoScale', 'input', firstTag, 'saveas', 'scaledEvents.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            expected = 2 .* cat(4, testCase.InputData(:, :, 1:2), testCase.InputData(:, :, 3:4));
            testCase.verifyAxes('scaledEvents.dat', {'Y', 'X', 'T', 'E'}, [2 3 2 2]);
            testCase.verifyHeader('scaledEvents.dat', 25, 3, expected);
        end

        function metaDataOutputSetsAxes(testCase)
            pm = testCase.newManager();
            pm.addStep('pmSourceInfoAxes', 'input', 'input.dat', 'saveas', 'relabelled.dat');

            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            testCase.verifyAxes('relabelled.dat', {'Y', 'X', 'E'}, [2 3 4]);
            testCase.verifyEqual(loadData(fullfile(testCase.SaveFolder, 'relabelled.dat')), ...
                testCase.InputData);
        end

        function generatedScriptWritesSameAxes(testCase)
            pm = testCase.newManager();
            pm.addStep('pmSourceInfoMean', 'input', 'input.dat', 'saveas', 'mean.dat');
            pm.addStep('pmSourceInfoEvents', 'input', 'input.dat', 'saveas', 'events.dat');

            testCase.runGeneratedScript(pm);

            testCase.verifyAxes('mean.dat', {'Y', 'X'}, [2 3]);
            testCase.verifyAxes('events.dat', {'Y', 'X', 'T', 'E'}, [2 3 2 2]);
            testCase.verifyHeader('events.dat', 25, 3, ...
                cat(4, testCase.InputData(:, :, 1:2), testCase.InputData(:, :, 3:4)));
        end

        % ------------------------------------------------ parallel branches (Phase 6b-2)
        function parallelBranchesKeepTheirOwnInfo(testCase)
            % Two sources with different rates, exposures, and axes feed
            % parallel branches whose values are held in RAM at the same
            % time; a fork of one branch changes its rate through metaData.
            % Every saved file carries its own branch's Info, in auto and
            % ramsafe modes and from a generated script.
            second = single(reshape(201:236, [2 3 3 2]));
            saveData(fullfile(testCase.SaveFolder, 'second.dat'), second, ...
                'DimNames', {'Y', 'X', 'T', 'E'}, 'Info', struct('frameRateHz', 40, 'exposureMsec', 7));

            for mode = ["auto", "ramsafe", "script"]
                testCase.clearOutputs({'s1', 's2', 'r1', 'after'});
                pm = testCase.newManager();
                firstTag = pm.addStep('pmSourceInfoScale', 'input', 'input.dat', 'saveas', 's1.dat');
                pm.addStep('pmSourceInfoScaleB', 'input', 'second.dat', 'saveas', 's2.dat');
                rateTag = pm.addStep('pmSourceInfoRate', 'input', 'input.dat', 'saveas', 'r1.dat');
                pm.addStep('pmSourceInfoEvents', 'input', rateTag, 'saveas', 'after.dat');
                if mode == "script"
                    testCase.runGeneratedScript(pm);
                else
                    pm.ramMode = char(mode);
                    result = pm.executePipeline('PrintSummary', false);
                    testCase.assertEqual(result.status, "completed", char(mode));
                end
                testCase.verifyHeader('s1.dat', 25, 3, 2 .* testCase.InputData);
                testCase.verifyHeader('s2.dat', 40, 7, 2 .* second);
                testCase.verifyAxes('s2.dat', {'Y', 'X', 'T', 'E'}, [2 3 3 2]);
                testCase.verifyHeader('r1.dat', 5, 9, testCase.InputData + 1);
                rated = testCase.InputData + 1;
                testCase.verifyHeader('after.dat', 5, 9, cat(4, rated(:, :, 1:2), rated(:, :, 3:4)));
                testCase.verifyAxes('after.dat', {'Y', 'X', 'T', 'E'}, [2 3 2 2]);
                testCase.verifyNotEmpty(firstTag);
            end
        end

        % ------------------------------------------------ resolveDimNames
        function resolveDimNamesUsesFittingSourceAxes(testCase)
            info = struct('dimNames', {{'Y', 'X', 'T', 'E'}});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3, 4, 5), info, ''), ...
                {'Y', 'X', 'T', 'E'});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3, 4), info, 'ImageTimeSeries'), ...
                {'Y', 'X', 'T', 'E'}, 'E = 1 still fits');
            yxe = struct('dimNames', {{'Y', 'X', 'E'}});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3, 4), yxe, 'Image'), ...
                {'Y', 'X', 'E'});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3, 4), yxe, ''), ...
                {'Y', 'X', 'E'});
        end

        function resolveDimNamesIgnoresDisagreeingOrUnfitSource(testCase)
            yxt = struct('dimNames', {{'Y', 'X', 'T'}});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3), yxt, 'Image'), ...
                {'Y', 'X'}, 'Image excludes a T axis');
            yx = struct('dimNames', {{'Y', 'X'}});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3), yx, 'ImageTimeSeries'), ...
                {'Y', 'X', 'T'}, 'ImageTimeSeries needs a T axis');
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3, 4, 5), yxt, ''), ...
                {'Y', 'X', 'T', 'E'}, '4-D does not fit Y-X-T');
            bad = struct('dimNames', {{'X', 'Y', 'T'}});
            testCase.verifyEqual(PipelineManager.resolveDimNames(zeros(2, 3, 4), bad, 'Image'), ...
                {'Y', 'X', 'E'}, 'an invalid layout is ignored');
        end

        function resolveDimNamesDefaults(testCase)
            cases = { ...
                [2 3],     'Image',           {'Y', 'X'}; ...
                [2 3 4],   'Image',           {'Y', 'X', 'E'}; ...
                [2 3],     'ImageTimeSeries', {'Y', 'X', 'T'}; ...
                [2 3 4],   'ImageTimeSeries', {'Y', 'X', 'T'}; ...
                [2 3 4 5], 'ImageTimeSeries', {'Y', 'X', 'T', 'E'}; ...
                [2 3],     '',                {'Y', 'X'}; ...
                [2 3 4],   '',                {'Y', 'X', 'T'}; ...
                [2 3 4 5], '',                {'Y', 'X', 'T', 'E'}};
            for k = 1:size(cases, 1)
                names = PipelineManager.resolveDimNames(zeros(cases{k, 1}), [], cases{k, 2});
                testCase.verifyEqual(names, cases{k, 3}, ...
                    sprintf('%d-D, type ''%s''', numel(cases{k, 1}), cases{k, 2}));
            end
        end

        function resolveDimNamesErrors(testCase)
            id = 'Umitoolbox:PipelineManager:unresolvedDimNames';
            testCase.verifyError(@() PipelineManager.resolveDimNames(zeros(2, 3, 4, 5), [], 'Image'), id);
            testCase.verifyError(@() PipelineManager.resolveDimNames(zeros(2, 3, 4, 5, 6), [], ''), id);
        end

        function generatedScriptWritesSameHeaders(testCase)
            pm = testCase.newManager();
            firstTag = pm.addStep('pmSourceInfoRate', 'input', 'input.dat', 'saveas', 'rated.dat');
            pm.addStep('pmSourceInfoScale', 'input', firstTag, 'saveas', 'after.dat');

            testCase.runGeneratedScript(pm);

            testCase.verifyHeader('rated.dat', 5, 9, testCase.InputData + 1);
            testCase.verifyHeader('after.dat', 5, 9, 2 .* (testCase.InputData + 1));
        end
    end

    methods (Access = private)
        function pm = newManager(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
        end

        function clearOutputs(testCase, bases)
            for k = 1:numel(bases)
                f = fullfile(testCase.SaveFolder, [bases{k} '.dat']);
                if isfile(f)
                    delete(f);
                end
            end
        end

        function runGeneratedScript(testCase, pm)
            scriptFolder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            scriptFile = fullfile(scriptFolder, 'sourceInfoScript.m');
            pm.generateScript(scriptFile);

            text = fileread(scriptFile);
            testCase.assertTrue(contains(text, "SaveFolder = '/path/to/save_folder';"), ...
                'expected SaveFolder placeholder in the generated script');
            text = strrep(text, "SaveFolder = '/path/to/save_folder';", ...
                "SaveFolder = '" + strrep(testCase.SaveFolder, '''', '''''') + "';");
            fid = fopen(scriptFile, 'w');
            fwrite(fid, text, 'char');
            fclose(fid);

            iRunScript(scriptFile);
        end

        function verifyAxes(testCase, fileName, names, sizes)
            f = fullfile(testCase.SaveFolder, fileName);
            testCase.assertTrue(isfile(f), sprintf('%s was not written', fileName));
            hdr = readDatHeader(f);
            testCase.verifyEqual(hdr.dimNames, names, sprintf('%s axes', fileName));
            testCase.verifyEqual(hdr.dimSizes, sizes, sprintf('%s sizes', fileName));
        end

        function verifyHeader(testCase, fileName, rate, exposure, expectedData)
            f = fullfile(testCase.SaveFolder, fileName);
            testCase.assertTrue(isfile(f), sprintf('%s was not written', fileName));
            testCase.assertTrue(isDatWithHeader(f), sprintf('%s must be headered', fileName));
            hdr = readDatHeader(f);
            [~, base] = fileparts(fileName);
            testCase.verifyEqual(hdr.frameRateHz, rate, sprintf('%s frame rate', fileName));
            if isnan(exposure)
                testCase.verifyTrue(isnan(hdr.exposureMsec), sprintf('%s exposure', fileName));
            else
                testCase.verifyEqual(hdr.exposureMsec, exposure, sprintf('%s exposure', fileName));
            end
            testCase.verifyEqual(hdr.channelName, base);
            testCase.verifyEqual(loadData(f), expectedData);
        end
    end
end

% =========================================================================
function iRunScript(scriptFile)
% Run the generated script in its own workspace.
run(scriptFile);
end

function parserFolder = iCreateFixtureCategory(projectRoot)
% PipelineManager discovers functions under Analysis/: copy the fixtures
% into a temporary category there (same approach as
% TestPipelineManagerLegacyMetaData).
analysisFolder = fullfile(projectRoot, 'Analysis');
parserFolder = tempname(analysisFolder);
mkdir(parserFolder);
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
names = {'pmSourceInfoScale', 'pmSourceInfoRate', 'pmSourceInfoBadMeta', 'pmSourceInfoMake', ...
    'pmSourceInfoMean', 'pmSourceInfoEvents', 'pmSourceInfoAxes'};
for k = 1:numel(names)
    copyfile(fullfile(fixtureFolder, [names{k} '.m']), fullfile(parserFolder, [names{k} '.m']));
end
% A second copy of the scaling step: the same function cannot appear twice
% in one branch.
source = fileread(fullfile(fixtureFolder, 'pmSourceInfoScale.m'));
source = strrep(source, "function outData = pmSourceInfoScale(", "function outData = pmSourceInfoScaleB(");
source = strrep(source, "'pmSourceInfoScale', 'Test", "'pmSourceInfoScaleB', 'Test");
fid = fopen(fullfile(parserFolder, 'pmSourceInfoScaleB.m'), 'w');
fwrite(fid, source, 'char');
fclose(fid);
addpath(parserFolder, '-begin');
clear pmSourceInfoScale pmSourceInfoScaleB pmSourceInfoRate pmSourceInfoBadMeta pmSourceInfoMake pmSourceInfoMean pmSourceInfoEvents pmSourceInfoAxes
rehash;
end

function iRemoveFixtureCategory(parserFolder, projectRoot)
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(parserFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerSourceInfo:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmSourceInfoScale pmSourceInfoScaleB pmSourceInfoRate pmSourceInfoBadMeta pmSourceInfoMake pmSourceInfoMean pmSourceInfoEvents pmSourceInfoAxes
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
end
