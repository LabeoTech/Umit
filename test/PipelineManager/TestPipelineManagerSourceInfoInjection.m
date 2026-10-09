classdef TestPipelineManagerSourceInfoInjection < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERSOURCEINFOINJECTION 'sourceInfo' inputs (.dat header Phase 6b-2).
    %
    %   Per-data metadata (frame rate, axes, exposure) reaches functions as
    %   explicit Name-Value parameters that PipelineManager injects from the
    %   data flowing into each step. The fixtures record what they receive
    %   in getappdata(0, 'pmInjLog').
    %
    %   Source files (AcqInfos.mat says 10 Hz, so a wrong source is visible):
    %       a.dat  Y-X-T   2x3x4    25 Hz, exposure 3
    %       b.dat  Y-X-T-E 2x3x2x2  40 Hz, exposure 7
    %       c.dat  Y-X-T   2x3x4    12 Hz, no exposure (NaN)
    %
    %   The complex pipeline holds several branches with different metadata
    %   in RAM at once and checks that every step receives its own branch's
    %   values, identically in auto mode, ramsafe mode, and a generated
    %   script.

    properties
        ProjectRoot char
        SaveFolder char
        DataA single
        DataB single
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
            testCase.DataA = single(reshape(1:24, [2 3 4]));
            testCase.DataB = single(reshape(101:124, [2 3 2 2]));

            AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                'FrameRateHz', 10, 'ExposureMsec', 1);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, [ ...
                struct('DatFile', 'a.dat', 'Length', 4, 'FrameRateHz', 10), ...
                struct('DatFile', 'b.dat', 'Length', 2, 'FrameRateHz', 10), ...
                struct('DatFile', 'c.dat', 'Length', 4, 'FrameRateHz', 10)]);
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            saveData(fullfile(testCase.SaveFolder, 'a.dat'), testCase.DataA, ...
                'DimNames', {'Y', 'X', 'T'}, 'Info', struct('frameRateHz', 25, 'exposureMsec', 3));
            saveData(fullfile(testCase.SaveFolder, 'b.dat'), testCase.DataB, ...
                'DimNames', {'Y', 'X', 'T', 'E'}, 'Info', struct('frameRateHz', 40, 'exposureMsec', 7));
            saveData(fullfile(testCase.SaveFolder, 'c.dat'), testCase.DataA, ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 12);

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));
            iResetLog();
            testCase.addTeardown(@iResetLog);
        end
    end

    methods (Test)
        % ------------------------------------------------ declaration
        function declarationsStayOutOfEditableLists(testCase)
            info = pmInjEcho('pipelineInfo');
            testCase.verifyEqual({info.sourceInfo.name}, {'FrameRateHz', 'DimNames', 'ExposureMsec'});
            testCase.verifyEqual([info.sourceInfo.required], [true false false]);
            testCase.verifyFalse(any(ismember({'FrameRateHz', 'DimNames', 'ExposureMsec'}, ...
                {info.inputs.name})), 'sourceInfo inputs must not be data or folder inputs');
            testCase.verifyTrue(isempty(info.parameters) || ~any(ismember( ...
                {'FrameRateHz', 'DimNames', 'ExposureMsec'}, {info.parameters.name})), ...
                'the config dialog lists info.parameters only');
            args = info.arguments(strcmp({info.arguments.kind}, 'sourceInfo'));
            testCase.verifyEqual({args.callType}, {'namevalue', 'namevalue', 'namevalue'});

            pm = testCase.newManager();
            pm.addStep('pmInjEcho', 'input', 'a.dat');
            node = pm.nodes(strcmp({pm.nodes.kind}, 'stream'));
            testCase.verifyTrue(isempty(node.info.parameters) || ...
                ~any(strcmp({node.info.parameters.name}, 'FrameRateHz')));
        end

        function invalidDeclarationsAreRejected(testCase)
            id = 'Umitoolbox:PipelineManager:invalidSourceInfoDeclaration';
            base = PipelineManager.createPipelineInfo('fx', 'fixture');
            base = PipelineManager.addInput(base, 'data', 'ImageTimeSeries', 'd', ...
                'position', 1, 'callType', 'positional', 'isData', true);

            testCase.verifyError(@() PipelineManager.addInput(base, 'FrameRateHz', 'sourceInfo', '', ...
                'kind', 'sourceInfo', 'sourceField', 'frameRate'), id);
            testCase.verifyError(@() PipelineManager.addInput(base, 'Rate', 'sourceInfo', '', ...
                'kind', 'sourceInfo', 'sourceField', 'frameRateHz'), id);

            bad = PipelineManager.addInput(base, 'FrameRateHz', 'sourceInfo', '', ...
                'kind', 'sourceInfo', 'sourceField', 'frameRateHz', 'sourceInput', 'ref');
            testCase.verifyError(@() PipelineManager.validateSourceInfoDecls(bad, 'fx'), id);

            noData = PipelineManager.createPipelineInfo('fx', 'fixture');
            noData = PipelineManager.addInput(noData, 'FrameRateHz', 'sourceInfo', '', ...
                'kind', 'sourceInfo', 'sourceField', 'frameRateHz');
            testCase.verifyError(@() PipelineManager.validateSourceInfoDecls(noData, 'fx'), id);

            good = PipelineManager.addInput(base, 'FrameRateHz', 'sourceInfo', '', ...
                'kind', 'sourceInfo', 'sourceField', 'frameRateHz', 'sourceInput', 'data');
            testCase.verifyWarningFree(@() PipelineManager.validateSourceInfoDecls(good, 'fx'));
        end

        % ------------------------------------------------ resolution
        function fileInputIsInjected(testCase)
            pm = testCase.newManager();
            pm.addStep('pmInjEcho', 'input', 'a.dat');
            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            rec = testCase.recordOf('pmInjEcho');
            testCase.verifyEqual(rec.FrameRateHz, 25);
            testCase.verifyEqual(rec.DimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(rec.ExposureMsec, 3);
        end

        function optionalUnresolvedFieldIsOmitted(testCase)
            pm = testCase.newManager();
            pm.addStep('pmInjEcho', 'input', 'c.dat');
            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            rec = testCase.recordOf('pmInjEcho');
            testCase.verifyEqual(rec.FrameRateHz, 12);
            testCase.verifyEmpty(rec.ExposureMsec, 'a NaN exposure must not be passed');
        end

        function upstreamMetaDataIsDeferredAndInjectedAtRunTime(testCase)
            pm = testCase.newManager();
            rateTag = pm.addStep('pmSourceInfoRate', 'input', 'a.dat');
            pm.addStep('pmInjEcho', 'input', rateTag);

            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyTrue(isValid, strjoin(report.errors, ' | '));
            testCase.verifyTrue(any(contains(string(report.warnings), 'resolved at run time')));

            result = pm.executePipeline('PrintSummary', false);
            testCase.assertEqual(result.status, "completed");
            rec = testCase.recordOf('pmInjEcho');
            testCase.verifyEqual(rec.FrameRateHz, 5, 'the metaData update must reach the step');
            testCase.verifyEqual(rec.ExposureMsec, 9);
        end

        function requiredUnresolvedFieldFailsPreflight(testCase)
            pm = testCase.newManager();
            makeTag = pm.addStep('pmSourceInfoMake');
            pm.addStep('pmInjEcho', 'input', makeTag);

            [isValid, report] = pm.validateNodePreflight('throwOnError', false);
            testCase.verifyFalse(isValid);
            msg = strjoin(report.errors, ' | ');
            testCase.verifySubstring(msg, 'pmInjEcho');
            testCase.verifySubstring(msg, '''FrameRateHz''');
            testCase.verifySubstring(msg, 'input "data"');
            testCase.verifySubstring(msg, 'frameRateHz');
            testCase.verifySubstring(msg, 'no source file');

            testCase.verifyError(@() pm.executePipeline('PrintSummary', false), ...
                'PipelineManager:validateNodePreflight:Failed');
            testCase.verifyEmpty(getappdata(0, 'pmInjLog'), 'the step must not run');
        end

        function generatedScriptRaisesForUnresolvedRequiredField(testCase)
            pm = testCase.newManager();
            makeTag = pm.addStep('pmSourceInfoMake');
            pm.addStep('pmInjEcho', 'input', makeTag);
            testCase.verifyError(@() testCase.runGeneratedScript(pm), ...
                'Umitoolbox:PipelineManager:unresolvedSourceInfo');
        end

        function nonPrimarySourceInputIsUsed(testCase)
            pm = testCase.newManager();
            pm.addStep('pmInjTwoInputs', 'input', {'a.dat', 'b.dat'});
            result = pm.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            rec = testCase.recordOf('pmInjTwoInputs');
            testCase.verifyEqual(rec.FrameRateHz, 40, 'FrameRateHz comes from input "ref"');
            testCase.verifyEqual(rec.ExposureMsec, 3, 'ExposureMsec comes from input "data"');
        end

        % ------------------------------------------------ complex pipelines
        function branchesKeepTheirOwnMetadataInAutoMode(testCase)
            pm = testCase.buildComplexPipeline();
            result = pm.executePipeline('PrintSummary', false);
            testCase.assertEqual(result.status, "completed");
            testCase.verifyComplexLog();
        end

        function branchesKeepTheirOwnMetadataInRamSafeMode(testCase)
            pm = testCase.buildComplexPipeline();
            pm.ramMode = 'ramsafe';
            result = pm.executePipeline('PrintSummary', false);
            testCase.assertEqual(result.status, "completed");
            testCase.verifyComplexLog();
        end

        function generatedScriptInjectsTheSameValues(testCase)
            pm = testCase.buildComplexPipeline();
            testCase.runGeneratedScript(pm);
            testCase.verifyComplexLog();
        end

        % ------------------------------------------------ saved pipelines
        function olderSavedPipelineGetsCurrentDeclarations(testCase)
            pm = testCase.newManager();
            pm.addStep('pmInjEcho', 'input', 'a.dat');
            pipeFile = fullfile(fileparts(testCase.SaveFolder), 'older.pipe');
            pm.savePipe(pipeFile);

            loaded = load(pipeFile, 'pipeStruct', '-mat');
            pipeStruct = loaded.pipeStruct;
            idx = find(strcmpi({pipeStruct.nodes.kind}, 'stream'), 1);
            info = pipeStruct.nodes(idx).info;
            info = rmfield(info, 'sourceInfo');
            info.arguments = info.arguments(~strcmp({info.arguments.kind}, 'sourceInfo'));
            pipeStruct.nodes(idx).info = info;
            save(pipeFile, 'pipeStruct', '-mat');

            restored = testCase.newManager();
            restored.loadPipe(pipeFile);
            result = restored.executePipeline('PrintSummary', false);

            testCase.assertEqual(result.status, "completed");
            testCase.verifyEqual(testCase.recordOf('pmInjEcho').FrameRateHz, 25);
        end
    end

    methods (Access = private)
        function pm = newManager(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
        end

        function pm = buildComplexPipeline(testCase)
            % Branch A: a.dat -> echo -> scale -> echoB          (25 Hz, 3, Y-X-T)
            % Branch B: b.dat -> echoC                           (40 Hz, 7, Y-X-T-E)
            % Fork of A: a.dat -> rate (metaData 5 Hz, 9) -> echoD
            % Join:  twoInputs(echo output, echoC output): rate from B, exposure from A
            pm = testCase.newManager();
            echoA = pm.addStep('pmInjEcho', 'input', 'a.dat');
            scaled = pm.addStep('pmSourceInfoScale', 'input', echoA);
            pm.addStep('pmInjEchoB', 'input', scaled, 'saveas', 'branchA.dat');
            echoB = pm.addStep('pmInjEchoC', 'input', 'b.dat');
            rated = pm.addStep('pmSourceInfoRate', 'input', 'a.dat');
            pm.addStep('pmInjEchoD', 'input', rated, 'saveas', 'fork.dat');
            pm.addStep('pmInjTwoInputs', 'input', {echoA, echoB}, 'saveas', 'joined.dat');
        end

        function verifyComplexLog(testCase)
            log = getappdata(0, 'pmInjLog');
            testCase.assertEqual(sort({log.fcn}), sort({'pmInjEcho', 'pmInjEchoB', 'pmInjEchoC', ...
                'pmInjEchoD', 'pmInjTwoInputs'}), 'every step runs once');
            expected = { ...
                'pmInjEcho',      25, {'Y', 'X', 'T'},      3; ...
                'pmInjEchoB',     25, {'Y', 'X', 'T'},      3; ...
                'pmInjEchoC',     40, {'Y', 'X', 'T', 'E'}, 7; ...
                'pmInjEchoD',      5, {'Y', 'X', 'T'},      9};
            for k = 1:size(expected, 1)
                rec = testCase.recordOf(expected{k, 1});
                testCase.verifyEqual(rec.FrameRateHz, expected{k, 2}, [expected{k, 1} ' FrameRateHz']);
                testCase.verifyEqual(rec.DimNames, expected{k, 3}, [expected{k, 1} ' DimNames']);
                testCase.verifyEqual(rec.ExposureMsec, expected{k, 4}, [expected{k, 1} ' ExposureMsec']);
            end
            rec = testCase.recordOf('pmInjTwoInputs');
            testCase.verifyEqual(rec.FrameRateHz, 40, 'join: rate of the ref branch (b.dat)');
            testCase.verifyEqual(rec.ExposureMsec, 3, 'join: exposure of the data branch (a.dat)');
            testCase.verifyEqual(readDatHeader(fullfile(testCase.SaveFolder, 'fork.dat')).frameRateHz, 5);
            testCase.verifyEqual(readDatHeader(fullfile(testCase.SaveFolder, 'branchA.dat')).frameRateHz, 25);
            testCase.verifyEqual(readDatHeader(fullfile(testCase.SaveFolder, 'joined.dat')).frameRateHz, 25);
        end

        function rec = recordOf(testCase, fcn)
            log = getappdata(0, 'pmInjLog');
            testCase.assertNotEmpty(log, 'no step recorded anything');
            idx = find(strcmp({log.fcn}, fcn));
            testCase.assertNumElements(idx, 1, sprintf('%s must run exactly once', fcn));
            rec = log(idx);
        end

        function runGeneratedScript(testCase, pm)
            scriptFolder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            scriptFile = fullfile(scriptFolder, 'injectionScript.m');
            pm.generateScript(scriptFile);
            text = fileread(scriptFile);
            testCase.assertTrue(contains(text, "SaveFolder = '/path/to/save_folder';"));
            text = strrep(text, "SaveFolder = '/path/to/save_folder';", ...
                "SaveFolder = '" + strrep(testCase.SaveFolder, '''', '''''') + "';");
            fid = fopen(scriptFile, 'w');
            fwrite(fid, text, 'char');
            fclose(fid);
            iRunScript(scriptFile);
        end
    end
end

% =========================================================================
function iRunScript(scriptFile)
run(scriptFile);
end

function iResetLog()
if isappdata(0, 'pmInjLog')
    rmappdata(0, 'pmInjLog');
end
end

function parserFolder = iCreateFixtureCategory(projectRoot)
% PipelineManager discovers functions under Analysis/: copy the fixtures
% into a temporary category there. Copies of pmInjEcho let the same step
% appear on several branches.
analysisFolder = fullfile(projectRoot, 'Analysis');
parserFolder = tempname(analysisFolder);
mkdir(parserFolder);
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
names = {'pmInjEcho', 'pmInjTwoInputs', 'pmSourceInfoScale', 'pmSourceInfoRate', 'pmSourceInfoMake'};
for k = 1:numel(names)
    copyfile(fullfile(fixtureFolder, [names{k} '.m']), fullfile(parserFolder, [names{k} '.m']));
end
source = fileread(fullfile(fixtureFolder, 'pmInjEcho.m'));
for suffix = {'B', 'C', 'D'}
    newName = ['pmInjEcho' suffix{1}];
    copyText = strrep(source, 'function outData = pmInjEcho(', ['function outData = ' newName '(']);
    fid = fopen(fullfile(parserFolder, [newName '.m']), 'w');
    fwrite(fid, copyText, 'char');
    fclose(fid);
end
addpath(parserFolder, '-begin');
clear pmInjEcho pmInjEchoB pmInjEchoC pmInjEchoD pmInjTwoInputs pmSourceInfoScale pmSourceInfoRate pmSourceInfoMake
rehash;
end

function iRemoveFixtureCategory(parserFolder, projectRoot)
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(parserFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerSourceInfoInjection:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmInjEcho pmInjEchoB pmInjEchoC pmInjEchoD pmInjTwoInputs pmSourceInfoScale pmSourceInfoRate pmSourceInfoMake
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
end
