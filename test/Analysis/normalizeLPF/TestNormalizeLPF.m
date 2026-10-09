classdef TestNormalizeLPF < matlab.unittest.TestCase
    %TESTNORMALIZELPF Unit tests for normalizeLPF.
    %
    % Coverage:
    %   1) pipelineInfo generation
    %   2) Raw YXT array inputs, with and without NaNs
    %   3) Raw YXTE (event-split) array inputs: each trial filtered on its own
    %   4) Raw .dat filename inputs (YXT and YXTE)
    %   5) Streamed .dat results across forced multi-slab runs
    %   6) Rejection of UMT structs, .umt/.mat files, and unsupported layouts
    %   7) Output representation preservation:
    %        - array in  -> array out, same size
    %        - .dat in   -> .dat out, same axes and sizes
    %
    % Fixture policy:
    %   - Sample files are copied from:
    %         <projectRoot>\test\Analysis\TestingData_with_events
    %   - A fresh temporary SaveFolder is used per test.
    %   - events.mat is recreated before each test.
    %   - The real recording is retained because this suite verifies
    %     numerical equivalence across array and DAT representations through
    %     the legacy filter core.

    properties
        SaveFolder char
        ProjectRoot char
        SampleDataFolder char
    end

    methods (TestClassSetup)
        function addProjectPathAndLocateFixture(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);

            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            cfg = [];
            if isappdata(0, 'NormalizeLPFTestConfig')
                cfg = getappdata(0, 'NormalizeLPFTestConfig');
            end

            if ~isempty(cfg) && isfield(cfg, 'sampleDataFolder') && ...
                    isfolder(cfg.sampleDataFolder)
                testCase.SampleDataFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.SampleDataFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', ...
                    'Analysis', ...
                    'TestingData_with_events');
            end

            testCase.assertTrue( ...
                isfile(fullfile(testCase.SampleDataFolder, 'green.dat')), ...
                ['Missing fixture file "green.dat". Put it in: ' ...
                 testCase.SampleDataFolder]);

            testCase.assertTrue( ...
                isfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat')), ...
                ['Missing fixture file "AcqInfos.mat". Put it in: ' ...
                 testCase.SampleDataFolder]);
        end
    end

    methods (TestMethodSetup)
        function createFreshWorkspace(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.SaveFolder = fx.Folder;

            % The exponential-fit path seeds fminsearch through rng('shuffle')
            % (DFR-20260819-001), so two independent calls on the same input
            % can converge to different optima. Shadow rng with the fixed-seed
            % test double so every comparison in this class is reproducible.
            projectRoot = extractBefore(mfilename('fullpath'), [filesep 'test' filesep]);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                projectRoot, 'test', 'IOIAnalysis', 'NormalisationFiltering', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_RNG_MODE', 'fixed'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'green.dat'), ...
                fullfile(testCase.SaveFolder, 'green.dat'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.SaveFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'normLPF.dat'));

            testCase.createFreshEvents(0.8);
        end
    end

    methods (TestMethodTeardown)
        function cleanupCreatedFiles(testCase)
            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'normLPF.dat'));
        end
    end

    methods (Test)
        function testPipelineInfoConstrainsBaselineCutoff(testCase)
            info = normalizeLPF('pipelineInfo');

            baselineCutoff = info.parameters(strcmp( ...
                {info.parameters.name}, 'BaselineCutoffHz'));
            testCase.verifyEqual(baselineCutoff.allowed, [0 Inf]);
        end

        function testArrayInputRecording(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));

            expected = testCase.expectedFilteredArray(rawData, ...
                0.0083, 1, true, false, testCase.getFrameRateHz());

            out = normalizeLPF(rawData, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, ...
                'SignalCutoffHz', 1, ...
                'Normalize', true, ...
                'bApplyExpFit', false, ...
                'FrameRateHz', testCase.acqFrameRateHz());

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testArrayInputRecordingWithNaNs(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            rawData(1,1,1) = NaN;
            rawData(2,3,5) = NaN;

            expected = testCase.expectedFilteredArray(rawData, ...
                0.01, 0.8, true, true, testCase.getFrameRateHz());

            out = normalizeLPF(rawData, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.01, ...
                'SignalCutoffHz', 0.8, ...
                'Normalize', true, ...
                'bApplyExpFit', true, ...
                'FrameRateHz', testCase.acqFrameRateHz());

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyTrue(isnan(out(1,1,1)));
            testCase.verifyTrue(isnan(out(2,3,5)));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testArrayInputYXTEFiltersEachTrialAlongT(testCase)
            trials = testCase.buildEventSplitData();
            expected = testCase.expectedFilteredTrials(trials, 0.0083, 1, true, false);

            out = normalizeLPF(trials, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, ...
                'SignalCutoffHz', 1, ...
                'Normalize', true, ...
                'bApplyExpFit', false, ...
                'FrameRateHz', testCase.acqFrameRateHz());

            testCase.verifyClass(out, 'single');
            testCase.verifyEqual(size(out), size(trials));
            testCase.verifyNumericEquivalent(out, expected);
        end

        function testDatInputRecording(testCase)
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));

            expected = testCase.expectedFilteredArray(rawData, ...
                0.0083, 1, true, false, testCase.getFrameRateHz());

            out = normalizeLPF(inFile, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, ...
                'SignalCutoffHz', 1, ...
                'Normalize', true, ...
                'bApplyExpFit', false);

            testCase.verifyTrue(ischar(out) || (isstring(out) && isscalar(out)));
            outFile = char(string(out));
            testCase.verifyTrue(isfile(outFile));

            actual = single(loadData(outFile));
            testCase.verifyEqual(size(actual), size(rawData));
            testCase.verifyNumericEquivalent(actual, single(expected));
        end

        function testDatInputYXTEKeepsTheAxesAndFiltersEachTrial(testCase)
            trials = testCase.buildEventSplitData();
            inFile = testCase.writeEventSplitDat(trials);

            out = normalizeLPF(inFile, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, ...
                'SignalCutoffHz', 1, ...
                'Normalize', true, ...
                'bApplyExpFit', false);

            outFile = char(string(out));
            outInfo = loadMetaData(outFile);
            testCase.verifyEqual(outInfo.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(outInfo.dimSizes(:)).', size(trials));
            testCase.verifyNumericEquivalent(single(loadData(outFile)), ...
                testCase.expectedFilteredTrials(trials, 0.0083, 1, true, false));
        end

        function testForcedMultiSlabYXTEMatchesInRam(testCase)
            % The memory mock makes the core ask for many X slabs: the
            % streamed result must equal the in-RAM one, without and with the
            % per-trial exponential fit (whose random start is pinned by the
            % rng mock so the two runs are comparable).
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'IOIAnalysis', 'NormalisationFiltering', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_RNG_MODE', 'fixed'));

            trials = testCase.buildEventSplitData();
            inFile = testCase.writeEventSplitDat(trials); %#ok<NASGU> used inside evalc
            rate = testCase.acqFrameRateHz();

            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            for bExpFit = [false true]
                inRam = normalizeLPF(trials, testCase.SaveFolder, ...
                    'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                    'Normalize', true, 'bApplyExpFit', bExpFit, 'FrameRateHz', rate);

                progress = evalc(['outFile = normalizeLPF(inFile, testCase.SaveFolder, ' ...
                    '''BaselineCutoffHz'', 0.0083, ''SignalCutoffHz'', 1, ' ...
                    '''Normalize'', true, ''bApplyExpFit'', bExpFit);']);

                slabCounts = cellfun(@(c) str2double(c{1}), ...
                    regexp(progress, 'Filtering event-split data \((\d+) chunk', 'tokens'));
                testCase.verifyGreaterThan(max(slabCounts), 1, ...
                    'The fixture must force more than one X slab.');
                testCase.verifyNumericEquivalent( ...
                    single(loadData(char(string(outFile)))), inRam);
            end
        end

        function testRejectsUMTAndUnsupportedInputs(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            rate = testCase.acqFrameRateHz();
            umt = genUMTStruct(rawData, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.SaveFolder, 'input.umt');
            saveData(umtFile, umt);
            matFile = fullfile(testCase.SaveFolder, 'input.mat');
            umtInput = umt;
            save(matFile, 'umtInput');
            yxeFile = fullfile(testCase.SaveFolder, 'yxe.dat');
            saveData(yxeFile, rawData(:, :, 1:3), 'DimNames', {'Y','X','E'}, ...
                'FrameRateHz', rate);

            testCase.verifyError(@() normalizeLPF(umt, testCase.SaveFolder), ...
                'normalizeLPF:UnsupportedInputType');
            testCase.verifyError(@() normalizeLPF(umtFile, testCase.SaveFolder), ...
                'normalizeLPF:UnsupportedInputFile');
            testCase.verifyError(@() normalizeLPF(matFile, testCase.SaveFolder), ...
                'normalizeLPF:UnsupportedInputFile');
            testCase.verifyError(@() normalizeLPF(yxeFile, testCase.SaveFolder), ...
                'Umitoolbox:normalizeLPF:unsupportedLayout');
            testCase.verifyError(@() normalizeLPF(rawData(:, :, 1), testCase.SaveFolder, ...
                'FrameRateHz', rate), 'normalizeLPF:InvalidArrayInput');
        end

        function testInRamInputWithoutFrameRateErrors(testCase)
            % Numeric input needs 'FrameRateHz'; AcqInfos.mat in the folder
            % is not used for it.
            testCase.assertTrue(isfile(fullfile(testCase.SaveFolder, 'AcqInfos.mat')));
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));

            testCase.verifyError(@() normalizeLPF(rawData, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, ...
                'SignalCutoffHz', 1, ...
                'Normalize', true, ...
                'bApplyExpFit', false), ...
                'Umitoolbox:normalizeLPF:missingFrameRateHz');
        end

        function testExplicitFrameRateDifferentFromDatHeaderWarns(testCase)
            % An explicit 'FrameRateHz' that differs from the .dat header
            % warns sourceInfoConflict (the explicit value is used).
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            headerRateHz = testCase.getFrameRateHz();

            testCase.verifyWarning(@() normalizeLPF(inFile, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, ...
                'SignalCutoffHz', 1, ...
                'Normalize', true, ...
                'bApplyExpFit', false, ...
                'FrameRateHz', 2 * headerRateHz), ...
                'Umitoolbox:normalizeLPF:sourceInfoConflict');
        end
    end

    methods (Access = private)

        function frameRateHz = acqFrameRateHz(testCase)
            % Frame rate of the fixture dataset (from its AcqInfos.mat).
            acq = load(fullfile(testCase.SaveFolder, 'AcqInfos.mat'));
            if isfield(acq, 'AcqInfoStream')
                acqInfo = acq.AcqInfoStream;
            else
                fn = fieldnames(acq);
                acqInfo = acq.(fn{1});
            end
            frameRateHz = double(acqInfo.FrameRateHz);
        end

        function createFreshEvents(testCase, baselinePeriodSec)
            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));

            acq = load(fullfile(testCase.SaveFolder, 'AcqInfos.mat'));
            if isfield(acq, 'AcqInfoStream')
                acqInfo = acq.AcqInfoStream;
            else
                fn = fieldnames(acq);
                acqInfo = acq.(fn{1});
            end

            frameRate = double(acqInfo.FrameRateHz);
            nFrames = double(acqInfo.Length);

            signalSR = 1000;
            durationS = max(nFrames / frameRate, 10);
            nSamples = ceil((durationS + 1) * signalSR);
            triggerSignal = zeros(nSamples, 1, 'single');

            onsetFrames = round(linspace( ...
                max(3, 0.15 * nFrames), ...
                max(4, 0.75 * nFrames), ...
                4));

            onsetFrames = unique(onsetFrames(:)', 'stable');
            onsetTimes = (onsetFrames - 1) ./ frameRate;
            pulseWidthS = 0.40;

            for iOn = 1:numel(onsetTimes)
                s0 = max(1, round(onsetTimes(iOn) * signalSR) + 1);
                s1 = min(nSamples, s0 + round(pulseWidthS * signalSR) - 1);
                triggerSignal(s0:s1) = 1;
            end

            evObj = EventsManager(testCase.SaveFolder, '', 'csv');
            evObj.getTriggersFromSignal(triggerSignal, signalSR, false);
            evObj.setBaselinePeriod(single(baselinePeriodSec));
            evObj.saveEvents(testCase.SaveFolder);

            testCase.verifyTrue(isfile(fullfile(testCase.SaveFolder, 'events.mat')));
        end

        function freqHz = getFrameRateHz(testCase)
            meta = loadMetaData(fullfile(testCase.SaveFolder, 'green.dat'));
            freqHz = double(meta.frameRateHz);
        end

        function trials = buildEventSplitData(testCase)
            % Realistic event-split data: the fixture recording split by the
            % events of events.mat (trials cropped to the shortest one).
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            trials = single(EventsManager(testCase.SaveFolder).splitDataByEvents( ...
                rawData, 'FrameRateHz', testCase.acqFrameRateHz(), 'IncludeIgnored', true));
            testCase.assertEqual(ndims(trials), 4);
        end

        function inFile = writeEventSplitDat(testCase, trials)
            inFile = fullfile(testCase.SaveFolder, 'byEvent.dat');
            saveData(inFile, trials, 'DimNames', {'Y','X','T','E'}, ...
                'FrameRateHz', testCase.acqFrameRateHz());
        end

        function out = expectedFilteredTrials(testCase, trials, baselineCutoffHz, signalCutoffHz, bNormalize, bApplyExpFit)
            % Reference: the legacy core applied to every trial on its own.
            out = zeros(size(trials), 'single');
            for iTrial = 1:size(trials, 4)
                out(:,:,:,iTrial) = testCase.expectedFilteredArray( ...
                    trials(:,:,:,iTrial), baselineCutoffHz, signalCutoffHz, ...
                    bNormalize, bApplyExpFit, testCase.getFrameRateHz());
            end
        end

        function out = expectedFilteredArray(~, inArray, baselineCutoffHz, signalCutoffHz, bNormalize, bApplyExpFit, Fs)
            if ~isa(inArray, 'single')
                workArray = single(inArray);
            else
                workArray = inArray;
            end

            idxNaN = isnan(workArray);
            if any(idxNaN(:))
                workArray(idxNaN) = 0;
            end

            out = NormalisationFiltering( ...
                pwd, workArray, ...
                baselineCutoffHz, ...
                signalCutoffHz, ...
                bNormalize, ...
                bApplyExpFit, ...
                Fs);

            out = single(out);

            if any(idxNaN(:))
                out(idxNaN) = NaN;
            end
        end

        function verifyNumericEquivalent(testCase, a, b)
            testCase.verifyEqual(size(a), size(b));

            if isequaln(a, b)
                return
            end

            diffVals = double(a(:)) - double(b(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                testCase.verifyEqual(a, b);
                return
            end

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 2e-3);
        end
    end
end

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
