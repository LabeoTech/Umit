classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestNormalizeLPFModes < matlab.unittest.TestCase
    %TESTNORMALIZELPFMODES normalizeLPF input modes against the legacy filter.
    %
    % Formerly TestNormalisationFiltering, the suite of the IOIAnalysis
    % core absorbed into normalizeLPF (2026-10-13).
    %
    % Coverage:
    %   1) In-memory YXT input with explicit frame rate
    %   2) File mode with relative filename
    %   3) File mode with absolute filename
    %   4) File mode with Normalize=false and bApplyExpFit=true
    %   5) Legacy sidecar metadata precedence over AcqInfos.mat
    %   6) Event-split file mode
    %   7) F-13: chunked file-mode bExpFit output is invariant to nChunks
    %   8) F-14: chunked file-mode NaN masking preserves NaN positions and
    %      roughly agrees with the in-RAM path
    %
    % Notes:
    %   - The goal of this suite is to validate input combinations and
    %     compatibility behavior, not to re-derive the filtering algorithm.
    %   - Expected outputs come from legacyNormFiltDirect (private), a frozen
    %     copy of the legacy in-memory filter, and are compared to the other
    %     input combinations.

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
            if isappdata(0, 'NormalisationFilteringTestConfig')
                cfg = getappdata(0, 'NormalisationFilteringTestConfig');
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
                projectRoot, 'test', 'Analysis', 'normalizeLPF', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_RNG_MODE', 'fixed'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'green.dat'), ...
                fullfile(testCase.SaveFolder, 'green.dat'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.SaveFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_info.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_NormFilt.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_NormFilt.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'customNorm.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'customNorm.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_events.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_events.mat'));

            testCase.createFreshEvents(0.8);
        end
    end

    methods (TestMethodTeardown)
        function cleanupCreatedFiles(testCase)
            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_info.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_NormFilt.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_NormFilt.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'customNorm.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'customNorm.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_events.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'green_events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'greenNaN.dat'));
        end
    end

    methods (Test)

        function testDirectArrayWithExplicitFreq(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            Fs = testCase.getFrameRateHz();

            out = normalizeLPF(rawData, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', false, 'FrameRateHz', Fs);

            expected = testCase.expectedDirect(rawData, 0.0083, 1, true, false, Fs);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testFileModeRelativeName(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            Fs = testCase.getFrameRateHz();

            out = loadData(normalizeLPF('green.dat', testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', false));

            expected = testCase.expectedDirect(rawData, 0.0083, 1, true, false, Fs);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testFileModeAbsolutePath(testCase)
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));
            Fs = testCase.getFrameRateHz();

            out = loadData(normalizeLPF(inFile, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', false));

            expected = testCase.expectedDirect(rawData, 0.0083, 1, true, false, Fs);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testFileModeSubtractWithExpFit(testCase)
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));
            Fs = testCase.getFrameRateHz();

            outFile = normalizeLPF('green.dat', testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.02, 'SignalCutoffHz', 0.9, ...
                'Normalize', false, 'bApplyExpFit', true);

            testCase.verifyTrue(isfile(outFile));

            actual = single(loadData(outFile));
            expected = testCase.expectedDirect(rawData, 0.02, 0.9, false, true, Fs);

            testCase.verifyEqual(size(actual), size(rawData));
            testCase.verifyNumericEquivalent(actual, single(expected));
        end

        function testLegacySidecarFreqTakesPrecedence(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));

            legacyFreq = 1.25; % deliberately different from AcqInfos.mat
            legacyMeta = struct();
            legacyMeta.dim_names = {'Y','X','T'};
            legacyMeta.datSize = [size(rawData,1), size(rawData,2)];
            legacyMeta.datLength = size(rawData,3);
            legacyMeta.Freq = legacyFreq;
            legacyMeta.Datatype = 'single';

            % A legacy (headerless) file: the fixture is headered since .dat
            % header Phase 5a, and a header would take precedence.
            writeTestDat(fullfile(testCase.SaveFolder, 'green.dat'), rawData, legacyFreq, ...
                'Format', 'legacySidecar');
            save(fullfile(testCase.SaveFolder, 'green.mat'), '-struct', 'legacyMeta');

            out = loadData(normalizeLPF('green.dat', testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 0.5, ...
                'Normalize', true, 'bApplyExpFit', false));

            expected = testCase.expectedDirect(rawData, 0.0083, 0.5, true, false, legacyFreq);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testEventSplitFileMode(testCase)
            [eventDatFile, eventData, legacyFreq] = testCase.buildEventSplitDatFixture();

            out = loadData(normalizeLPF(eventDatFile, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', false));

            expected = zeros(size(eventData), 'single');
            for iTrial = 1:size(eventData, 4)
                expected(:,:,:,iTrial) = testCase.expectedDirect( ...
                    eventData(:,:,:,iTrial), ...
                    0.0083, 1, true, false, legacyFreq);
            end

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(eventData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testEventSplitFileModeWritesHeaderedYXTE(testCase)
            % .dat header Phase 4c-2b: the file-mode output keeps the
            % input's Y,X,T,E axes and frame rate.
            [eventDatFile, eventData, Fs] = testCase.buildEventSplitDatFixture();
            outFile = fullfile(testCase.SaveFolder, 'normLPF.dat');
            testCase.addTeardown(@() deleteIfExists(outFile));

            normalizeLPF(eventDatFile, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', false);

            testCase.assertTrue(isDatWithHeader(outFile));
            hdr = readDatHeader(outFile);
            testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T', 'E'});
            testCase.verifyEqual(hdr.dimSizes, size(eventData));
            testCase.verifyEqual(hdr.frameRateHz, double(single(Fs)));
            testCase.verifyEqual(hdr.channelName, 'normLPF');
            testCase.verifyTrue(hdr.writeComplete);
        end

        function testEventSplitStreamedResultIsInvariantToChunkCount(testCase)
            % The Y,X,T,E file mode streams X slabs. Pixels and trials are
            % independent along T, and the exponential fit is one whole-image
            % fit per trial (a pre-pass over the slabs), so one slab and many
            % slabs must give the same output, which must also equal the
            % core's own per-trial array mode. The rng mock is pinned so the
            % randomized fminsearch start is reproducible.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk-count forcing relies on shadowing the PCWIN64 memory() built-in.');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'Analysis', 'normalizeLPF', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_RNG_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));

            [eventDatFile, eventData, Fs] = testCase.buildEventSplitDatFixture();

            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(1e12, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(1e12, '%d')));
            lpfArgs = {'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', true};
            outSingleChunk = loadData(normalizeLPF( ...
                eventDatFile, testCase.SaveFolder, lpfArgs{:}));

            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));
            progress = evalc(['outManyChunks = loadData(normalizeLPF(' ...
                'eventDatFile, testCase.SaveFolder, lpfArgs{:}));']);

            nChunks = str2double(regexp(progress, ...
                'Filtering event-split data \((\d+) chunk', 'tokens', 'once'));
            testCase.verifyGreaterThan(nChunks, 1, ...
                'The fixture must force more than one X slab.');
            % Same fit, same per-pixel filtering: slabs of different widths
            % may flip the last bit of a single value, nothing more.
            testCase.verifyEqual(double(outManyChunks), double(outSingleChunk), 'AbsTol', 1e-5);

            expected = zeros(size(eventData), 'single');
            for iTrial = 1:size(eventData, 4)
                expected(:,:,:,iTrial) = testCase.expectedDirect( ...
                    eventData(:,:,:,iTrial), 0.0083, 1, true, true, Fs);
            end
            testCase.verifyNumericEquivalent(single(outSingleChunk), expected);
        end

        function testChunkCountInvariantForGlobalExpFit(testCase)
            % F-13 regression: the double-exponential detrend fit used to be
            % recomputed per X-chunk on the chunked .dat path, so the output
            % depended on nChunks, which itself depends on available RAM via
            % calculateMaxChunkSize. Force two very different chunk counts
            % (1 vs. many) via the memory() mock and confirm identical
            % output now that the fit is computed once, globally.
            %
            % fminsearch's initial guess is drawn via rng('shuffle')/rand,
            % so two independent invocations are not bit-comparable on their
            % own (tracked separately as DFR-20260819-001); shadow rng() to
            % force the same deterministic seed on both calls so this test
            % isolates the chunk-count effect alone.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk-count forcing relies on shadowing the PCWIN64 memory() built-in.');

            rngMocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'Analysis', 'normalizeLPF', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rngMocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_RNG_MODE', 'fixed'));

            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));

            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(1e12, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(1e12, '%d')));
            lpfArgs = {'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', true};
            outSingleChunk = loadData(normalizeLPF( ...
                'green.dat', testCase.SaveFolder, lpfArgs{:}));

            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(2e9, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(6.13e8, '%d')));
            outManyChunks = loadData(normalizeLPF( ...
                'green.dat', testCase.SaveFolder, lpfArgs{:}));

            testCase.verifyEqual(outManyChunks, outSingleChunk);
        end

        function testDatFileWithNaNPreservesPositionsAndMatchesInRAMPath(testCase)
            % F-14 regression: the chunked .dat path used to send NaNs
            % straight into filtfilt/backslash, poisoning the whole temporal
            % trace of every affected pixel. Mask a few scattered (Y,X,T)
            % elements as NaN -- the same style TestNormalizeLPF already
            % uses for the in-RAM path -- run the file-mode (chunked) path
            % with bExpFit on, and confirm NaN positions in the output match
            % the input exactly, with no leakage into pixels that were valid
            % on input. Also cross-check against the in-RAM path
            % (normalizeLPF's array case), which already applies the same
            % NaN policy.
            %
            % Note: a pixel masked across *every* T sample would zero out
            % its whole temporal trace, and the exponential-fit division
            % (Signal ./ Approx) resolves that all-zero trace to 0/0 = NaN
            % on both paths alike -- a pre-existing, out-of-scope limitation
            % of the bExpFit division shared by the array and file paths,
            % not something F-13/F-14 touches. Scattered single-sample NaNs
            % avoid that edge case while still exercising the masking policy.
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            Fs = testCase.getFrameRateHz();

            maskedData = rawData;
            maskedData(1,1,1) = NaN;
            maskedData(2,3,5) = NaN;
            maskedData(10,20,500) = NaN;
            idxNaN = isnan(maskedData);

            naNFile = fullfile(testCase.SaveFolder, 'greenNaN.dat');
            % Headered input (.dat header Phase 5a).
            writeTestDat(naNFile, maskedData, Fs);

            outFile = loadData(normalizeLPF('greenNaN.dat', testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', true, 'FrameRateHz', Fs));

            testCase.verifyEqual(isnan(outFile), idxNaN);

            % In-RAM input takes its frame rate explicitly (.dat header
            % Phase 6b-2): the same Fs the file path uses.
            inRAM = normalizeLPF(maskedData, testCase.SaveFolder, ...
                'BaselineCutoffHz', 0.0083, 'SignalCutoffHz', 1, ...
                'Normalize', true, 'bApplyExpFit', true, 'FrameRateHz', Fs);

            testCase.verifyEqual(isnan(inRAM), idxNaN);
            testCase.verifyNumericEquivalent(single(outFile), single(inRAM));

            deleteIfExists(naNFile);
        end
    end

    methods (Access = private)

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

        function Fs = getFrameRateHz(testCase)
            S = load(fullfile(testCase.SaveFolder, 'AcqInfos.mat'));
            if isfield(S, 'AcqInfoStream')
                acqInfo = S.AcqInfoStream;
            else
                fn = fieldnames(S);
                acqInfo = S.(fn{1});
            end
            Fs = double(acqInfo.FrameRateHz);
        end

        function expected = expectedDirect(~, inArray, lowFreq, highFreq, bDivide, bExpFit, Fs)
            expected = legacyNormFiltDirect(single(inArray), ...
                lowFreq, highFreq, bDivide, bExpFit, Fs);
            expected = single(expected);
        end

        function [eventDatFile, eventData, Fs] = buildEventSplitDatFixture(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            evObj = EventsManager(testCase.SaveFolder, '', 'csv');
            [frMat, ~, ~] = evObj.getFrameMatrix(size(rawData, 3), 'FrameRateHz', evObj.AcqInfo.FrameRateHz);

            testCase.assertNotEmpty(frMat, ...
                'The event frame matrix is empty.');

            validLens = sum(~isnan(frMat), 2);
            trialLen = min(validLens);
            testCase.assertGreaterThan(trialLen, 0);

            nTrials = size(frMat, 1);
            nY = size(rawData, 1);
            nX = size(rawData, 2);

            % Store event data on disk in Y,X,T,E order
            eventData = zeros(nY, nX, trialLen, nTrials, 'single');

            for iTrial = 1:nTrials
                frameIdx = frMat(iTrial, 1:trialLen);
                eventData(:,:,:,iTrial) = rawData(:,:,frameIdx);
            end

            Fs = testCase.getFrameRateHz();

            % Event-split data are written as a headered .dat (the layout
            % event-split outputs use); legacy sidecars never described
            % Y,X,T,E with a three-element datSize.
            hdr = struct('dataClass', 'single', 'frameRateHz', Fs, ...
                'exposureMsec', NaN, 'channelName', 'green_events', ...
                'dimNames', {{'Y','X','T','E'}}, ...
                'dimSizes', [nY, nX, trialLen, nTrials], 'writeComplete', true);

            eventDatFile = fullfile(testCase.SaveFolder, 'green_events.dat');
            fid = fopen(eventDatFile, 'w', 'ieee-le');
            testCase.assertNotEqual(fid, -1, 'Failed to create green_events.dat.');
            cleaner = onCleanup(@() fclose(fid));
            fwrite(fid, encodeDatHeader(hdr), 'uint8');
            fwrite(fid, eventData, 'single');
            clear cleaner
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