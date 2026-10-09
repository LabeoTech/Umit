classdef TestNormalizeBSLN < matlab.unittest.TestCase
    %TESTNORMALIZEBSLN Unit tests for normalizeBSLN.
    %
    % Contract: .dat input only (Y-X-T or event-split Y-X-T-E), a .dat
    % output with the same axes, streamed in X slabs.
    %
    % Coverage:
    %   1) Y-X-T recording normalization (auto and numeric baseline,
    %      centered at one) against values computed outside the function
    %   2) Recording-mode clamping for too-long and tiny baseline windows
    %   3) Y-X-T-E files: per-trial recording and trial (events.mat baseline)
    %      normalization, ignored instances kept on the E axis
    %   4) Forced multi-slab streaming against the single-slab result
    %   5) Rejection of arrays, UMT inputs, other files, unsupported layouts,
    %      trial mode on continuous data, and a missing baseline period

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
            if isappdata(0, 'NormalizeBSLNTestConfig')
                cfg = getappdata(0, 'NormalizeBSLNTestConfig');
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

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'green.dat'), ...
                fullfile(testCase.SaveFolder, 'green.dat'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.SaveFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));

            testCase.createFreshEvents(0.8);
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = normalizeBSLN('pipelineInfo');

            testCase.verifyEqual(info.name, 'normalizeBSLN');
            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');
            testCase.verifyNumElements(info.outputs, 1);
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'normBSLN.dat');
        end

        function testDatRecordingAuto(testCase)
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(rawFile));
            freqHz = testCase.getFrameRateHz();
            expected = testCase.expectedRecording(rawData, freqHz, 'auto', false);

            out = normalizeBSLN(rawFile, testCase.SaveFolder, ...
                'normalizationMode', 'recording', ...
                'baselineMode', 'auto', ...
                'b_centerAtOne', false);

            testCase.verifyEqual(out, fullfile(testCase.SaveFolder, 'normBSLN.dat'));
            testCase.verifyEqual(iAxes(out), {'Y','X','T'});
            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testDatRecordingNumeric(testCase)
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(rawFile));
            freqHz = testCase.getFrameRateHz();
            baselineSec = 1.2;
            expected = testCase.expectedRecording(rawData, freqHz, baselineSec, true);

            out = normalizeBSLN(rawFile, testCase.SaveFolder, ...
                'normalizationMode', 'recording', ...
                'baselineMode', baselineSec, ...
                'b_centerAtOne', true);

            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testOutputKeepsTheInputRate(testCase)
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');

            out = normalizeBSLN(rawFile, testCase.SaveFolder);

            testCase.verifyEqual(loadMetaData(out).frameRateHz, loadMetaData(rawFile).frameRateHz);
            testCase.verifyEqual(loadMetaData(out).dataClass, 'single');
        end

        function testRecordingBaselineTooLongClamps(testCase)
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(rawFile));
            freqHz = testCase.getFrameRateHz();

            baselineSec = 1e6;
            expected = testCase.expectedRecording(rawData, freqHz, baselineSec, false);

            out = normalizeBSLN(rawFile, testCase.SaveFolder, ...
                'normalizationMode', 'recording', ...
                'baselineMode', baselineSec);

            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testRecordingTinyBaselineClampsToOneFrame(testCase)
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(rawFile));
            freqHz = testCase.getFrameRateHz();

            baselineSec = 1e-12;
            expected = testCase.expectedRecording(rawData, freqHz, baselineSec, false);

            out = normalizeBSLN(rawFile, testCase.SaveFolder, ...
                'normalizationMode', 'recording', ...
                'baselineMode', baselineSec);

            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testEventSplitRecordingNormalizesEveryTrial(testCase)
            % Y-X-T-E, recording mode: the baseline of every trial is its
            % own first 20% of T.
            [splitFile, byEv] = testCase.writeEventSplitFile();
            nT = size(byEv, 3);

            out = normalizeBSLN(splitFile, testCase.SaveFolder, ...
                'normalizationMode', 'recording', 'b_centerAtOne', true);

            testCase.verifyEqual(iAxes(out), {'Y','X','T','E'});
            expected = testCase.expectedPerTrial(byEv, max(1, round(0.2 * nT)), true);
            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testEventSplitTrialModeUsesTheEventsBaseline(testCase, centerAtOne)
            [splitFile, byEv] = testCase.writeEventSplitFile();
            ev = EventsManager(testCase.SaveFolder);
            nBase = min(size(byEv, 3), ...
                max(1, round(double(ev.baselinePeriod) * testCase.getFrameRateHz())));

            out = normalizeBSLN(splitFile, testCase.SaveFolder, ...
                'normalizationMode', 'trial', 'baselineMode', 'auto', ...
                'b_centerAtOne', centerAtOne);

            testCase.verifyEqual(iAxes(out), {'Y','X','T','E'});
            expected = testCase.expectedPerTrial(byEv, nBase, centerAtOne);
            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testTrialKeepsIgnoredInstancesFlagged(testCase)
            % .dat header Phase 8c: trial normalization keeps every instance
            % on the E axis, ignored ones included and flagged by the mapping.
            [splitFile, byEv] = testCase.writeEventSplitFile();
            ev = EventsManager(testCase.SaveFolder);
            inst = ev.getEventInstances();
            ev.removeRepetition(char(inst.eventName(1)), inst.repetitionIndex(1));
            ev.saveEvents(testCase.SaveFolder);

            out = normalizeBSLN(splitFile, testCase.SaveFolder, ...
                'normalizationMode', 'trial');

            testCase.verifyEqual(size(loadData(out), 4), size(byEv, 4));
            m = resolveDatEventMapping(loadMetaData(out), testCase.SaveFolder);
            testCase.verifyEqual(m.status, 'matched');
            expectedSelected = inst.selected;
            expectedSelected(1) = false;
            testCase.verifyEqual(m.eventInfo.selected, expectedSelected);
        end

        function testRerunCanOverwriteItsOwnInput(testCase)
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');
            first = normalizeBSLN(rawFile, testCase.SaveFolder);

            second = normalizeBSLN(first, testCase.SaveFolder);

            testCase.verifyEqual(second, first);
            testCase.verifyEmpty(dir(fullfile(testCase.SaveFolder, '*_writing.dat')));
        end

        function testForcedMultiSlabMatchesSingleSlab(testCase)
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            rawFile = fullfile(testCase.SaveFolder, 'green.dat');
            [splitFile, ~] = testCase.writeEventSplitFile();
            outYXT1 = single(loadData(normalizeBSLN(rawFile, testCase.SaveFolder)));
            outYXTE1 = single(loadData(normalizeBSLN(splitFile, testCase.SaveFolder, ...
                'normalizationMode', 'trial')));

            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            out = '';
            progress = evalc("out = normalizeBSLN(rawFile, testCase.SaveFolder);");
            testCase.verifySubstring(progress, 'Chunk 2/');
            testCase.verifyEqual(single(loadData(out)), outYXT1, 'AbsTol', 1e-6);

            outE = '';
            progress = evalc(['outE = normalizeBSLN(splitFile, testCase.SaveFolder, ' ...
                '''normalizationMode'', ''trial'');']);
            testCase.verifySubstring(progress, 'Chunk 2/');
            testCase.verifyEqual(single(loadData(outE)), outYXTE1, 'AbsTol', 1e-6);
        end

        function testTrialModeOnContinuousDataIsRejected(testCase)
            testCase.verifyError(@() normalizeBSLN( ...
                fullfile(testCase.SaveFolder, 'green.dat'), testCase.SaveFolder, ...
                'normalizationMode', 'trial'), ...
                'normalizeBSLN:TrialModeRequiresEventSplit');
        end

        function testTrialModeWithoutBaselinePeriodIsRejected(testCase)
            [splitFile, ~] = testCase.writeEventSplitFile();
            delete(fullfile(testCase.SaveFolder, 'events.mat'));

            testCase.verifyError(@() normalizeBSLN(splitFile, testCase.SaveFolder, ...
                'normalizationMode', 'trial'), ...
                'normalizeBSLN:MissingBaselinePeriod');
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.SaveFolder;

            % Arrays and UMT structs.
            testCase.verifyError(@() normalizeBSLN(rand(6, 5, 20, 'single'), sv), ...
                'normalizeBSLN:UnsupportedInputType');
            umt = genUMTStruct(rand(6, 5, 20, 'single'), 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y','X','T'});
            testCase.verifyError(@() normalizeBSLN(umt, sv), ...
                'normalizeBSLN:UnsupportedInputType');

            % .umt file, missing file.
            umtFile = fullfile(sv, 'x.umt');
            saveData(umtFile, umt);
            testCase.verifyError(@() normalizeBSLN(umtFile, sv), ...
                'normalizeBSLN:UnsupportedInputFile');
            testCase.verifyError(@() normalizeBSLN('missing.dat', sv), ...
                'normalizeBSLN:InputFileNotFound');

            % Layouts without T.
            yxe = fullfile(sv, 'yxe.dat');
            writeTestDat(yxe, rand(6, 5, 3, 'single'), 10, 'DimNames', {'Y','X','E'});
            testCase.verifyError(@() normalizeBSLN(yxe, sv), ...
                'Umitoolbox:normalizeBSLN:unsupportedLayout');
            yx = fullfile(sv, 'yx.dat');
            writeTestDat(yx, rand(6, 5, 'single'), 10, 'DimNames', {'Y','X'});
            testCase.verifyError(@() normalizeBSLN(yx, sv), ...
                'Umitoolbox:normalizeBSLN:unsupportedLayout');
        end

        function testTrialBaselineTooLongRejectedByEventsManager(testCase)
            tooLongSec = 1e6;

            try
                testCase.createFreshEvents(tooLongSec);
                testCase.verifyFail(['Expected EventsManager.setBaselinePeriod to ' ...
                    'reject a baseline period longer than the allowed inter-trigger window.']);
            catch ME
                testCase.verifyNotEmpty(ME.message);
                testCase.verifyTrue(contains(ME.message, ...
                    'Baseline time period is too long'), ...
                    ['Unexpected error message: ' ME.message]);
            end
        end

        function testTrialTinyBaselinePeriodIsAccepted(testCase)
            % .dat header Phase 7a: EventsManager no longer bounds the
            % baseline by one frame period from AcqInfos.mat; a tiny positive
            % baseline is accepted and normalization uses at least one frame.
            tinyBaselineSec = 1e-12;
            testCase.createFreshEvents(tinyBaselineSec);
            ev = EventsManager(testCase.SaveFolder);
            testCase.verifyEqual(double(ev.baselinePeriod), tinyBaselineSec, 'RelTol', 1e-5);
        end

    end

    properties (TestParameter)
        centerAtOne = {false, true}
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

        function [splitFile, byEv] = writeEventSplitFile(testCase)
            %WRITEEVENTSPLITFILE Event-split Y-X-T-E .dat of green.dat and its values.
            splitFile = char(string(split_data_by_event( ...
                fullfile(testCase.SaveFolder, 'green.dat'), testCase.SaveFolder)));
            byEv = single(loadData(splitFile));
        end

        function freqHz = getFrameRateHz(testCase)
            meta = loadMetaData(fullfile(testCase.SaveFolder, 'green.dat'));
            freqHz = double(meta.frameRateHz);
        end

        function expected = expectedRecording(testCase, dataIn, freqHz, baselineMode, b_centerAtOne)
            nT = size(dataIn, 3);
            nBaseFrames = testCase.resolveRecordingBaselineFrames(nT, freqHz, baselineMode);

            bsln = median(single(dataIn(:,:,1:nBaseFrames)), 3, 'omitnan');
            bsln(bsln == 0) = 1;
            expected = (single(dataIn) - bsln) ./ bsln;

            if b_centerAtOne
                expected = expected + 1;
            end
        end

        function expected = expectedPerTrial(~, byEv, nBaseFrames, b_centerAtOne)
            %EXPECTEDPERTRIAL Every trial normalized by the median of its first frames.
            expected = single(byEv);
            for iTrial = 1:size(byEv, 4)
                bsln = median(expected(:,:,1:nBaseFrames,iTrial), 3, 'omitnan');
                bsln(bsln == 0) = 1;
                expected(:,:,:,iTrial) = (expected(:,:,:,iTrial) - bsln) ./ bsln;
                if b_centerAtOne
                    expected(:,:,:,iTrial) = expected(:,:,:,iTrial) + 1;
                end
            end
        end

        function nBaseFrames = resolveRecordingBaselineFrames(~, nT, freqHz, baselineMode)
            if ischar(baselineMode) || (isstring(baselineMode) && isscalar(baselineMode))
                nBaseFrames = round(0.2 * nT);
            else
                nBaseFrames = round(double(baselineMode) * freqHz);
            end

            nBaseFrames = max(1, nBaseFrames);
            nBaseFrames = min(nBaseFrames, nT);
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

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end
    end
end

function names = iAxes(datFile)
%IAXES Axis names of a .dat file, as a row cell.
names = cellstr(string(loadMetaData(datFile).dimNames(:).'));
end

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
