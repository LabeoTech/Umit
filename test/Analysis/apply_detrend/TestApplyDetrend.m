classdef TestApplyDetrend < matlab.unittest.TestCase
    %TESTAPPLYDETREND Unit tests for apply_detrend.
    %
    % Coverage:
    %   1) pipelineInfo generation
    %   2) Raw YXT array input
    %   3) Raw YXTE array input
    %   4) Raw .dat filename input (YXT and event-split YXTE)
    %   5) Standard vs low-RAM comparison, including forced multi-slab runs
    %   6) Rejection of UMT struct, .umt file, and unsupported layouts
    %
    % Ground truth:
    %   Uses the same detrend algorithm outside the function to compare
    %   implementations without changing the algorithm.

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
            if isappdata(0, 'ApplyDetrendTestConfig')
                cfg = getappdata(0, 'ApplyDetrendTestConfig');
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
            deleteIfExists(fullfile(testCase.SaveFolder, 'data_detrended_PREALLOCATED_FILE.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'input.umt'));

            testCase.createFreshEvents(0.8);
        end
    end

    methods (TestMethodTeardown)
        function cleanupCreatedFiles(testCase)
            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'data_detrended_PREALLOCATED_FILE.dat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'input.umt'));
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = apply_detrend('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));

            reqFields = {'name','description','version','inputs','outputs'};
            testCase.verifyTrue(all(ismember(reqFields, fieldnames(info))), ...
                'pipelineInfo is missing one or more required top-level fields.');

            testCase.verifyEqual(info.name, 'apply_detrend');
            testCase.verifyEqual(numel(info.inputs), 2);
            testCase.verifyEqual(numel(info.outputs), 1);

            inputNames = {info.inputs.name};
            outputNames = {info.outputs.name};

            testCase.verifyEqual(inputNames, {'data','SaveFolder'});
            testCase.verifyEqual(outputNames, {'outData'});

            dataInput = info.inputs(strcmp(inputNames, 'data'));
            outDataOutput = info.outputs(strcmp(outputNames, 'outData'));

            dataTypes = dataInput.type;
            outTypes = outDataOutput.type;

            if ischar(dataTypes)
                dataTypes = {dataTypes};
            end
            if ischar(outTypes)
                outTypes = {outTypes};
            end

            % 'ImageTimeSeriesByEvents' is not part of the current type
            % taxonomy (no Analysis function declares it; siblings such as
            % spatialGaussFilt.m use the same {Image,ImageTimeSeries,
            % ProcessedData,UnknownDataType} set) -- only assert the types
            % apply_detrend actually declares.
            testCase.verifyTrue(ismember('ImageTimeSeries', dataTypes));
            testCase.verifyTrue(ismember('ProcessedData', dataTypes));
            testCase.verifyTrue(ismember('ImageTimeSeries', outTypes));
            testCase.verifyTrue(ismember('ProcessedData', outTypes));
        end

        function testArrayInputYXT(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            expected = testCase.expectedDetrendYXT(rawData);

            out = apply_detrend(rawData, testCase.SaveFolder, ...
                'FrameRateHz', testCase.acqFrameRateHz());

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testArrayInputYXTE(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            eventData = testCase.buildSyntheticYXTE(rawData);
            expected = testCase.expectedDetrendYXTE(eventData);

            out = apply_detrend(eventData, testCase.SaveFolder, ...
                'FrameRateHz', testCase.acqFrameRateHz());

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(eventData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
        end

        function testDatInput(testCase)
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));
            expected = testCase.expectedDetrendYXT(rawData);

            out = apply_detrend(inFile, testCase.SaveFolder);

            testCase.verifyTrue(ischar(out) || (isstring(out) && isscalar(out)));
            outFile = char(string(out));
            testCase.verifyTrue(isfile(outFile));

            actual = single(loadData(outFile));
            testCase.verifyEqual(size(actual), size(rawData));
            testCase.verifyNumericEquivalent(actual, single(expected));
        end

        function testStandardVsLowRAM(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));

            outStandard = apply_detrend(rawData, testCase.SaveFolder, ...
                'FrameRateHz', testCase.acqFrameRateHz());
            outLowRAMFile = apply_detrend(fullfile(testCase.SaveFolder, 'green.dat'), ...
                testCase.SaveFolder);

            testCase.verifyTrue(ischar(outLowRAMFile) || (isstring(outLowRAMFile) && isscalar(outLowRAMFile)));
            testCase.verifyTrue(isfile(char(string(outLowRAMFile))));

            outLowRAM = single(loadData(char(string(outLowRAMFile))));

            testCase.verifyEqual(size(outStandard), size(outLowRAM));

            if isequaln(outStandard, outLowRAM)
                return
            end

            diffVals = double(outStandard(:)) - double(outLowRAM(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                testCase.verifyEqual(outStandard, outLowRAM);
                return
            end

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end

        function testDatInputYXTEKeepsTheAxesAndDetrendsEachTrial(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            eventData = testCase.buildSyntheticYXTE(rawData);
            inFile = testCase.writeEventSplitDat(eventData, 'byEvent.dat');

            outFile = apply_detrend(inFile, testCase.SaveFolder);

            outInfo = loadMetaData(char(string(outFile)));
            testCase.verifyEqual(outInfo.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(outInfo.dimSizes(:)).', size(eventData));
            testCase.verifyNumericEquivalent(single(loadData(char(string(outFile)))), ...
                single(testCase.expectedDetrendYXTE(eventData)));
        end

        function testStandardVsLowRAMYXTE(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            eventData = testCase.buildSyntheticYXTE(rawData);
            inFile = testCase.writeEventSplitDat(eventData, 'byEvent.dat');

            outStandard = apply_detrend(eventData, testCase.SaveFolder, ...
                'FrameRateHz', testCase.acqFrameRateHz());
            outFile = apply_detrend(inFile, testCase.SaveFolder);

            testCase.verifyNumericEquivalent(single(loadData(char(string(outFile)))), ...
                single(outStandard));
        end

        function testForcedMultiSlabMatchesInRamForYXTAndYXTE(testCase)
            % The memory mock makes calculateMaxChunkSize ask for many X
            % slabs: the streamed result must still equal the in-RAM one.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            eventData = testCase.buildSyntheticYXTE(rawData);
            byEventFile = testCase.writeEventSplitDat(eventData, 'byEvent.dat');
            rate = testCase.acqFrameRateHz();

            cases = struct( ...
                'file', {fullfile(testCase.SaveFolder, 'green.dat'), byEventFile}, ...
                'inRam', {apply_detrend(rawData, testCase.SaveFolder, 'FrameRateHz', rate), ...
                          apply_detrend(eventData, testCase.SaveFolder, 'FrameRateHz', rate)});
            for iCase = 1:numel(cases)
                progress = evalc('outFile = apply_detrend(cases(iCase).file, testCase.SaveFolder);');
                slabCounts = cellfun(@(c) str2double(c{1}), ...
                    regexp(progress, 'Chunk \d+/(\d+) \[Reading', 'tokens'));
                testCase.verifyGreaterThan(max(slabCounts), 1, ...
                    'The fixture must force more than one X slab.');
                testCase.verifyNumericEquivalent( ...
                    single(loadData(char(string(outFile)))), single(cases(iCase).inRam));
            end
        end

        function testRejectsUMTAndUnsupportedLayouts(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            umt = genUMTStruct(rawData, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.SaveFolder, 'input.umt');
            saveData(umtFile, umt);
            yxeFile = fullfile(testCase.SaveFolder, 'yxe.dat');
            saveData(yxeFile, rawData(:, :, 1:3), 'DimNames', {'Y','X','E'}, ...
                'FrameRateHz', testCase.acqFrameRateHz());
            rate = testCase.acqFrameRateHz();

            testCase.verifyError(@() apply_detrend(umt, testCase.SaveFolder, ...
                'FrameRateHz', rate), 'apply_detrend:UnsupportedInputType');
            testCase.verifyError(@() apply_detrend(umtFile, testCase.SaveFolder), ...
                'apply_detrend:UnsupportedInputFile');
            testCase.verifyError(@() apply_detrend(yxeFile, testCase.SaveFolder), ...
                'Umitoolbox:apply_detrend:unsupportedLayout');
            testCase.verifyError(@() apply_detrend(rawData(:, :, 1), testCase.SaveFolder, ...
                'FrameRateHz', rate), 'apply_detrend:InvalidArrayInput');
        end

        function testInRamInputWithoutFrameRateErrors(testCase)
            % events.mat defines a baseline period, so the in-RAM path needs
            % a frame rate; AcqInfos.mat in the folder is not used for it.
            testCase.assertTrue(isfile(fullfile(testCase.SaveFolder, 'AcqInfos.mat')));
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));

            testCase.verifyError(@() apply_detrend(rawData, testCase.SaveFolder), ...
                'Umitoolbox:apply_detrend:missingFrameRateHz');
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

        function out = expectedDetrendYXT(testCase, inData)
            frames = testCase.getDetrendFrameCount(size(inData, 3));

            slab = reshape(single(inData), [], size(inData, 3));
            slab = testCase.expectedDetrend2D(slab, frames);
            out = reshape(slab, size(inData));
        end

        function out = expectedDetrendYXTE(testCase, inData)
            out = zeros(size(inData), 'like', inData);
            frames = testCase.getDetrendFrameCount(size(inData, 3));

            for iEvent = 1:size(inData, 4)
                slab = reshape(single(inData(:,:,:,iEvent)), [], size(inData, 3));
                slab = testCase.expectedDetrend2D(slab, frames);
                out(:,:,:,iEvent) = reshape(slab, size(inData,1), size(inData,2), size(inData,3));
            end
        end

        function out2D = expectedDetrend2D(~, in2D, frames)
            Nt = size(in2D, 2);

            frames = min(frames, Nt);
            frames = max(frames, 3);

            if mod(frames, 2) == 0
                frames = frames - 1;
            end

            if frames >= Nt
                frames = max(3, Nt - 1);
                if mod(frames, 2) == 0
                    frames = frames - 1;
                end
            end

            if frames < 3 || Nt < 3
                out2D = in2D;
                return
            end

            delta_y = median(in2D(:, end-frames+1:end), 2, 'omitnan') - ...
                      median(in2D(:, 1:frames), 2, 'omitnan');

            delta_x = Nt - frames;
            M = delta_y ./ delta_x;
            b = median(in2D(:, 1:frames), 2, 'omitnan');

            trend = M .* linspace(-2, Nt-3, Nt) + b;
            out2D = in2D - trend + b;
        end

        function frames = getDetrendFrameCount(testCase, Nt)
            frames = 7;

            freqHz = [];
            meta = loadMetaData(fullfile(testCase.SaveFolder, 'green.dat'));
            if isfield(meta, 'frameRateHz') && isfinite(meta.frameRateHz)
                freqHz = double(meta.frameRateHz);
            end

            baselineSec = [];
            evObj = EventsManager(testCase.SaveFolder, '', 'csv');
            if ~isempty(evObj.baselinePeriod)
                baselineSec = double(evObj.baselinePeriod);
            end

            if ~isempty(freqHz) && ~isempty(baselineSec)
                frames = round(baselineSec * freqHz);
            end

            frames = max(frames, 3);

            if mod(frames, 2) == 0
                frames = frames + 1;
            end

            frames = min(frames, Nt);
            if mod(frames, 2) == 0 && frames > 3
                frames = frames - 1;
            end
            frames = max(min(frames, Nt), 3);
        end

        function eventData = buildSyntheticYXTE(~, rawData)
            nT = size(rawData, 3);
            trialLen = max(5, floor(nT / 4));
            nEvents = 3;

            eventData = zeros(size(rawData,1), size(rawData,2), trialLen, nEvents, 'single');

            for iEvent = 1:nEvents
                startIdx = 1 + (iEvent-1) * trialLen;
                stopIdx = min(startIdx + trialLen - 1, nT);
                thisLen = stopIdx - startIdx + 1;
                eventData(:,:,1:thisLen,iEvent) = rawData(:,:,startIdx:stopIdx);
            end
        end

        function inFile = writeEventSplitDat(testCase, eventData, fileName)
            % Save a Y-X-T-E array as a headered .dat at the fixture's rate.
            inFile = fullfile(testCase.SaveFolder, fileName);
            saveData(inFile, eventData, 'DimNames', {'Y','X','T','E'}, ...
                'FrameRateHz', testCase.acqFrameRateHz());
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

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end