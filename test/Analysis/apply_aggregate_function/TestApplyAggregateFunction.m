classdef TestApplyAggregateFunction < matlab.unittest.TestCase
    %TESTAPPLY_AGGREGATE_FUNCTION Unit tests for apply_aggregate_function.
    %
    % Coverage:
    %   1) .dat files (Y-X-T, Y-X-T-E, Y-X-E) in, .dat files out
    %   2) UMT structs and .umt files (single entry) in, UMT out
    %   3) All supported aggregations for both T and E
    %   4) Streamed (forced multi-slab) vs single-slab results
    %   5) Rejection of arrays, other files, missing axes, multi-entry UMTs
    %
    % Fixture policy:
    %   - green.dat and AcqInfos.mat are copied into a fresh temporary
    %     SaveFolder before each test.
    %   - Any pre-existing events.mat is deleted.
    %   - A fresh events.mat is then created using EventsManager with an
    %     external synthetic trigger signal.
    %
    % Notes:
    %   - A .dat input returns the path of a .dat output; a UMT input
    %     returns a valid UMT struct.
    %   - Numeric comparisons first try exact equality. If exact equality
    %     fails, the fallback criterion is:
    %         std(diff(:), 'omitnan') <= 1e-4

    properties
        SaveFolder char
        ProjectRoot char
        SampleDataFolder char
    end

    properties (TestParameter)
        aggregateFcn = {'mean', 'median', 'std', 'max', 'min', 'sum'}
    end

    methods (TestClassSetup)
        function addProjectPathAndLocateFixture(testCase)
            %ADDPROJECTPATHANDLOCATEFIXTURE Add project path and locate sample data.

            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);

            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);

            % Add the project recursively so analysis functions, classes,
            % and helper utilities are all available during the tests.
            addpath(genpath(testCase.ProjectRoot));

            cfg = [];
            if isappdata(0, 'ApplyAggregateFunctionTestConfig')
                cfg = getappdata(0, 'ApplyAggregateFunctionTestConfig');
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
            %CREATEFRESHWORKSPACE Create a fresh temporary SaveFolder.

            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fx.Folder;

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'green.dat'), ...
                fullfile(testCase.SaveFolder, 'green.dat'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.SaveFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));

            % Create a fresh events.mat for this test method.
            testCase.createFreshEvents();
        end
    end

    methods (TestMethodTeardown)
        function cleanupCreatedFiles(testCase)
            %CLEANUPCREATEDFILES Remove files explicitly created by the tests.

            deleteIfExists(fullfile(testCase.SaveFolder, 'events.mat'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'fixture_input.umt'));
            deleteIfExists(fullfile(testCase.SaveFolder, 'fixture_input.mat'));
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = apply_aggregate_function('pipelineInfo');

            testCase.verifyEqual(info.name, 'apply_aggregate_function');
            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');
            testCase.verifyNumElements(info.outputs, 1);
        end

        function testDatInputAggregateT(testCase, aggregateFcn)
            %TESTDATINPUTAGGREGATET Y-X-T .dat in -> Y-X .dat out.

            dataFile = fullfile(testCase.SaveFolder, 'green.dat');

            out = apply_aggregate_function( ...
                dataFile, ...
                testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, ...
                'dimensionName', 'T');

            testCase.verifyEqual(out, fullfile(testCase.SaveFolder, 'aggFcn_applied.dat'));
            testCase.verifyEqual(iAxes(out), {'Y','X'});

            rawData = single(loadData(dataFile));
            expected = testCase.buildExpectedRawTAggregation(rawData, aggregateFcn);
            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testDatInputAggregateEFromContinuousFile(testCase, aggregateFcn)
            %TESTDATINPUTAGGREGATEEFROMCONTINUOUSFILE Y-X-T .dat in -> Y-X-T-E .dat out.

            dataFile = fullfile(testCase.SaveFolder, 'green.dat');

            out = apply_aggregate_function( ...
                dataFile, ...
                testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, ...
                'dimensionName', 'E');

            testCase.verifyEqual(iAxes(out), {'Y','X','T','E'});
            rawData = single(loadData(dataFile));
            expected = testCase.buildExpectedRawEAggregation(rawData, aggregateFcn);
            testCase.verifyEqual(single(loadData(out)), expected, 'AbsTol', 1e-5);

            % The T axis is kept, so the output keeps the input's rate.
            testCase.verifyEqual(loadMetaData(out).frameRateHz, ...
                loadMetaData(dataFile).frameRateHz);
        end

        function testEventSplitDatAggregateT(testCase, aggregateFcn)
            %TESTEVENTSPLITDATAGGREGATET Y-X-T-E .dat in -> Y-X-E .dat out.

            [splitFile, byEv] = testCase.writeEventSplitFile();

            out = apply_aggregate_function(splitFile, testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, 'dimensionName', 'T');

            testCase.verifyEqual(iAxes(out), {'Y','X','E'});
            expected = testCase.aggregateAlongDim(byEv, 3, aggregateFcn);
            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testEventSplitDatAggregateE(testCase)
            %TESTEVENTSPLITDATAGGREGATEE Y-X-T-E .dat in -> Y-X-T-E .dat per condition.
            %   The file is matched to events.mat and reduced per condition,
            %   ignored instances excluded; an aggregated file is refused.

            sv = testCase.SaveFolder;
            [splitFile, byEv] = testCase.writeEventSplitFile();

            ev = EventsManager(sv);
            inst = ev.getEventInstances();
            ev.removeRepetition(char(inst.eventName(1)), 1);
            ev.saveEvents(sv);
            rateHz = loadMetaData(splitFile).frameRateHz;
            plan = EventsManager.conditionAggregationPlan( ...
                ev.exportEventInfo('FrameRateHz', rateHz));

            out = apply_aggregate_function(splitFile, sv, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E');

            expected = EventsManager.reduceByCondition(byEv, plan, ...
                @(x) mean(x, 4, 'omitnan'), 4);
            testCase.verifyEqual(iAxes(out), {'Y','X','T','E'});
            testCase.verifyEqual(single(loadData(out)), expected, 'AbsTol', 1e-5);
            testCase.verifySize(loadData(out), ...
                [size(byEv, 1:3), numel(plan.conditionID)]);

            % Direct aggregation of the continuous file gives the same values.
            direct = apply_aggregate_function(fullfile(sv, 'green.dat'), sv, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E');
            testCase.verifyEqual(single(loadData(direct)), single(loadData(out)), 'AbsTol', 1e-5);

            % An already aggregated file is refused.
            testCase.assumeLessThan(numel(plan.conditionID), numel(inst.eventID), ...
                'refusal needs fewer conditions than instances');
            testCase.verifyError(@() apply_aggregate_function(out, sv, ...
                'dimensionName', 'E'), 'Umitoolbox:EventsManager:alreadyAggregated');
        end

        function testYXEDatAggregateE(testCase, aggregateFcn)
            %TESTYXEDATAGGREGATEE Y-X-E .dat (one frame per instance) per condition.

            sv = testCase.SaveFolder;
            [splitFile, byEv] = testCase.writeEventSplitFile();
            perInstance = apply_aggregate_function(splitFile, sv, ...
                'aggregateFcn', 'mean', 'dimensionName', 'T');
            yxeFile = fullfile(sv, 'yxe.dat');
            movefile(perInstance, yxeFile);

            out = apply_aggregate_function(yxeFile, sv, ...
                'aggregateFcn', aggregateFcn, 'dimensionName', 'E');

            ev = EventsManager(sv);
            rateHz = loadMetaData(splitFile).frameRateHz;
            plan = EventsManager.conditionAggregationPlan( ...
                ev.exportEventInfo('FrameRateHz', rateHz, 'IncludeIgnored', true));
            perInstanceMean = testCase.aggregateAlongDim(byEv, 3, 'mean');
            expected = EventsManager.reduceByCondition(perInstanceMean, plan, ...
                @(x) reshape(testCase.aggregateAlongDim(x, 3, aggregateFcn), ...
                [size(x, 1), size(x, 2), 1]), 3);
            testCase.verifyEqual(iAxes(out), {'Y','X','E'});
            testCase.verifyNumericEquivalent(single(loadData(out)), expected);
        end

        function testIgnoredInstancesExcludedAndAllIgnoredIsNaN(testCase)
            % .dat header Phase 8c: an ignored repetition is left out of its
            % condition's aggregate; a condition whose instances are all
            % ignored is kept as a NaN slice; slices follow first appearance.
            dataFile = fullfile(testCase.SaveFolder, 'green.dat');
            data = single(loadData(dataFile));
            rate = testCase.acqFrameRateHz();
            ev = EventsManager(testCase.SaveFolder);
            inst = ev.getEventInstances();
            condIDs = EventsManager.conditionOrder(inst.eventID);
            nRep1 = nnz(inst.eventID == condIDs(1));
            testCase.assumeGreaterThanOrEqual(nRep1, 2, 'needs a condition with 2+ repetitions');

            full = loadData(apply_aggregate_function(dataFile, testCase.SaveFolder, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E'));

            ev.removeRepetition(char(inst.eventName(find(inst.eventID == condIDs(1), 1))), 1);
            if numel(condIDs) >= 2
                ev.removeCondition(char(inst.eventName(find(inst.eventID == condIDs(2), 1))));
            end
            ev.saveEvents(testCase.SaveFolder);

            out = loadData(apply_aggregate_function(dataFile, testCase.SaveFolder, ...
                'aggregateFcn', 'mean', 'dimensionName', 'E'));
            testCase.verifySize(out, size(full), 'every condition keeps its slice');

            % Condition 1 without its first repetition.
            [frMat, condList] = ev.getFrameMatrix(size(data, 3), 'FrameRateHz', rate, ...
                'IncludeIgnored', true);
            frMat = iCropLikeSplit(frMat);
            rows = find(double(condList) == condIDs(1));
            rows = rows(2:end);
            trials = nan(size(data, 1), size(data, 2), size(frMat, 2), numel(rows), 'single');
            for k = 1:numel(rows)
                valid = ~isnan(frMat(rows(k), :));
                trials(:, :, valid, k) = data(:, :, frMat(rows(k), valid));
            end
            testCase.verifyEqual(out(:, :, :, 1), mean(trials, 4, 'omitnan'), 'AbsTol', 1e-5);
            if numel(condIDs) >= 2
                testCase.verifyTrue(all(isnan(out(:, :, :, 2)), 'all'), 'all-ignored condition is NaN');
            end
        end

        function testForcedMultiChunkMatchesSingleChunk(testCase, aggregateFcn)
            %TESTFORCEDMULTICHUNKMATCHESSINGLECHUNK The memory mock forces many X
            %   slabs; the streamed result must equal the one-slab result for
            %   T aggregation and for E aggregation of a continuous file.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            dataFile = fullfile(testCase.SaveFolder, 'green.dat');
            outT1 = loadData(apply_aggregate_function(dataFile, testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, 'dimensionName', 'T'));
            outE1 = loadData(apply_aggregate_function(dataFile, testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, 'dimensionName', 'E'));

            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            outT = '';
            progressT = evalc(['outT = apply_aggregate_function(dataFile, testCase.SaveFolder, ' ...
                '''aggregateFcn'', aggregateFcn, ''dimensionName'', ''T'');']);
            testCase.verifySubstring(progressT, 'Chunk 2/');
            testCase.verifyEqual(loadData(outT), outT1, 'AbsTol', 1e-5);

            outE = '';
            progressE = evalc(['outE = apply_aggregate_function(dataFile, testCase.SaveFolder, ' ...
                '''aggregateFcn'', aggregateFcn, ''dimensionName'', ''E'');']);
            testCase.verifySubstring(progressE, 'Chunk 2/');
            testCase.verifyEqual(loadData(outE), outE1, 'AbsTol', 1e-5);
        end

        function testRerunOverwritesTheOutputWhenInputIsTheOutput(testCase)
            %TESTRERUNOVERWRITESTHEOUTPUTWHENINPUTISTHEOUTPUT The input can be
            %   the file that the output replaces (a pipeline re-run).
            [splitFile, byEv] = testCase.writeEventSplitFile();
            first = apply_aggregate_function(splitFile, testCase.SaveFolder, ...
                'dimensionName', 'T');
            second = apply_aggregate_function(first, testCase.SaveFolder, ...
                'dimensionName', 'E');

            testCase.verifyEqual(second, first);
            testCase.verifyEqual(iAxes(second), {'Y','X','E'});
            testCase.verifyNotEmpty(byEv);
            testCase.verifyEmpty(dir(fullfile(testCase.SaveFolder, '*_writing.dat')));
        end

        function testUMTInputAggregateT(testCase, aggregateFcn)
            %TESTUMTINPUTAGGREGATET Validate T aggregation for UMT input.

            umt = testCase.buildInputUMTFixture();

            out = apply_aggregate_function( ...
                umt, ...
                testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, ...
                'dimensionName', 'T');

            inVal = single(umt.data.main.value);
            testCase.verifyUMTOutput( ...
                out, ...
                {'main'}, ...
                {{'Y', 'X', 'E'}}, ...
                {testCase.aggregateAlongDim(inVal, 3, aggregateFcn)}, ...
                umt.eventInfo, ...
                true);
        end

        function testUMTInputAggregateE(testCase, aggregateFcn)
            %TESTUMTINPUTAGGREGATEE Validate E aggregation for UMT input.

            umt = testCase.buildInputUMTFixture();

            out = apply_aggregate_function( ...
                umt, ...
                testCase.SaveFolder, ...
                'aggregateFcn', aggregateFcn, ...
                'dimensionName', 'E');

            condIDs = unique(umt.eventInfo.eventID(:), 'stable');
            eventNames = cell(numel(condIDs), 1);

            for iCond = 1:numel(condIDs)
                idxFirst = find(umt.eventInfo.eventID(:) == condIDs(iCond), 1, 'first');
                eventNames{iCond} = umt.eventInfo.eventName{idxFirst};
            end

            inVal = single(umt.data.main.value);
            expectedData = testCase.aggregateGroupedDim( ...
                inVal, 4, umt.eventInfo.eventID(:), aggregateFcn);

            expectedEventInfo = struct();
            expectedEventInfo.eventID = condIDs(:);
            expectedEventInfo.repetitionIndex = zeros(numel(condIDs), 1);
            expectedEventInfo.eventName = eventNames;
            expectedEventInfo.eventAxisMode = 'aggregated_repetitions';

            testCase.verifyUMTOutput( ...
                out, ...
                {'main'}, ...
                {{'Y', 'X', 'T', 'E'}}, ...
                {expectedData}, ...
                expectedEventInfo, ...
                true);
        end

        function testUMTFileInputGivesTheSameUMT(testCase)
            %TESTUMTFILEINPUTGIVESTHESAMEUMT A .umt file gives the UMT of the struct.

            umt = testCase.buildInputUMTFixture();
            umtFile = fullfile(testCase.SaveFolder, 'fixture_input.umt');
            saveData(umtFile, umt);

            fromStruct = apply_aggregate_function(umt, testCase.SaveFolder, ...
                'dimensionName', 'E');
            ws = warning('off', 'apply_aggregate_function:UMTFileLoadsInRAM');
            restoreWarning = onCleanup(@() warning(ws));
            fromFile = apply_aggregate_function(umtFile, testCase.SaveFolder, ...
                'dimensionName', 'E');
            clear restoreWarning

            testCase.verifyUMTEquivalent(fromStruct, fromFile);
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.SaveFolder;

            % Arrays are not an input form.
            testCase.verifyError(@() apply_aggregate_function( ...
                rand(6, 5, 3, 'single'), sv), ...
                'apply_aggregate_function:UnsupportedInputType');

            % Other file types and missing files.
            matFile = fullfile(sv, 'x.mat');
            fclose(fopen(matFile, 'w'));
            testCase.verifyError(@() apply_aggregate_function(matFile, sv), ...
                'apply_aggregate_function:UnsupportedInputFile');
            testCase.verifyError(@() apply_aggregate_function('missing.dat', sv), ...
                'apply_aggregate_function:InputFileNotFound');

            % Layouts without T or E.
            yx = fullfile(sv, 'yx.dat');
            writeTestDat(yx, rand(6, 5, 'single'), 10, 'DimNames', {'Y','X'});
            testCase.verifyError(@() apply_aggregate_function(yx, sv), ...
                'Umitoolbox:apply_aggregate_function:unsupportedLayout');

            % The requested axis must exist: no T in Y-X-E.
            yxe = fullfile(sv, 'yxe_only.dat');
            writeTestDat(yxe, rand(6, 5, 3, 'single'), 10, 'DimNames', {'Y','X','E'});
            testCase.verifyError(@() apply_aggregate_function(yxe, sv, 'dimensionName', 'T'), ...
                'apply_aggregate_function:MissingDimension');

            % A UMT with several entries.
            umt = testCase.buildInputUMTFixture();
            umt = genUMTStruct(umt, 'value', 2 * umt.data.main.value, ...
                'entryName', 'scaled', 'dimNames', {'Y', 'X', 'T', 'E'});
            testCase.verifyError(@() apply_aggregate_function(umt, sv), ...
                'apply_aggregate_function:multipleCompatibleUMTEntries');
        end
    end

    methods (Access = private)

        function frameRateHz = acqFrameRateHz(testCase)
            %ACQFRAMERATEHZ Frame rate of the fixture dataset (from its AcqInfos.mat).

            acq = load(fullfile(testCase.SaveFolder, 'AcqInfos.mat'));
            if isfield(acq, 'AcqInfoStream')
                acqInfo = acq.AcqInfoStream;
            else
                fn = fieldnames(acq);
                acqInfo = acq.(fn{1});
            end
            frameRateHz = double(acqInfo.FrameRateHz);
        end

        function createFreshEvents(testCase)
            %CREATEFRESHEVENTS Create a fresh events.mat using EventsManager.

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

            % Synthetic trigger signal with several repetitions of one
            % generic condition. This is sufficient to test event-based
            % aggregation and keeps the event fixture self-contained.
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

            evObj = EventsManager(testCase.SaveFolder);
            evObj.getTriggersFromSignal(triggerSignal, signalSR, false);
            evObj.saveEvents(testCase.SaveFolder);

            testCase.verifyTrue(isfile(fullfile(testCase.SaveFolder, 'events.mat')));
        end

        function [splitFile, byEv] = writeEventSplitFile(testCase)
            %WRITEEVENTSPLITFILE Event-split Y-X-T-E .dat of green.dat and its values.

            sv = testCase.SaveFolder;
            splitFile = char(string(split_data_by_event(fullfile(sv, 'green.dat'), sv)));
            byEv = single(loadData(splitFile));
        end

        function umt = buildInputUMTFixture(testCase)
            %BUILDINPUTUMTFIXTURE Build one valid image-only UMT fixture from green.dat.

            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            evObj = EventsManager(testCase.SaveFolder);

            [frMat, conditionIDlist, repetitionList] = evObj.getFrameMatrix(size(rawData, 3), 'FrameRateHz', evObj.AcqInfo.FrameRateHz);

            testCase.assertNotEmpty(frMat, ...
                'The event frame matrix is empty. Failed to build the UMT fixture.');

            nTrial = size(frMat, 1);
            nY = size(rawData, 1);
            nX = size(rawData, 2);
            trialLen = size(frMat, 2);

            dataByEv = nan(nY, nX, trialLen, nTrial, 'single');

            for iTrial = 1:nTrial
                validMask = ~isnan(frMat(iTrial, :));
                if any(validMask)
                    frameIdx = frMat(iTrial, validMask);
                    dataByEv(:, :, validMask, iTrial) = rawData(:, :, frameIdx);
                end
            end

            % Match the splitDataByEvents behavior: crop all trials to the
            % shortest valid length if NaNs are present.
            if any(isnan(frMat(:)))
                firstNaNCol = find(any(isnan(frMat), 1), 1, 'first');
                if ~isempty(firstNaNCol)
                    dataByEv(:, :, firstNaNCol:end, :) = [];
                end
            end

            eventNames = evObj.eventNameList(conditionIDlist);

            umt = genUMTStruct( ...
                dataByEv, ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y', 'X', 'T', 'E'});

            umt = appendUMTEventInfo( ...
                umt, ...
                'eventID', conditionIDlist(:), ...
                'repetitionIndex', repetitionList(:), ...
                'eventName', eventNames(:), ...
                'eventAxisMode', 'instances');

            validateUMTStruct(umt, 'requireEventInfo', true);
        end

        function expected = buildExpectedRawTAggregation(testCase, data, aggregateFcn)
            %BUILDEXPECTEDRAWTAGGREGATION Build expected YX output for raw T aggregation.

            expected = testCase.aggregateAlongDim(single(data), 3, aggregateFcn);
        end

        function [expected, eventInfo] = buildExpectedRawEAggregation(testCase, data, aggregateFcn)
            %BUILDEXPECTEDRAWEAGGREGATION Build expected YXTE output for raw E aggregation.

            evObj = EventsManager(testCase.SaveFolder);
            [frMat, conditionIDlist, ~] = evObj.getFrameMatrix(size(data, 3), 'FrameRateHz', evObj.AcqInfo.FrameRateHz);

            testCase.assertNotEmpty(frMat, ...
                'The event frame matrix is empty.');
            frMat = iCropLikeSplit(frMat);

            condIDs = unique(conditionIDlist(:), 'stable');
            nCond = numel(condIDs);
            nY = size(data, 1);
            nX = size(data, 2);
            trialLen = size(frMat, 2);

            expected = zeros(nY, nX, trialLen, nCond, 'single');
            eventNames = cell(nCond, 1);

            for iCond = 1:nCond
                rowIdx = find(conditionIDlist(:) == condIDs(iCond));
                nRep = numel(rowIdx);

                condBlock = nan(nY, nX, trialLen, nRep, 'single');

                for iRep = 1:nRep
                    validMask = ~isnan(frMat(rowIdx(iRep), :));
                    if any(validMask)
                        frameIdx = frMat(rowIdx(iRep), validMask);
                        condBlock(:, :, validMask, iRep) = data(:, :, frameIdx);
                    end
                end

                expected(:, :, :, iCond) = testCase.aggregateAlongDim(condBlock, 4, aggregateFcn);
                eventNames{iCond} = evObj.eventNameList{condIDs(iCond)};
            end

            eventInfo = struct();
            eventInfo.eventID = condIDs(:);
            eventInfo.repetitionIndex = zeros(nCond, 1);
            eventInfo.eventName = eventNames;
            eventInfo.eventAxisMode = 'aggregated_repetitions';
            % Raw/.dat E outputs are .dat since Phase 8c; eventInfo is kept
            % for callers that still build a UMT expectation.
        end

        function outData = aggregateAlongDim(testCase, dataIn, dimIdx, aggregateFcn)
            %AGGREGATEALONGDIM Aggregate along one dimension and remove it.

            permOrder = [dimIdx, setdiff(1:ndims(dataIn), dimIdx)];
            dataP = permute(dataIn, permOrder);
            szP = size(dataP);
            dataP = reshape(dataP, szP(1), []);

            aggFlat = testCase.calcAgg(dataP, aggregateFcn);
            outData = reshape(single(aggFlat), szP(2:end));

            if isvector(outData) && ~isscalar(outData)
                outData = outData(:);
            end
        end

        function outData = aggregateGroupedDim(testCase, dataIn, dimIdx, groupID, aggregateFcn)
            %AGGREGATEGROUPEDDIM Aggregate one dimension by stable group IDs.

            groupID = groupID(:);
            uniqueID = unique(groupID, 'stable');

            permOrder = [dimIdx, setdiff(1:ndims(dataIn), dimIdx)];
            dataP = permute(dataIn, permOrder);
            szP = size(dataP);
            dataP = reshape(dataP, szP(1), []);

            outFlat = zeros(numel(uniqueID), size(dataP, 2), 'single');

            for iGroup = 1:numel(uniqueID)
                idx = groupID == uniqueID(iGroup);
                outFlat(iGroup, :) = single(testCase.calcAgg(dataP(idx, :), aggregateFcn));
            end

            outP = reshape(outFlat, [numel(uniqueID), szP(2:end)]);
            outData = ipermute(outP, permOrder);
        end

        function out = calcAgg(~, vals, aggregateFcn)
            %CALCAGG Apply one supported aggregation along the first dimension.

            switch lower(aggregateFcn)
                case 'mean'
                    out = mean(vals, 1, 'omitnan');

                case 'median'
                    out = median(vals, 1, 'omitnan');

                case 'std'
                    out = std(vals, 0, 1, 'omitnan');

                case 'max'
                    out = max(vals, [], 1, 'omitnan');

                case 'min'
                    out = min(vals, [], 1, 'omitnan');

                case 'sum'
                    out = sum(vals, 1, 'omitnan');

                otherwise
                    error('TestApplyAggregateFunction:UnsupportedAggregateFcn', ...
                        'Unsupported aggregate function "%s".', aggregateFcn);
            end
        end

        function verifyUMTEquivalent(testCase, a, b)
            %VERIFYUMTEQUIVALENT Compare two UMT outputs with tolerant numeric checks.

            validateUMTStruct(a, 'requireEventInfo', false);
            validateUMTStruct(b, 'requireEventInfo', false);

            testCase.verifyEqual(lower(char(string(a.kind))), lower(char(string(b.kind))));
            testCase.verifyEqual(a.version, b.version);

            aNames = fieldnames(a.data);
            bNames = fieldnames(b.data);
            testCase.verifyEqual(aNames, bNames);

            for iEntry = 1:numel(aNames)
                nm = aNames{iEntry};
                testCase.verifyEqual( ...
                    cellstr(string(a.data.(nm).dimNames)), ...
                    cellstr(string(b.data.(nm).dimNames)));

                testCase.verifyNumericEquivalent( ...
                    single(a.data.(nm).value), ...
                    single(b.data.(nm).value));
            end

            if isfield(a, 'eventInfo') || isfield(b, 'eventInfo')
                testCase.verifyTrue(isfield(a, 'eventInfo'));
                testCase.verifyTrue(isfield(b, 'eventInfo'));
                testCase.verifyEqual(a.eventInfo.eventID, b.eventInfo.eventID);
                testCase.verifyEqual(a.eventInfo.repetitionIndex, b.eventInfo.repetitionIndex);
                testCase.verifyEqual(a.eventInfo.eventName, b.eventInfo.eventName);
                testCase.verifyEqual( ...
                    char(string(a.eventInfo.eventAxisMode)), ...
                    char(string(b.eventInfo.eventAxisMode)));
            end
        end

        function verifyUMTOutput(testCase, umt, entryNames, entryDims, entryData, eventInfo, expectEventInfo)
            %VERIFYUMTOUTPUT Validate a UMT output against expected content.

            validateUMTStruct(umt, 'requireEventInfo', false);

            testCase.verifyEqual(lower(char(string(umt.kind))), 'image');

            outNames = fieldnames(umt.data);
            testCase.verifyEqual(outNames, entryNames(:));

            for iEntry = 1:numel(entryNames)
                nm = entryNames{iEntry};
                testCase.verifyEqual( ...
                    cellstr(string(umt.data.(nm).dimNames)), ...
                    entryDims{iEntry});
                testCase.verifyNumericEquivalent( ...
                    single(umt.data.(nm).value), ...
                    single(entryData{iEntry}));
            end

            if expectEventInfo
                testCase.verifyTrue(isfield(umt, 'eventInfo'));
                testCase.verifyEqual(umt.eventInfo.eventID, eventInfo.eventID);
                testCase.verifyEqual(umt.eventInfo.repetitionIndex, eventInfo.repetitionIndex);
                testCase.verifyEqual(umt.eventInfo.eventName, eventInfo.eventName);
                testCase.verifyEqual( ...
                    char(string(umt.eventInfo.eventAxisMode)), ...
                    char(string(eventInfo.eventAxisMode)));
            else
                testCase.verifyFalse(isfield(umt, 'eventInfo'));
            end
        end

        function verifyNumericEquivalent(testCase, a, b)
            %VERIFYNUMERICEQUIVALENT Compare numeric arrays with exact-or-small-diff rule.

            testCase.verifyEqual(size(a), size(b));

            if isequaln(a, b)
                return
            end

            diffVals = double(a(:)) - double(b(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                % Arrays differ only by NaN placement or both are effectively empty.
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
%DELETEIFEXISTS Delete file if it exists.

if isfile(filePath)
    delete(filePath);
end
end
function frMat = iCropLikeSplit(frMat)
%ICROPLIKESPLIT Trial length of EventsManager.splitDataByEvents: crop from
% the first frame column that any instance lacks.
firstNaNCol = find(any(isnan(frMat), 1), 1, 'first');
if ~isempty(firstNaNCol)
    frMat(:, firstNaNCol:end) = [];
end
end
