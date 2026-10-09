classdef TestGenAmplitudeMaps < matlab.unittest.TestCase
    %TESTGENAMPLITUDEMAPS Unit tests for genAmplitudeMaps.
    %
    % Contract: event-split Y-X-T-E .dat input only (no internal split), a
    % .dat Y-X-E output with one map per condition, streamed in X slabs.
    %
    % These tests validate:
    %   - the amplitude values against a direct in-RAM computation, for
    %     every baseline/response measure and a time window
    %   - ignored instances, aggregated (one-slice-per-condition) files, and
    %     a missing events.mat
    %   - forced multi-slab streaming against the single-slab result
    %   - pipelineInfo consistency
    %   - rejection of continuous Y-X-T data, arrays, UMT inputs, unmatched
    %     E axes, and invalid time windows

    properties
        ProjectRoot
        FixtureFolder
        TempFolder
        SourceInfo
        SplitFile
    end

    properties (TestParameter)
        responseMeasure = {'mean', 'median', 'min', 'max'}
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
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
                testCase.FixtureFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.FixtureFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', ...
                    'Analysis', ...
                    'TestingData_with_events');
            end

            testCase.verifyTrue(isfolder(testCase.FixtureFolder), ...
                'Fixture folder was not found.');

            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;

            acqInfoSrc = fullfile(testCase.FixtureFolder, 'AcqInfos.mat');
            testCase.verifyTrue(isfile(acqInfoSrc), ...
                'Fixture AcqInfos.mat was not found.');
            copyfile(acqInfoSrc, fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            sourceDatFile = iFindFirstDatFile(testCase.FixtureFolder);
            testCase.verifyTrue(~isempty(sourceDatFile) && isfile(sourceDatFile), ...
                'No .dat fixture file was found.');

            [dataYXT, info] = loadData(sourceDatFile);
            testCase.SourceInfo = info;

            saveData(fullfile(testCase.TempFolder, 'inputData.dat'), single(dataYXT), ...
                'DimNames', {'Y', 'X', 'T'}, 'Info', info);
            iCreateSyntheticEvents(testCase.TempFolder, size(dataYXT, 3), info.frameRateHz);

            % The event-split input: one slice per event instance.
            testCase.SplitFile = char(string(split_data_by_event( ...
                fullfile(testCase.TempFolder, 'inputData.dat'), testCase.TempFolder)));
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = genAmplitudeMaps('pipelineInfo');

            testCase.verifyEqual(info.name, 'genAmplitudeMaps');
            testCase.verifyFalse(isempty(info.inputs));
            testCase.verifyNumElements(info.outputs, 1);

            dataInputIdx = find(strcmp({info.inputs.name}, 'data'), 1, 'first');
            saveFolderIdx = find(strcmp({info.inputs.name}, 'SaveFolder'), 1, 'first');
            outIdx = find(strcmp({info.outputs.name}, 'outData'), 1, 'first');

            testCase.verifyNotEmpty(dataInputIdx);
            testCase.verifyNotEmpty(saveFolderIdx);
            testCase.verifyNotEmpty(outIdx);

            testCase.verifyTrue(info.inputs(dataInputIdx).isData);
            testCase.verifyTrue(info.inputs(dataInputIdx).supportsFile);
            testCase.verifyEqual(info.inputs(dataInputIdx).dataMode, 'file');
            testCase.verifyFalse(info.inputs(saveFolderIdx).isData);
            testCase.verifyEqual(info.outputs(outIdx).type, {'ProcessedData'});
            testCase.verifyEqual(info.outputs(outIdx).defOutfilename, 'amplitudeMap.dat');
        end

        function testEventSplitDatGivesOneMapPerCondition(testCase)
            out = genAmplitudeMaps(testCase.SplitFile, testCase.TempFolder);

            testCase.verifyEqual(out, fullfile(testCase.TempFolder, 'amplitudeMap.dat'));
            info = loadMetaData(out);
            testCase.verifyEqual(cellstr(string(info.dimNames(:).')), {'Y','X','E'});
            testCase.verifyEqual(info.dataClass, 'single');
            byEv = loadData(testCase.SplitFile);
            testCase.verifyEqual(double(info.dimSizes(:).'), ...
                [size(byEv, 1), size(byEv, 2), 2]);
        end

        function testMatchesDirectComputation(testCase, responseMeasure)
            out = genAmplitudeMaps(testCase.SplitFile, testCase.TempFolder, ...
                'BaselineMeasure', 'mean', 'ResponseMeasure', responseMeasure);

            expected = testCase.expectedMaps('mean', responseMeasure, []);
            testCase.verifyEqual(single(loadData(out)), expected, 'AbsTol', 1e-5);
        end

        function testTimeWindowStartingAtZeroIsAccepted(testCase)
            % 'allowed' advertises [0 Inf] and "0 to N seconds after onset" is
            % the natural request (P1-8).
            out = genAmplitudeMaps(testCase.SplitFile, testCase.TempFolder, ...
                'TimeWindow_sec', [0 0.1]);

            expected = testCase.expectedMaps('median', 'max', [0 0.1]);
            testCase.verifyEqual(single(loadData(out)), expected, 'AbsTol', 1e-5);
        end

        function testIgnoredEventIsExcluded(testCase)
            % .dat header Phase 8c: an ignored instance stays on the E axis of
            % the file but no longer contributes to its condition's map.
            full = loadData(genAmplitudeMaps(testCase.SplitFile, testCase.TempFolder, ...
                'ResponseMeasure', 'mean', 'BaselineMeasure', 'mean'));

            ev = EventsManager(testCase.TempFolder);
            inst = ev.getEventInstances();
            ev.removeRepetition(char(inst.eventName(1)), inst.repetitionIndex(1));
            ev.saveEvents(testCase.TempFolder);

            out = loadData(genAmplitudeMaps(testCase.SplitFile, testCase.TempFolder, ...
                'ResponseMeasure', 'mean', 'BaselineMeasure', 'mean'));
            testCase.verifySize(out, size(full));
            testCase.verifyNotEqual(out(:, :, 1), full(:, :, 1), ...
                'the ignored repetition no longer contributes');
            testCase.verifyEqual(single(out), ...
                testCase.expectedMaps('mean', 'mean', []), 'AbsTol', 1e-5);
        end

        function testAggregatedFileGivesOnePerSlice(testCase)
            % A file with one slice per condition (apply_aggregate_function
            % over E) is matched as aggregated: one map per slice.
            aggFile = apply_aggregate_function(testCase.SplitFile, testCase.TempFolder, ...
                'dimensionName', 'E', 'aggregateFcn', 'mean');
            aggData = single(loadData(aggFile));

            out = genAmplitudeMaps(aggFile, testCase.TempFolder, ...
                'BaselineMeasure', 'mean', 'ResponseMeasure', 'max');

            expected = zeros(size(aggData, 1), size(aggData, 2), size(aggData, 4), 'single');
            for k = 1:size(aggData, 4)
                resp = max(aggData(:, :, 8:end, k), [], 3, 'omitnan');
                base = mean(aggData(:, :, 1:7, k), 3, 'omitnan');
                expected(:, :, k) = resp - base;
            end
            testCase.verifyEqual(single(loadData(out)), expected, 'AbsTol', 1e-5);
        end

        function testWithoutEventsFileUsesOneCondition(testCase)
            delete(fullfile(testCase.TempFolder, 'events.mat'));

            testCase.verifyWarning(@() genAmplitudeMaps(testCase.SplitFile, ...
                testCase.TempFolder), 'Umitoolbox:genAmplitudeMaps:eventsNotMatched');
            out = loadData(fullfile(testCase.TempFolder, 'amplitudeMap.dat'));
            testCase.verifyEqual(size(out, 3), 1);
        end

        function testForcedMultiSlabMatchesSingleSlab(testCase)
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            single1 = loadData(genAmplitudeMaps(testCase.SplitFile, testCase.TempFolder));

            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            out = '';
            progress = evalc(['out = genAmplitudeMaps(testCase.SplitFile, ' ...
                'testCase.TempFolder);']);
            testCase.verifySubstring(progress, 'Chunk 2/');
            testCase.verifyEqual(loadData(out), single1, 'AbsTol', 1e-5);
        end

        function testInvalidTimeWindowStartGreaterThanEnd(testCase)
            testCase.verifyError(@() genAmplitudeMaps( ...
                testCase.SplitFile, testCase.TempFolder, 'TimeWindow_sec', [2 1]), ...
                'Umitoolbox:genAmplitudeMaps:invalidTimeWindow');
        end

        function testInvalidTimeWindowUnknownString(testCase)
            testCase.verifyError(@() genAmplitudeMaps( ...
                testCase.SplitFile, testCase.TempFolder, 'TimeWindow_sec', 'badString'), ...
                'Umitoolbox:genAmplitudeMaps:invalidTimeWindow');
        end

        function testInvalidTimeWindowBeyondTrialDuration(testCase)
            testCase.verifyError(@() genAmplitudeMaps( ...
                testCase.SplitFile, testCase.TempFolder, 'TimeWindow_sec', [0 1e6]), ...
                'Umitoolbox:genAmplitudeMaps:invalidTimeWindow');
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            % Continuous Y-X-T data must be split first.
            testCase.verifyError(@() genAmplitudeMaps( ...
                fullfile(sv, 'inputData.dat'), sv), ...
                'Umitoolbox:genAmplitudeMaps:unsupportedLayout');

            % Arrays and UMT structs.
            testCase.verifyError(@() genAmplitudeMaps(rand(6, 5, 20, 3, 'single'), sv), ...
                'Umitoolbox:genAmplitudeMaps:unsupportedInput');
            umt = genUMTStruct(rand(6, 5, 20, 'single'), 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y','X','T'});
            testCase.verifyError(@() genAmplitudeMaps(umt, sv), ...
                'Umitoolbox:genAmplitudeMaps:unsupportedInput');

            % .umt file and missing file.
            umtFile = fullfile(sv, 'x.umt');
            saveData(umtFile, umt);
            testCase.verifyError(@() genAmplitudeMaps(umtFile, sv), ...
                'Umitoolbox:genAmplitudeMaps:unsupportedExtension');
            testCase.verifyError(@() genAmplitudeMaps('missing.dat', sv), ...
                'Umitoolbox:genAmplitudeMaps:inputFileNotFound');

            % An E axis that matches neither instances (4) nor conditions (2).
            bad = fullfile(sv, 'bad.dat');
            writeTestDat(bad, rand(6, 5, 20, 3, 'single'), 10, 'DimNames', {'Y','X','T','E'});
            testCase.verifyError(@() genAmplitudeMaps(bad, sv), ...
                'Umitoolbox:genAmplitudeMaps:eventMappingMismatch');
        end
    end

    methods (Access = private)
        function expected = expectedMaps(testCase, baselineMeasure, responseMeasure, timeWindowSec)
            %EXPECTEDMAPS Direct in-RAM amplitude maps of the event-split file.

            byEv = single(loadData(testCase.SplitFile));
            rate = double(testCase.SourceInfo.frameRateHz);
            ev = EventsManager(testCase.TempFolder);
            plan = EventsManager.conditionAggregationPlan( ...
                ev.exportEventInfo('FrameRateHz', rate, 'IncludeIgnored', true));

            baselineFrames = 1:round(double(ev.baselinePeriod) * rate);
            if isempty(timeWindowSec)
                responseFrames = (baselineFrames(end) + 1):size(byEv, 3);
            else
                responseFrames = (baselineFrames(end) + 1 + round(timeWindowSec(1) * rate)): ...
                    (baselineFrames(end) + 1 + round(timeWindowSec(2) * rate));
            end

            reducer = @(x) iAgg(reshape(x(:, :, responseFrames, :), size(x, 1), size(x, 2), []), ...
                responseMeasure) - iAgg(reshape(x(:, :, baselineFrames, :), ...
                size(x, 1), size(x, 2), []), baselineMeasure);
            perCondition = EventsManager.reduceByCondition(byEv, plan, reducer, 4);
            expected = reshape(single(perCondition), size(perCondition, 1), size(perCondition, 2), []);
        end
    end
end

function out = iAgg(vals, name)
switch name
    case 'mean'
        out = mean(vals, 3, 'omitnan');
    case 'median'
        out = median(vals, 3, 'omitnan');
    case 'max'
        out = max(vals, [], 3, 'omitnan');
    case 'min'
        out = min(vals, [], 3, 'omitnan');
end
end

function datFile = iFindFirstDatFile(rootFolder)
files = dir(fullfile(rootFolder, '**', '*.dat'));
if isempty(files)
    datFile = '';
else
    datFile = fullfile(files(1).folder, files(1).name);
end
end

function iCreateSyntheticEvents(saveFolder, datLen, frameRateHz)
% Create a synthetic events.mat file aligned with the current EventsManager
% loading conventions.

baselinePeriod = single(7 / frameRateHz);

onsets = [15 45 75 105];
dur = 10;
offsets = min(onsets + dur - 1, datLen);

timestamps = reshape([onsets; offsets], [], 1);
timestamps = single((timestamps - 1) ./ frameRateHz);

state = logical(repmat([1; 0], numel(onsets), 1));
eventID = uint16([1;1;2;2;1;1;2;2]);
eventNameList = {'A','B'};
selectedEvents = true(size(eventID));

save(fullfile(saveFolder, 'events.mat'), ...
    'timestamps', 'state', 'eventID', 'eventNameList', ...
    'baselinePeriod', 'selectedEvents', '-mat');
end
