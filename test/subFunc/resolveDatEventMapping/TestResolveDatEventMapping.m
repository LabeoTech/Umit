classdef TestResolveDatEventMapping < matlab.unittest.TestCase
    %TESTRESOLVEDATEVENTMAPPING E axis of an event-split .dat onto events.mat.
    %
    %   .dat header Phase 8b: the mapping is used only when the E size equals
    %   the number of event instances in events.mat (ignored ones included);
    %   otherwise every slice is a repetition of one condition.

    properties
        Folder char
    end

    methods (TestMethodSetup)
        function createFolderWithEvents(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            writeMinimalAcqInfosMat(testCase.Folder, 'FrameRateHz', 10, 'AISampleRate', 100);
            ev = EventsManager(testCase.Folder, '', 'csv');
            ev.EventFileParseMethod = 'csv';
            ev.getTriggersFromSignal(makePulseSignal(400, [50 150 250], 10, 'Amplitude', 5), 100, false);
            csvFile = writeCSVConditionFile(testCase.Folder, 'labels.csv', {'A', 'B', 'A'});
            testCase.assertTrue(ev.readConditionFile(csvFile));
            ev.removeRepetition('A', 2);
            ev.saveEvents(testCase.Folder);
        end
    end

    methods (Test)
        function matchingCountMapsEveryInstance(testCase)
            info = iWriteSplit(testCase.Folder, 3);
            m = resolveDatEventMapping(info, testCase.Folder);
            testCase.verifyEqual(m.status, 'matched');
            testCase.verifyEqual(m.nE, 3);
            testCase.verifyEqual(m.nEvents, 3);
            testCase.verifyEqual(m.eventInfo.eventID, [1; 2; 1]);
            testCase.verifyEqual(m.eventInfo.repetitionIndex, [1; 1; 2]);
            testCase.verifyEqual(m.eventInfo.eventName, ["A"; "B"; "A"]);
            testCase.verifyEqual(m.eventInfo.selected, logical([1; 1; 0]));
            testCase.verifyEqual(m.eventInfo.durationSec, [0.1; 0.1; 0.1], 'AbsTol', 1e-6);
            testCase.verifyTrue(isfield(m.eventInfo, 'baselinePeriod'));
            testCase.verifyEmpty(m.message);
        end

        function mismatchedCountFallsBackToOneCondition(testCase)
            % 4 slices: neither 3 instances nor 2 conditions.
            info = iWriteSplit(testCase.Folder, 4);
            m = resolveDatEventMapping(info, testCase.Folder);
            testCase.verifyEqual(m.status, 'mismatch');
            testCase.verifyEqual(m.eventInfo.eventID, ones(4, 1));
            testCase.verifyEqual(m.eventInfo.repetitionIndex, (1:4).');
            testCase.verifyEqual(m.eventInfo.selected, true(4, 1));
            testCase.verifyTrue(all(isnan(m.eventInfo.durationSec)));
            testCase.verifySubstring(m.message, '4 event slices');
            testCase.verifySubstring(m.message, '3 event instances');
        end

        function conditionCountMapsAggregatedSlices(testCase)
            % .dat header Phase 8c: E = number of conditions (2 < 3
            % instances) is an aggregated output, one slice per condition in
            % first-appearance order; counts follow the current flags.
            info = iWriteSplit(testCase.Folder, 2);
            m = resolveDatEventMapping(info, testCase.Folder);
            testCase.verifyEqual(m.status, 'aggregated');
            testCase.verifyEqual(m.eventInfo.eventAxisMode, 'aggregated_repetitions');
            testCase.verifyEqual(double(m.eventInfo.eventID), [1; 2]);
            testCase.verifyEqual(m.eventInfo.eventName, ["A"; "B"]);
            testCase.verifyEqual(m.eventInfo.repetitionIndex, [0; 0]);
            testCase.verifyEqual(m.eventInfo.nInstances, [1; 1], 'A has one selected repetition');
            testCase.verifyEqual(m.eventInfo.selected, true(2, 1));
            testCase.verifyEqual(m.eventInfo.durationSec, [0.1; 0.1], 'AbsTol', 1e-6);
        end

        function oneRepetitionPerConditionReadsPerInstance(testCase)
            % Equal counts: per-instance and aggregated data are the same;
            % the file is read per instance.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            writeMinimalAcqInfosMat(folder, 'FrameRateHz', 10, 'AISampleRate', 100);
            ev = EventsManager(folder, '', 'csv');
            ev.EventFileParseMethod = 'csv';
            ev.getTriggersFromSignal(makePulseSignal(400, [50 250], 10, 'Amplitude', 5), 100, false);
            testCase.assertTrue(ev.readConditionFile(writeCSVConditionFile(folder, 'l.csv', {'B', 'A'})));
            ev.saveEvents(folder);

            m = resolveDatEventMapping(iWriteSplit(folder, 2), folder);
            testCase.verifyEqual(m.status, 'matched');
            testCase.verifyEqual(m.eventInfo.eventName, ["B"; "A"]);
        end

        function noEventsFileFallsBackToOneCondition(testCase)
            delete(fullfile(testCase.Folder, 'events.mat'));
            info = iWriteSplit(testCase.Folder, 3);
            m = resolveDatEventMapping(info, testCase.Folder);
            testCase.verifyEqual(m.status, 'noEvents');
            testCase.verifyEqual(m.eventInfo.repetitionIndex, (1:3).');
            testCase.verifyFalse(isfield(m.eventInfo, 'baselinePeriod'));
        end

        function continuousFileIsRejected(testCase)
            f = fullfile(testCase.Folder, 'cont.dat');
            saveData(f, zeros(2, 2, 3, 'single'), 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);
            testCase.verifyError(@() resolveDatEventMapping(loadMetaData(f), testCase.Folder), ...
                'Umitoolbox:resolveDatEventMapping:noEventAxis');
        end
    end
end

function info = iWriteSplit(folder, nE)
f = fullfile(folder, sprintf('byEv%d.dat', nE));
saveData(f, zeros(2, 2, 4, nE, 'single'), 'DimNames', {'Y', 'X', 'T', 'E'}, 'FrameRateHz', 10);
info = loadMetaData(f);
end
