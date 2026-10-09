classdef TestEventsManagerSaveLoad < matlab.unittest.TestCase
    %TESTEVENTSMANAGERSAVELOAD Round-trip tests for events.mat persistence.

    properties
        TempFolder char
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
            addpath(fullfile(fileparts(mfilename('fullpath')), 'EMTestHelpers'));
        end
    end

    methods (Test)
        function testSaveLoadRoundTripAnalog(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj1 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            obj1.removeRepetition('Main', 2);
            obj1.removeRepetition('Main', 3, 'Purge', true);
            obj1.saveEvents(testCase.TempFolder);

            obj2 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifyFalse(obj2.b_hasExternalSignal);
            testCase.verifyEqual(obj2.timestamps, obj1.timestamps);
            testCase.verifyEqual(obj2.state, obj1.state);
            testCase.verifyEqual(obj2.eventID, obj1.eventID);
            testCase.verifyEqual(obj2.eventNameList, obj1.eventNameList);
            testCase.verifyEqual(obj2.selectedEvents, obj1.selectedEvents);
            testCase.verifyEqual(obj2.PurgedEvents, obj1.PurgedEvents);
            % .dat header Phase 8b: the purged repetition was removed.
            testCase.verifyEqual(numel(obj2.eventID), 4);
            testCase.verifyFalse(any(obj2.PurgedEvents));
            testCase.verifyEqual(obj2.baselinePeriod, obj1.baselinePeriod, 'AbsTol', 1e-6);
            testCase.verifyEqual(obj2.trigChanName, obj1.trigChanName);
        end

        function testLegacyPurgeFlagsRemovedOnSave(testCase)
            % .dat header Phase 8b: events flagged by an older purge are
            % removed when events.mat is saved and no event-split file exists.
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], 'PulseWidth', 80, 'StimName', 'Main');
            obj1 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            obj1.saveEvents(testCase.TempFolder);
            iFlagLegacyPurge(testCase.TempFolder, [3 4]);

            obj2 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            testCase.verifyEqual(nnz(obj2.PurgedEvents), 2, 'loading never modifies the list');
            obj2.saveEvents(testCase.TempFolder);
            testCase.verifyEqual(numel(obj2.eventID), 4);
            testCase.verifyFalse(any(obj2.PurgedEvents));
            S = load(fullfile(testCase.TempFolder, 'events.mat'));
            testCase.verifyEqual(numel(S.eventID), 4);
        end

        function testLegacyPurgeFlagsKeptWithEventSplitFiles(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], 'PulseWidth', 80, 'StimName', 'Main');
            obj1 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            obj1.saveEvents(testCase.TempFolder);
            iFlagLegacyPurge(testCase.TempFolder, [3 4]);
            saveData(fullfile(testCase.TempFolder, 'byEv.dat'), zeros(2, 2, 3, 2, 'single'), ...
                'DimNames', {'Y', 'X', 'T', 'E'}, 'FrameRateHz', 10);

            obj2 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            testCase.verifyWarning(@() obj2.saveEvents(testCase.TempFolder), ...
                'Umitoolbox:EventsManager:legacyPurgeKept');
            testCase.verifyEqual(numel(obj2.eventID), 6);
            testCase.verifyEqual(nnz(obj2.PurgedEvents), 2);
        end

        function testLoadBackwardCompatWithoutPurgedEvents(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj1 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            obj1.saveEvents(testCase.TempFolder);

            % Simulate an events.mat saved before PurgedEvents existed.
            evFile = fullfile(testCase.TempFolder, 'events.mat');
            evData = load(evFile);
            evData = rmfield(evData, 'PurgedEvents');
            save(evFile, '-struct', 'evData');

            obj2 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifyEqual(obj2.PurgedEvents, false(size(obj2.eventID)));
            [frMat1] = obj1.getFrameMatrix(40, 'FrameRateHz', obj1.AcqInfo.FrameRateHz);
            [frMat2] = obj2.getFrameMatrix(40, 'FrameRateHz', obj2.AcqInfo.FrameRateHz);
            testCase.verifyEqual(frMat2, frMat1);
        end

        function testSaveLoadRoundTripExternalSignal(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj1 = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            obj1.getTriggersFromSignal(signal, 10000, false);
            obj1.removeRepetition('extSignal', 2);
            obj1.saveEvents(testCase.TempFolder);

            obj2 = EventsManager(testCase.TempFolder, '', 'csv');

            testCase.verifyTrue(obj2.b_hasExternalSignal);
            testCase.verifyEqual(obj2.AnalogIN, single(signal(:)));
            testCase.verifyEqual(obj2.timestamps, obj1.timestamps);
            testCase.verifyEqual(obj2.state, obj1.state);
            testCase.verifyEqual(obj2.eventID, obj1.eventID);
            testCase.verifyEqual(obj2.eventNameList, {'extSignal'});
            testCase.verifyEqual(obj2.selectedEvents, obj1.selectedEvents);
            testCase.verifyEqual(obj2.sr, single(10000));
        end

        function testExportEventInfoMatchesOnEvents(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            obj.getTriggersFromSignal(signal, 10000, false);
            obj.removeRepetition('extSignal', 2);

            evInfo = obj.exportEventInfo('FrameRateHz', obj.AcqInfo.FrameRateHz);

            testCase.verifyEqual(evInfo.eventNameList, {'extSignal'});
            testCase.verifyEqual(evInfo.FrameRateHz, obj.AcqInfo.FrameRateHz);
            testCase.verifyEqual(evInfo.eventID, obj.eventID(obj.state));
            testCase.verifyEqual(evInfo.selectedEvents, obj.selectedEvents(obj.state));
            testCase.verifyFalse(isfield(evInfo, 'PurgedEvents'));
            testCase.verifyEqual(evInfo.selected, obj.selectedEvents(obj.state));
            testCase.verifyEqual(evInfo.repetitionIndex, obj.repetitionID(obj.state));
            testCase.verifyEqual(evInfo.baselinePeriod, obj.baselinePeriod, 'AbsTol', 1e-6);
        end

        function testPlotUsesDarkThemeAndBaselinePatch(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000], 10000, 'Amplitude', 5);
            obj.getTriggersFromSignal(signal, 10000, false);
            obj.setBaselinePeriod(0.1);

            fig = figure('Visible', 'off', 'Color', [0.1 0.1 0.1]);
            cleanup = onCleanup(@() delete(fig));
            ax = axes(fig);
            obj.plot([], ax);

            trace = findobj(ax, 'Type', 'line', 'Tag', '');
            baseline = findobj(ax, 'Tag', 'BaselinePatch');
            testCase.verifyNotEmpty(trace);
            testCase.verifyEqual(trace(1).Color, [1 1 1], 'AbsTol', 1e-12);
            testCase.verifyNotEmpty(baseline);

            legendObject = legend(ax);
            testCase.verifyTrue(any(contains(string(legendObject.String), 'Baseline (')));
        end
    end
end

function iFlagLegacyPurge(folder, transitions)
%IFLAGLEGACYPURGE Simulate an events.mat written by an older purge.
S = load(fullfile(folder, 'events.mat'));
S.PurgedEvents = false(size(S.eventID));
S.PurgedEvents(transitions) = true;
S.selectedEvents(transitions) = false;
save(fullfile(folder, 'events.mat'), '-struct', 'S');
end
