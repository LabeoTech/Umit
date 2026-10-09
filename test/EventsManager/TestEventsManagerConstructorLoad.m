classdef TestEventsManagerConstructorLoad < matlab.unittest.TestCase
    %TESTEVENTSMANAGERCONSTRUCTORLOAD Constructor and load-path behavior tests.

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
        function testConstructorLoadsExternalEventsMatWithoutMetadata(testCase)
            % Build and save external-signal events in a folder that
            % initially has minimal metadata.
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj1 = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            obj1.getTriggersFromSignal(signal, 10000, false);
            obj1.saveEvents(testCase.TempFolder);

            % Remove metadata so constructor must rely on events.mat only.
            delete(fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            obj2 = EventsManager(testCase.TempFolder, '', 'csv');

            testCase.verifyTrue(obj2.b_hasExternalSignal);
            testCase.verifyEqual(obj2.AnalogIN, single(signal(:)));
            testCase.verifyEqual(obj2.timestamps, obj1.timestamps);
            testCase.verifyEqual(obj2.eventID, obj1.eventID);
            testCase.verifyEqual(obj2.eventNameList, {'extSignal'});
            testCase.verifyEqual(obj2.selectedEvents, obj1.selectedEvents);
            testCase.verifyEqual(obj2.sr, single(10000));
        end

        function testConstructorAutoInitializesFromRawFolderAnalog(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifyFalse(isempty(obj.AcqInfo));
            testCase.verifyFalse(isempty(obj.AnalogIN));
            testCase.verifyEqual(obj.trigChanName, {'StimAna1'});
            testCase.verifyEqual(obj.eventNameList(:)', {'Main'});
            testCase.verifyEqual(nnz(obj.state), 3);
        end

        function testLoadEventsRepairsInvalidSelectedEvents(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj1 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            obj1.saveEvents(testCase.TempFolder);

            ev = load(fullfile(testCase.TempFolder, 'events.mat'));
            ev.selectedEvents = true(1, 3); % intentionally wrong size
            save(fullfile(testCase.TempFolder, 'events.mat'), '-struct', 'ev');

            obj2 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifySize(obj2.selectedEvents, size(obj2.eventID));
            testCase.verifyTrue(all(obj2.selectedEvents));
        end

        function testConstructorUsesExistingEventsMatBeforeRawProcessing(testCase)
            % Save external events.mat in a folder that also contains an
            % invalid raw acquisition layout. Constructor should still load
            % from events.mat first.
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj1 = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000], 60, 'Amplitude', 5);
            obj1.getTriggersFromSignal(signal, 10000, false);
            obj1.saveEvents(testCase.TempFolder);

            % Add an invalid raw file sequence that would fail if processed.
            writeInfoTxtStimAna1(testCase.TempFolder);
            data = zeros(12000, 12, 'single');
            data(1001:11000,1) = 5;
            data(5000:5079,2) = 5;
            writeSyntheticAIFile(testCase.TempFolder, data, 'FileName', 'ai_00001.bin');

            obj2 = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifyTrue(obj2.b_hasExternalSignal);
            testCase.verifyEqual(obj2.eventNameList, {'extSignal'});
            testCase.verifyEqual(obj2.timestamps, obj1.timestamps);
        end
    end
end
