classdef TestEventsManagerTrigPolarity < matlab.unittest.TestCase
    %TESTEVENTSMANAGERTRIGPOLARITY Tests for trigPolarity behavior.

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
        function testTrigPolaritySetterAcceptsMixedCase(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');

            obj.trigPolarity = 'NEGATIVE';
            testCase.verifyEqual(obj.trigPolarity, 'negative');

            obj.trigPolarity = 'Positive';
            testCase.verifyEqual(obj.trigPolarity, 'positive');
        end

        function testTrigPolaritySetterRejectsInvalidValue(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');

            didError = false;
            try
                obj.trigPolarity = 'sideways';
            catch
                didError = true;
            end
            testCase.verifyTrue(didError);
        end

        function testNegativePolarityExternalSignalSucceedsWherePositiveFails(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');

            signal = makePulseSignal(20000, [3000 8000 13000], 80, ...
                'Amplitude', 5, 'Baseline', 0, 'Polarity', -1);

            obj.trigPolarity = 'positive';
            outPositive = obj.getTriggersFromSignal(signal, 10000, false);

            obj.trigPolarity = 'negative';
            outNegative = obj.getTriggersFromSignal(signal, 10000, false);

            testCase.verifyTrue(isempty(outPositive.timestamps) || nnz(outPositive.state) ~= 3, ...
                ['Positive polarity should not interpret a purely negative-going trigger train ' ...
                 'as three clean onset events.']);

            testCase.verifyEqual(nnz(outNegative.state), 3);
            testCase.verifyEqual(numel(outNegative.timestamps), 6);
            testCase.verifyEqual(outNegative.eventNameList(:)', {'extSignal'});
            testCase.verifyEqual(outNegative.repetitionID, uint16([1;1;2;2;3;3]));
        end

        function testNegativePolarityFileBackedAnalogDetection(testCase)
            folderPath = testCase.TempFolder;
            if ~isfolder(folderPath)
                mkdir(folderPath);
            end

            nSamples = 25000;
            nChan = 12;
            data = zeros(nSamples, nChan, 'single');
            data(1001:24000, 1) = 5; % Camera trigger
            data(:, 2) = makePulseSignal(nSamples, [5000 10000 15000], 80, ...
                'Amplitude', 5, 'Baseline', 0, 'Polarity', -1);

            writeInfoTxtStimAna1(folderPath, ...
                'FrameRateHz', 60, ...
                'AISampleRate', 10000, ...
                'StimName', 'Main', ...
                'StimNRepeat', 3);
            writeSyntheticAIFile(folderPath, data, 'FileName', 'ai_00000.bin');

            obj = EventsManager(folderPath, folderPath, 'csv');
            obj.trigPolarity = 'negative';
            obj.getTriggers('',false);

            testCase.verifyFalse(obj.b_isDigital);
            testCase.verifyEqual(obj.trigChanName, {'StimAna1'});
            testCase.verifyEqual(obj.eventNameList(:)', {'Main'});
            testCase.verifyEqual(nnz(obj.state), 3);
            testCase.verifyEqual(numel(obj.timestamps), 6);
            testCase.verifyTrue(all(obj.selectedEvents));
        end

        function testSaveLoadRoundTripPreservesTrigPolarity(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj1 = EventsManager(testCase.TempFolder, '', 'csv');
            obj1.trigPolarity = 'negative';
            signal = makePulseSignal(22000, [3000 7000 11000], 60, ...
                'Amplitude', 5, 'Baseline', 0, 'Polarity', -1);
            obj1.getTriggersFromSignal(signal, 10000, false);
            obj1.saveEvents(testCase.TempFolder);

            obj2 = EventsManager(testCase.TempFolder, '', 'csv');

            testCase.verifyEqual(obj2.trigPolarity, 'negative');
            testCase.verifyTrue(obj2.b_hasExternalSignal);
            testCase.verifyEqual(obj2.timestamps, obj1.timestamps);
            testCase.verifyEqual(obj2.state, obj1.state);
            testCase.verifyEqual(obj2.eventID, obj1.eventID);
        end
    end
end
