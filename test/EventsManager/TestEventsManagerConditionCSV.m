classdef TestEventsManagerConditionCSV < matlab.unittest.TestCase
    %TESTEVENTSMANAGERCONDITIONCSV Unit tests for CSV condition-file parsing.

    properties
        TempFolder char
        Obj
    end

    methods (TestMethodSetup)
        function createTestObject(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
            addpath(testCase.TempFolder);
            addpath(fullfile(fileparts(mfilename('fullpath')), 'EMTestHelpers'));
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            testCase.Obj = EventsManager(testCase.TempFolder, '', 'csv');
            testCase.Obj.EventFileParseMethod = 'csv';
        end
    end

    methods (Test)
        function testSingleColumnCSVRelabeling(testCase)
            sr = 10000;
            signal = makePulseSignal(22000, [3000 7000 11000 15000], 60, 'Amplitude', 5);
            testCase.Obj.getTriggersFromSignal(signal, sr, false);

            csvFile = writeCSVConditionFile(testCase.TempFolder, 'events.csv', {'A','B','A','C'});
            status = testCase.Obj.readConditionFile(csvFile);

            testCase.verifyTrue(status);
            testCase.verifyEqual(testCase.Obj.eventNameList(:)', {'A','B','C'});
            testCase.verifyEqual(testCase.Obj.eventID, uint16([1;1;2;2;1;1;3;3]));
        end

        function testMultiColumnCSVRelabeling(testCase)
            sr = 10000;
            signal = makePulseSignal(22000, [3000 7000 11000 15000], 60, 'Amplitude', 5);
            testCase.Obj.getTriggersFromSignal(signal, sr, false);

            spec = struct();
            spec.Type = {'A','B','A','C'};
            spec.Side = {'L','R','L','R'};
            csvFile = writeCSVConditionFile(testCase.TempFolder, 'events_multi.csv', spec);
            status = testCase.Obj.readConditionFile(csvFile, 'CSVcols', {'Type','Side'});

            testCase.verifyTrue(status);
            expNames = {'Type-A-Side-L','Type-B-Side-R','Type-C-Side-R'};
            testCase.verifyEqual(testCase.Obj.eventNameList(:)', expNames);
            testCase.verifyEqual(testCase.Obj.eventID, uint16([1;1;2;2;1;1;3;3]));
        end

        function testCSVTooShortFails(testCase)
            sr = 10000;
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            testCase.Obj.getTriggersFromSignal(signal, sr, false);

            csvFile = writeCSVConditionFile(testCase.TempFolder, 'events_short.csv', {'A','B'});
            didError = false;
            try
                testCase.Obj.readConditionFile(csvFile);
            catch
                didError = true;
            end
            testCase.verifyTrue(didError);
        end

        function testCSVTrimAfterLeadingPartialTrigger(testCase)
            sr = 10000;
            signal = makePulseSignal(24000, [5000 9000 13000], 60, 'Amplitude', 5);
            signal(1:1200) = 5;
            testCase.Obj.getTriggersFromSignal(signal, sr, false);

            csvFile = writeCSVConditionFile(testCase.TempFolder, 'events_trim.csv', ...
                {'DropMe','Keep1','Keep2','Keep3'});
            status = testCase.Obj.readConditionFile(csvFile);

            testCase.verifyTrue(status);
            testCase.verifyEqual(testCase.Obj.eventNameList(:)', {'Keep1','Keep2','Keep3'});
            testCase.verifyEqual(testCase.Obj.eventID, uint16([1;1;2;2;3;3]));
        end
    end
end
