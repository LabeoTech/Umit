classdef TestEventsManagerMalformedInputs < matlab.unittest.TestCase
    %TESTEVENTSMANAGERMALFORMEDINPUTS Malformed-input handling tests.

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
        function testReadConditionFileInvalidCSVColumnFails(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            obj.getTriggersFromSignal(signal, 10000, false);

            csvFile = writeCSVConditionFile(testCase.TempFolder, 'events.csv', {'A','B','C'});

            didError = false;
            try
                obj.readConditionFile(csvFile, 'CSVcols', {'NoSuchColumn'});
            catch
                didError = true;
            end
            testCase.verifyTrue(didError);
        end

        function testCorruptedAnalogPayloadLengthFails(testCase)
            writeInfoTxtStimAna1(testCase.TempFolder);
            fid = fopen(fullfile(testCase.TempFolder, 'ai_00000.bin'), 'w');
            cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid, zeros(5,1,'int32'), 'int32');
            fwrite(fid, [1 2 3], 'double'); % invalid payload length

            didError = false;
            try
                obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv'); %#ok<NASGU>
            catch
                didError = true;
            end
            testCase.verifyTrue(didError);
        end

        function testMissingAnalogFileGapFails(testCase)
            writeInfoTxtStimAna1(testCase.TempFolder);
            data = zeros(12000, 12, 'single');
            data(1001:11000,1) = 5;
            data(5000:5079,2) = 5;
            writeSyntheticAIFile(testCase.TempFolder, data, 'FileName', 'ai_00000.bin');
            writeSyntheticAIFile(testCase.TempFolder, data, 'FileName', 'ai_00002.bin');

            didError = false;
            try
                obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv'); %#ok<NASGU>
            catch
                didError = true;
            end
            testCase.verifyTrue(didError);
        end

        function testMalformedInfoTxtLeavesObjectUninitialized(testCase)
            fid = fopen(fullfile(testCase.TempFolder, 'info.txt'), 'w');
            cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fprintf(fid, 'this is not a valid acquisition metadata file\n');

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifyTrue(isempty(fieldnames(obj.AcqInfo)) || isempty(obj.AcqInfo));
            testCase.verifyTrue(isempty(obj.AnalogIN));
            testCase.verifyTrue(isempty(obj.eventID));
        end
    end
end
