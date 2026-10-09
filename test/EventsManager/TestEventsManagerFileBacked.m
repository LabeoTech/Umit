classdef TestEventsManagerFileBacked < matlab.unittest.TestCase
    %TESTEVENTSMANAGERFILEBACKED Integration tests for file-backed trigger detection.

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
        function testAnalogStimAna1Detection(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');

            testCase.verifyFalse(obj.b_isDigital);
            testCase.verifyEqual(obj.trigChanName, {'StimAna1'});
            testCase.verifyEqual(obj.eventNameList(:)', {'Main'});
            testCase.verifyEqual(nnz(obj.state), 3);
            testCase.verifyEqual(numel(obj.timestamps), 6);
            testCase.verifyTrue(all(obj.selectedEvents));
        end

        function testDigitalStimDigDetection(testCase)
            buildDigitalAcquisitionFolder(testCase.TempFolder, ...
                'EventOrder', [3 1 2], ...
                'EventNames', {'block1','block2','block3'}, ...
                'Durations', [3 5 1], ...
                'PulseStarts', [5000 9000 13000]);

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'none');

            testCase.verifyTrue(obj.b_isDigital);
            testCase.verifyEqual(obj.trigChanName, {'StimDig'});
            testCase.verifyEqual(obj.eventNameList(:)', {'block1','block2','block3'});
            testCase.verifyEqual(obj.eventID(obj.state), uint16([3;1;2]));
            testCase.verifyEqual(nnz(obj.state), 3);

            onIDs = obj.eventID(obj.state);
            offIDs = obj.eventID(~obj.state);
            onTimes = obj.timestamps(obj.state);
            offTimes = obj.timestamps(~obj.state);
            actualDur = zeros(numel(onIDs), 1, 'single');
            for ii = 1:numel(onIDs)
                actualDur(ii) = offTimes(find(offIDs == onIDs(ii), 1, 'first')) - onTimes(ii);
            end
            expDur = single([3;5;1]);
            expDur = expDur(double(onIDs));
            testCase.verifyEqual(actualDur, expDur, 'AbsTol', 1e-4);
        end
    end
end
