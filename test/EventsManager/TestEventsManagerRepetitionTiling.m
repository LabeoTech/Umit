classdef TestEventsManagerRepetitionTiling < matlab.unittest.TestCase
    %TESTEVENTSMANAGERREPETITIONTILING Regression for DFR-20260921-001.
    %
    % Reproduces two digital conditions whose declared durations tile
    % continuous time with no gap (one condition's ON registers before the
    % other condition's own OFF), which previously broke the positional
    % ON/OFF adjacency assumed by get.repetitionID's repelem-based doubling.

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
        function testRepetitionIDWithNoGapTiling(testCase)
            % ON/OFF chronology this produces (id1 dur=1s, id2 dur=3s):
            %   t=2.0 ON id1 rep1   t=2.5 ON id2 rep1   t=3.0 OFF id1 rep1
            %   t=5.0 ON id1 rep2   t=5.5 OFF id2 rep1  t=6.0 OFF id1 rep2
            %   t=6.5 ON id2 rep2   t=9.5 OFF id2 rep2
            % id2's rep1 ON (t=2.5) precedes id1's rep1 OFF (t=3.0), and id1's
            % rep2 ON (t=5.0) precedes id2's rep1 OFF (t=5.5): both break the
            % "next array slot is this event's own OFF" assumption.
            buildDigitalAcquisitionFolder(testCase.TempFolder, ...
                'AISampleRate', 1000, ...
                'NSamples', 12000, ...
                'EventOrder', [1 2 1 2], ...
                'EventNames', {'PuffON','baseline'}, ...
                'Durations', [1 3], ...
                'PulseStarts', [2000 2500 5000 6500]);

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'none');

            testCase.assertTrue(obj.b_isDigital);
            testCase.assertEqual(obj.eventNameList(:)', {'PuffON','baseline'});

            repID = obj.repetitionID;
            testCase.verifyEqual(repID, uint16([1;1;1;2;1;2;2;2]));

            for condID = 1:numel(obj.eventNameList)
                onRep = repID(obj.state & obj.eventID == condID);
                offRep = repID(~obj.state & obj.eventID == condID);
                testCase.verifyEqual(onRep, uint16((1:numel(onRep))'), ...
                    sprintf('ON repetition labels for condition %d must be sequential.', condID));
                testCase.verifyEqual(offRep, onRep, ...
                    sprintf('OFF repetition labels for condition %d must match their own ON.', condID));
            end
        end
    end
end
