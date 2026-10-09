classdef TestEventsManagerExternalSignal < matlab.unittest.TestCase
    %TESTEVENTSMANAGEREXTERNALSIGNAL Unit tests for getTriggersFromSignal.

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
        function testEdgeSetCleanPulses(testCase)
            sr = 10000;
            signal = makePulseSignal(20000, [2000 6000 10000], 80, 'Amplitude', 5);
            out = testCase.Obj.getTriggersFromSignal(signal, sr, false);

            testCase.verifyEqual(nnz(out.state), 3);
            testCase.verifyEqual(numel(out.timestamps), 6);
            testCase.verifyEqual(out.eventNameList(:)', {'extSignal'});
            testCase.verifyEqual(unique(out.eventID), uint16(1));
            testCase.verifyTrue(all(out.selectedEvents));
            testCase.verifyEqual(out.repetitionID, uint16([1;1;2;2;3;3]));
            testCase.verifyTrue(testCase.Obj.b_hasExternalSignal);
            testCase.verifyEqual(testCase.Obj.eventNameList(:)', {'extSignal'});
        end

        function testEdgeToggleInterpretation(testCase)
            sr = 10000;
            testCase.Obj.trigType = 'EdgeToggle';
            signal = makePulseSignal(16000, [2000 6000 10000 14000], 50, 'Amplitude', 5);
            out = testCase.Obj.getTriggersFromSignal(signal, sr, false);

            testCase.verifyEqual(numel(out.timestamps), 4);
            testCase.verifyEqual(out.state, logical([1;0;1;0]));
            testCase.verifyEqual(out.repetitionID, uint16([1;1;2;2]));
        end

        function testLeadingPartialEventTrim(testCase)
            sr = 10000;
            signal = makePulseSignal(20000, [5000 9000], 100, 'Amplitude', 5);
            signal(1:1500) = 5; % Start while already HIGH.
            out = testCase.Obj.getTriggersFromSignal(signal, sr, false);

            testCase.verifyEqual(out.triggerTrimInfo.nLeadDropped, 1);
            testCase.verifyEqual(nnz(out.state), 2);
            testCase.verifyEqual(numel(out.timestamps), 4);
        end

        function testTrailingPartialEventTrim(testCase)
            sr = 10000;
            signal = makePulseSignal(20000, [3000 8000], 100, 'Amplitude', 5);
            signal(15000:end) = 5; % End while still HIGH.
            out = testCase.Obj.getTriggersFromSignal(signal, sr, false);

            testCase.verifyEqual(out.triggerTrimInfo.nTrailDropped, 1);
            testCase.verifyEqual(nnz(out.state), 2);
            testCase.verifyEqual(numel(out.timestamps), 4);
        end

        function testNoTriggerSignalLeavesBaselineEmpty(testCase)
            sr = 10000;
            signal = zeros(20000, 1, 'single');
            out = testCase.Obj.getTriggersFromSignal(signal, sr, false);

            testCase.verifyEmpty(out.timestamps);
            testCase.verifyEmpty(out.eventID);
            testCase.verifyEmpty(out.baselinePeriod);
            testCase.verifyEmpty(testCase.Obj.baselinePeriod);
        end

        function testAutomaticBaselineUsesShortestOnsetInterval(testCase)
            sr = 1000;
            signal = makePulseSignal(14000, [2000 6000 12000], 50, ...
                'Amplitude', 5);

            out = testCase.Obj.getTriggersFromSignal(signal, sr, false);

            % Onset intervals are 4 s and 6 s. EventsManager deliberately
            % uses 20% of the shortest interval (0.8 s), not the previous
            % split-data fallback's rounded-mean result (1.0 s).
            testCase.verifyEqual(double(out.baselinePeriod), 0.8, ...
                'AbsTol', 1e-6);
        end
    end
end
