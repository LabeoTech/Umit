classdef TestEventTimeVector < matlab.unittest.TestCase
    %TESTEVENTTIMEVECTOR Event-locked time axis (.dat header Phase 8b).
    %
    %   Time 0 is column round(baselinePeriod * frameRateHz), the onset
    %   column of EventsManager.getFrameMatrix (checked against EventsManager
    %   in TestEventsManagerSelectionAndSplit).

    methods (Test)
        function zeroAtTheRoundedBaselineColumn(testCase)
            % The StimAna1 case: 3.1133 s at 5 Hz -> onset column 16 (not 17).
            t = eventTimeVector(78, 5, 3.1133);
            testCase.verifySize(t, [1 78]);
            testCase.verifyEqual(t(16), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(t(1), -15 / 5, 'AbsTol', 1e-12);
            testCase.verifyEqual(diff(t), repmat(0.2, 1, 77), 'AbsTol', 1e-12);
        end

        function noBaselineStartsAtZero(testCase)
            testCase.verifyEqual(eventTimeVector(4, 10, []), (0:3) / 10, 'AbsTol', 1e-12);
            testCase.verifyEqual(eventTimeVector(4, 10, NaN), (0:3) / 10, 'AbsTol', 1e-12);
        end

        function invalidInputsError(testCase)
            testCase.verifyError(@() eventTimeVector(0, 10, 1), 'MATLAB:eventTimeVector:expectedPositive');
            testCase.verifyError(@() eventTimeVector(5, 0, 1), 'MATLAB:eventTimeVector:expectedPositive');
            testCase.verifyError(@() eventTimeVector(5, 10, -1), 'MATLAB:eventTimeVector:expectedNonnegative');
        end
    end
end
