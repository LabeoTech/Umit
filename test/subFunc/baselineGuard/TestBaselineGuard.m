classdef TestBaselineGuard < matlab.unittest.TestCase
    %TESTBASELINEGUARD countLowBaseline and warnIfBaselineNearZero.
    %
    %   The two helpers let a step that divides by a per-trace baseline notice
    %   zero-centered input (already normalized, z-scored or baseline-subtracted
    %   data), whose baseline is within the noise.

    properties (Constant)
        WarningId = 'Umitoolbox:testStep:ZeroCenteredInput'
    end

    methods (Test)
        function intensityLikeTracesAreNotLow(testCase)
            baseline = [1000 800 1200];
            fluctuation = [20 25 30];

            [nChecked, nLow] = countLowBaseline(baseline, fluctuation);

            testCase.verifyEqual(nChecked, 3);
            testCase.verifyEqual(nLow, 0);
        end

        function zeroCenteredTracesAreLow(testCase)
            baseline = [0 0.01 -0.2 0.5];
            fluctuation = [1 0.5 0.4 0.3];

            [nChecked, nLow] = countLowBaseline(baseline, fluctuation);

            testCase.verifyEqual(nChecked, 4);
            testCase.verifyEqual(nLow, 4);
        end

        function thresholdIsThreeTimesTheFluctuation(testCase)
            % Strictly below 3 x fluctuation counts as low.
            [~, nBelow] = countLowBaseline(2.99, 1);
            [~, nAt] = countLowBaseline(3, 1);
            [~, nAbove] = countLowBaseline(-3.01, 1);

            testCase.verifyEqual(nBelow, 1);
            testCase.verifyEqual(nAt, 0);
            testCase.verifyEqual(nAbove, 0);
        end

        function maskedAndNonFiniteTracesAreNotExamined(testCase)
            baseline = [0 NaN 5 Inf 1000];
            fluctuation = [0 1 NaN 2 20];

            [nChecked, nLow] = countLowBaseline(baseline, fluctuation);

            testCase.verifyEqual(nChecked, 1, ...
                'Only the last trace is finite and not masked.');
            testCase.verifyEqual(nLow, 0);
        end

        function constantNonZeroTraceIsNotLow(testCase)
            % A flat trace at a large level has fluctuation 0 and a usable baseline.
            [nChecked, nLow] = countLowBaseline(1000, 0);

            testCase.verifyEqual(nChecked, 1);
            testCase.verifyEqual(nLow, 0);
        end

        function arraysOfAnyShapeAreCounted(testCase)
            baseline = ones(4, 5, 1, 2) * 1000;
            fluctuation = ones(4, 5, 1, 2) * 10;
            baseline(1:2, :, 1, 1) = 0.1;

            [nChecked, nLow] = countLowBaseline(baseline, fluctuation);

            testCase.verifyEqual(nChecked, 40);
            testCase.verifyEqual(nLow, 10);
        end

        function sizeMismatchIsRejected(testCase)
            testCase.verifyError(@() countLowBaseline(ones(2, 3), ones(3, 2)), ...
                'Umitoolbox:countLowBaseline:sizeMismatch');
        end

        function noWarningBelowTheShareThreshold(testCase)
            testCase.verifyWarningFree(@() warnIfBaselineNearZero( ...
                testCase.WarningId, 'testStep', 100, 19, 'Do something.'));
        end

        function warnsFromTheShareThresholdUp(testCase)
            testCase.verifyWarning(@() warnIfBaselineNearZero( ...
                testCase.WarningId, 'testStep', 100, 20, 'Do something.'), ...
                testCase.WarningId);
            testCase.verifyWarning(@() warnIfBaselineNearZero( ...
                testCase.WarningId, 'testStep', 100, 100, 'Do something.'), ...
                testCase.WarningId);
        end

        function noWarningWhenNothingWasExamined(testCase)
            testCase.verifyWarningFree(@() warnIfBaselineNearZero( ...
                testCase.WarningId, 'testStep', 0, 0, 'Do something.'));
        end

        function messageNamesTheStepTheShareAndTheAdvice(testCase)
            lastwarn('');
            warning('off', testCase.WarningId);
            cleanup = onCleanup(@() warning('on', testCase.WarningId));
            warnIfBaselineNearZero(testCase.WarningId, 'testStep', 80, 60, 'Use raw data.');
            [message, id] = lastwarn;
            clear cleanup

            testCase.verifyEqual(id, testCase.WarningId);
            testCase.verifySubstring(message, 'testStep');
            testCase.verifySubstring(message, '75%');
            testCase.verifySubstring(message, '80 traces');
            testCase.verifySubstring(message, 'Use raw data.');
        end
    end
end
