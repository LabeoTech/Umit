classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestSplitTrialFramePlan < matlab.unittest.TestCase
    %TESTSPLITTRIALFRAMEPLAN Unit tests of splitTrialFramePlan and padLeadingBaseline.
    %
    %   Frame matrices are written as getFrameMatrix returns them: one row per
    %   trial, NaN for frames outside the recording. The baseline is
    %   evOnsetFrame columns, the onset is the next column.

    methods (TestMethodSetup)
        function addPaths(testCase) %#ok<MANU>
            testFolder = fileparts(mfilename('fullpath'));
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(fileparts(testFolder)));
            end
            addpath(genpath(char(projectRoot)));
        end
    end

    methods (Test)
        function testCompleteMatrixIsUnchanged(testCase)
            fr = [1:6; 7:12];

            [out, nLead, ok, msg] = splitTrialFramePlan(fr, 3, 10);

            testCase.verifyEqual(out, fr);
            testCase.verifyEqual(nLead, 0);
            testCase.verifyTrue(ok);
            testCase.verifyEmpty(msg);
        end

        function testTrailingNaNCropsAllTrials(testCase)
            fr = [1:6; 7:12; 13:17 NaN];

            [out, nLead, ok] = splitTrialFramePlan(fr, 3, 10);

            testCase.verifyEqual(out, fr(:, 1:5));
            testCase.verifyEqual(nLead, 0);
            testCase.verifyTrue(ok);
        end

        function testLeadingNaNOfTheFirstTrialIsKeptAndCounted(testCase)
            % 3 baseline columns, the first trial lacks its first 2.
            fr = [NaN NaN 1 2 3 4; 7:12];

            [out, nLead, ok, msg] = splitTrialFramePlan(fr, 3, 10);

            testCase.verifyEqual(out, fr);
            testCase.verifyEqual(nLead, 2);
            testCase.verifyTrue(ok);
            testCase.verifySubstring(msg, '2 frames (0.2 s) before the recording');
            testCase.verifySubstring(msg, 'mean of its 1 valid baseline');
            testCase.verifySubstring(msg, 'reduce the baselinePeriod parameter in the EventsManager');
            testCase.verifySubstring(msg, 'run this function again');
        end

        function testLeadingNaNAndTrailingNaNAreTreatedSeparately(testCase)
            fr = [NaN 1 2 3 4 5; 7:12; 13:17 NaN];

            [out, nLead, ok] = splitTrialFramePlan(fr, 3, 10);

            testCase.verifyEqual(out, fr(:, 1:5));
            testCase.verifyEqual(nLead, 1);
            testCase.verifyTrue(ok);
        end

        function testFirstTrialWithoutAnyBaselineFrameIsNotOk(testCase)
            % All 3 baseline columns are outside the recording.
            fr = [NaN NaN NaN 1 2 3; 7:12];

            [~, ~, ok] = splitTrialFramePlan(fr, 3, 10);

            testCase.verifyFalse(ok);
        end

        function testFirstTrialEntirelyOutsideIsNotOk(testCase)
            [~, ~, ok] = splitTrialFramePlan([NaN NaN NaN; 1 2 3], 1, 10);

            testCase.verifyFalse(ok);
        end

        function testLaterTrialStartingOutsideIsNotOk(testCase)
            [~, ~, ok] = splitTrialFramePlan([1:6; NaN 2:6], 3, 10);

            testCase.verifyFalse(ok);
        end

        function testEmptyMatrixIsNotOk(testCase)
            [~, ~, ok] = splitTrialFramePlan([], 3, 10);

            testCase.verifyFalse(ok);
        end

        function testPadFillsWithThePerPixelMeanOfTheValidBaseline(testCase)
            % 4 baseline columns, the first trial lacks 2: they take the
            % mean of columns 3-4, pixel by pixel (non-square, so axes mix up).
            trial = nan(2, 3, 6, 'single');
            trial(:, :, 3:6) = reshape(single(1:2*3*4), 2, 3, 4);

            out = padLeadingBaseline(trial, 2, 4);

            expectedFill = mean(trial(:, :, 3:4), 3);
            testCase.verifyEqual(out(:, :, 1), expectedFill);
            testCase.verifyEqual(out(:, :, 2), expectedFill);
            testCase.verifyEqual(out(:, :, 3:6), trial(:, :, 3:6));
            testCase.verifyClass(out, 'single');
        end

        function testPadWithNoMissingFramesIsANoOp(testCase)
            trial = rand(2, 3, 6, 'single');

            testCase.verifyEqual(padLeadingBaseline(trial, 0, 4), trial);
        end

        function testPadWithoutAValidBaselineFrameErrors(testCase)
            trial = nan(2, 3, 6, 'single');

            testCase.verifyError(@() padLeadingBaseline(trial, 4, 4), ...
                'Umitoolbox:padLeadingBaseline:noBaselineFrame');
        end
    end
end
