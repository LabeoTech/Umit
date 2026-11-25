classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestEventViewFirstTrialFill < matlab.unittest.TestCase
    %TESTEVENTVIEWFIRSTTRIALFILL Event view of a first trial that starts before the recording.
    %
    %   Covers planEventViewFrames, padEventTrialTrace and eventAverageFrame,
    %   the backend of DataViewer's .dat event view. Frame matrices are written
    %   as EventsManager.getFrameMatrix returns them (NaN = outside the
    %   recording); the baseline is evOnsetFrame columns.

    methods (TestMethodSetup)
        function addPaths(testCase) %#ok<MANU>
            testFolder = fileparts(mfilename('fullpath'));
            projectRoot = fileparts(fileparts(testFolder));
            addpath(fullfile(projectRoot, 'GUI', 'DataViewer'));
            addpath(fullfile(projectRoot, 'subFunc'));
        end
    end

    methods (Test)
        function testAverageKeepsLeadingNaNOfFirstTrial(testCase)
            % The old crop at the first NaN column emptied the matrix.
            fr = [NaN NaN 1 2 3 4; 5 6 7 8 9 10];

            [out, nLead, evOnsetFrame] = planEventViewFrames(fr, 0.4, 10, true);

            testCase.verifyEqual(out, fr);
            testCase.verifyEqual(nLead, 2);
            testCase.verifyEqual(evOnsetFrame, 4);
        end

        function testAverageCropsAtFirstMissingFrameOfAnyTrial(testCase)
            fr = [NaN 1 2 3 4 5; 6 7 8 9 10 NaN];

            [out, nLead] = planEventViewFrames(fr, 0.3, 10, true);

            testCase.verifyEqual(out, fr(:, 1:5));
            testCase.verifyEqual(nLead, 1);
        end

        function testAverageFallsBackToOldCropWhenPlanFails(testCase)
            % A later trial starting outside the recording is not fillable.
            fr = [1:5; NaN 2:5];

            [out, nLead] = planEventViewFrames(fr, 0.3, 10, true);

            testCase.verifyEmpty(out);
            testCase.verifyEqual(nLead, 0);
        end

        function testSingleRepetitionIsNotCropped(testCase)
            fr = [NaN NaN 1 2 3 NaN];

            [out, nLead] = planEventViewFrames(fr, 0.4, 10, false);

            testCase.verifyEqual(out, fr);
            testCase.verifyEqual(nLead, 2);
        end

        function testSingleRepetitionWithoutBaselineFrameIsNotFilled(testCase)
            % The onset (column 3) is the first recording frame: no baseline to fill from.
            fr = [NaN NaN 1 2 3];

            [~, nLead] = planEventViewFrames(fr, 0.2, 10, false);

            testCase.verifyEqual(nLead, 0);
        end

        function testPadTraceUsesBaselineMean(testCase)
            trace = [NaN NaN 2 4 10 20];

            out = padEventTrialTrace(trace, 2, 4);

            testCase.verifyEqual(out, [3 3 2 4 10 20]);
        end

        function testPadTraceKeepsColumnOrientation(testCase)
            trace = [NaN; 2; 4; 9];

            out = padEventTrialTrace(trace, 1, 3);

            testCase.verifyEqual(out, [3; 2; 4; 9]);
        end

        function testPadTraceWithoutLeadIsUnchanged(testCase)
            trace = [1 2 3];

            testCase.verifyEqual(padEventTrialTrace(trace, 0, 2), trace);
        end

        function testAverageFrameFillsFirstTrialWithBaselineMean(testCase)
            % Non-square frames (Y ~= X); frame k is the constant image k.
            getFrame = @(k) k * ones(3, 5, 'single');
            fr = [NaN NaN 1 2 3 4; 11 12 13 14 15 16];

            % Column 1: trial 2 = 11, trial 1 filled with mean(frames 1,2) = 1.5.
            [frame, padFrame] = eventAverageFrame(getFrame, fr, 1, 2, 4, []);

            testCase.verifyEqual(size(frame), [3 5]);
            testCase.verifyEqual(frame, single((11 + 1.5) / 2) * ones(3, 5, 'single'), 'AbsTol', 1e-5);
            testCase.verifyEqual(padFrame, 1.5 * ones(3, 5, 'single'), 'AbsTol', 1e-5);

            % A cached pad frame is reused: baseline frames 1-2 are not read again.
            frame2 = eventAverageFrame(@iReadAbove10, [NaN NaN; 11 12], 2, 2, 4, padFrame);
            testCase.verifyEqual(frame2, single((12 + 1.5) / 2) * ones(3, 5, 'single'), 'AbsTol', 1e-5);
        end

        function testAverageFrameAfterLeadAveragesAllTrials(testCase)
            getFrame = @(k) k * ones(2, 4, 'single');
            fr = [NaN NaN 1 2; 11 12 13 14];

            frame = eventAverageFrame(getFrame, fr, 3, 2, 3, []);

            testCase.verifyEqual(frame, single((1 + 13) / 2) * ones(2, 4, 'single'), 'AbsTol', 1e-5);
        end

        function testAverageFrameSingleRepetitionFillsLeadingColumns(testCase)
            getFrame = @(k) k * ones(2, 3, 'single');
            fr = [NaN 4 6 8];

            frame = eventAverageFrame(getFrame, fr, 1, 1, 2, []);

            testCase.verifyEqual(frame, 4 * ones(2, 3, 'single'), 'AbsTol', 1e-5);
        end

        function testAverageFrameWithoutAnyFrameIsEmpty(testCase)
            frame = eventAverageFrame(@(k) k, [NaN 1; NaN 2], 1, 0, 1, []);

            testCase.verifyEmpty(frame);
        end
    end
end

function img = iReadAbove10(k)
%IREADABOVE10 Frame reader that fails on a baseline frame (1-2) of the first trial.
assert(k > 10, 'TestEventViewFirstTrialFill:unexpectedRead', ...
    'Frame %d should not be read: the baseline mean is cached.', k);
img = k * ones(3, 5, 'single');
end
