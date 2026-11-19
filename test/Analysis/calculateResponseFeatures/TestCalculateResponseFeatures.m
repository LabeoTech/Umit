classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestCalculateResponseFeatures < matlab.unittest.TestCase
    %TESTCALCULATERESPONSEFEATURES Unit tests for calculateResponseFeatures.

    methods (Test)

        function testPipelineInfo(testCase)
            info = calculateResponseFeatures('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyEqual(info.name, 'calculateResponseFeatures');

            timeWindow = info.parameters(strcmp( ...
                {info.parameters.name}, 'TimeWindow_sec'));
            testCase.verifyEqual(timeWindow.allowed, {'all',[0 Inf]});

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');
        end

        function testRejectsUnsupportedInputs(testCase)
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            umt = iBuildValidFixture();

            % UMT structs and arrays are not input forms.
            testCase.verifyError(@() calculateResponseFeatures(umt), ...
                'Umitoolbox:calculateResponseFeatures:InvalidInputType');
            testCase.verifyError(@() calculateResponseFeatures(rand(2, 12, 2)), ...
                'Umitoolbox:calculateResponseFeatures:InvalidInputType');

            % Missing file, other extensions.
            testCase.verifyError(@() calculateResponseFeatures( ...
                fullfile(folder, 'missing.umt')), ...
                'Umitoolbox:calculateResponseFeatures:InputFileNotFound');
            matFile = fullfile(folder, 'x.mat');
            fclose(fopen(matFile, 'w'));
            testCase.verifyError(@() calculateResponseFeatures(matFile), ...
                'Umitoolbox:calculateResponseFeatures:UnsupportedInputFile');

            % An image UMT is not a ROI UMT.
            imageUMT = genUMTStruct(rand(4, 4, 12, 'single'), 'kind', 'image', ...
                'entryName', 'main', 'dimNames', {'Y','X','T'});
            imageFile = fullfile(folder, 'image.umt');
            saveData(imageFile, imageUMT);
            testCase.verifyError(@() calculateResponseFeatures(imageFile), ...
                'Umitoolbox:calculateResponseFeatures:InvalidUMTKind');
        end

        function testValidInput(testCase)
            umt = iBuildValidFixture();
            out = iRun(testCase, umt);

            validateUMTStruct(out, 'requireEventInfo', true);
            testCase.verifyEqual(char(string(out.kind)), 'roi');

            entryName = fieldnames(out.data);
            entryName = entryName{1};
            testCase.verifyEqual(cellstr(string(out.data.(entryName).dimNames)), {'ROI','Measure','E'});
            testCase.verifyEqual(size(out.data.(entryName).value), [2 6 2]);
        end

        function testRejectPixelInput(testCase)
            umt = iBuildPixelFixture();
            testCase.verifyError(@() iRun(testCase, umt), ...
                'Umitoolbox:calculateResponseFeatures:PixelDimensionNotSupported');
        end

        function testRejectMissingFrameRate(testCase)
            umt = iBuildValidFixture();
            entryName = fieldnames(umt.data);
            entryName = entryName{1};
            umt.data.(entryName).meta = struct();

            testCase.verifyError(@() iRun(testCase, umt), ...
                'Umitoolbox:calculateResponseFeatures:MissingFrameRate');
        end

        function testRejectInvalidTimeWindow(testCase)
            umt = iBuildValidFixture();
            testCase.verifyError(@() iRun(testCase, umt, 'TimeWindow_sec', [3 1]), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end

        function testTimeWindowStartingAtZeroIsAccepted(testCase)
            % 'allowed' advertises [0 Inf], and "0 to N seconds after onset"
            % is the natural request, but the window used to be anchored one
            % frame early and the guard rejected it (P1-8).
            umt = iBuildValidFixture();

            out = iRun(testCase, umt, 'TimeWindow_sec', [0 0.6]);

            entryName = fieldnames(out.data);
            vals = out.data.(entryName{1}).value;

            % Frame 6 is t = 0, so [0 0.6] must cover the whole response and
            % agree with the 'all' window.
            outAll = iRun(testCase, umt, 'TimeWindow_sec', 'all');
            testCase.verifyEqual(vals, outAll.data.(entryName{1}).value, 'AbsTol', 1e-6);
        end

        function testLatenciesAreMeasuredFromFirstResponseFrame(testCase)
            % A response peaking on the first frame after the baseline is at
            % t = 0, not t = 1/Fs (P1-8).
            umt = iBuildValidFixture();

            out = iRun(testCase, umt);
            entryName = fieldnames(out.data);
            vals = out.data.(entryName{1}).value;

            % ROI1/event1 rises on frame 6 (t = 0) and peaks on frame 8.
            onsetLatency = vals(1, 6, 1);
            peakLatency = vals(1, 2, 1);

            testCase.verifyEqual(double(onsetLatency), 0, 'AbsTol', 1e-6);
            testCase.verifyEqual(double(peakLatency), 0.2, 'AbsTol', 1e-6);
        end

        function testAUCIntegratesBaselineCorrectedWindow(testCase)
            % The old implementation subtracted a baseline integral taken
            % over a different number of samples, so the correction scaled
            % with the ratio of the two window lengths (P1-9).
            umt = iBuildOffsetBaselineFixture();

            out = iRun(testCase, umt);
            entryName = fieldnames(out.data);
            auc = double(out.data.(entryName{1}).value(1, 4, 1));

            % Window is frames 6:12 with values baseline+[1 2 3 2 1 0 0] and a
            % constant baseline, so the correct area is trapz of the response
            % shape alone.
            testCase.verifyEqual(auc, trapz([1 2 3 2 1 0 0]), 'AbsTol', 1e-5);
        end

        function testAUCOmitsNaNSamples(testCase)
            % Every neighbouring reduction omits NaN; TRAPZ did not, so one
            % masked sample turned AUCamplitude into NaN while the other
            % features stayed finite (P1-9).
            umt = iBuildOffsetBaselineFixture();
            entryName = fieldnames(umt.data);
            umt.data.(entryName{1}).value(1, 9, 1) = NaN;

            out = iRun(testCase, umt);
            outName = fieldnames(out.data);
            auc = double(out.data.(outName{1}).value(1, 4, 1));

            testCase.verifyTrue(isfinite(auc));
        end

        function testRejectBaselineThatConsumesRecording(testCase)
            umt = iBuildValidFixture();
            umt.eventInfo.baselinePeriod = 2;

            testCase.verifyError(@() iRun(testCase, umt), ...
                'Umitoolbox:calculateResponseFeatures:InvalidTimeWindow');
        end

        % -----------------------------------------------------------------
        % First trial whose baseline was filled by split_data_by_event
        % -----------------------------------------------------------------
        function testFilledFirstTrialTakesItsStdFromTheOtherTrialsOfTheCondition(testCase)
            % Trial 1's first 6 baseline frames are identical (filled). Its own
            % std is about 0.105, the mean std of trials 2-3 is about 0.234, so
            % the onset (first response frame above the threshold) moves from
            % frame 3 (0.2 s) to frame 5 (0.4 s).
            umt = iBuildFilledFirstTrialFixture([1;1;1], []);

            [out, warnId, warnMsg] = iRunWarn(testCase, umt);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');
            testCase.verifySubstring(warnMsg, 'replaced by the mean of the 2 other selected trial(s)');
            testCase.verifyEqual(iOnsetLatency(out, 1), single([0.4; 0.4]), 'AbsTol', single(1e-5));
        end

        function testFilledFirstTrialWithoutAnotherTrialKeepsItsStdAndWarns(testCase)
            % Trials 2-3 belong to another condition: nothing to take the std from.
            umt = iBuildFilledFirstTrialFixture([1;2;2], []);

            [out, warnId, warnMsg] = iRunWarn(testCase, umt);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');
            testCase.verifySubstring(warnMsg, 'at least 2 are needed');
            testCase.verifyEqual(iOnsetLatency(out, 1), single([0.2; 0.2]), 'AbsTol', single(1e-5));
        end

        function testIgnoredTrialsAreNotUsedForTheFilledTrialStd(testCase)
            % Trial 2 is ignored: only trial 3 (std about 0.313) counts, so the
            % threshold is higher still and the onset moves to frame 7 (0.6 s).
            umt = iBuildFilledFirstTrialFixture([1;1;1], [true; false; true]);

            [out, warnId, warnMsg] = iRunWarn(testCase, umt);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');
            testCase.verifySubstring(warnMsg, 'the 1 other selected trial(s)');
            testCase.verifyEqual(iOnsetLatency(out, 1), single([0.6; 0.6]), 'AbsTol', single(1e-5));
        end

        function testOnsetSecGivesTheExactNumberOfFilledFrames(testCase)
            % First onset at 0.9 s = frame 9 with a 10-frame baseline: exactly
            % 1 frame was filled, whatever the leading values look like.
            umt = iBuildFilledFirstTrialFixture([1;1;1], []);
            umt.eventInfo.onsetSec = [0.9; 3; 6];

            [out, warnId, warnMsg] = iRunWarn(testCase, umt);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');
            testCase.verifySubstring(warnMsg, 'The first 1 baseline frames');
            testCase.verifyEqual(iOnsetLatency(out, 1), single([0.4; 0.4]), 'AbsTol', single(1e-5));
        end

        function testOnsetSecAtTheBaselineMeansNothingWasFilled(testCase)
            % First onset at 1 s with a 1 s baseline: the trial is complete, so
            % a constant run at the start of the baseline is real data. It is
            % not corrected (own std 0.105: onset at frame 3, 0.2 s).
            umt = iBuildFilledFirstTrialFixture([1;1;1], []);
            umt.eventInfo.onsetSec = [1; 3; 6];

            [out, warnId] = iRunWarn(testCase, umt);

            testCase.verifyNotEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');
            testCase.verifyEqual(iOnsetLatency(out, 1), single([0.2; 0.2]), 'AbsTol', single(1e-5));
        end

        function testTrialsOtherThanTheFirstAreNotAffectedByTheCorrection(testCase)
            umt = iBuildFilledFirstTrialFixture([1;1;1], []);

            out = iRunWarn(testCase, umt);

            % Trials 2 and 3 use their own std (0.156 and 0.313): the first
            % response frame above 0.156 is frame 4 (0.3 s), above 0.313 is 7 (0.6 s).
            testCase.verifyEqual(iOnsetLatency(out, 2), single([0.3; 0.3]), 'AbsTol', single(1e-5));
            testCase.verifyEqual(iOnsetLatency(out, 3), single([0.6; 0.6]), 'AbsTol', single(1e-5));
        end

        function testBaselineWithoutAFilledRunRaisesNoWarning(testCase)
            % A noisy first baseline and an all-constant baseline are both left alone.
            umt = iBuildFilledFirstTrialFixture([1;1;1], []);
            name = fieldnames(umt.data);
            umt.data.(name{1}).value(:, 1:6, 1) = repmat(single([0.3 -0.1 0.2 -0.3 0.1 -0.2]), 2, 1);
            [~, warnId] = iRunWarn(testCase, umt);
            testCase.verifyNotEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');

            [~, warnId] = iRunWarn(testCase, iBuildValidFixture());
            testCase.verifyNotEqual(warnId, 'Umitoolbox:calculateResponseFeatures:firstTrialBaselineFilled');
        end

        function testRepetitionsRemainSeparateEventInstances(testCase)
            umt = iBuildValidFixture();
            umt = appendUMTEventInfo(umt, ...
                'eventID', [1;1], ...
                'repetitionIndex', [1;2], ...
                'eventName', {'A';'A'}, ...
                'eventAxisMode', 'instances', ...
                'overwrite', true);
            umt.eventInfo.baselinePeriod = 0.5;

            out = iRun(testCase, umt);
            entryName = fieldnames(out.data);
            values = out.data.(entryName{1}).value;

            testCase.verifyEqual(size(values, 3), 2);
            testCase.verifyEqual( ...
                squeeze(values(1, 1, :)), single([3;4]), ...
                'AbsTol', single(10 * eps('single')));
        end
    end
end

function out = iRun(testCase, umt, varargin)
%IRUN Save a UMT struct as a .umt file and run the function on that file.
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
umtFile = fullfile(folder, 'responseInput.umt');
saveData(umtFile, umt);
out = calculateResponseFeatures(umtFile, varargin{:});
end

function [out, warnId, warnMsg] = iRunWarn(testCase, umt)
%IRUNWARN iRun that also returns the last warning raised by the call.
lastwarn('', '');
out = iRun(testCase, umt);
[warnMsg, warnId] = lastwarn;
end

function lat = iOnsetLatency(out, iE)
%IONSETLATENCY OnsetLatency (measure 6) of every ROI for event instance IE.
name = fieldnames(out.data);
lat = out.data.(name{1}).value(:, 6, iE);
end

function umt = iBuildFilledFirstTrialFixture(eventID, selected)
%IBUILDFILLEDFIRSTTRIALFIXTURE 2 ROI x 30 T x 3 E at 10 Hz, baseline 1 s (10 frames).
%
%   Trial 1 has its first 6 baseline frames filled with the mean (0) of its
%   last 4 real ones, as split_data_by_event does. Trials 2 and 3 have noisy
%   baselines (std 0.156 and 0.313). The response is a ramp of 0.05 per frame.
%   ROI 2 is twice ROI 1.

frameRateHz = 10;
base2 = single([0.1 -0.1 0.2 -0.2 0.1 -0.1 0.2 -0.2 0.1 -0.1]);
base1 = [zeros(1, 6, 'single'), single([0.1 -0.1 0.2 -0.2])];
resp = single(0.05 * (1:20));

row = zeros(1, 30, 3, 'single');
row(1, :, 1) = [base1, resp];
row(1, :, 2) = [base2, resp];
row(1, :, 3) = [2 * base2, resp];
vals = [row; 2 * row];

umt = genUMTStruct(vals, ...
    'kind', 'roi', ...
    'entryName', 'main', ...
    'dimNames', {'ROI','T','E'}, ...
    'labels', struct('ROI', {{'ROI1','ROI2'}}));
umt.data.main.meta = struct('FrameRateHz', frameRateHz);

names = repmat({'A'}, 3, 1);
names(eventID == 2) = {'B'};
umt = appendUMTEventInfo(umt, ...
    'eventID', eventID(:), ...
    'repetitionIndex', [1; 2; 3], ...
    'eventName', names, ...
    'eventAxisMode', 'instances', ...
    'overwrite', true);
umt.eventInfo.baselinePeriod = 1;
if ~isempty(selected)
    umt.eventInfo.selected = logical(selected(:));
end
end

function umt = iBuildValidFixture()
frameRateHz = 10;
baselinePeriod = 0.5;

t = 12;
e = 2;
roi = 2;
vals = zeros(roi, t, e, 'single');
vals(1,:,1) = [0 0 0 0 0 1 2 3 2 1 0 0];
vals(1,:,2) = [0 0 0 0 0 2 3 4 3 2 0 0];
vals(2,:,1) = [0 0 0 0 0 0.5 1 1.5 1 0.5 0 0];
vals(2,:,2) = [0 0 0 0 0 1 1.5 2 1.5 1 0 0];

labels = struct();
labels.ROI = {'ROI1','ROI2'};

umt = genUMTStruct(vals, ...
    'kind', 'roi', ...
    'entryName', 'main', ...
    'dimNames', {'ROI','T','E'}, ...
    'labels', labels);

umt.data.main.meta = struct('FrameRateHz', frameRateHz);

umt = appendUMTEventInfo(umt, ...
    'eventID', [1;2], ...
    'repetitionIndex', [1;1], ...
    'eventName', {'A';'B'}, ...
    'eventAxisMode', 'instances', ...
    'overwrite', true);

umt.eventInfo.baselinePeriod = baselinePeriod;
end

function umt = iBuildOffsetBaselineFixture()
%IBUILDOFFSETBASELINEFIXTURE Same shape as the valid fixture, baseline != 0.
%
%   A non-zero baseline is what separates the correct AUC from the old
%   difference-of-integrals: with an all-zero baseline both formulas agree.

umt = iBuildValidFixture();
entryName = fieldnames(umt.data);
umt.data.(entryName{1}).value = umt.data.(entryName{1}).value + 2;
end

function umt = iBuildPixelFixture()
frameRateHz = 10;
baselinePeriod = 0.5;
vals = rand(2, 3, 12, 2, 'single');
labels = struct();
labels.ROI = {'ROI1','ROI2'};
labels.Pixel = {'1','2','3'};

umt = genUMTStruct(vals, ...
    'kind', 'roi', ...
    'entryName', 'main', ...
    'dimNames', {'ROI','Pixel','T','E'}, ...
    'labels', labels);

umt.data.main.meta = struct('FrameRateHz', frameRateHz);

umt = appendUMTEventInfo(umt, ...
    'eventID', [1;2], ...
    'repetitionIndex', [1;1], ...
    'eventName', {'A';'B'}, ...
    'eventAxisMode', 'instances', ...
    'overwrite', true);

umt.eventInfo.baselinePeriod = baselinePeriod;
end
