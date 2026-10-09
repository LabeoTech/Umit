classdef TestCalculateResponseFeatures < matlab.unittest.TestCase
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
