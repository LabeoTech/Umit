classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestHemoCorrection < matlab.unittest.TestCase
    %TESTHEMOCORRECTION Unit tests for HemoCorrection.
    %
    %   These tests use a synthetic AcqInfos-driven dataset with unequal
    %   imported-channel timelines. The fluorescence channel follows the
    %   base timeline, while the red reference channel has twice the number
    %   of frames and twice the frame rate. This exercises the reference
    %   resampling path used by the current metadata model.

    properties
        ProjectRoot
        TempFolder
        AcqInfo
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestHemoCorrection_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            testCase.AcqInfo = iCreateRepeatedTimelineDataset(testCase.TempFolder);
        end
    end

    methods (TestMethodTeardown)
        function removeTempFolder(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                try %#ok<TRYNC>
                    rmdir(testCase.TempFolder, 's');
                end
            end
        end
    end

    methods (Test)
        function testNumericInputResamplesReferenceChannel(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            out = HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifySize(out, size(dataYXT));
        end

        function testRawDatInputResamplesReferenceChannel(testCase)
            out = HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false);

            testCase.verifyEqual(out, 'hemoCorr_fluo.dat');
            outPath = fullfile(testCase.TempFolder, out);
            testCase.verifyTrue(isfile(outPath));

            info = loadMetaData(outPath);
            testCase.verifyEqual(datAxisSize(info, 'T'), testCase.AcqInfo.Length);
            testCase.verifyEqual(info.frameRateHz, testCase.AcqInfo.FrameRateHz);
        end

        function testLowPassStillAcceptsValidCutoff(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            out = HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'LowPassFreq', testCase.AcqInfo.FrameRateHz / 10, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifySize(out, size(dataYXT));
        end

        function testInvalidLowPassFreq(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'LowPassFreq', testCase.AcqInfo.FrameRateHz, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:invalidLowPassFreq');
        end

        function testLowPassFreqFiltersTheReference(testCase)
            % The reference and the fluorescence share a 4 Hz component.
            % A 1 Hz cutoff removes it from the reference before the
            % regression, so it is no longer regressed out of the
            % fluorescence.
            folder = fullfile(testCase.TempFolder, 'lowPass');
            mkdir(folder);
            [dataYXT, frameRateHz] = iCreateLowPassDataset(folder);
            args = {'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', frameRateHz};

            outNoLP = HemoCorrection(dataYXT, folder, args{:});
            outLP = HemoCorrection(dataYXT, folder, args{:}, 'LowPassFreq', 1);

            amp4Hz = @(d) iComponentAmplitude(d, 4, frameRateHz);
            testCase.verifyGreaterThan(amp4Hz(outLP), 0.9 * amp4Hz(dataYXT), ...
                'The 4 Hz component must survive a 1 Hz reference cutoff.');
            testCase.verifyLessThan(amp4Hz(outNoLP), amp4Hz(outLP), ...
                'Without LowPassFreq the 4 Hz component is regressed out.');
        end

        function testLowPassFreqGivesSameResultForArrayAndDat(testCase)
            folder = fullfile(testCase.TempFolder, 'lowPassModes');
            mkdir(folder);
            [dataYXT, frameRateHz] = iCreateLowPassDataset(folder);

            fromArray = HemoCorrection(dataYXT, folder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', frameRateHz, 'LowPassFreq', 1);
            outFile = HemoCorrection('yellow.dat', folder, ...
                'Red', true, 'Green', false, 'Amber', false, 'LowPassFreq', 1);
            fromDat = loadData(fullfile(folder, outFile));

            testCase.verifyEqual(fromDat, fromArray, 'AbsTol', single(1e-6));
        end

        function testLowPassFreqAtReferenceNyquistWarnsAndSkips(testCase)
            % A 5 Hz reference (Nyquist 2.5 Hz) cannot be filtered at 3 Hz,
            % although 3 Hz is below the 10 Hz fluorescence Nyquist (5 Hz):
            % the filter is skipped with a warning, as with LowPassFreq 0.
            folder = fullfile(testCase.TempFolder, 'slowReference');
            mkdir(folder);
            [dataYXT, frameRateHz] = iCreateLowPassDataset(folder, 0.5);
            args = {'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', frameRateHz};

            outLP = testCase.verifyWarning( ...
                @() HemoCorrection(dataYXT, folder, args{:}, 'LowPassFreq', 3), ...
                'Umitoolbox:HemoCorrection:LowPassFreqSkipped');
            outNoLP = HemoCorrection(dataYXT, folder, args{:});

            testCase.verifyEqual(outLP, outNoLP);
        end

        function testLowPassFreqBelowReferenceNyquistDoesNotWarn(testCase)
            folder = fullfile(testCase.TempFolder, 'fastReference');
            mkdir(folder);
            [dataYXT, frameRateHz] = iCreateLowPassDataset(folder);

            testCase.verifyWarningFree(@() HemoCorrection(dataYXT, folder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', frameRateHz, 'LowPassFreq', 1));
        end

        function testLowPassFreqRejectedForRatiometric(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'LowPassFreq', 1, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:UnsupportedLowPassFreq');
        end

        function testLowPassFreqRejectedForEventSplitData(testCase)
            dataYXTE = repmat(iLoadNumericData(testCase.TempFolder), 1, 1, 1, 2);

            testCase.verifyError(@() HemoCorrection(dataYXTE, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'LowPassFreq', 1, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:UnsupportedLowPassFreq');
        end

        function testNegativeLowPassFreqRejected(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'LowPassFreq', -1, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end

        function testLowPassFreqDeclaredInPipelineInfo(testCase)
            info = HemoCorrection('pipelineInfo');
            lp = info.parameters(strcmp({info.parameters.name}, 'LowPassFreq'));

            testCase.assertNumElements(lp, 1);
            testCase.verifyEqual(lp.default, 0);
            testCase.verifySubstring(char(string(lp.description)), 'LinearRegression');
            testCase.verifySubstring(char(string(lp.description)), 'Ratiometric');
        end

        function testNoReferenceChannelErrors(testCase)
            % LinearRegression without any reference channel raises a clear
            % error for both input modes (no channel-selection dialog).
            dataYXT = iLoadNumericData(testCase.TempFolder);
            noRef = {'Red', false, 'Green', false, 'Amber', false};

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                noRef{:}, 'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:InvalidInput');
            testCase.verifyError(@() HemoCorrection('yellow.dat', testCase.TempFolder, ...
                noRef{:}), ...
                'Umitoolbox:HemoCorrection:InvalidInput');
        end

        function testMissingReferenceChannel(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', false, 'Green', true, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:missingReferenceChannel');
        end

        function testRunsWithoutAcqInfos(testCase)
            % Nothing reads AcqInfos.mat: the frame rates come from the
            % .dat headers (or FrameRateHz), so a folder without it works.
            dataYXT = iLoadNumericData(testCase.TempFolder);
            delete(fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            out = HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);
            testCase.verifySize(out, size(dataYXT));

            outFile = HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false);
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, outFile)));
        end

        function testInRamInputWithoutFrameRateErrors(testCase)
            % Numeric input needs 'FrameRateHz'; AcqInfos.mat in the folder
            % is not used for the fluorescence frame rate.
            testCase.assertTrue(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Red', true, 'Green', false, 'Amber', false), ...
                'Umitoolbox:HemoCorrection:missingFrameRateHz');
        end


        function testNumericInputRejectsReferenceDurationMismatch(testCase)
            mismatchFolder = fullfile(testCase.TempFolder, 'durationMismatchNumeric');
            mkdir(mismatchFolder);
            mismatchAcqInfo = iCreateDurationMismatchDataset(mismatchFolder);
            dataYXT = iLoadNumericData(mismatchFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, mismatchFolder, ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', mismatchAcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:DurationMismatch');
        end

        function testRawDatInputRejectsReferenceDurationMismatch(testCase)
            mismatchFolder = fullfile(testCase.TempFolder, 'durationMismatchRaw');
            mkdir(mismatchFolder);
            iCreateDurationMismatchDataset(mismatchFolder);

            testCase.verifyError(@() HemoCorrection('yellow.dat', mismatchFolder, ...
                'Red', true, 'Green', false, 'Amber', false), ...
                'Umitoolbox:HemoCorrection:DurationMismatch');
        end
    end
end


function AcqInfoStream = iCreateDurationMismatchDataset(saveFolder)
%ICREATEDURATIONMISMATCHDATASET Create valid files with incompatible durations.
%
% The red channel is written with metadata that matches the file on disk,
% but its Length/FrameRateHz duration is shorter than the fluorescence
% channel. This should trigger the duration/span assertion before
% resampling.

Ny = 12;
Nx = 10;
Nt = 16;
baseFreq = 10;
redNt = 2 * Nt - 2;
redFreq = 2 * baseFreq;

[y, x, t] = ndgrid(single(1:Ny), single(1:Nx), single(1:Nt));
yellow = 1 + 0.01 * sin(2*pi*t/Nt) + 0.001 * y + 0.0015 * x;
yellow = single(yellow);

[y2, x2, t2] = ndgrid(single(1:Ny), single(1:Nx), single(1:redNt));
red = 1 + 0.012 * cos(2*pi*t2/max(redNt, 1)) + 0.001 * y2 + 0.0015 * x2;
red = single(red);

iWriteDat(fullfile(saveFolder, 'yellow.dat'), yellow, baseFreq, 5);
iWriteDat(fullfile(saveFolder, 'red.dat'), red, redFreq, 5);

AcqInfoStream = struct();
AcqInfoStream.Height = Ny;
AcqInfoStream.Width = Nx;
AcqInfoStream.Length = Nt;
AcqInfoStream.FrameRateHz = baseFreq;
AcqInfoStream.Datatype = 'single';
AcqInfoStream.BinningSpatial = 1;
AcqInfoStream.BinningTemp = 1;
AcqInfoStream.MultiCam = false;
AcqInfoStream.ExposureMsec = 5;
AcqInfoStream.ImportedChannels = struct( ...
    'DatFile', {'yellow.dat', 'red.dat'}, ...
    'Length', {Nt, redNt}, ...
    'FrameRateHz', {baseFreq, redFreq}, ...
    'ExposureMsec', {5, 5}, ...
    'CamIdx', {1, 1});

save(fullfile(saveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
end


function [yellow, frameRateHz] = iCreateLowPassDataset(saveFolder, refRateFactor)
%ICREATELOWPASSDATASET Fluorescence and a reference sharing a 4 Hz component.
%
% The reference runs at REFRATEFACTOR x the 10 Hz fluorescence rate
% (default 2). The fluorescence also carries a 0.5 Hz signal that the
% reference lacks. Both channels span the same 4 s recording.
if nargin < 2
    refRateFactor = 2;
end
Ny = 6;
Nx = 5;
Nt = 40;
frameRateHz = 10;
refRateHz = refRateFactor * frameRateHz;
refNt = round(refRateFactor * Nt);

[y, x, t] = ndgrid(single(1:Ny), single(1:Nx), single(0:Nt-1) / frameRateHz);
yellow = single(1 + 0.02 * sin(2*pi*4*t) + 0.01 * sin(2*pi*0.5*t) ...
    + 0.001 * y + 0.001 * x);

[y2, x2, t2] = ndgrid(single(1:Ny), single(1:Nx), single(0:refNt-1) / refRateHz);
red = single(1 + 0.02 * sin(2*pi*4*t2) + 0.001 * y2 + 0.001 * x2);

iWriteDat(fullfile(saveFolder, 'yellow.dat'), yellow, frameRateHz, 5);
iWriteDat(fullfile(saveFolder, 'red.dat'), red, refRateHz, 5);
end

function amp = iComponentAmplitude(data, freqHz, frameRateHz)
%ICOMPONENTAMPLITUDE Amplitude of one frequency in the spatially averaged trace.
nT = size(data, 3);
trace = double(mean(reshape(data, [], nT), 1));
t = (0:nT-1) / frameRateHz;
amp = abs(trace * exp(-1i * 2 * pi * freqHz * t).') * 2 / nT;
end

function dataYXT = iLoadNumericData(saveFolder)
loaded = loadData(fullfile(saveFolder, 'yellow.dat'));
assert(isnumeric(loaded) && ndims(loaded) == 3, ...
    'Fixture data must resolve to numeric YXT data.');
dataYXT = single(loaded);
end

function AcqInfoStream = iCreateRepeatedTimelineDataset(saveFolder)
Ny = 12;
Nx = 10;
Nt = 16;
baseFreq = 10;

[y, x, t] = ndgrid(single(1:Ny), single(1:Nx), single(1:Nt));
yellow = 1 + 0.01 * sin(2*pi*t/Nt) + 0.001 * y + 0.0015 * x;
yellow = single(yellow);

[y2, x2, t2] = ndgrid(single(1:Ny), single(1:Nx), single(1:(2*Nt)));
red = 1 + 0.012 * cos(2*pi*t2/(2*Nt)) + 0.001 * y2 + 0.0015 * x2;
red = single(red);

iWriteDat(fullfile(saveFolder, 'yellow.dat'), yellow, baseFreq, 5);
iWriteDat(fullfile(saveFolder, 'red.dat'), red, 2*baseFreq, 5);

AcqInfoStream = struct();
AcqInfoStream.Height = Ny;
AcqInfoStream.Width = Nx;
AcqInfoStream.Length = Nt;
AcqInfoStream.FrameRateHz = baseFreq;
AcqInfoStream.Datatype = 'single';
AcqInfoStream.BinningSpatial = 1;
AcqInfoStream.BinningTemp = 1;
AcqInfoStream.MultiCam = false;
AcqInfoStream.ExposureMsec = 5;
AcqInfoStream.ImportedChannels = struct( ...
    'DatFile', {'yellow.dat', 'red.dat'}, ...
    'Length', {Nt, 2*Nt}, ...
    'FrameRateHz', {baseFreq, 2*baseFreq}, ...
    'ExposureMsec', {5, 5}, ...
    'CamIdx', {1, 1});

save(fullfile(saveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
end

function iWriteDat(filePath, data, frameRateHz, exposureMsec)
% Headered input (.dat header Phase 5a).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
end
