classdef TestHemoCorrection < matlab.unittest.TestCase
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
                'ChannelList', {'red'}, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifySize(out, size(dataYXT));
        end

        function testRawDatInputResamplesReferenceChannel(testCase)
            out = HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'ChannelList', {'red'});

            testCase.verifyEqual(out, 'fluoHemoCorr.dat');
            outPath = fullfile(testCase.TempFolder, out);
            testCase.verifyTrue(isfile(outPath));

            info = loadMetaData(outPath);
            testCase.verifyEqual(datAxisSize(info, 'T'), testCase.AcqInfo.Length);
            testCase.verifyEqual(info.frameRateHz, testCase.AcqInfo.FrameRateHz);
        end

        function testRawDatInputOverwritesDeclaredOutputOnRerun(testCase)
            % A pre-existing output must be replaced, not sidestepped under a
            % different name: the declared output is fluoHemoCorr.dat, so a
            % re-run has to keep writing there and must not leave a stale
            % original or scratch file behind (DFR-20260819-008 sibling).
            staleFile = fullfile(testCase.TempFolder, 'fluoHemoCorr.dat');
            fid = fopen(staleFile, 'w');
            fclose(fid);
            staleBytes = dir(staleFile).bytes;

            out = HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'ChannelList', {'red'});

            testCase.verifyEqual(out, 'fluoHemoCorr.dat');
            testCase.verifyTrue(isfile(staleFile));
            testCase.verifyGreaterThan(dir(staleFile).bytes, staleBytes);
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'fluoHemoCorr_preallocData.dat')));
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'fluoHemoCorr_writing.dat')));
        end

        function testLowPassStillAcceptsValidCutoff(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            out = HemoCorrection(dataYXT, testCase.TempFolder, ...
                'ChannelList', {'red'}, ...
                'LowPassFreq', testCase.AcqInfo.FrameRateHz / 10, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifySize(out, size(dataYXT));
        end

        function testInvalidLowPassFreq(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'ChannelList', {'red'}, ...
                'LowPassFreq', testCase.AcqInfo.FrameRateHz, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:invalidLowPassFreq');
        end

        function testMissingReferenceChannel(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'ChannelList', {'green'}, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:missingReferenceChannel');
        end

        function testRunsWithoutAcqInfos(testCase)
            % Nothing reads AcqInfos.mat: the frame rates come from the
            % .dat headers (or FrameRateHz), so a folder without it works.
            dataYXT = iLoadNumericData(testCase.TempFolder);
            delete(fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            out = HemoCorrection(dataYXT, testCase.TempFolder, ...
                'ChannelList', {'red'}, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);
            testCase.verifySize(out, size(dataYXT));

            outFile = HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'ChannelList', {'red'});
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, outFile)));
        end

        function testInRamInputWithoutFrameRateErrors(testCase)
            % Numeric input needs 'FrameRateHz'; AcqInfos.mat in the folder
            % is not used for the fluorescence frame rate.
            testCase.assertTrue(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, testCase.TempFolder, ...
                'ChannelList', {'red'}), ...
                'Umitoolbox:HemoCorrection:missingFrameRateHz');
        end


        function testNumericInputRejectsReferenceDurationMismatch(testCase)
            mismatchFolder = fullfile(testCase.TempFolder, 'durationMismatchNumeric');
            mkdir(mismatchFolder);
            mismatchAcqInfo = iCreateDurationMismatchDataset(mismatchFolder);
            dataYXT = iLoadNumericData(mismatchFolder);

            testCase.verifyError(@() HemoCorrection(dataYXT, mismatchFolder, ...
                'ChannelList', {'red'}, ...
                'FrameRateHz', mismatchAcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:DurationMismatch');
        end

        function testRawDatInputRejectsReferenceDurationMismatch(testCase)
            mismatchFolder = fullfile(testCase.TempFolder, 'durationMismatchRaw');
            mkdir(mismatchFolder);
            iCreateDurationMismatchDataset(mismatchFolder);

            testCase.verifyError(@() HemoCorrection('yellow.dat', mismatchFolder, ...
                'ChannelList', {'red'}), ...
                'Umitoolbox:HemoCorrection:DurationMismatch');
        end

        function testPipelineInfo(testCase)
            %#ok<NASGU>
            info = HemoCorrection('pipelineInfo');
            testCase.verifyEqual(info.name, 'HemoCorrection');
            testCase.verifyFalse(info.legacyOpts);
            testCase.verifyEqual(info.outputs(1).type, {'ImageTimeSeries'});
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'fluoHemoCorr.dat');
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
