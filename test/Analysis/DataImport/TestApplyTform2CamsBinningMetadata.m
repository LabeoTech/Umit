classdef TestApplyTform2CamsBinningMetadata < matlab.unittest.TestCase
    %TESTAPPLYTFORM2CAMSBINNINGMETADATA Effective-binning regression tests.

    properties
        TempFolder char
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fixture = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fixture.Folder;
        end
    end

    methods (Test)
        function testUsesHardwareTimesClassificationBinning(testCase)
            testCase.writeDualCameraFolder(4, 1);

            tform = affine2d([1 0 0; 0 1 0; 4 -2 1]);
            tformInfo = struct( ...
                'Binning', 2, ...
                'BinningSpatial', 1, ...
                'Rotation', 0, ...
                'X_Offset', 0, ...
                'Y_Offset', 0);

            [status, warnmsg] = applyTform2Cams( ...
                testCase.TempFolder, tform, tformInfo);

            testCase.verifyTrue(status, warnmsg);
            saved = load(fullfile(testCase.TempFolder, 'tformDualCam.mat'), 'tform');
            expected = [1 0 0; 0 1 0; 2 -1 1];
            testCase.verifyEqual(saved.tform.T, expected, 'AbsTol', 1e-12);
        end

        function testProcessesHeaderedCameraTwoFile(testCase)
            % .dat header Phase 4e-1: a headered Camera-2 file is
            % coregistered and stays headered.
            testCase.writeDualCameraFolder(4, 1);
            iRewriteWithHeader(fullfile(testCase.TempFolder, 'green.dat'), [8 8 1], 1);
            tform = affine2d([1 0 0; 0 1 0; 4 -2 1]);
            tformInfo = struct('Binning', 2, 'BinningSpatial', 1, 'Rotation', 0, ...
                'X_Offset', 0, 'Y_Offset', 0);

            [status, warnmsg] = applyTform2Cams(testCase.TempFolder, tform, tformInfo);

            testCase.verifyTrue(status, warnmsg);
            f = fullfile(testCase.TempFolder, 'green.dat');
            testCase.verifyTrue(isDatWithHeader(f));
            hdr = readDatHeader(f);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.dimSizes, [8 8 1]);
            testCase.verifyEqual(hdr.channelName, 'green');
            testCase.verifyFalse(isfile([f '.tmp']));
        end

        function testRejectsLegacyAcquisitionMetadata(testCase)
            testCase.writeDualCameraFolder(2, []);
            tform = affine2d(eye(3));
            tformInfo = struct('Binning', 2, 'BinningSpatial', 1);

            testCase.verifyError( ...
                @() applyTform2Cams(testCase.TempFolder, tform, tformInfo), ...
                'Umitoolbox:applyTform2Cams:MissingBinningMetadata');
        end

        function testRejectsLegacyTransformMetadata(testCase)
            testCase.writeDualCameraFolder(2, 1);
            tform = affine2d(eye(3));
            tformInfo = struct('Binning', 2);

            testCase.verifyError( ...
                @() applyTform2Cams(testCase.TempFolder, tform, tformInfo), ...
                'Umitoolbox:applyTform2Cams:MissingBinningMetadata');
        end

        function testGenerationRejectsLegacyAcquisitionMetadata(testCase)
            testCase.writeDualCameraFolder(2, []);

            testCase.verifyError( ...
                @() genTform2Cams(testCase.TempFolder, false, false), ...
                'Umitoolbox:genTform2Cams:MissingBinningMetadata');
        end
    end

    methods (Access = private)
        function writeDualCameraFolder(testCase, hardwareBinning, softwareBinning)
            frameSize = [8, 8];
            AcqInfoStream = struct(); %#ok<NASGU>
            AcqInfoStream.Width = frameSize(2);
            AcqInfoStream.Height = frameSize(1);
            AcqInfoStream.Length = 1;
            AcqInfoStream.FrameRateHz = 1;
            AcqInfoStream.Datatype = 'single';
            AcqInfoStream.MultiCam = 1;
            AcqInfoStream.Binning = hardwareBinning;
            if ~isempty(softwareBinning)
                AcqInfoStream.BinningSpatial = softwareBinning;
            end
            AcqInfoStream.Illumination1 = struct( ...
                'Color', 'red', 'CamIdx', 1, 'FrameIdx', 1);
            AcqInfoStream.Illumination2 = struct( ...
                'Color', 'green', 'CamIdx', 2, 'FrameIdx', 1);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            testCase.writeDatFile('red.dat', zeros(frameSize, 'single'), AcqInfoStream.FrameRateHz);
            testCase.writeDatFile('green.dat', zeros(frameSize, 'single'), AcqInfoStream.FrameRateHz);
        end

        function writeDatFile(testCase, fileName, data, frameRateHz)
            % Headered input (.dat header Phase 5a).
            writeTestDat(fullfile(testCase.TempFolder, fileName), single(data), frameRateHz);
        end
    end
end

function iRewriteWithHeader(filePath, dimSizes, rate)
% Put a v1 header in front of the file's existing single values
% (.dat header Phase 4c-1 guard test).
values = loadData(filePath);
[~, base] = fileparts(filePath);
hdr = struct('dataClass', 'single', 'frameRateHz', rate, 'exposureMsec', NaN, ...
    'channelName', base, 'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', dimSizes, ...
    'writeComplete', true);
fid = fopen(filePath, 'w', 'ieee-le');
fwrite(fid, encodeDatHeader(hdr), 'uint8');
fwrite(fid, values, 'single');
fclose(fid);
end
