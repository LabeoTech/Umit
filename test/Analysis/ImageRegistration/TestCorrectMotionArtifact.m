classdef TestCorrectMotionArtifact < matlab.unittest.TestCase
    %TESTCORRECTMOTIONARTIFACT Unit tests for correctMotionArtifact.
    %
    % The function returns the corrected data (array in -> array out, .dat in
    % -> motionCorrected.dat). The estimated shifts are only available in the
    % <name>_MotionCorrection.mat that SaveShifts=true writes, so the shift
    % tests read them from there.

    properties
        TempFolder
        Ny = 64
        Nx = 64
        Nt = 4
        RefPattern
        Shifts    % Nt-by-2 [rowShift colShift] estimated for green.dat
        RedShifts % Nt-by-2 [rowShift colShift] estimated for red.dat
    end

    methods (TestMethodSetup)
        function setup(testCase)
            testCase.TempFolder = fullfile(tempdir, ...
                ['TestCorrectMotionArtifact_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            % Textured pattern: the registration high-pass filter removes
            % large smooth structure, so the test image needs fine detail.
            stream = RandStream('mt19937ar', 'Seed', 1);
            texture = imgaussfilt(randn(stream, testCase.Ny + 40, testCase.Nx + 40), 1.5);
            texture = texture(21:end-20, 21:end-20);
            texture = (texture - min(texture(:))) / (max(texture(:)) - min(texture(:)));
            testCase.RefPattern = single(texture);

            testCase.Shifts = [0 0; 2 -3; -1 4; 3 1];
            testCase.RedShifts = [0 0; -2 1; 3 2; 1 -4]; % different motion

            greenData = iBuildShiftedStack(testCase.RefPattern, testCase.Shifts);
            redData = iBuildShiftedStack(testCase.RefPattern * 0.6, testCase.RedShifts);

            iWriteSingleDat(fullfile(testCase.TempFolder, 'green.dat'), greenData, 10, 20);
            iWriteSingleDat(fullfile(testCase.TempFolder, 'red.dat'), redData, 10, 20);

            AcqInfoStream = struct('Width', testCase.Nx, 'Height', testCase.Ny, ...
                'Length', testCase.Nt, 'FrameRateHz', 10, 'ExposureMsec', 20);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testEstimatesKnownShifts(testCase)
            [~, params] = testCase.runAndReadShifts('green.dat', 'UpsamplingFactor', 10);

            testCase.verifyEqual(size(params), [testCase.Nt, 4]);
            testCase.verifyEqual(params(1,:), [0 0 0 1]);
            % Columns are [tx ty rotationDeg scale]; tx is the column shift.
            testCase.verifyEqual(params(:, [2 1]), testCase.Shifts, 'AbsTol', 0.15);
            testCase.verifyEqual(params(:, 3), zeros(testCase.Nt, 1));
            testCase.verifyEqual(params(:, 4), ones(testCase.Nt, 1));
        end

        function testChannelsAreCorrectedIndependently(testCase)
            % Each recording is registered against its own first frame, so two
            % channels with different motion must get different transforms.
            [~, greenParams] = testCase.runAndReadShifts('green.dat', 'UpsamplingFactor', 10);
            [~, redParams] = testCase.runAndReadShifts('red.dat', 'UpsamplingFactor', 10);

            testCase.verifyEqual(greenParams(:, [2 1]), testCase.Shifts, 'AbsTol', 0.15);
            testCase.verifyEqual(redParams(:, [2 1]), testCase.RedShifts, 'AbsTol', 0.15);
        end

        function testDatOutputIsTheDeclaredFileAndSourcesAreUntouched(testCase)
            greenBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat'));
            redBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat'));

            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'motionCorrected.dat'));
            testCase.verifyTrue(isfile(outFile));
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat')), greenBefore);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat')), redBefore);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));

            info = loadMetaData(outFile);
            testCase.verifyEqual(info.dimNames, {'Y','X','T'});
            testCase.verifyEqual(double(info.dimSizes(:)).', [testCase.Ny testCase.Nx testCase.Nt]);
            testCase.verifyEqual(info.dataClass, 'single');
        end

        function testRerunReplacesTheDeclaredOutput(testCase)
            staleFile = fullfile(testCase.TempFolder, 'motionCorrected.dat');
            fclose(fopen(staleFile, 'w'));
            staleBytes = dir(staleFile).bytes;

            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(outFile, staleFile);
            testCase.verifyGreaterThan(dir(staleFile).bytes, staleBytes);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));
        end

        function testCorrectionRealignsToFirstFrame(testCase)
            greenCorrected = loadData(correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'UpsamplingFactor', 10, 'SaveShifts', false, 'ShowPlot', false));
            redCorrected = loadData(correctMotionArtifact('red.dat', testCase.TempFolder, ...
                'UpsamplingFactor', 10, 'SaveShifts', false, 'ShowPlot', false));

            % Interior region only: a shift zero-fills pixels warped in from
            % outside the original frame near the border.
            interior = 10:(testCase.Ny-10);
            for t = 1:testCase.Nt
                testCase.verifyEqual( ...
                    greenCorrected(interior, interior, t), ...
                    testCase.RefPattern(interior, interior), 'AbsTol', 0.1);
                testCase.verifyEqual( ...
                    redCorrected(interior, interior, t), ...
                    testCase.RefPattern(interior, interior) * 0.6, 'AbsTol', 0.1);
            end
        end

        function testArrayInputReturnsASameSizeArrayEqualToTheDatResult(testCase)
            dataYXT = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));

            outArray = correctMotionArtifact(dataYXT, testCase.TempFolder, ...
                'UpsamplingFactor', 10, 'SaveShifts', false, 'ShowPlot', false);
            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'UpsamplingFactor', 10, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyClass(outArray, 'single');
            testCase.verifyEqual(size(outArray), size(dataYXT));
            testCase.verifyEqual(outArray, single(loadData(outFile)));
        end

        function testSimilarityTransformOnTranslationOnlyMotion(testCase)
            [~, params] = testCase.runAndReadShifts('green.dat', 'TransformType', 'similarity');

            testCase.verifyEqual(size(params), [testCase.Nt, 4]);
            testCase.verifyEqual(params(1,:), [0 0 0 1]);
            % Motion is pure translation: no meaningful rotation or scale.
            testCase.verifyLessThan(max(abs(params(:, 3))), 1);
            testCase.verifyEqual(params(:, 4), ones(testCase.Nt, 1), 'AbsTol', 0.05);
            testCase.verifyEqual(params(:, [2 1]), testCase.Shifts, 'AbsTol', 1);
        end

        function testSaveShiftsWritesProvenanceAndQC(testCase)
            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, 'ShowPlot', false);

            provPath = fullfile(testCase.TempFolder, 'green_MotionCorrection.mat');
            testCase.verifyTrue(isfile(provPath));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'green_MotionCorrectionQC.png')));
            S = load(provPath, 'MotionCorrection');
            testCase.verifyEqual(S.MotionCorrection.inputFile, 'green.dat');
            testCase.verifyEqual(S.MotionCorrection.transformType, 'translation');
            testCase.verifyEqual(S.MotionCorrection.interpolation, 'cubic');
            testCase.verifyEqual(S.MotionCorrection.correctedFile, outFile);
            testCase.verifyEqual(size(S.MotionCorrection.params), [testCase.Nt, 4]);
            testCase.verifyEqual(size(S.MotionCorrection.corrBefore), [testCase.Nt, 1]);
            % Correction must raise the correlation with the first frame.
            testCase.verifyGreaterThan(mean(S.MotionCorrection.corrAfter(2:end)), ...
                mean(S.MotionCorrection.corrBefore(2:end)));
        end

        function testNothingIsSavedWhenSaveShiftsIsFalse(testCase)
            correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_MotionCorrection.mat')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_MotionCorrectionQC.png')));
        end

        function testArrayInputSavesShiftsUnderAFixedName(testCase)
            dataYXT = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));

            correctMotionArtifact(dataYXT, testCase.TempFolder, 'ShowPlot', false);

            S = load(fullfile(testCase.TempFolder, 'MotionCorrection.mat'), 'MotionCorrection');
            testCase.verifyEqual(size(S.MotionCorrection.params), [testCase.Nt, 4]);
            testCase.verifyEqual(S.MotionCorrection.correctedFile, '');
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'MotionCorrectionQC.png')));
        end

        function testProvenanceIsPerFile(testCase)
            correctMotionArtifact('green.dat', testCase.TempFolder, 'ShowPlot', false);
            correctMotionArtifact('red.dat', testCase.TempFolder, 'ShowPlot', false);

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'green_MotionCorrection.mat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'red_MotionCorrection.mat')));
        end

        function testProcessesHeaderedDat(testCase)
            % .dat header Phase 4e-1: a headered file is corrected like a
            % headerless one, and the output is headered.
            iRewriteWithHeader(fullfile(testCase.TempFolder, 'red.dat'), ...
                [testCase.Ny testCase.Nx testCase.Nt], 10);

            [outFile, params] = testCase.runAndReadShifts('red.dat', 'UpsamplingFactor', 10);

            testCase.verifyEqual(params(:, [2 1]), testCase.RedShifts, 'AbsTol', 0.15);
            testCase.verifyTrue(isDatWithHeader(outFile));
            hdr = readDatHeader(outFile);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.frameRateHz, 10);
            testCase.verifyEqual(hdr.channelName, 'motionCorrected');
        end

        function testInputMayLiveOutsideSaveFolder(testCase)
            % The corrected file goes to SaveFolder, so the source no longer
            % has to be there.
            externalFolder = fullfile(testCase.TempFolder, 'external');
            mkdir(externalFolder);
            externalFile = fullfile(externalFolder, 'ext.dat');
            copyfile(fullfile(testCase.TempFolder, 'green.dat'), externalFile);

            outFile = correctMotionArtifact(externalFile, testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'motionCorrected.dat'));
            testCase.verifyTrue(isfile(outFile));
        end

        function testSingleFrameDatIsRejected(testCase)
            iWriteSingleDat(fullfile(testCase.TempFolder, 'one.dat'), ...
                zeros(testCase.Ny, testCase.Nx, 1, 'single'), 10, 20);

            testCase.verifyError(@() correctMotionArtifact('one.dat', testCase.TempFolder), ...
                'Umitoolbox:correctMotionArtifact:TooFewFrames');
        end

        function testRejectsUnsupportedInputs(testCase)
            dataYXT = single(loadData(fullfile(testCase.TempFolder, 'green.dat')));
            umt = genUMTStruct(dataYXT, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.TempFolder, 'green.umt');
            saveData(umtFile, umt);
            byEvent = fullfile(testCase.TempFolder, 'byEvent.dat');
            saveData(byEvent, cat(4, dataYXT, dataYXT), 'DimNames', {'Y','X','T','E'}, ...
                'FrameRateHz', 10);
            folder = testCase.TempFolder;

            testCase.verifyError(@() correctMotionArtifact(umt, folder), ...
                'Umitoolbox:correctMotionArtifact:UnsupportedInputType');
            testCase.verifyError(@() correctMotionArtifact(umtFile, folder), ...
                'Umitoolbox:correctMotionArtifact:UnsupportedInputFile');
            testCase.verifyError(@() correctMotionArtifact(byEvent, folder), ...
                'Umitoolbox:correctMotionArtifact:UnsupportedLayout');
            testCase.verifyError(@() correctMotionArtifact(cat(4, dataYXT, dataYXT), folder), ...
                'Umitoolbox:correctMotionArtifact:UnsupportedLayout');
            testCase.verifyError(@() correctMotionArtifact(dataYXT(:, :, 1), folder), ...
                'Umitoolbox:correctMotionArtifact:UnsupportedLayout');
        end

        function testPipelineInfoDeclaresContract(testCase)
            info = correctMotionArtifact('pipelineInfo');

            testCase.verifyEqual({info.inputs.name}, {'data', 'SaveFolder'});
            testCase.verifyEqual({info.parameters.name}, ...
                {'TransformType', 'UpsamplingFactor', 'SaveShifts', 'ShowPlot'});
            testCase.verifyEqual({info.outputs.name}, {'outData'});
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'motionCorrected.dat');
            testCase.verifyEqual(info.outputs(1).type, {'ImageTimeSeries', 'ProcessedData'});
        end
    end

    methods (Access = private)
        function [outData, params] = runAndReadShifts(testCase, dataFile, varargin)
            % Run with SaveShifts on and read the shift matrix back from the
            % provenance file (the only place it is available).
            outData = correctMotionArtifact(dataFile, testCase.TempFolder, ...
                varargin{:}, 'SaveShifts', true, 'ShowPlot', false);
            [~, stem] = fileparts(dataFile);
            S = load(fullfile(testCase.TempFolder, [stem '_MotionCorrection.mat']), ...
                'MotionCorrection');
            params = S.MotionCorrection.params;
        end
    end
end

function stack = iBuildShiftedStack(pattern, expectedShifts)
%IBUILDSHIFTEDSTACK Build a YXT stack matching known correctMotionArtifact shifts.
%
% correctMotionArtifact corrects frame t by translating it forward by
% shifts(t,:) = [rowShift, colShift] so it lands back on frame 1. For the
% estimate to come out equal to expectedShifts, frame t must therefore be
% frame 1 translated forward by the negated vector.

nt = size(expectedShifts, 1);
stack = zeros([size(pattern), nt], 'single');
for t = 1:nt
    rowShift = expectedShifts(t, 1);
    colShift = expectedShifts(t, 2);
    stack(:,:,t) = imtranslate(pattern, [-colShift, -rowShift], 'cubic', 'FillValues', 0);
end
end

function iWriteSingleDat(filePath, data, frameRateHz, exposureMsec)
%IWRITESINGLEDAT Write YXT single data to a .dat file.

% Headered input (.dat header Phase 5a).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
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
