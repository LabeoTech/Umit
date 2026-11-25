classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestCorrectMotionArtifact < matlab.unittest.TestCase
    %TESTCORRECTMOTIONARTIFACT Unit tests for correctMotionArtifact.
    %
    % The function returns the corrected data (array in -> array out, .dat in
    % -> <name>_MotionCorrected.dat). The estimated shifts are only available in the
    % <name>_MotionCorrection.mat that SaveShifts=true writes, so the shift
    % tests read them from there.

    properties
        TempFolder
        Ny = 48  % non-square on purpose (DFR-20261109-002)
        Nx = 80
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
            for channel = {'green.dat', 'red.dat'}
                AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                    'DatFile', channel{1}, 'Length', testCase.Nt, ...
                    'FrameRateHz', 10, 'ExposureMsec', 20));
            end
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

        function testDatOutputIsSuffixedAndSourcesAreUntouched(testCase)
            greenBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat'));
            redBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat'));

            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'green_MotionCorrected.dat'));
            testCase.verifyTrue(isfile(outFile));
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat')), greenBefore);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat')), redBefore);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));

            info = loadMetaData(outFile);
            testCase.verifyEqual(info.dimNames, {'Y','X','T'});
            testCase.verifyEqual(double(info.dimSizes(:)).', [testCase.Ny testCase.Nx testCase.Nt]);
            testCase.verifyEqual(info.dataClass, 'single');
        end

        function testRerunReplacesItsOwnPreviousOutput(testCase)
            staleFile = fullfile(testCase.TempFolder, 'green_MotionCorrected.dat');
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
            rowsIn = 10:(testCase.Ny-10);
            colsIn = 10:(testCase.Nx-10);
            % AbsTol is ignored for single arrays, so compare in double.
            for t = 1:testCase.Nt
                testCase.verifyEqual( ...
                    double(greenCorrected(rowsIn, colsIn, t)), ...
                    double(testCase.RefPattern(rowsIn, colsIn)), 'AbsTol', 0.1);
                testCase.verifyEqual( ...
                    double(redCorrected(rowsIn, colsIn, t)), ...
                    double(testCase.RefPattern(rowsIn, colsIn)) * 0.6, 'AbsTol', 0.1);
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
            testCase.verifyEqual(hdr.channelName, 'red_MotionCorrected');
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

            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'ext_MotionCorrected.dat'));
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

        function testOverwriteReplacesTheSourceInPlace(testCase)
            sourceFile = fullfile(testCase.TempFolder, 'green.dat');
            before = computeFileChecksum(sourceFile);

            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'Overwrite', true, 'UpsamplingFactor', 10, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(outFile, sourceFile);
            testCase.verifyNotEqual(computeFileChecksum(sourceFile), before);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_MotionCorrected.dat')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));
            testCase.verifyTrue(readDatHeader(sourceFile).writeComplete);
            rowsIn = 10:(testCase.Ny-10);
            colsIn = 10:(testCase.Nx-10);
            corrected = double(loadData(sourceFile));
            for t = 1:testCase.Nt
                testCase.verifyEqual(corrected(rowsIn, colsIn, t), ...
                    double(testCase.RefPattern(rowsIn, colsIn)), 'AbsTol', 0.1);
            end
        end

        function testRunLogsEveryWrittenFileEvenWithoutSaveShifts(testCase)
            correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            S = load(fullfile(testCase.TempFolder, 'MotionCorrectionLog.mat'), 'MotionCorrectionLog');
            testCase.verifyEqual({S.MotionCorrectionLog.file}, {'green_MotionCorrected.dat'});
            testCase.verifyEqual(S.MotionCorrectionLog.inputFile, 'green.dat');
            testCase.verifyEqual(S.MotionCorrectionLog.bytes, ...
                dir(fullfile(testCase.TempFolder, 'green_MotionCorrected.dat')).bytes);
        end

        function testLogStoresTheModifiedTimeAsReadableText(testCase)
            correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            S = load(fullfile(testCase.TempFolder, 'MotionCorrectionLog.mat'), 'MotionCorrectionLog');
            entry = S.MotionCorrectionLog;
            testCase.verifyClass(entry.modifiedText, 'char');
            testCase.verifyNotEmpty(regexp(entry.modifiedText, ...
                '^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}$', 'once'));
            % The text is the same instant as the numeric value the guard compares.
            readBack = datetime(entry.modifiedText, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', ...
                'TimeZone', 'local');
            testCase.verifyEqual(posixtime(readBack) * 1000, entry.modified, 'AbsTol', 1);
        end

        function testOlderLogWithoutReadableTimeIsUpgradedAndStillRecognised(testCase)
            % green.dat is corrected in place: only the log can recognise it afterwards.
            correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'Overwrite', true, 'SaveShifts', false, 'ShowPlot', false);
            logFile = fullfile(testCase.TempFolder, 'MotionCorrectionLog.mat');
            S = load(logFile, 'MotionCorrectionLog');
            MotionCorrectionLog = rmfield(S.MotionCorrectionLog, 'modifiedText');
            save(logFile, 'MotionCorrectionLog');

            correctMotionArtifact('red.dat', testCase.TempFolder, ...
                'Overwrite', true, 'SaveShifts', false, 'ShowPlot', false);

            S = load(logFile, 'MotionCorrectionLog');
            log = S.MotionCorrectionLog;
            testCase.verifyEqual(sort({log.file}), {'green.dat', 'red.dat'});
            testCase.verifyTrue(all(arrayfun(@(e) ~isempty(e.modifiedText), log)));
            % The upgraded entry still carries its numeric time, so it stays recognised.
            testCase.verifyError(@() correctMotionArtifact('green.dat', ...
                testCase.TempFolder, 'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
        end

        function testSuffixedCorrectedFileIsRefused(testCase)
            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);
            before = computeFileChecksum(outFile);

            testCase.verifyError(@() correctMotionArtifact(outFile, testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
            testCase.verifyEqual(computeFileChecksum(outFile), before);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_MotionCorrected_MotionCorrected.dat')));
        end

        function testSuffixIsRefusedEvenWithoutALog(testCase)
            % A corrected copy brought in from elsewhere has no log entry.
            copyfile(fullfile(testCase.TempFolder, 'green.dat'), ...
                fullfile(testCase.TempFolder, 'other_MotionCorrected.dat'));

            testCase.verifyError(@() correctMotionArtifact('other_MotionCorrected.dat', ...
                testCase.TempFolder, 'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
        end

        function testFileOverwrittenInPlaceIsRefused(testCase)
            % No suffix to tell it apart: only the log can.
            sourceFile = fullfile(testCase.TempFolder, 'green.dat');
            correctMotionArtifact('green.dat', testCase.TempFolder, 'Overwrite', true, ...
                'SaveShifts', false, 'ShowPlot', false);
            before = computeFileChecksum(sourceFile);

            testCase.verifyError(@() correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'Overwrite', true, 'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
            testCase.verifyError(@() correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
            testCase.verifyEqual(computeFileChecksum(sourceFile), before);
        end

        function testAllowRecorrectionWarnsAndProceeds(testCase)
            correctMotionArtifact('green.dat', testCase.TempFolder, 'Overwrite', true, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyWarning(@() correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'AllowRecorrection', true, 'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'green_MotionCorrected.dat')));
        end

        function testReplacedFileIsNoLongerTreatedAsCorrected(testCase)
            % Importing again produces a new file under the same name.
            correctMotionArtifact('green.dat', testCase.TempFolder, 'Overwrite', true, ...
                'SaveShifts', false, 'ShowPlot', false);
            pause(0.2);
            iWriteSingleDat(fullfile(testCase.TempFolder, 'green.dat'), ...
                iBuildShiftedStack(testCase.RefPattern, testCase.Shifts), 10, 20);

            outFile = correctMotionArtifact('green.dat', testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyTrue(isfile(outFile));
        end

        function testBatchCorrectsEveryListedChannelIndependently(testCase)
            greenBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat'));
            redBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat'));

            outFiles = correctMotionArtifact(testCase.TempFolder, ...
                'UpsamplingFactor', 10, 'ShowPlot', false);

            testCase.verifyEqual(outFiles, { ...
                fullfile(testCase.TempFolder, 'green_MotionCorrected.dat'), ...
                fullfile(testCase.TempFolder, 'red_MotionCorrected.dat')});
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat')), greenBefore);
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat')), redBefore);
            for k = 1:numel(outFiles)
                testCase.verifyTrue(isfile(outFiles{k}));
            end
            greenParams = iReadParams(testCase.TempFolder, 'green');
            redParams = iReadParams(testCase.TempFolder, 'red');
            testCase.verifyEqual(greenParams(:, [2 1]), testCase.Shifts, 'AbsTol', 0.15);
            testCase.verifyEqual(redParams(:, [2 1]), testCase.RedShifts, 'AbsTol', 0.15);
            S = load(fullfile(testCase.TempFolder, 'MotionCorrectionLog.mat'), 'MotionCorrectionLog');
            testCase.verifyEqual(sort({S.MotionCorrectionLog.file}), ...
                {'green_MotionCorrected.dat', 'red_MotionCorrected.dat'});
        end

        function testBatchOverwriteReplacesEveryChannelInPlace(testCase)
            greenBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat'));
            redBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'red.dat'));

            outFiles = correctMotionArtifact(testCase.TempFolder, 'Overwrite', true, ...
                'UpsamplingFactor', 10, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(outFiles, { ...
                fullfile(testCase.TempFolder, 'green.dat'), ...
                fullfile(testCase.TempFolder, 'red.dat')});
            testCase.verifyNotEqual(computeFileChecksum(outFiles{1}), greenBefore);
            testCase.verifyNotEqual(computeFileChecksum(outFiles{2}), redBefore);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_MotionCorrected.dat')));

            % A second batch finds nothing left to correct and changes nothing.
            after = computeFileChecksum(outFiles{1});
            testCase.verifyError(@() correctMotionArtifact(testCase.TempFolder, ...
                'Overwrite', true, 'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');
            testCase.verifyEqual(computeFileChecksum(outFiles{1}), after);
        end

        function testBatchSkipsAnAlreadyCorrectedChannelWithAWarning(testCase)
            correctMotionArtifact('green.dat', testCase.TempFolder, 'Overwrite', true, ...
                'SaveShifts', false, 'ShowPlot', false);
            greenAfterFirst = computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat'));

            testCase.verifyWarning(@() correctMotionArtifact(testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:AlreadyCorrected');

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'red_MotionCorrected.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'green_MotionCorrected.dat')));
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat')), greenAfterFirst);
        end

        function testBatchSkipsAListedChannelWhoseFileIsMissing(testCase)
            iAddImportedChannel(testCase.TempFolder, 'yellow.dat', testCase.Nt);

            testCase.verifyWarning(@() correctMotionArtifact(testCase.TempFolder, ...
                'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:MissingChannelFile');

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'green_MotionCorrected.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'red_MotionCorrected.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'yellow_MotionCorrected.dat')));
        end

        function testBatchChecksEveryChannelBeforeWritingAnything(testCase)
            iWriteSingleDat(fullfile(testCase.TempFolder, 'yellow.dat'), ...
                zeros(testCase.Ny, testCase.Nx, 1, 'single'), 10, 20);
            iAddImportedChannel(testCase.TempFolder, 'yellow.dat', 1);

            testCase.verifyError(@() correctMotionArtifact(testCase.TempFolder, ...
                'Overwrite', true, 'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:TooFewFrames');

            % green and red come first in the list but must be untouched.
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, 'MotionCorrectionLog.mat')));
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_writing.dat')));
        end

        function testBatchErrorsWhenNoChannelFileExists(testCase)
            delete(fullfile(testCase.TempFolder, 'green.dat'));
            delete(fullfile(testCase.TempFolder, 'red.dat'));

            testCase.verifyError(@() correctMotionArtifact(testCase.TempFolder, ...
                'ShowPlot', false), 'Umitoolbox:correctMotionArtifact:NoFiles');
        end

        function testBatchOnALegacyFolderIncludesSpeckleAndFluo(testCase)
            % No ImportedChannels: the channels come from the declared
            % illuminations. The shared fallback only knows red, green and
            % yellow, so speckle and fluo must not be dropped.
            folder = iMakeLegacyFolder(testCase, 'legacyAll', { ...
                {'green.dat', 'Green'}, {'red.dat', 'Red'}, {'speckle.dat', 'Speckle'}, ...
                {'fluo_475.dat', 'Fluo #1 475nm'}, {'fluo.dat', 'Fluo'}});

            outFiles = correctMotionArtifact(folder, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(sort(cellfun(@(f) iFileName(f), outFiles, 'UniformOutput', false)), ...
                sort({'green_MotionCorrected.dat', 'red_MotionCorrected.dat', ...
                'speckle_MotionCorrected.dat', 'fluo_475_MotionCorrected.dat', ...
                'fluo_MotionCorrected.dat'}));
            S = load(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyFalse(isfield(S.AcqInfoStream, 'ImportedChannels'));
        end

        function testBatchOnALegacyFolderWithOnlySpeckleAndFluo(testCase)
            % Here the shared fallback finds nothing at all.
            folder = iMakeLegacyFolder(testCase, 'legacyNoRgy', { ...
                {'speckle.dat', 'Speckle'}, {'fluo_567.dat', 'Fluo #1 567nm'}});

            outFiles = correctMotionArtifact(folder, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(sort(cellfun(@(f) iFileName(f), outFiles, 'UniformOutput', false)), ...
                sort({'speckle_MotionCorrected.dat', 'fluo_567_MotionCorrected.dat'}));
        end

        function testLegacyDeclaredChannelWithAMissingFileIsSkippedNotAnError(testCase)
            folder = iMakeLegacyFolder(testCase, 'legacyMissing', { ...
                {'green.dat', 'Green'}, {'red.dat', 'Red'}, {'speckle.dat', 'Speckle'}});
            delete(fullfile(folder, 'red.dat'));
            delete(fullfile(folder, 'red.mat'));

            testCase.verifyWarning(@() correctMotionArtifact(folder, ...
                'SaveShifts', false, 'ShowPlot', false), ...
                'Umitoolbox:correctMotionArtifact:MissingChannelFile');

            testCase.verifyTrue(isfile(fullfile(folder, 'green_MotionCorrected.dat')));
            testCase.verifyTrue(isfile(fullfile(folder, 'speckle_MotionCorrected.dat')));
            testCase.verifyFalse(isfile(fullfile(folder, 'red_MotionCorrected.dat')));
        end

        function testLegacyFolderWithoutIlluminationsUsesTheFilesInTheFolder(testCase)
            folder = iMakeLegacyFolder(testCase, 'legacyBare', { ...
                {'speckle.dat', ''}, {'fluo.dat', ''}});

            outFiles = correctMotionArtifact(folder, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(sort(cellfun(@(f) iFileName(f), outFiles, 'UniformOutput', false)), ...
                sort({'speckle_MotionCorrected.dat', 'fluo_MotionCorrected.dat'}));
        end

        function testImportedChannelsStaysAuthoritative(testCase)
            % A speckle.dat that the registry does not list (for example a
            % derived file) is not an imported channel and is left alone.
            writeTestDat(fullfile(testCase.TempFolder, 'speckle.dat'), ...
                zeros(testCase.Ny, testCase.Nx, testCase.Nt, 'single') + 1, 10, 20);

            outFiles = correctMotionArtifact(testCase.TempFolder, 'SaveShifts', false, 'ShowPlot', false);

            testCase.verifyEqual(sort(cellfun(@(f) iFileName(f), outFiles, 'UniformOutput', false)), ...
                sort({'green_MotionCorrected.dat', 'red_MotionCorrected.dat'}));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'speckle_MotionCorrected.dat')));
        end

        function testBatchNeedsAcqInfos(testCase)
            delete(fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            % One argument: the folder alone is the batch form.
            testCase.verifyError(@() correctMotionArtifact(testCase.TempFolder), ...
                'Umitoolbox:correctMotionArtifact:InvalidInput');
        end

        function testPipelineInfoDeclaresContract(testCase)
            info = correctMotionArtifact('pipelineInfo');

            % A setup step: a folder input, no DATA input or output.
            testCase.verifyEqual({info.inputs.name}, {'SaveFolder'});
            testCase.verifyFalse(info.inputs(1).isData);
            testCase.verifyEqual({info.parameters.name}, {'TransformType', 'UpsamplingFactor', ...
                'Overwrite', 'AllowRecorrection', 'SaveShifts', 'ShowPlot'});
            testCase.verifyEqual({info.outputs.name}, {'correctedDatFiles'});
            testCase.verifyFalse(info.outputs(1).isData);
            testCase.verifyFalse(info.outputs(1).returnsValue);
            testCase.verifyEqual(info.outputs(1).defOutfilename, '*.dat');
        end

        function testPipelineRunsItAsASetupStepThatOverwritesTheChannels(testCase)
            % The point of the setup role: the channel files that
            % ImportedChannels lists are corrected in place before any later
            % step reads them.
            greenFile = fullfile(testCase.TempFolder, 'green.dat');
            redFile = fullfile(testCase.TempFolder, 'red.dat');
            greenBefore = computeFileChecksum(greenFile);
            redBefore = computeFileChecksum(redFile);

            pm = buildPMForScenario(testCase.TempFolder, 'correctMotionArtifact', 'auto', ...
                'Parameters', struct('Overwrite', true, 'UpsamplingFactor', 10, ...
                'SaveShifts', false, 'ShowPlot', false));
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyNotEqual(computeFileChecksum(greenFile), greenBefore);
            testCase.verifyNotEqual(computeFileChecksum(redFile), redBefore);
            testCase.verifyEmpty(dir(fullfile(testCase.TempFolder, '*_MotionCorrected.dat')));
            rowsIn = 10:(testCase.Ny-10);
            colsIn = 10:(testCase.Nx-10);
            corrected = double(loadData(greenFile));
            for t = 1:testCase.Nt
                testCase.verifyEqual(corrected(rowsIn, colsIn, t), ...
                    double(testCase.RefPattern(rowsIn, colsIn)), 'AbsTol', 0.1);
            end
        end

        function testPipelineDefaultKeepsTheChannelsUntouched(testCase)
            greenBefore = computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat'));

            pm = buildPMForScenario(testCase.TempFolder, 'correctMotionArtifact', 'auto', ...
                'Parameters', struct('SaveShifts', false, 'ShowPlot', false));
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyEqual(computeFileChecksum(fullfile(testCase.TempFolder, 'green.dat')), greenBefore);
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'green_MotionCorrected.dat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'red_MotionCorrected.dat')));
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

function params = iReadParams(folder, stem)
%IREADPARAMS Shift matrix saved by a run on <stem>.dat.
S = load(fullfile(folder, [stem '_MotionCorrection.mat']), 'MotionCorrection');
params = S.MotionCorrection.params;
end

function iAddImportedChannel(folder, datFile, nt)
%IADDIMPORTEDCHANNEL List one more .dat file in the folder's AcqInfos.mat.
S = load(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
AcqInfoStream = appendImportedChannelInfo(S.AcqInfoStream, struct( ...
    'DatFile', datFile, 'Length', nt, 'FrameRateHz', 10, 'ExposureMsec', 20));
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
end

function name = iFileName(filePath)
%IFILENAME Name with extension of a path.
[~, base, ext] = fileparts(filePath);
name = [base ext];
end

function folder = iMakeLegacyFolder(testCase, folderName, channels)
%IMAKELEGACYFOLDER Legacy folder: sidecar .dat files, Illumination<N> in
% AcqInfos.mat and no ImportedChannels. channels is {{datName, color}, ...};
% an empty color declares no illumination for that file.
folder = fullfile(testCase.TempFolder, folderName);
mkdir(folder);
stack = iBuildShiftedStack(testCase.RefPattern, testCase.Shifts);
AcqInfoStream = struct('Width', testCase.Nx, 'Height', testCase.Ny, ...
    'Length', testCase.Nt, 'FrameRateHz', 10, 'ExposureMsec', 20);
nIllum = 0;
for k = 1:numel(channels)
    writeTestDat(fullfile(folder, channels{k}{1}), stack, 10, 20, 'Format', 'legacySidecar');
    if ~isempty(channels{k}{2})
        nIllum = nIllum + 1;
        AcqInfoStream.(sprintf('Illumination%d', nIllum)) = ...
            struct('Color', channels{k}{2}, 'FrameRateHz', 10, 'ExposureMsec', 20);
    end
end
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
end
