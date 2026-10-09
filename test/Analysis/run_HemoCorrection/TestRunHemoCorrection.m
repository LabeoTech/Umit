classdef TestRunHemoCorrection < matlab.unittest.TestCase
    %TESTRUNHEMOCORRECTION Unit tests for run_HemoCorrection.
    %
    %   These tests use a synthetic AcqInfos-driven dataset with unequal
    %   imported-channel timelines. The yellow fluorescence channel follows
    %   the base timeline, while the red reference channel has twice the
    %   number of frames and twice the frame rate.

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
                ['TestRunHemoCorrection_' char(java.util.UUID.randomUUID)]);
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
        function testLinearRegressionNumericResamplesReference(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            out = run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'LinearRegression', ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifySize(out, size(dataYXT));
        end

        function testLinearRegressionRawDatResamplesReference(testCase)
            out = run_HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'Algorithm', 'LinearRegression', ...
                'Red', true, 'Green', false, 'Amber', false);

            testCase.verifyEqual(out, 'hemoCorr_fluo.dat');
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'hemoCorr_fluo.dat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'fluoHemoCorr.dat')));
        end

        function testRatiometricNumericResamplesReference(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            out = run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifySize(out, size(dataYXT));
        end

        function testNumericInputIgnoresAcqInfosFrameSize(testCase)
            % .dat header Phase 4e-2: after an alignment AcqInfos.mat keeps
            % the raw frame size while the (headered) files have the new
            % one. In-memory input takes Y/X from the array itself.
            dataYXT = iLoadNumericData(testCase.TempFolder);
            iHeaderizeFolder(testCase.TempFolder);
            iSetAcqInfosFrameSize(testCase.TempFolder, size(dataYXT, 1) + 10, size(dataYXT, 2) + 6);

            for algorithm = {'LinearRegression', 'Ratiometric'}
                out = run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                    'Algorithm', algorithm{1}, 'Red', true, 'Green', false, 'Amber', false, ...
                    'FrameRateHz', testCase.AcqInfo.FrameRateHz);
                testCase.verifySize(out, size(dataYXT), algorithm{1});
            end
        end

        function testNumericInputSpatialMismatchStillRefused(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);
            iHeaderizeFolder(testCase.TempFolder);
            cropped = dataYXT(1:end-1, :, :);

            testCase.verifyError(@() run_HemoCorrection(cropped, testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', 'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:run_HemoCorrection:SpatialMismatch');
            testCase.verifyError(@() run_HemoCorrection(cropped, testCase.TempFolder, ...
                'Algorithm', 'LinearRegression', 'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:SpatialMismatch');
        end

        function testRatiometricRawDatResamplesReference(testCase)
            out = run_HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false);

            testCase.verifyEqual(out, 'hemoCorr_fluo.dat');
            outPath = fullfile(testCase.TempFolder, out);
            testCase.verifyTrue(isfile(outPath));

            info = loadMetaData(outPath);
            testCase.verifyEqual(datAxisSize(info, 'T'), testCase.AcqInfo.Length);
            testCase.verifyEqual(info.frameRateHz, testCase.AcqInfo.FrameRateHz);
        end

        function testInvalidAlgorithm(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'BadAlgo'), ...
                'Umitoolbox:run_HemoCorrection:InvalidInput');
        end

        function testRatiometricRequiresSingleReference(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', true, 'Amber', false), ...
                'Umitoolbox:run_HemoCorrection:InvalidInput');
        end

        function testMissingReferenceFile(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', false, 'Green', true, 'Amber', false), ...
                'Umitoolbox:run_HemoCorrection:FileNotFound');
        end

        function testInRamInputWithoutFrameRateErrors(testCase)
            % Numeric fluorescence input needs 'FrameRateHz'; AcqInfos.mat
            % in the folder is not used for it.
            testCase.assertTrue(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));
            dataYXT = iLoadNumericData(testCase.TempFolder);

            testCase.verifyError(@() run_HemoCorrection(dataYXT, testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false), ...
                'Umitoolbox:run_HemoCorrection:missingFrameRateHz');
        end

        function testRerunOverwritesDeclaredOutput(testCase)
            % A pre-existing output must be replaced, not sidestepped under a
            % different name: the declared pipeline output is hemoCorr_fluo.dat,
            % so a re-run has to keep writing there and must not leave a stale
            % original behind (P1-5).
            staleFile = fullfile(testCase.TempFolder, 'hemoCorr_fluo.dat');
            fid = fopen(staleFile, 'w');
            fclose(fid);
            staleBytes = dir(staleFile).bytes;

            out = run_HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false);

            testCase.verifyEqual(out, 'hemoCorr_fluo.dat');
            testCase.verifyTrue(isfile(staleFile));
            testCase.verifyGreaterThan(dir(staleFile).bytes, staleBytes);
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'hemoCorr_fluo_preallocData.dat')));
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'hemoCorr_fluo_writing.dat')));
        end

        function testLinearRegressionRerunOverwritesDeclaredOutput(testCase)
            % The wrapper owns hemoCorr_fluo.dat even though the legacy
            % LinearRegression core initially writes fluoHemoCorr.dat.
            staleFile = fullfile(testCase.TempFolder, 'hemoCorr_fluo.dat');
            fid = fopen(staleFile, 'w');
            fclose(fid);
            staleBytes = dir(staleFile).bytes;

            out = run_HemoCorrection('yellow.dat', testCase.TempFolder, ...
                'Algorithm', 'LinearRegression', ...
                'Red', true, 'Green', false, 'Amber', false);

            testCase.verifyEqual(out, 'hemoCorr_fluo.dat');
            testCase.verifyTrue(isfile(staleFile));
            testCase.verifyGreaterThan(dir(staleFile).bytes, staleBytes);
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'fluoHemoCorr_preallocData.dat')));
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'fluoHemoCorr_writing.dat')));
            testCase.verifyFalse( ...
                isfile(fullfile(testCase.TempFolder, 'fluoHemoCorr.dat')));
        end

        function testLinearRegressionRejectsReferenceDurationMismatch(testCase)
            mismatchFolder = fullfile(testCase.TempFolder, 'durationMismatchLinear');
            mkdir(mismatchFolder);
            mismatchAcqInfo = iCreateDurationMismatchDataset(mismatchFolder);
            dataYXT = iLoadNumericData(mismatchFolder);

            testCase.verifyError(@() run_HemoCorrection(dataYXT, mismatchFolder, ...
                'Algorithm', 'LinearRegression', ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', mismatchAcqInfo.FrameRateHz), ...
                'Umitoolbox:HemoCorrection:DurationMismatch');
        end

        function testRatiometricRejectsReferenceDurationMismatch(testCase)
            mismatchFolder = fullfile(testCase.TempFolder, 'durationMismatchRatiometric');
            mkdir(mismatchFolder);
            mismatchAcqInfo = iCreateDurationMismatchDataset(mismatchFolder);
            dataYXT = iLoadNumericData(mismatchFolder);

            testCase.verifyError(@() run_HemoCorrection(dataYXT, mismatchFolder, ...
                'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', mismatchAcqInfo.FrameRateHz), ...
                'Umitoolbox:run_HemoCorrection:DurationMismatch');
        end

        function testPipelineInfo(testCase)
            info = run_HemoCorrection('pipelineInfo');
            testCase.verifyEqual(info.name, 'run_HemoCorrection');
            testCase.verifyFalse(info.legacyOpts);
            testCase.verifyEqual(info.outputs(1).type, {'ImageTimeSeries', 'ProcessedData'});
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'hemoCorr_fluo.dat');

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(all(ismember({'ImageTimeSeries', 'ProcessedData'}, ...
                dataInput.type)));

            % The injected frame rate is the fluorescence one: it is taken
            % from the "data" input, never from a reference channel.
            rateDecl = info.sourceInfo(strcmp({info.sourceInfo.name}, 'FrameRateHz'));
            testCase.verifyEqual(char(string(rateDecl.sourceInput)), 'data');
            testCase.verifyEqual(char(string(rateDecl.sourceField)), 'frameRateHz');
        end

        function testPipelineManagerInjectsTheFluorescenceFrameRate(testCase)
            % yellow.dat (fluorescence) is 10 Hz while red.dat (reference) is
            % 20 Hz. An injected reference rate would differ from the .dat
            % header of the fluorescence and raise the rate-conflict warning.
            workFolder = iPreparePMHemoCorrectionFolder(testCase.TempFolder);
            parameters = struct('Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false);
            pm = buildPMForScenario(workFolder, 'run_HemoCorrection', 'auto', ...
                'Input', 'yellow.dat', 'Parameters', parameters);

            testCase.verifyWarningFree(@() pm.executePipeline('PrintSummary', false));

            outInfo = loadMetaData(fullfile(workFolder, 'hemoCorr_fluo.dat'));
            testCase.verifyEqual(outInfo.frameRateHz, testCase.AcqInfo.FrameRateHz);
        end

        function testRunsWithoutAcqInfos(testCase)
            % Nothing reads AcqInfos.mat: rates come from the .dat headers (or
            % FrameRateHz), so a folder without it must work for both algorithms.
            dataYXT = iLoadNumericData(testCase.TempFolder);
            delete(fullfile(testCase.TempFolder, 'AcqInfos.mat'));

            for algorithm = {'LinearRegression', 'Ratiometric'}
                args = {'Algorithm', algorithm{1}, ...
                    'Red', true, 'Green', false, 'Amber', false};
                outArray = run_HemoCorrection(dataYXT, testCase.TempFolder, args{:}, ...
                    'FrameRateHz', testCase.AcqInfo.FrameRateHz);
                testCase.verifySize(outArray, size(dataYXT), algorithm{1});

                outFile = run_HemoCorrection('yellow.dat', testCase.TempFolder, args{:});
                testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, outFile)), algorithm{1});
            end
        end

        function testRejectsUMTAndUnsupportedLayouts(testCase)
            dataYXT = iLoadNumericData(testCase.TempFolder);
            umt = genUMTStruct(dataYXT, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.TempFolder, 'fluo.umt');
            saveData(umtFile, umt);
            yxeFile = fullfile(testCase.TempFolder, 'fluo_yxe.dat');
            saveData(yxeFile, dataYXT(:, :, 1:3), 'DimNames', {'Y','X','E'}, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz);
            args = {'Red', true, 'Green', false, 'Amber', false, ...
                'FrameRateHz', testCase.AcqInfo.FrameRateHz};

            testCase.verifyError(@() run_HemoCorrection(umt, testCase.TempFolder, args{:}), ...
                'Umitoolbox:run_HemoCorrection:UnsupportedInputType');
            testCase.verifyError(@() run_HemoCorrection(umtFile, testCase.TempFolder, args{:}), ...
                'Umitoolbox:run_HemoCorrection:UnsupportedInputFile');
            testCase.verifyError(@() run_HemoCorrection(yxeFile, testCase.TempFolder, args{:}), ...
                'Umitoolbox:run_HemoCorrection:unsupportedLayout');
            testCase.verifyError(@() run_HemoCorrection(dataYXT(:, :, 1), testCase.TempFolder, args{:}), ...
                'Umitoolbox:run_HemoCorrection:InvalidInput');
        end

        % -----------------------------------------------------------------
        % Event-split (Y-X-T-E) fluorescence
        % -----------------------------------------------------------------
        function testEventSplitRatiometricEqualRatesMatchesPerTrialReference(testCase)
            s = iBuildEventFolder(testCase);

            out = run_HemoCorrection('fluoByEv.dat', s.folder, 'Algorithm', 'Ratiometric', ...
                'Red', false, 'Green', false, 'Amber', false, 'Other', 'ref.dat');

            testCase.verifyEqual(out, 'hemoCorr_fluo.dat');
            outInfo = loadMetaData(fullfile(s.folder, out));
            testCase.verifyEqual(outInfo.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(double(outInfo.dimSizes(:)).', size(s.fluoTrials));
            testCase.verifyFalse(isfile(fullfile(s.folder, 'hemoCorr_fluo_writing.dat')));

            % The reference trials come from split_data_by_event itself: at
            % equal rates the function must reproduce them exactly.
            expected = iExpectedRatiometric(s.fluoTrials, s.refTrials);
            testCase.verifyLessThan(max(abs(single(loadData(fullfile(s.folder, out))) - expected), [], 'all'), 1e-5);
        end

        function testEventSplitLinearRegressionEqualRatesMatchesPerTrialRef(testCase)
            s = iBuildEventFolder(testCase);

            out = run_HemoCorrection('fluoByEv.dat', s.folder, 'Algorithm', 'LinearRegression', ...
                'Red', false, 'Green', false, 'Amber', false, 'Other', 'ref.dat');

            expected = iExpectedLinear(s.fluoTrials, s.refTrials);
            actual = single(loadData(fullfile(s.folder, out)));
            testCase.verifyEqual(size(actual), size(expected));
            testCase.verifyLessThan(max(abs(actual - expected), [], 'all'), 1e-4);
        end

        function testEventSplitUpsamplesASlowerReference(testCase)
            % refHalf.dat holds every second frame at half the rate: every
            % trial frame must be interpolated between its reference frames.
            s = iBuildEventFolder(testCase);

            out = run_HemoCorrection('fluoByEv.dat', s.folder, 'Algorithm', 'Ratiometric', ...
                'Red', false, 'Green', false, 'Amber', false, 'Other', 'refHalf.dat');

            refAtTrials = iInterpolateReference(s.refHalf, s.rate / 2, s.trialTimes);
            expected = iExpectedRatiometric(s.fluoTrials, refAtTrials);
            testCase.verifyLessThan(max(abs(single(loadData(fullfile(s.folder, out))) - expected), [], 'all'), 1e-5);
        end

        function testEventSplitDownsamplesAFasterReference(testCase)
            % refDouble.dat is the same smooth reference at twice the rate:
            % after the anti-alias filter and the interpolation the corrected
            % trials must agree with the equal-rate result.
            s = iBuildEventFolder(testCase);
            args = {'Algorithm', 'Ratiometric', 'Red', false, 'Green', false, 'Amber', false};

            run_HemoCorrection('fluoByEv.dat', s.folder, args{:}, 'Other', 'ref.dat');
            equalRate = single(loadData(fullfile(s.folder, 'hemoCorr_fluo.dat')));
            run_HemoCorrection('fluoByEv.dat', s.folder, args{:}, 'Other', 'refDouble.dat');
            fromFaster = single(loadData(fullfile(s.folder, 'hemoCorr_fluo.dat')));

            testCase.verifyLessThan(max(abs(fromFaster - equalRate), [], 'all'), 2e-3);
        end

        function testEventSplitArrayMatchesDatWithTwoReferencesOfDifferentRates(testCase)
            s = iBuildEventFolder(testCase);
            copyfile(fullfile(s.folder, 'ref.dat'), fullfile(s.folder, 'red.dat'));
            copyfile(fullfile(s.folder, 'refHalf.dat'), fullfile(s.folder, 'green.dat'));
            args = {'Algorithm', 'LinearRegression', 'Red', true, 'Green', true, 'Amber', false};

            outDat = run_HemoCorrection('fluoByEv.dat', s.folder, args{:});
            outArray = run_HemoCorrection(s.fluoTrials, s.folder, args{:}, ...
                'FrameRateHz', s.rate);

            testCase.verifyTrue(isnumeric(outArray));
            testCase.verifySize(outArray, size(s.fluoTrials));
            testCase.verifyLessThan(max(abs(outArray - single(loadData(fullfile(s.folder, outDat)))), [], 'all'), 1e-5);
        end

        function testEventSplitRejectsUnmatchedEventsAndReferenceDurations(testCase)
            s = iBuildEventFolder(testCase);
            args = {'Algorithm', 'Ratiometric', 'Red', false, 'Green', false, ...
                'Amber', false, 'Other', 'ref.dat'};

            % One E slice fewer than the event instances of events.mat.
            testCase.verifyError(@() run_HemoCorrection( ...
                s.fluoTrials(:, :, :, 1:end-1), s.folder, args{:}, 'FrameRateHz', s.rate), ...
                'Umitoolbox:run_HemoCorrection:EventsNotMatched');

            % A reference that does not span the same recording.
            saveData(fullfile(s.folder, 'red.dat'), s.ref(:, :, 1:end-4), ...
                'DimNames', {'Y','X','T'}, 'FrameRateHz', s.rate);
            copyfile(fullfile(s.folder, 'ref.dat'), fullfile(s.folder, 'green.dat'));
            testCase.verifyError(@() run_HemoCorrection('fluoByEv.dat', s.folder, ...
                'Algorithm', 'LinearRegression', 'Red', true, 'Green', true, 'Amber', false), ...
                'Umitoolbox:run_HemoCorrection:DurationMismatch');
        end

        function testPipelineManagerExecutesForcedChunkAllRamScenarios(testCase)
            mocksFolder = fullfile(testCase.ProjectRoot, 'test', 'subFunc', ...
                'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '1900'));

            prepareFcn = @() iPreparePMHemoCorrectionFolder(testCase.TempFolder);
            parameters = struct('Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false);
            outputs = pmCollectScenarioOutputs(prepareFcn, 'run_HemoCorrection', ...
                'Input', 'yellow.dat', ...
                'PMOptions', {'Parameters', parameters});

            for iScenario = 1:numel(outputs)
                testCase.verifyEqual(outputs(iScenario).files, ...
                    {'hemoCorr_fluo.dat'}, ...
                    sprintf('run_HemoCorrection output mismatch under %s.', ...
                    outputs(iScenario).scenario));
            end
        end
    end
end

function folder = iPreparePMHemoCorrectionFolder(rootFolder)
folder = tempname(rootFolder);
mkdir(folder);
iCreateRepeatedTimelineDataset(folder);
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

function iHeaderizeFolder(folder)
% Give every .dat a v1 header (same values; rate and exposure from
% loadMetaData of the headerless file).
listing = dir(fullfile(folder, '*.dat'));
for k = 1:numel(listing)
    f = fullfile(folder, listing(k).name);
    if isDatWithHeader(f)
        continue
    end
    info = loadMetaData(f);
    values = loadData(f);
    [~, base] = fileparts(f);
    hdr = datHeaderFromInfo(info, base);
    hdr.writeComplete = true;
    fid = fopen(f, 'w', 'ieee-le');
    fwrite(fid, encodeDatHeader(hdr), 'uint8');
    fwrite(fid, values, info.dataClass);
    fclose(fid);
end
end

function iSetAcqInfosFrameSize(folder, height, width)
S = load(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
AcqInfoStream = S.AcqInfoStream;
AcqInfoStream.Height = height;
AcqInfoStream.Width = width;
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
end

% =========================================================================
% Event-split fixtures and per-trial references
% =========================================================================
function s = iBuildEventFolder(testCase)
%IBUILDEVENTFOLDER Synthetic recordings on the clock of a real events.mat.
%
% A 6 x 6 fluorescence recording and smooth reference channels at the same,
% half, and double frame rate share the duration of the events fixture, so
% the real event timing applies. The fluorescence is split by
% split_data_by_event and saved as fluoByEv.dat (Y-X-T-E).

fixture = fullfile(testCase.ProjectRoot, 'test', 'Analysis', 'TestingData_with_events');
s.folder = fullfile(testCase.TempFolder, 'events');
mkdir(s.folder);
copyfile(fullfile(fixture, 'AcqInfos.mat'), s.folder);
copyfile(fullfile(fixture, 'events.mat'), s.folder);

acq = load(fullfile(s.folder, 'AcqInfos.mat'), 'AcqInfoStream');
s.rate = double(acq.AcqInfoStream.FrameRateHz);
nT = datAxisSize(loadMetaData(fullfile(fixture, 'green.dat')), 'T');
nT = 2 * floor(nT / 2);

Ny = 6;
Nx = 6;
t = (0:nT-1) / s.rate;
stream = RandStream('mt19937ar', 'Seed', 3);
slow = reshape(sin(2 * pi * t / (t(end) / 3)), 1, 1, nT);
other = reshape(cos(2 * pi * t / (t(end) / 7)), 1, 1, nT);
pixRef = 1 + 0.05 * rand(stream, Ny, Nx);
pixFluo = 1 + 0.05 * rand(stream, Ny, Nx);

s.ref = single(pixRef .* (1 + 0.10 * slow));
fluo = single(pixFluo .* (1 + 0.05 * slow + 0.03 * other) + 0.002 * randn(stream, Ny, Nx, nT));

s.refHalf = s.ref(:, :, 1:2:end);
refDouble = zeros(Ny, Nx, 2 * nT, 'single');
refDouble(:, :, 1:2:end) = s.ref;
refDouble(:, :, 2:2:end) = (s.ref + s.ref(:, :, [2:end end])) / 2;

dims = {'Y','X','T'};
saveData(fullfile(s.folder, 'ref.dat'), s.ref, 'DimNames', dims, 'FrameRateHz', s.rate);
saveData(fullfile(s.folder, 'refHalf.dat'), s.refHalf, 'DimNames', dims, 'FrameRateHz', s.rate / 2);
saveData(fullfile(s.folder, 'refDouble.dat'), refDouble, 'DimNames', dims, 'FrameRateHz', 2 * s.rate);

evSplit = EventsManager(s.folder);
s.fluoTrials = single(evSplit.splitDataByEvents(fluo, 'FrameRateHz', s.rate, 'IncludeIgnored', true));
s.refTrials = single(evSplit.splitDataByEvents(s.ref, 'FrameRateHz', s.rate, 'IncludeIgnored', true));
saveData(fullfile(s.folder, 'fluoByEv.dat'), s.fluoTrials, ...
    'DimNames', {'Y','X','T','E'}, 'FrameRateHz', s.rate);

% Frame times (s) of every trial frame, from the same events.mat.
frMat = EventsManager(s.folder).getFrameMatrix(nT + 2, '', [], ...
    'FrameRateHz', s.rate, 'IncludeIgnored', true);
s.trialTimes = (frMat(:, 1:size(s.fluoTrials, 3)) - 1) / s.rate;
end

function trials = iInterpolateReference(ref, refRate, trialTimes)
%IINTERPOLATEREFERENCE Linear interpolation of a continuous reference at trial times.
[Ny, Nx, M] = size(ref);
[Ne, Nt] = size(trialTimes);
pos = min(max(trialTimes * refRate + 1, 1), M);
series = double(reshape(ref, [], M).');                 % M x pixels
trials = zeros(Ny, Nx, Nt, Ne, 'single');
for e = 1:Ne
    vals = interp1((1:M).', series, pos(e, :).', 'linear');   % Nt x pixels
    trials(:, :, :, e) = reshape(single(vals).', Ny, Nx, Nt);
end
end

function out = iExpectedRatiometric(fluoTrials, refTrials)
%IEXPECTEDRATIOMETRIC Ratiometric correction of every trial on its own.
mF = mean(fluoTrials, 3);
mR = mean(refTrials, 3);
out = (((fluoTrials - mF) ./ mF) - ((refTrials - mR) ./ mR)) .* mF + mF;
end

function out = iExpectedLinear(fluoTrials, refTrials)
%IEXPECTEDLINEAR Per-trial, per-pixel regression on constant, drift and reference.
[Ny, Nx, Nt, Ne] = size(fluoTrials);
out = zeros(size(fluoTrials), 'single');
base = [ones(Nt, 1), linspace(0, 1, Nt).'];
for e = 1:Ne
    refSmooth = imgaussfilt(refTrials(:, :, :, e), 1, 'Padding', 'symmetric');
    mR = mean(refSmooth, 3);
    refNorm = (refSmooth - mR) ./ mR;
    f = fluoTrials(:, :, :, e);
    mF = mean(f, 3);
    fNorm = (f - mF) ./ mF;
    for y = 1:Ny
        for x = 1:Nx
            X = [base, double(reshape(refNorm(y, x, :), Nt, 1))];
            yy = double(reshape(fNorm(y, x, :), Nt, 1));
            residual = yy - X * (X \ yy);
            out(y, x, :, e) = single(reshape( ...
                residual * double(mF(y, x)) + double(mF(y, x)), 1, 1, Nt));
        end
    end
end
end
