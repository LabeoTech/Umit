classdef TestRunBloodFlow < matlab.unittest.TestCase
    %TESTRUNBLOODFLOW Unit tests for run_BloodFlow.
    %
    % Uses the real fixture folder:
    %   Analysis/TestingData_speckle

    properties
        FixtureFolder = ''
        TempFolder = ''
    end

    properties (TestParameter)
        sType = {'Spatial', 'Temporal'}
        kernelSize = {3, 5, 7}
    end

    methods (TestMethodSetup)
        function setup(testCase)
            thisFile = which(class(testCase));
            testRoot = fileparts(thisFile);
            analysisRoot = fileparts(testRoot);
            testCase.FixtureFolder = fullfile(analysisRoot, 'TestingData_speckle');

            testCase.assumeTrue(isfolder(testCase.FixtureFolder), ...
                sprintf('Fixture folder not found: %s', testCase.FixtureFolder));

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestRunBloodFlow_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            copyfile(fullfile(testCase.FixtureFolder, 'speckle.dat'), ...
                fullfile(testCase.TempFolder, 'speckle.dat'));
            copyfile(fullfile(testCase.FixtureFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.TempFolder, 'AcqInfos.mat'));
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = run_BloodFlow('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyTrue(isfield(info, 'inputs'));
            testCase.verifyTrue(isfield(info, 'outputs'));

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');
        end

        function testPipelineInfoParameterTypesAreDialogSupported(testCase)
            % PipelineManager.buildParameterDialog only renders these types;
            % any other type (e.g. 'double') errors and hides the control.
            info = run_BloodFlow('pipelineInfo');
            names = {info.parameters.name};
            testCase.verifyEqual(sort(names), ...
                sort({'sType', 'ExposureMsec', 'KernelSize', 'bNormalize'}));
            for k = 1:numel(info.parameters)
                testCase.verifyTrue(ismember(lower(info.parameters(k).type), ...
                    {'logical', 'numeric', 'char'}), ...
                    sprintf('Parameter "%s" has type "%s".', ...
                    info.parameters(k).name, info.parameters(k).type));
            end
            kernel = info.parameters(strcmp(names, 'KernelSize'));
            testCase.verifyEqual(kernel.type, 'numeric');
            testCase.verifyEqual(kernel.default, 5);
        end

        function testExposureAcceptsNumericText(testCase)
            % The dialog edits the mixed 'auto'/number ExposureMsec in a text
            % field, so a typed number arrives as char.
            outText = iFlow(testCase, 'ExposureMsec', '2');
            outNum = iFlow(testCase, 'ExposureMsec', 2);
            testCase.verifyEqual(outText, outNum);

            testCase.verifyError( ...
                @() run_BloodFlow(testCase.TempFolder, 'speckle.dat', 'ExposureMsec', 'abc'), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end

        function testMatchesFormula(testCase, sType, kernelSize)
            datFile = fullfile(testCase.TempFolder, 'speckle.dat');
            dat = loadData(datFile);
            md = loadMetaData(datFile);
            T = double(md.exposureMsec) / 1000;   % the file's own (speckle) exposure

            out = iFlow(testCase, 'sType', sType, 'KernelSize', kernelSize);

            testCase.verifyClass(out, 'single');
            testCase.verifySize(out, size(dat));

            % Exact reference: explicit std over each symmetric-padded
            % window at sampled voxels, including first/last frames and
            % image borders.
            [ny, nx, nt] = size(dat);
            I = double(dat) ./ mean(double(dat), 3, 'omitnan');
            h = (kernelSize - 1) / 2;
            disk = fspecial('disk', h) > 0;
            % Interior, border, and first/last-frame voxels (coordinates scaled
            % by 1/4 when the speckle fixture was binned to 128 x 128 in .dat
            % header Phase 5a).
            samples = [1 1 1; ny nx nt; 50 50 300; 96 5 1; ...
                88 101 nt; 2 nx-1 2; ny 3 nt-1];
            for k = 1:size(samples, 1)
                y = samples(k, 1); x = samples(k, 2); t = samples(k, 3);
                if strcmpi(sType, 'Spatial')
                    frame = padarray(I(:, :, t), [h h], 'symmetric');
                    win = frame(y:y+2*h, x:x+2*h);
                    kExact = std(win(disk));
                else
                    trace = padarray(squeeze(I(y, x, :)), h, 'symmetric');
                    kExact = std(trace(t:t+2*h));
                end
                expected = single(1 / (T * kExact^2));
                testCase.verifyEqual(out(y, x, t), expected, 'RelTol', single(1e-5), ...
                    sprintf('Voxel (%d,%d,%d)', y, x, t));
            end
        end

        function testWritesAHeaderedYXTFile(testCase, sType)
            datFile = fullfile(testCase.TempFolder, 'speckle.dat');

            outFile = run_BloodFlow(testCase.TempFolder, 'speckle.dat', 'sType', sType, ...
                'KernelSize', 7);

            testCase.verifyEqual(outFile, fullfile(testCase.TempFolder, 'BloodFlow.dat'));
            testCase.verifyTrue(isfile(outFile));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'BloodFlow_compute.dat')));
            outInfo = loadMetaData(outFile);
            inInfo = loadMetaData(datFile);
            testCase.verifyEqual(outInfo.dimNames, {'Y','X','T'});
            testCase.verifyEqual(outInfo.dimSizes, inInfo.dimSizes);
            testCase.verifyEqual(outInfo.dataClass, 'single');
            testCase.verifyEqual(outInfo.frameRateHz, inInfo.frameRateHz);
        end

        function testForcedMultiChunkMatchesSingleChunk(testCase, sType)
            % The memory mock forces many chunks (temporal slabs for Spatial,
            % X slabs for Temporal); the result must equal the one-chunk run,
            % with a non-default window to exercise chunk borders.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            expected = iFlow(testCase, 'sType', sType, 'KernelSize', 7);
            expectedNorm = iFlow(testCase, 'sType', sType, 'KernelSize', 7, ...
                'bNormalize', true);

            projectRoot = extractBefore(mfilename('fullpath'), [filesep 'test' filesep]);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                projectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000000'));
            info = loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat'));
            totalBytes = prod(double(info.dimSizes)) * 4;
            testCase.assertGreaterThan(calculateMaxChunkSize(totalBytes, 24, .15), 1, ...
                'The fixture must force more than one chunk.');

            actual = iFlow(testCase, 'sType', sType, 'KernelSize', 7);
            actualNorm = iFlow(testCase, 'sType', sType, 'KernelSize', 7, 'bNormalize', true);

            testCase.verifyEqual(actual, expected, 'RelTol', single(1e-4));
            testCase.verifyEqual(actualNorm, expectedNorm, 'RelTol', single(1e-4));
        end

        function testNormalizeDividesByTemporalMean(testCase, sType)
            dat = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));

            raw = iFlow(testCase, 'sType', sType);
            out = iFlow(testCase, 'sType', sType, 'bNormalize', true);

            expected = single(double(raw) ./ mean(double(raw), 3, 'omitnan'));
            testCase.verifyClass(out, 'single');
            testCase.verifySize(out, size(dat));
            testCase.verifyEqual(out, expected, 'RelTol', single(1e-6));

            pixelMean = mean(double(out), 3, 'omitnan');
            pixelMean = pixelMean(isfinite(pixelMean));
            testCase.verifyEqual(pixelMean, ones(size(pixelMean)), 'AbsTol', 1e-5);
        end

        function testExposureOverrideScalesFlow(testCase)
            out1 = iFlow(testCase, 'ExposureMsec', 1);
            out2 = iFlow(testCase, 'ExposureMsec', 2);

            testCase.verifyEqual(out2, out1 / 2, 'RelTol', single(1e-5));
        end

        function testDefaultSpatialKernelMatchesSpeckleMapping(testCase)
            % The default spatial window must be SpeckleMapping's 21-pixel
            % 5 x 5 disk (corners excluded), not a full square.
            dat = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));
            out = iFlow(testCase, 'sType', 'Spatial');
            ref = iFlow(testCase, 'sType', 'Spatial', 'KernelSize', 5);
            testCase.verifyEqual(out, ref);

            I = double(dat(:, :, 1)) ./ mean(double(dat), 3, 'omitnan') - 1;
            K = stdfilt(I, fspecial('disk', 2) > 0);
            T = double(loadMetaData(fullfile(testCase.TempFolder, 'speckle.dat')) ...
                .exposureMsec) / 1000;   % the file's own (speckle) exposure
            testCase.verifyEqual(out(:, :, 1), single(1 ./ (T .* K.^2)), ...
                'RelTol', single(1e-5));
        end

        function testInvalidKernelSizeErrors(testCase)
            for badSize = {1, 4, 0, -3, 2.5, [3 3], '5'}
                testCase.verifyError( ...
                    @() run_BloodFlow(testCase.TempFolder, 'speckle.dat', ...
                    'KernelSize', badSize{1}), ...
                    'MATLAB:InputParser:ArgumentFailedValidation');
            end
        end

        function testMissingFileInputErrors(testCase)
            delete(fullfile(testCase.TempFolder, 'speckle.dat'));

            testCase.verifyError( ...
                @() run_BloodFlow(testCase.TempFolder, 'speckle.dat'), ...
                'Umitoolbox:run_BloodFlow:FileNotFound');
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            % Arrays are not an input form.
            testCase.verifyError(@() run_BloodFlow(sv, rand(6, 5, 20, 'single')), ...
                'Umitoolbox:run_BloodFlow:UnsupportedInputType');

            % Other extensions.
            testCase.verifyError(@() run_BloodFlow(sv, 'speckle.tif'), ...
                'Umitoolbox:run_BloodFlow:UnsupportedInputFile');

            % Layouts other than Y-X-T.
            layouts = {{'Y','X','T','E'}, [6 5 4 2]; {'Y','X','E'}, [6 5 3]; {'Y','X'}, [6 5]};
            for k = 1:size(layouts, 1)
                f = fullfile(sv, 'layout.dat');
                writeTestDat(f, rand([layouts{k, 2}, 1], 'single'), 10, 5, ...
                    'DimNames', layouts{k, 1});
                before = dir(sv);
                testCase.verifyError(@() run_BloodFlow(sv, f), ...
                    'Umitoolbox:run_BloodFlow:unsupportedLayout');
                testCase.verifyEqual(numel(dir(sv)), numel(before), ...
                    'no output may be written for a refused input');
            end
        end
    end
end

function out = iFlow(testCase, varargin)
%IFLOW Run run_BloodFlow on the fixture file and return its output values.
outFile = run_BloodFlow(testCase.TempFolder, 'speckle.dat', varargin{:});
out = loadData(outFile);
end
