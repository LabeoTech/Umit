classdef TestSpeckleMappingChunkOverrun < matlab.unittest.TestCase
    %TESTSPECKLEMAPPINGCHUNKOVERRUN Regression: RAM-safe chunk loops must not
    % overrun Nt/Nx when calculateMaxChunkSize's returned nChunks, combined
    % with a ceil()-rounded chunk size, needs fewer chunks than requested.
    %
    % Real system RAM cannot be controlled deterministically in a unit test,
    % so the overshoot is forced via the calculateMaxChunkSize RAM-query mock
    % (test/subFunc/calculateMaxChunkSize/mocks), following the same pattern
    % as TestHemoCompute's RAM-exhaustion/chunk-count-forcing tests.
    %
    % With TOTAL=10000, AVAILABLE=1900 against a 12x12x10 single-precision
    % fixture, all three RAM-safe chunking calls in SpeckleMapping.m
    % overshoot (nChunks vs. ceil(N/chunkSize)):
    %   Pass 1:            nChunks=7,   chunkT=2, actually needed=5
    %   'spatial' Pass 2:  nChunks=64,  chunkT=1, actually needed=10
    %   'temporal' Pass 2: nChunks=173, chunkX=1, actually needed=12
    % Before the fix, the last loop iteration(s) read/index past the data
    % (negative frame count on the Pass-1/spatial arithmetic paths, or an
    % empty xIdx reaching spatialSlabIO's x0 = xIdx(1) on the temporal
    % path), crashing with MATLAB:badsize_mx or MATLAB:badsubscript.

    properties
        TempFolder = ''
        ProjectRoot = ''
    end

    methods (TestClassSetup)
        function configurePath(testCase)
            thisFile = which(class(testCase));
            testRoot = fileparts(thisFile);
            testCase.ProjectRoot = fileparts(fileparts(fileparts(testRoot)));
        end
    end

    methods (TestMethodSetup)
        function setup(testCase)
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk-count forcing relies on shadowing the PCWIN64 memory() built-in.');

            testCase.TempFolder = fullfile(tempdir, ...
                ['TestSpeckleMappingChunkOverrun_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            Ny = 12; Nx = 12; Nt = 10;

            rng(7);
            data = 0.5 + rand(Ny, Nx, Nt, 'single');

            % Headered input (.dat header Phase 5a).
            writeTestDat(fullfile(testCase.TempFolder, 'speckle.dat'), data, 5);

            AcqInfoStream = struct( ...
                'Height', Ny, ...
                'Width', Nx, ...
                'Length', Nt, ...
                'FrameRateHz', 5, ...
                'Datatype', 'single');
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            mocksFolder = fullfile(testCase.ProjectRoot, ...
                'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocksFolder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(10000, '%d')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(1900, '%d')));
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
        function testSpatialModeDoesNotOverrunWithForcedExcessChunks(testCase)
            out = SpeckleMapping('speckle.dat', testCase.TempFolder, 'spatial', false, false);

            testCase.verifyTrue(isstruct(out) && isscalar(out));
            map = out.data.SpeckleMap.value;
            testCase.verifyFalse(any(isnan(map(:))), ...
                'A NaN-free synthetic fixture must not produce a NaN-populated map.');
        end

        function testTemporalModeDoesNotOverrunWithForcedExcessChunks(testCase)
            out = SpeckleMapping('speckle.dat', testCase.TempFolder, 'temporal', false, false);

            testCase.verifyTrue(isstruct(out) && isscalar(out));
            map = out.data.SpeckleMap.value;
            testCase.verifyFalse(any(isnan(map(:))), ...
                'A NaN-free synthetic fixture must not produce a NaN-populated map.');
        end

        function testSpatialModeMatchesStandardUnderForcedExcessChunks(testCase)
            % Not asserted bit-for-bit here: summing across many more
            % chunks than the 'temporal' path reassociates floating-point
            % addition differently than a single whole-array sum, which is
            % expected IEEE 754 noise (~1e-7 scale), not a regression. The
            % NaN-denominator fix's own bit-for-bit guarantee is covered by
            % TestSpeckleMappingNaNDenominator.m under a single chunk.
            outRAM = SpeckleMapping('speckle.dat', testCase.TempFolder, 'spatial', false, false);
            outStd = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'spatial', false, false);

            testCase.verifyEqual(outRAM.data.SpeckleMap.value, outStd.data.SpeckleMap.value, ...
                'AbsTol', single(1e-4));
        end

        function testTemporalModeMatchesStandardUnderForcedExcessChunks(testCase)
            outRAM = SpeckleMapping('speckle.dat', testCase.TempFolder, 'temporal', false, false);
            outStd = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'temporal', false, false);

            testCase.verifyEqual(outRAM.data.SpeckleMap.value, outStd.data.SpeckleMap.value);
        end

        function testChunkedArrayMatchesWholeArrayAlgorithm(testCase)
            % Under the forced RAM shortage an array no longer fits whole and
            % is chunked; the map must still equal the whole-array STDFILT
            % algorithm (computed here directly).
            arr = iArray(testCase);
            dat = arr ./ mean(arr, 3, 'omitnan');
            kernels = struct('spatial', single(fspecial('disk', 2) > 0), ...
                'temporal', ones(1, 1, 5, 'single'));

            for sType = {'spatial', 'temporal'}
                expected = single(mean(stdfilt(dat, kernels.(sType{1})), 3, 'omitnan'));
                out = SpeckleMapping(arr, testCase.TempFolder, sType{1}, false, false);

                testCase.verifyEqual(out.data.SpeckleMap.value, expected, ...
                    'AbsTol', single(1e-4), ...
                    sprintf('Chunked array map differs for sType=%s.', sType{1}));
            end
        end

        function testPipelineManagerForcedChunkAllRamScenarios(testCase)
            prepareFcn = @() iPreparePMSpeckleFolder( ...
                testCase.TempFolder, testCase.ProjectRoot);
            outputs = pmCollectScenarioOutputs(prepareFcn, ...
                'run_SpeckleMapping', 'Input', 'speckle.dat');

            for iScenario = 1:numel(outputs)
                testCase.verifyEqual(outputs(iScenario).files, ...
                    {'SpeckleMap.umt'}, ...
                    sprintf('run_SpeckleMapping output mismatch under %s.', ...
                    outputs(iScenario).scenario));
            end
        end
    end
end

function folder = iPreparePMSpeckleFolder(sourceFolder, projectRoot)
folder = tempname(sourceFolder);
mkdir(folder);
copyfile(fullfile(sourceFolder, 'speckle.dat'), folder);
copyfile(fullfile(sourceFolder, 'AcqInfos.mat'), folder);
ensurePMReadyAcqInfos(folder, 'speckle.dat');
addpath(genpath(fullfile(projectRoot, 'test', 'Analysis', 'PMTestHelpers')));
end

function arr = iArray(testCase)
%IARRAY The fixture recording as an in-RAM Y-X-T array (Standard mode).
arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));
end
