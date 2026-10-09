classdef TestSpatialGaussFilt < matlab.unittest.TestCase
    %TESTSPATIALGAUSSFILT Unit tests for spatialGaussFilt.
    %
    % Coverage:
    %   1) pipelineInfo generation
    %   2) Raw YXT array input
    %   3) Raw YXTE array input (every frame filtered on its own)
    %   4) Raw .dat filename input (YXT and YXTE)
    %   5) Standard vs low-RAM comparison, including forced multi-chunk runs
    %      and the command-window progress of the low-RAM mode
    %   6) NaN handling and integer-class input
    %   7) Rejection of UMT structs, .umt files, and unsupported layouts
    %
    % Ground truth:
    %   Uses direct calls to IMGAUSSFILT with the same NaN handling as the
    %   function to compare algorithms.

    properties
        SaveFolder char
        ProjectRoot char
        SampleDataFolder char
    end

    methods (TestClassSetup)
        function addProjectPathAndLocateFixture(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);

            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end

            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            cfg = [];
            if isappdata(0, 'SpatialGaussFiltTestConfig')
                cfg = getappdata(0, 'SpatialGaussFiltTestConfig');
            end

            if ~isempty(cfg) && isfield(cfg, 'sampleDataFolder') && ...
                    isfolder(cfg.sampleDataFolder)
                testCase.SampleDataFolder = char(string(cfg.sampleDataFolder));
            else
                testCase.SampleDataFolder = fullfile( ...
                    testCase.ProjectRoot, ...
                    'test', ...
                    'Analysis', ...
                    'TestingData_with_events');
            end

            testCase.assertTrue( ...
                isfile(fullfile(testCase.SampleDataFolder, 'green.dat')), ...
                ['Missing fixture file "green.dat". Put it in: ' ...
                testCase.SampleDataFolder]);

            testCase.assertTrue( ...
                isfile(fullfile(testCase.SampleDataFolder, 'AcqInfos.mat')), ...
                ['Missing fixture file "AcqInfos.mat". Put it in: ' ...
                testCase.SampleDataFolder]);
        end
    end

    methods (TestMethodSetup)
        function createFreshWorkspace(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.SaveFolder = fx.Folder;

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'green.dat'), ...
                fullfile(testCase.SaveFolder, 'green.dat'));

            copyfile( ...
                fullfile(testCase.SampleDataFolder, 'AcqInfos.mat'), ...
                fullfile(testCase.SaveFolder, 'AcqInfos.mat'));

            deleteIfExists(fullfile(testCase.SaveFolder, 'spatialGaussFilt.dat'));
        end
    end

    methods (TestMethodTeardown)
        function cleanupCreatedFiles(testCase)
            deleteIfExists(fullfile(testCase.SaveFolder, 'spatialGaussFilt.dat'));
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = spatialGaussFilt('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));

            reqFields = {'name','description','version','inputs','outputs'};
            testCase.verifyTrue(all(ismember(reqFields, fieldnames(info))), ...
                'pipelineInfo is missing one or more required top-level fields.');

            testCase.verifyEqual(info.name, 'spatialGaussFilt');
            testCase.verifyEqual(numel(info.inputs), 2);
            testCase.verifyEqual(numel(info.parameters), 1);
            testCase.verifyEqual(numel(info.outputs), 1);

            inputNames = {info.inputs.name};
            outputNames = {info.outputs.name};
            paramNames = {info.parameters.name};
            testCase.verifyEqual(inputNames, {'data','SaveFolder'});
            testCase.verifyEqual(paramNames, {'Sigma'});
            testCase.verifyEqual(outputNames, {'outData'});
            testCase.verifyEqual(info.outputs(1).defOutfilename, 'spatialGaussFilt.dat');
        end

        function testArrayInput(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            rawData(1,1,:) = NaN;
            rawData(5,3,:) = NaN;

            sigma = 1.25;
            expected = testCase.expectedFrames(rawData, sigma);

            out = spatialGaussFilt(rawData, testCase.SaveFolder, 'Sigma', sigma);

            testCase.verifyTrue(isnumeric(out));
            testCase.verifyEqual(size(out), size(rawData));
            testCase.verifyNumericEquivalent(single(out), single(expected));
            testCase.verifyTrue(isnan(out(1,1,1)));
            testCase.verifyTrue(isnan(out(5,3,7)));
        end

        function testDatInput(testCase)
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));
            sigma = 1.1;

            expected = testCase.expectedFrames(rawData, sigma);

            out = spatialGaussFilt(inFile, testCase.SaveFolder, 'Sigma', sigma);

            testCase.verifyTrue(ischar(out) || (isstring(out) && isscalar(out)));
            outFile = char(string(out));
            testCase.verifyTrue(isfile(outFile));

            actual = single(loadData(outFile));
            testCase.verifyEqual(size(actual), size(rawData));
            testCase.verifyNumericEquivalent(actual, single(expected));
        end

        function testStandardVsLowRAM(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            sigma = 1.4;

            outStandard = spatialGaussFilt(rawData, testCase.SaveFolder, 'Sigma', sigma);
            outLowRAMFile = spatialGaussFilt(fullfile(testCase.SaveFolder, 'green.dat'), ...
                testCase.SaveFolder, 'Sigma', sigma);

            testCase.verifyTrue(ischar(outLowRAMFile) || (isstring(outLowRAMFile) && isscalar(outLowRAMFile)));
            testCase.verifyTrue(isfile(char(string(outLowRAMFile))));

            outLowRAM = single(loadData(char(string(outLowRAMFile))));

            testCase.verifyEqual(size(outStandard), size(outLowRAM));
            testCase.verifyNumericEquivalent(outLowRAM, outStandard);
        end

        function testLowRAMPreservesNaNLikeStandardMode(testCase)
            % The in-RAM path zero-fills NaN, filters, then restores NaN. The
            % .dat path used to call IMGAUSSFILT on the raw slab, which lets
            % NaN bleed outwards by the kernel radius, so the same input gave
            % materially different results in the two modes (P1-4).
            inFile = fullfile(testCase.SaveFolder, 'green.dat');
            rawData = single(loadData(inFile));

            % Mask a spatial border for the whole time series, which is the
            % normal state after GSR or normalization in this toolbox.
            rawData(1:3, :, :) = NaN;
            rawData(:, end-2:end, :) = NaN;
            inInfo = loadMetaData(inFile);
            writeSingleDatFile(inFile, rawData, inInfo.frameRateHz, inInfo.exposureMsec);

            sigma = 1.4;
            outStandard = spatialGaussFilt(rawData, testCase.SaveFolder, 'Sigma', sigma);
            outLowRAMFile = spatialGaussFilt(inFile, testCase.SaveFolder, 'Sigma', sigma);
            outLowRAM = single(loadData(char(string(outLowRAMFile))));

            testCase.verifyEqual(size(outLowRAM), size(outStandard));

            % The NaN mask itself must not grow in file mode.
            testCase.verifyEqual(isnan(outLowRAM), isnan(outStandard), ...
                'Low-RAM mode changed the NaN mask.');

            finiteIdx = isfinite(outStandard) & isfinite(outLowRAM);
            testCase.verifyTrue(any(finiteIdx(:)));
            testCase.verifyLessThanOrEqual( ...
                max(abs(double(outStandard(finiteIdx)) - double(outLowRAM(finiteIdx)))), 1e-4);
        end

        function testLowRAMWritesDeclaredOutputName(testCase)
            % The declared pipeline output must be the file actually written,
            % not a scaffolding "_PREALLOC" name (P1-5).
            info = spatialGaussFilt('pipelineInfo');
            declaredName = info.outputs(1).defOutfilename;
            if iscell(declaredName)
                declaredName = declaredName{1};
            end

            out = spatialGaussFilt(fullfile(testCase.SaveFolder, 'green.dat'), ...
                testCase.SaveFolder, 'Sigma', 1.1);

            [~, actualStem, actualExt] = fileparts(char(string(out)));
            testCase.verifyEqual([actualStem actualExt], char(string(declaredName)));
            testCase.verifyTrue(isfile(char(string(out))));
            testCase.verifyEmpty(dir(fullfile(testCase.SaveFolder, '*_writing.dat')));
        end

        function testYXTEArrayAndDatFilterEachFrame(testCase)
            % Event-split data: every Y-X frame is filtered on its own and the
            % output keeps the input's dimensions and axes, for both modes.
            sigma = 1.1;
            trials = testCase.buildSyntheticYXTE();
            expected = testCase.expectedFrames(trials, sigma);

            outArray = spatialGaussFilt(trials, testCase.SaveFolder, 'Sigma', sigma);
            testCase.verifyEqual(size(outArray), size(trials));
            testCase.verifyNumericEquivalent(single(outArray), single(expected));

            inFile = testCase.writeYXTEDat(trials);
            outFile = char(string(spatialGaussFilt(inFile, testCase.SaveFolder, 'Sigma', sigma)));
            hdr = readDatHeader(outFile);
            testCase.verifyEqual(hdr.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(hdr.dimSizes, size(trials));
            testCase.verifyNumericEquivalent(single(loadData(outFile)), single(expected));
        end

        function testForcedMultiChunkMatchesInRamAndPrintsProgress(testCase)
            % The memory mock forces many frame blocks. The result must equal
            % the in-RAM one for YXT and YXTE (also with NaN pixels, whose
            % mask is per frame), and the low-RAM mode reports its progress in
            % the command window.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            sigma = 1.3;
            stream = RandStream('mt19937ar', 'Seed', 4);
            movie = randn(stream, 12, 10, 20, 'single');
            movie(2, 3, 5:9) = NaN;          % a pixel that is NaN in some frames only
            movie(8, 6, :) = NaN;            % and one that is NaN throughout
            trials = testCase.buildSyntheticYXTE();
            trials(4, 5, :, 2) = NaN;

            cases = {movie, trials};
            for iCase = 1:numel(cases)
                data = cases{iCase};
                axesNames = {'Y','X','T','E'};
                inFile = fullfile(testCase.SaveFolder, 'blocked.dat');
                saveData(inFile, data, 'DimNames', axesNames(1:ndims(data)), 'FrameRateHz', 10);

                progress = evalc('outFile = spatialGaussFilt(inFile, testCase.SaveFolder, ''Sigma'', sigma);');

                nChunks = str2double(regexp(progress, ...
                    'frame\(s\) in (\d+) chunk', 'tokens', 'once'));
                testCase.verifyGreaterThan(nChunks, 1, ...
                    'The fixture must force more than one frame block.');
                testCase.verifyNotEmpty(regexp(progress, 'Chunk 1/\d+ \[Reading', 'once'));
                testCase.verifyNotEmpty(regexp(progress, ...
                    sprintf('Chunk %d/%d \\[Writing', nChunks, nChunks), 'once'));

                inRam = spatialGaussFilt(data, testCase.SaveFolder, 'Sigma', sigma);
                actual = single(loadData(char(string(outFile))));
                testCase.verifyEqual(isnan(actual), isnan(inRam));
                testCase.verifyNumericEquivalent(actual, single(inRam));
            end
        end

        function testIntegerAndLogicalArraysAreFilteredAsSingle(testCase)
            % imgaussfilt returns its input class: an integer array would be
            % rounded back to integers. The function converts to single first.
            raw = uint16(mod(reshape(1:12*10*4, 12, 10, 4), 251));

            out = spatialGaussFilt(raw, testCase.SaveFolder, 'Sigma', 1.2);

            testCase.verifyClass(out, 'single');
            testCase.verifyNumericEquivalent(out, ...
                single(testCase.expectedFrames(single(raw), 1.2)));
            testCase.verifyClass(spatialGaussFilt(raw > 100, testCase.SaveFolder), 'single');
        end

        function testRejectsUMTAndUnsupportedLayouts(testCase)
            rawData = single(loadData(fullfile(testCase.SaveFolder, 'green.dat')));
            umt = genUMTStruct(rawData, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});
            umtFile = fullfile(testCase.SaveFolder, 'input.umt');
            saveData(umtFile, umt);
            yxFile = fullfile(testCase.SaveFolder, 'yx.dat');
            saveData(yxFile, rawData(:, :, 1), 'DimNames', {'Y','X'}, 'FrameRateHz', 10);
            yxeFile = fullfile(testCase.SaveFolder, 'yxe.dat');
            saveData(yxeFile, rawData(:, :, 1:3), 'DimNames', {'Y','X','E'}, 'FrameRateHz', 10);

            testCase.verifyError(@() spatialGaussFilt(umt, testCase.SaveFolder), ...
                'spatialGaussFilt:UnsupportedInputType');
            testCase.verifyError(@() spatialGaussFilt(umtFile, testCase.SaveFolder), ...
                'spatialGaussFilt:UnsupportedInputFile');
            testCase.verifyError(@() spatialGaussFilt(yxFile, testCase.SaveFolder), ...
                'Umitoolbox:spatialGaussFilt:unsupportedLayout');
            testCase.verifyError(@() spatialGaussFilt(yxeFile, testCase.SaveFolder), ...
                'Umitoolbox:spatialGaussFilt:unsupportedLayout');
            testCase.verifyError(@() spatialGaussFilt(rawData(:, :, 1), testCase.SaveFolder), ...
                'spatialGaussFilt:InvalidArrayInput');
        end
    end

    methods (Access = private)

        function expected = expectedFrames(~, block, sigma)
            % Reference: every frame filtered with its own NaN mask (zero-fill,
            % filter, restore).
            frames = reshape(block, size(block, 1), size(block, 2), []);
            spatialMask = isnan(frames);
            work = frames;
            work(spatialMask) = 0;
            filtered = imgaussfilt(work, sigma, 'FilterDomain', 'spatial');
            filtered(spatialMask) = NaN;
            expected = reshape(filtered, size(block));
        end

        function trials = buildSyntheticYXTE(~)
            stream = RandStream('mt19937ar', 'Seed', 9);
            trials = randn(stream, 12, 10, 3, 2, 'single');
        end

        function inFile = writeYXTEDat(testCase, trials)
            inFile = fullfile(testCase.SaveFolder, 'byEvent.dat');
            saveData(inFile, trials, 'DimNames', {'Y','X','T','E'}, 'FrameRateHz', 10);
        end

        function verifyNumericEquivalent(testCase, a, b)
            testCase.verifyEqual(size(a), size(b));

            if isequaln(a, b)
                return
            end

            diffVals = double(a(:)) - double(b(:));
            diffVals = diffVals(isfinite(diffVals));

            if isempty(diffVals)
                testCase.verifyEqual(a, b);
                return
            end

            testCase.verifyLessThanOrEqual(std(diffVals, 0, 'omitnan'), 1e-4);
        end
    end
end

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function writeSingleDatFile(filePath, data, frameRateHz, exposureMsec)
%WRITESINGLEDATFILE Overwrite a single-precision .dat file in place.

% Headered input (.dat header Phase 5a).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
end
