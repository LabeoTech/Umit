classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestNonSquareEquivariance < matlab.unittest.TestCase
    %TESTNONSQUAREEQUIVARIANCE Real recordings, cropped to non-square frames.
    %
    %   DFR-20261109-002, phase 2. The shared real fixtures
    %   (TestingData_with_events 56x56, TestingData_speckle 128x128,
    %   TestingData_retinotopy 64x64) are square and several tests compare them
    %   against stored hashes, so they cannot change. These tests crop them,
    %   inside the test, to non-square regions and check a property that holds
    %   for every pixelwise or isotropic function: the result for the
    %   Y/X-transposed recording equals the Y/X-transposed result for the
    %   original one. A transposed read, a swapped frame size or a wrong slab
    %   offset breaks it only when the frame is not square.
    %
    %   Every case runs twice: with the normal chunking and with a RAM mock
    %   that forces many chunks (the low-RAM loops over X or T).
    %
    %   The last test builds a small non-square raw acquisition and checks
    %   that ImagesClassification keeps the pixel orientation.

    properties
        TestRoot char
        ProjectRoot char
    end

    properties (TestParameter)
        chunking = {'normal', 'manyChunks'}
        filterCase = struct( ...
            'spatialGaussFilt', @(f, d) spatialGaussFilt(f, d, 'Sigma', 2), ...
            'normalizeZScore', @(f, d) normalizeZScore(f, d), ...
            'apply_detrend', @(f, d) apply_detrend(f, d), ...
            'GSR', @(f, d) GSR(f, d), ...
            'normalizeLPF', @(f, d) normalizeLPF(f, d, ...
                'BaselineCutoffHz', 0.05, 'SignalCutoffHz', 1), ...
            'normalizeBSLN', @(f, d) normalizeBSLN(f, d), ...
            'aggregateOverT', @(f, d) apply_aggregate_function(f, d, ...
                'aggregateFcn', 'mean', 'dimensionName', 'T'))
        speckleCase = struct( ...
            'SpeckleMappingTemporal', @(d, f) SpeckleMapping(d, f, 'sType', 'Temporal'), ...
            'SpeckleMappingSpatial', @(d, f) SpeckleMapping(d, f, 'sType', 'Spatial'), ...
            'AnaSpeckle', @(d, f) Ana_Speckle(d, f), ...
            'BloodFlowTemporal', @(d, f) run_BloodFlow(d, f, 'sType', 'Temporal'), ...
            'BloodFlowSpatial', @(d, f) run_BloodFlow(d, f, 'sType', 'Spatial'))
    end

    methods (TestClassSetup)
        function locateFixtures(testCase)
            testCase.TestRoot = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = fileparts(fileparts(testCase.TestRoot));
        end
    end

    methods (Test)
        function filtersOnCroppedEventsRecording(testCase, filterCase, chunking)
            % A 40 x 48 crop of the 56 x 56 recording.
            testCase.forceChunking(chunking);
            fcn = filterCase;
            results = cell(1, 2);
            for k = 1:2
                folder = testCase.prepareFolder('TestingData_with_events', 'green.dat', ...
                    5:44, 3:50, k == 2, {'AcqInfos.mat', 'events.mat'});
                results{k} = iLoadOutput(fcn('data.dat', folder), folder);
            end
            iVerifyTransposed(testCase, results{1}, results{2});
        end

        function eventPipelineOnCroppedRecording(testCase, chunking)
            % split -> amplitude maps / aggregation over E / trial baseline.
            testCase.forceChunking(chunking);
            stages = {'dataByEv', 'amplitude', 'aggregateOverE', 'trialBaseline'};
            results = struct();
            for k = 1:2
                folder = testCase.prepareFolder('TestingData_with_events', 'green.dat', ...
                    5:44, 3:50, k == 2, {'AcqInfos.mat', 'events.mat'});
                split = char(string(split_data_by_event(fullfile(folder, 'data.dat'), folder)));
                results(k).dataByEv = iLoadOutput(split, folder);
                results(k).amplitude = iLoadOutput(genAmplitudeMaps(split, folder), folder);
                results(k).aggregateOverE = iLoadOutput(apply_aggregate_function(split, folder, ...
                    'aggregateFcn', 'mean', 'dimensionName', 'E'), folder);
                results(k).trialBaseline = iLoadOutput(normalizeBSLN(split, folder, ...
                    'normalizationMode', 'trial'), folder);
            end
            for s = stages
                iVerifyTransposed(testCase, results(1).(s{1}), results(2).(s{1}), s{1});
            end
            testCase.verifyEqual(size(results(1).dataByEv, 1:2), [40 48]);
        end

        function speckleFamilyOnCroppedRecording(testCase, speckleCase, chunking)
            % A 40 x 56 crop of the 128 x 128 speckle recording.
            testCase.forceChunking(chunking);
            fcn = speckleCase;
            results = cell(1, 2);
            for k = 1:2
                folder = testCase.prepareFolder('TestingData_speckle', 'speckle.dat', ...
                    1:40, 1:56, k == 2, {'AcqInfos.mat'}, 1:200);
                results{k} = iLoadOutput(fcn(folder, 'data.dat'), folder);
            end
            iVerifyTransposed(testCase, results{1}, results{2});
        end

        function retinotopyOnCroppedRecording(testCase, chunking)
            % A 40 x 56 crop of the 64 x 64 retinotopy recording.
            testCase.forceChunking(chunking);
            maps = cell(1, 2);
            for k = 1:2
                folder = testCase.prepareFolder('TestingData_retinotopy', 'fluo.dat', ...
                    1:40, 1:56, k == 2, {'AcqInfos.mat', 'events.mat'});
                maps{k} = iCollectArrays(genRetinotopyMaps(fullfile(folder, 'data.dat'), folder));
            end
            testCase.assertNotEmpty(maps{1});
            testCase.assertEqual(numel(maps{1}), numel(maps{2}));
            for m = 1:numel(maps{1})
                iVerifyTransposed(testCase, maps{1}{m}, maps{2}{m}, sprintf('map %d', m));
            end
        end

        function roiMeansOnCroppedRecording(testCase)
            % ROI means must equal a direct extraction, for the recording and
            % for its transpose with transposed masks.
            wide = false(40, 48); wide(3:10, 25:44) = true;
            tall = false(40, 48); tall(20:36, 6:11) = true;
            for transposed = [false true]
                folder = testCase.prepareFolder('TestingData_with_events', 'green.dat', ...
                    5:44, 3:50, transposed, {'AcqInfos.mat', 'events.mat'}, 1:300);
                masks = {wide, tall};
                if transposed
                    masks = cellfun(@(m) m.', masks, 'UniformOutput', false);
                end
                roiFile = createROIFile(size(masks{1}), 'ROIs', ...
                    [iRectROI('Wide', masks{1}), iRectROI('Tall', masks{2})]);
                saveROIFile(fullfile(folder, 'myROI.roi'), roiFile);

                out = getDataFromROI('data.dat', folder, 'SpatialAggFcn', 'mean');
                values = double(out.data.(char(fieldnames(out.data))).value);

                data = double(loadData(fullfile(folder, 'data.dat')));
                flat = reshape(data, [], size(data, 3));
                expected = [mean(flat(masks{1}(:), :), 1); mean(flat(masks{2}(:), :), 1)];
                % Values are about 1e4 in single precision: compare relatively.
                testCase.verifyEqual(values, expected, 'RelTol', 1e-6, ...
                    sprintf('transposed = %d', transposed));
            end
        end

        function rawNonSquareAcquisitionKeepsPixelOrientation(testCase)
            % ImagesClassification on a synthetic 48 x 80 raw acquisition
            % (the shared raw fixture is 112 x 112). Raw layout: 5 int32
            % [-, nx, ny, -, -], then per frame 3 uint64 and nx*ny uint16
            % stored x-fastest.
            nx = 80; ny = 48; nColors = 3; cycles = 10; nFrames = nColors * cycles;
            fixtureInfo = fullfile(testCase.TestRoot, 'TestingData_with_events', 'info.txt');
            testCase.assumeTrue(isfile(fixtureInfo), 'Raw fixture info.txt is missing.');
            txt = fileread(fixtureInfo);
            txt = regexprep(txt, 'Width: \d+', sprintf('Width: %d', nx));
            txt = regexprep(txt, 'Height: \d+', sprintf('Height: %d', ny));

            rawFolder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            fid = fopen(fullfile(rawFolder, 'info.txt'), 'w'); fwrite(fid, txt); fclose(fid);
            [xx, yy] = meshgrid(1:nx, 1:ny);
            pixels = uint16((yy - 1) * nx + xx);     % unique per pixel, Y-first
            fid = fopen(fullfile(rawFolder, 'img_00000.bin'), 'w');
            fwrite(fid, int32([0 nx ny 0 0]), 'int32');
            for k = 1:nFrames
                fwrite(fid, uint64([k; 0; k]), 'uint64');
                fwrite(fid, pixels.' + uint16(4000 * mod(k, 12)), 'uint16');
            end
            fclose(fid);

            for binning = [1 2]
                saveFolder = testCase.applyFixture( ...
                    matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
                ImagesClassification(rawFolder, saveFolder, 'ApplyRig', false, ...
                    'backupOpts', 'ERASE', 'BinningSpatial', binning);
                names = {'red.dat', 'yellow.dat', 'green.dat'};
                for c = 1:3
                    frames = loadData(fullfile(saveFolder, names{c}));
                    testCase.verifySize(frames, [[ny nx] / binning, cycles]);
                end
                if binning == 1
                    red = loadData(fullfile(saveFolder, 'red.dat'));
                    for j = [1 2 5 cycles]
                        k = (j - 1) * nColors + 1;
                        testCase.verifyEqual(red(:, :, j), ...
                            single(pixels + uint16(4000 * mod(k, 12))), ...
                            sprintf('red frame %d', j));
                    end
                end
            end
        end
    end

    methods (Access = private)
        function folder = prepareFolder(testCase, fixtureName, datName, rows, cols, ...
                transposed, copyFiles, frames)
            %PREPAREFOLDER Temporary folder with a cropped (and transposed) data.dat.
            fixtureFolder = fullfile(testCase.TestRoot, fixtureName);
            testCase.assumeTrue(isfile(fullfile(fixtureFolder, datName)), ...
                sprintf('Fixture file %s/%s is missing.', fixtureName, datName));
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            for name = copyFiles
                copyfile(fullfile(fixtureFolder, name{1}), fullfile(folder, name{1}));
            end
            [data, info] = loadData(fullfile(fixtureFolder, datName));
            if nargin >= 8
                data = data(:, :, frames);
            end
            data = data(rows, cols, :);
            if transposed
                data = permute(data, [2 1 3]);
            end
            writeTestDat(fullfile(folder, 'data.dat'), single(data), info.frameRateHz, ...
                info.exposureMsec);
        end

        function forceChunking(testCase, chunking)
            %FORCECHUNKING Shadow memory() so the low-RAM loops use many chunks.
            if strcmp(chunking, 'normal')
                return
            end
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            mocks = fullfile(testCase.ProjectRoot, 'test', 'subFunc', ...
                'calculateMaxChunkSize', 'mocks');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mocks));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));
        end
    end
end

function A = iLoadOutput(out, folder)
%ILOADOUTPUT Array of a function output: a numeric array, a .dat path/name, or a UMT.
if isnumeric(out)
    A = out;
elseif ischar(out) || (isstring(out) && isscalar(out))
    out = char(out);
    if ~isfile(out)
        out = fullfile(folder, out);
    end
    A = loadData(out);
else
    arrays = iCollectArrays(out);
    A = arrays{1};
end
A = single(A);
while ~ismatrix(A) && size(A, ndims(A)) == 1
    A = reshape(A, size(A, 1:ndims(A) - 1));
end
end

function arrays = iCollectArrays(s)
%ICOLLECTARRAYS Numeric arrays (more than 4 values) inside a UMT struct, in order.
arrays = {};
if isnumeric(s)
    if numel(s) > 4
        arrays = {s};
    end
elseif iscell(s)
    for i = 1:numel(s)
        arrays = [arrays, iCollectArrays(s{i})]; %#ok<AGROW>
    end
elseif isstruct(s)
    for e = 1:numel(s)
        names = fieldnames(s);
        for i = 1:numel(names)
            if ismember(names{i}, {'meta', 'info', 'eventInfo'})
                continue
            end
            arrays = [arrays, iCollectArrays(s(e).(names{i}))]; %#ok<AGROW>
        end
    end
end
end

function iVerifyTransposed(testCase, original, transposed, label)
%IVERIFYTRANSPOSED permute(result(original), swap Y/X) == result(transposed).
if nargin < 4
    label = '';
end
expected = permute(original, [2 1 3:ndims(original)]);
testCase.verifyEqual(size(expected), size(transposed), label);
if ~isequal(size(expected), size(transposed))
    return
end
testCase.verifyEqual(isnan(expected), isnan(transposed), label);
expected(isnan(expected)) = 0;
transposed(isnan(transposed)) = 0;
testCase.verifyEqual(double(transposed), double(expected), ...
    'AbsTol', 1e-5, 'RelTol', 1e-5, label);
end

function roi = iRectROI(name, mask)
%IRECTROI Schema-valid polygon ROI covering the bounding box of MASK.
[r, c] = find(mask);
rows = min(r):max(r);
cols = min(c):max(c);
vertices = [min(cols) - 0.5, min(rows) - 0.5; max(cols) + 0.5, min(rows) - 0.5; ...
    max(cols) + 0.5, max(rows) + 0.5; min(cols) - 0.5, max(rows) + 0.5];
stats = struct('computedOn', datetime('now'), 'NPixels', nnz(mask), ...
    'areaPx2', nnz(mask), 'areaMM2', [], 'centroidXY_px', [], 'centroidXY_mm', [], ...
    'distanceFromOrigin_px', [], 'distanceFromOrigin_mm', [], 'spatialMean', [], ...
    'spatialStd', [], 'spatialMedian', [], 'spatialMin', [], 'spatialMax', []);
roi = struct('name', name, 'type', 'polygon', 'DOC', datetime('now'), ...
    'modifiedOn', datetime('now'), 'color', [1 0 0], 'notes', '', ...
    'geometry', struct('polyshape', polyshape(vertices(:, 1), vertices(:, 2)), ...
        'verticesXY_px', vertices, 'ROIType', 'polygon', ...
        'ROIParameters', struct('ROIType', 'polygon', 'Position', vertices, ...
            'Vertices', vertices)), ...
    'mask', logical(mask), 'stats', stats);
end
