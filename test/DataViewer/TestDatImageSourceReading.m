classdef TestDatImageSourceReading < matlab.unittest.TestCase
    %TESTDATIMAGESOURCEREADING Every DatImageSource read path against loadData.
    %
    %   Runs on copies of the Y-X-T single reference files (headered, legacy
    %   sidecar; the AcqInfos-bound one is rejected) in a temporary folder: the constructor may
    %   write rig metadata into AcqInfos.mat, so reference files are never
    %   opened in place. A small cache budget forces a partial spatial cache
    %   so both cached and uncached reads are exercised.

    properties (Constant)
        RigWarning = 'Umitoolbox:DatImageSource:rigAssociationSkipped'
    end

    properties
        FixtureRoot
    end

    properties (TestParameter)
        kind = struct( ...
            'header', struct('rel', {{'v1', 'yxt_single.dat'}}, 'extra', {{}}), ...
            'legacySidecar', struct('rel', {{'legacy_sidecar', 'green.dat'}}, 'extra', {{'green.mat'}}))
        rejected = struct( ...
            'yx_uint16', struct('file', 'yx_uint16.dat', 'id', 'DatImageSource:UnsupportedDatatype'))
        eventSplit = struct( ...
            'yxte_single', struct('file', 'yxte_single.dat'), ...
            'yxe_single', struct('file', 'yxe_single.dat'))
    end

    methods (TestClassSetup)
        function setupPaths(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(testFolder));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(repoRoot, 'GUI', 'DataViewer')));
            testCase.FixtureRoot = fullfile(repoRoot, 'test', 'subFunc', 'datHeader', 'fixtures');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat')), ...
                'Reference files missing: run makeDatFixturesV1 once.');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'green.dat')), ...
                'Legacy reference file missing: run makeLegacySidecarFixture once.');
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture(testCase.RigWarning));
        end
    end

    methods (Test)
        function sizeAndMetadata(testCase, kind)
            [src, ref] = testCase.openCopy(kind, []);
            testCase.verifyEqual(src.getSize(), [size(ref), 1]);
            testCase.verifyEqual(src.Precision, 'single');
            testCase.verifyEqual(src.FrameRateHz, 10);
        end

        function framesMatchLoadData(testCase, kind)
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            testCase.assertTrue(src.hasPartialTemporalCache(), 'cache should be partial');
            for t = [1 2 size(ref, 3)]
                testCase.verifyEqual(src.getFrame(t), ref(:, :, t), sprintf('frame %d', t));
            end
        end

        function framesFromFullCacheMatchLoadData(testCase, kind)
            [src, ref] = testCase.openCopy(kind, []);
            testCase.assertFalse(src.hasPartialTemporalCache(), 'cache should hold the whole image');
            testCase.verifyEqual(src.getFrame(3), ref(:, :, 3));
        end

        function blocksMatchLoadData(testCase, kind)
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            yIn = src.CacheYRange;
            xIn = src.CacheXRange;
            testCase.verifyEqual(src.getFrameBlock(2, yIn, xIn), ref(yIn, xIn, 2), 'inside cache');
            testCase.verifyEqual(src.getFrameBlock(4, 1:size(ref, 1), 1:size(ref, 2)), ...
                ref(:, :, 4), 'whole frame (outside cache)');
            testCase.verifyEqual(src.getFrameBlock(1, 2:5, 3:5), ref(2:5, 3:5, 1), 'partial block');
        end

        function pixelTracesInsideAndOutsideCache(testCase, kind)
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            yIn = src.CacheYRange(1);
            xIn = src.CacheXRange(1);
            [trace, status] = src.getPixelTrace(yIn, xIn);
            testCase.verifyEqual(status, 'ok');
            testCase.verifyEqual(trace, squeeze(ref(yIn, xIn, :)));

            [yOut, xOut] = testCase.pixelOutsideCache(src, size(ref));
            [trace, status] = src.getPixelTrace(yOut, xOut);
            testCase.verifyEqual(status, 'cache_rebuilt');
            testCase.verifyEqual(trace, squeeze(ref(yOut, xOut, :)));
        end

        function cacheFillMatchesLoadData(testCase, kind)
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            src.updateCacheAround(size(ref, 1), size(ref, 2));
            testCase.verifyEqual(src.CacheData, ref(src.CacheYRange, src.CacheXRange, :));
            src.updateCacheAround(1, 1);
            testCase.verifyEqual(src.CacheData, ref(src.CacheYRange, src.CacheXRange, :));
        end

        function roiTracesFromCacheAndDirect(testCase, kind)
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            [Ny, Nx, Nt] = size(ref);

            cacheMask = false(Ny, Nx);
            cacheMask(src.CacheYRange, src.CacheXRange(1)) = true;
            fullMask = false(Ny, Nx);
            fullMask([1 Ny], [1 Nx]) = true;

            [traces, mode] = src.getROIMeanTraceMatrix(cacheMask, 1:Nt);
            testCase.verifyEqual(mode, 'cache');
            testCase.verifyEqual(traces, testCase.roiMeans(ref, cacheMask, 1:Nt), 'AbsTol', 1e-5);

            [traces, mode] = src.getROIMeanTraceMatrix({cacheMask, fullMask}, 1:Nt);
            testCase.verifyEqual(mode, 'direct_full_frame');
            expected = [testCase.roiMeans(ref, cacheMask, 1:Nt); testCase.roiMeans(ref, fullMask, 1:Nt)];
            testCase.verifyEqual(traces, expected, 'AbsTol', 1e-5);
        end

        function roiEventMatrixWithRepeatedAndMissingFrames(testCase, kind)
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            [Ny, Nx, ~] = size(ref);
            fullMask = false(Ny, Nx);
            fullMask(:, [1 Nx]) = true;
            frameIdx = [1 2 3; 3 3 NaN];

            [traces, mode] = src.getROIMeanTraceMatrix(fullMask, frameIdx);
            testCase.verifyEqual(mode, 'direct_full_frame');
            testCase.verifySize(traces, [1 2 3]);
            testCase.verifyEqual(squeeze(traces(1, 1, :)).', testCase.roiMeans(ref, fullMask, [1 2 3]), 'AbsTol', 1e-5);
            testCase.verifyEqual(squeeze(traces(1, 2, 1:2)).', testCase.roiMeans(ref, fullMask, [3 3]), 'AbsTol', 1e-5);
            testCase.verifyTrue(isnan(traces(1, 2, 3)));
        end

        function noFileLeftOpen(testCase, kind)
            before = fopen('all');
            [src, ref] = testCase.openCopy(kind, testCase.partialCacheBytes());
            src.getFrame(2);
            src.getFrameBlock(3, 1:2, 1:2);
            [yOut, xOut] = testCase.pixelOutsideCache(src, size(ref));
            src.getPixelTrace(yOut, xOut);
            src.getROIMeanTraceMatrix(true(size(ref, 1), size(ref, 2)), 1:size(ref, 3));
            testCase.verifyEqual(sort(fopen('all')), sort(before));
        end

        function unsupportedLayoutsAreRejected(testCase, rejected)
            folder = testCase.tempFolder();
            f = fullfile(folder, rejected.file);
            copyfile(fullfile(testCase.FixtureRoot, 'v1', rejected.file), f);
            testCase.verifyError(@() DatImageSource(f), rejected.id);
        end

        function eventSplitFilesReadLikeLoadData(testCase, eventSplit)
            % .dat header Phase 8b: Y-X-T-E and Y-X-E files open; (t, e)
            % reads, blocks, and pixel traces match loadData, from a partial
            % and from a full cache.
            folder = testCase.tempFolder();
            f = fullfile(folder, eventSplit.file);
            copyfile(fullfile(testCase.FixtureRoot, 'v1', eventSplit.file), f);
            ref = loadData(f);
            md = loadMetaData(f);
            nT = max(1, datAxisSize(md, 'T'));
            nE = datAxisSize(md, 'E');
            ref = reshape(ref, size(ref, 1), size(ref, 2), nT, nE);

            for cacheBytes = {testCase.partialCacheBytes(), []}
                src = DatImageSource(f, 'maxCacheBytes', cacheBytes{1});
                testCase.verifyEqual(src.getSize(), [size(ref, 1), size(ref, 2), nT, nE]);
                testCase.verifyTrue(src.HasEventAxis);
                for e = 1:nE
                    for t = 1:nT
                        testCase.verifyEqual(src.getFrame(t, e), ref(:, :, t, e));
                        testCase.verifyEqual(src.getFrameBlock(t, 1:2, 2:3, e), ref(1:2, 2:3, t, e));
                    end
                    testCase.verifyEqual(src.getPixelTrace(2, 3, e), squeeze(ref(2, 3, :, e)));
                end
                testCase.verifyEqual(src.getFrame(1), ref(:, :, 1, 1), 'event index defaults to 1');
                testCase.verifyError(@() src.getFrame(1, nE + 1), 'MATLAB:DatImageSource:notLessEqual');
                clear src
            end
        end

        function eventSplitFileWithoutEventsShowsOneCondition(testCase, eventSplit)
            % No events.mat: the E slices are repetitions of one condition.
            folder = testCase.tempFolder();
            f = fullfile(folder, eventSplit.file);
            copyfile(fullfile(testCase.FixtureRoot, 'v1', eventSplit.file), f);
            nE = datAxisSize(loadMetaData(f), 'E');

            src = DatImageSource(f);
            testCase.verifyEqual(src.EventMapping.status, 'noEvents');
            info = src.getEventInfo();
            testCase.verifyEqual(info.eventID, ones(nE, 1));
            testCase.verifyEqual(info.repetitionIndex, (1:nE).');
            testCase.verifyEqual(info.eventAxisMode, 'instances');
            testCase.verifyEqual(info.selected, true(nE, 1));
            clear src
        end

        function continuousFileHasNoEventAxis(testCase)
            [src, ~] = testCase.openCopy(struct('rel', {{'v1', 'yxt_single.dat'}}, 'extra', {{}}), []);
            testCase.verifyFalse(src.HasEventAxis);
            testCase.verifyEqual(src.getEventInfo(), struct());
        end

        function singleFrameFileOpensAsOneFrame(testCase)
            % .dat header Phase 6b-1: a Y-X file is shown as one frame.
            folder = testCase.tempFolder();
            f = fullfile(folder, 'frame.dat');
            img = reshape(single(1:30), 6, 5);
            writeTestDat(f, img, 10, 'DimNames', {'Y', 'X'});

            src = DatImageSource(f);
            testCase.verifyEqual(src.getSize(), [6 5 1 1]);
            testCase.verifyEqual(src.getFrame(1), img);
            testCase.verifyTrue(isnan(src.FrameRateHz));
            clear src
        end

        function acqInfosBoundFileIsRejected(testCase)
            folder = testCase.tempFolder();
            srcFolder = fullfile(testCase.FixtureRoot, 'v1', 'legacy');
            copyfile(fullfile(srcFolder, 'green.dat'), folder);
            copyfile(fullfile(srcFolder, 'AcqInfos.mat'), folder);
            testCase.verifyError(@() DatImageSource(fullfile(folder, 'green.dat')), ...
                'Umitoolbox:loadMetaData:acqInfosBoundUnsupported');
        end
    end

    methods (Access = private)
        function [src, ref] = openCopy(testCase, kind, maxCacheBytes)
            folder = testCase.tempFolder();
            srcFolder = fullfile(testCase.FixtureRoot, kind.rel{1:end-1});
            copyfile(fullfile(srcFolder, kind.rel{end}), folder);
            for k = 1:numel(kind.extra)
                copyfile(fullfile(srcFolder, kind.extra{k}), folder);
            end
            f = fullfile(folder, kind.rel{end});
            ref = loadData(f);
            if isempty(maxCacheBytes)
                src = DatImageSource(f);
            else
                src = DatImageSource(f, 'maxCacheBytes', maxCacheBytes);
            end
        end

        function bytes = partialCacheBytes(~)
            % 6 pixel traces of 4 single frames: a partial cache of the 6x5 image.
            bytes = 6 * 4 * 4;
        end

        function [y, x] = pixelOutsideCache(~, src, sz)
            [Y, X] = ndgrid(1:sz(1), 1:sz(2));
            outside = ~(ismember(Y, src.CacheYRange) & ismember(X, src.CacheXRange));
            idx = find(outside, 1);
            [y, x] = ind2sub(sz(1:2), idx);
        end

        function m = roiMeans(~, ref, mask, frames)
            m = zeros(1, numel(frames));
            for k = 1:numel(frames)
                frame = double(ref(:, :, frames(k)));
                m(k) = mean(frame(mask));
            end
        end

        function folder = tempFolder(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end
end
