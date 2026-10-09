classdef TestGenCorrelationMatrix < matlab.unittest.TestCase
    %TESTGENCORRELATIONMATRIX Unit tests for genCorrelationMatrix.
    %
    % Coverage (DFR-20260819-003):
    %   - pipelineInfo contract
    %   - the three CorrAlgorithm modes (centroid_vs_centroid, avg_vs_avg,
    %     centroid_vs_agg, the latter across all four SpatialAggFcn values)
    %   - the b_FisherZ_transform option, including the +-0.998 truncation
    %   - the SPC-map file output (b_genSPCMaps) and the UMT-struct output
    %   - streaming a .dat across several X slabs against a direct reference
    %   - rejection of array, .umt, and event-split inputs (.dat-only input)
    %   - .roi loading through loadROIFile and rejection of pre-.roi files
    %   - the degenerate-ROI guard (an all-false mask no longer aborts with a
    %     bare index error)
    %   - confirmation of audit finding P1-7: correlations now omit NaN
    %     (fully-masked traces are reported as NaN without poisoning other
    %     ROIs; partially-masked traces use the pairwise-complete estimator)
    %     instead of corrcoef's whole-matrix-NaN behavior.
    %
    % Fixture policy:
    %   - All tests use small, fully synthetic single-pixel or few-pixel ROIs
    %     with deterministic, non-orthogonal signals so expected correlation
    %     values can be computed independently with the builtin CORR and
    %     compared directly. No external fixture data is required.

    properties
        TempFolder
    end

    properties (TestParameter)
        aggFcn = {'mean', 'median', 'min', 'max'}
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = genCorrelationMatrix('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyEqual(info.name, 'genCorrelationMatrix');

            algoParam = info.parameters(strcmp({info.parameters.name}, 'CorrAlgorithm'));
            testCase.verifyEqual(algoParam.allowed, ...
                {'centroid_vs_centroid','centroid_vs_agg','avg_vs_avg'});

            aggParam = info.parameters(strcmp({info.parameters.name}, 'SpatialAggFcn'));
            testCase.verifyEqual(aggParam.allowed, {'mean','max','min','median'});

            testCase.verifyEqual({info.outputs.name}, {'outData', 'spcFile'});
            out = info.outputs(strcmp({info.outputs.name}, 'outData'));
            testCase.verifyTrue(out.isData);
            spcOut = info.outputs(strcmp({info.outputs.name}, 'spcFile'));
            testCase.verifyFalse(spcOut.isData);
            testCase.verifyFalse(spcOut.isRequired);
        end

        function testCentroidVsCentroidAndAvgVsAvgAgreeForSinglePixelROIs(testCase)
            [imgData, roiA, roiB, sigA, sigB] = iBuildTwoSinglePixelROIs();
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [4 4]);

            expected = corr(sigA(:), sigB(:));

            for algo = {'centroid_vs_centroid', 'avg_vs_avg'}
                [rho, labels] = iRunAndReadCorr(saveFolder, imgData, ...
                    'CorrAlgorithm', algo{1});
                ia = find(strcmp(labels, 'ROI_A'), 1);
                ib = find(strcmp(labels, 'ROI_B'), 1);

                testCase.verifyEqual(double(rho(ia, ib)), double(expected), 'AbsTol', 1e-4);
                testCase.verifyEqual(double(rho(ia, ia)), double(1), 'AbsTol', 1e-4);
            end
        end

        function testCentroidVsAggMatchesManualAggregation(testCase, aggFcn)
            [imgData, roiSeed, roiTarget, sigS, pixelSigs] = iBuildSeedAndThreePixelTarget();
            saveFolder = iSaveROIs(testCase, [roiSeed, roiTarget], [5 5]);

            perPixelRho = cellfun(@(s) corr(sigS(:), s(:)), pixelSigs);
            expected = feval(aggFcn, perPixelRho);

            [rho, labels] = iRunAndReadCorr(saveFolder, imgData, ...
                'CorrAlgorithm', 'centroid_vs_agg', 'SpatialAggFcn', aggFcn);
            iSeed = find(strcmp(labels, 'Seed'), 1);
            iTarget = find(strcmp(labels, 'Target'), 1);

            testCase.verifyEqual(double(rho(iSeed, iTarget)), double(expected), 'AbsTol', 1e-4);
        end

        function testFisherZTransformTruncatesBeforeAtanh(testCase)
            [imgData, roiA, roiB] = iBuildTwoIdenticalSinglePixelROIs();
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [3 3]);

            [rhoNoZ, labelsNoZ] = iRunAndReadCorr(saveFolder, imgData);
            ia = find(strcmp(labelsNoZ, 'ROI_A'), 1);
            ib = find(strcmp(labelsNoZ, 'ROI_B'), 1);
            testCase.verifyEqual(double(rhoNoZ(ia, ib)), double(1), 'AbsTol', 1e-6);

            [rhoZ, ~] = iRunAndReadCorr(saveFolder, imgData, ...
                'b_FisherZ_transform', true);
            testCase.verifyEqual(double(rhoZ(ia, ib)), double(atanh(0.998)), 'AbsTol', 1e-3);
        end

        function testReturnsUMTStructAndSPCMapsAsSecondOutput(testCase)
            [imgData, roiA, roiB] = iBuildTwoIdenticalSinglePixelROIs();
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [3 3]);
            datFile = iWriteDat(saveFolder, imgData);

            [outData, spcFile] = genCorrelationMatrix(datFile, saveFolder, ...
                'b_genSPCMaps', true);

            validateUMTStruct(outData, 'requireEventInfo', false);
            testCase.verifyEqual(lower(char(string(outData.kind))), 'roi');
            testCase.verifyEqual(fieldnames(outData.data), {'CorrMatrix'});
            testCase.verifyEqual(outData.data.CorrMatrix.dimNames, {'ROI','ROI'});
            testCase.verifyEqual(cellstr(string(outData.labels.ROI(:))).', {'ROI_A','ROI_B'});
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'corrMatrix.umt')), ...
                'The correlation matrix is returned, not saved by the function.');

            testCase.verifyEqual(spcFile, 'corrMatrix_SPCMaps.umt');
            testCase.verifyTrue(isfile(fullfile(saveFolder, spcFile)));

            spc = load(fullfile(saveFolder, spcFile), '-mat');
            testCase.verifyEqual(char(string(spc.kind)), 'image');
            testCase.verifyEqual(sort(fieldnames(spc.data)), sort({'ROI_A'; 'ROI_B'}));

            entryA = spc.data.ROI_A;
            testCase.verifyEqual(cellstr(string(entryA.dimNames)), {'Y','X'});
            testCase.verifyEqual(size(entryA.value), [3 3]);

            [~, spcFileOff] = genCorrelationMatrix(datFile, saveFolder);
            testCase.verifyEmpty(spcFileOff);
        end

        function testStreamingAcrossSlabsMatchesReference(testCase)
            % 32 x 64 x 4096 single (34 MB) exceeds the 128 MB slab budget
            % once the working-copy factor is applied, so the SPC run reads
            % two slabs, and the T ROI straddles the slab boundary (columns
            % 42|43). Runs without SPC read only the ROI columns, which are
            % not contiguous.
            imageSizeYX = [32 64];
            Nt = 4096;
            stream = RandStream('mt19937ar', 'Seed', 11);
            imgData = randn(stream, [imageSizeYX, Nt], 'single') + ...
                0.5 * randn(stream, [1 1 Nt], 'single');

            maskT = false(imageSizeYX); maskT(10:11, 42:43) = true;
            rois = [iMakeSinglePixelROIStruct('S1', imageSizeYX, [3 5]), ...
                    iMakeSinglePixelROIStruct('S2', imageSizeYX, [20 60]), ...
                    iMakeROIStructFromMask('T', maskT)];
            saveFolder = iSaveROIs(testCase, rois, imageSizeYX);
            datFile = iWriteDat(saveFolder, imgData);

            trace = @(y, x) double(squeeze(imgData(y, x, :)));
            sigS1 = trace(3, 5);
            sigS2 = trace(20, 60);
            tPixels = [trace(10, 42), trace(11, 42), trace(10, 43), trace(11, 43)];

            iS1 = 1; iS2 = 2; iT = 3;

            outCC = genCorrelationMatrix(datFile, saveFolder);
            testCase.verifyEqual(double(outCC.data.CorrMatrix.value(iS1, iS2)), ...
                corr(sigS1, sigS2), 'AbsTol', 1e-4);

            outAvg = genCorrelationMatrix(datFile, saveFolder, 'CorrAlgorithm', 'avg_vs_avg');
            testCase.verifyEqual(double(outAvg.data.CorrMatrix.value(iS1, iT)), ...
                corr(sigS1, mean(tPixels, 2)), 'AbsTol', 1e-4);

            expectedAgg = mean(arrayfun(@(k) corr(sigS1, tPixels(:, k)), 1:4));
            outAgg = genCorrelationMatrix(datFile, saveFolder, ...
                'CorrAlgorithm', 'centroid_vs_agg');
            testCase.verifyEqual(double(outAgg.data.CorrMatrix.value(iS1, iT)), ...
                expectedAgg, 'AbsTol', 1e-4);

            % Both at once: the aggregate and the SPC maps share one pass.
            [outBoth, spcFile] = genCorrelationMatrix(datFile, saveFolder, ...
                'CorrAlgorithm', 'centroid_vs_agg', 'b_genSPCMaps', true);
            testCase.verifyEqual(double(outBoth.data.CorrMatrix.value(iS1, iT)), ...
                expectedAgg, 'AbsTol', 1e-4);

            spc = load(fullfile(saveFolder, spcFile), '-mat');
            testCase.verifyEqual(double(spc.data.S1.value(3, 5)), 1, 'AbsTol', 1e-4);
            testCase.verifyEqual(double(spc.data.S1.value(10, 43)), ...
                corr(sigS1, trace(10, 43)), 'AbsTol', 1e-4);
            testCase.verifyEqual(double(spc.data.S2.value(32, 64)), ...
                corr(sigS2, trace(32, 64)), 'AbsTol', 1e-4);
            testCase.verifyEqual(double(spc.data.S2.value(1, 1)), ...
                corr(sigS2, trace(1, 1)), 'AbsTol', 1e-4);
        end

        function testRejectsArrayUMTFileAndEventSplitInput(testCase)
            [imgData, roiA, roiB] = iBuildTwoIdenticalSinglePixelROIs();
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [3 3]);

            testCase.verifyError(@() genCorrelationMatrix(imgData, saveFolder), ...
                'Umitoolbox:genCorrelationMatrix:UnsupportedInputType');

            umtFile = fullfile(saveFolder, 'input.umt');
            fclose(fopen(umtFile, 'w'));
            testCase.verifyError(@() genCorrelationMatrix(umtFile, saveFolder), ...
                'Umitoolbox:genCorrelationMatrix:UnsupportedInputFile');

            byEvent = fullfile(saveFolder, 'byEvent.dat');
            saveData(byEvent, cat(4, imgData, imgData), ...
                'DimNames', {'Y','X','T','E'}, 'FrameRateHz', 10);
            testCase.verifyError(@() genCorrelationMatrix(byEvent, saveFolder), ...
                'Umitoolbox:genCorrelationMatrix:unsupportedLayout');
        end

        function testRoiFileLoadingRejectsPreRoiFiles(testCase)
            [imgData, roiA, roiB] = iBuildTwoIdenticalSinglePixelROIs();
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [3 3]);

            datFile = iWriteDat(saveFolder, imgData);

            legacyFile = fullfile(saveFolder, 'legacy.mat');
            fid = fopen(legacyFile, 'w');
            testCase.assertNotEqual(fid, -1);
            fclose(fid);

            testCase.verifyError(@() genCorrelationMatrix(datFile, saveFolder, ...
                'ROImasks_filename', legacyFile), ...
                'Umitoolbox:genCorrelationMatrix:UnsupportedROIFile');
        end

        function testMissingRoiFileThrows(testCase)
            [imgData, roiA, roiB] = iBuildTwoIdenticalSinglePixelROIs();
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [3 3]);
            datFile = iWriteDat(saveFolder, imgData);

            testCase.verifyError(@() genCorrelationMatrix(datFile, saveFolder, ...
                'ROImasks_filename', 'doesNotExist.roi'), ...
                'Umitoolbox:genCorrelationMatrix:MissingROIFile');
        end

        function testDegenerateROIGuardOnEmptyMask(testCase)
            imageSizeYX = [4 4];
            Nt = 10;
            imgData = single(reshape(1:(imageSizeYX(1)*imageSizeYX(2)*Nt), ...
                [imageSizeYX, Nt]));

            emptyMask = false(imageSizeYX);
            roiEmpty = iMakeSinglePixelROIStruct('EmptyROI', imageSizeYX, [1 1]);
            roiEmpty.mask = emptyMask;
            roiEmpty.stats.NPixels = 0;
            roiEmpty.stats.areaPx2 = 0;

            saveFolder = iSaveROIs(testCase, roiEmpty, imageSizeYX);
            datFile = iWriteDat(saveFolder, imgData);

            % An empty mask must be rejected with a clear error instead of a
            % bare "index exceeds array bounds" failure from cIdx(1).
            testCase.verifyError(@() genCorrelationMatrix(datFile, saveFolder), ...
                'Umitoolbox:genCorrelationMatrix:EmptyROIMask');
        end

        function testNaNHandlingUsesPairwiseCorrelationForPartiallyMaskedTraces(testCase)
            Nt = 30;
            t = 0:(Nt-1);
            sigA = sin(2*pi*t/9);
            sigB = cos(2*pi*t/9) + 0.3*sin(2*pi*t/5);
            sigA(1, [3 10 20]) = NaN;

            [imgData, roiA, roiB] = iBuildTwoSinglePixelROIsFromSignals(sigA, sigB);
            saveFolder = iSaveROIs(testCase, [roiA, roiB], [4 4]);

            [rho, labels] = iRunAndReadCorr(saveFolder, imgData);
            ia = find(strcmp(labels, 'ROI_A'), 1);
            ib = find(strcmp(labels, 'ROI_B'), 1);

            expected = corr(sigA(:), sigB(:), 'rows', 'pairwise');

            % Before the P1-7 fix, bare CORRCOEF returned NaN for any pair
            % touching a NaN sample -- confirm that no longer happens and the
            % pairwise-complete estimator is used instead.
            testCase.verifyFalse(isnan(rho(ia, ib)));
            testCase.verifyEqual(double(rho(ia, ib)), double(expected), 'AbsTol', 1e-4);
        end

        function testNaNHandlingOmitsFullyMaskedPixelFromAggregate(testCase)
            Nt = 30;
            t = 0:(Nt-1);
            sigS = sin(2*pi*t/8);
            sigT = cos(2*pi*t/8) + 0.2*sin(2*pi*t/3);

            imageSizeYX = [4 4];
            imgData = zeros([imageSizeYX, Nt], 'single');
            imgData(1, 1, :) = sigS;
            imgData(2, 2, :) = NaN;
            imgData(3, 3, :) = sigT;

            maskSeed = false(imageSizeYX); maskSeed(1, 1) = true;
            maskTarget = false(imageSizeYX); maskTarget(2, 2) = true; maskTarget(3, 3) = true;
            roiSeed = iMakeROIStructFromMask('Seed', maskSeed);
            roiTarget = iMakeROIStructFromMask('Target', maskTarget);

            saveFolder = iSaveROIs(testCase, [roiSeed, roiTarget], imageSizeYX);

            [rho, labels] = iRunAndReadCorr(saveFolder, imgData, ...
                'CorrAlgorithm', 'centroid_vs_agg', 'SpatialAggFcn', 'mean');
            iSeed = find(strcmp(labels, 'Seed'), 1);
            iTarget = find(strcmp(labels, 'Target'), 1);

            % The fully-NaN target pixel must be excluded from the omitnan
            % mean instead of making the whole aggregate NaN.
            expected = corr(sigS(:), sigT(:));
            testCase.verifyEqual(double(rho(iSeed, iTarget)), double(expected), 'AbsTol', 1e-4);
            testCase.verifyEqual(double(rho(iSeed, iSeed)), double(1), 'AbsTol', 1e-4);
        end
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function saveFolder = iSaveROIs(testCase, rois, imageSizeYX)
%ISAVEROIS Save an ROI struct array as myROI.roi in a fresh subfolder.
saveFolder = testCase.TempFolder;
ROIFile = createROIFile(imageSizeYX, 'ROIs', rois);
saveROIFile(fullfile(saveFolder, 'myROI.roi'), ROIFile);
end

function datFile = iWriteDat(saveFolder, imgData)
%IWRITEDAT Write a Y-X-T array as the .dat input of genCorrelationMatrix.
datFile = fullfile(saveFolder, 'input.dat');
saveData(datFile, single(imgData), 'DimNames', {'Y','X','T'}, 'FrameRateHz', 10);
end

function [rho, labels] = iRunAndReadCorr(saveFolder, imgData, varargin)
%IRUNANDREADCORR Run genCorrelationMatrix and return the matrix and ROI labels.
outData = genCorrelationMatrix(iWriteDat(saveFolder, imgData), saveFolder, varargin{:});
rho = outData.data.CorrMatrix.value;
labels = outData.labels.ROI;
end

function roi = iMakeSinglePixelROIStruct(name, imageSizeYX, rc)
%IMAKESINGLEPIXELROISTRUCT Build one schema-valid single-pixel ROI.
mask = false(imageSizeYX);
mask(rc(1), rc(2)) = true;
roi = iMakeROIStructFromMask(name, mask);
end

function roi = iMakeROIStructFromMask(name, mask)
%IMAKEROISTRUCTFROMMASK Build one schema-valid ROI struct from a mask.
%
% The bounding box below is only used to produce a schema-valid polyshape;
% the mask field (not the polygon) is authoritative for pixel membership.
[rows, cols] = find(mask);
if isempty(rows)
    rows = 1; cols = 1;
end
x0 = min(cols) - 0.5; x1 = max(cols) + 0.5;
y0 = min(rows) - 0.5; y1 = max(rows) + 0.5;
vertices = [x0 y0; x1 y0; x1 y1; x0 y1];
pgon = polyshape(vertices(:,1), vertices(:,2));

stats = struct( ...
    'computedOn', datetime('now'), ...
    'NPixels', nnz(mask), ...
    'areaPx2', nnz(mask), ...
    'areaMM2', [], ...
    'centroidXY_px', [], ...
    'centroidXY_mm', [], ...
    'distanceFromOrigin_px', [], ...
    'distanceFromOrigin_mm', [], ...
    'spatialMean', [], ...
    'spatialStd', [], ...
    'spatialMedian', [], ...
    'spatialMin', [], ...
    'spatialMax', []);

roi = struct( ...
    'name', name, ...
    'type', 'polygon', ...
    'DOC', datetime('now'), ...
    'modifiedOn', datetime('now'), ...
    'color', [1 0 0], ...
    'notes', '', ...
    'geometry', struct( ...
        'polyshape', pgon, ...
        'verticesXY_px', vertices, ...
        'ROIType', 'polygon', ...
        'ROIParameters', struct('ROIType', 'polygon', ...
            'Position', vertices, 'Vertices', vertices)), ...
    'mask', logical(mask), ...
    'stats', stats);
end

function [imgData, roiA, roiB, sigA, sigB] = iBuildTwoSinglePixelROIs()
%IBUILDTWOSINGLEPIXELROIS Two 1-pixel ROIs with distinct, non-orthogonal signals.
Nt = 25;
t = 0:(Nt-1);
sigA = sin(2*pi*t/9);
sigB = cos(2*pi*t/9) + 0.15*t/Nt;
[imgData, roiA, roiB] = iBuildTwoSinglePixelROIsFromSignals(sigA, sigB);
end

function [imgData, roiA, roiB] = iBuildTwoSinglePixelROIsFromSignals(sigA, sigB)
%IBUILDTWOSINGLEPIXELROISFROMSIGNALS Place two given traces at two pixels.
imageSizeYX = [4 4];
Nt = numel(sigA);
imgData = zeros([imageSizeYX, Nt], 'single');
imgData(1, 1, :) = sigA;
imgData(3, 3, :) = sigB;

roiA = iMakeSinglePixelROIStruct('ROI_A', imageSizeYX, [1 1]);
roiB = iMakeSinglePixelROIStruct('ROI_B', imageSizeYX, [3 3]);
end

function [imgData, roiA, roiB] = iBuildTwoIdenticalSinglePixelROIs()
%IBUILDTWOIDENTICALSINGLEPIXELROIS Two 1-pixel ROIs sharing the same trace.
Nt = 30;
t = 0:(Nt-1);
sig = sin(2*pi*t/7);
imageSizeYX = [3 3];
imgData = zeros([imageSizeYX, Nt], 'single');
imgData(1, 1, :) = sig;
imgData(2, 2, :) = sig;

roiA = iMakeSinglePixelROIStruct('ROI_A', imageSizeYX, [1 1]);
roiB = iMakeSinglePixelROIStruct('ROI_B', imageSizeYX, [2 2]);
end

function [imgData, roiSeed, roiTarget, sigS, pixelSigs] = iBuildSeedAndThreePixelTarget()
%IBUILDSEEDANDTHREEPIXELTARGET One 1-pixel seed ROI vs. a 3-pixel target ROI.
imageSizeYX = [5 5];
Nt = 25;
t = 0:(Nt-1);
sigS = sin(2*pi*t/6);
sig1 = cos(2*pi*t/6);
sig2 = cos(2*pi*t/6 + 0.3) + 0.1*t/Nt;
sig3 = -sin(2*pi*t/4);
pixelSigs = {sig1, sig2, sig3};

imgData = zeros([imageSizeYX, Nt], 'single');
imgData(1, 1, :) = sigS;
imgData(3, 1, :) = sig1;
imgData(3, 2, :) = sig2;
imgData(3, 3, :) = sig3;

maskTarget = false(imageSizeYX);
maskTarget(3, 1) = true; maskTarget(3, 2) = true; maskTarget(3, 3) = true;

roiSeed = iMakeSinglePixelROIStruct('Seed', imageSizeYX, [1 1]);
roiTarget = iMakeROIStructFromMask('Target', maskTarget);
end
