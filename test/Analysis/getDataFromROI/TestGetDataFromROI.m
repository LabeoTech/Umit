classdef TestGetDataFromROI < matlab.unittest.TestCase
    %TESTGETDATAFROMROI Unit tests for getDataFromROI.
    %
    % Coverage (DFR-20260819-002):
    %   - pipelineInfo contract
    %   - top-level eventInfo propagation, including baselinePeriod
    %   - per-entry meta.FrameRateHz propagation (UMT meta and .dat header)
    %   - SpatialAggFcn='none' {ROI,Pixel,...} output and its NaN padding to
    %     the largest ROI
    %   - .roi loading through loadROIFile and rejection of pre-.roi files
    %   - ROI filename resolution against SaveFolder (bare name vs. full path)
    %   - the round trip split_data_by_event -> getDataFromROI ->
    %     calculateResponseFeatures
    %   - the input-form matrix: .dat Y-X/Y-X-T/Y-X-F and UMT Y-X/Y-X-T/Y-X-E/
    %     Y-X-T-E/Y-X-F (each with and without eventInfo where the schema
    %     permits), as a standing guard against re-introducing a coupling to
    %     event-split producers
    %   - the streamed .dat path equals the in-RAM UMT path for every
    %     aggregation, also when the memory mock forces many frame blocks
    %   - Y-X-F input gives {ROI,F} with the F labels carried over
    %   - rejection of numeric arrays, other file types, ROI-axis images, size
    %     mismatches, and 'none' on an F axis
    %
    % Fixture policy:
    %   - Most tests use a small fully synthetic 6x6 image with two
    %     differently-sized ROIs (9 px and 12 px), which needs no external
    %     data.
    %   - The round-trip test needs a real EventsManager-compatible
    %     AcqInfos.mat (EventsManager reads camera/AIN fields that a minimal
    %     synthetic struct does not provide) and reuses the existing
    %     test/Analysis/TestingData_with_events fixture, matching sibling
    %     tests such as TestApplyAggregateFunction and TestGenAmplitudeMaps.

    properties
        TempFolder
        ImageSizeYX = [6 6]
        FrameRateHz = 10
        NFrames = 12
    end

    properties (TestParameter)
        inputForm = {'datYX', 'datYXT', 'datYXF', 'umtYX', 'umtYXT', ...
            'umtYXF', 'umtYXEnoEvt', 'umtYXEwithEvt', 'umtYXTEnoEvt', ...
            'umtYXTEwithEvt'}
        aggFcn = {'mean','max','min','median','mode','sum','std'}
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;

            AcqInfoStream = struct( ...
                'Height', testCase.ImageSizeYX(1), ...
                'Width', testCase.ImageSizeYX(2), ...
                'Length', testCase.NFrames, ...
                'FrameRateHz', testCase.FrameRateHz, ...
                'Datatype', 'single', ...
                'dim_names', {{'Y','X','T'}}, ...
                'ExposureMsec', 5, ...
                'MultiCam', false);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            roi1 = iMakeRectROI('ROI1', testCase.ImageSizeYX, 1:3, 1:3, [1 0 0]);
            roi2 = iMakeRectROI('ROI2', testCase.ImageSizeYX, 4:6, 3:6, [0 1 0]);
            ROIFile = createROIFile(testCase.ImageSizeYX, 'ROIs', [roi1, roi2]);
            saveROIFile(fullfile(testCase.TempFolder, 'myROI.roi'), ROIFile);

            soloROI = iMakeRectROI('SoloROI', testCase.ImageSizeYX, 1:3, 1:3, [1 0 0]);
            soloROIFile = createROIFile(testCase.ImageSizeYX, 'ROIs', soloROI);
            saveROIFile(fullfile(testCase.TempFolder, 'soloROI.roi'), soloROIFile);
        end
    end

    methods (Test)

        function testPipelineInfo(testCase)
            info = getDataFromROI('pipelineInfo');

            testCase.verifyTrue(isstruct(info) && isscalar(info));
            testCase.verifyEqual(info.name, 'getDataFromROI');

            dataInput = info.inputs(strcmp({info.inputs.name}, 'data'));
            testCase.verifyTrue(dataInput.isData);
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');

            aggParam = info.parameters(strcmp({info.parameters.name}, 'SpatialAggFcn'));
            testCase.verifyEqual(aggParam.allowed, ...
                {'none','mean','max','min','median','mode','sum','std'});

            out = info.outputs(strcmp({info.outputs.name}, 'outData'));
            testCase.verifyTrue(out.isData);
        end

        function testEventInfoAndBaselinePeriodPropagation(testCase)
            umt = iBuildYXEImage(testCase.ImageSizeYX, 3);
            umt = appendUMTEventInfo(umt, ...
                'eventID', [1;2;1], ...
                'repetitionIndex', [1;1;2], ...
                'eventName', {'A';'B';'A'}, ...
                'eventAxisMode', 'instances', ...
                'overwrite', true);
            umt.eventInfo.baselinePeriod = 0.75;
            validateUMTStruct(umt, 'requireEventInfo', true);

            out = getDataFromROI(umt, testCase.TempFolder);

            testCase.verifyTrue(isfield(out, 'eventInfo'));
            testCase.verifyEqual(out.eventInfo, umt.eventInfo);
        end

        function testEntryMetaFrameRatePropagation(testCase)
            umt = iBuildYXTImage(testCase.ImageSizeYX, testCase.NFrames);
            umt.data.main.meta = struct('FrameRateHz', 37);

            out = getDataFromROI(umt, testCase.TempFolder);

            entryName = fieldnames(out.data);
            entryName = entryName{1};
            testCase.verifyEqual(out.data.(entryName).meta.FrameRateHz, 37);
        end

        function testFrameRateComesFromTheDataNotAcqInfos(testCase)
            % .dat header Phase 7a: genUMTStruct no longer fills
            % meta.FrameRateHz from AcqInfos.mat (present here). A UMT input
            % without its own rate gives none; a .dat input carries the rate
            % of its own header (not AcqInfos.mat's 10 Hz).
            umt = iBuildYXTImage(testCase.ImageSizeYX, testCase.NFrames);
            testCase.verifyFalse(isfield(umt.data.main, 'meta'));
            out = getDataFromROI(umt, testCase.TempFolder);
            entryName = fieldnames(out.data);
            entry = out.data.(entryName{1});
            testCase.verifyFalse(isfield(entry, 'meta') && isfield(entry.meta, 'FrameRateHz'));

            datFile = iWriteDat(testCase.TempFolder, 'rate.dat', ...
                iSyntheticYXT(testCase.ImageSizeYX, testCase.NFrames), {'Y','X','T'}, 12.5);
            out = getDataFromROI(datFile, testCase.TempFolder);
            entryName = fieldnames(out.data);
            testCase.verifyEqual(out.data.(entryName{1}).meta.FrameRateHz, 12.5);
        end

        function testSpatialAggNonePadsToLargestROIWithNaN(testCase)
            data = iSyntheticYXT(testCase.ImageSizeYX, testCase.NFrames);
            datFile = iWriteDat(testCase.TempFolder, 'in.dat', data, {'Y','X','T'}, 10);

            out = getDataFromROI(datFile, testCase.TempFolder, 'SpatialAggFcn', 'none');

            entryName = fieldnames(out.data);
            entryName = entryName{1};
            entry = out.data.(entryName);

            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'ROI','Pixel','T'});
            testCase.verifyEqual(size(entry.value), [2, 12, testCase.NFrames]);

            % ROI1 has 9 pixels: the trailing 3 pixel slots must be NaN.
            roi1Slice = squeeze(entry.value(1, :, :));
            testCase.verifyEqual(nnz(all(~isnan(roi1Slice), 2)), 9);
            testCase.verifyEqual(nnz(all(isnan(roi1Slice), 2)), 3);

            % ROI2 has 12 pixels: fully populated, no NaN padding.
            roi2Slice = squeeze(entry.value(2, :, :));
            testCase.verifyFalse(any(isnan(roi2Slice), 'all'));

            % Cross-check actual values against a direct mask extraction.
            mask1 = false(testCase.ImageSizeYX); mask1(1:3, 1:3) = true;
            data2D = reshape(data, prod(testCase.ImageSizeYX), testCase.NFrames);
            expectedROI1 = data2D(mask1(:), :);
            testCase.verifyEqual(double(roi1Slice(1:9, :)), double(expectedROI1));
        end

        function testSingleROIAggregatedYXTDoesNotCrash(testCase)
            % Regression for DFR-20260819-007: a one-ROI file paired with an
            % aggregated {Y,X,T} entry used to throw
            % Umitoolbox:genUMTStruct:invalidInput because the resulting
            % [1, NFrames] value was misidentified as a mis-oriented 1-D
            % vector.
            data = iSyntheticYXT(testCase.ImageSizeYX, testCase.NFrames);
            datFile = iWriteDat(testCase.TempFolder, 'in.dat', data, {'Y','X','T'}, 10);

            out = getDataFromROI(datFile, testCase.TempFolder, ...
                'ROImasks_filename', 'soloROI.roi', 'SpatialAggFcn', 'mean');

            entryName = fieldnames(out.data);
            entryName = entryName{1};
            entry = out.data.(entryName);

            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'ROI','T'});
            testCase.verifyEqual(size(entry.value), [1, testCase.NFrames]);

            mask = false(testCase.ImageSizeYX); mask(1:3, 1:3) = true;
            data2D = reshape(data, prod(testCase.ImageSizeYX), testCase.NFrames);
            expected = mean(double(data2D(mask(:), :)), 1);
            testCase.verifyEqual(double(entry.value), expected, 'AbsTol', 1e-4);
        end

        function testSingleROIAggregatedYXEDoesNotCrash(testCase)
            % Regression for DFR-20260819-007, {Y,X,E} entry variant.
            umt = iBuildYXEImage(testCase.ImageSizeYX, 3);

            out = getDataFromROI(umt, testCase.TempFolder, ...
                'ROImasks_filename', 'soloROI.roi', 'SpatialAggFcn', 'mean');

            entryName = fieldnames(out.data);
            entryName = entryName{1};
            entry = out.data.(entryName);

            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'ROI','E'});
            testCase.verifyEqual(size(entry.value), [1, 3]);
        end

        function testSingleROISingleEventPreservesNamedDimensions(testCase)
            % A one-ROI, one-event result is a MATLAB scalar, but it still
            % semantically represents the declared {ROI,E} dimensions.
            umt = iBuildYXEImage(testCase.ImageSizeYX, 1);
            umt = appendUMTEventInfo(umt, ...
                'eventID', 1, ...
                'repetitionIndex', 0, ...
                'eventName', {'OnlyEvent'}, ...
                'eventAxisMode', 'aggregated_repetitions');

            out = getDataFromROI(umt, testCase.TempFolder, ...
                'ROImasks_filename', 'soloROI.roi', 'SpatialAggFcn', 'mean');

            entry = out.data.main;
            testCase.verifyTrue(isscalar(entry.value));
            testCase.verifyEqual(entry.dimNames, {'ROI','E'});
            testCase.verifyEqual(out.eventInfo, umt.eventInfo);
            testCase.verifyWarningFree(@() validateUMTStruct(out));
        end

        function testRoiFilenameResolvesAgainstSaveFolderAndAbsolutePath(testCase)
            umt = iBuildYXTImage(testCase.ImageSizeYX, testCase.NFrames);

            outBare = getDataFromROI(umt, testCase.TempFolder, ...
                'ROImasks_filename', 'myROI.roi');
            outAbs = getDataFromROI(umt, testCase.TempFolder, ...
                'ROImasks_filename', fullfile(testCase.TempFolder, 'myROI.roi'));

            testCase.verifyEqual(outBare, outAbs);
        end

        function testMissingRoiFileThrows(testCase)
            umt = iBuildYXTImage(testCase.ImageSizeYX, testCase.NFrames);

            testCase.verifyError(@() getDataFromROI(umt, testCase.TempFolder, ...
                'ROImasks_filename', 'doesNotExist.roi'), ...
                'Umitoolbox:getDataFromROI:MissingROIFile');
        end

        function testPreRoiFileIsRejected(testCase)
            umt = iBuildYXTImage(testCase.ImageSizeYX, testCase.NFrames);

            legacyFile = fullfile(testCase.TempFolder, 'legacy.mat');
            fid = fopen(legacyFile, 'w');
            testCase.assertNotEqual(fid, -1);
            fclose(fid);

            testCase.verifyError(@() getDataFromROI(umt, testCase.TempFolder, ...
                'ROImasks_filename', legacyFile), ...
                'Umitoolbox:getDataFromROI:UnsupportedROIFile');
        end

        function testStreamedDatMatchesInRamUmtForEveryAggregation(testCase, aggFcn)
            % The .dat path is streamed in frame blocks and aggregated per
            % block; the UMT path aggregates in RAM. Both must agree.
            stream = RandStream('mt19937ar', 'Seed', 2);
            data = round(10 * randn(stream, [testCase.ImageSizeYX, 40], 'single'));
            data(2, 2, 5:9) = NaN;
            datFile = iWriteDat(testCase.TempFolder, 'in.dat', data, {'Y','X','T'}, 10);
            umt = genUMTStruct(data, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});

            outDat = getDataFromROI(datFile, testCase.TempFolder, 'SpatialAggFcn', aggFcn);
            outUmt = getDataFromROI(umt, testCase.TempFolder, 'SpatialAggFcn', aggFcn);

            testCase.verifyEqual(outDat.data.main.dimNames, outUmt.data.main.dimNames);
            testCase.verifyEqual(double(outDat.data.main.value), double(outUmt.data.main.value), ...
                'RelTol', 1e-6, 'AbsTol', 1e-6);
        end

        function testForcedMultiBlockStreamingMatchesInRam(testCase)
            % The memory mock forces many frame blocks. The streamed result
            % must still equal the in-RAM one for every aggregation,
            % 'none' included.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                testCase.ProjectRoot(), 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));

            nFrames = 600;
            stream = RandStream('mt19937ar', 'Seed', 5);
            data = round(10 * randn(stream, [testCase.ImageSizeYX, nFrames], 'single'));
            datFile = iWriteDat(testCase.TempFolder, 'in.dat', data, {'Y','X','T'}, 10);
            umt = genUMTStruct(data, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T'});

            % The same sizing call the function makes: ROI columns 1:6.
            nBlocks = calculateMaxChunkSize(prod(testCase.ImageSizeYX) * nFrames * 4, 2, 0.1);
            testCase.assertGreaterThan(nBlocks, 1, ...
                'The fixture must force more than one frame block.');

            for agg = {'none','mean','median','mode','std','max','min','sum'}
                outDat = getDataFromROI(datFile, testCase.TempFolder, 'SpatialAggFcn', agg{1});
                outUmt = getDataFromROI(umt, testCase.TempFolder, 'SpatialAggFcn', agg{1});
                testCase.verifyEqual(double(outDat.data.main.value), ...
                    double(outUmt.data.main.value), 'RelTol', 1e-6, 'AbsTol', 1e-6, agg{1});
            end
        end

        function testYXFInputGivesRoiFWithTheFLabels(testCase)
            maps = iSyntheticYXT(testCase.ImageSizeYX, 2);
            umt = genUMTStruct(maps, 'kind', 'image', 'entryName', 'AzimuthMap', ...
                'dimNames', {'Y','X','F'}, 'labels', struct('F', {{'Amplitude','Phase'}}));

            out = getDataFromROI(umt, testCase.TempFolder, 'SpatialAggFcn', 'mean');

            entry = out.data.AzimuthMap;
            testCase.verifyEqual(cellstr(string(entry.dimNames)), {'ROI','F'});
            testCase.verifyEqual(size(entry.value), [2, 2]);
            testCase.verifyEqual(cellstr(string(out.labels.F(:))), {'Amplitude';'Phase'});
            testCase.verifyEqual(cellstr(string(out.labels.ROI(:))), {'ROI1';'ROI2'});

            mask1 = false(testCase.ImageSizeYX); mask1(1:3, 1:3) = true;
            frames = reshape(maps, [], 2);
            expected1 = mean(double(frames(mask1(:), :)), 1);
            testCase.verifyEqual(double(entry.value(1, :)), expected1, 'AbsTol', 1e-4);
        end

        function testRoundTripSplitEventsToROIToResponseFeatures(testCase)
            fixtureFolder = fullfile(fileparts(fileparts(fileparts( ...
                mfilename('fullpath')))), 'Analysis', 'TestingData_with_events');
            testCase.assumeTrue(isfile(fullfile(fixtureFolder, 'AcqInfos.mat')) && ...
                ~isempty(dir(fullfile(fixtureFolder, '*.dat'))), ...
                'TestingData_with_events fixture is not available locally.');

            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fx.Folder;

            copyfile(fullfile(fixtureFolder, 'AcqInfos.mat'), fullfile(saveFolder, 'AcqInfos.mat'));
            datFiles = dir(fullfile(fixtureFolder, '*.dat'));
            [dataYXT, info] = loadData(fullfile(datFiles(1).folder, datFiles(1).name));
            saveData(fullfile(saveFolder, 'inputData.dat'), single(dataYXT), ...
                'DimNames', {'Y', 'X', 'T'}, 'Info', info);

            iCreateSyntheticEvents(saveFolder, size(dataYXT, 3), info.frameRateHz);

            imageSizeYX = [size(dataYXT, 1), size(dataYXT, 2)];
            roi1 = iMakeRectROI('ROI1', imageSizeYX, ...
                round(imageSizeYX(1)*0.1):round(imageSizeYX(1)*0.3), ...
                round(imageSizeYX(2)*0.1):round(imageSizeYX(2)*0.3), [1 0 0]);
            roi2 = iMakeRectROI('ROI2', imageSizeYX, ...
                round(imageSizeYX(1)*0.6):round(imageSizeYX(1)*0.9), ...
                round(imageSizeYX(2)*0.6):round(imageSizeYX(2)*0.9), [0 1 0]);
            ROIFile = createROIFile(imageSizeYX, 'ROIs', [roi1, roi2]);
            saveROIFile(fullfile(saveFolder, 'myROI.roi'), ROIFile);

            % .dat header Phase 8c: the split is event-split image data saved
            % as .dat (its path is returned); getDataFromROI takes its labels
            % from events.mat.
            byEvFile = char(string(split_data_by_event( ...
                fullfile(saveFolder, 'inputData.dat'), saveFolder)));
            mapping = resolveDatEventMapping(loadMetaData(byEvFile), saveFolder);
            testCase.assertEqual(mapping.status, 'matched');

            roiUMT = getDataFromROI(byEvFile, saveFolder, 'SpatialAggFcn', 'mean');
            entryName = fieldnames(roiUMT.data);
            entryName = entryName{1};
            testCase.verifyEqual(cellstr(string(roiUMT.data.(entryName).dimNames)), {'ROI','T','E'});
            testCase.verifyEqual(roiUMT.eventInfo.eventID, double(mapping.eventInfo.eventID));
            testCase.verifyEqual(roiUMT.eventInfo.selected, mapping.eventInfo.selected);
            testCase.verifyEqual(roiUMT.eventInfo.baselinePeriod, mapping.eventInfo.baselinePeriod);

            % The same event-split data as an in-RAM UMT (with the mapped
            % eventInfo) gives the same result.
            byEv = single(loadData(byEvFile));
            umt = genUMTStruct(byEv, 'kind', 'image', 'entryName', 'main', ...
                'dimNames', {'Y','X','T','E'});
            umt = appendUMTEventInfo(umt, 'eventInfo', mapping.eventInfo);
            roiRam = getDataFromROI(umt, saveFolder, 'SpatialAggFcn', 'mean');
            testCase.verifyEqual(double(roiRam.data.(entryName).value), ...
                double(roiUMT.data.(entryName).value), 'RelTol', 1e-6, 'AbsTol', 1e-6);
            testCase.verifyEqual(roiRam.eventInfo, roiUMT.eventInfo);

            roiFile = fullfile(saveFolder, 'roiData.umt');
            saveData(roiFile, roiUMT);
            featOut = calculateResponseFeatures(roiFile);
            testCase.verifyEqual(char(string(featOut.kind)), 'roi');
            featEntryName = fieldnames(featOut.data);
            featEntryName = featEntryName{1};
            testCase.verifyEqual(cellstr(string(featOut.data.(featEntryName).dimNames)), ...
                {'ROI','Measure','E'});
            testCase.verifyEqual(size(featOut.data.(featEntryName).value, 1), 2);
        end

        function testRejectsUnsupportedInputs(testCase)
            folder = testCase.TempFolder;

            % Numeric arrays are not an input form.
            testCase.verifyError(@() getDataFromROI( ...
                iSyntheticYXT(testCase.ImageSizeYX, 3), folder), ...
                'Umitoolbox:getDataFromROI:UnsupportedInputType');

            % Other file types.
            matFile = fullfile(folder, 'x.mat');
            fclose(fopen(matFile, 'w'));
            testCase.verifyError(@() getDataFromROI(matFile, folder), ...
                'Umitoolbox:getDataFromROI:UnsupportedInputFile');

            % An image whose third axis is ROI is not a spatial recording.
            roiAxis = genUMTStruct(iSyntheticYXT(testCase.ImageSizeYX, 2), ...
                'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','ROI'});
            testCase.verifyError(@() getDataFromROI(roiAxis, folder), ...
                'Umitoolbox:getDataFromROI:InvalidUMTEntryDims');

            % A frame size that differs from the ROI file's, in both forms.
            bigUMT = iBuildYXTImage([8 8], 3);
            testCase.verifyError(@() getDataFromROI(bigUMT, folder), ...
                'Umitoolbox:getDataFromROI:IncompatibleSizes');
            bigDat = iWriteDat(folder, 'big.dat', iSyntheticYXT([8 8], 3), {'Y','X','T'}, 10);
            testCase.verifyError(@() getDataFromROI(bigDat, folder), ...
                'Umitoolbox:getDataFromROI:IncompatibleSizes');

            % 'none' has no {ROI,Pixel,F} layout to produce.
            yxfUMT = genUMTStruct(iSyntheticYXT(testCase.ImageSizeYX, 2), ...
                'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','F'});
            testCase.verifyError(@() getDataFromROI(yxfUMT, folder, 'SpatialAggFcn', 'none'), ...
                'Umitoolbox:getDataFromROI:NoneNotSupportedForF');
            yxfDat = iWriteDat(folder, 'maps.dat', iSyntheticYXT(testCase.ImageSizeYX, 2), ...
                {'Y','X','F'}, 10);
            testCase.verifyError(@() getDataFromROI(yxfDat, folder, 'SpatialAggFcn', 'none'), ...
                'Umitoolbox:getDataFromROI:NoneNotSupportedForF');
        end

        function testInputFormMatrix(testCase, inputForm)
            imgYXT = iSyntheticYXT(testCase.ImageSizeYX, testCase.NFrames);
            imgYXE = iSyntheticYXT(testCase.ImageSizeYX, 3);
            imgYXF = iSyntheticYXT(testCase.ImageSizeYX, 2);
            folder = testCase.TempFolder;

            switch inputForm
                case 'datYX'
                    data = iWriteDat(folder, 'raw.dat', imgYXT(:,:,1), {'Y','X'}, 10);
                    expectDims = {'ROI'};
                    expectEvt = false;

                case 'datYXT'
                    data = iWriteDat(folder, 'raw.dat', imgYXT, {'Y','X','T'}, 10);
                    expectDims = {'ROI','T'};
                    expectEvt = false;

                case 'datYXF'
                    data = iWriteDat(folder, 'raw.dat', imgYXF, {'Y','X','F'}, 10);
                    expectDims = {'ROI','F'};
                    expectEvt = false;

                case 'umtYX'
                    data = genUMTStruct(imgYXT(:,:,1), ...
                        'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X'});
                    expectDims = {'ROI'};
                    expectEvt = false;

                case 'umtYXT'
                    data = iBuildYXTImage(testCase.ImageSizeYX, testCase.NFrames);
                    expectDims = {'ROI','T'};
                    expectEvt = false;

                case 'umtYXF'
                    data = genUMTStruct(imgYXF, ...
                        'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','F'});
                    expectDims = {'ROI','F'};
                    expectEvt = false;

                case 'umtYXEnoEvt'
                    data = genUMTStruct(imgYXE, ...
                        'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','E'});
                    expectDims = {'ROI','E'};
                    expectEvt = false;

                case 'umtYXEwithEvt'
                    umt = genUMTStruct(imgYXE, ...
                        'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','E'});
                    data = appendUMTEventInfo(umt, ...
                        'eventID', [1;2;1], 'repetitionIndex', [1;1;2], ...
                        'eventName', {'A';'B';'A'}, 'eventAxisMode', 'instances', ...
                        'overwrite', true);
                    expectDims = {'ROI','E'};
                    expectEvt = true;

                case 'umtYXTEnoEvt'
                    imgYXTE = cat(4, imgYXT, imgYXT*2, imgYXT*3);
                    data = genUMTStruct(imgYXTE, ...
                        'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','T','E'});
                    expectDims = {'ROI','T','E'};
                    expectEvt = false;

                case 'umtYXTEwithEvt'
                    imgYXTE = cat(4, imgYXT, imgYXT*2, imgYXT*3);
                    umt = genUMTStruct(imgYXTE, ...
                        'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','T','E'});
                    data = appendUMTEventInfo(umt, ...
                        'eventID', [1;2;1], 'repetitionIndex', [1;1;2], ...
                        'eventName', {'A';'B';'A'}, 'eventAxisMode', 'instances', ...
                        'overwrite', true);
                    expectDims = {'ROI','T','E'};
                    expectEvt = true;
            end

            out = getDataFromROI(data, folder, 'SpatialAggFcn', 'mean');

            entryName = fieldnames(out.data);
            entryName = entryName{1};
            testCase.verifyEqual(cellstr(string(out.data.(entryName).dimNames)), expectDims);
            testCase.verifyEqual(isfield(out, 'eventInfo'), expectEvt);
        end
    end

    methods
        function root = ProjectRoot(~)
            thisFile = mfilename('fullpath');
            root = extractBefore(thisFile, [filesep 'test' filesep]);
        end
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function datFile = iWriteDat(folder, fileName, data, dimNames, frameRateHz)
%IWRITEDAT Write an image array as a headered .dat file.
datFile = fullfile(folder, fileName);
saveData(datFile, single(data), 'DimNames', dimNames, 'FrameRateHz', frameRateHz);
end

function data = iSyntheticYXT(imageSizeYX, nFrames)
%ISYNTHETICYXT Deterministic single YXT array.
data = single(reshape(1:(imageSizeYX(1)*imageSizeYX(2)*nFrames), ...
    [imageSizeYX(1), imageSizeYX(2), nFrames]));
end

function umt = iBuildYXTImage(imageSizeYX, nFrames)
%IBUILDYXTIMAGE Build a continuous {Y,X,T} image UMT with no meta.
data = iSyntheticYXT(imageSizeYX, nFrames);
umt = genUMTStruct(data, 'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','T'});
end

function umt = iBuildYXEImage(imageSizeYX, nEvents)
%IBUILDYXEIMAGE Build an event-only {Y,X,E} image UMT with no meta.
data = iSyntheticYXT(imageSizeYX, nEvents);
umt = genUMTStruct(data, 'kind', 'image', 'entryName', 'main', 'dimNames', {'Y','X','E'});
end

function roi = iMakeRectROI(name, imageSizeYX, rowRange, colRange, color)
%IMAKERECTROI Build one deterministic schema-valid rectangle ROI.

height = imageSizeYX(1);
width = imageSizeYX(2);
mask = false(height, width);
mask(rowRange, colRange) = true;

x0 = min(colRange) - 0.5; x1 = max(colRange) + 0.5;
y0 = min(rowRange) - 0.5; y1 = max(rowRange) + 0.5;
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
    'color', color, ...
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

function iCreateSyntheticEvents(saveFolder, datLen, frameRateHz)
%ICREATESYNTHETICEVENTS Create a synthetic events.mat aligned with the
%current EventsManager loading conventions (mirrors the pattern used by
%TestGenAmplitudeMaps).

onsets = round(linspace(max(3, 0.15*datLen), max(4, 0.7*datLen), 3));
dur = max(3, round(0.05*datLen));
offsets = min(onsets + dur - 1, datLen);

timestamps = reshape([onsets; offsets], [], 1);
timestamps = single((timestamps - 1) ./ frameRateHz);

state = logical(repmat([1; 0], numel(onsets), 1));
eventID = uint16(ones(numel(onsets)*2, 1));
eventNameList = {'CondA'};
selectedEvents = true(size(eventID));
baselinePeriod = single(max(2, round(0.03*datLen)) / frameRateHz);

save(fullfile(saveFolder, 'events.mat'), ...
    'timestamps', 'state', 'eventID', 'eventNameList', ...
    'baselinePeriod', 'selectedEvents', '-mat');
end
