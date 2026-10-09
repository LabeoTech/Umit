classdef TestImportFromTif < matlab.unittest.TestCase
    % TESTIMPORTFROMTIF Unit tests for importFromTif using synthetic TIFF data.
    %
    % This suite generates synthetic 3-D arrays, writes them as multi-page
    % TIFF files, creates a matching info.json file, runs importFromTif,
    % and validates the generated .dat files and AcqInfos.mat content.

    properties (Access = private)
        TempFolder char = ''
        RigRoot char = ''
        RigActivationFixture = struct()
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            fx = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;

            suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
            rigStore = UMITRigStore.create(struct('rigID', ['TifImport_' suffix(1:8)]));
            testCase.RigRoot = rigStore.RigRoot;
            testCase.RigActivationFixture = ...
                activateRigTemporarily(rigStore.getRigInfo().uuid);
        end
    end

    methods (TestMethodTeardown)
        function removeRigFixture(testCase)
            deactivateRigTemporarily(testCase.RigActivationFixture);
            if ~isempty(testCase.RigRoot) && isfolder(testCase.RigRoot)
                rmdir(testCase.RigRoot, 's');
            end
        end
    end

    methods (Test)
        function testPipelineManagerExecutesImporter(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_pm');
            saveFolder = fullfile(testCase.TempFolder, 'Save_pm');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:5*6*4), [5, 6, 4]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), data);
            tiffiles = struct('filename', 'img_red.tif', ...
                'FrameRateHz', 30.0, 'ExposureMsec', 0.1, ...
                'IlluminationColor', 'Red');
            globalMeta = struct('DateTime', '20260411_120000', ...
                'Camera_Model', 'CS2100M');
            localWriteInfoJson(rawFolder, globalMeta, tiffiles);

            pm = buildPMForScenario(saveFolder, 'importFromTif', 'auto', ...
                'RawFolder', rawFolder);
            pm.executePipeline('PrintSummary', false);

            testCase.verifyTrue(isfile(fullfile(saveFolder, 'red.dat')));
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
        end

        function testPipelineInfoContract(testCase)
            info = importFromTif('pipelineInfo');

            testCase.verifyTrue( ...
                isstruct(info), ...
                'importFromTif(''pipelineInfo'') must return a struct.');
            testCase.verifyEqual(info.freshSaveFolderRole, ...
                'acquisition-initializer');

            testCase.verifyTrue( ...
                isfield(info, 'outputs') && ~isempty(info.outputs), ...
                'pipelineInfo must declare at least one output.');

            outIdx = find(strcmp({info.outputs.name}, 'outFile'), 1, 'first');
            testCase.assertNotEmpty( ...
                outIdx, ...
                'pipelineInfo must declare the output "outFile".');

            testCase.verifyEqual(info.outputs(outIdx).outputMode, 'file');
            testCase.verifyTrue(info.outputs(outIdx).isData);

            paramNames = string.empty(0,1);

            if isfield(info, 'parameters') && ~isempty(info.parameters)
                tmpNames = {info.parameters.name};
                paramNames = [paramNames; string(tmpNames(:))];
            end

            if isfield(info, 'arguments') && ~isempty(info.arguments)
                argKinds = {info.arguments.kind};
                argMask = strcmpi(argKinds, 'parameter');
                if any(argMask)
                    tmpNames = {info.arguments(argMask).name};
                    paramNames = [paramNames; string(tmpNames(:))];
                end
            end

            paramNames = unique(paramNames, 'stable');

            testCase.verifyTrue( ...
                any(strcmpi(paramNames, "BinningSpatial")), ...
                'pipelineInfo must expose "BinningSpatial" as a parameter.');

            testCase.verifyTrue( ...
                any(strcmpi(paramNames, "BinningTemp")), ...
                'pipelineInfo must expose "BinningTemp" as a parameter.');
        end

        function testSingleChannelImportNoBinning(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_single');
            saveFolder = fullfile(testCase.TempFolder, 'Save_single');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:5*6*4), [5, 6, 4]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), data);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 30.0, ...
                'ExposureMsec', 0.1, ...
                'IlluminationColor', 'Red');

            globalMeta = struct( ...
                'DateTime', '20260411_120000', ...
                'Camera_Model', 'CS2100M');

            localWriteInfoJson(rawFolder, globalMeta, tiffiles);

            [outFile, rigResolution] = importFromTif(rawFolder, saveFolder);

            testCase.verifyEqual(outFile, {'red.dat'});
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'red.dat')));
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));

            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            testCase.verifyEqual(AcqInfoStream.rigUUID, rigResolution.rigUUID);
            testCase.verifyEqual(AcqInfoStream.rigID, rigResolution.rigID);
            testCase.verifyFalse(rigResolution.wasCreated);

            testCase.verifyEqual(AcqInfoStream.Width, size(data,2));
            testCase.verifyEqual(AcqInfoStream.Height, size(data,1));
            localVerifyRawScope(testCase, AcqInfoStream);
            testCase.verifyEqual(AcqInfoStream.FrameRateHz, 30.0);
            testCase.verifyEqual(AcqInfoStream.ExposureMsec, 0.1);
            testCase.verifyEqual(char(string(AcqInfoStream.DateTime)), '20260411_120000');
            testCase.verifyEqual(char(string(AcqInfoStream.Camera_Model)), 'CS2100M');

            imported = localReadDatFile(fullfile(saveFolder, 'red.dat'));

            testCase.verifyEqual(imported, single(data));

            % .dat header Phase 4d: the channel file carries its own rate,
            % exposure, and tag.
            hdr = readDatHeader(fullfile(saveFolder, 'red.dat'));
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.dimSizes, size(data));
            testCase.verifyEqual(hdr.frameRateHz, 30);
            testCase.verifyEqual(hdr.exposureMsec, double(single(0.1)));
            testCase.verifyEqual(hdr.channelName, 'red');
        end

        function testMultiChannelImportCreatesMultipleDatFiles(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_multi');
            saveFolder = fullfile(testCase.TempFolder, 'Save_multi');
            mkdir(rawFolder);
            mkdir(saveFolder);

            redData = reshape(uint16(1:4*5*3), [4, 5, 3]);
            fluoData = reshape(uint16(1000 + (1:4*5*3)), [4, 5, 3]);

            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), redData);
            localWriteMultipageTif(fullfile(rawFolder, 'img_fluo.tif'), fluoData);

            tiffiles(1) = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');
            tiffiles(2) = struct( ...
                'filename', 'img_fluo.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Fluo');

            globalMeta = struct( ...
                'DateTime', '20260411_120100', ...
                'Camera_Model', 'CS2100M');

            localWriteInfoJson(rawFolder, globalMeta, tiffiles);

            outFile = importFromTif(rawFolder, saveFolder);

            testCase.verifyEqual(outFile, {'red.dat', 'fluo.dat'});
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'red.dat')));
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'fluo.dat')));
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));

            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            testCase.verifyEqual(AcqInfoStream.Width, size(redData,2));
            testCase.verifyEqual(AcqInfoStream.Height, size(redData,1));
            localVerifyRawScope(testCase, AcqInfoStream);

            redImported = localReadDatFile(fullfile(saveFolder, 'red.dat'));

            fluoImported = localReadDatFile(fullfile(saveFolder, 'fluo.dat'));

            testCase.verifyEqual(redImported, single(redData));
            testCase.verifyEqual(fluoImported, single(fluoData));
        end

        function testTemporalBinning(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_tempBin');
            saveFolder = fullfile(testCase.TempFolder, 'Save_tempBin');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:6*4*6), [6, 4, 6]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_green.tif'), data);

            tiffiles = struct( ...
                'filename', 'img_green.tif', ...
                'FrameRateHz', 24.0, ...
                'ExposureMsec', 0.5, ...
                'IlluminationColor', 'Green');

            globalMeta = struct('DateTime', '20260411_120200');
            localWriteInfoJson(rawFolder, globalMeta, tiffiles);

            outFile = importFromTif(rawFolder, saveFolder, 'BinningTemp', 2);

            testCase.verifyEqual(outFile, {'green.dat'});

            expected = imresize3(single(data), [size(data,1), size(data,2), size(data,3)/2], 'linear');
            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            imported = localReadDatFile(fullfile(saveFolder, 'green.dat'));

            % .dat header Phase 4d: the header rate is the binned rate.
            hdr = readDatHeader(fullfile(saveFolder, 'green.dat'));
            testCase.verifyEqual(hdr.frameRateHz, 12);
            testCase.verifyEqual(hdr.dimSizes, size(expected));
            % .dat header Phase 7b: AcqInfos.mat keeps the raw rate.
            testCase.verifyEqual(AcqInfoStream.FrameRateHz, 24.0);
            testCase.verifyEqual(AcqInfoStream.BinningTemp, 2);
            localVerifyRawScope(testCase, AcqInfoStream);
            testCase.verifyEqual(imported, expected, 'AbsTol', 1e-5);
        end

        function testSpatialBinning(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_spatBin');
            saveFolder = fullfile(testCase.TempFolder, 'Save_spatBin');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:8*6*4), [8, 6, 4]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_yellow.tif'), data);

            tiffiles = struct( ...
                'filename', 'img_yellow.tif', ...
                'FrameRateHz', 15.0, ...
                'ExposureMsec', 1.0, ...
                'IlluminationColor', 'Yellow');

            globalMeta = struct('DateTime', '20260411_120300');
            localWriteInfoJson(rawFolder, globalMeta, tiffiles);

            outFile = importFromTif(rawFolder, saveFolder, 'BinningSpatial', 2);

            testCase.verifyEqual(outFile, {'yellow.dat'});

            expected = imresize(single(data), 1/2);
            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            imported = localReadDatFile(fullfile(saveFolder, 'yellow.dat'));

            % .dat header Phase 7b: raw size in AcqInfos.mat, binned in the header.
            testCase.verifyEqual(AcqInfoStream.Height, size(data,1));
            testCase.verifyEqual(AcqInfoStream.Width, size(data,2));
            testCase.verifyEqual(AcqInfoStream.BinningSpatial, 2);
            testCase.verifyEqual(readDatHeader(fullfile(saveFolder, 'yellow.dat')).dimSizes, ...
                size(expected));
            localVerifyRawScope(testCase, AcqInfoStream);
            testCase.verifyEqual(imported, expected, 'AbsTol', 1e-5);
        end

        function testFileSequenceConcatenation(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_seq');
            saveFolder = fullfile(testCase.TempFolder, 'Save_seq');
            mkdir(rawFolder);
            mkdir(saveFolder);

            dataA = reshape(uint16(1:5*4*2), [5, 4, 2]);
            dataB = reshape(uint16(100 + (1:5*4*3)), [5, 4, 3]);

            localWriteMultipageTif(fullfile(rawFolder, 'img_red_001.tif'), dataA);
            localWriteMultipageTif(fullfile(rawFolder, 'img_red_002.tif'), dataB);

            tiffiles = struct( ...
                'filename', 'img_red_001.tif', ...
                'FrameRateHz', 12.0, ...
                'ExposureMsec', 0.3, ...
                'IlluminationColor', 'Red');

            globalMeta = struct('DateTime', '20260411_120400');
            localWriteInfoJson(rawFolder, globalMeta, tiffiles);

            outFile = importFromTif(rawFolder, saveFolder);

            testCase.verifyEqual(outFile, {'red.dat'});

            expected = cat(3, single(dataA), single(dataB));
            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            imported = localReadDatFile(fullfile(saveFolder, 'red.dat'));

            localVerifyRawScope(testCase, AcqInfoStream);
            testCase.verifyEqual(imported, expected);
        end

        function testUnnumberedBaseFileIsIncludedInSequence(testCase)
            % "img.tif" alongside "img_1.tif"/"img_2.tif" matches no numeric
            % suffix, so it used to be dropped from the sequence and its
            % frames were silently never imported (P1-11).
            rawFolder = fullfile(testCase.TempFolder, 'Raw_baseSeq');
            saveFolder = fullfile(testCase.TempFolder, 'Save_baseSeq');
            mkdir(rawFolder);
            mkdir(saveFolder);

            dataBase = reshape(uint16(1:5*4*2), [5, 4, 2]);
            data1 = reshape(uint16(100 + (1:5*4*3)), [5, 4, 3]);
            data2 = reshape(uint16(200 + (1:5*4*2)), [5, 4, 2]);

            localWriteMultipageTif(fullfile(rawFolder, 'img.tif'), dataBase);
            localWriteMultipageTif(fullfile(rawFolder, 'img_1.tif'), data1);
            localWriteMultipageTif(fullfile(rawFolder, 'img_2.tif'), data2);

            tiffiles = struct( ...
                'filename', 'img.tif', ...
                'FrameRateHz', 12.0, ...
                'ExposureMsec', 0.3, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_120400'), tiffiles);

            outFile = importFromTif(rawFolder, saveFolder);
            testCase.verifyEqual(outFile, {'red.dat'});

            % The base file supplies the first frames of the sequence.
            expected = cat(3, single(dataBase), single(data1), single(data2));
            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            imported = localReadDatFile(fullfile(saveFolder, 'red.dat'));

            localVerifyRawScope(testCase, AcqInfoStream);
            testCase.verifyEqual(imported, expected);
        end

        function testExistingAcqInfosAllowsMatchingSecondImport(testCase)
            saveFolder = fullfile(testCase.TempFolder, 'Save_existingMatch');
            mkdir(saveFolder);

            % First import
            rawFolder1 = fullfile(testCase.TempFolder, 'Raw_existingMatch_1');
            mkdir(rawFolder1);

            data1 = reshape(uint16(1:4*5*3), [4, 5, 3]);
            localWriteMultipageTif(fullfile(rawFolder1, 'img_red.tif'), data1);

            tiffiles1 = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder1, struct('DateTime', '20260411_120500'), tiffiles1);
            outFile1 = importFromTif(rawFolder1, saveFolder);

            testCase.verifyEqual(outFile1, {'red.dat'});
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));

            % Second import with matching dimensions/timing
            rawFolder2 = fullfile(testCase.TempFolder, 'Raw_existingMatch_2');
            mkdir(rawFolder2);

            data2 = reshape(uint16(1000 + (1:4*5*3)), [4, 5, 3]);
            localWriteMultipageTif(fullfile(rawFolder2, 'img_green.tif'), data2);

            tiffiles2 = struct( ...
                'filename', 'img_green.tif', ...
                'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, ...
                'IlluminationColor', 'Green');

            localWriteInfoJson(rawFolder2, struct('DateTime', '20260411_120501'), tiffiles2);
            outFile2 = importFromTif(rawFolder2, saveFolder);

            testCase.verifyEqual(outFile2, {'green.dat'});
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'red.dat')));
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'green.dat')));

            AcqInfoStream = localLoadAcqInfo(fullfile(saveFolder, 'AcqInfos.mat'));
            testCase.verifyEqual(AcqInfoStream.Width, 5);
            testCase.verifyEqual(AcqInfoStream.Height, 4);
            localVerifyRawScope(testCase, AcqInfoStream);
            testCase.verifyEqual(AcqInfoStream.FrameRateHz, 18.0);
        end

        function testExistingAcqInfosMismatchErrors(testCase)
            saveFolder = fullfile(testCase.TempFolder, 'Save_existingMismatch');
            mkdir(saveFolder);

            % First import creates AcqInfos.mat
            rawFolder1 = fullfile(testCase.TempFolder, 'Raw_existingMismatch_1');
            mkdir(rawFolder1);

            data1 = reshape(uint16(1:4*5*3), [4, 5, 3]);
            localWriteMultipageTif(fullfile(rawFolder1, 'img_red.tif'), data1);

            tiffiles1 = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder1, struct('DateTime', '20260411_120600'), tiffiles1);
            importFromTif(rawFolder1, saveFolder);

            % Second import has mismatched dimensions
            rawFolder2 = fullfile(testCase.TempFolder, 'Raw_existingMismatch_2');
            mkdir(rawFolder2);

            data2 = reshape(uint16(1000 + (1:6*5*3)), [6, 5, 3]);
            localWriteMultipageTif(fullfile(rawFolder2, 'img_green.tif'), data2);

            tiffiles2 = struct( ...
                'filename', 'img_green.tif', ...
                'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, ...
                'IlluminationColor', 'Green');

            localWriteInfoJson(rawFolder2, struct('DateTime', '20260411_120601'), tiffiles2);

            didThrow = false;
            try
                importFromTif(rawFolder2, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'must have identical raw Width, Height, FrameRateHz, and ExposureMsec'), ...
                    'Mismatch import should fail with the shared-dimension validation message.');
            end

            testCase.verifyTrue(didThrow, ...
                'A mismatched second import should raise an error.');
        end

        function testExistingAcqInfosBinningMismatchErrors(testCase)
            % .dat header Phase 7b: raw values match, but a different
            % binning would give .dat files of another size.
            saveFolder = fullfile(testCase.TempFolder, 'Save_binMismatch');
            mkdir(saveFolder);

            rawFolder1 = fullfile(testCase.TempFolder, 'Raw_binMismatch_1');
            mkdir(rawFolder1);
            localWriteMultipageTif(fullfile(rawFolder1, 'img_red.tif'), ...
                reshape(uint16(1:4*4*2), [4, 4, 2]));
            localWriteInfoJson(rawFolder1, struct('DateTime', '20260411_120700'), ...
                struct('filename', 'img_red.tif', 'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, 'IlluminationColor', 'Red'));
            importFromTif(rawFolder1, saveFolder);

            rawFolder2 = fullfile(testCase.TempFolder, 'Raw_binMismatch_2');
            mkdir(rawFolder2);
            localWriteMultipageTif(fullfile(rawFolder2, 'img_green.tif'), ...
                reshape(uint16(1:4*4*2), [4, 4, 2]));
            localWriteInfoJson(rawFolder2, struct('DateTime', '20260411_120701'), ...
                struct('filename', 'img_green.tif', 'FrameRateHz', 18.0, ...
                'ExposureMsec', 0.4, 'IlluminationColor', 'Green'));

            didThrow = false;
            try
                importFromTif(rawFolder2, saveFolder, 'BinningSpatial', 2);
            catch ME
                didThrow = true;
                testCase.verifyTrue(contains(ME.message, 'same BinningSpatial and BinningTemp'), ...
                    'A binning mismatch should fail with the shared-dimension message.');
            end
            testCase.verifyTrue(didThrow, ...
                'A second import with a different binning should raise an error.');
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'green.dat')));
        end

        function testFrameCountMismatchWithinImportErrors(testCase)
            % .dat header Phase 7b: the frame count is compared in memory.
            rawFolder = fullfile(testCase.TempFolder, 'Raw_lenMismatch');
            saveFolder = fullfile(testCase.TempFolder, 'Save_lenMismatch');
            mkdir(rawFolder);
            mkdir(saveFolder);

            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), ...
                reshape(uint16(1:4*5*3), [4, 5, 3]));
            localWriteMultipageTif(fullfile(rawFolder, 'img_fluo.tif'), ...
                reshape(uint16(1:4*5*2), [4, 5, 2]));
            tiffiles(1) = struct('filename', 'img_red.tif', 'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, 'IlluminationColor', 'Red');
            tiffiles(2) = struct('filename', 'img_fluo.tif', 'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, 'IlluminationColor', 'Fluo');
            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_120800'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue(contains(ME.message, 'same number of frames'), ...
                    'A frame-count mismatch should fail with an informative message.');
            end
            testCase.verifyTrue(didThrow, ...
                'TIFF sequences with different frame counts should raise an error.');
        end

        function testMissingInfoJsonErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingInfo');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingInfo');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:4*5*2), [4, 5, 2]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), data);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'info.json'), ...
                    'Missing info.json should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif without info.json should raise an error.');
        end

        function testMissingTiffilesFieldErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingTiffiles');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingTiffiles');
            mkdir(rawFolder);
            mkdir(saveFolder);

            meta = struct('DateTime', '20260411_130000');
            infoPath = fullfile(rawFolder, 'info.json');

            fid = fopen(infoPath, 'w');
            testCase.assertNotEqual(fid, -1, 'Failed to create info.json.');
            c = onCleanup(@() fclose(fid));
            fwrite(fid, jsonencode(meta), 'char');
            clear c

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'Tiffiles'), ...
                    'Missing Tiffiles field should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif without Tiffiles should raise an error.');
        end

        function testEmptyTiffilesErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_emptyTiffiles');
            saveFolder = fullfile(testCase.TempFolder, 'Save_emptyTiffiles');
            mkdir(rawFolder);
            mkdir(saveFolder);

            meta = struct( ...
                'DateTime', '20260411_130100', ...
                'Tiffiles', {{}});
            infoPath = fullfile(rawFolder, 'info.json');

            fid = fopen(infoPath, 'w');
            testCase.assertNotEqual(fid, -1, 'Failed to create info.json.');
            c = onCleanup(@() fclose(fid));
            fwrite(fid, jsonencode(meta), 'char');
            clear c

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'Missing information about TIFF files'), ...
                    'Empty Tiffiles should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with empty Tiffiles should raise an error.');
        end

        function testMissingReferencedTifErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingTif');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingTif');
            mkdir(rawFolder);
            mkdir(saveFolder);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130200'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'No TIFF file found'), ...
                    'Missing referenced TIFF should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with a missing TIFF file should raise an error.');
        end

        function testMissingFilenameFieldErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingFilename');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingFilename');
            mkdir(rawFolder);
            mkdir(saveFolder);

            tiffiles = struct( ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130300'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'Missing field "filename"'), ...
                    'Missing filename should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with missing filename should raise an error.');
        end

        function testMissingFrameRateFieldErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingFrameRate');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingFrameRate');
            mkdir(rawFolder);
            mkdir(saveFolder);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130400'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'Missing field "FrameRateHz"'), ...
                    'Missing FrameRateHz should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with missing FrameRateHz should raise an error.');
        end

        function testMissingExposureFieldErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingExposure');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingExposure');
            mkdir(rawFolder);
            mkdir(saveFolder);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130500'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'Missing field "ExposureMsec"'), ...
                    'Missing ExposureMsec should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with missing ExposureMsec should raise an error.');
        end
        function testMissingIlluminationColorFieldErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingIllumination');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingIllumination');
            mkdir(rawFolder);
            mkdir(saveFolder);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2);

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130600'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'Missing field "IlluminationColor"'), ...
                    'Missing IlluminationColor should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with missing IlluminationColor should raise an error.');
        end

        function testTemporalBinningNonDivisibleFrameCountErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_badTempBin');
            saveFolder = fullfile(testCase.TempFolder, 'Save_badTempBin');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:4*5*5), [4, 5, 5]); % 5 frames, not divisible by 2
            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), data);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130700'), tiffiles);

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder, 'BinningTemp', 2);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'divisible by BinningTemp'), ...
                    'Non-divisible temporal binning should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with non-divisible BinningTemp should raise an error.');
        end

        function testMalformedExistingAcqInfosErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_badAcqInfo');
            saveFolder = fullfile(testCase.TempFolder, 'Save_badAcqInfo');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:4*5*3), [4, 5, 3]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), data);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130800'), tiffiles);

            junk = 42;
            save(fullfile(saveFolder, 'AcqInfos.mat'), 'junk');

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'does not contain a valid acquisition-info structure'), ...
                    'Malformed AcqInfos.mat should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with malformed AcqInfos.mat should raise an error.');
        end

        function testExistingAcqInfosMissingRequiredFieldsErrors(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_missingAcqFields');
            saveFolder = fullfile(testCase.TempFolder, 'Save_missingAcqFields');
            mkdir(rawFolder);
            mkdir(saveFolder);

            data = reshape(uint16(1:4*5*3), [4, 5, 3]);
            localWriteMultipageTif(fullfile(rawFolder, 'img_red.tif'), data);

            tiffiles = struct( ...
                'filename', 'img_red.tif', ...
                'FrameRateHz', 20.0, ...
                'ExposureMsec', 0.2, ...
                'IlluminationColor', 'Red');

            localWriteInfoJson(rawFolder, struct('DateTime', '20260411_130900'), tiffiles);

            AcqInfoStream = struct('Width', 5);
            save(fullfile(saveFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            didThrow = false;
            try
                importFromTif(rawFolder, saveFolder);
            catch ME
                didThrow = true;
                testCase.verifyTrue( ...
                    contains(ME.message, 'does not contain the required fields'), ...
                    'AcqInfos.mat missing required fields should raise an informative error.');
            end

            testCase.verifyTrue(didThrow, ...
                'Calling importFromTif with incomplete AcqInfos.mat should raise an error.');
        end
        function testInvalidBinningSpatialThrows(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_badSpatialParam');
            saveFolder = fullfile(testCase.TempFolder, 'Save_badSpatialParam');
            mkdir(rawFolder);
            mkdir(saveFolder);

            testCase.verifyError( ...
                @() importFromTif(rawFolder, saveFolder, 'BinningSpatial', 3), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end
        function testInvalidBinningTempThrows(testCase)
            rawFolder = fullfile(testCase.TempFolder, 'Raw_badTempParam');
            saveFolder = fullfile(testCase.TempFolder, 'Save_badTempParam');
            mkdir(rawFolder);
            mkdir(saveFolder);

            testCase.verifyError( ...
                @() importFromTif(rawFolder, saveFolder, 'BinningTemp', 9), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end
    end
end

% =========================================================================
% Local helpers
% =========================================================================
function localWriteMultipageTif(filePath, data)
%LOCALWRITEMULTIPAGETIF Write a 3-D array as a multi-page TIFF.

validateattributes(data, {'numeric'}, {'3d', 'nonempty', 'real', 'nonsparse'}, ...
    'localWriteMultipageTif', 'data');

data = data(:,:,1:end);

for ii = 1:size(data,3)
    if ii == 1
        imwrite(data(:,:,ii), filePath, 'tif', 'Compression', 'none');
    else
        imwrite(data(:,:,ii), filePath, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end
end
end

function localWriteInfoJson(rawFolder, globalMeta, tiffiles)
%LOCALWRITEINFOJSON Create an info.json file that follows the documented template.

meta = globalMeta;
meta.Tiffiles = tiffiles;

infoPath = fullfile(rawFolder, 'info.json');
fid = fopen(infoPath, 'w');
assert(fid ~= -1, 'Failed to create info.json in "%s".', rawFolder);
c = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(meta), 'char');
clear c
end

function AcqInfoStream = localLoadAcqInfo(acqInfoPath)
%LOCALLOADACQINFO Load AcqInfoStream from AcqInfos.mat.

tmp = load(acqInfoPath);

if isfield(tmp, 'AcqInfoStream')
    AcqInfoStream = tmp.AcqInfoStream;
elseif isfield(tmp, 'AcqInfos')
    AcqInfoStream = tmp.AcqInfos;
else
    fn = fieldnames(tmp);
    assert(~isempty(fn), 'The MAT file "%s" does not contain any variables.', acqInfoPath);
    AcqInfoStream = tmp.(fn{1});
end
end

function data = localReadDatFile(datPath)
%LOCALREADDATFILE Read an imported .dat file back into YxXxT.
%
% Imported files are headered since .dat header Phase 4d, and their header
% is the only description of the imported size (.dat header Phase 7b).
% Callers compare the whole array, size included.

data = loadData(datPath);
end

function localVerifyRawScope(testCase, AcqInfoStream)
%LOCALVERIFYRAWSCOPE AcqInfos.mat holds no imported-data facts (Phase 7b).

testCase.verifyFalse(isfield(AcqInfoStream, 'Length'), ...
    'AcqInfos.mat must not store the imported Length.');
testCase.verifyFalse(isfield(AcqInfoStream, 'Datatype'), ...
    'AcqInfos.mat must not store the imported Datatype.');
testCase.verifyTrue(isfield(AcqInfoStream, 'BinningSpatial'));
testCase.verifyTrue(isfield(AcqInfoStream, 'BinningTemp'));
testCase.verifyTrue(isfield(AcqInfoStream, 'ImportedChannels'));
end
