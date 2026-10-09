classdef TestImagesClassification < matlab.unittest.TestCase
    %TESTIMAGESCLASSIFICATION Unit tests for ImagesClassification.
    %
    %   This suite verifies file generation, backup handling, the shared
    %   AcqInfos.mat metadata organization, and repeated-illumination
    %   channel grouping using a synthetic info.txt mutation of the fixture
    %   dataset.

    properties
        ProjectRoot
        FixtureFolder
        TempSaveFolder
        TempRawFolder
    end

    methods (TestMethodSetup)
        function setup(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            if isempty(projectRoot)
                projectRoot = fileparts(fileparts(testFolder));
            end
            testCase.ProjectRoot = char(projectRoot);
            addpath(genpath(testCase.ProjectRoot));

            testCase.FixtureFolder = fullfile(testCase.ProjectRoot, 'test', 'Analysis', 'TestingData_with_events');
            testCase.TempSaveFolder = fullfile(tempdir, ['TestImagesClassification_' char(java.util.UUID.randomUUID)]);
            testCase.TempRawFolder = fullfile(tempdir, ['TestImagesClassificationRaw_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempSaveFolder);
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if isfolder(testCase.TempSaveFolder)
                rmdir(testCase.TempSaveFolder, 's');
            end
            if isfolder(testCase.TempRawFolder)
                rmdir(testCase.TempRawFolder, 's');
            end
        end
    end

    methods (Test)
        function testGeneratesManifestAndAcqInfos(testCase)
            %TESTGENERATESMANIFESTANDACQINFOS Generate .dat files and AcqInfos.mat.

            rigsBefore = UMITRigStore.listRigs();
            outFile = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'ERASE');
            rigsAfter = UMITRigStore.listRigs();

            testCase.verifyTrue(iscell(outFile));
            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
            testCase.verifyEqual(rigsAfter.RigUUID, rigsBefore.RigUUID);
            info = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyFalse(isfield(info.AcqInfoStream, 'rigUUID'));
            testCase.verifyFalse(isfield(info.AcqInfoStream, 'rigID'));

            datFiles = outFile(endsWith(outFile, '.dat'));
            testCase.verifyEqual(numel(datFiles), numel(outFile), ...
                'ImagesClassification should return imported .dat outputs only.');
            testCase.verifyNotEmpty(datFiles);
            for iFile = 1:numel(datFiles)
                testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, datFiles{iFile})));
            end
        end

        function testNoOutputArgumentStillWritesFiles(testCase)
            %TESTNOOUTPUTARGUMENTSTILLWRITESFILES Write output files without manifest request.

            ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'ERASE');

            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
            listing = dir(fullfile(testCase.TempSaveFolder, '*.dat'));
            testCase.verifyNotEmpty(listing);
        end

        function testBackupEraseRemovesManagedFiles(testCase)
            %TESTBACKUPERASEREMOVESMANAGEDFILES Remove existing files before import.

            dummyFile = fullfile(testCase.TempSaveFolder, 'old_result.dat');
            writeTestDat(dummyFile, single(1), 1);

            testCase.verifyTrue(isfile(dummyFile));

            ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'ERASE');

            testCase.verifyFalse(isfile(dummyFile));
            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
        end

        function testBackupGenBackupCreatesZip(testCase)
            %TESTBACKUPGENBACKUPCREATESZIP Back up existing files before import.

            dummyFile = fullfile(testCase.TempSaveFolder, 'old_result.dat');
            writeTestDat(dummyFile, single(1), 1);

            ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'GENBACKUP');

            zipList = dir(fullfile(testCase.TempSaveFolder, 'bkp_*.zip'));
            testCase.verifyNotEmpty(zipList);
            testCase.verifyFalse(isfile(dummyFile));
            testCase.verifyTrue(isfile(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat')));
        end

        function testAcqInfosContainsImportedChannels(testCase)
            %TESTACQINFOSCONTAINSIMPORTEDCHANNELS Validate imported-channel metadata.

            outFile = ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'ERASE');

            S = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyTrue(isfield(S, 'AcqInfoStream'));
            AcqInfoStream = S.AcqInfoStream;

            testCase.verifyTrue(isfield(AcqInfoStream, 'ImportedChannels'));
            importedChannels = AcqInfoStream.ImportedChannels;
            testCase.verifyNotEmpty(importedChannels);

            requiredFields = {'DatFile', 'Length', 'FrameRateHz', 'ExposureMsec', 'CamIdx'};
            for iField = 1:numel(requiredFields)
                testCase.verifyTrue(isfield(importedChannels, requiredFields{iField}), ...
                    ['Missing ImportedChannels field: ' requiredFields{iField}]);
            end

            testCase.verifyFalse(isfield(importedChannels, 'RepeatCount'));
            testCase.verifyFalse(isfield(importedChannels, 'SequenceIdx'));

            datFiles = outFile(endsWith(outFile, '.dat'));
            testCase.verifyEqual(sort({importedChannels.DatFile}), sort(datFiles));

            for iChan = 1:numel(importedChannels)
                datPath = fullfile(testCase.TempSaveFolder, importedChannels(iChan).DatFile);
                testCase.verifyTrue(isfile(datPath));

                actualLength = testCase.getDatLengthFromFile(datPath, AcqInfoStream.Height, AcqInfoStream.Width);
                testCase.verifyEqual(actualLength, double(importedChannels(iChan).Length));
                testCase.verifyEqual(double(importedChannels(iChan).CamIdx), 1);

                % .dat header Phase 4d: each channel file carries its own
                % rate, exposure, and tag.
                testCase.assertTrue(isDatWithHeader(datPath));
                hdr = readDatHeader(datPath);
                [~, tag] = fileparts(importedChannels(iChan).DatFile);
                testCase.verifyTrue(hdr.writeComplete);
                testCase.verifyEqual(hdr.dataClass, 'single');
                testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T'});
                testCase.verifyEqual(hdr.dimSizes, [double(AcqInfoStream.Height), ...
                    double(AcqInfoStream.Width), double(importedChannels(iChan).Length)]);
                testCase.verifyEqual(hdr.frameRateHz, ...
                    double(single(importedChannels(iChan).FrameRateHz)));
                testCase.verifyEqual(hdr.exposureMsec, ...
                    double(single(importedChannels(iChan).ExposureMsec)));
                testCase.verifyEqual(hdr.channelName, tag);
            end
        end

        function testPixelValuesMatchHeaderlessReference(testCase)
            %TESTPIXELVALUESMATCHHEADERLESSREFERENCE Values unchanged by the header.
            %
            %   .dat header Phase 4d: fixtures/imagesClassificationReference.mat
            %   holds, per binning, each channel's size and the SHA-256 of
            %   its values as written by the headerless ImagesClassification
            %   (commit 7fe2f77) from the TestingData_with_events fixture.

            S = load(fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
                'imagesClassificationReference.mat'), 'ref');
            for iRef = 1:numel(S.ref)
                ref = S.ref(iRef);
                outFolder = fullfile(testCase.TempSaveFolder, sprintf('bin%d', iRef));
                mkdir(outFolder);
                ImagesClassification(testCase.FixtureFolder, outFolder, ...
                    ref.binning(1), ref.binning(2), 0, 'backupOpts', 'ERASE');

                for iChan = 1:numel(ref.channels)
                    datPath = fullfile(outFolder, ref.channels(iChan).DatFile);
                    values = loadData(datPath);
                    testCase.verifyEqual(size(values), ref.channels(iChan).dimSizes, ...
                        sprintf('%s size (binning %s)', ref.channels(iChan).DatFile, mat2str(ref.binning)));
                    testCase.verifyEqual(iSha256(values), ref.channels(iChan).sha256, ...
                        sprintf('%s values (binning %s)', ref.channels(iChan).DatFile, mat2str(ref.binning)));
                end
            end
        end

        function testLoadMetaDataReturnsFileFacingFields(testCase)
            %TESTLOADMETADATARETURNSFILEFACINGFIELDS Keep metadata focused on the .dat file.

            ImagesClassification( ...
                testCase.FixtureFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'ERASE');

            listing = dir(fullfile(testCase.TempSaveFolder, '*.dat'));
            testCase.assumeNotEmpty(listing);

            Info = loadMetaData(fullfile(listing(1).folder, listing(1).name));

            % .dat header Phase 8a: channel files are headered, so Info is
            % exactly the .dat schema. The former compatibility names
            % (Height, Length, FrameRateHz, MetadataSource, ...) are gone.
            expectedFields = {'filePath', 'format', 'dataOffset', 'dataClass', ...
                'dimNames', 'dimSizes', 'frameRateHz', 'exposureMsec', ...
                'channelName', 'writeComplete'};
            testCase.verifyEqual(Info.format, 'header');
            testCase.verifyEqual(sort(fieldnames(Info)), sort(expectedFields(:)), ...
                'Info of a headered file must hold exactly the .dat schema fields.');

            rejectedFields = {'fileName', 'ImportedChannels', 'Illumination1', ...
                'AISampleRate', 'Acquisition_Duration', 'OriginalLength', ...
                'Tag', 'Color', 'MetadataSource', 'Length', 'FrameRateHz', ...
                'Height', 'Width', 'Datatype', 'datLength', 'datSize', ...
                'dim_names', 'Freq', 'ExposureMsec'};

            for iField = 1:numel(rejectedFields)
                testCase.verifyFalse(isfield(Info, rejectedFields{iField}), ...
                    ['Unexpected Info field: ' rejectedFields{iField}]);
            end
        end

        function testRepeatingIlluminationCreatesGroupedChannelLengths(testCase)
            %TESTREPEATINGILLUMINATIONCREATESGROUPEDCHANNELLENGTHS Group repeated entries.

            [rawFolder, repeatedDatFile] = testCase.makeRepeatedIlluminationFixture();

            ImagesClassification( ...
                rawFolder, ...
                testCase.TempSaveFolder, ...
                1, 1, 0, ...
                'backupOpts', 'ERASE');

            S = load(fullfile(testCase.TempSaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            AcqInfoStream = S.AcqInfoStream;
            importedChannels = AcqInfoStream.ImportedChannels;

            testCase.verifyTrue(any(strcmp({importedChannels.DatFile}, repeatedDatFile)), ...
                ['Repeated channel file not found: ' repeatedDatFile]);

            idxRep = find(strcmp({importedChannels.DatFile}, repeatedDatFile), 1, 'first');
            % .dat header Phase 7b: AcqInfos.mat keeps the raw info.txt
            % values and no Length; the base timeline is a non-repeated
            % channel's header.
            rawInfo = ReadInfoFile(rawFolder);
            testCase.verifyEqual(double(AcqInfoStream.FrameRateHz), double(rawInfo.FrameRateHz));
            testCase.verifyFalse(isfield(AcqInfoStream, 'Length'));
            testCase.verifyFalse(isfield(AcqInfoStream, 'Datatype'));
            idxBase = find(~strcmp({importedChannels.DatFile}, repeatedDatFile), 1, 'first');
            baseInfo = loadMetaData(fullfile(testCase.TempSaveFolder, importedChannels(idxBase).DatFile));
            baseLength = datAxisSize(baseInfo, 'T');
            baseFreq = double(baseInfo.frameRateHz);

            testCase.verifyEqual(double(importedChannels(idxRep).Length), 2 * baseLength);
            testCase.verifyEqual(double(importedChannels(idxRep).FrameRateHz), 2 * baseFreq, ...
                'AbsTol', max(1e-9, abs(2 * baseFreq) * 1e-9));

            for iChan = 1:numel(importedChannels)
                if iChan == idxRep
                    continue
                end
                testCase.verifyEqual(double(importedChannels(iChan).Length), baseLength);
                testCase.verifyEqual(double(importedChannels(iChan).FrameRateHz), baseFreq, ...
                    'AbsTol', max(1e-9, abs(baseFreq) * 1e-9));
            end

            repPath = fullfile(testCase.TempSaveFolder, repeatedDatFile);
            Info = loadMetaData(repPath);
            testCase.verifyEqual(datAxisSize(Info, 'T'), 2 * baseLength);
            testCase.verifyEqual(double(Info.frameRateHz), 2 * baseFreq, ...
                'AbsTol', max(1e-9, abs(2 * baseFreq) * 1e-9));
        end

        function testMultiCamMissingCam2FilesErrors(testCase)
            %TESTMULTICAMMISSINGCAM2FILESERRORS Error clearly when a MultiCam
            %acquisition has no imgCam2_*.bin files (F-19).

            rawFolder = testCase.TempRawFolder;
            if isfolder(rawFolder)
                rmdir(rawFolder, 's');
            end
            mkdir(rawFolder);
            copyfile(fullfile(testCase.FixtureFolder, '*'), rawFolder);

            infoPath = fullfile(rawFolder, 'info.txt');
            testCase.assertTrue(isfile(infoPath), 'Fixture folder must contain info.txt.');
            txt = fileread(infoPath);
            % Any "Illumination<n>CameraIdx" entry forces AcqInfoStream.MultiCam
            % to true (see ReadInfoFile), regardless of actual Cam2 files.
            txt = sprintf('%s\nIllumination1CameraIdx: 1\n', txt);
            fid = fopen(infoPath, 'w');
            testCase.assertNotEqual(fid, -1, 'Failed to open synthetic info.txt for writing.');
            fwrite(fid, txt, 'char');
            fclose(fid);

            didThrow = false;
            try
                ImagesClassification( ...
                    rawFolder, ...
                    testCase.TempSaveFolder, ...
                    1, 1, 0, ...
                    'backupOpts', 'ERASE');
            catch ME
                didThrow = true;
                testCase.verifyTrue(contains(ME.message, 'Camera #2'), ...
                    'MultiCam acquisitions missing Cam2 binaries should error with a Cam2-specific message.');
            end
            testCase.verifyTrue(didThrow, ...
                'A MultiCam acquisition missing Cam2 binaries should raise an error.');
        end

        function testNonNumericFrameIndexErrors(testCase)
            %TESTNONNUMERICFRAMEINDEXERRORS Reject a non-numeric img_*.bin index (F-19).

            rawFolder = testCase.TempRawFolder;
            if isfolder(rawFolder)
                rmdir(rawFolder, 's');
            end
            mkdir(rawFolder);
            copyfile(fullfile(testCase.FixtureFolder, '*'), rawFolder);

            seqFiles = dir(fullfile(rawFolder, 'img_*.bin'));
            testCase.assumeNotEmpty(seqFiles, 'Fixture folder must contain img_*.bin files.');
            names = sort({seqFiles.name})';
            targetName = names{1};
            corruptedName = regexprep(targetName, '(\d)(?=\.bin$)', 'X');
            testCase.assertNotEqual(corruptedName, targetName, ...
                'Fixture filename must contain a trailing digit to corrupt.');
            movefile(fullfile(rawFolder, targetName), fullfile(rawFolder, corruptedName));

            didThrow = false;
            try
                ImagesClassification( ...
                    rawFolder, ...
                    testCase.TempSaveFolder, ...
                    1, 1, 0, ...
                    'backupOpts', 'ERASE');
            catch ME
                didThrow = true;
                testCase.verifyTrue(contains(ME.message, 'Image binary files missing'), ...
                    'A non-numeric frame index should raise the file-set validation error.');
            end
            testCase.verifyTrue(didThrow, ...
                'A non-numeric frame index should raise an error.');
        end

        function testEmbeddedImgTokenIsRejectedByAnchoredPattern(testCase)
            %TESTEMBEDDEDIMGTOKENISREJECTEDBYANCHOREDPATTERN Anchored parsing
            %rejects filenames with an "img_" token embedded in the middle
            %(not as the true prefix/suffix), matching master's behavior (F-19).
            %
            %   Under the previous unanchored erase(), a filename like
            %   "img_00000img_.bin" would have every "img_"/".bin" occurrence
            %   stripped wherever it appears, silently recovering the correct
            %   numeric index ("00000") and letting the corrupted filename
            %   pass validation. The restored anchored regexprep only strips
            %   the true leading/trailing occurrences, leaving the embedded
            %   token in place, so the value is non-numeric and rejected by
            %   the NaN check.

            rawFolder = testCase.TempRawFolder;
            if isfolder(rawFolder)
                rmdir(rawFolder, 's');
            end
            mkdir(rawFolder);
            copyfile(fullfile(testCase.FixtureFolder, '*'), rawFolder);

            seqFiles = dir(fullfile(rawFolder, 'img_*.bin'));
            testCase.assumeNotEmpty(seqFiles, 'Fixture folder must contain img_*.bin files.');
            names = sort({seqFiles.name})';
            targetName = names{1};
            numeral = regexp(targetName, '^img_(\d+)\.bin$', 'tokens', 'once');
            testCase.assumeNotEmpty(numeral, ...
                'Fixture filenames must match the img_<digits>.bin naming scheme.');
            corruptedName = ['img_' numeral{1} 'img_.bin'];
            movefile(fullfile(rawFolder, targetName), fullfile(rawFolder, corruptedName));

            didThrow = false;
            try
                ImagesClassification( ...
                    rawFolder, ...
                    testCase.TempSaveFolder, ...
                    1, 1, 0, ...
                    'backupOpts', 'ERASE');
            catch ME
                didThrow = true;
                testCase.verifyTrue(contains(ME.message, 'Image binary files missing'), ...
                    'An embedded img_ token should raise the file-set validation error.');
            end
            testCase.verifyTrue(didThrow, ...
                ['A filename with an embedded "img_" token should be rejected, not silently ' ...
                 'accepted by stripping the substring wherever it occurs.']);
        end
    end

    methods (Access = private)
        function [rawFolder, repeatedDatFile] = makeRepeatedIlluminationFixture(testCase)
            %MAKEREPEATEDILLUMINATIONFIXTURE Copy fixture and repeat one illumination.

            rawFolder = testCase.TempRawFolder;
            if isfolder(rawFolder)
                rmdir(rawFolder, 's');
            end
            mkdir(rawFolder);
            copyfile(fullfile(testCase.FixtureFolder, '*'), rawFolder);

            infoPath = fullfile(rawFolder, 'info.txt');
            testCase.assertTrue(isfile(infoPath), 'Fixture folder must contain info.txt.');

            txt = fileread(infoPath);
            expr = '(?m)^Illumination(\d+):\s*([^\r\n]*)';
            tok = regexp(txt, expr, 'tokens');
            testCase.assumeGreaterThanOrEqual(numel(tok), 3, ...
                'Repeating-illumination test requires at least three illumination entries.');

            illumIdx = cellfun(@(x) str2double(x{1}), tok);
            illumValues = cellfun(@(x) strtrim(x{2}), tok, 'UniformOutput', false);
            datTags = cellfun(@(x) testCase.getDatTagFromColor(x), illumValues, 'UniformOutput', false);

            candidateIdx = find(cellfun(@(x) sum(strcmp(datTags, x)) == 1, datTags), 1, 'first');
            testCase.assumeNotEmpty(candidateIdx, ...
                'Repeating-illumination test requires at least one uniquely tagged illumination entry.');

            targetIdx = numel(illumIdx);
            if targetIdx == candidateIdx
                targetIdx = numel(illumIdx) - 1;
            end

            repeatedColor = illumValues{candidateIdx};
            repeatedDatFile = datTags{candidateIdx};

            pattern = sprintf('(?m)^(Illumination%d:\\s*)[^\\r\\n]*', illumIdx(targetIdx));
            txt = regexprep(txt, pattern, ['$1' repeatedColor], 'once');

            fid = fopen(infoPath, 'w');
            testCase.assertNotEqual(fid, -1, 'Failed to open synthetic info.txt for writing.');
            cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid, txt, 'char');
        end
    end

    methods (Static, Access = private)
        function dTag = getDatTagFromColor(colorName)
            %GETDATTAGFROMCOLOR Match ImagesClassification output naming.

            if contains(colorName, {'red', 'green'}, 'IgnoreCase', true)
                dTag = [lower(char(string(colorName))) '.dat'];
            elseif contains(colorName, 'amber', 'IgnoreCase', true)
                dTag = 'yellow.dat';
            elseif contains(colorName, 'fluo', 'IgnoreCase', true)
                waveTag = regexp(colorName, '[0-9]{3}', 'match');
                if ~isempty(waveTag)
                    dTag = ['fluo_' waveTag{1} '.dat'];
                else
                    dTag = 'fluo.dat';
                end
            else
                dTag = 'speckle.dat';
            end
        end

        function datLength = getDatLengthFromFile(datPath, ~, ~)
            %GETDATLENGTHFROMFILE T of a headered Y-X-T .dat file.
            %
            %   Channel files are headered since .dat header Phase 4d, so T
            %   comes from the header instead of the file size.

            hdr = readDatHeader(datPath);
            datLength = hdr.dimSizes(strcmp(hdr.dimNames, 'T'));
        end
    end
end

function h = iSha256(values)
% SHA-256 of the values' bytes in memory (column-major) order.
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(values(:), 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
