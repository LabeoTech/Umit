classdef TestDatFormatLoading < matlab.unittest.TestCase
    %TESTDATFORMATLOADING loadMetaData, loadData, and mapDat on every readable .dat kind.
    %
    %   Kinds covered (.dat header Phase 2, updated in Phase 5b):
    %     - headered files, described from the header alone;
    %     - legacy files with a per-file sidecar .mat (Astrocyte format);
    %     - AcqInfos-bound headerless files (no header, no sidecar), which
    %       every reader rejects since Phase 5b.
    %
    %   Reference files are read from test/subFunc/datHeader/fixtures and
    %   are never written by these tests; files that must be altered are
    %   copied into a temporary folder first.

    properties (Constant)
        AcqInfosBoundId = 'Umitoolbox:loadMetaData:acqInfosBoundUnsupported'
        SchemaFields = {'filePath', 'format', 'dataOffset', 'dataClass', 'dimNames', ...
            'dimSizes', 'frameRateHz', 'exposureMsec', 'channelName', 'writeComplete'}
        IncompleteId = 'Umitoolbox:loadMetaData:incompleteFile'
    end

    properties
        FixtureRoot
    end

    properties (TestParameter)
        headerCase = struct( ...
            'yxt_single', struct('file', 'yxt_single.dat', 'dataClass', 'single', ...
                'frameRateHz', 10, 'exposureMsec', 5, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [6 5 4], ...
                'writeComplete', true, 'values', single(1:120)), ...
            'yx_uint16', struct('file', 'yx_uint16.dat', 'dataClass', 'uint16', ...
                'frameRateHz', NaN, 'exposureMsec', NaN, 'channelName', 'amplitude', ...
                'dimNames', {{'Y', 'X'}}, 'dimSizes', [6 5], ...
                'writeComplete', true, 'values', uint16(1:30)), ...
            'yxte_single', struct('file', 'yxte_single.dat', 'dataClass', 'single', ...
                'frameRateHz', 10, 'exposureMsec', 5, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T', 'E'}}, 'dimSizes', [4 3 5 2], ...
                'writeComplete', true, 'values', single(1:120)), ...
            'yxe_single', struct('file', 'yxe_single.dat', 'dataClass', 'single', ...
                'frameRateHz', NaN, 'exposureMsec', NaN, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'E'}}, 'dimSizes', [4 3 2], ...
                'writeComplete', true, 'values', single(1:24)), ...
            'yxt_incomplete', struct('file', 'yxt_incomplete.dat', 'dataClass', 'single', ...
                'frameRateHz', 10, 'exposureMsec', 5, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [6 5 4], ...
                'writeComplete', false, 'values', single(1:120)))
    end

    methods (TestClassSetup)
        function locateFixtures(testCase)
            testCase.FixtureRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'datHeader', 'fixtures');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat')), ...
                'Header reference files missing: run makeDatFixturesV1 once.');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'green.dat')), ...
                'Legacy reference file missing: run makeLegacySidecarFixture once.');
        end
    end

    methods (Test)
        % ------------------------------------------------------ headered files
        function headeredMetadataComesFromHeader(testCase, headerCase)
            f = testCase.headerFile(headerCase);
            Info = testCase.callExpectingIncompleteWarning(headerCase, @() loadMetaData(f));

            testCase.verifyEqual(Info.filePath, f);
            testCase.verifyEqual(Info.format, 'header');
            testCase.verifyEqual(Info.dataOffset, 512);
            testCase.verifyEqual(Info.dataClass, headerCase.dataClass);
            testCase.verifyEqual(Info.dimNames, headerCase.dimNames);
            testCase.verifyEqual(Info.dimSizes, headerCase.dimSizes);
            testCase.verifyEqual(Info.frameRateHz, headerCase.frameRateHz);
            testCase.verifyEqual(Info.exposureMsec, headerCase.exposureMsec);
            testCase.verifyEqual(Info.channelName, headerCase.channelName);
            testCase.verifyEqual(Info.writeComplete, headerCase.writeComplete);
        end

        function headeredInfoHasOnlySchemaFields(testCase, headerCase)
            % .dat header Phase 8a: the compatibility block of old names is
            % removed; a headered file's Info is exactly the schema.
            f = testCase.headerFile(headerCase);
            Info = testCase.callExpectingIncompleteWarning(headerCase, @() loadMetaData(f));
            testCase.verifyEqual(sort(fieldnames(Info)), sort({'filePath', 'format', 'dataOffset', 'dataClass', 'dimNames', 'dimSizes', 'frameRateHz', 'exposureMsec', 'channelName', 'writeComplete'}'));
        end

        function loadDataReadsHeaderedValues(testCase, headerCase)
            f = testCase.headerFile(headerCase);
            data = testCase.callExpectingIncompleteWarning(headerCase, @() loadData(f));
            testCase.verifyClass(data, headerCase.dataClass);
            testCase.verifyEqual(data, reshape(headerCase.values, headerCase.dimSizes));
        end

        function mapDatMapsHeaderedValues(testCase, headerCase)
            f = testCase.headerFile(headerCase);
            mm = testCase.callExpectingIncompleteWarning(headerCase, @() mapDat(f));
            testCase.verifyEqual(mm.Offset, 512);
            testCase.verifyEqual(mm.Data.data, reshape(headerCase.values, headerCase.dimSizes));
            clear mm
        end

        function headeredFileIsSelfContained(testCase)
            folder = testCase.tempFolder();
            copyfile(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat'), fullfile(folder, 'green.dat'));
            testCase.assertFalse(isfile(fullfile(folder, 'AcqInfos.mat')));

            data = loadData(fullfile(folder, 'green.dat'));
            testCase.verifyEqual(data, reshape(single(1:120), 6, 5, 4));
        end

        function headerWinsOverAcqInfosAndSidecar(testCase)
            folder = testCase.tempFolder();
            copyfile(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat'), fullfile(folder, 'green.dat'));
            testCase.writeConflictingMetadata(folder);

            Info = loadMetaData(fullfile(folder, 'green.dat'));
            testCase.verifyEqual(Info.format, 'header');
            testCase.verifyEqual(Info.dimSizes, [6 5 4]);
            testCase.verifyEqual(Info.frameRateHz, 10);
            testCase.verifyEqual(datAxisSize(Info, 'Y'), 6);
        end

        function incompleteFileWarnsOnceAndLoads(testCase)
            f = fullfile(testCase.FixtureRoot, 'v1', 'yxt_incomplete.dat');
            testCase.verifyThat(@() loadData(f), matlab.unittest.constraints.IssuesWarnings( ...
                {testCase.IncompleteId}, 'RespectingCount', true));
            data = testCase.verifyWarning(@() loadData(f), testCase.IncompleteId);
            testCase.verifyEqual(data, reshape(single(1:120), 6, 5, 4));
        end

        function sizeMismatchOfCompleteFileErrors(testCase)
            folder = testCase.tempFolder();
            f = fullfile(folder, 'green.dat');
            bytes = testCase.readBytes(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat'));
            testCase.writeBytes(f, bytes(1:end-4));

            id = 'Umitoolbox:validateDatHeader:fileSizeMismatch';
            testCase.verifyError(@() loadMetaData(f), id);
            testCase.verifyError(@() loadData(f), id);
            testCase.verifyError(@() mapDat(f), id);
        end

        function corruptHeaderNeverFallsBack(testCase)
            folder = testCase.tempFolder();
            f = fullfile(folder, 'green.dat');
            bytes = testCase.readBytes(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat'));
            bytes(9:10) = uint8([2 0]);   % headerVersion 2
            testCase.writeBytes(f, bytes);
            testCase.writeUsableHeaderlessMetadata(folder);

            id = 'Umitoolbox:decodeDatHeader:corruptHeader';
            testCase.verifyError(@() loadMetaData(f), id);
            testCase.verifyError(@() loadData(f), id);
            testCase.verifyError(@() mapDat(f), id);
        end

        % ------------------------------------------------------ headerless files
        function legacySidecarReferenceOpens(testCase)
            f = fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'green.dat');
            testCase.assertFalse(isfile(fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'AcqInfos.mat')));

            Info = loadMetaData(f);
            testCase.verifyEqual(Info.format, 'legacySidecar');
            testCase.verifyEqual(Info.dataOffset, 0);
            testCase.verifyEqual(Info.dataClass, 'single');
            testCase.verifyEqual(Info.dimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(Info.dimSizes, [6 5 4]);
            testCase.verifyEqual(Info.frameRateHz, 10);
            testCase.verifyEqual(Info.exposureMsec, NaN);
            testCase.verifyEqual(Info.channelName, '');
            testCase.verifyTrue(Info.writeComplete);
            testCase.verifyEqual(sort(fieldnames(Info)), sort({'filePath', 'format', 'dataOffset', 'dataClass', 'dimNames', 'dimSizes', 'frameRateHz', 'exposureMsec', 'channelName', 'writeComplete'}'), ...
                'a legacy sidecar file''s Info is exactly the schema (Phase 8a)');

            testCase.verifyEqual(loadData(f), reshape(single(1:120), 6, 5, 4));
            mm = mapDat(f);
            testCase.verifyEqual(mm.Offset, 0);
            testCase.verifyEqual(mm.Data.data, reshape(single(1:120), 6, 5, 4));
            clear mm
        end

        function acqInfosBoundReferenceIsRejected(testCase)
            % The frozen Phase 1 reference: headerless, no sidecar, next to
            % an AcqInfos.mat that fully describes it.
            f = fullfile(testCase.FixtureRoot, 'v1', 'legacy', 'green.dat');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'v1', 'legacy', 'AcqInfos.mat')));

            testCase.verifyError(@() loadMetaData(f), testCase.AcqInfosBoundId);
            try
                loadMetaData(f);
                message = '';
            catch ME
                message = ME.message;
            end
            testCase.verifySubstring(message, f);
            testCase.verifySubstring(message, 'Re-import the raw data');
        end

        function rejectionIsTheSameFromEveryReader(testCase)
            % A copy without AcqInfos.mat, so DatImageSource does not touch
            % the rig store; the rejection does not depend on AcqInfos.mat.
            folder = testCase.tempFolder();
            f = fullfile(folder, 'green.dat');
            copyfile(fullfile(testCase.FixtureRoot, 'v1', 'legacy', 'green.dat'), f);

            id = testCase.AcqInfosBoundId;
            testCase.verifyError(@() loadMetaData(f), id);
            testCase.verifyError(@() loadData(f), id);
            testCase.verifyError(@() mapDat(f), id);
            testCase.verifyError(@() spatialSlabIO('open', f), id);
            testCase.verifyError(@() DatImageSource(f), id);
        end

        function everyKindReturnsTheFullSchema(testCase)
            files = {fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat'), ...
                fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'green.dat')};
            for k = 1:numel(files)
                Info = loadMetaData(files{k});
                missing = setdiff(testCase.SchemaFields, fieldnames(Info));
                testCase.verifyEmpty(missing, sprintf('%s lacks schema fields', files{k}));
            end
        end

        function speckleFilesUseTheirOwnExposure(testCase)
            % Legacy sidecar files next to a folder-level AcqInfos.mat with
            % a general and a speckle exposure: they still open, and
            % exposureMsec is each file's own exposure. The old names are
            % gone since Phase 8a.
            folder = testCase.tempFolder();
            channel = @(name) struct('DatFile', [name '.dat'], 'Tag', name, 'Color', name, ...
                'Length', 4, 'FrameRateHz', 10, 'ExposureMsec', 5, 'CamIdx', 1);
            AcqInfoStream = struct('Width', 5, 'Height', 6, 'Length', 4, 'FrameRateHz', 10, ...
                'ExposureMsec', 5, 'ExposureSpeckleMsec', 2);
            AcqInfoStream.ImportedChannels = [channel('green'), channel('speckle')];
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            for name = {'green', 'speckle'}
                writeTestDat(fullfile(folder, [name{1} '.dat']), ...
                    reshape(single(1:120), 6, 5, 4), 10, 'Format', 'legacySidecar');
            end

            speckle = loadMetaData(fullfile(folder, 'speckle.dat'));
            testCase.verifyEqual(speckle.format, 'legacySidecar');
            testCase.verifyEqual(speckle.dimSizes, [6 5 4]);
            testCase.verifyEqual(speckle.frameRateHz, 10);
            testCase.verifyEqual(speckle.exposureMsec, 2);
            testCase.verifyFalse(any(isfield(speckle, {'ExposureSpeckleMsec', 'ExposureMsec', ...
                'MetadataSource', 'CamIdx'})));

            green = loadMetaData(fullfile(folder, 'green.dat'));
            testCase.verifyEqual(green.format, 'legacySidecar');
            testCase.verifyEqual(green.exposureMsec, 5);
            testCase.verifyEqual(loadData(fullfile(folder, 'green.dat')), ...
                reshape(single(1:120), 6, 5, 4));
        end

        % ------------------------------------------------------ datAxisSize
        function datAxisSizeReturnsSizeOrZero(testCase)
            Info = struct('dimNames', {{'Y', 'X', 'E'}}, 'dimSizes', [4 3 2]);
            testCase.verifyEqual(datAxisSize(Info, 'X'), 3);
            testCase.verifyEqual(datAxisSize(Info, "E"), 2);
            testCase.verifyEqual(datAxisSize(Info, 'T'), 0);

            legacy = struct('dimNames', {{'E', 'Y', 'X', 'T'}}, 'dimSizes', [3 112 112 320]);
            testCase.verifyEqual(datAxisSize(legacy, 'Y'), 112);
            testCase.verifyEqual(datAxisSize(legacy, 'T'), 320);
        end

        function datAxisSizeRejectsInvalidInput(testCase)
            id = 'Umitoolbox:datAxisSize:invalidInput';
            testCase.verifyError(@() datAxisSize(struct('Height', 6), 'Y'), id);
            testCase.verifyError(@() datAxisSize(struct('dimNames', {{'Y'}}, 'dimSizes', 6), 1), id);
        end
    end

    methods (Access = private)
        function f = headerFile(testCase, c)
            f = fullfile(testCase.FixtureRoot, 'v1', c.file);
        end

        function out = callExpectingIncompleteWarning(testCase, c, fcn)
            % Complete files must load silently; incomplete ones warn once.
            if c.writeComplete
                out = testCase.verifyWarningFree(fcn);
            else
                out = testCase.verifyWarning(fcn, testCase.IncompleteId);
            end
        end

        function folder = tempFolder(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end

        function writeConflictingMetadata(~, folder)
            % AcqInfos.mat and sidecar that disagree with the header.
            AcqInfoStream = struct('Width', 99, 'Height', 99, 'Length', 1, ...
                'FrameRateHz', 1, 'ExposureMsec', 1);
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            datSize = [99 99];
            datLength = 1;
            Freq = 1;
            Datatype = 'single';
            dim_names = {'Y', 'X', 'T'};
            save(fullfile(folder, 'green.mat'), 'datSize', 'datLength', 'Freq', 'Datatype', 'dim_names');
        end

        function writeUsableHeaderlessMetadata(~, folder)
            % Metadata that would let the file open as headerless if the
            % reader (wrongly) fell back.
            AcqInfoStream = struct('Width', 5, 'Height', 6, 'Length', 4, ...
                'FrameRateHz', 10, 'ExposureMsec', 5);
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
            datSize = [6 5];
            datLength = 4;
            Freq = 10;
            Datatype = 'single';
            dim_names = {'Y', 'X', 'T'};
            save(fullfile(folder, 'green.mat'), 'datSize', 'datLength', 'Freq', 'Datatype', 'dim_names');
        end

        function bytes = readBytes(~, f)
            fid = fopen(f, 'r');
            bytes = fread(fid, inf, '*uint8').';
            fclose(fid);
        end

        function writeBytes(~, f, bytes)
            fid = fopen(f, 'w');
            fwrite(fid, bytes, 'uint8');
            fclose(fid);
        end
    end
end
