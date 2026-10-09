classdef TestMapDat < matlab.unittest.TestCase
    %TESTMAPDAT Characterization tests for mapDat on headerless .dat files.
    %
    %   Records mapDat behavior on headerless files so that the .dat header
    %   work changes it only deliberately:
    %     - the Phase 1 reference file legacy/green.dat + AcqInfos.mat is an
    %       AcqInfos-bound file (no header, no sidecar), not a legacy file;
    %       it is rejected since .dat header Phase 5b;
    %     - a per-file sidecar .mat (the real legacy format) takes
    %       precedence over AcqInfos.mat.
    %   Headered files are covered by TestDatFormatLoading. Behaviors
    %   recorded as bugs in the Phase 0 impact inventory are intentionally
    %   not asserted.

    properties
        LegacyFolder
    end

    methods (TestClassSetup)
        function locateFixtures(testCase)
            testCase.LegacyFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'datHeader', 'fixtures', 'v1', 'legacy');
            testCase.assertTrue(isfile(fullfile(testCase.LegacyFolder, 'green.dat')), ...
                'AcqInfos-bound reference file missing: run makeDatFixturesV1 once.');
        end
    end

    methods (Test)
        function rejectsAcqInfosBoundFile(testCase)
            testCase.verifyError(@() mapDat(fullfile(testCase.LegacyFolder, 'green.dat')), ...
                'Umitoolbox:loadMetaData:acqInfosBoundUnsupported');
        end

        function legacySidecarOverridesAcqInfos(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;

            % Folder-level AcqInfos deliberately disagrees with the file.
            AcqInfoStream = struct('Width', 99, 'Height', 99, 'Length', 1, ...
                'FrameRateHz', 1, 'ExposureMsec', 5);
            save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');

            datFile = fullfile(folder, 'green.dat');
            fid = fopen(datFile, 'w');
            fwrite(fid, single(1:120), 'single');
            fclose(fid);

            datSize = [6 5]; datLength = 4; Freq = 20;
            save(fullfile(folder, 'green.mat'), 'datSize', 'datLength', 'Freq');

            [mmFile, info] = mapDat(datFile);
            testCase.verifyEqual(mmFile.Data.data, reshape(single(1:120), 6, 5, 4));
            testCase.verifyEqual(datAxisSize(info, 'Y'), 6);
            testCase.verifyEqual(datAxisSize(info, 'X'), 5);
            testCase.verifyEqual(datAxisSize(info, 'T'), 4);
            testCase.verifyEqual(info.frameRateHz, 20);
            clear mmFile
        end

        function rejectsPathThatIsNotAnExistingDatFile(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            txtFile = fullfile(folder, 'notes.txt');
            fclose(fopen(txtFile, 'w'));

            testCase.verifyError(@() mapDat(txtFile), 'MATLAB:InputParser:ArgumentFailedValidation');
            testCase.verifyError(@() mapDat(fullfile(folder, 'missing.dat')), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end
    end
end
