classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestDatImageSource < matlab.unittest.TestCase
    %TESTDATIMAGESOURCE Focused regression tests for DatImageSource.
    %
    %   Covers the rig behaviour of loading (rig-ledger redesign, Phase 5):
    %   opening a .dat never writes rig data (no AcqInfos.mat stamping, no
    %   ledger creation), and an invalid folder rig ledger is reported on
    %   RigAssociationIssue without blocking the load.

    properties
        TempRoot
    end

    methods (TestMethodSetup)
        function createTemporaryRoot(testCase)
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
            testCase.addTeardown(@() rmdir(testCase.TempRoot, 's'));
        end
    end

    methods (Test)
        function testConstructionSucceedsWithoutAcqInfos(testCase)
            % Baseline: no AcqInfos.mat means Rig association is never
            % attempted, and RigAssociationIssue stays empty.
            [datFile, ~] = testCase.writeContinuousDatFixture(3, 4, 5, 10);

            src = DatImageSource(datFile);

            testCase.verifyEqual(src.getSize(), [3, 4, 5, 1]);
            testCase.verifyEmpty(src.RigAssociationIssue);
        end

        function testLoadingNeverWritesRigData(testCase)
            % A folder with AcqInfos.mat but no rig association and no
            % ledger used to be stamped with the Active rig on open.
            [datFile, folderPath] = testCase.writeContinuousDatFixture(3, 4, 5, 10);
            AcqInfoStream = struct('Width', 4, 'Height', 3);
            acqPath = fullfile(folderPath, 'AcqInfos.mat');
            save(acqPath, 'AcqInfoStream', '-mat');
            before = dir(folderPath);

            src = testCase.verifyWarningFree(@() testCase.constructSource(datFile));

            testCase.verifyEqual(src.getSize(), [3, 4, 5, 1]);
            testCase.verifyEmpty(src.RigAssociationIssue);
            loaded = load(acqPath, 'AcqInfoStream');
            testCase.verifyFalse(any(isfield(loaded.AcqInfoStream, {'rigUUID', 'rigID'})));
            testCase.verifyFalse(UMITRigLedger.exists(folderPath));
            after = dir(folderPath);
            testCase.verifyEqual({after.name}, {before.name});
            testCase.verifyEqual([after.datenum], [before.datenum]);
        end

        function testInvalidRigLedgerIsReportedWithoutBlockingLoad(testCase)
            [datFile, folderPath] = testCase.writeContinuousDatFixture(3, 4, 5, 10);
            fid = fopen(UMITRigLedger.filePath(folderPath), 'w');
            fprintf(fid, '{not json');
            fclose(fid);

            testCase.verifyWarning(@() testCase.constructSource(datFile), ...
                'Umitoolbox:DatImageSource:invalidRigLedger');

            src = testCase.constructSource(datFile);
            testCase.verifyEqual(src.getSize(), [3, 4, 5, 1]);
            testCase.verifySubstring(src.RigAssociationIssue, 'rig ledger');
        end
    end

    methods (Access = private)
        function src = constructSource(~, datFile)
            src = DatImageSource(datFile);
        end

        function [datFile, folderPath] = writeContinuousDatFixture( ...
                testCase, ny, nx, nt, freqHz)
            %WRITECONTINUOUSDATFIXTURE Write a minimal continuous [Y,X,T]
            %.dat file with a legacy sidecar, bypassing AcqInfoStream-based
            %timeline resolution entirely.

            folderPath = fullfile(testCase.TempRoot, ...
                char(java.util.UUID.randomUUID()));
            mkdir(folderPath);

            datFile = fullfile(folderPath, 'data.dat');
            fid = fopen(datFile, 'w');
            fwrite(fid, zeros(ny * nx * nt, 1, 'single'), 'single');
            fclose(fid);

            dim_names = {'Y', 'X', 'T'};
            datSize = [ny, nx];
            datLength = nt;
            Freq = freqHz;
            Datatype = 'single';
            save(fullfile(folderPath, 'data.mat'), ...
                'dim_names', 'datSize', 'datLength', 'Freq', 'Datatype', '-mat');
        end
    end
end
