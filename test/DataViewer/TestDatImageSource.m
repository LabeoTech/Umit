classdef TestDatImageSource < matlab.unittest.TestCase
    %TESTDATIMAGESOURCE Focused regression tests for DatImageSource.
    %
    %   Covers the non-blocking Rig-association behavior added so an
    %   unresolvable stored Rig UUID/ID (deleted Rig, dataset moved from
    %   another install, etc.) never aborts a .dat load.

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

        function testUnresolvableRigUUIDDoesNotBlockLoad(testCase)
            % An AcqInfos.mat pointing at a Rig UUID that UMITRigStore
            % cannot resolve must not abort the load: the constructor
            % catches the failure, records it on RigAssociationIssue, and
            % still finishes building a usable DatImageSource.
            [datFile, folderPath] = testCase.writeContinuousDatFixture(3, 4, 5, 10);
            AcqInfoStream = struct('rigUUID', char(java.util.UUID.randomUUID()));
            save(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream', '-mat');

            testCase.verifyWarning(@() testCase.constructSource(datFile), ...
                'Umitoolbox:DatImageSource:rigAssociationSkipped');

            src = testCase.constructSource(datFile);
            testCase.verifyEqual(src.getSize(), [3, 4, 5, 1]);
            testCase.verifyNotEmpty(src.RigAssociationIssue);
            testCase.verifySubstring(src.RigAssociationIssue, 'could not be resolved');
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
