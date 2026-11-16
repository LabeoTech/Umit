classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestAssertSingleDat < matlab.unittest.TestCase
    %TESTASSERTSINGLEDAT The shared single-precision .dat gate (DFR-20261019-006).
    %
    %   umIT supports only single-precision .dat files. assertSingleDat is the
    %   one place that says so; these tests guard the policy and the cases it
    %   must leave alone (non-.dat paths, missing files, legacy sidecar files).

    properties
        Folder
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            addpath(genpath(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))));
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function singleHeaderedFilePasses(testCase)
            f = fullfile(testCase.Folder, 'ok.dat');
            writeTestDat(f, single(rand(4, 5, 3)), 10);
            assertSingleDat(f);
            assertSingleDat(f, 'Context', 'Anything');
        end

        function singleLegacySidecarFilePasses(testCase)
            f = fullfile(testCase.Folder, 'legacy.dat');
            writeTestDat(f, single(rand(4, 5, 3)), 10, 'Format', 'legacySidecar');
            assertSingleDat(f);
        end

        function everyOtherClassIsRefused(testCase)
            for cls = {'uint8', 'uint16', 'int16', 'uint32', 'int32', 'double'}
                f = fullfile(testCase.Folder, ['v_' cls{1} '.dat']);
                writeTestDat(f, cast(rand(4, 5, 3) * 100, cls{1}), 10);
                testCase.verifyError(@() assertSingleDat(f), 'Umitoolbox:Dat:nonSingle', cls{1});
            end
        end

        function messageNamesFileClassAndPolicy(testCase)
            f = fullfile(testCase.Folder, 'counts.dat');
            writeTestDat(f, uint16(1:60), 10, 'DimNames', {'Y', 'X', 'T'});
            try
                assertSingleDat(f);
                testCase.verifyFail('uint16 must be refused.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:Dat:nonSingle');
                testCase.verifySubstring(ME.message, '"counts.dat"');
                testCase.verifySubstring(ME.message, 'uint16');
                testCase.verifySubstring(ME.message, 'single-precision');
            end
            try
                assertSingleDat(f, 'Context', 'Step 2 output "outData"');
            catch ME
                testCase.verifySubstring(ME.message, 'Step 2 output "outData"');
            end
        end

        function infoStructIsChecked(testCase)
            assertSingleDat(struct('dataClass', 'single'));
            testCase.verifyError(@() assertSingleDat(struct('dataClass', 'int16')), ...
                'Umitoolbox:Dat:nonSingle');
            testCase.verifyError(@() assertSingleDat(struct('other', 1)), ...
                'Umitoolbox:assertSingleDat:invalidInput');
            testCase.verifyError(@() assertSingleDat(42), ...
                'Umitoolbox:assertSingleDat:invalidInput');
        end

        function filesThatAreNotDatOrMissingAreIgnored(testCase)
            notDat = fullfile(testCase.Folder, 'events.mat');
            x = 1;
            save(notDat, 'x');
            assertSingleDat(notDat);
            assertSingleDat(fullfile(testCase.Folder, 'missing.dat'));
            assertSingleDat(fullfile(testCase.Folder, 'map.umt'));
        end
    end
end
