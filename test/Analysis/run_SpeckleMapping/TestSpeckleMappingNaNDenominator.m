classdef TestSpeckleMappingNaNDenominator < matlab.unittest.TestCase
    %TESTSPECKLEMAPPINGNANDENOMINATOR Regression for TASK_F1_F2: RAMSafe
    % 'spatial' must use a non-NaN-count denominator, not a fixed Nt, so its
    % output matches Standard mode bit-for-bit on NaN-containing input.
    %
    % Synthetic fixture (built here, not the shared TestingData_speckle
    % fixture) with three pixel categories along T:
    %   - all-valid
    %   - all-NaN
    %   - partially-NaN

    properties
        TempFolder = ''
        AllNaNPixel = [3, 3]
        PartialNaNPixel = [8, 8]
    end

    methods (TestMethodSetup)
        function setup(testCase)
            testCase.TempFolder = fullfile(tempdir, ...
                ['TestSpeckleMappingNaNDenominator_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);

            Ny = 12; Nx = 12; Nt = 10;

            rng(42);
            data = 0.5 + rand(Ny, Nx, Nt, 'single');

            data(testCase.AllNaNPixel(1), testCase.AllNaNPixel(2), :) = NaN;
            data(testCase.PartialNaNPixel(1), testCase.PartialNaNPixel(2), 1:4) = NaN;

            % Headered input (.dat header Phase 5a).
            writeTestDat(fullfile(testCase.TempFolder, 'speckle.dat'), data, 5);

            AcqInfoStream = struct( ...
                'Height', Ny, ...
                'Width', Nx, ...
                'Length', Nt, ...
                'FrameRateHz', 5, ...
                'Datatype', 'single');
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            if ~isempty(testCase.TempFolder) && isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testSpatialModeMatchesAcrossStandardAndRAMSafe(testCase)
            outStd = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'spatial', false, false);
            outRAM = SpeckleMapping('speckle.dat', testCase.TempFolder, 'spatial', false, false);

            mapStd = outStd.data.SpeckleMap.value;
            mapRAM = outRAM.data.SpeckleMap.value;

            testCase.verifyTrue(isequaln(mapStd, mapRAM), ...
                ['RAMSafe ''spatial'' output must match Standard output ' ...
                 'bit-for-bit given identical NaN-containing input.']);
        end

        function testTemporalModeMatchesAcrossStandardAndRAMSafe(testCase)
            outStd = SpeckleMapping(iArray(testCase), testCase.TempFolder, 'temporal', false, false);
            outRAM = SpeckleMapping('speckle.dat', testCase.TempFolder, 'temporal', false, false);

            mapStd = outStd.data.SpeckleMap.value;
            mapRAM = outRAM.data.SpeckleMap.value;

            testCase.verifyTrue(isequaln(mapStd, mapRAM), ...
                ['RAMSafe ''temporal'' output must match Standard output ' ...
                 'bit-for-bit given identical NaN-containing input.']);
        end

    end
end

function arr = iArray(testCase)
%IARRAY The fixture recording as an in-RAM Y-X-T array (Standard mode).
arr = loadData(fullfile(testCase.TempFolder, 'speckle.dat'));
end
