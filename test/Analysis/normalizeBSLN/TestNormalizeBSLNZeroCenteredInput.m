classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestNormalizeBSLNZeroCenteredInput < matlab.unittest.TestCase
    %TESTNORMALIZEBSLNZEROCENTEREDINPUT normalizeBSLN warns on zero-centered input.
    %
    %   normalizeBSLN divides by the baseline of every pixel. Input that is
    %   already zero-centered (baseline-normalized without b_centerAtOne,
    %   z-scored, or baseline-subtracted) has a baseline within the noise, so
    %   the function warns once, with the id
    %   Umitoolbox:normalizeBSLN:ZeroCenteredInput, and does not change its
    %   result. Synthetic non-square data (Y ~= X): no dataset.

    properties (Constant)
        WarningId = 'Umitoolbox:normalizeBSLN:ZeroCenteredInput'
        FrameRateHz = 10
    end

    properties
        SaveFolder char
    end

    methods (TestMethodSetup)
        function createSaveFolder(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fixture.Folder;
        end
    end

    methods (Test)
        function rawArrayDoesNotWarn(testCase)
            raw = makeIntensityData(6, 9, 400);

            testCase.verifyWarningFree(@() normalizeBSLN(raw, testCase.SaveFolder));
        end

        function zeroCenteredArrayWarns(testCase)
            raw = makeIntensityData(6, 9, 400);
            centered = raw - mean(raw, 3);

            testCase.verifyWarning(@() normalizeBSLN(centered, testCase.SaveFolder), ...
                testCase.WarningId);
        end

        function alreadyNormalizedArrayWarns(testCase)
            raw = makeIntensityData(6, 9, 400);
            normalized = normalizeBSLN(raw, testCase.SaveFolder);

            testCase.verifyWarning(@() normalizeBSLN(normalized, testCase.SaveFolder), ...
                testCase.WarningId);
        end

        function dataCenteredAtOneAreNotFlagged(testCase)
            % b_centerAtOne gives a baseline of 1 with a fluctuation of a few
            % percent: the ratio is large, so nothing warns.
            raw = makeIntensityData(6, 9, 400);
            centeredAtOne = normalizeBSLN(raw, testCase.SaveFolder, 'b_centerAtOne', true);

            testCase.verifyWarningFree( ...
                @() normalizeBSLN(centeredAtOne, testCase.SaveFolder));
        end

        function maskedPixelsDoNotTriggerTheWarning(testCase)
            raw = makeIntensityData(6, 9, 400);
            raw(1:2, :, :) = 0;
            raw(3, :, :) = NaN;

            testCase.verifyWarningFree(@() normalizeBSLN(raw, testCase.SaveFolder));
        end

        function tooFewBaselineFramesSkipTheCheck(testCase)
            % 0.2 s at 10 Hz is 2 baseline frames: no meaningful spread.
            raw = makeIntensityData(6, 9, 400);
            centered = raw - mean(raw, 3);

            testCase.verifyWarningFree(@() normalizeBSLN(centered, testCase.SaveFolder, ...
                'baselineMode', 0.2, 'FrameRateHz', testCase.FrameRateHz));
        end

        function datFileRawDoesNotWarn(testCase)
            raw = makeIntensityData(6, 9, 400);
            file = testCase.writeDat('raw.dat', raw);

            testCase.verifyWarningFree(@() normalizeBSLN(file, testCase.SaveFolder));
        end

        function datFileZeroCenteredWarnsOnce(testCase)
            raw = makeIntensityData(6, 9, 400);
            file = testCase.writeDat('centered.dat', raw - mean(raw, 3)); %#ok<NASGU> used inside the evalc string

            output = evalc('normalizeBSLN(file, testCase.SaveFolder);');

            testCase.verifyEqual(numel(strfind(output, 'looks zero-centered')), 1);
        end

        function multiSlabFileWarnsOnceForTheWholeFile(testCase)
            % The counts of every slab are summed and the warning is raised
            % at the end: a file read in many slabs still gives one warning.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            raw = makeIntensityData(6, 9, 400);
            file = testCase.writeDat('centered.dat', raw - mean(raw, 3)); %#ok<NASGU> used inside the evalc string
            testCase.forceManySlabs();

            output = evalc('normalizeBSLN(file, testCase.SaveFolder);');

            testCase.verifySubstring(output, 'Chunk 3/', ...
                'The file must have been read in several slabs.');
            testCase.verifyEqual(numel(strfind(output, 'looks zero-centered')), 1);
        end

        function warningDoesNotChangeTheResult(testCase)
            raw = makeIntensityData(6, 9, 400);
            centered = raw - mean(raw, 3);

            warning('off', testCase.WarningId);
            cleanup = onCleanup(@() warning('on', testCase.WarningId));
            silenced = normalizeBSLN(centered, testCase.SaveFolder);
            clear cleanup
            warned = testCase.verifyWarning( ...
                @() normalizeBSLN(centered, testCase.SaveFolder), testCase.WarningId);

            testCase.verifyEqual(warned, silenced);
        end
    end

    methods (Access = private)
        function file = writeDat(testCase, name, data)
            file = fullfile(testCase.SaveFolder, name);
            saveData(file, data, 'DimNames', {'Y', 'X', 'T'}, ...
                'FrameRateHz', testCase.FrameRateHz);
        end

        function forceManySlabs(testCase)
            %FORCEMANYSLABS Make calculateMaxChunkSize see very little free RAM.
            projectRoot = extractBefore(mfilename('fullpath'), [filesep 'test' filesep]);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile( ...
                projectRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', '10000'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', '5000'));
        end
    end
end

function data = makeIntensityData(ny, nx, nt)
%MAKEINTENSITYDATA Raw-like single data (large positive level, small fluctuation).
%
%   Y ~= X, so a transposition shows up.
[yy, xx, tt] = ndgrid(1:ny, 1:nx, 1:nt);
level = 1000 + 10 * yy + 5 * xx;
wave = 20 * sin(2 * pi * tt / 17 + 0.3 * yy + 0.2 * xx);
data = single(level + wave);
end
