classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestNormalizeLPFZeroCenteredInput < matlab.unittest.TestCase
    %TESTNORMALIZELPFZEROCENTEREDINPUT normalizeLPF warns on zero-centered input.
    %
    %   With Normalize = true the filtered signal is divided by its low-pass
    %   baseline. Zero-centered input (already baseline-normalized, z-scored
    %   or baseline-subtracted) has a baseline near 0, so the division gives
    %   spikes of arbitrary size. The function warns once, with the id
    %   Umitoolbox:normalizeLPF:ZeroCenteredInput, and does not change its
    %   result. Synthetic non-square data (Y ~= X, distinct E): no dataset.

    properties (Constant)
        WarningId = 'Umitoolbox:normalizeLPF:ZeroCenteredInput'
        FrameRateHz = 10
        BaselineCutoffHz = 0.05   % period 20 s, shorter than the 40 s recordings
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

            testCase.verifyWarningFree(@() testCase.runLPF(raw));
        end

        function zeroCenteredArrayWarns(testCase)
            raw = makeIntensityData(6, 9, 400);
            centered = raw - mean(raw, 3);

            testCase.verifyWarning(@() testCase.runLPF(centered), testCase.WarningId);
        end

        function baselineNormalizedArrayWarns(testCase)
            % DeltaR/R0 of raw data: the baseline is about 0.
            raw = makeIntensityData(6, 9, 400);
            baseline = median(raw(:, :, 1:80), 3);
            normalized = (raw - baseline) ./ baseline;

            testCase.verifyWarning(@() testCase.runLPF(normalized), testCase.WarningId);
        end

        function eventSplitArrayWarnsOnceForAllTrials(testCase)
            raw = makeIntensityData(6, 9, 400, 3);
            centered = raw - mean(raw, 3); %#ok<NASGU> used inside the evalc string

            output = evalc('testCase.runLPF(centered);');

            testCase.verifyEqual(numel(strfind(output, 'looks zero-centered')), 1, ...
                'One warning for the whole call, not one per trial.');
        end

        function eventSplitRawArrayDoesNotWarn(testCase)
            raw = makeIntensityData(6, 9, 400, 3);

            testCase.verifyWarningFree(@() testCase.runLPF(raw));
        end

        function withoutNormalizeItNeverWarns(testCase)
            raw = makeIntensityData(6, 9, 400);
            centered = raw - mean(raw, 3);

            testCase.verifyWarningFree(@() testCase.runLPF(centered, 'Normalize', false));
        end

        function maskedPixelsDoNotTriggerTheWarning(testCase)
            % Half of the image is masked (zeros and NaN); the rest is raw.
            raw = makeIntensityData(6, 9, 400);
            raw(1:2, :, :) = 0;
            raw(3, :, :) = NaN;

            testCase.verifyWarningFree(@() testCase.runLPF(raw));
        end

        function mostlyRawDataBelowTheShareDoesNotWarn(testCase)
            % 1 of 6 rows zero-centered: 16 % of the traces, below the 20 % share.
            raw = makeIntensityData(6, 9, 400);
            raw(1, :, :) = raw(1, :, :) - mean(raw(1, :, :), 3);

            testCase.verifyWarningFree(@() testCase.runLPF(raw));
        end

        function datFileRawDoesNotWarn(testCase)
            raw = makeIntensityData(6, 9, 400);
            file = testCase.writeDat('raw.dat', raw);

            testCase.verifyWarningFree(@() testCase.runLPF(file));
        end

        function datFileZeroCenteredWarns(testCase)
            raw = makeIntensityData(6, 9, 400);
            file = testCase.writeDat('centered.dat', raw - mean(raw, 3));

            testCase.verifyWarning(@() testCase.runLPF(file), testCase.WarningId);
        end

        function multiSlabFileWarnsOnceForTheWholeFile(testCase)
            % The counts of every slab are summed and the warning is raised
            % at the end: a file read in many slabs still gives one warning.
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'Chunk forcing relies on shadowing the PCWIN64 memory() built-in.');
            raw = makeIntensityData(6, 9, 400);
            file = testCase.writeDat('centered.dat', raw - mean(raw, 3)); %#ok<NASGU> used inside the evalc string
            testCase.forceManySlabs();

            output = evalc('testCase.runLPF(file);');

            testCase.verifySubstring(output, 'Chunk #3', ...
                'The file must have been read in several slabs.');
            testCase.verifyEqual(numel(strfind(output, 'looks zero-centered')), 1);
        end

        function warningDoesNotChangeTheResult(testCase)
            raw = makeIntensityData(6, 9, 400);
            centered = raw - mean(raw, 3);

            warning('off', testCase.WarningId);
            cleanup = onCleanup(@() warning('on', testCase.WarningId));
            silenced = testCase.runLPF(centered);
            clear cleanup
            warned = testCase.verifyWarning(@() testCase.runLPF(centered), testCase.WarningId);

            testCase.verifyEqual(warned, silenced);
        end
    end

    methods (Access = private)
        function out = runLPF(testCase, data, varargin)
            %RUNLPF normalizeLPF with fixed rate and a baseline cutoff that fits.
            out = normalizeLPF(data, testCase.SaveFolder, ...
                'FrameRateHz', testCase.FrameRateHz, ...
                'BaselineCutoffHz', testCase.BaselineCutoffHz, varargin{:});
        end

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

function data = makeIntensityData(ny, nx, nt, ne)
%MAKEINTENSITYDATA Raw-like single data (large positive level, small fluctuation).
%
%   Y ~= X; the optional fourth argument adds an E axis of distinct size.
if nargin < 4
    ne = 1;
end
[yy, xx, tt, ee] = ndgrid(1:ny, 1:nx, 1:nt, 1:ne);
level = 1000 + 10 * yy + 5 * xx + 2 * ee;
wave = 20 * sin(2 * pi * tt / 17 + 0.3 * yy + 0.2 * xx + ee);
data = single(level + wave);
if ne == 1
    data = reshape(data, ny, nx, nt);
end
end
