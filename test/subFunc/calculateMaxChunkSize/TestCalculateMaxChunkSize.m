classdef TestCalculateMaxChunkSize < matlab.unittest.TestCase
    %TESTCALCULATEMAXCHUNKSIZE Fallback + warning behavior on chunk-calc failure (TASK 5.1).
    %
    %   Covers both failure branches of calculateMaxChunkSize: a failed RAM
    %   query and a successful query that leaves no usable RAM. Both must
    %   fall back to the fixed 128 MB conservative budget, warn with a
    %   distinct ID, and never return [].
    %
    %   RAM is mocked by shadowing the PCWIN64 `memory` built-in with
    %   mocks/memory.m, controlled via environment variables. These tests
    %   only run on PCWIN64.

    properties (Constant)
        FallbackBudgetBytes = 128 * 1024 * 1024
    end

    properties
        RepoRoot
        MocksFolder
    end

    methods (TestClassSetup)
        function configurePath(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.RepoRoot = fileparts(fileparts(fileparts(testFolder)));
            testCase.MocksFolder = fullfile(testFolder, 'mocks');

            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.RepoRoot, 'subFunc')));
        end
    end

    methods (TestMethodSetup)
        function requireMemoryShadowingPlatform(testCase)
            testCase.assumeEqual(computer, 'PCWIN64', ...
                'RAM-query mocking relies on shadowing the PCWIN64 memory() built-in.');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testCase.MocksFolder));
        end
    end

    methods (Test)
        function testRAMQueryFailureFallsBackToProportionateBudget(testCase)
            testCase.mockMemoryFails();
            requiredBytes = 5 * testCase.FallbackBudgetBytes;
            expectedChunks = 5;

            [nChunks, warnId, warnMsg] = testCase.callWithLastWarning(requiredBytes);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateMaxChunkSize:RAMQueryFailed');
            testCase.verifySubstring(warnMsg, sprintf('%d', testCase.FallbackBudgetBytes / (1024*1024)));
            testCase.verifySubstring(warnMsg, sprintf('%d', expectedChunks));
            testCase.verifyEqual(nChunks, expectedChunks);
        end

        function testRAMQueryFailureSmallInputStaysSingleChunk(testCase)
            testCase.mockMemoryFails();
            requiredBytes = testCase.FallbackBudgetBytes / 4;

            [nChunks, warnId] = testCase.callWithLastWarning(requiredBytes);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateMaxChunkSize:RAMQueryFailed');
            testCase.verifyEqual(nChunks, 1);
        end

        function testNoUsableRAMFallsBackToProportionateBudget(testCase)
            testCase.mockMemoryFixed(1e9, 0);
            requiredBytes = 5 * testCase.FallbackBudgetBytes;
            expectedChunks = 5;

            [nChunks, warnId, warnMsg] = testCase.callWithLastWarning(requiredBytes);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateMaxChunkSize:NoUsableRAM');
            testCase.verifySubstring(warnMsg, sprintf('%d', testCase.FallbackBudgetBytes / (1024*1024)));
            testCase.verifySubstring(warnMsg, sprintf('%d', expectedChunks));
            testCase.verifyEqual(nChunks, expectedChunks);
        end

        function testNoUsableRAMSmallInputStaysSingleChunk(testCase)
            testCase.mockMemoryFixed(1e9, 0);
            requiredBytes = testCase.FallbackBudgetBytes / 4;

            [nChunks, warnId] = testCase.callWithLastWarning(requiredBytes);

            testCase.verifyEqual(warnId, 'Umitoolbox:calculateMaxChunkSize:NoUsableRAM');
            testCase.verifyEqual(nChunks, 1);
        end

        function testNoUsableRAMNeverReturnsEmpty(testCase)
            % Regression guard for the original defect: usable RAM at/below
            % zero used to return [] with no warning, silently zero-filling
            % chunked output.
            testCase.mockMemoryFixed(1e9, 0);

            nChunks = calculateMaxChunkSize(1024, 1, 0.2);

            testCase.verifyNotEmpty(nChunks);
            testCase.verifyGreaterThanOrEqual(nChunks, 1);
        end
    end

    methods (Access = private)
        function mockMemoryFails(testCase)
            import matlab.unittest.fixtures.EnvironmentVariableFixture
            testCase.applyFixture(EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fail'));
        end

        function mockMemoryFixed(testCase, totalRAM, availableRAM)
            import matlab.unittest.fixtures.EnvironmentVariableFixture
            testCase.applyFixture(EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_MODE', 'fixed'));
            testCase.applyFixture(EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_TOTAL', num2str(totalRAM, '%d')));
            testCase.applyFixture(EnvironmentVariableFixture( ...
                'UMIT_TEST_MOCK_MEMORY_AVAILABLE', num2str(availableRAM, '%d')));
        end

        function [nChunks, warnId, warnMsg] = callWithLastWarning(~, requiredBytes)
            lastwarn('');
            nChunks = calculateMaxChunkSize(requiredBytes, 1, 0.2);
            [warnMsg, warnId] = lastwarn();
        end
    end
end
