classdef TestValidateTaskChanges < matlab.unittest.TestCase
%TESTVALIDATETASKCHANGES Focused tests for the consolidated task validator.

    properties
        RepoRoot
        ToolsDirectory
    end

    methods (TestClassSetup)
        function configurePath(testCase)
            testCase.RepoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.ToolsDirectory = fullfile(testCase.RepoRoot, '_tools');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testCase.ToolsDirectory));
        end
    end

    methods (Test)
        function testSuccessfulFocusedValidation(testCase)
            report = validateTaskChanges( ...
                'ChangedFiles', "test/Tools/fixtures/validatorSampleFunction.m", ...
                'TestTargets', "test/Tools/fixtures/TestValidatorPassing.m", ...
                'Mode', "targeted");

            testCase.verifyEqual(report.status, "passed");
            testCase.verifyEqual(report.counts.testsRun, 1);
            testCase.verifyTrue(isfile(report.diagnosticsPath));
        end

        function testMockFoldersStayOffThePath(testCase)
            % DFR-20260928-001: mock folders (e.g. a memory.m shadowing the
            % built-in) must not be put on the path by the validator, and a
            % mock folder left there by an earlier run is removed.
            originalPath = path;
            testCase.addTeardown(@() path(originalPath));
            staleMock = fullfile(testCase.RepoRoot, 'test', 'subFunc', 'calculateMaxChunkSize', 'mocks');
            addpath(staleMock);

            report = validateTaskChanges( ...
                'TestTargets', "test/Tools/fixtures/TestValidatorPassing.m", ...
                'Mode', "targeted");

            testCase.verifyEqual(report.status, "passed");
            entries = strsplit(path, pathsep);
            isMock = cellfun(@(folder) any(strcmpi(strsplit(folder, filesep), 'mocks')), entries);
            testCase.verifyEmpty(entries(isMock), 'no mock folder may remain on the path');
            testCase.verifySubstring(which('memory'), matlabroot, ...
                'memory must resolve to the MATLAB built-in');
        end

        function testMissingTargetIsInsufficientValidation(testCase)
            report = validateTaskChanges( ...
                'TestTargets', "test/Tools/fixtures/DoesNotExist.m");

            testCase.verifyEqual(report.status, "insufficient_validation");
            testCase.verifyNotEmpty(report.failures);
        end

        function testNoApplicableTestsIsInsufficientValidation(testCase)
            report = validateTaskChanges( ...
                'TestMethods', ...
                "test/Tools/fixtures/TestValidatorPassing.m::testDoesNotExist");

            testCase.verifyEqual(report.status, "insufficient_validation");
            testCase.verifyEqual(report.counts.testsRun, 0);
        end

        function testFailingTestIsReported(testCase)
            report = validateTaskChanges( ...
                'TestTargets', "test/Tools/fixtures/TestValidatorFailing.m");

            testCase.verifyEqual(report.status, "failed");
            testCase.verifyEqual(report.counts.testsFailed, 1);
            testCase.verifyNotEmpty(report.failures);
        end

        function testExplicitMethodTargeting(testCase)
            report = validateTaskChanges( ...
                'TestMethods', ...
                "test/Tools/fixtures/TestValidatorPassing.m::testSampleFunction");

            testCase.verifyEqual(report.status, "passed");
            testCase.verifyEqual(report.counts.testsRun, 1);
            testCase.verifyEqual(report.tests.targets, ...
                replace(string(fullfile(testCase.RepoRoot, ...
                'test/Tools/fixtures/TestValidatorPassing.m')), '\', '/') + ...
                "::testSampleFunction");
        end

        function testDuplicateTargetsAreRejected(testCase)
            target = "test/Tools/fixtures/TestValidatorPassing.m";
            report = validateTaskChanges('TestTargets', [target; target]);

            testCase.verifyEqual(report.status, "insufficient_validation");
            testCase.verifyEqual(report.counts.testsRun, 0);
            testCase.verifyTrue(contains(report.failures, "Duplicate validation target"));
        end

        function testOverlappingTargetsAreRejected(testCase)
            report = validateTaskChanges( ...
                'TestTargets', ["test/Tools/fixtures"; ...
                "test/Tools/fixtures/TestValidatorPassing.m"]);

            testCase.verifyEqual(report.status, "insufficient_validation");
            testCase.verifyEqual(report.counts.testsRun, 0);
            testCase.verifyTrue(contains(report.failures, "Overlapping validation targets"));
        end

        function testExecutionExceptionIsInfrastructureFailure(testCase)
            report = validateTaskChanges( ...
                'TestTargets', "test/Tools/fixtures/TestValidatorThrows.m");

            testCase.verifyEqual(report.status, "infrastructure_failure");
            testCase.verifyNotEmpty(report.tests.infrastructureFailures);
        end

        function testEmptyRequestSetIsInsufficientValidation(testCase)
            report = validateTaskChanges();

            testCase.verifyEqual(report.status, "insufficient_validation");
            testCase.verifyEqual(report.counts.testsRun, 0);
            testCase.verifyNotEmpty(report.deferredValidation);
        end
    end
end
