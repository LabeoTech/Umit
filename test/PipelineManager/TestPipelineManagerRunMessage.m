classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerRunMessage < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERRUNMESSAGE One-sentence run outcome for the end-of-run dialog.
    %
    %   PipelineManager.describeExecutionResult(result) turns an executePipeline
    %   result (or an MException) into a severity and a short message. It reads
    %   result.status and the Status / Messages_short columns of
    %   result.globalPipeLog, so these tests build minimal synthetic results.

    methods (Test)
        function completedIsSuccess(testCase)
            info = PipelineManager.describeExecutionResult(iResult("completed", {}, {}));
            testCase.verifyEqual(info.severity, 'success');
            testCase.verifyEqual(info.message, 'Pipeline completed successfully.');
        end

        function failedCarriesFirstErrorWithoutFunctionDetail(testCase)
            msg = 'Step "GSR" failed: bad input (Function: runGSR, Line: 42) | second problem';
            info = PipelineManager.describeExecutionResult( ...
                iResult("failed", {"Failed"}, {msg}));
            testCase.verifyEqual(info.severity, 'error');
            testCase.verifyEqual(info.message, 'Pipeline failed: Step "GSR" failed: bad input');
        end

        function failedWithoutUsableMessageStillReads(testCase)
            info = PipelineManager.describeExecutionResult( ...
                iResult("failed", {"Failed"}, {'No Errors'}));
            testCase.verifyEqual(info.severity, 'error');
            testCase.verifyEqual(info.message, 'Pipeline failed.');
        end

        function severalFailedFoldersAreCountedAndFirstMessageUsed(testCase)
            info = PipelineManager.describeExecutionResult(iResult("failed", ...
                {"Completed", "Failed", "Failed"}, {'No Errors', 'first error', 'other error'}));
            testCase.verifyEqual(info.message, 'Pipeline failed: first error (2 folders failed)');
        end

        function completedWithErrorsIsWarning(testCase)
            info = PipelineManager.describeExecutionResult(iResult("completed_with_errors", ...
                {"Completed", "Partially executed"}, {'No Errors', 'one step broke'}));
            testCase.verifyEqual(info.severity, 'warning');
            testCase.verifyEqual(info.message, 'Pipeline finished with errors: one step broke');
        end

        function cancelledAndSkippedAreNotErrors(testCase)
            cancelled = PipelineManager.describeExecutionResult(iResult("cancelled", {}, {}));
            testCase.verifyEqual(cancelled.severity, 'warning');
            testCase.verifyEqual(cancelled.message, 'Pipeline was cancelled.');

            skipped = PipelineManager.describeExecutionResult(iResult("skipped", {}, {}));
            testCase.verifyEqual(skipped.severity, 'info');
            testCase.verifyTrue(contains(skipped.message, 'skipped'));
        end

        function longMultilineMessageIsCollapsedAndTruncated(testCase)
            raw = ['first line' newline repmat('x', 1, 400)];
            info = PipelineManager.describeExecutionResult(iResult("failed", {"Failed"}, {raw}));
            prefix = 'Pipeline failed: ';
            detail = info.message(numel(prefix) + 1:end);
            testCase.verifyFalse(contains(info.message, newline));
            testCase.verifyLessThanOrEqual(numel(detail), 200);
            testCase.verifyTrue(endsWith(detail, '...'));
            testCase.verifyTrue(startsWith(detail, 'first line xxx'));
        end

        function thrownErrorBecomesFailure(testCase)
            ME = MException('PipelineManager:test', 'Preflight failed: missing input.');
            info = PipelineManager.describeExecutionResult(ME);
            testCase.verifyEqual(info.severity, 'error');
            testCase.verifyEqual(info.message, 'Pipeline failed: Preflight failed: missing input.');
        end

        function incompleteResultsNeverThrow(testCase)
            info = PipelineManager.describeExecutionResult(struct());
            testCase.verifyEqual(info.severity, 'info');
            testCase.verifyEqual(info.message, 'Pipeline finished (status: unknown).');

            info = PipelineManager.describeExecutionResult([]);
            testCase.verifyEqual(info.severity, 'info');

            noLog = struct('status', "failed");
            info = PipelineManager.describeExecutionResult(noLog);
            testCase.verifyEqual(info.message, 'Pipeline failed.');

            info = PipelineManager.describeExecutionResult(struct('status', "weird_status"));
            testCase.verifyEqual(info.message, 'Pipeline finished (status: weird_status).');
        end
    end
end

function result = iResult(status, folderStatus, shortMessages)
result = struct('status', status);
n = numel(folderStatus);
result.globalPipeLog = table(repmat({'save'}, n, 1), string(folderStatus(:)), ...
    shortMessages(:), 'VariableNames', {'SaveFolder', 'Status', 'Messages_short'});
end
