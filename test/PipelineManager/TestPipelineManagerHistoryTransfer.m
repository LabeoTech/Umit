classdef TestPipelineManagerHistoryTransfer < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERHISTORYTRANSFER Validate PipelineManager.transferDataHistoryEntry.
    %
    %   dataHistory.mat is keyed by file name, so provenance must follow a
    %   file that DataViewer's Save as... copies or renames.

    methods (Test)
        function testMoveRekeysEntryAndKeepsProvenance(testCase)
            folder = testCase.newFolder();
            testCase.writeHistory(folder, {'PMTMP_a.dat', 'keep.dat'});

            PipelineManager.transferDataHistoryEntry( ...
                folder, 'PMTMP_a.dat', 'final.dat', 'move');

            files = testCase.readHistoryFiles(folder);
            testCase.verifyEqual(files, {'final.dat'; 'keep.dat'});
            history = load(fullfile(folder, 'dataHistory.mat'));
            testCase.verifyEqual(history.dataHistory(1).info, 'info-PMTMP_a.dat');
        end

        function testCopyKeepsSourceEntry(testCase)
            folder = testCase.newFolder();
            testCase.writeHistory(folder, {'a.dat'});

            PipelineManager.transferDataHistoryEntry( ...
                folder, 'a.dat', 'b.dat', 'copy');

            testCase.verifyEqual(testCase.readHistoryFiles(folder), ...
                {'a.dat'; 'b.dat'});
        end

        function testStaleTargetEntryIsReplaced(testCase)
            folder = testCase.newFolder();
            testCase.writeHistory(folder, {'a.dat', 'b.dat'});

            PipelineManager.transferDataHistoryEntry( ...
                folder, 'a.dat', 'b.dat', 'move');

            files = testCase.readHistoryFiles(folder);
            testCase.verifyEqual(files, {'b.dat'});
            history = load(fullfile(folder, 'dataHistory.mat'));
            testCase.verifyEqual(history.dataHistory(1).info, 'info-a.dat');
        end

        function testMissingHistoryOrEntryIsNoOp(testCase)
            folder = testCase.newFolder();
            PipelineManager.transferDataHistoryEntry( ...
                folder, 'a.dat', 'b.dat', 'move');
            testCase.verifyFalse(isfile(fullfile(folder, 'dataHistory.mat')));

            testCase.writeHistory(folder, {'x.dat'});
            PipelineManager.transferDataHistoryEntry( ...
                folder, 'a.dat', 'b.dat', 'copy');
            testCase.verifyEqual(testCase.readHistoryFiles(folder), {'x.dat'});
        end

        function testInvalidModeErrors(testCase)
            folder = testCase.newFolder();
            testCase.verifyError(@() PipelineManager.transferDataHistoryEntry( ...
                folder, 'a.dat', 'b.dat', 'link'), ...
                'Umitoolbox:PipelineManager:invalidHistoryTransferMode');
        end
    end

    methods (Access = private)
        function folder = newFolder(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            folder = fixture.Folder;
        end

        function writeHistory(~, folder, fileNames)
            dataHistory = struct('file', fileNames, ...
                'info', cellfun(@(f) ['info-' f], fileNames, 'UniformOutput', false)); %#ok<NASGU>
            save(fullfile(folder, 'dataHistory.mat'), 'dataHistory');
        end

        function files = readHistoryFiles(~, folder)
            history = load(fullfile(folder, 'dataHistory.mat'));
            files = arrayfun(@(e) char(string(e.file)), history.dataHistory(:), ...
                'UniformOutput', false);
        end
    end
end
