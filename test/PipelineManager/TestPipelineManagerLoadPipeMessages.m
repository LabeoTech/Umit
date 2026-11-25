classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerLoadPipeMessages < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERLOADPIPEMESSAGES What loadPipe says about a file that is not a pipeline.
    %
    %   A corrupt or empty .pipe used to surface MATLAB's own load error ("Not a binary
    %   MAT-file. Try load -ASCII to read as text."), which means nothing for a pipeline
    %   file. loadPipe now raises PipelineManager:loadPipe:InvalidFile with a plain
    %   message and keeps MATLAB's error as the cause. The graph stays untouched.

    properties
        ProjectRoot char
        Folder char
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.Folder);
        end
    end

    methods (Test)
        function aGarbageFileGetsAPlainMessageAndKeepsTheRawErrorAsCause(testCase)
            pipeFile = fullfile(testCase.Folder, 'broken.pipe');
            fid = fopen(pipeFile, 'w');
            fprintf(fid, 'this is not a mat file');
            fclose(fid);
            pm = PipelineManager(testCase.Folder, {'Missing'}, testCase.ProjectRoot);

            try
                pm.loadPipe(pipeFile);
                testCase.verifyFail('Expected the load to be refused.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'PipelineManager:loadPipe:InvalidFile');
                testCase.verifySubstring(ME.message, 'broken.pipe');
                testCase.verifySubstring(ME.message, 'not a valid pipeline file');
                testCase.verifyFalse(contains(ME.message, 'load -ASCII'));
                testCase.verifyNotEmpty(ME.cause);
                testCase.verifyEqual(ME.cause{1}.identifier, 'MATLAB:load:notBinaryFile');
            end
        end

        function anEmptyFileIsRefusedTheSameWay(testCase)
            pipeFile = fullfile(testCase.Folder, 'empty.pipe');
            fclose(fopen(pipeFile, 'w'));
            pm = PipelineManager(testCase.Folder, {'Missing'}, testCase.ProjectRoot);

            testCase.verifyError(@() pm.loadPipe(pipeFile), 'PipelineManager:loadPipe:InvalidFile');
        end

        function aMatFileWithoutAPipelineKeepsItsOwnMessage(testCase)
            pipeFile = fullfile(testCase.Folder, 'other.pipe');
            x = 1;
            save(pipeFile, 'x');
            pm = PipelineManager(testCase.Folder, {'Missing'}, testCase.ProjectRoot);

            try
                pm.loadPipe(pipeFile);
                testCase.verifyFail('Expected the load to be refused.');
            catch ME
                testCase.verifyEqual(ME.identifier, 'PipelineManager:loadPipe:InvalidFile');
                testCase.verifySubstring(ME.message, 'does not contain a valid pipeStruct');
            end
        end

        function aRefusedFileLeavesTheCurrentPipelineAlone(testCase)
            pipeFile = fullfile(testCase.Folder, 'broken.pipe');
            fid = fopen(pipeFile, 'w');
            fprintf(fid, 'garbage');
            fclose(fid);
            pm = PipelineManager(testCase.Folder, {'Missing'}, testCase.ProjectRoot);
            nBefore = numel(pm.nodes);

            testCase.verifyError(@() pm.loadPipe(pipeFile), 'PipelineManager:loadPipe:InvalidFile');

            testCase.verifyEqual(numel(pm.nodes), nBefore);
        end
    end
end
