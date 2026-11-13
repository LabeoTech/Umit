classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerFileSourceType < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERFILESOURCETYPE FileSource nodes accept only .dat and .umt files.
    %
    %   Companion artifacts such as events.mat or AcqInfos.mat are owned by other
    %   managers and are not data streams, so they must be rejected as sources.

    properties
        ProjectRoot char
        SaveFolder char
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
            sampleFolder = fullfile(testCase.ProjectRoot, 'test', 'Analysis', 'TestingData_with_events');

            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);

            for name = {'AcqInfos.mat', 'events.mat', 'green.dat'}
                copyfile(fullfile(sampleFolder, name{1}), testCase.SaveFolder);
            end

            sourceInfo = loadMetaData(fullfile(testCase.SaveFolder, 'green.dat'));
            appendImportedChannelInfo(testCase.SaveFolder, struct('DatFile', 'green.dat', ...
                'Length', sourceInfo.dimSizes(3), 'FrameRateHz', sourceInfo.frameRateHz));
        end
    end

    methods (Test)
        function datFileIsAccepted(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            [stepTag, srcNodeID] = pm.ensureFileSourceStep('green.dat');
            testCase.verifyEqual(stepTag, 'green.dat');
            testCase.verifyNotEmpty(srcNodeID);
        end

        function eventsMatIsRejectedAsSource(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            testCase.verifyError(@() pm.ensureFileSourceStep('events.mat'), ...
                'PipelineManager:ensureFileSourceNode:UnsupportedFileType');
            testCase.verifyError( ...
                @() pm.ensureFileSourceStep(fullfile(testCase.SaveFolder, 'AcqInfos.mat')), ...
                'PipelineManager:ensureFileSourceNode:UnsupportedFileType');
            testCase.verifyEmpty(pm.nodes);
        end

        function eventsMatIsRejectedAsStepInput(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            testCase.verifyError( ...
                @() pm.addStep('split_data_by_event', 'input', 'events.mat'), ...
                'PipelineManager:ensureFileSourceNode:UnsupportedFileType');
            testCase.verifyEmpty(pm.nodes);
        end
    end
end
