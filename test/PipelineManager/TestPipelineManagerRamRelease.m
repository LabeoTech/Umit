classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerRamRelease < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERRAMRELEASE A finished run must not keep arrays in RAM.
    %
    %   DataViewer keeps its PipelineManager alive after a run. A leaf output
    %   has no downstream consumer, so decrementConsumers never released it and
    %   its ramValue stayed in the run registry until the next run, starving
    %   later data loads (DatImageSource:InsufficientRAM). Every record may keep
    %   its file backing, but none may keep an array.

    properties
        ProjectRoot char
        SaveFolder char
        SourceFile char
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            projectRoot = extractBefore(testFolder, [filesep 'test']);
            testCase.ProjectRoot = projectRoot;
            sampleFolder = fullfile(projectRoot, 'test', 'Analysis', 'TestingData_with_events');

            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);

            for name = {'AcqInfos.mat', 'events.mat', 'green.dat'}
                copyfile(fullfile(sampleFolder, name{1}), testCase.SaveFolder);
            end
            testCase.SourceFile = fullfile(testCase.SaveFolder, 'green.dat');

            sourceInfo = loadMetaData(testCase.SourceFile);
            appendImportedChannelInfo(testCase.SaveFolder, struct('DatFile', 'green.dat', ...
                'Length', sourceInfo.dimSizes(3), 'FrameRateHz', sourceInfo.frameRateHz));
        end
    end

    methods (Test)
        function leafArrayIsReleasedAfterRun(testCase)
            for policy = {'viewerTemp', 'saveLeaves'}
                folder = fullfile(testCase.SaveFolder, ['run_' policy{1}]);
                mkdir(folder);
                for name = {'AcqInfos.mat', 'events.mat', 'green.dat'}
                    copyfile(fullfile(testCase.SaveFolder, name{1}), folder);
                end
                sourceInfo = loadMetaData(fullfile(folder, 'green.dat'));
                appendImportedChannelInfo(folder, struct('DatFile', 'green.dat', ...
                    'Length', sourceInfo.dimSizes(3), 'FrameRateHz', sourceInfo.frameRateHz));

                pm = PipelineManager(folder, '', testCase.ProjectRoot);
                idSplit = pm.addStep('split_data_by_event', 'input', fullfile(folder, 'green.dat'));
                pm.addStep('normalizeBSLN', 'input', idSplit);

                result = pm.executePipeline('LeafOutputPolicy', policy{1}, 'PrintSummary', false);

                testCase.verifyEqual(result.status, "completed", policy{1});
                [nStored, nRamBytes] = TestPipelineManagerRamRelease.inspectDataStore(pm);
                testCase.verifyGreaterThan(nStored, 0, ...
                    [policy{1} ': file-backed records must remain after the run.']);
                testCase.verifyEqual(nRamBytes, 0, ...
                    [policy{1} ': no ramValue may stay in the run registry.']);
            end
        end
    end

    methods (Static, Access = private)
        function [nStored, nRamBytes] = inspectDataStore(pm)
            %INSPECTDATASTORE Count dataStore records and the bytes of arrays they hold.
            warnState = warning('off', 'MATLAB:structOnObject');
            restoreWarning = onCleanup(@() warning(warnState));
            snapshot = struct(pm);
            store = snapshot.dataStore;
            keys = store.keys;
            nStored = numel(keys);
            nRamBytes = 0;
            for iKey = 1:nStored
                rec = store(keys{iKey});
                if isfield(rec, 'ramValue') && ~isempty(rec.ramValue)
                    value = rec.ramValue; %#ok<NASGU>
                    info = whos('value');
                    nRamBytes = nRamBytes + info.bytes;
                end
            end
        end
    end
end
