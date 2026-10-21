classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerViewerTempOwnedFile < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERVIEWERTEMPOWNEDFILE viewerTemp promotes function-written files.
    %
    %   split_data_by_event writes dataByEv.dat itself for a .dat input. Under the
    %   'viewerTemp' policy the leaf must end up as a single PMTMP_* file: the
    %   function-named dataByEv.dat may not stay beside it. Under 'saveLeaves' the
    %   function-named file is the permanent output and is kept.

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

            % The sample AcqInfos.mat predates the ImportedChannels registry, which
            % PipelineManager requires.
            sourceInfo = loadMetaData(testCase.SourceFile);
            appendImportedChannelInfo(testCase.SaveFolder, struct('DatFile', 'green.dat', ...
                'Length', sourceInfo.dimSizes(3), 'FrameRateHz', sourceInfo.frameRateHz));
        end
    end

    methods (Test)
        function viewerTempLeavesNoFunctionNamedDatFile(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.addStep('split_data_by_event', 'input', testCase.SourceFile);

            result = pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyFalse(isfile(fullfile(testCase.SaveFolder, 'dataByEv.dat')), ...
                'dataByEv.dat must not stay beside the PMTMP_* file.');

            tempListing = dir(fullfile(testCase.SaveFolder, 'PMTMP_*.dat'));
            testCase.verifyEqual(numel(tempListing), 1);

            manifest = result.outputManifest(result.outputManifest.OutputName == "outData", :);
            testCase.verifyEqual(height(manifest), 1);
            testCase.verifyEqual(char(manifest.ActualFileName(1)), tempListing(1).name);
            testCase.verifyEqual(manifest.ActualPersistence(1), "saved_temp");
            testCase.verifyTrue(manifest.FileExists(1));
            testCase.verifyEqual(cellstr(result.createdFiles.FileName), {tempListing(1).name}, ...
                'The run must report only the PMTMP_* file as created.');
            testCase.verifyTrue(isfile(fullfile(testCase.SaveFolder, 'green.dat')), ...
                'The input file must be untouched.');

            stored = load(fullfile(testCase.SaveFolder, 'dataHistory.mat'), 'dataHistory');
            historyFiles = string({stored.dataHistory.file});
            testCase.verifyFalse(any(strcmpi(historyFiles, 'dataByEv.dat')));
            testCase.verifyTrue(any(strcmpi(historyFiles, tempListing(1).name)));
        end

        function saveLeavesKeepsFunctionNamedDatFile(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.addStep('split_data_by_event', 'input', testCase.SourceFile);

            result = pm.executePipeline('LeafOutputPolicy', 'saveLeaves', 'PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(isfile(fullfile(testCase.SaveFolder, 'dataByEv.dat')));
            testCase.verifyEmpty(dir(fullfile(testCase.SaveFolder, 'PMTMP_*')));
        end
    end
end
