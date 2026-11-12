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

        function intermediateFunctionWrittenFilesAreRemoved(testCase)
            % split_data_by_event and normalizeBSLN both write default-named .dat
            % files. The split output only feeds normalizeBSLN, so dataByEv.dat
            % must not survive the run; the leaf is the only product.
            for policy = {'saveLeaves', 'viewerTemp'}
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
                testCase.verifyFalse(isfile(fullfile(folder, 'dataByEv.dat')), ...
                    [policy{1} ': intermediate dataByEv.dat must be deleted.']);
                testCase.verifyEmpty(dir(fullfile(folder, 'PMTMP_*split_data_by_event*')), ...
                    [policy{1} ': renamed intermediate must be deleted after its consumer ran.']);
                if strcmp(policy{1}, 'saveLeaves')
                    testCase.verifyTrue(isfile(fullfile(folder, 'normBSLN.dat')));
                else
                    testCase.verifyFalse(isfile(fullfile(folder, 'normBSLN.dat')));
                    testCase.verifyNumElements(dir(fullfile(folder, 'PMTMP_*.dat')), 1);
                end
            end
        end

        function explicitSaveOfIntermediateKeepsItPermanent(testCase)
            % The user asked to keep the intermediate: saving it under its own
            % default name, or under a custom name, must leave exactly one
            % permanent copy and no PMTMP leftovers.
            cases = {'dataByEv.dat', 'mySplit.dat'};
            for iCase = 1:numel(cases)
                saveName = cases{iCase};
                folder = fullfile(testCase.SaveFolder, sprintf('explicit_%d', iCase));
                mkdir(folder);
                for name = {'AcqInfos.mat', 'events.mat', 'green.dat'}
                    copyfile(fullfile(testCase.SaveFolder, name{1}), folder);
                end
                sourceInfo = loadMetaData(fullfile(folder, 'green.dat'));
                appendImportedChannelInfo(folder, struct('DatFile', 'green.dat', ...
                    'Length', sourceInfo.dimSizes(3), 'FrameRateHz', sourceInfo.frameRateHz));

                pm = PipelineManager(folder, '', testCase.ProjectRoot);
                idSplit = pm.addStep('split_data_by_event', 'input', ...
                    fullfile(folder, 'green.dat'), 'saveas', saveName);
                pm.addStep('normalizeBSLN', 'input', idSplit, 'saveas', 'final.dat');

                result = pm.executePipeline('LeafOutputPolicy', 'viewerTemp', 'PrintSummary', false);

                testCase.verifyEqual(result.status, "completed", saveName);
                testCase.verifyTrue(isfile(fullfile(folder, saveName)), ...
                    [saveName ': explicitly saved intermediate must persist.']);
                testCase.verifyTrue(isfile(fullfile(folder, 'final.dat')));
                testCase.verifyEmpty(dir(fullfile(folder, 'PMTMP_*')), ...
                    [saveName ': no PMTMP leftovers.']);
                if ~strcmp(saveName, 'dataByEv.dat')
                    testCase.verifyFalse(isfile(fullfile(folder, 'dataByEv.dat')), ...
                        [saveName ': function-named copy must not stay beside the saved file.']);
                end
            end
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
