classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestGetEvents < matlab.unittest.TestCase
    %TESTGETEVENTS PipelineManager integration coverage for getEvents.

    properties
        TempFolder
        RawFolder
        RigRoot
    end

    methods (TestMethodSetup)
        function setup(testCase)
            projectRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            testCase.RawFolder = fullfile(projectRoot, 'test', 'Analysis', ...
                'TestingData_with_events');
            testCase.TempFolder = fullfile(tempdir, ...
                ['TestGetEvents_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempFolder);
            % ImagesClassification resolves the Default Rig, so make a
            % throwaway rig the Default one in this class's temporary root.
            rigStore = UMITRigStore.create(struct('rigID', ...
                ['TestGetEvents_' strrep(char(java.util.UUID.randomUUID()), '-', '')]));
            UMITRigStore.setDefaultRig(rigStore.getRigInfo().uuid);
            testCase.RigRoot = rigStore.RigRoot;
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            schema = getUMITRigSchema();
            defaultFile = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.store.internalFolder, schema.store.defaultFile);
            if isfile(defaultFile)
                delete(defaultFile);
            end
            if ~isempty(testCase.RigRoot) && isfolder(testCase.RigRoot)
                rmdir(testCase.RigRoot, 's');
            end
            if isfolder(testCase.TempFolder)
                rmdir(testCase.TempFolder, 's');
            end
        end
    end

    methods (Test)
        function testPipelineManagerExecutesEventImporter(testCase)
            pm = buildPMForScenario(testCase.TempFolder, ...
                'ImagesClassification', 'auto', ...
                'RawFolder', testCase.RawFolder);
            pm.addStep('getEvents');
            pm.executePipeline('PrintSummary', false);

            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'events.mat')));
        end

        function testImportOverLegacyAcqInfosIsBackedUpAndReplaced(testCase)
            % Like DataViewer's import: ImagesClassification plus getEvents
            % over a SaveFolder whose AcqInfos.mat is old-schema. The fresh
            % import zips the old metadata, then replaces it (DFR-20261002-001).
            AcqInfoStream = struct('Width', 1, 'Height', 1);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream', '-mat');
            testCase.assertTrue(isLegacySchemaFolder(testCase.TempFolder));

            pm = buildPMForScenario(testCase.TempFolder, ...
                'ImagesClassification', 'auto', ...
                'RawFolder', testCase.RawFolder);
            pm.addStep('getEvents');
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyFalse(isLegacySchemaFolder(testCase.TempFolder));
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, 'events.mat')));
            zips = dir(fullfile(testCase.TempFolder, 'bkp_*.zip'));
            testCase.assertNumElements(zips, 1, 'The fresh import must zip the old files.');
            extractFolder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            unzip(fullfile(zips(1).folder, zips(1).name), extractFolder.Folder);
            old = load(fullfile(extractFolder.Folder, 'AcqInfos.mat'), 'AcqInfoStream');
            testCase.verifyEqual(old.AcqInfoStream, AcqInfoStream);
        end

        function testPipelineManagerRejectsEventImporterWithoutAcquisition(testCase)
            pm = buildPMForScenario(testCase.TempFolder, 'getEvents', 'auto', ...
                'RawFolder', testCase.RawFolder);

            testCase.verifyError(@() pm.executePipeline('PrintSummary', false), ...
                'PipelineManager:executePipeline:FreshSaveFolderNotInitializable');
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'events.mat')));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));
        end

        function testPipelineInfoDeclaresCompanionRole(testCase)
            info = getEvents('pipelineInfo');
            testCase.verifyEqual(info.freshSaveFolderRole, ...
                'acquisition-companion');
        end
    end
end
