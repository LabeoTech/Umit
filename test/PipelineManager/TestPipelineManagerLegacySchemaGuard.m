classdef TestPipelineManagerLegacySchemaGuard < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERLEGACYSCHEMAGUARD Validate the legacy-schema processing block.
    %
    %   Covers isLegacySchemaFolder and PipelineManager.executePipeline's
    %   refusal path (task-legacy-readonly-mode.md), using a synthetic
    %   legacy-schema SaveFolder (AcqInfos.mat without ImportedChannels,
    %   mirroring an Astrocyte-era AcqInfos.mat) and a current-schema
    %   SaveFolder (AcqInfos.mat with ImportedChannels) as the two cases.

    properties
        ProjectRoot char
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            thisFile = mfilename('fullpath');
            testFolder = fileparts(thisFile);
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (Test)
        function testLegacySidecarSchemaIsFlagged(testCase)
            saveFolder = testCase.createSaveFolder('legacy');
            [isLegacy, message] = isLegacySchemaFolder(saveFolder);

            testCase.verifyTrue(isLegacy);
            testCase.verifySubstring(message, 'legacy metadata schema');
            testCase.verifySubstring(message, 'Reprocess this dataset from raw/continuous data');
        end

        function testMissingAcqInfosIsFlagged(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(saveFolder);

            [isLegacy, message] = isLegacySchemaFolder(saveFolder);

            testCase.verifyTrue(isLegacy);
            testCase.verifyNotEmpty(message);
        end

        function testCurrentSchemaIsNotFlagged(testCase)
            saveFolder = testCase.createSaveFolder('current');
            [isLegacy, message] = isLegacySchemaFolder(saveFolder);

            testCase.verifyFalse(isLegacy);
            testCase.verifyEmpty(message);
        end

        function testExecutePipelineRefusesLegacySaveFolder(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            saveFolder = testCase.createSaveFolder('legacy');

            pm = PipelineManager(saveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'pipelineTemplate.dat');

            testCase.verifyError(@() pm.executePipeline(), ...
                'PipelineManager:executePipeline:LegacySchema');
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'pipelineTemplate.dat')));
        end

        function testExecutePipelineStillRunsForCurrentSaveFolder(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            saveFolder = testCase.createSaveFolder('current');

            pm = PipelineManager(saveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'pipelineTemplate.dat');
            pm.executePipeline();

            testCase.verifyTrue(isfile(fullfile(saveFolder, 'pipelineTemplate.dat')));
        end

        function testExecutePipelineRefusesAmbiguousNonEmptySaveFolder(testCase)
            parserFolder = createTemporaryTemplateCategory(testCase.ProjectRoot);
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(saveFolder);
            fid = fopen(fullfile(saveFolder, 'input.dat'), 'w');
            testCase.assertNotEqual(fid, -1);
            fwrite(fid, zeros(1, 1, 2, 'single'), 'single');
            fclose(fid);

            pm = PipelineManager(saveFolder, '', testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('funcTemplate', ...
                'input', 'input.dat', ...
                'saveas', 'pipelineTemplate.dat');

            % Without AcqInfos.mat the folder is fresh (only AcqInfos.mat
            % counts, decided 2026-10-02), but a non-importer pipeline cannot
            % initialize it: still refused, nothing written.
            testCase.verifyError(@() pm.executePipeline(), ...
                'PipelineManager:executePipeline:FreshSaveFolderNotInitializable');
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'pipelineTemplate.dat')));
        end

        function testRawFolderCanAlsoBeFreshSaveFolder(testCase)
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshValidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            sharedFolder = fullfile(fixture.Folder, 'RawAndSaveFolder');
            mkdir(sharedFolder);
            createEmptyFile(fullfile(sharedFolder, 'img_00001.bin'));
            createEmptyFile(fullfile(sharedFolder, 'ai_00001.bin'));
            createEmptyFile(fullfile(sharedFolder, 'info.json'));
            mkdir(fullfile(sharedFolder, 'raw_auxiliary'));

            pm = PipelineManager( ...
                sharedFolder, sharedFolder, testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmFreshValidMetadataInitializer');
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(isfile(fullfile( ...
                sharedFolder, 'AcqInfos.mat')));
            testCase.verifyTrue(isfile(fullfile( ...
                sharedFolder, 'imported.dat')));
        end

        function testExistingAcqInfosBlocksInitialization(testCase)
            % Only AcqInfos.mat makes a SaveFolder "not fresh": an existing
            % but invalid one keeps the legacy-schema block for importers.
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshValidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            [saveFolder, rawFolder] = testCase.newFolderPair();
            unexpectedVariable = 1; %#ok<NASGU>
            save(fullfile(saveFolder, 'AcqInfos.mat'), 'unexpectedVariable', '-mat');

            pm = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmFreshValidMetadataInitializer');

            testCase.verifyError(@() pm.executePipeline('PrintSummary', false), ...
                'PipelineManager:executePipeline:LegacySchema');
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'imported.dat')));
        end

        function testLeftoverArtifactsDoNotBlockInitialization(testCase)
            % Leftovers of a failed import (stale, corrupt, or valid .dat
            % files and other artifacts) are ignored without AcqInfos.mat:
            % the importer runs and may overwrite them (decided 2026-10-02).
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshValidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            artifactNames = { ...
                'existing.dat', ...
                'existing.umt', ...
                'existing.roi', ...
                'DataParams.mat', ...
                'events.mat', ...
                'Text_events.mat', ...
                'dataHistory.mat', ...
                'pipeLog.mat', ...
                'LogBook.mat'};

            for iArtifact = 1:numel(artifactNames)
                [saveFolder, rawFolder] = testCase.newFolderPair();
                createEmptyFile(fullfile(saveFolder, artifactNames{iArtifact}));

                pm = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
                pm.b_skipSteps = false;
                pm.addStep('pmFreshValidMetadataInitializer');
                result = pm.executePipeline('PrintSummary', false);

                testCase.verifyEqual(result.status, "completed", ...
                    sprintf('Leftover must not block initialization: %s', artifactNames{iArtifact}));
                testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')), artifactNames{iArtifact});
                testCase.verifyTrue(isfile(fullfile(saveFolder, 'imported.dat')), artifactNames{iArtifact});
            end
        end

        function testStaleChannelFromFailedImportDoesNotBlock(testCase)
            % The DataViewer case: a complete red.dat written by an import
            % that failed before AcqInfos.mat.
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshValidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            [saveFolder, rawFolder] = testCase.newFolderPair();
            saveData(fullfile(saveFolder, 'red.dat'), single(ones(4, 3, 5)), ...
                'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20);

            pm = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmFreshValidMetadataInitializer');
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "completed");
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
        end

        function testInvalidInitializerCannotUnlockCompanion(testCase)
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshInvalidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            rawFolder = fullfile(fixture.Folder, 'RawFolder');
            mkdir(saveFolder);
            mkdir(rawFolder);

            pm = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmFreshInvalidMetadataInitializer');
            pm.addStep('getEvents');
            result = pm.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "failed");
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'events.mat')));
            testCase.verifySubstring( ...
                char(pm.globalPipeLog.Messages_long{1}), ...
                'Fresh SaveFolder initialization did not complete successfully');
        end

        function testFreshFolderRejectsMultipleInitializers(testCase)
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshInvalidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            rawFolder = fullfile(fixture.Folder, 'RawFolder');
            mkdir(saveFolder);
            mkdir(rawFolder);

            pm = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmFreshInvalidMetadataInitializer');
            pm.addStep('pmFreshInvalidMetadataInitializer');

            testCase.verifyError(@() pm.executePipeline('PrintSummary', false), ...
                'PipelineManager:executePipeline:FreshSaveFolderNotInitializable');
            testCase.verifyFalse(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
        end

        function testOlderSavedPipelineResolvesCurrentInitializerRole(testCase)
            parserFolder = createTemporaryFixtureCategory( ...
                testCase.ProjectRoot, 'pmFreshInvalidMetadataInitializer');
            cleanup = onCleanup(@() cleanupTemplateCategory( ...
                parserFolder, testCase.ProjectRoot));

            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            rawFolder = fullfile(fixture.Folder, 'RawFolder');
            pipeFile = fullfile(fixture.Folder, 'olderFreshImporter.pipe');
            mkdir(saveFolder);
            mkdir(rawFolder);

            pm = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
            pm.b_skipSteps = false;
            pm.addStep('pmFreshInvalidMetadataInitializer');
            pm.savePipe(pipeFile);

            loaded = load(pipeFile, 'pipeStruct', '-mat');
            pipeStruct = loaded.pipeStruct;
            streamIdx = find(strcmpi({pipeStruct.nodes.kind}, 'stream'), 1, 'first');
            pipeStruct.nodes(streamIdx).info = rmfield( ...
                pipeStruct.nodes(streamIdx).info, 'freshSaveFolderRole');
            save(pipeFile, 'pipeStruct', '-mat');

            restored = PipelineManager(saveFolder, rawFolder, testCase.ProjectRoot);
            restored.b_skipSteps = false;
            restored.loadPipe(pipeFile);
            result = restored.executePipeline('PrintSummary', false);

            testCase.verifyEqual(result.status, "failed");
            testCase.verifyTrue(isfile(fullfile(saveFolder, 'AcqInfos.mat')));
            testCase.verifySubstring( ...
                char(restored.globalPipeLog.Messages_long{1}), ...
                'Acquisition importer');
        end
    end

    methods (Access = private)
        function [saveFolder, rawFolder] = newFolderPair(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            rawFolder = fullfile(fixture.Folder, 'RawFolder');
            mkdir(saveFolder);
            mkdir(rawFolder);
        end

        function saveFolder = createSaveFolder(testCase, schema)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            saveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(saveFolder);

            inputData = single(reshape(-12:11, [2 3 4]));

            AcqInfoStream = struct();
            AcqInfoStream.Width = 3;
            AcqInfoStream.Height = 2;
            AcqInfoStream.Length = 4;
            AcqInfoStream.FrameRateHz = 20;
            AcqInfoStream.ExposureMsec = 1;

            if strcmpi(schema, 'current')
                AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                    'DatFile', 'input.dat', 'Length', 4, 'FrameRateHz', 20));
            end

            save(fullfile(saveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            saveData(fullfile(saveFolder, 'input.dat'), inputData, 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20);
        end
    end
end

function createEmptyFile(filePath)
%CREATEEMPTYFILE Create one zero-byte fixture file.

fid = fopen(filePath, 'w');
if fid == -1
    error('Umitoolbox:TestPipelineManagerLegacySchemaGuard:FileOpenFailed', ...
        'Could not create fixture file "%s".', filePath);
end
fclose(fid);
end

function parserFolder = createTemporaryTemplateCategory(projectRoot)
%CREATETEMPORARYTEMPLATECATEGORY Copy funcTemplate into a scanned category.
%
%   PipelineManager.createFcnList only scans Analysis/*/*.m, so the
%   root-level funcTemplate.m must be copied into a category subfolder to
%   be discoverable via addStep (mirrors
%   TestPipelineManagerFunctionTemplate's identical helper).

parserFolder = createTemporaryFixtureCategory(projectRoot, 'funcTemplate');
end

function parserFolder = createTemporaryFixtureCategory(projectRoot, fixtureName)
%CREATETEMPORARYFIXTURECATEGORY Copy one fixture into a scanned category.

analysisFolder = fullfile(projectRoot, 'Analysis');
parserFolder = tempname(analysisFolder);
mkdir(parserFolder);

if strcmp(fixtureName, 'funcTemplate')
    sourceFile = fullfile(analysisFolder, 'funcTemplate.m');
else
    sourceFile = fullfile(projectRoot, 'test', 'PipelineManager', ...
        'fixtures', [fixtureName '.m']);
end
copyfile(sourceFile, fullfile(parserFolder, [fixtureName '.m']));

addpath(parserFolder, '-begin');
clear(fixtureName);
rehash;
end

function cleanupTemplateCategory(parserFolder, projectRoot)
%CLEANUPTEMPLATECATEGORY Remove the temporary scanned category.

analysisFolder = fullfile(projectRoot, 'Analysis');
expectedPrefix = [analysisFolder filesep];
if ~startsWith(parserFolder, expectedPrefix, 'IgnoreCase', true)
    error('umIToolbox:TestPipelineManagerLegacySchemaGuard:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end

fixtureFiles = dir(fullfile(parserFolder, '*.m'));
for iFile = 1:numel(fixtureFiles)
    [~, fixtureName] = fileparts(fixtureFiles(iFile).name);
    clear(fixtureName);
end
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
rehash;
end
