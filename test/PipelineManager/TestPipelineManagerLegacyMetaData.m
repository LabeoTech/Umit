classdef TestPipelineManagerLegacyMetaData < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERLEGACYMETADATA PipelineManager supports only pipelineInfo functions.
    %
    %   .dat header Phase 8a removed legacy pipeline support (Bruno,
    %   2026-10-05). This suite replaces the Phase 3c-2/7a checks of the
    %   legacy metaData input:
    %     - a function that does not answer fcn('pipelineInfo') is skipped
    %       with PipelineManager:createFcnList:NoPipelineInfo;
    %     - a pipelineInfo that declares a metaData input is rejected with
    %       PipelineManager:createFcnList:MetaDataInputUnsupported;
    %     - generated scripts no longer contain localGetMetaData.
    %   Fixtures are copied into a temporary Analysis/ category, as
    %   PipelineManager discovers functions there.

    properties
        ProjectRoot char
        SaveFolder char
        FixtureCategory char = ''
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createSaveFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.SaveFolder = fullfile(fixture.Folder, 'SaveFolder');
            mkdir(testCase.SaveFolder);
            AcqInfoStream = struct('Width', 3, 'Height', 2, 'FrameRateHz', 20);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
                'DatFile', 'input.dat', 'Length', 4, 'FrameRateHz', 20));
            save(fullfile(testCase.SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
            saveData(fullfile(testCase.SaveFolder, 'input.dat'), ...
                single(reshape(-12:11, [2 3 4])), 'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20);
            testCase.addTeardown(@() iRemoveFixtureCategory(testCase.FixtureCategory, testCase.ProjectRoot));
        end
    end

    methods (Test)
        function functionWithoutPipelineInfoIsSkipped(testCase)
            testCase.FixtureCategory = iCreateFixtureCategory(testCase.ProjectRoot, ...
                'pmNoPipelineInfoLegacy');
            pm = testCase.verifyWarning(@() PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot), ...
                'PipelineManager:createFcnList:NoPipelineInfo');
            testCase.verifyFalse(any(strcmp({pm.funcList.name}, 'pmNoPipelineInfoLegacy')), ...
                'a function without pipelineInfo must not be listed');
            testCase.verifyTrue(any(strcmp({pm.funcList.name}, 'spatialGaussFilt')), ...
                'pipelineInfo functions are still listed');
        end

        function metaDataInputIsRejected(testCase)
            testCase.FixtureCategory = iCreateFixtureCategory(testCase.ProjectRoot, ...
                'pmLegacyMetaDataEcho');
            testCase.verifyError(@() PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot), ...
                'PipelineManager:createFcnList:MetaDataInputUnsupported');
        end

        function generatedScriptHasNoMetaDataHelper(testCase)
            pm = PipelineManager(testCase.SaveFolder, '', testCase.ProjectRoot);
            pm.addStep('spatialGaussFilt', 'input', 'input.dat');
            scriptFile = fullfile(testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder, 'pipelineScript.m');
            pm.generateScript(scriptFile);
            text = fileread(scriptFile);
            testCase.verifyFalse(contains(text, 'localGetMetaData'));
            testCase.verifyFalse(contains(text, 'legacyMetaDataFromSourceInfo'));
        end
    end
end

% =========================================================================
function categoryFolder = iCreateFixtureCategory(projectRoot, fixtureName)
% Copy one fixture into a temporary Analysis/ category.
categoryFolder = tempname(fullfile(projectRoot, 'Analysis'));
mkdir(categoryFolder);
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
copyfile(fullfile(fixtureFolder, [fixtureName '.m']), fullfile(categoryFolder, [fixtureName '.m']));
addpath(categoryFolder, '-begin');
clear(fixtureName);
rehash;
end

function iRemoveFixtureCategory(categoryFolder, projectRoot)
if isempty(categoryFolder)
    return
end
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(categoryFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerLegacyMetaData:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmLegacyMetaDataEcho pmNoPipelineInfoLegacy
if any(strcmp(strsplit(path, pathsep), categoryFolder))
    rmpath(categoryFolder);
end
if isfolder(categoryFolder)
    rmdir(categoryFolder, 's');
end
end
