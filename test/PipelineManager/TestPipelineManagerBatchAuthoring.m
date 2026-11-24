classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestPipelineManagerBatchAuthoring < matlab.unittest.TestCase
    %TESTPIPELINEMANAGERBATCHAUTHORING Adding steps to a pipeline that spans several folders.
    %
    %   addStep checks a file source only in the first SaveFolder, so the window needs
    %   findFoldersMissingFile to warn about the other folders before a run. A refused add
    %   (duplicate function in a branch) must return an empty tag and leave the graph as it
    %   was, because the window reports "not added" for that case, not "cancelled".
    %
    %   Data is 1:24 as 2 x 3 x 4 (Y, X, T): non-square.

    properties
        ProjectRoot char
        Folders cell
    end

    methods (TestClassSetup)
        function resolveProjectRoot(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            testCase.ProjectRoot = extractBefore(testFolder, [filesep 'test']);
        end
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folders = cell(1, 3);
            for iFolder = 1:3
                folder = fullfile(fixture.Folder, sprintf('SaveFolder%d', iFolder));
                mkdir(folder);
                testCase.Folders{iFolder} = folder;

                AcqInfoStream = struct('Width', 3, 'Height', 2, 'Length', 4, ...
                    'FrameRateHz', 10, 'ExposureMsec', 1);
                AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, ...
                    struct('DatFile', 'green.dat', 'Length', 4, 'FrameRateHz', 10));
                save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
                if iFolder ~= 3
                    saveData(fullfile(folder, 'green.dat'), single(reshape(1:24, [2 3 4])), ...
                        'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 10);
                end
            end

            parserFolder = iCreateFixtureCategory(testCase.ProjectRoot);
            testCase.addTeardown(@() iRemoveFixtureCategory(parserFolder, testCase.ProjectRoot));
        end
    end

    methods (Test)
        function reportsTheFoldersThatLackTheSourceFile(testCase)
            pm = PipelineManager(testCase.Folders, {'Missing', 'Missing', 'Missing'}, testCase.ProjectRoot);

            missing = pm.findFoldersMissingFile('green.dat');

            testCase.verifyEqual(numel(missing), 1);
            testCase.verifyTrue(strcmpi(missing{1}, iNormalized(testCase.Folders{3})));
        end

        function usesOnlyTheFileNameOfAFullPath(testCase)
            pm = PipelineManager(testCase.Folders, {'Missing', 'Missing', 'Missing'}, testCase.ProjectRoot);

            missing = pm.findFoldersMissingFile(fullfile(testCase.Folders{1}, 'green.dat'));

            testCase.verifyEqual(numel(missing), 1);
        end

        function returnsNothingWhenEveryFolderHasTheFile(testCase)
            pm = PipelineManager(testCase.Folders(1:2), {'Missing', 'Missing'}, testCase.ProjectRoot);

            testCase.verifyEmpty(pm.findFoldersMissingFile('green.dat'));
        end

        function reportsEveryFolderForAFileNoneHas(testCase)
            pm = PipelineManager(testCase.Folders, {'Missing', 'Missing', 'Missing'}, testCase.ProjectRoot);

            testCase.verifyEqual(numel(pm.findFoldersMissingFile('nothere.dat')), 3);
        end

        function addingTheSourceStillChecksOnlyTheFirstFolder(testCase)
            pm = PipelineManager(testCase.Folders, {'Missing', 'Missing', 'Missing'}, testCase.ProjectRoot);

            tag = pm.addStep('pmSkipStepA', 'input', 'green.dat');

            testCase.verifyNotEmpty(tag, ...
                'The add succeeds, which is why the window must check the other folders itself.');
        end

        function refusedDuplicateReturnsEmptyTagAndKeepsTheGraph(testCase)
            pm = PipelineManager(testCase.Folders(1:2), {'Missing', 'Missing'}, testCase.ProjectRoot);
            tagA = pm.addStep('pmSkipStepA', 'input', 'green.dat');
            nNodes = numel(pm.nodes);
            nConnections = numel(pm.connections);

            lastwarn('');
            tagAgain = pm.addStep('pmSkipStepA', 'input', tagA);
            [~, warnID] = lastwarn;

            testCase.verifyEmpty(tagAgain);
            testCase.verifyEqual(warnID, 'PipelineManager:addStep:DuplicateFunctionInBranch');
            testCase.verifyEqual(numel(pm.nodes), nNodes);
            testCase.verifyEqual(numel(pm.connections), nConnections);
        end
    end
end

function folder = iNormalized(folder)
[~, info] = fileattrib(folder);
folder = info.Name;
end

function parserFolder = iCreateFixtureCategory(projectRoot)
% PipelineManager discovers functions under Analysis/: copy the fixtures into a
% temporary category there (createAnalysisFixtureCategory).
fixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
names = {'pmSkipStepA', 'pmSkipStepB'};
files = struct('name', names, 'source', fullfile(fixtureFolder, strcat(names, '.m')));
parserFolder = createAnalysisFixtureCategory(projectRoot, files);
end

function iRemoveFixtureCategory(parserFolder, projectRoot)
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~startsWith(parserFolder, [analysisFolder filesep], 'IgnoreCase', true)
    error('Umitoolbox:TestPipelineManagerBatchAuthoring:UnsafeCleanupPath', ...
        'Refusing to remove a folder outside Analysis.');
end
clear pmSkipStepA pmSkipStepB
if any(strcmp(strsplit(path, pathsep), parserFolder))
    rmpath(parserFolder);
end
if isfolder(parserFolder)
    rmdir(parserFolder, 's');
end
end
