classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestEnsureUMITPath < matlab.unittest.TestCase
    %TESTENSUREUMITPATH Which checkout folders ensureUMITPath keeps on the MATLAB path.
    %
    %   The saved path of a developer's MATLAB held hundreds of .git and docs
    %   entries (printed as "Name is nonexistent or not a directory" at every
    %   startup) and local folders such as _tasks. ensureUMITPath now adds only
    %   production folders (plus the _tools and _mlapp_tools development tools)
    %   and removes every excluded entry it finds. All tests use Persist = false
    %   and restore the path afterwards.

    properties
        Root char
        SavedPath char
    end

    methods (TestMethodSetup)
        function rememberPath(testCase)
            testCase.Root = fileparts(fileparts(mfilename('fullpath')));
            testCase.SavedPath = path;
            testCase.addTeardown(@() path(testCase.SavedPath));
        end
    end

    methods (Test)
        function excludedFoldersAreRemovedFromThePath(testCase)
            excluded = testCase.existingExcludedFolders();
            testCase.assertNotEmpty(excluded, 'The checkout has no excluded folders to test with.');
            addpath(excluded{:});

            ensureUMITPath('Persist', false);

            onPath = strsplit(path, pathsep);
            for k = 1:numel(excluded)
                testCase.verifyFalse(any(strcmpi(onPath, excluded{k})), ...
                    sprintf('%s must not stay on the path.', excluded{k}));
            end
        end

        function noGitDocsHiddenOrLocalFolderStaysOnThePath(testCase)
            excluded = testCase.existingExcludedFolders();
            addpath(excluded{:});

            ensureUMITPath('Persist', false);

            entries = strsplit(path, pathsep);
            inRepo = entries(startsWith(entries, [testCase.Root filesep], 'IgnoreCase', true));
            relative = extractAfter(inRepo, [testCase.Root filesep]);
            for k = 1:numel(relative)
                segments = strsplit(relative{k}, {'/', '\'});
                testCase.verifyFalse(any(startsWith(segments, '.')), relative{k});
                testCase.verifyFalse(any(strcmpi(segments{1}, {'docs', 'claude-context', 'test'})), relative{k});
                if startsWith(segments{1}, '_')
                    testCase.verifyTrue(any(strcmpi(segments{1}, {'_tools', '_mlapp_tools'})), relative{k});
                end
            end
        end

        function productionFoldersAndDevelopmentToolsResolve(testCase)
            ensureUMITPath('Persist', false);

            for name = {'PipelineManager', 'loadData', 'normalizeBSLN', 'validateTaskChanges', ...
                    'exportMlappToM'}
                located = which(name{1});
                testCase.verifyNotEmpty(located, sprintf('%s must resolve.', name{1}));
                testCase.verifyTrue(startsWith(located, testCase.Root, 'IgnoreCase', true), ...
                    sprintf('%s resolves outside the checkout: %s', name{1}, located));
            end
        end

        function testMocksDoNotShadowMatlabBuiltIns(testCase)
            ensureUMITPath('Persist', false);

            memoryLocation = which('memory');
            testCase.verifyFalse(contains(memoryLocation, [filesep 'test' filesep], 'IgnoreCase', true), ...
                memoryLocation);
        end

        function aCleanPathIsLeftAlone(testCase)
            ensureUMITPath('Persist', false);
            pathAfterFirstCall = path;

            result = ensureUMITPath('Persist', false);

            testCase.verifyFalse(result.Changed);
            testCase.verifyEmpty(result.MissingFolders);
            testCase.verifyEmpty(result.RemovedFolders);
            testCase.verifyFalse(result.SaveAttempted);
            testCase.verifyEqual(path, pathAfterFirstCall);
        end

        function aMissingProductionFolderIsAddedBack(testCase)
            ensureUMITPath('Persist', false);
            rmpath(fullfile(testCase.Root, 'subFunc'));

            result = ensureUMITPath('Persist', false);

            testCase.verifyTrue(result.Changed);
            testCase.verifyTrue(any(strcmpi(result.MissingFolders, fullfile(testCase.Root, 'subFunc'))));
            testCase.verifyTrue(any(strcmpi(strsplit(path, pathsep), fullfile(testCase.Root, 'subFunc'))));
        end

        function theSavedPathListsExcludedAndVanishedCheckoutEntries(testCase)
            % MATLAB drops vanished folders from the session path at startup, so only the
            % saved path file still names them (and prints the startup warnings).
            root = testCase.Root;
            entries = { ...
                fullfile(root, '.git', 'objects', '08'), ...        % excluded (and gone)
                fullfile(root, 'docs', 'html'), ...                  % excluded (and gone)
                fullfile(root, 'GUI'), ...                           % production, exists
                fullfile(root, 'GUI', 'noLongerThere'), ...          % vanished
                fullfile(root, '_tools'), ...                        % development tool, exists
                fullfile(root, '_tasks'), ...                        % excluded
                fullfile(matlabroot, 'toolbox', 'matlab', 'general')}; % not the checkout
            pathDefFile = testCase.writePathDef(entries);

            result = ensureUMITPath('Persist', false, 'PathDefFile', pathDefFile);

            testCase.verifyEqual(result.StaleSavedEntries, ...
                {entries{1}, entries{2}, entries{4}, entries{6}});
            testCase.verifyFalse(result.SaveAttempted);
        end

        function aCleanSavedPathReportsNothingStale(testCase)
            root = testCase.Root;
            pathDefFile = testCase.writePathDef({root, fullfile(root, 'GUI'), ...
                fullfile(root, 'classes'), fullfile(root, '_mlapp_tools')});

            result = ensureUMITPath('Persist', false, 'PathDefFile', pathDefFile);

            testCase.verifyEmpty(result.StaleSavedEntries);
        end

        function persistFalseNeverSavesThePath(testCase)
            excluded = testCase.existingExcludedFolders();
            addpath(excluded{:});

            result = ensureUMITPath('Persist', false);

            testCase.verifyTrue(result.Changed);
            testCase.verifyFalse(result.SaveAttempted);
            testCase.verifyTrue(isnan(result.SaveStatus));
        end
    end

    methods (Access = private)
        function file = writePathDef(testCase, entries)
            %WRITEPATHDEF Write a pathdef.m look-alike listing the given path entries.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(folder.Folder, 'pathdef.m');
            fid = fopen(file, 'w');
            assert(fid ~= -1, 'Could not create the pathdef test file.');
            fprintf(fid, 'function p = pathdef\n%%%% BEGIN ENTRIES %%%%\np = [...\n');
            for k = 1:numel(entries)
                fprintf(fid, '''%s;'', ...\n', entries{k});
            end
            fprintf(fid, '];\n%%%% END ENTRIES %%%%\nend\n');
            fclose(fid);
        end

        function folders = existingExcludedFolders(testCase)
            %EXISTINGEXCLUDEDFOLDERS Excluded folders of this checkout that exist.
            candidates = {'.github', '.claude', 'docs', 'claude-context', '_tasks', '_logs', ...
                fullfile('test', 'subFunc')};
            folders = cell(1, 0);
            for k = 1:numel(candidates)
                candidate = fullfile(testCase.Root, candidates{k});
                if isfolder(candidate)
                    folders{end + 1} = candidate; %#ok<AGROW>
                end
            end
        end
    end
end
