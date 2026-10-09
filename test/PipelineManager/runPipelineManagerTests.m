function results = runPipelineManagerTests(saveFolder, rawFolder)
%RUNPIPELINEMANAGERTESTS Run the PipelineManager unit test suite.
%
%   results = runPipelineManagerTests()
%   results = runPipelineManagerTests(saveFolder, rawFolder)
%
%   If saveFolder/rawFolder are provided, the test suite creates an isolated
%   temporary test-run folder inside saveFolder. Each test then creates its own
%   SaveFolder/RawFolder pair inside that run folder. The complete test-run
%   folder is deleted automatically when this function exits.
%
%   Output:
%       results - matlab.unittest.TestResult array.

    if nargin < 1 || isempty(saveFolder)
        saveFolder = '';
    end
    if nargin < 2 || isempty(rawFolder)
        rawFolder = saveFolder;
    end

    thisFile = mfilename('fullpath');
    testFolder = fileparts(thisFile);

    assert(isfolder(testFolder), ...
        'runPipelineManagerTests:MissingTestFolder', ...
        'Test folder not found: %s', testFolder);

    if ~isempty(saveFolder)
        assert(isfolder(saveFolder), ...
            'runPipelineManagerTests:InvalidSaveFolder', ...
            'Provided saveFolder does not exist: %s', saveFolder);
    end

    if ~isempty(rawFolder)
        assert(isfolder(rawFolder), ...
            'runPipelineManagerTests:InvalidRawFolder', ...
            'Provided rawFolder does not exist: %s', rawFolder);
    end

    testRunRoot = '';

    if ~isempty(saveFolder)
        timeTag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
        testRunRoot = fullfile(char(string(saveFolder)), ['PipelineManagerTestRun_' timeTag]);
        mkdir(testRunRoot);
        cleanupRunRoot = onCleanup(@() removeFolderIfExists(testRunRoot));
    end

    cfg = struct();
    cfg.useExternalFolders = ~isempty(saveFolder);

    if cfg.useExternalFolders
        cfg.saveFolder = testRunRoot;
    else
        cfg.saveFolder = '';
    end

    % Stored for traceability. The current synthetic unit tests create their own
    % per-test RawFolder, but preserving this field keeps the config extensible.
    cfg.rawFolder = char(string(rawFolder));

    setappdata(0, 'PipelineManagerTestConfig', cfg);
    cleanupCfg = onCleanup(@() rmappdataIfExists('PipelineManagerTestConfig'));

    suite = testsuite(testFolder, 'IncludeSubfolders', true);
    suiteNames = string({suite.Name});
    keepMask = contains(suiteNames, "TestPipelineManager");

    assert(any(keepMask), ...
        'runPipelineManagerTests:MissingTestClass', ...
        'Could not find TestPipelineManager tests under: %s', testFolder);

    results = run(suite(keepMask));

    disp(table(results))
end


function rmappdataIfExists(name)
%RMAPPDATAIFEXISTS Remove appdata only when it exists.

    if isappdata(0, name)
        rmappdata(0, name);
    end
end


function removeFolderIfExists(folderPath)
%REMOVEFOLDERIFEXISTS Delete one temporary folder tree if it still exists.

    folderPath = char(string(folderPath));

    if isempty(folderPath) || ~isfolder(folderPath)
        return
    end

    try
        rmdir(folderPath, 's');
    catch ME
        warning('runPipelineManagerTests:CleanupFailed', ...
            'Failed to delete temporary test folder "%s". Error: %s', ...
            folderPath, ME.message);
    end
end
