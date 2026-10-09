function results = run_TestApplyAggregateFunction(sampleDataFolder)
%RUN_TESTAPPLYAGGREGATEFUNCTION Run the apply_aggregate_function unit tests.
%
%   results = run_TestApplyAggregateFunction()
%   results = run_TestApplyAggregateFunction(sampleDataFolder)
%
%   This wrapper configures the fixture folder used by the
%   TestApplyAggregateFunction test class, runs the full suite, and then
%   clears the temporary configuration.
%
%   Inputs:
%       sampleDataFolder - Optional path to the folder containing the test
%                          fixture files:
%                              green.dat
%                              AcqInfos.mat
%
%                          If omitted, the default location is:
%                              <projectRoot>/test/TestData/ApplyAggregateFunction
%
%   Output:
%       results          - matlab.unittest.TestResult array returned by runtests.
%
%   Behavior:
%       - Adds the project root to the MATLAB path recursively.
%       - Stores the fixture folder in root appdata under:
%             'ApplyAggregateFunctionTestConfig'
%       - Runs TestApplyAggregateFunction
%       - Clears the appdata configuration after the run
%
%   Notes:
%       - The test class itself copies fixture files into a fresh temporary
%         SaveFolder for each test method.
%       - Only green.dat and AcqInfos.mat are required in the fixture
%         folder. events.mat is recreated automatically by the tests.

    % ---------------------------------------------------------------------
    % Resolve project root from this wrapper location
    % ---------------------------------------------------------------------
    thisFile = mfilename('fullpath');
    thisFolder = fileparts(thisFile);

    projectRoot = extractBefore(thisFolder, [filesep 'test']);
    if isempty(projectRoot)
        projectRoot = fileparts(fileparts(thisFolder));
    end
    projectRoot = char(projectRoot);

    addpath(genpath(projectRoot));

    % ---------------------------------------------------------------------
    % Resolve fixture folder
    % ---------------------------------------------------------------------
    if nargin < 1 || isempty(sampleDataFolder)
        sampleDataFolder = fullfile( ...
            projectRoot, ...
            'test', ...
            'Analysis', ...
            'TestingData_with_events');
    end

    sampleDataFolder = char(string(sampleDataFolder));

    if ~isfolder(sampleDataFolder)
        error('run_TestApplyAggregateFunction:MissingFixtureFolder', ...
            'Fixture folder not found: %s', sampleDataFolder);
    end

    if ~isfile(fullfile(sampleDataFolder, 'green.dat'))
        error('run_TestApplyAggregateFunction:MissingGreenDat', ...
            'Fixture file "green.dat" was not found in: %s', sampleDataFolder);
    end

    if ~isfile(fullfile(sampleDataFolder, 'AcqInfos.mat'))
        error('run_TestApplyAggregateFunction:MissingAcqInfos', ...
            'Fixture file "AcqInfos.mat" was not found in: %s', sampleDataFolder);
    end

    % ---------------------------------------------------------------------
    % Store test configuration for the test class
    % ---------------------------------------------------------------------
    cfg = struct();
    cfg.sampleDataFolder = sampleDataFolder;
    setappdata(0, 'ApplyAggregateFunctionTestConfig', cfg);

    cleanupObj = onCleanup(@() rmappdataIfExists(0, 'ApplyAggregateFunctionTestConfig')); %#ok<NASGU>

    fprintf('Running TestApplyAggregateFunction\n');
    fprintf('Fixture folder: %s\n\n', sampleDataFolder);

    % ---------------------------------------------------------------------
    % Run tests
    % ---------------------------------------------------------------------
    results = runtests('TestApplyAggregateFunction');

    % ---------------------------------------------------------------------
    % Print short summary
    % ---------------------------------------------------------------------
    nPassed = sum([results.Passed]);
    nFailed = sum([results.Failed]);
    nIncomplete = sum([results.Incomplete]);

    fprintf('\nTestApplyAggregateFunction summary\n');
    fprintf('  Passed     : %d\n', nPassed);
    fprintf('  Failed     : %d\n', nFailed);
    fprintf('  Incomplete : %d\n', nIncomplete);

    if nFailed > 0 || nIncomplete > 0
        disp(table({results.Name}', [results.Passed]', [results.Failed]', [results.Incomplete]', ...
            'VariableNames', {'Name','Passed','Failed','Incomplete'}));
    end
end

function rmappdataIfExists(h, key)
%RMAPPDATAIFEXISTS Remove appdata key if it exists.

    if isappdata(h, key)
        rmappdata(h, key);
    end
end