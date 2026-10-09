function results = runEventsManagerTests()
%RUNEVENTSMANAGERTESTS Run the EventsManager unit-test suite.
%
% Output:
%   results - matlab.unittest.TestResult array.

    rootFolder = fileparts(mfilename('fullpath'));
    addpath(rootFolder);
    addpath(fullfile(rootFolder, 'EMTestHelpers'));

    import matlab.unittest.TestSuite
    suite = TestSuite.fromFolder(rootFolder, 'IncludingSubfolders', false);
    results = run(suite);
end
