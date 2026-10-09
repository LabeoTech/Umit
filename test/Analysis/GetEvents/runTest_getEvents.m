function results = runTest_getEvents(parentFolder)
% RUNTEST_GETEVENTS Run the TestGetEvents unit-test class.
%
% Syntax:
%   results = runTest_getEvents
%   results = runTest_getEvents(parentFolder)
%
% Input:
%   parentFolder - Optional parent folder containing the getEvents test
%                  case subfolders. If provided, it is stored as the
%                  preference used by TestGetEvents.
%
% Output:
%   results      - matlab.unittest.TestResult array.
%
% Notes:
%   - If "parentFolder" is omitted, TestGetEvents resolves the case root
%     using its internal priority:
%         1) environment variable UMIT_GETEVENTS_TEST_PARENT
%         2) MATLAB preference
%         3) test class folder
%   - This wrapper prints a compact summary table after execution.

    if nargin >= 1 && ~isempty(parentFolder)
        validateattributes(parentFolder, {'char','string'}, {'scalartext'}, ...
            mfilename, 'parentFolder');
        parentFolder = char(string(parentFolder));
        assert(isfolder(parentFolder), ...
            'The parent folder "%s" does not exist.', parentFolder);

        setpref('umIToolbox', 'GetEventsTestParentFolder', parentFolder);
    end

    results = runtests('TestGetEvents');
    disp(table(results))
end