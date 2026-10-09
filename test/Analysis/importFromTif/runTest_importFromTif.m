function results = runTest_importFromTif()
% RUNTEST_IMPORTFROMTIF Run the TestImportFromTif unit-test class.
%
% Syntax:
%   results = runTest_importFromTif
%
% Output:
%   results - matlab.unittest.TestResult array.
%
% Notes:
%   - This wrapper runs the full TestImportFromTif class and prints a
%     compact summary table.

    results = runtests('TestImportFromTif');
    disp(table(results))
end