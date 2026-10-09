function writeFailedTests(folder, res)
%writeFailedTests Save failed tests report to failed_tests.txt in folder
%   folder - path to folder where failed_tests.txt will be saved
%   res    - array of test result structures (as returned by runPipelineManagerTests)

if nargin < 2
    error('writeFailedTests requires two inputs: folder and res');
end
if ~ischar(folder) && ~isstring(folder)
    error('folder must be a character vector or string');
end
folder = char(folder);

% Ensure folder exists
if ~isfolder(folder)
    error('Specified folder does not exist: %s', folder);
end

failed = res([res.Failed] == 1);

fid = fopen(fullfile(folder, 'failed_tests.txt'), 'w');
if fid == -1
    error('Unable to open file for writing: %s', fullfile(folder, 'failed_tests.txt'));
end

for ii = 1:numel(failed)
    fprintf(fid, 'Test %s failed: %s\n', ...
        failed(ii).Name, failed(ii).Details.DiagnosticRecord.Report);
end

fclose(fid);
end