function result = exportReportOrLogToJson(parentFigure, startFolder)
%EXPORTREPORTORLOGTOJSON Ask for a report or log MAT file and write a JSON copy of it.
%
%   result = exportReportOrLogToJson(parentFigure, startFolder)
%
%   The DataViewer menu action "Report or log to JSON..." (compiled
%   application only) calls this helper. It lets the user pick one of the
%   report and log files DataViewer and its tools write (the alignment
%   report manualAlignmentReport_*.mat, pipeLog.mat, dataHistory.mat,
%   MotionCorrectionLog.mat), asks where to save the JSON copy (default: the
%   same folder and name with the .json extension), converts it with
%   exportMatToJson and reports where the file is and whether any large array
%   was summarised.
%
%   Inputs:
%       parentFigure - uifigure that owns the dialogs.
%       startFolder  - folder the file dialog opens in (default: current folder).
%
%   Output:
%       result - the exportMatToJson result, or [] when the user cancelled or
%                the conversion failed (the failure is shown in a dialog).
%
%   See also exportMatToJson.

if nargin < 2 || isempty(startFolder) || ~isfolder(startFolder)
    startFolder = pwd;
end

result = [];

filterSpec = { ...
    'manualAlignmentReport_*.mat;pipeLog.mat;dataHistory.mat;MotionCorrectionLog.mat', ...
    'Reports and logs (manualAlignmentReport_*.mat, pipeLog.mat, dataHistory.mat, MotionCorrectionLog.mat)'};
[fileName, folderName] = uigetfile(filterSpec, 'Select a report or log file', startFolder);
if isequal(fileName, 0)
    return
end

[~, baseName] = fileparts(fileName);
[jsonName, jsonFolder] = uiputfile({'*.json', 'JSON file (*.json)'}, ...
    'Save the JSON copy as', fullfile(folderName, [baseName '.json']));
if isequal(jsonName, 0)
    return
end

try
    result = exportMatToJson(fullfile(folderName, fileName), fullfile(jsonFolder, jsonName));
catch ME
    result = [];
    uialert(parentFigure, ME.message, 'JSON export failed', 'Icon', 'error');
    return
end

message = sprintf('JSON copy written:\n%s', result.JsonFile);
if ~isempty(result.SummarisedVariables)
    message = sprintf('%s\n\nLarge arrays were summarised (not written in full) in: %s', ...
        message, strjoin(cellstr(result.SummarisedVariables), ', '));
end
uialert(parentFigure, message, 'JSON export', 'Icon', 'success');
end
