function outFiles = funcTemplateFileManifest(data, SaveFolder)
%FUNCTEMPLATEFILEMANIFEST Template for a file-manifest pipeline output.
%
%   OUTFILES = FUNCTEMPLATEFILEMANIFEST(DATA, SAVEFOLDER) calculates a few
%   descriptive values for a numeric array, writes a MAT file and a text
%   report in SAVEFOLDER, and returns their relative filenames as one file
%   manifest. A file manifest is the list of files actually created by one
%   logical PipelineManager output port.
%
%   INFO = FUNCTEMPLATEFILEMANIFEST('pipelineInfo') returns the metadata
%   used by PipelineManager to discover and run the function.
%
% Inputs
%   data       - Nonempty, real numeric or logical array held in RAM.
%   SaveFolder - Existing save-folder path supplied by PipelineManager.
%                All reported files must be created inside this folder.
%
% Output
%   outFiles - Cell array containing the relative filenames
%              {'funcTemplateSummary.mat','funcTemplateSummary.txt'}.
%              Return only files that were successfully created.
%
% Files Created or Modified
%   funcTemplateSummary.mat - Contains scalar structure "summary" with
%                             input size, class, minimum, maximum, and mean.
%   funcTemplateSummary.txt - Human-readable form of the same summary.
%
% Example
%   data = reshape(1:24, [2 3 4]);
%   saveFolder = tempdir;
%   outFiles = funcTemplateFileManifest(data, saveFolder);
%
% PipelineManager Compatibility
%   outputMode='file' tells PipelineManager that the function creates its
%   own files and returns their filenames. defOutfilename lists all files
%   the default execution may create, while OUTFILES reports the files
%   actually created. This report is declared isData=false because it is a
%   side artifact, not a downstream pipeline data source. A copied
%   function should use isData=true only when every returned file is valid
%   data that PipelineManager can route to another analysis function.
%
%   Copy this file to Analysis/<Category>/<YourFunction>.m, then rename the
%   function, pipelineInfo name, filenames, and summary calculation.
%   Preserve the one-argument 'pipelineInfo' query.
%
% See also PIPELINEMANAGER, PIPELINEMANAGER.ADDOUTPUT, FULLFILE.

defaultOutput = {'funcTemplateSummary.mat', 'funcTemplateSummary.txt'};

% KEEP: PipelineManager discovers modern functions through this query.
if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outFiles = localGetPipelineInfo(defaultOutput);
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder');
parse(p, data, SaveFolder);

data = p.Results.data;
SaveFolder = p.Results.SaveFolder;
localValidateData(data);
localValidateSaveFolder(SaveFolder);
SaveFolder = char(string(SaveFolder));

% EDIT: Replace this summary with the desired file-producing calculation.
numericData = double(data(:));
summary = struct();
summary.Size = size(data);
summary.MatlabClass = class(data);
summary.Minimum = min(numericData, [], 'omitnan');
summary.Maximum = max(numericData, [], 'omitnan');
summary.Mean = mean(numericData, 'omitnan');

matPath = fullfile(SaveFolder, defaultOutput{1});
textPath = fullfile(SaveFolder, defaultOutput{2});
save(matPath, 'summary', '-mat');

fid = fopen(textPath, 'w');
if fid == -1
    error('Umitoolbox:funcTemplateFileManifest:FileOpenFailed', ...
        'Could not create report file "%s".', textPath);
end
fileCleanup = onCleanup(@() safeFclose(fid));

fprintf(fid, 'Function template summary\n');
fprintf(fid, 'Size: %s\n', mat2str(summary.Size));
fprintf(fid, 'MATLAB class: %s\n', summary.MatlabClass);
fprintf(fid, 'Minimum: %.17g\n', summary.Minimum);
fprintf(fid, 'Maximum: %.17g\n', summary.Maximum);
fprintf(fid, 'Mean: %.17g\n', summary.Mean);

% Return SaveFolder-relative names so the manifest remains relocatable.
outFiles = defaultOutput;

end

function info = localGetPipelineInfo(defaultOutput)
%LOCALGETPIPELINEINFO Declare the PipelineManager integration contract.

info = PipelineManager.createPipelineInfo( ...
    'funcTemplateFileManifest', ...
    'Write summary artifacts and return their file manifest.');

info = PipelineManager.addInput( ...
    info, ...
    'data', ...
    'UnknownDataType', ...
    'Numeric or logical data summarized in RAM.', ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', true, ...
    'supportsFile', false, ...
    'dataMode', 'ram');

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    'Folder where the declared report files are created.', ...
    'kind', 'input', ...
    'position', 2, ...
    'callType', 'positional', ...
    'isData', false);

% A file manifest is one logical output even when it contains many files.
info = PipelineManager.addOutput( ...
    info, ...
    'outFiles', ...
    'UnknownDataType', ...
    'file', ...
    'Manifest of report artifacts created in SaveFolder.', ...
    defaultOutput, ...
    1, ...
    'isData', false, ...
    'saveFileName', '');

info.notes = { ...
    'Return only files that the function actually created.'; ...
    'Use isData=true only for files that are valid downstream pipeline data.'};

end

function localValidateData(data)
%LOCALVALIDATEDATA Require nonempty real numeric or logical input.

if ~((isnumeric(data) || islogical(data)) && isreal(data) && ~isempty(data))
    error('Umitoolbox:funcTemplateFileManifest:InvalidData', ...
        'DATA must be a nonempty, real numeric or logical array.');
end
end

function localValidateSaveFolder(SaveFolder)
%LOCALVALIDATESAVEFOLDER Require an existing dataset save folder.

isTextScalar = ischar(SaveFolder) || ...
    (isstring(SaveFolder) && isscalar(SaveFolder));
if ~isTextScalar || ~isfolder(SaveFolder)
    error('Umitoolbox:funcTemplateFileManifest:InvalidSaveFolder', ...
        'SAVEFOLDER must be the path to an existing folder.');
end
end
