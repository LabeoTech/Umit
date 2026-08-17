function saveROIFile(filePath, ROIFile)
%SAVEROIFILE Validate and save a UMIT ROI file structure.
%
%   saveROIFile(filePath, ROIFile)
%
%   Saves ROIFile as a MAT file with custom .roi extension. Runtime-only
%   fields such as ID and runtime are removed before validation/saving.

if ~(ischar(filePath) || isstring(filePath))
    error('saveROIFile:InvalidFilePath', ...
        'filePath must be a char vector or string scalar.');
end

filePath = char(string(filePath));
[folderPath, fileName, ext] = fileparts(filePath);

if isempty(folderPath)
    folderPath = pwd;
end

if isempty(fileName)
    error('saveROIFile:InvalidFilePath', ...
        'filePath must include a file name.');
end

if isempty(ext)
    ext = '.roi';
elseif ~strcmpi(ext, '.roi')
    error('saveROIFile:InvalidExtension', ...
        'ROI files must use the .roi extension.');
end

if ~isfolder(folderPath)
    error('saveROIFile:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

ROIFile = iStripRuntimeFields(ROIFile);
ROIFile.modifiedOn = datetime('now');

validateROIFile(ROIFile);

filePath = fullfile(folderPath, [fileName, ext]);
save(filePath, 'ROIFile', '-mat');

end

function ROIFile = iStripRuntimeFields(ROIFile)
%ISTRIPRUNTIMEFIELDS Remove app-only fields before saving.

if ~isfield(ROIFile, 'ROIs') || isempty(ROIFile.ROIs)
    return
end

runtimeFields = intersect({'ID', 'runtime'}, fieldnames(ROIFile.ROIs));
if ~isempty(runtimeFields)
    % Remove fields from the complete struct array in one operation. Removing
    % fields one element at a time can create heterogeneous-struct assignment
    % errors when the array contains more than one ROI.
    ROIFile.ROIs = rmfield(ROIFile.ROIs, runtimeFields);
end

end
