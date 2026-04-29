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

for iROI = 1:numel(ROIFile.ROIs)
    if isfield(ROIFile.ROIs(iROI), 'ID')
        ROIFile.ROIs(iROI) = rmfield(ROIFile.ROIs(iROI), 'ID');
    end

    if isfield(ROIFile.ROIs(iROI), 'runtime')
        ROIFile.ROIs(iROI) = rmfield(ROIFile.ROIs(iROI), 'runtime');
    end
end

end
