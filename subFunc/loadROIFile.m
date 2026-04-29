function ROIFile = loadROIFile(filePath)
%LOADROIFILE Load, migrate, and validate a UMIT .roi file.
%
%   ROIFile = loadROIFile(filePath)
%
%   Loads variable ROIFile from a MAT-backed .roi file, applies supported
%   schema migration, and validates the result.

if ~(ischar(filePath) || isstring(filePath))
    error('loadROIFile:InvalidFilePath', ...
        'filePath must be a char vector or string scalar.');
end

filePath = char(string(filePath));

if isempty(fileparts(filePath))
    filePath = fullfile(pwd, filePath);
end

if ~isfile(filePath)
    error('loadROIFile:FileNotFound', ...
        'File not found: %s', filePath);
end

[~, ~, ext] = fileparts(filePath);
if ~strcmpi(ext, '.roi')
    error('loadROIFile:InvalidExtension', ...
        'ROI files must use the .roi extension.');
end

S = load(filePath, 'ROIFile');

if ~isfield(S, 'ROIFile')
    error('loadROIFile:MissingROIFileVariable', ...
        'File does not contain variable ROIFile: %s', filePath);
end

ROIFile = migrateROIFile(S.ROIFile);
validateROIFile(ROIFile);

end
