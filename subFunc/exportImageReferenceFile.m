function filePath = exportImageReferenceFile(ImageReference, filePath, opts)
%EXPORTIMAGEREFERENCEFILE Write a standalone image reference file (source kind 'file').
%
%   filePath = exportImageReferenceFile(ImageReference, filePath)
%   filePath = exportImageReferenceFile(ImageReference, filePath, 'ClearIdentity', true)
%
%   Writes ImageReference (version 2; a version-1 struct is upgraded first)
%   as a MAT file holding the single variable 'ImageReference', atomically.
%   Works with or without a project: a reference created without one simply
%   has empty identity fields. The file can be selected as an alignment
%   reference (resolveImageReferenceSource kind 'file') or imported into a
%   subject library (importImageReferenceFile).
%
%   Name-Value options:
%       ClearIdentity - true empties the project, subject and session
%                       identity fields before writing (e.g. to share a
%                       reference outside its project). Default: false.
%
%   See also importImageReferenceFile, genImageReferenceStruct,
%   resolveImageReferenceSource.

arguments
    ImageReference (1,1) struct
    filePath (1,1) string
    opts.ClearIdentity (1,1) logical = false
end

ImageReference = normalizeImageReferenceStruct(ImageReference);
if opts.ClearIdentity
    ImageReference = clearImageReferenceIdentity(ImageReference);
end

filePath = char(filePath);
[folder, ~, ext] = fileparts(filePath);
if ~strcmpi(ext, '.mat')
    error('Umitoolbox:ImageRegistration:invalidExportPath', ...
        'Image reference files use the .mat extension: %s', filePath);
end
if ~isempty(folder) && ~isfolder(folder)
    error('Umitoolbox:ImageRegistration:invalidExportPath', ...
        'Folder does not exist: %s', folder);
end
saveMatAtomic(filePath, 'ImageReference', ImageReference);
end
