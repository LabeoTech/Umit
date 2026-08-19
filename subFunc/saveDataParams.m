function saveDataParams(folderPath, DataParams)
%SAVEDATAPARAMS Validate and save DataParams.mat to a folder.
%
%   saveDataParams(folderPath, DataParams)
%
%   Normalizes the supplied DataParams structure for the current schema,
%   validates it, updates the "lastModified" timestamp, and saves it as
%   DataParams.mat in the specified folder.
%
%   Inputs:
%       folderPath  - Path to the folder.
%       DataParams  - Folder-global DataParams structure.
%
%   Notes:
%       - The folder must already exist.
%       - The variable is saved with the name "DataParams".
%       - SaveFolder-level RawFolder context is stored under
%         DataParams.folders.

if ~(ischar(folderPath) || isstring(folderPath))
    error('saveDataParams:InvalidFolderPath', ...
        'folderPath must be a char vector or string scalar.');
end

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('saveDataParams:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

DataParams = iNormalizeDataParams(DataParams);
validateDataParams(DataParams);
DataParams.lastModified = datetime('now');

filePath = fullfile(folderPath, 'DataParams.mat');
saveMatAtomic(filePath, 'DataParams', DataParams);
end

function DataParams = iNormalizeDataParams(DataParams)
%INORMALIZEDATAPARAMS Add backward-compatible fields before validation.

if ~isstruct(DataParams) || ~isscalar(DataParams)
    return
end

if ~isfield(DataParams, 'folders') || ~isstruct(DataParams.folders) || ...
        ~isscalar(DataParams.folders)
    DataParams.folders = struct();
end

if ~isfield(DataParams.folders, 'RawFolder') || isempty(DataParams.folders.RawFolder)
    legacyRawFolder = '';

    if isfield(DataParams, 'RawFolder') && ~isempty(DataParams.RawFolder)
        legacyRawFolder = DataParams.RawFolder;
    elseif isfield(DataParams, 'rawFolder') && ~isempty(DataParams.rawFolder)
        legacyRawFolder = DataParams.rawFolder;
    end

    if isempty(legacyRawFolder)
        DataParams.folders.RawFolder = 'Missing';
    else
        DataParams.folders.RawFolder = char(string(legacyRawFolder));
    end
end

DataParams.folders.RawFolder = char(string(DataParams.folders.RawFolder));

if ~isfield(DataParams.folders, 'RawFolderStatus') || ...
        isempty(DataParams.folders.RawFolderStatus)
    DataParams.folders.RawFolderStatus = iInferRawFolderStatus(DataParams.folders.RawFolder);
else
    DataParams.folders.RawFolderStatus = char(string(DataParams.folders.RawFolderStatus));
end

if ~isfield(DataParams.folders, 'RawFolderSetOn')
    DataParams.folders.RawFolderSetOn = [];
end

if ~isfield(DataParams.folders, 'RawFolderSetBy') || isempty(DataParams.folders.RawFolderSetBy)
    DataParams.folders.RawFolderSetBy = '';
else
    DataParams.folders.RawFolderSetBy = char(string(DataParams.folders.RawFolderSetBy));
end

if ~isfield(DataParams, 'registration') || ...
        ~isstruct(DataParams.registration) || ...
        ~isscalar(DataParams.registration)
    DataParams.registration = struct();
end

registrationProvenanceDefaults = struct( ...
    'imageReferenceUUID', '', ...
    'imageReferenceChecksum', '');
provenanceFields = fieldnames(registrationProvenanceDefaults);
for k = 1:numel(provenanceFields)
    thisField = provenanceFields{k};
    if ~isfield(DataParams.registration, thisField)
        DataParams.registration.(thisField) = ...
            registrationProvenanceDefaults.(thisField);
    else
        DataParams.registration.(thisField) = char(string( ...
            DataParams.registration.(thisField)));
    end
end

end

function statusText = iInferRawFolderStatus(rawFolder)
%IINFERRAWFOLDERSTATUS Infer status from a RawFolder value.

rawFolder = char(string(rawFolder));

if isempty(rawFolder) || strcmpi(rawFolder, 'Missing')
    statusText = 'missing';
elseif isfolder(rawFolder)
    statusText = 'valid';
else
    statusText = 'invalid';
end

end
