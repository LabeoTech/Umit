function DataParams = loadDataParams(folderPath)
%LOADDATAPARAMS Load DataParams.mat from a folder.
%
%   DataParams = loadDataParams(folderPath)
%
%   Loads the folder-global DataParams structure from DataParams.mat,
%   normalizes missing schema fields for backward compatibility, and
%   validates the result.
%
%   Input:
%       folderPath - Path to the save folder.
%
%   Output:
%       DataParams - Loaded DataParams structure.
%
%   Error:
%       Throws an error if the folder or DataParams.mat file does not exist,
%       or if the file does not contain the variable "DataParams".

if ~(ischar(folderPath) || isstring(folderPath))
    error('loadDataParams:InvalidFolderPath', ...
        'folderPath must be a char vector or string scalar.');
end

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('loadDataParams:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

filePath = fullfile(folderPath, 'DataParams.mat');

if ~exist(filePath, 'file')
    error('loadDataParams:MissingFile', ...
        'DataParams.mat was not found in folder: %s', folderPath);
end

S = load(filePath, 'DataParams');

if ~isfield(S, 'DataParams')
    error('loadDataParams:MissingVariable', ...
        'File "%s" does not contain variable "DataParams".', filePath);
end

DataParams = iNormalizeDataParams(S.DataParams);
validateDataParams(DataParams);
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

if ~isfield(DataParams, 'registration') || ~isstruct(DataParams.registration) || ...
        ~isscalar(DataParams.registration)
    DataParams.registration = struct();
end

if ~isfield(DataParams.registration, 'resourceUUID')
    DataParams.registration.resourceUUID = '';
else
    DataParams.registration.resourceUUID = char(string(DataParams.registration.resourceUUID));
end

if ~isfield(DataParams, 'cameraCoregistration') || ...
        ~isstruct(DataParams.cameraCoregistration) || ...
        ~isscalar(DataParams.cameraCoregistration)
    DataParams.cameraCoregistration = struct();
end

ccFieldDefaults = struct( ...
    'isCoregistered', false, ...
    'isReviewed', false, ...
    'tform', [], ...
    'resourceUUID', '', ...
    'rigID', '', ...
    'transformType', '', ...
    'method', '', ...
    'sourceFile', '', ...
    'sourceFileTimestamp', '', ...
    'createdOn', '', ...
    'source', '', ...
    'qcStatus', '', ...
    'qcWarning', '', ...
    'qcFigureFile', '', ...
    'qcPreviewImageFile', '', ...
    'appliedOn', '', ...
    'appliedBy', '', ...
    'confirmationMode', '', ...
    'notes', '');

ccFieldNames = fieldnames(ccFieldDefaults);
for k = 1:numel(ccFieldNames)
    thisField = ccFieldNames{k};
    if ~isfield(DataParams.cameraCoregistration, thisField)
        DataParams.cameraCoregistration.(thisField) = ccFieldDefaults.(thisField);
    end
end

if ~isfield(DataParams.cameraCoregistration, 'qcMetrics') || ...
        ~isstruct(DataParams.cameraCoregistration.qcMetrics) || ...
        ~isscalar(DataParams.cameraCoregistration.qcMetrics)
    DataParams.cameraCoregistration.qcMetrics = struct();
end

ccQcMetricsFields = {'MIBefore', 'MIAfter', 'MIDelta', 'translationXY_px', ...
    'rotationDeg', 'scaleXY', 'determinant'};
for k = 1:numel(ccQcMetricsFields)
    thisField = ccQcMetricsFields{k};
    if ~isfield(DataParams.cameraCoregistration.qcMetrics, thisField)
        DataParams.cameraCoregistration.qcMetrics.(thisField) = [];
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
