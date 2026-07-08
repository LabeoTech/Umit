function DataParams = createDataParams(folderPath, varargin)
%CREATEDATAPARAMS Create folder-global DataParams.mat in a save folder.
%
%   DataParams = createDataParams(folderPath)
%   DataParams = createDataParams(folderPath, 'overwrite', true)
%   DataParams = createDataParams(folderPath, ...
%       'axisConvention', 'imageXY_topLeft')
%
%   Creates a default folder-global DataParams structure and saves it as
%   DataParams.mat in the specified folder.
%
%   If AcqInfos.mat exists in the folder, this function loads
%   AcqInfoStream and uses its Height and Width fields to:
%       - populate DataParams.view.imageSizeYX
%       - create a default full TRUE mask
%       - set the default image origin to [1 1]
%
%   If AcqInfos.mat does not exist, the related spatial fields are left
%   empty.
%
%   If DataParams.mat already exists:
%       - By default, a warning is issued and the existing file is loaded,
%         schema-normalized in memory, validated, and returned.
%       - If overwrite is true, the existing file is replaced.
%
%   Inputs:
%       folderPath - Target folder where DataParams.mat will be created.
%
%   Name-Value options:
%       overwrite      - Logical scalar. If true, overwrite an existing
%                        DataParams.mat file. Default: false.
%       axisConvention - Coordinate convention string. Default:
%                        'imageXY_topLeft'.
%
%   Output:
%       DataParams     - Created or existing DataParams structure.
%
%   Notes:
%       - DataParams stores folder-global spatial/view parameters, masks,
%         registration metadata, folder context, and user-defined custom
%         metadata.
%       - Binary reconstruction metadata should remain in AcqInfos.mat.
%       - RawFolder is SaveFolder-level context and is stored under
%         DataParams.folders, not in AcqInfos.mat.

p = inputParser;
p.FunctionName = 'createDataParams';

addRequired(p, 'folderPath', @(x) ischar(x) || isstring(x));
addParameter(p, 'overwrite', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'axisConvention', 'imageXY_topLeft', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

parse(p, folderPath, varargin{:});
R = p.Results;

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('createDataParams:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

dataParamsPath = fullfile(folderPath, 'DataParams.mat');

% -------------------------------------------------------------------------
% Preserve existing DataParams.mat unless overwrite is requested.
% -------------------------------------------------------------------------
if exist(dataParamsPath, 'file') && ~R.overwrite
    warning('createDataParams:FileExists', ...
        'DataParams.mat already exists in folder. Existing file was kept.');

    DataParams = iLoadExistingDataParams(dataParamsPath);
    DataParams = iEnsureDataParamsFolderFields(DataParams);
    validateDataParams(DataParams);
    return
end

tNow = datetime('now');

% -------------------------------------------------------------------------
% Default spatial/view values
% -------------------------------------------------------------------------
defaultImageSizeYX = [];
defaultMask = [];
defaultOrigin = [];
defaultPixelSize = [];

% -------------------------------------------------------------------------
% Load AcqInfos.mat if available.
% -------------------------------------------------------------------------
acqInfoPath = fullfile(folderPath, 'AcqInfos.mat');

if exist(acqInfoPath, 'file')
    S = load(acqInfoPath, 'AcqInfoStream');

    if ~isfield(S, 'AcqInfoStream')
        error('createDataParams:MissingAcqInfoStream', ...
            'File "%s" does not contain variable "AcqInfoStream".', ...
            acqInfoPath);
    end

    AcqInfoStream = S.AcqInfoStream;

    if ~isstruct(AcqInfoStream) || ~isscalar(AcqInfoStream)
        error('createDataParams:InvalidAcqInfoStream', ...
            'Variable "AcqInfoStream" in "%s" must be a scalar struct.', ...
            acqInfoPath);
    end

    if ~isfield(AcqInfoStream, 'Height') || ...
            ~isfield(AcqInfoStream, 'Width')
        error('createDataParams:MissingAcqFields', ...
            'AcqInfoStream must contain fields "Height" and "Width".');
    end

    Ny = double(AcqInfoStream.Height);
    Nx = double(AcqInfoStream.Width);

    validateattributes(Ny, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
        'createDataParams', 'AcqInfoStream.Height');

    validateattributes(Nx, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
        'createDataParams', 'AcqInfoStream.Width');

    defaultImageSizeYX = [Ny, Nx];
    defaultMask = true(Ny, Nx);
    defaultOrigin = [1, 1];
    defaultPixelSize = [];
end

% -------------------------------------------------------------------------
% Build DataParams structure.
% -------------------------------------------------------------------------
DataParams = struct();

DataParams.schemaVersion = '1.1';
DataParams.createdOn = tNow;
DataParams.lastModified = tNow;
DataParams.notes = "";

% -------------------------------------------------------------------------
% View metadata
% -------------------------------------------------------------------------
DataParams.view = struct();
DataParams.view.pixelSize_px_per_mm = defaultPixelSize;
DataParams.view.origin_xy_px = defaultOrigin;
DataParams.view.imageSizeYX = defaultImageSizeYX;
DataParams.view.axisConvention = char(string(R.axisConvention));

% -------------------------------------------------------------------------
% Mask metadata
% -------------------------------------------------------------------------
DataParams.mask = struct();
DataParams.mask.logical = defaultMask;
DataParams.mask.name = '';
DataParams.mask.description = '';
DataParams.mask.space = 'native';
DataParams.mask.createdOn = [];
DataParams.mask.source = '';

if ~isempty(defaultMask)
    DataParams.mask.name = 'default';
    DataParams.mask.description = 'Default full-field logical mask.';
    DataParams.mask.createdOn = tNow;
    DataParams.mask.source = 'createDataParams';
end

% -------------------------------------------------------------------------
% Folder-context metadata
% -------------------------------------------------------------------------
DataParams.folders = struct();
DataParams.folders.RawFolder = 'Missing';
DataParams.folders.RawFolderStatus = 'missing';
DataParams.folders.RawFolderSetOn = [];
DataParams.folders.RawFolderSetBy = '';

% -------------------------------------------------------------------------
% Registration metadata
% -------------------------------------------------------------------------
DataParams.registration = struct();

DataParams.registration.isRegistered = false;
DataParams.registration.isReviewed = false;
DataParams.registration.tform = [];

DataParams.registration.transformType = '';
DataParams.registration.method = '';
DataParams.registration.referenceDescription = '';
DataParams.registration.referenceFile = '';
DataParams.registration.referenceImage = [];
DataParams.registration.createdOn = '';
DataParams.registration.source = '';

DataParams.registration.qcStatus = '';
DataParams.registration.qcWarning = '';
DataParams.registration.qcFigureFile = '';
DataParams.registration.qcPreviewImageFile = '';

DataParams.registration.appliedOn = '';
DataParams.registration.appliedBy = '';
DataParams.registration.confirmationMode = '';
DataParams.registration.notes = '';

DataParams.registration.qcMetrics = struct();
DataParams.registration.qcMetrics.MIBefore = [];
DataParams.registration.qcMetrics.MIAfter = [];
DataParams.registration.qcMetrics.MIDelta = [];
DataParams.registration.qcMetrics.translationXY_px = [];
DataParams.registration.qcMetrics.rotationDeg = [];
DataParams.registration.qcMetrics.scaleXY = [];
DataParams.registration.qcMetrics.determinant = [];

% -------------------------------------------------------------------------
% Custom metadata
% -------------------------------------------------------------------------
DataParams.custom = struct();

% -------------------------------------------------------------------------
% Validate and save.
% -------------------------------------------------------------------------
validateDataParams(DataParams);

save(dataParamsPath, 'DataParams', '-mat');

end

% =========================================================================
% Local helpers
% =========================================================================

function DataParams = iLoadExistingDataParams(dataParamsPath)
%ILOADEXISTINGDATAPARAMS Load DataParams from an existing MAT file.

S = load(dataParamsPath, 'DataParams');

if ~isfield(S, 'DataParams')
    error('createDataParams:MissingDataParamsVariable', ...
        'Existing file "%s" does not contain variable "DataParams".', ...
        dataParamsPath);
end

DataParams = S.DataParams;

end

function DataParams = iEnsureDataParamsFolderFields(DataParams)
%IENSUREDATAPARAMSFOLDERFIELDS Add folder-context fields for older files.

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

if ~isfield(DataParams.folders, 'RawFolderStatus') || ...
        isempty(DataParams.folders.RawFolderStatus)
    DataParams.folders.RawFolderStatus = iInferRawFolderStatus(DataParams.folders.RawFolder);
end

if ~isfield(DataParams.folders, 'RawFolderSetOn')
    DataParams.folders.RawFolderSetOn = [];
end

if ~isfield(DataParams.folders, 'RawFolderSetBy') || isempty(DataParams.folders.RawFolderSetBy)
    DataParams.folders.RawFolderSetBy = '';
end

DataParams.folders.RawFolder = char(string(DataParams.folders.RawFolder));
DataParams.folders.RawFolderStatus = char(string(DataParams.folders.RawFolderStatus));
DataParams.folders.RawFolderSetBy = char(string(DataParams.folders.RawFolderSetBy));

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
