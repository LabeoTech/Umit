function DataParams = createDataParams(folderPath, varargin)
%CREATEDATAPARAMS Create folder-global DataParams.mat in a save folder.
%
%   DataParams = createDataParams(folderPath)
%   DataParams = createDataParams(folderPath, 'overwrite', true)
%   DataParams = createDataParams(folderPath, 'axisConvention', 'imageXY_topLeft')
%
%   Creates a default folder-global DataParams structure and saves it as
%   DataParams.mat in the specified folder.
%
%   If "AcqInfos.mat" exists in the folder, this function loads the
%   structure "AcqInfoStream" and uses its "Height" and "Width" fields to:
%       - populate DataParams.view.imageSizeYX
%       - create a default full TRUE mask
%       - assign default mask metadata
%
%   If "AcqInfos.mat" does not exist, the related fields are left empty.
%
%   If DataParams.mat already exists:
%       - By default, a warning is issued and the existing file is loaded
%         and returned unchanged.
%       - If 'overwrite' is true, the existing file is replaced.
%
%   Input:
%       folderPath - Path to the target folder.
%
%   Name-Value options:
%       overwrite      - Logical scalar. If true, overwrite an existing
%                        DataParams.mat file. Default: false
%       axisConvention - Coordinate convention string.
%
%   Output:
%       DataParams - Created or existing DataParams structure.
%
%   Notes:
%       - This file stores folder-global spatial and shared analysis
%         parameters for all .dat files in the folder.
%       - Binary reconstruction metadata should remain in AcqInfos.mat.

p = inputParser;
p.FunctionName = 'createDataParams';

addRequired(p, 'folderPath', @(x) ischar(x) || isstring(x));
addParameter(p, 'overwrite', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'axisConvention', 'imageXY_topLeft', @(x) ischar(x) || isstring(x));

parse(p, folderPath, varargin{:});
R = p.Results;

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('createDataParams:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

dataParamsPath = fullfile(folderPath, 'DataParams.mat');

if exist(dataParamsPath, 'file') && ~R.overwrite
    warning('createDataParams:FileExists', ...
        'DataParams.mat already exists in folder. Existing file was kept.');
    DataParams = loadDataParams(folderPath);
    return
end

tNow = datetime('now');

% -------------------------------------------------------------
% Default values
% -------------------------------------------------------------
defaultImageSizeYX = [];
defaultMask = [];
defaultOrigin = [];
defaultPixelSize = [];

% -------------------------------------------------------------
% Load AcqInfos.mat if available
% -------------------------------------------------------------
acqInfoPath = fullfile(folderPath, 'AcqInfos.mat');

if exist(acqInfoPath, 'file')
    S = load(acqInfoPath, 'AcqInfoStream');

    if ~isfield(S, 'AcqInfoStream')
        error('createDataParams:MissingAcqInfoStream', ...
            'File "%s" does not contain variable "AcqInfoStream".', acqInfoPath);
    end

    AcqInfoStream = S.AcqInfoStream;

    if ~isstruct(AcqInfoStream) || ~isscalar(AcqInfoStream)
        error('createDataParams:InvalidAcqInfoStream', ...
            'Variable "AcqInfoStream" in "%s" must be a scalar struct.', acqInfoPath);
    end

    if ~isfield(AcqInfoStream, 'Height') || ~isfield(AcqInfoStream, 'Width')
        error('createDataParams:MissingAcqFields', ...
            'AcqInfoStream must contain fields "Height" and "Width".');
    end

    Ny = AcqInfoStream.Height;
    Nx = AcqInfoStream.Width;

    if ~isnumeric(Ny) || ~isscalar(Ny) || ~isfinite(Ny) || Ny <= 0 || mod(Ny, 1) ~= 0
        error('createDataParams:InvalidHeight', ...
            'AcqInfoStream.Height must be a positive integer scalar.');
    end

    if ~isnumeric(Nx) || ~isscalar(Nx) || ~isfinite(Nx) || Nx <= 0 || mod(Nx, 1) ~= 0
        error('createDataParams:InvalidWidth', ...
            'AcqInfoStream.Width must be a positive integer scalar.');
    end

    defaultImageSizeYX = [Ny, Nx];
    defaultMask = true(Ny, Nx);
    defaultOrigin = [1, 1];
    defaultPixelSize = [];
end

% -------------------------------------------------------------
% Build DataParams structure
% -------------------------------------------------------------
DataParams = struct();
DataParams.schemaVersion = '1.0';
DataParams.createdOn = tNow;
DataParams.lastModified = tNow;
DataParams.notes = "";

DataParams.view = struct();
DataParams.view.pixelSize_px_per_mm = defaultPixelSize;
DataParams.view.origin_xy_px = defaultOrigin;
DataParams.view.imageSizeYX = defaultImageSizeYX;
DataParams.view.axisConvention = char(string(R.axisConvention));

DataParams.mask = struct();
DataParams.mask.logical = defaultMask;
DataParams.mask.name = '';
DataParams.mask.description = '';
DataParams.mask.space = 'native';
DataParams.mask.createdOn = [];
DataParams.mask.source = '';

if ~isempty(defaultMask)
    DataParams.mask.name = 'default';
    DataParams.mask.description = 'Default mask';
    DataParams.mask.createdOn = tNow;
    DataParams.mask.source = 'createDataParams';
end

DataParams.registration = struct();
DataParams.registration.isRegistered = false;
DataParams.registration.tform = [];
DataParams.registration.transformType = '';
DataParams.registration.method = '';
DataParams.registration.referenceDescription = '';
DataParams.registration.referenceFile = '';
DataParams.registration.referenceImage = [];
DataParams.registration.createdOn = [];
DataParams.registration.source = '';

DataParams.custom = struct();

validateDataParams(DataParams);
save(dataParamsPath, 'DataParams', '-mat');
end