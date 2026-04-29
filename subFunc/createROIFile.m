function ROIFile = createROIFile(imageSizeYX, varargin)
%CREATEROIFILE Create an empty UMIT ROI file structure.
%
%   ROIFile = createROIFile(imageSizeYX)
%   ROIFile = createROIFile(imageSizeYX, Name, Value, ...)
%
%   Creates a versioned ROIFile structure following getROISchema(1).
%   The structure can be saved later using saveROIFile.
%
%   Inputs:
%       imageSizeYX - [Height Width] of the image used to define ROI masks.
%
%   Name-Value options:
%       dataFile             - Source image file path or name. Default: ''.
%       entryName            - UMT entry name. Empty for .dat. Default: ''.
%       coordinateSystem     - Coordinate convention. Default: 'imageXY_px'.
%       originXY_px          - [x y] origin in image pixels. Default: [1 1].
%       pixelSize_px_per_mm  - [] or positive scalar. Default: [].
%       statsImage           - Single image used for ROI statistics. Default: [].
%       statsDescription     - Description of statsImage. Default: ''.
%       frameIndex           - Source/display frame index. Default: [].
%       viewMode             - '', 'normal', or 'event'. Default: ''.
%       eventCondition       - Event condition name. Default: ''.
%       eventRepetition      - Event repetition value. Default: ''.
%       ROIs                 - ROI struct array. Default: struct([]).
%
%   Output:
%       ROIFile              - Versioned ROI file structure.

p = inputParser;
p.FunctionName = 'createROIFile';

addRequired(p, 'imageSizeYX', @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'dataFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'entryName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'coordinateSystem', 'imageXY_px', @(x) ischar(x) || isstring(x));
addParameter(p, 'originXY_px', [1 1], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'pixelSize_px_per_mm', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'statsImage', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'statsDescription', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'frameIndex', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'viewMode', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'eventCondition', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'eventRepetition', '', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'ROIs', struct([]), @(x) isstruct(x));

parse(p, imageSizeYX, varargin{:});
R = p.Results;

schema = getROISchema(1);

imageSizeYX = double(imageSizeYX(:).');
if numel(imageSizeYX) ~= 2 || any(~isfinite(imageSizeYX)) || ...
        any(imageSizeYX <= 0) || any(imageSizeYX ~= round(imageSizeYX))
    error('createROIFile:InvalidImageSize', ...
        'imageSizeYX must be a finite positive integer [Y X] vector.');
end

originXY = double(R.originXY_px(:).');
if numel(originXY) ~= 2 || any(~isfinite(originXY))
    error('createROIFile:InvalidOrigin', ...
        'originXY_px must be a finite [x y] coordinate.');
end

Ny = imageSizeYX(1);
Nx = imageSizeYX(2);

if originXY(1) < 1 || originXY(1) > Nx || ...
        originXY(2) < 1 || originXY(2) > Ny
    error('createROIFile:OriginOutOfBounds', ...
        'originXY_px must be inside image bounds: X [1 %d], Y [1 %d].', Nx, Ny);
end

pixelSize = iNormalizePixelSize(R.pixelSize_px_per_mm);

coordSystem = char(string(R.coordinateSystem));
if ~ismember(coordSystem, schema.allowedCoordinateSystems)
    error('createROIFile:InvalidCoordinateSystem', ...
        'Unsupported coordinate system: %s.', coordSystem);
end

viewMode = char(string(R.viewMode));
if ~ismember(viewMode, schema.allowedViewModes)
    error('createROIFile:InvalidViewMode', ...
        'Unsupported view mode: %s.', viewMode);
end

statsImage = R.statsImage;
if ~isempty(statsImage)
    if ~isequal(size(statsImage, 1), Ny) || ~isequal(size(statsImage, 2), Nx)
        error('createROIFile:StatsImageSizeMismatch', ...
            'statsImage size must match imageSizeYX.');
    end
end

tNow = datetime('now');

ROIFile = struct();
ROIFile.schemaName = schema.schemaName;
ROIFile.version = schema.version;
ROIFile.createdOn = tNow;
ROIFile.modifiedOn = tNow;

ROIFile.imageInfo = struct( ...
    'dataFile', char(string(R.dataFile)), ...
    'entryName', char(string(R.entryName)), ...
    'imageSizeYX', imageSizeYX, ...
    'coordinateSystem', coordSystem);

ROIFile.view = struct( ...
    'originXY_px', round(originXY), ...
    'pixelSize_px_per_mm', pixelSize);

if isempty(statsImage)
    statsCreatedOn = datetime.empty;
else
    statsCreatedOn = tNow;
end

ROIFile.statsImage = struct( ...
    'image', statsImage, ...
    'description', char(string(R.statsDescription)), ...
    'frameIndex', R.frameIndex, ...
    'viewMode', viewMode, ...
    'eventCondition', char(string(R.eventCondition)), ...
    'eventRepetition', char(string(R.eventRepetition)), ...
    'createdOn', statsCreatedOn);

ROIFile.ROIs = R.ROIs;

validateROIFile(ROIFile);

end

function pixelSize = iNormalizePixelSize(pixelSize)
%INORMALIZEPIXELSIZE Normalize []/0/positive scalar pixel size.

if isempty(pixelSize)
    pixelSize = [];
    return
end

pixelSize = double(pixelSize);

if ~isscalar(pixelSize) || ~isfinite(pixelSize) || pixelSize < 0
    error('createROIFile:InvalidPixelSize', ...
        'pixelSize_px_per_mm must be [] or a non-negative scalar.');
end

if pixelSize == 0
    pixelSize = [];
end

end
