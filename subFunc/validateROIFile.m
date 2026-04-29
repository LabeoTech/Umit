function validateROIFile(ROIFile)
%VALIDATEROIFILE Validate a UMIT ROI file structure.
%
%   validateROIFile(ROIFile)
%
%   Validates a versioned ROIFile structure against getROISchema(version).
%   This function intentionally does not support pre-version legacy ROI
%   structures. Retrocompatibility can be added later through migrateROIFile.

if ~isstruct(ROIFile) || ~isscalar(ROIFile)
    error('validateROIFile:InvalidROIFile', ...
        'ROIFile must be a scalar structure.');
end

if ~isfield(ROIFile, 'version')
    error('validateROIFile:MissingVersion', ...
        'ROIFile must contain a version field.');
end

schema = getROISchema(ROIFile.version);

iRequireFields(ROIFile, schema.requiredTopFields, 'ROIFile');

if ~strcmp(char(string(ROIFile.schemaName)), schema.schemaName)
    error('validateROIFile:InvalidSchemaName', ...
        'ROIFile.schemaName must be "%s".', schema.schemaName);
end

if ROIFile.version ~= schema.version
    error('validateROIFile:InvalidVersion', ...
        'ROIFile.version must be %g.', schema.version);
end

iValidateImageInfo(ROIFile.imageInfo, schema);
iValidateView(ROIFile.view);
iValidateStatsImage(ROIFile.statsImage, ROIFile.imageInfo, schema);
iValidateROIs(ROIFile.ROIs, ROIFile.imageInfo, schema);

end

function iValidateImageInfo(imageInfo, schema)
%IVALIDATEIMAGEINFO Validate ROIFile.imageInfo.

iRequireScalarStruct(imageInfo, 'ROIFile.imageInfo');
iRequireFields(imageInfo, schema.requiredImageInfoFields, 'ROIFile.imageInfo');

imageSizeYX = double(imageInfo.imageSizeYX(:).');
if numel(imageSizeYX) ~= 2 || any(~isfinite(imageSizeYX)) || ...
        any(imageSizeYX <= 0) || any(imageSizeYX ~= round(imageSizeYX))
    error('validateROIFile:InvalidImageSize', ...
        'ROIFile.imageInfo.imageSizeYX must be a finite positive integer [Y X] vector.');
end

coordSystem = char(string(imageInfo.coordinateSystem));
if ~ismember(coordSystem, schema.allowedCoordinateSystems)
    error('validateROIFile:InvalidCoordinateSystem', ...
        'Unsupported coordinate system: %s.', coordSystem);
end

end

function iValidateView(viewInfo)
%IVALIDATEVIEW Validate ROIFile.view.

iRequireScalarStruct(viewInfo, 'ROIFile.view');
iRequireFields(viewInfo, {'originXY_px', 'pixelSize_px_per_mm'}, 'ROIFile.view');

originXY = double(viewInfo.originXY_px(:).');
if numel(originXY) ~= 2 || any(~isfinite(originXY))
    error('validateROIFile:InvalidOrigin', ...
        'ROIFile.view.originXY_px must be a finite [x y] coordinate.');
end

px = viewInfo.pixelSize_px_per_mm;
if ~isempty(px)
    px = double(px);
    if ~isscalar(px) || ~isfinite(px) || px <= 0
        error('validateROIFile:InvalidPixelSize', ...
            'ROIFile.view.pixelSize_px_per_mm must be [] or a positive scalar.');
    end
end

end

function iValidateStatsImage(statsImage, imageInfo, schema)
%IVALIDATESTATSIMAGE Validate ROIFile.statsImage.

iRequireScalarStruct(statsImage, 'ROIFile.statsImage');
iRequireFields(statsImage, schema.requiredStatsImageFields, 'ROIFile.statsImage');

viewMode = char(string(statsImage.viewMode));
if ~ismember(viewMode, schema.allowedViewModes)
    error('validateROIFile:InvalidStatsImageViewMode', ...
        'Unsupported ROIFile.statsImage.viewMode: %s.', viewMode);
end

if ~isempty(statsImage.image)
    if ~isnumeric(statsImage.image) || ndims(statsImage.image) ~= 2 %#ok<ISMAT>
        error('validateROIFile:InvalidStatsImage', ...
            'ROIFile.statsImage.image must be empty or a 2D numeric array.');
    end

    imageSizeYX = double(imageInfo.imageSizeYX(:).');
    if ~isequal(size(statsImage.image, 1), imageSizeYX(1)) || ...
            ~isequal(size(statsImage.image, 2), imageSizeYX(2))
        error('validateROIFile:StatsImageSizeMismatch', ...
            'ROIFile.statsImage.image size must match imageInfo.imageSizeYX.');
    end
end

end

function iValidateROIs(ROIs, imageInfo, schema)
%IVALIDATEROIS Validate ROIFile.ROIs.

if ~isstruct(ROIs)
    error('validateROIFile:InvalidROIs', ...
        'ROIFile.ROIs must be a structure array.');
end

if isempty(ROIs)
    return
end

names = cell(1, numel(ROIs));
imageSizeYX = double(imageInfo.imageSizeYX(:).');

for iROI = 1:numel(ROIs)
    roi = ROIs(iROI);
    roiPath = sprintf('ROIFile.ROIs(%d)', iROI);

    iRequireFields(roi, schema.requiredROIFields, roiPath);

    roiName = char(string(roi.name));
    if isempty(strtrim(roiName))
        error('validateROIFile:InvalidROIName', ...
            '%s.name must be non-empty.', roiPath);
    end
    names{iROI} = roiName;

    roiType = lower(char(string(roi.type)));
    if ~ismember(roiType, schema.allowedROITypes)
        error('validateROIFile:InvalidROIType', ...
            '%s.type is unsupported: %s.', roiPath, roiType);
    end

    iValidateRGB(roi.color, [roiPath, '.color']);
    iValidateGeometry(roi.geometry, schema, roiPath);
    iValidateMask(roi.mask, imageSizeYX, roiPath);
    iValidateStats(roi.stats, schema, roiPath);
end

if numel(unique(names, 'stable')) ~= numel(names)
    error('validateROIFile:DuplicateROINames', ...
        'ROI names must be unique within ROIFile.ROIs.');
end

end

function iValidateGeometry(geometry, schema, roiPath)
%IVALIDATEGEOMETRY Validate ROI geometry.

iRequireScalarStruct(geometry, [roiPath, '.geometry']);
iRequireFields(geometry, schema.requiredGeometryFields, [roiPath, '.geometry']);

if ~isa(geometry.polyshape, 'polyshape')
    error('validateROIFile:InvalidPolyshape', ...
        '%s.geometry.polyshape must be a polyshape object.', roiPath);
end

if ~isempty(geometry.verticesXY_px)
    if ~isnumeric(geometry.verticesXY_px) || size(geometry.verticesXY_px, 2) ~= 2 || ...
            any(~isfinite(geometry.verticesXY_px(:)))
        error('validateROIFile:InvalidVertices', ...
            '%s.geometry.verticesXY_px must be empty or finite [N x 2].', roiPath);
    end
end

roiType = lower(char(string(geometry.ROIType)));
if ~ismember(roiType, schema.allowedROITypes)
    error('validateROIFile:InvalidGeometryROIType', ...
        '%s.geometry.ROIType is unsupported: %s.', roiPath, roiType);
end

if ~isstruct(geometry.ROIParameters)
    error('validateROIFile:InvalidROIParameters', ...
        '%s.geometry.ROIParameters must be a structure.', roiPath);
end

end

function iValidateMask(mask, imageSizeYX, roiPath)
%IVALIDATEMASK Validate ROI mask.

if ~(islogical(mask) || isnumeric(mask)) || ndims(mask) ~= 2 %#ok<ISMAT>
    error('validateROIFile:InvalidMask', ...
        '%s.mask must be a 2D logical or numeric array.', roiPath);
end

if ~isequal(size(mask, 1), imageSizeYX(1)) || ~isequal(size(mask, 2), imageSizeYX(2))
    error('validateROIFile:MaskSizeMismatch', ...
        '%s.mask size must match ROIFile.imageInfo.imageSizeYX.', roiPath);
end

end

function iValidateStats(stats, schema, roiPath)
%IVALIDATESTATS Validate ROI stats structure.

iRequireScalarStruct(stats, [roiPath, '.stats']);
iRequireFields(stats, schema.requiredStatsFields, [roiPath, '.stats']);

numericFields = setdiff(schema.requiredStatsFields, {'computedOn'});
for iField = 1:numel(numericFields)
    fieldName = numericFields{iField};
    value = stats.(fieldName);

    if ~isempty(value) && ~isnumeric(value)
        error('validateROIFile:InvalidStatsField', ...
            '%s.stats.%s must be numeric or empty.', roiPath, fieldName);
    end
end

end

function iValidateRGB(value, fieldPath)
%IVALIDATERGB Validate RGB triplet.

value = double(value(:).');
if numel(value) ~= 3 || any(~isfinite(value)) || any(value < 0) || any(value > 1)
    error('validateROIFile:InvalidRGB', ...
        '%s must be an RGB triplet with values between 0 and 1.', fieldPath);
end

end

function iRequireScalarStruct(value, fieldPath)
%IREQUIRESCALARSTRUCT Require scalar struct.

if ~isstruct(value) || ~isscalar(value)
    error('validateROIFile:InvalidStruct', ...
        '%s must be a scalar structure.', fieldPath);
end

end

function iRequireFields(S, fields, structName)
%IREQUIREFIELDS Require fields in struct.

missing = fields(~isfield(S, fields));

if ~isempty(missing)
    error('validateROIFile:MissingField', ...
        '%s is missing required field(s): %s.', ...
        structName, strjoin(missing, ', '));
end

end
