function validateDataParams(DataParams)
%VALIDATEDATAPARAMS Validate structure and field consistency of DataParams.
%
%   validateDataParams(DataParams)
%
%   Checks whether DataParams follows the expected schema and verifies
%   consistency of folder-global fields such as image size, origin, mask,
%   and registration metadata.
%
%   Required top-level fields:
%       schemaVersion
%       createdOn
%       lastModified
%       notes
%       view
%       mask
%       registration
%       custom
%
%   Notes:
%       - This function throws an error when validation fails.
%       - Validation is intentionally permissive for optional fields.

if ~isstruct(DataParams) || ~isscalar(DataParams)
    error('validateDataParams:InvalidType', ...
        'DataParams must be a scalar struct.');
end

reqTopFields = { ...
    'schemaVersion', ...
    'createdOn', ...
    'lastModified', ...
    'notes', ...
    'view', ...
    'mask', ...
    'registration', ...
    'custom'};

for k = 1:numel(reqTopFields)
    if ~isfield(DataParams, reqTopFields{k})
        error('validateDataParams:MissingField', ...
            'Missing top-level field: DataParams.%s', reqTopFields{k});
    end
end

if ~isstruct(DataParams.view) || ~isscalar(DataParams.view)
    error('validateDataParams:InvalidView', ...
        'DataParams.view must be a scalar struct.');
end

reqViewFields = { ...
    'pixelSize_px_per_mm', ...
    'origin_xy_px', ...
    'imageSizeYX', ...
    'axisConvention'};

for k = 1:numel(reqViewFields)
    if ~isfield(DataParams.view, reqViewFields{k})
        error('validateDataParams:MissingViewField', ...
            'Missing field: DataParams.view.%s', reqViewFields{k});
    end
end

if ~isempty(DataParams.view.pixelSize_px_per_mm)
    x = DataParams.view.pixelSize_px_per_mm;
    if ~isnumeric(x) || ~(isscalar(x) || numel(x) == 2) || any(~isfinite(x(:))) || any(x(:) <= 0)
        error('validateDataParams:InvalidPixelSize', ...
            'DataParams.view.pixelSize_px_per_mm must be a positive scalar or 2-element numeric vector.');
    end
end

if ~isempty(DataParams.view.origin_xy_px)
    x = DataParams.view.origin_xy_px;
    if ~isnumeric(x) || numel(x) ~= 2 || any(~isfinite(x(:)))
        error('validateDataParams:InvalidOrigin', ...
            'DataParams.view.origin_xy_px must be a 2-element numeric vector [x y].');
    end
end

if ~isempty(DataParams.view.imageSizeYX)
    x = DataParams.view.imageSizeYX;
    if ~isnumeric(x) || numel(x) ~= 2 || any(~isfinite(x(:))) || any(x(:) <= 0) || any(mod(x(:),1) ~= 0)
        error('validateDataParams:InvalidImageSize', ...
            'DataParams.view.imageSizeYX must be a positive integer 2-element vector [Ny Nx].');
    end
end

if ~(ischar(DataParams.view.axisConvention) || isstring(DataParams.view.axisConvention))
    error('validateDataParams:InvalidAxisConvention', ...
        'DataParams.view.axisConvention must be a char vector or string scalar.');
end

if ~isstruct(DataParams.mask) || ~isscalar(DataParams.mask)
    error('validateDataParams:InvalidMask', ...
        'DataParams.mask must be a scalar struct.');
end

reqMaskFields = { ...
    'logical', ...
    'name', ...
    'description', ...
    'space', ...
    'createdOn', ...
    'source'};

for k = 1:numel(reqMaskFields)
    if ~isfield(DataParams.mask, reqMaskFields{k})
        error('validateDataParams:MissingMaskField', ...
            'Missing field: DataParams.mask.%s', reqMaskFields{k});
    end
end

if ~isempty(DataParams.mask.logical)
    if ~islogical(DataParams.mask.logical) || ndims(DataParams.mask.logical) ~= 2
        error('validateDataParams:InvalidMaskLogical', ...
            'DataParams.mask.logical must be a 2D logical array.');
    end

    if ~isempty(DataParams.view.imageSizeYX)
        if ~isequal(size(DataParams.mask.logical), DataParams.view.imageSizeYX)
            error('validateDataParams:MaskSizeMismatch', ...
                'DataParams.mask.logical size must match DataParams.view.imageSizeYX.');
        end
    end
end

if ~isstruct(DataParams.registration) || ~isscalar(DataParams.registration)
    error('validateDataParams:InvalidRegistration', ...
        'DataParams.registration must be a scalar struct.');
end

reqRegFields = { ...
    'isRegistered', ...
    'tform', ...
    'transformType', ...
    'method', ...
    'referenceDescription', ...
    'referenceFile', ...
    'referenceImage', ...
    'createdOn', ...
    'source'};

for k = 1:numel(reqRegFields)
    if ~isfield(DataParams.registration, reqRegFields{k})
        error('validateDataParams:MissingRegistrationField', ...
            'Missing field: DataParams.registration.%s', reqRegFields{k});
    end
end

if ~islogical(DataParams.registration.isRegistered) || ~isscalar(DataParams.registration.isRegistered)
    error('validateDataParams:InvalidIsRegistered', ...
        'DataParams.registration.isRegistered must be a scalar logical.');
end

if DataParams.registration.isRegistered && isempty(DataParams.registration.tform)
    error('validateDataParams:MissingTform', ...
        'DataParams.registration.tform cannot be empty when isRegistered is true.');
end

if ~(isempty(DataParams.registration.referenceFile) || ...
        ischar(DataParams.registration.referenceFile) || ...
        isstring(DataParams.registration.referenceFile))
    error('validateDataParams:InvalidReferenceFile', ...
        'DataParams.registration.referenceFile must be empty, a char vector, or a string scalar.');
end

if ~isempty(DataParams.registration.referenceImage)
    if ~isnumeric(DataParams.registration.referenceImage) && ...
            ~islogical(DataParams.registration.referenceImage)
        error('validateDataParams:InvalidReferenceImage', ...
            'DataParams.registration.referenceImage must be numeric or logical.');
    end

    if ndims(DataParams.registration.referenceImage) ~= 2
        error('validateDataParams:InvalidReferenceImageDims', ...
            'DataParams.registration.referenceImage must be a 2D array.');
    end

    if ~isempty(DataParams.view.imageSizeYX)
        if ~isequal(size(DataParams.registration.referenceImage), DataParams.view.imageSizeYX)
            error('validateDataParams:ReferenceImageSizeMismatch', ...
                ['DataParams.registration.referenceImage size must match ' ...
                 'DataParams.view.imageSizeYX.']);
        end
    end
end

if ~isempty(DataParams.registration.tform)
    tf = DataParams.registration.tform;
    isValidTform = isa(tf, 'affine2d') || ...
                   isa(tf, 'projective2d') || ...
                   isa(tf, 'rigidtform2d') || ...
                   isa(tf, 'simtform2d') || ...
                   isa(tf, 'affinetform2d') || ...
                   isa(tf, 'images.geotrans.GeometricTransformation2D');

    if ~isValidTform
        error('validateDataParams:InvalidTform', ...
            'DataParams.registration.tform is not a recognized 2D geometric transform object.');
    end
end
end