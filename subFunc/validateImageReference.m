function validateImageReference(ImageReference)
%VALIDATEIMAGEREFERENCE Validate UMIT ImageReference structure.
%
%   validateImageReference(ImageReference)
%
%   Checks whether ImageReference follows the expected flat schema.
%
%   Required fields:
%       schemaVersion
%       createdOn
%       lastModified
%       notes
%       image
%       imageSizeYX
%       axisConvention
%       name
%       description
%       projectName
%       groupName
%       mouseID
%       sessionID
%       sourceFolder
%       sourceFile
%       sourceFrame
%       sourceType
%       createdBy
%       wasManuallyTransformed
%       processingNotes
%       relativePath
%       fileName
%       custom
%
%   Optional provenance fields:
%       manualTransformAppliedOn
%       manualTransformType
%       manualTransformNote
%       manualTransformTform
%       lastModifiedBy

if ~isstruct(ImageReference) || ~isscalar(ImageReference)
    error('validateImageReference:InvalidType', ...
        'ImageReference must be a scalar struct.');
end

reqFields = { ...
    'schemaVersion', ...
    'createdOn', ...
    'lastModified', ...
    'notes', ...
    'image', ...
    'imageSizeYX', ...
    'axisConvention', ...
    'name', ...
    'description', ...
    'projectName', ...
    'groupName', ...
    'mouseID', ...
    'sessionID', ...
    'sourceFolder', ...
    'sourceFile', ...
    'sourceFrame', ...
    'sourceType', ...
    'createdBy', ...
    'wasManuallyTransformed', ...
    'processingNotes', ...
    'relativePath', ...
    'fileName', ...
    'custom'};

for k = 1:numel(reqFields)
    if ~isfield(ImageReference, reqFields{k})
        error('validateImageReference:MissingField', ...
            'Missing field: ImageReference.%s', reqFields{k});
    end
end

% -------------------------------------------------------------------------
% Image and spatial metadata
% -------------------------------------------------------------------------
refImg = ImageReference.image;

if ~(isnumeric(refImg) || islogical(refImg)) || ndims(refImg) ~= 2 || isempty(refImg)
    error('validateImageReference:InvalidImage', ...
        'ImageReference.image must be a non-empty 2D numeric or logical array.');
end

imageSizeYX = ImageReference.imageSizeYX;

if ~isnumeric(imageSizeYX) || numel(imageSizeYX) ~= 2 || ...
        any(~isfinite(imageSizeYX(:))) || any(imageSizeYX(:) <= 0) || ...
        any(mod(imageSizeYX(:), 1) ~= 0)
    error('validateImageReference:InvalidImageSize', ...
        'ImageReference.imageSizeYX must be a positive integer [Ny Nx] vector.');
end

if ~isequal(size(refImg), double(imageSizeYX(:).'))
    error('validateImageReference:ImageSizeMismatch', ...
        'ImageReference.image size must match ImageReference.imageSizeYX.');
end

if ~iIsTextScalar(ImageReference.axisConvention)
    error('validateImageReference:InvalidAxisConvention', ...
        'ImageReference.axisConvention must be a char vector or string scalar.');
end

% -------------------------------------------------------------------------
% Required non-empty identifying metadata
% -------------------------------------------------------------------------
nonEmptyTextFields = { ...
    'projectName', ...
    'groupName', ...
    'mouseID', ...
    'sessionID'};

for k = 1:numel(nonEmptyTextFields)
    fieldName = nonEmptyTextFields{k};
    val = ImageReference.(fieldName);

    if ~iIsTextScalar(val) || strlength(strtrim(string(val))) == 0
        error('validateImageReference:InvalidRequiredText', ...
            'ImageReference.%s must be non-empty text.', fieldName);
    end
end

% -------------------------------------------------------------------------
% Text metadata
% -------------------------------------------------------------------------
textFields = { ...
    'schemaVersion', ...
    'notes', ...
    'name', ...
    'description', ...
    'sourceFolder', ...
    'sourceFile', ...
    'sourceType', ...
    'createdBy', ...
    'processingNotes', ...
    'relativePath', ...
    'fileName'};

optionalTextFields = { ...
    'manualTransformType', ...
    'manualTransformNote', ...
    'lastModifiedBy'};

textFields = [textFields, optionalTextFields(isfield(ImageReference, optionalTextFields))];

for k = 1:numel(textFields)
    fieldName = textFields{k};
    val = ImageReference.(fieldName);

    if ~(isempty(val) || iIsTextScalar(val))
        error('validateImageReference:InvalidTextField', ...
            'ImageReference.%s must be empty, char, or string scalar.', fieldName);
    end
end

if strlength(strtrim(string(ImageReference.relativePath))) == 0
    error('validateImageReference:MissingRelativePath', ...
        'ImageReference.relativePath cannot be empty.');
end

if strlength(strtrim(string(ImageReference.fileName))) == 0
    error('validateImageReference:MissingFileName', ...
        'ImageReference.fileName cannot be empty.');
end

% -------------------------------------------------------------------------
% Dates and flags
% -------------------------------------------------------------------------
if ~(isdatetime(ImageReference.createdOn) && isscalar(ImageReference.createdOn))
    error('validateImageReference:InvalidCreatedOn', ...
        'ImageReference.createdOn must be a datetime scalar.');
end

if ~(isdatetime(ImageReference.lastModified) && isscalar(ImageReference.lastModified))
    error('validateImageReference:InvalidLastModified', ...
        'ImageReference.lastModified must be a datetime scalar.');
end

if isfield(ImageReference, 'manualTransformAppliedOn') && ...
        ~isempty(ImageReference.manualTransformAppliedOn)
    if ~(isdatetime(ImageReference.manualTransformAppliedOn) && ...
            isscalar(ImageReference.manualTransformAppliedOn))
        error('validateImageReference:InvalidManualTransformAppliedOn', ...
            'ImageReference.manualTransformAppliedOn must be empty or a datetime scalar.');
    end
end

if ~islogical(ImageReference.wasManuallyTransformed) || ...
        ~isscalar(ImageReference.wasManuallyTransformed)
    error('validateImageReference:InvalidManualTransformFlag', ...
        'ImageReference.wasManuallyTransformed must be a scalar logical.');
end

if ~isempty(ImageReference.sourceFrame)
    if ~(isnumeric(ImageReference.sourceFrame) && isscalar(ImageReference.sourceFrame) && ...
            isfinite(ImageReference.sourceFrame) && ImageReference.sourceFrame >= 1)
        error('validateImageReference:InvalidSourceFrame', ...
            'ImageReference.sourceFrame must be empty or a positive numeric scalar.');
    end
end

if isfield(ImageReference, 'manualTransformTform') && ...
        ~isempty(ImageReference.manualTransformTform)
    if ~iIsSupportedTform2D(ImageReference.manualTransformTform)
        error('validateImageReference:InvalidManualTransformTform', ...
            'ImageReference.manualTransformTform must be empty or a supported 2D geometric transform object.');
    end
end

if ~isstruct(ImageReference.custom) || ~isscalar(ImageReference.custom)
    error('validateImageReference:InvalidCustom', ...
        'ImageReference.custom must be a scalar struct.');
end

end

function tf = iIsTextScalar(x)
%IISTEXTSCALAR Return true for char vectors or string scalars.

tf = ischar(x) || (isstring(x) && isscalar(x));
end

function tf = iIsSupportedTform2D(tform)
%IISSUPPORTEDTFORM2D Return true for common MATLAB 2D geometric transforms.

tf = false;

if isempty(tform)
    tf = true;
    return
end

validClasses = { ...
    'affine2d', ...
    'projective2d', ...
    'rigidtform2d', ...
    'simtform2d', ...
    'affinetform2d', ...
    'projtfom2d'};

for iClass = 1:numel(validClasses)
    try
        if isa(tform, validClasses{iClass})
            tf = true;
            return
        end
    catch
    end
end

% Fallback for release/version-specific transform classes.
className = lower(class(tform));
tf = contains(className, 'tform2d') || contains(className, 'affine2d') || ...
    contains(className, 'projective2d');
end