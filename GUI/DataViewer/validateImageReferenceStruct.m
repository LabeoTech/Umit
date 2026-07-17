function validateImageReferenceStruct(ImageReference)
%VALIDATEIMAGEREFERENCESTRUCT Validate a canonical UMIT ImageReference.
%
%   validateImageReferenceStruct(ImageReference)
%
%   Validates the in-memory payload stored in managed Image Reference MAT
%   files. Managed resource identity, path, checksum, status, and filename
%   are intentionally excluded because they belong to UMITProjectStore.
%
%   Required fields:
%       version
%       image
%       name
%       description
%       projectUUID
%       projectName
%       subjectUUID
%       subjectID
%       sessionUUID
%       sessionID
%       createdOn
%       lastModified
%       imageSizeYX
%       axisConvention
%       sourceFolder
%       sourceFile
%       sourceFrame
%       sourceType
%       processingNotes
%       createdBy
%       wasManuallyTransformed
%       parentResourceUUID
%       manualTransformTform
%       manualTransformAppliedOn
%       manualTransformNote
%
%   Error:
%       Throws Umitoolbox:validateImageReferenceStruct:invalidInput when
%       validation fails.

errID = 'Umitoolbox:validateImageReferenceStruct:invalidInput';

if ~isstruct(ImageReference) || ~isscalar(ImageReference)
    error(errID, 'ImageReference must be a scalar struct.');
end

requiredFields = { ...
    'version', ...
    'image', ...
    'name', ...
    'description', ...
    'projectUUID', ...
    'projectName', ...
    'subjectUUID', ...
    'subjectID', ...
    'sessionUUID', ...
    'sessionID', ...
    'createdOn', ...
    'lastModified', ...
    'imageSizeYX', ...
    'axisConvention', ...
    'sourceFolder', ...
    'sourceFile', ...
    'sourceFrame', ...
    'sourceType', ...
    'processingNotes', ...
    'createdBy', ...
    'wasManuallyTransformed', ...
    'parentResourceUUID', ...
    'manualTransformTform', ...
    'manualTransformAppliedOn', ...
    'manualTransformNote'};

fields = fieldnames(ImageReference);
missingFields = setdiff(requiredFields, fields);
if ~isempty(missingFields)
    error(errID, ...
        'ImageReference is missing required field(s): %s.', ...
        strjoin(missingFields, ', '));
end

unknownFields = setdiff(fields, requiredFields);
if ~isempty(unknownFields)
    error(errID, ...
        'ImageReference contains unsupported field(s): %s.', ...
        strjoin(unknownFields, ', '));
end

if ~isnumeric(ImageReference.version) || ...
        ~isscalar(ImageReference.version) || ...
        ImageReference.version ~= 1
    error(errID, '"version" must be numeric scalar 1.');
end

imageData = ImageReference.image;
if ~(isnumeric(imageData) || islogical(imageData)) || ...
        isempty(imageData) || ~ismatrix(imageData)
    error(errID, ...
        '"image" must be a non-empty numeric or logical 2-D matrix.');
end
if any(~isfinite(imageData), 'all')
    error(errID, '"image" must contain only finite values.');
end

imageSizeYX = ImageReference.imageSizeYX;
if ~isnumeric(imageSizeYX) || ~isvector(imageSizeYX) || ...
        numel(imageSizeYX) ~= 2 || any(~isfinite(imageSizeYX)) || ...
        any(imageSizeYX < 1) || any(mod(imageSizeYX, 1) ~= 0)
    error(errID, ...
        '"imageSizeYX" must be a two-element vector of positive integers.');
end
if ~isequal(double(imageSizeYX(:).'), double(size(imageData)))
    error(errID, '"imageSizeYX" does not match the stored image size.');
end

requiredNonemptyText = { ...
    'name', 'projectUUID', 'projectName', 'subjectUUID', ...
    'subjectID', 'axisConvention', 'createdBy'};
optionalText = { ...
    'description', 'sessionUUID', 'sessionID', 'sourceFolder', ...
    'sourceFile', 'sourceType', 'processingNotes', ...
    'parentResourceUUID', 'manualTransformNote'};

for iField = 1:numel(requiredNonemptyText)
    fieldName = requiredNonemptyText{iField};
    value = ImageReference.(fieldName);
    if ~iIsTextScalar(value) || strlength(strtrim(string(value))) == 0
        error(errID, '"%s" must be a non-empty text scalar.', fieldName);
    end
end

for iField = 1:numel(optionalText)
    fieldName = optionalText{iField};
    if ~iIsTextScalar(ImageReference.(fieldName))
        error(errID, '"%s" must be a text scalar.', fieldName);
    end
end

uuidFields = {'projectUUID', 'subjectUUID'};
for iField = 1:numel(uuidFields)
    fieldName = uuidFields{iField};
    iValidateUUID(ImageReference.(fieldName), fieldName, false, errID);
end
iValidateUUID( ...
    ImageReference.sessionUUID, 'sessionUUID', true, errID);
iValidateUUID( ...
    ImageReference.parentResourceUUID, 'parentResourceUUID', true, errID);

sessionUUID = char(string(ImageReference.sessionUUID));
sessionID = char(string(ImageReference.sessionID));
if xor(isempty(sessionUUID), isempty(sessionID))
    error(errID, ...
        '"sessionUUID" and "sessionID" must either both be empty or both be populated.');
end

dateFields = {'createdOn', 'lastModified'};
for iField = 1:numel(dateFields)
    fieldName = dateFields{iField};
    value = ImageReference.(fieldName);
    if ~isdatetime(value) || ~isscalar(value) || isnat(value)
        error(errID, '"%s" must be a non-NaT datetime scalar.', fieldName);
    end
end
if ImageReference.lastModified < ImageReference.createdOn
    error(errID, '"lastModified" cannot precede "createdOn".');
end

sourceFrame = ImageReference.sourceFrame;
if ~(isempty(sourceFrame) || ...
        (isnumeric(sourceFrame) && isscalar(sourceFrame) && ...
        isfinite(sourceFrame) && sourceFrame >= 1 && ...
        mod(sourceFrame, 1) == 0))
    error(errID, ...
        '"sourceFrame" must be empty or a positive integer scalar.');
end

if ~islogical(ImageReference.wasManuallyTransformed) || ...
        ~isscalar(ImageReference.wasManuallyTransformed)
    error(errID, ...
        '"wasManuallyTransformed" must be a logical scalar.');
end

manualTform = ImageReference.manualTransformTform;
manualAppliedOn = ImageReference.manualTransformAppliedOn;

if ImageReference.wasManuallyTransformed
    if ~isa(manualTform, 'affine2d')
        error(errID, ...
            ['"manualTransformTform" must be an affine2d when ' ...
             'wasManuallyTransformed=true.']);
    end
    if ~isdatetime(manualAppliedOn) || ...
            ~isscalar(manualAppliedOn) || isnat(manualAppliedOn)
        error(errID, ...
            ['"manualTransformAppliedOn" must be a non-NaT datetime ' ...
             'when wasManuallyTransformed=true.']);
    end
else
    if ~isempty(manualTform)
        error(errID, ...
            ['"manualTransformTform" must be empty when ' ...
             'wasManuallyTransformed=false.']);
    end
    if ~isdatetime(manualAppliedOn) || ...
            ~isscalar(manualAppliedOn) || ~isnat(manualAppliedOn)
        error(errID, ...
            ['"manualTransformAppliedOn" must be scalar NaT when ' ...
             'wasManuallyTransformed=false.']);
    end
end

end

function tf = iIsTextScalar(x)
%IISTEXTSCALAR Return true for a character vector or string scalar.

tf = ischar(x) || (isstring(x) && isscalar(x));

end

function iValidateUUID(value, fieldName, allowEmpty, errID)
%IVALIDATEUUID Validate one canonical UUID text field.

if ~iIsTextScalar(value)
    error(errID, '"%s" must be a text scalar.', fieldName);
end

value = char(string(value));
if isempty(value)
    if allowEmpty
        return
    end
    error(errID, '"%s" cannot be empty.', fieldName);
end

uuidPattern = [ ...
    '^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-' ...
    '[1-5][0-9a-fA-F]{3}-' ...
    '[89abAB][0-9a-fA-F]{3}-' ...
    '[0-9a-fA-F]{12}$'];

if isempty(regexp(value, uuidPattern, 'once'))
    error(errID, '"%s" must contain a valid UUID.', fieldName);
end

end
