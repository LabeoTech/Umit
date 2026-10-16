function validateImageReferenceStruct(ImageReference)
%VALIDATEIMAGEREFERENCESTRUCT Validate a canonical UMIT ImageReference.
%
%   validateImageReferenceStruct(ImageReference)
%
%   Validates the in-memory payload stored in Image Reference MAT files
%   (subject libraries and standalone exported files). Managed resource
%   identity, path, checksum, status, and filename are intentionally
%   excluded because they belong to UMITProjectStore.
%
%   Versions:
%       2 (written today): identity is OPTIONAL. projectUUID/projectName and
%         subjectUUID/subjectID are both-or-neither pairs, and a subject
%         needs a project. contentHash must equal
%         computeImageContentHash(image).
%       1 (read only): the original contract; project and subject identity
%         are required and there is no contentHash. Use
%         normalizeImageReferenceStruct to upgrade a v1 struct in memory.
%
%   Required fields (v2 adds contentHash):
%       version
%       image
%       contentHash (v2)
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

if ~isfield(ImageReference, 'version') || ~isnumeric(ImageReference.version) || ...
        ~isscalar(ImageReference.version) || ~ismember(ImageReference.version, [1 2])
    error(errID, '"version" must be numeric scalar 1 or 2.');
end
version = double(ImageReference.version);
if version == 2
    requiredFields = [requiredFields(1:2), {'contentHash'}, requiredFields(3:end)];
end

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

identityText = {'projectUUID', 'projectName', 'subjectUUID', 'subjectID'};
requiredNonemptyText = {'name', 'axisConvention', 'createdBy'};
optionalText = { ...
    'description', 'sessionUUID', 'sessionID', 'sourceFolder', ...
    'sourceFile', 'sourceType', 'processingNotes', ...
    'parentResourceUUID', 'manualTransformNote'};
if version == 1
    requiredNonemptyText = [requiredNonemptyText, identityText];
else
    optionalText = [optionalText, identityText];
end

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
    iValidateUUID(ImageReference.(fieldName), fieldName, version == 2, errID);
end
if version == 2
    iValidatePair(ImageReference, 'projectUUID', 'projectName', errID);
    iValidatePair(ImageReference, 'subjectUUID', 'subjectID', errID);
    if ~isempty(char(string(ImageReference.subjectUUID))) && ...
            isempty(char(string(ImageReference.projectUUID)))
        error(errID, 'A subject identity requires a project identity.');
    end
    contentHash = ImageReference.contentHash;
    if ~iIsTextScalar(contentHash) || ...
            isempty(regexp(char(string(contentHash)), '^sha256:[0-9a-f]{64}$', 'once'))
        error(errID, '"contentHash" must be ''sha256:'' followed by 64 lowercase hex digits.');
    end
    if ~strcmp(char(string(contentHash)), computeImageContentHash(imageData))
        error(errID, '"contentHash" does not match the stored image.');
    end
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

function iValidatePair(ImageReference, firstField, secondField, errID)
%IVALIDATEPAIR Both-or-neither identity fields (version 2).

first = strtrim(char(string(ImageReference.(firstField))));
second = strtrim(char(string(ImageReference.(secondField))));
if xor(isempty(first), isempty(second))
    error(errID, '"%s" and "%s" must either both be empty or both be populated.', ...
        firstField, secondField);
end

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
