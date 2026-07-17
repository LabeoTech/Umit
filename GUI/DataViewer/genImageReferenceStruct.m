function ImageReference = genImageReferenceStruct(imageData, varargin)
%GENIMAGEREFERENCESTRUCT Create a canonical UMIT ImageReference structure.
%
%   ImageReference = genImageReferenceStruct(imageData, ...
%       'Name', name, ...
%       'ProjectUUID', projectUUID, ...
%       'ProjectName', projectName, ...
%       'SubjectUUID', subjectUUID, ...
%       'SubjectID', subjectID)
%
%   ImageReference = genImageReferenceStruct(..., ...
%       'SessionUUID', sessionUUID, ...
%       'SessionID', sessionID, ...
%       'Description', description, ...
%       'SourceFolder', sourceFolder, ...
%       'SourceFile', sourceFile, ...
%       'SourceFrame', sourceFrame, ...
%       'SourceType', sourceType, ...
%       'ProcessingNotes', processingNotes, ...
%       'CreatedBy', createdBy, ...
%       'ParentResourceUUID', parentResourceUUID, ...
%       'WasManuallyTransformed', true, ...
%       'ManualTransformTform', tform, ...
%       'ManualTransformNote', note)
%
%   Creates an in-memory ImageReference payload. This function does not
%   choose a destination folder or save a managed project resource. Use
%   UMITProjectStore.addImageReference to import the resulting MAT-file.
%
%   Inputs:
%       imageData - Non-empty numeric or logical 2-D image.
%
%   Required name-value arguments:
%       Name        - User-facing reference name.
%       ProjectUUID - Immutable project UUID.
%       ProjectName - Project display name.
%       SubjectUUID - Immutable subject UUID.
%       SubjectID   - Current subject ID.
%
%   Optional name-value arguments:
%       SessionUUID            - Source session UUID. Default: ''
%       SessionID              - Source session ID. Default: ''
%       Description            - User description. Default: ''
%       SourceFolder           - External source folder. Default: ''
%       SourceFile             - External source file. Default: ''
%       SourceFrame            - Positive frame index or []. Default: []
%       SourceType             - Source type label. Default: ''
%       ProcessingNotes        - Processing notes. Default: ''
%       CreatedBy              - Creator identifier. Default: 'unknown'
%       ParentResourceUUID     - Parent ImageReference UUID. Default: ''
%       WasManuallyTransformed - Logical scalar. Default: false
%       ManualTransformTform   - Empty, affine2d, or numeric 3-by-3 matrix.
%       ManualTransformNote    - Transform provenance note. Default: ''
%       AxisConvention         - Coordinate convention.
%                                Default: 'imageXY_topLeft'
%       CreatedOn              - Creation datetime. Default: datetime('now')
%
%   Output:
%       ImageReference - Canonical scalar structure validated by
%                        validateImageReferenceStruct.
%
%   Notes:
%       - Images are stored as finite single-precision matrices.
%       - Non-finite image values are replaced by zero.
%       - Managed filenames, paths, checksums, and resource UUIDs are not
%         embedded in this payload. Those belong to UMITProjectStore.

errID = 'Umitoolbox:genImageReferenceStruct:invalidInput';

if ~(isnumeric(imageData) || islogical(imageData)) || ...
        isempty(imageData) || ~ismatrix(imageData)
    error(errID, ...
        '"imageData" must be a non-empty numeric or logical 2-D matrix.');
end

p = inputParser;
p.FunctionName = 'genImageReferenceStruct';

addParameter(p, 'Name', '', @iIsTextScalar);
addParameter(p, 'ProjectUUID', '', @iIsTextScalar);
addParameter(p, 'ProjectName', '', @iIsTextScalar);
addParameter(p, 'SubjectUUID', '', @iIsTextScalar);
addParameter(p, 'SubjectID', '', @iIsTextScalar);
addParameter(p, 'SessionUUID', '', @iIsTextScalar);
addParameter(p, 'SessionID', '', @iIsTextScalar);
addParameter(p, 'Description', '', @iIsTextScalar);
addParameter(p, 'SourceFolder', '', @iIsTextScalar);
addParameter(p, 'SourceFile', '', @iIsTextScalar);
addParameter(p, 'SourceFrame', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 1 && mod(x, 1) == 0));
addParameter(p, 'SourceType', '', @iIsTextScalar);
addParameter(p, 'ProcessingNotes', '', @iIsTextScalar);
addParameter(p, 'CreatedBy', 'unknown', @iIsTextScalar);
addParameter(p, 'ParentResourceUUID', '', @iIsTextScalar);
addParameter(p, 'WasManuallyTransformed', false, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'ManualTransformTform', [], @iIsTransformInput);
addParameter(p, 'ManualTransformNote', '', @iIsTextScalar);
addParameter(p, 'AxisConvention', 'imageXY_topLeft', @iIsTextScalar);
addParameter(p, 'CreatedOn', datetime('now'), ...
    @(x) isdatetime(x) && isscalar(x) && ~isnat(x));

parse(p, varargin{:});
R = p.Results;

requiredText = {'Name', 'ProjectUUID', 'ProjectName', ...
    'SubjectUUID', 'SubjectID'};
for iField = 1:numel(requiredText)
    fieldName = requiredText{iField};
    if strlength(strtrim(string(R.(fieldName)))) == 0
        error(errID, '"%s" cannot be empty.', fieldName);
    end
end

sessionUUID = char(string(R.SessionUUID));
sessionID = char(string(R.SessionID));
if xor(isempty(sessionUUID), isempty(sessionID))
    error(errID, ...
        '"SessionUUID" and "SessionID" must either both be empty or both be provided.');
end

manualTform = iNormalizeTransform(R.ManualTransformTform, errID);
if R.WasManuallyTransformed
    if isempty(manualTform)
        error(errID, ...
            '"ManualTransformTform" is required when WasManuallyTransformed=true.');
    end
    manualAppliedOn = R.CreatedOn;
else
    if ~isempty(manualTform)
        error(errID, ...
            ['"ManualTransformTform" must be empty when ' ...
             'WasManuallyTransformed=false.']);
    end
    manualAppliedOn = NaT;
end

imageData = single(imageData);
imageData(~isfinite(imageData)) = 0;

ImageReference = struct();
ImageReference.version = 1;
ImageReference.image = imageData;
ImageReference.name = char(string(R.Name));
ImageReference.description = char(string(R.Description));
ImageReference.projectUUID = lower(char(string(R.ProjectUUID)));
ImageReference.projectName = char(string(R.ProjectName));
ImageReference.subjectUUID = lower(char(string(R.SubjectUUID)));
ImageReference.subjectID = char(string(R.SubjectID));
ImageReference.sessionUUID = lower(sessionUUID);
ImageReference.sessionID = sessionID;
ImageReference.createdOn = R.CreatedOn;
ImageReference.lastModified = R.CreatedOn;
ImageReference.imageSizeYX = double(size(imageData));
ImageReference.axisConvention = char(string(R.AxisConvention));
ImageReference.sourceFolder = char(string(R.SourceFolder));
ImageReference.sourceFile = char(string(R.SourceFile));
ImageReference.sourceFrame = R.SourceFrame;
ImageReference.sourceType = char(string(R.SourceType));
ImageReference.processingNotes = char(string(R.ProcessingNotes));
ImageReference.createdBy = char(string(R.CreatedBy));
ImageReference.wasManuallyTransformed = R.WasManuallyTransformed;
ImageReference.parentResourceUUID = lower(char(string(R.ParentResourceUUID)));
ImageReference.manualTransformTform = manualTform;
ImageReference.manualTransformAppliedOn = manualAppliedOn;
ImageReference.manualTransformNote = char(string(R.ManualTransformNote));

validateImageReferenceStruct(ImageReference);

end

function tf = iIsTextScalar(x)
%IISTEXTSCALAR Return true for a character vector or string scalar.

tf = ischar(x) || (isstring(x) && isscalar(x));

end

function tf = iIsTransformInput(x)
%IISTRANSFORMINPUT Return true for a supported transform representation.

tf = isempty(x) || isa(x, 'affine2d') || ...
    (isnumeric(x) && isequal(size(x), [3 3]) && all(isfinite(x), 'all'));

end

function tform = iNormalizeTransform(value, errID)
%INORMALIZETRANSFORM Normalize transform input to affine2d.

if isempty(value)
    tform = [];
elseif isa(value, 'affine2d')
    tform = value;
elseif isnumeric(value) && isequal(size(value), [3 3]) && ...
        all(isfinite(value), 'all')
    tform = affine2d(double(value));
else
    error(errID, ...
        ['"ManualTransformTform" must be empty, affine2d, or a finite ' ...
         'numeric 3-by-3 matrix.']);
end

end

