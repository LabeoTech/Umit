function [ImageReference, referenceFileAbs, referenceFileRel] = createImageReference(referenceImage, opts)
%CREATEIMAGEREFERENCE Create and save a reusable UMIT image reference file.
%
%   [ImageReference, referenceFileAbs, referenceFileRel] = ...
%       createImageReference(referenceImage, ...
%           'projectName', projectName, ...
%           'groupName', groupName, ...
%           'mouseID', mouseID, ...
%           'sessionID', sessionID)
%
%   Creates an ImageReference structure and saves it under:
%
%       getUmitFolder('referenceImages/ProjectName/GroupName/MouseID/SessionID')
%
%   The saved file is named:
%
%       ImageReference_<yyyyMMdd_HHmmss>.mat
%
%   Inputs:
%       referenceImage - 2D numeric/logical image used as target reference.
%
%   Name-Value options:
%       projectName            - Project folder/name.
%       groupName              - Group folder/name.
%       mouseID                - Mouse ID.
%       sessionID              - Session ID.
%       referenceName          - User-facing reference name.
%       description            - Free-text reference description.
%       sourceFolder           - Folder from which image was generated.
%       sourceFile             - Source file path/name, if applicable.
%       sourceFrame            - Source frame index, if applicable.
%       sourceType             - Source type: current frame, dat, workspace, etc.
%       wasManuallyTransformed - Whether the image was manually transformed
%                                before saving.
%       processingNotes        - Notes about preprocessing/manual transform.
%       axisConvention         - Coordinate convention. Default:
%                                'imageXY_topLeft'
%       createdBy              - Source app/function name.
%
%   Outputs:
%       ImageReference   - Saved ImageReference structure.
%       referenceFileAbs - Absolute path to saved .mat file.
%       referenceFileRel - Path relative to the UMIT root folder.
%
%   Notes:
%       ImageReference is intentionally flat. It is a reusable reference
%       artifact, not a folder-level DataParams object.
%
%       manualTransformTform is provenance only. It must not be confused
%       with DataParams.registration.tform, which is the transform used to
%       align a dataset to this reference.

arguments
    referenceImage {mustBeNumericOrLogical2D}
    opts.projectName (1,1) string = ""
    opts.groupName (1,1) string = ""
    opts.mouseID (1,1) string = ""
    opts.sessionID (1,1) string = ""
    opts.referenceName (1,1) string = ""
    opts.description (1,1) string = ""
    opts.sourceFolder (1,1) string = ""
    opts.sourceFile (1,1) string = ""
    opts.sourceFrame = []
    opts.sourceType (1,1) string = ""
    opts.wasManuallyTransformed (1,1) logical = false
    opts.processingNotes (1,1) string = ""
    opts.axisConvention (1,1) string = "imageXY_topLeft"
    opts.createdBy (1,1) string = "DataViewer Set Image Reference"
end

projectName = iRequireText(opts.projectName, 'projectName');
groupName = iRequireText(opts.groupName, 'groupName');
mouseID = iRequireText(opts.mouseID, 'mouseID');
sessionID = iRequireText(opts.sessionID, 'sessionID');

referenceName = strtrim(opts.referenceName);
if referenceName == ""
    referenceName = "ImageReference";
end

projectFolder = iSanitizePathPart(projectName);
groupFolder = iSanitizePathPart(groupName);
mouseFolder = iSanitizePathPart(mouseID);
sessionFolder = iSanitizePathPart(sessionID);

referenceSubFolder = fullfile( ...
    'referenceImages', ...
    projectFolder, ...
    groupFolder, ...
    mouseFolder, ...
    sessionFolder);

referenceFolderAbs = getUmitFolder(referenceSubFolder);

dateTag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
referenceFileName = ['ImageReference_' dateTag '.mat'];

referenceFileAbs = fullfile(referenceFolderAbs, referenceFileName);
referenceFileRel = fullfile(referenceSubFolder, referenceFileName);

if ~isa(referenceImage, 'single')
    referenceImage = single(referenceImage);
end

tNow = datetime('now');

ImageReference = struct();

ImageReference.schemaVersion = '1.0';
ImageReference.createdOn = tNow;
ImageReference.lastModified = tNow;
ImageReference.notes = "";

ImageReference.image = referenceImage;
ImageReference.imageSizeYX = size(referenceImage);
ImageReference.axisConvention = char(opts.axisConvention);

ImageReference.name = char(referenceName);
ImageReference.description = char(opts.description);

ImageReference.projectName = char(projectName);
ImageReference.groupName = char(groupName);
ImageReference.mouseID = char(mouseID);
ImageReference.sessionID = char(sessionID);

ImageReference.sourceFolder = char(opts.sourceFolder);
ImageReference.sourceFile = char(opts.sourceFile);
ImageReference.sourceFrame = opts.sourceFrame;
ImageReference.sourceType = char(opts.sourceType);
ImageReference.createdBy = char(opts.createdBy);
ImageReference.wasManuallyTransformed = opts.wasManuallyTransformed;
ImageReference.processingNotes = char(opts.processingNotes);

% Manual reference-image curation provenance.
% This transform is only a record of how ImageReference.image was curated.
% It is not the dataset-alignment transform.
ImageReference.manualTransformAppliedOn = [];
ImageReference.manualTransformType = '';
ImageReference.manualTransformNote = '';
ImageReference.manualTransformTform = [];
ImageReference.lastModifiedBy = char(opts.createdBy);

ImageReference.relativePath = referenceFileRel;
ImageReference.fileName = referenceFileName;

ImageReference.custom = struct();

validateImageReference(ImageReference);

save(referenceFileAbs, 'ImageReference', '-mat');

end

function mustBeNumericOrLogical2D(x)
%MUSTBENUMERICORLOGICAL2D Validate 2D numeric/logical image input.

if ~(isnumeric(x) || islogical(x)) || ndims(x) ~= 2 || isempty(x)
    error('createImageReference:InvalidReferenceImage', ...
        'referenceImage must be a non-empty 2D numeric or logical array.');
end
end

function value = iRequireText(value, name)
%IREQUIRETEXT Require non-empty text scalar.

value = strtrim(string(value));

if value == ""
    error('createImageReference:MissingMetadata', ...
        '"%s" must be provided.', name);
end
end

function part = iSanitizePathPart(part)
%ISANITIZEPATHPART Make one metadata token safe for folder names.

part = strtrim(char(string(part)));

invalidChars = '<>:"/\|?*';
for iChar = 1:numel(invalidChars)
    part = strrep(part, invalidChars(iChar), '_');
end

part = regexprep(part, '\s+', '_');

if isempty(part)
    part = 'Unknown';
end
end