function [ImageReference, referenceFileAbs] = updateImageReferenceImage(referenceFile, newReferenceImage, opts)
%UPDATEIMAGEREFERENCEIMAGE Update the stored image inside an ImageReference file.
%
%   ImageReference = updateImageReferenceImage(referenceFile, newReferenceImage)
%
%   Loads an existing ImageReference_*.mat file, replaces ImageReference.image
%   with newReferenceImage, updates provenance fields, validates the structure,
%   and saves it back to disk.
%
%   Inputs:
%       referenceFile      - Absolute path or UMIT-root-relative path to an
%                            ImageReference_*.mat file.
%       newReferenceImage  - New 2D numeric/logical reference image.
%
%   Name-Value options:
%       manualTransformTform - Geometric transform used to curate the
%                              reference image. Stored only as provenance.
%       manualTransformNote  - Human-readable note about the curation.
%       updatedBy            - App/function name applying the update.
%
%   Outputs:
%       ImageReference   - Updated ImageReference structure.
%       referenceFileAbs - Absolute path to saved file.
%
%   Notes:
%       The manualTransformTform stored here must not be treated as the
%       dataset-alignment transform. Dataset alignment uses
%       DataParams.registration.tform after data are aligned to this
%       ImageReference.

arguments
    referenceFile (1,1) string
    newReferenceImage {mustBeNumericOrLogical2D}
    opts.manualTransformTform = []
    opts.manualTransformNote (1,1) string = ""
    opts.updatedBy (1,1) string = "ImageReferenceManager"
end

[ImageReference, referenceFileAbs] = loadImageReference(referenceFile);

if ~isa(newReferenceImage, 'single')
    newReferenceImage = single(newReferenceImage);
end

newReferenceImage(~isfinite(newReferenceImage)) = 0;

ImageReference.image = newReferenceImage;
ImageReference.imageSizeYX = size(newReferenceImage);
ImageReference.lastModified = datetime('now');

ImageReference.wasManuallyTransformed = true;
ImageReference.manualTransformAppliedOn = datetime('now');
ImageReference.manualTransformType = 'referenceImageCuration';
ImageReference.manualTransformNote = char(opts.manualTransformNote);
ImageReference.manualTransformTform = opts.manualTransformTform;
ImageReference.lastModifiedBy = char(opts.updatedBy);

validateImageReference(ImageReference);

save(referenceFileAbs, 'ImageReference', '-mat');

end

function mustBeNumericOrLogical2D(x)
%MUSTBENUMERICORLOGICAL2D Validate 2D numeric/logical image input.

if ~(isnumeric(x) || islogical(x)) || ndims(x) ~= 2 || isempty(x)
    error('updateImageReferenceImage:InvalidReferenceImage', ...
        'newReferenceImage must be a non-empty 2D numeric or logical array.');
end
end