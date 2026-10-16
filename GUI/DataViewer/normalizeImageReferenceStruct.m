function ImageReference = normalizeImageReferenceStruct(ImageReference)
%NORMALIZEIMAGEREFERENCESTRUCT Validate an ImageReference and return it as v2.
%
%   ImageReference = normalizeImageReferenceStruct(ImageReference)
%
%   Validates the struct (validateImageReferenceStruct) and upgrades a
%   version-1 struct IN MEMORY to version 2: the image is normalized as
%   genImageReferenceStruct does and contentHash is added. Nothing is
%   written. Version-1 files therefore stay readable (existing subject
%   libraries), while every new payload is written as version 2.
%
%   See also genImageReferenceStruct, validateImageReferenceStruct,
%   computeImageContentHash.

validateImageReferenceStruct(ImageReference);
if ImageReference.version == 2
    return
end

[contentHash, image] = computeImageContentHash(ImageReference.image);
upgraded = struct('version', 2, 'image', image, 'contentHash', contentHash);
names = fieldnames(ImageReference);
for iName = 1:numel(names)
    if ~ismember(names{iName}, {'version', 'image'})
        upgraded.(names{iName}) = ImageReference.(names{iName});
    end
end
ImageReference = upgraded;
validateImageReferenceStruct(ImageReference);
end
