function resourceUUID = importImageReferenceFile(store, subjectID, filePath, opts)
%IMPORTIMAGEREFERENCEFILE Import a standalone reference file into a subject library.
%
%   resourceUUID = importImageReferenceFile(store, subjectID, filePath)
%   resourceUUID = importImageReferenceFile(..., 'ClearIdentity', true, ...
%       'ResourceInfo', struct('displayName', name, 'description', text))
%
%   Reads the ImageReference in FILEPATH (an exported reference file; a
%   version-1 payload is upgraded), optionally clears its identity, and
%   imports it with UMITProjectStore.addImageReference. The store accepts
%   only version 2 and verifies that any identity fields present name this
%   project and subject (Umitoolbox:UMITProjectStore:
%   imageReferenceIdentityMismatch). FILEPATH itself is never modified.
%
%   Name-Value options:
%       ClearIdentity - Empty the identity fields first (adopt a reference
%                       made elsewhere). Default: false.
%       ResourceInfo  - displayName/description of the library resource.
%                       Default: the reference name.
%
%   See also exportImageReferenceFile, UMITProjectStore.addImageReference.

arguments
    store (1,1) UMITProjectStore
    subjectID (1,1) string
    filePath (1,1) string
    opts.ClearIdentity (1,1) logical = false
    opts.ResourceInfo = []
end

reference = resolveImageReferenceSource(struct('kind', 'file', 'path', filePath), ...
    'Interactive', false);
if ~strcmp(reference.source.fileFormat, 'umitReference')
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        ['"%s" is not an exported image reference file. Create a reference from it ' ...
         'first (Image Reference Manager).'], filePath);
end

loaded = load(char(filePath), 'ImageReference', '-mat');
ImageReference = normalizeImageReferenceStruct(loaded.ImageReference);
if opts.ClearIdentity
    ImageReference = clearImageReferenceIdentity(ImageReference);
end

resourceInfo = opts.ResourceInfo;
if isempty(resourceInfo)
    resourceInfo = struct('displayName', ImageReference.name, ...
        'description', ImageReference.description);
end

stagedPath = [tempname '.mat'];
cleanup = onCleanup(@() iDelete(stagedPath));
save(stagedPath, 'ImageReference', '-mat');
resourceUUID = store.addImageReference(char(subjectID), stagedPath, resourceInfo);
end

function iDelete(path)
if isfile(path)
    delete(path);
end
end
