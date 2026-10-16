function ImageReference = clearImageReferenceIdentity(ImageReference)
%CLEARIMAGEREFERENCEIDENTITY Empty the identity fields of a version-2 ImageReference.
%
%   ImageReference = clearImageReferenceIdentity(ImageReference)
%
%   Clears projectUUID/projectName, subjectUUID/subjectID and
%   sessionUUID/sessionID. The image, its content hash and every other
%   field are kept, so the reference keeps its identity by content.
%
%   See also exportImageReferenceFile, importImageReferenceFile.

ImageReference = normalizeImageReferenceStruct(ImageReference);
for fieldName = {'projectUUID', 'projectName', 'subjectUUID', 'subjectID', ...
        'sessionUUID', 'sessionID'}
    ImageReference.(fieldName{1}) = '';
end
validateImageReferenceStruct(ImageReference);
end
