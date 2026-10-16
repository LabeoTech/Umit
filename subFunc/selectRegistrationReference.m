function reference = selectRegistrationReference(SaveFolder, source, opts)
%SELECTREGISTRATIONREFERENCE Record the folder's selected alignment reference.
%
%   reference = selectRegistrationReference(SaveFolder, source)
%   reference = selectRegistrationReference(SaveFolder, source, 'Interactive', false)
%
%   Resolves SOURCE (resolveImageReferenceSource; a source struct, or a
%   path: a .dat is a 'umitDat' source, anything else a 'file' source) and
%   stores it as DataParams.registration.reference (copy of the image,
%   content hash, descriptor, selectedOn), through validateDataParams +
%   atomic saveDataParams. The source folder is never written.
%
%   Before apply, selection is free: selecting another reference replaces
%   the previous one; an existing estimate is kept but no longer matches
%   (applying it is refused until the folder is estimated again, see
%   applyImageAlignmentToFolder). After apply the folder is committed:
%   selecting the committed reference again is a no-op, any other one is
%   refused (Umitoolbox:ImageRegistration:referenceCommitted).
%
%   See also resolveImageReferenceSource, createRegistrationTform,
%   assertFolderCanBeAligned.

arguments
    SaveFolder (1,1) string
    source
    opts.Interactive (1,1) logical = true
end

SaveFolder = char(SaveFolder);
reference = resolveImageReferenceSource(normalizeReferenceSourceInput(source), ...
    'Interactive', opts.Interactive);
DataParams = loadDataParams(SaveFolder);
registration = DataParams.registration;

if registration.apply.isRegistered
    if strcmp(reference.contentHash, registration.apply.appliedReferenceHash)
        return
    end
    why = sprintf(['the folder''s image data were aligned to %s and are committed to it; ' ...
        'the requested reference is %s.'], ...
        describeRegistrationReference(registration.reference, ...
        registration.apply.appliedReferenceHash), describeRegistrationReference(reference));
    message = UMITLedgerCore.buildErrorMessage('Selecting the alignment reference', ...
        SaveFolder, why, {['A fresh import (re-import of the raw data into the save ' ...
        'folder) is the only way to align this recording to a different reference.']});
    throwAsCaller(MException('Umitoolbox:ImageRegistration:referenceCommitted', '%s', message));
end

[~, referenceTemplate] = newRegistrationRecord();
record = referenceTemplate;
for fieldName = {'image', 'contentHash', 'name', 'description', 'axisConvention', ...
        'imageSizeYX', 'source'}
    record.(fieldName{1}) = reference.(fieldName{1});
end
record.selectedOn = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
DataParams.registration.reference = record;
validateDataParams(DataParams);
saveDataParams(SaveFolder, DataParams);
end
