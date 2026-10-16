function verifyRegistrationReference(folder, registration, operation)
%VERIFYREGISTRATIONREFERENCE Re-verify the stored reference image against its hash.
%
%   verifyRegistrationReference(folder, registration, operation)
%
%   The reference recorded in DataParams.registration is identified only by
%   the content hash of its stored image copy. Estimate and apply call this
%   before they rely on it; a stored image that no longer matches its hash
%   is refused. Read-only. A folder without a recorded reference passes.
%
%   Error:
%       Umitoolbox:ImageRegistration:referenceHashMismatch
%
%   See also computeImageContentHash, assertFolderCanBeAligned.

reference = registration.reference;
if isempty(reference)
    return
end
actual = computeImageContentHash(reference.image);
if strcmp(actual, reference.contentHash)
    return
end
why = sprintf(['the reference image stored in DataParams.mat (%s) no longer matches ' ...
    'its recorded content hash; it was modified outside UMIT.'], ...
    describeRegistrationReference(reference));
fixes = {'Select the reference again before estimating (it is free before apply).', ...
    ['If the folder is already aligned, restore DataParams.mat from a backup copy, ' ...
     'or re-import the raw data (fresh import).']};
message = UMITLedgerCore.buildErrorMessage(operation, folder, why, fixes);
throwAsCaller(MException('Umitoolbox:ImageRegistration:referenceHashMismatch', '%s', message));
end
