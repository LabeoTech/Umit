function text = describeRegistrationReference(reference, hash)
%DESCRIBEREGISTRATIONREFERENCE User-facing text for a reference: name, kind, short hash.
%
%   text = describeRegistrationReference(reference)
%   text = describeRegistrationReference(reference, hash)
%
%   REFERENCE is DataParams.registration.reference (or [] when unknown) or
%   a resolveImageReferenceSource result. HASH overrides the hash shown
%   (e.g. the committed appliedReferenceHash). Example output:
%       reference "Mouse12 vessels" (project, sha256:2f1e8f998893)
%
%   See also assertFolderCanBeAligned, UMITLedgerCore.shortHash.

if nargin < 2
    hash = '';
end
if isempty(reference)
    if isempty(hash)
        text = 'an unknown reference';
    else
        text = sprintf('an unrecorded reference (%s)', UMITLedgerCore.shortHash(hash));
    end
    return
end
if isempty(hash)
    hash = reference.contentHash;
end
kind = '';
if isfield(reference, 'source') && isstruct(reference.source) && isfield(reference.source, 'kind')
    kind = [char(string(reference.source.kind)) ', '];
end
text = sprintf('reference "%s" (%s%s)', char(string(reference.name)), kind, ...
    UMITLedgerCore.shortHash(hash));
end
