function assertFolderCanBeAligned(folder, DataParams, operation, requested)
%ASSERTFOLDERCANBEALIGNED Refuse alignment work on an already-aligned folder.
%
%   assertFolderCanBeAligned(folder, DataParams, operation)
%   assertFolderCanBeAligned(folder, DataParams, operation, requested)
%
%   Shared guard of every image-alignment entry point (estimate and apply,
%   GUI and pipeline). Applying an alignment rewrites the folder's image
%   data irreversibly, so an aligned folder is committed to the reference it
%   was aligned to: estimating or applying again, even against the same
%   reference, is refused (brief section 4 Q1: re-applying would transform
%   the data twice). The only way to align the data again is a fresh import
%   of the raw data.
%
%   Inputs:
%       folder     - Save folder (used in the message).
%       DataParams - The folder's DataParams struct (loadDataParams;
%                    registration version 2).
%       operation  - Refused operation, e.g. 'Estimating the image
%                    alignment' (used in the message).
%       requested  - Optional: the reference the caller asked for (a
%                    resolveImageReferenceSource result or a registration
%                    reference record), named in the message.
%
%   Error:
%       Umitoolbox:ImageRegistration:alreadyAligned when
%       DataParams.registration.apply.isRegistered is true. The message
%       names the operation, the folder, the committed reference (name,
%       source kind, short hash), the requested one, and the fix.
%
%   See also applyImageAlignmentToFolder, createRegistrationTform,
%   selectRegistrationReference, describeRegistrationReference.

if nargin < 4
    requested = [];
end
registration = DataParams.registration;
if ~registration.apply.isRegistered
    return
end

apply = registration.apply;
when = '';
if ~isempty(apply.appliedOn)
    when = sprintf(' on %s', char(string(apply.appliedOn)));
end
if ~isempty(apply.appliedBy)
    when = sprintf('%s (by %s)', when, char(string(apply.appliedBy)));
end
committed = describeRegistrationReference(registration.reference, apply.appliedReferenceHash);
why = sprintf(['the folder''s image data were already aligned to %s%s. Alignment ' ...
    'rewrites the data irreversibly, so the folder is committed to that reference; ' ...
    'aligning again, even to the same reference, would transform the data twice.'], ...
    committed, when);
if ~isempty(requested)
    if strcmp(requested.contentHash, apply.appliedReferenceHash)
        why = sprintf('%s The requested reference is that same %s.', why, ...
            describeRegistrationReference(requested));
    else
        why = sprintf('%s The requested reference is %s.', why, ...
            describeRegistrationReference(requested));
    end
end
fixes = {['To align this recording to another reference (or again), re-import ' ...
    'the raw data into the save folder (a fresh import), then align the new import.']};
message = UMITLedgerCore.buildErrorMessage(operation, folder, why, fixes);
throwAsCaller(MException('Umitoolbox:ImageRegistration:alreadyAligned', '%s', message));
end
