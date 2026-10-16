function plan = planReferenceAdoption(saveFolder)
%PLANREFERENCEADOPTION What binding an aligned folder means for its subject library.
%
%   plan = planReferenceAdoption(saveFolder)
%
%   Project redesign Phase 7 (section 4 Q7). A folder aligned without a
%   project carries its committed reference in DataParams.registration.
%   Once the folder is bound to a project session, this read-only plan
%   compares that reference with the subject's image-reference library by
%   content hash:
%     - same hash in the library: the match is shown only (nothing is
%       recorded anywhere);
%     - no match: the reference can be added to the library, with the
%       user's confirmation (adoptFolderReference).
%   DataParams.registration is never modified: it records what was applied.
%
%   PLAN fields:
%       folder, state ('notBound', 'notAligned', 'inLibrary', 'canAdd')
%       committedHash, referenceName, sourceKind  - the folder's commitment
%       subjectID, projectName                    - the bound subject
%       match     - [] or the matching library resource (uuid, displayName,
%                   status), active or archived
%       message   - one line for the GUI
%
%   See also adoptFolderReference, describeFolderRegistration.

folder = char(string(saveFolder));
plan = struct('folder', folder, 'state', 'notBound', 'committedHash', '', ...
    'referenceName', '', 'sourceKind', '', 'subjectID', '', 'projectName', '', ...
    'match', [], 'message', '');

if ~strcmp(describeProjectBinding(folder).state, 'bound')
    plan.message = 'The folder is not bound to a project session.';
    return
end
[context, store] = UMITProjectStore.resolveProjectBinding(folder);
plan.subjectID = context.subjectID;
plan.projectName = context.projectName;

registration = loadDataParams(folder).registration;
if ~registration.apply.isRegistered || isempty(registration.reference)
    plan.state = 'notAligned';
    plan.message = 'The folder is not aligned: there is no committed reference to adopt.';
    return
end
plan.committedHash = registration.apply.appliedReferenceHash;
if isempty(plan.committedHash)
    plan.committedHash = registration.reference.contentHash;
end
plan.referenceName = char(string(registration.reference.name));
plan.sourceKind = registration.reference.source.kind;

resources = store.listSubjectResources(context.subjectID, ...
    'Type', 'imageReference', 'Status', {});
for iResource = 1:numel(resources)
    hash = iLibraryHash(resources(iResource).absolutePath);
    if strcmp(hash, plan.committedHash)
        plan.state = 'inLibrary';
        plan.match = struct('uuid', resources(iResource).uuid, ...
            'displayName', resources(iResource).displayName, ...
            'status', resources(iResource).status);
        plan.message = sprintf(['The committed reference "%s" (%s) is the library ' ...
            'reference "%s" (%s) of subject %s.'], plan.referenceName, ...
            UMITLedgerCore.shortHash(plan.committedHash), plan.match.displayName, ...
            plan.match.status, plan.subjectID);
        return
    end
end
plan.state = 'canAdd';
plan.message = sprintf(['The committed reference "%s" (%s) is not in the library of ' ...
    'subject %s. It can be added so other recordings can use it.'], ...
    plan.referenceName, UMITLedgerCore.shortHash(plan.committedHash), plan.subjectID);
end

function hash = iLibraryHash(path)
% Content hash of one library payload; '' when it cannot be read.
hash = '';
try
    loaded = load(path, 'ImageReference', '-mat');
    hash = normalizeImageReferenceStruct(loaded.ImageReference).contentHash;
catch
end
end
