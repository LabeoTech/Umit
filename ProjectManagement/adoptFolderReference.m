function resourceUUID = adoptFolderReference(saveFolder, opts)
%ADOPTFOLDERREFERENCE Add an aligned folder's committed reference to its subject library.
%
%   resourceUUID = adoptFolderReference(saveFolder)
%   resourceUUID = adoptFolderReference(saveFolder, 'Name', name, ...
%       'Description', text, 'SetActive', true)
%
%   Project redesign Phase 7 (section 4 Q7): call it after the user
%   confirmed the offer of planReferenceAdoption (state 'canAdd'). The
%   committed reference image in DataParams.registration is wrapped in an
%   ImageReference v2 with the bound project and subject identity, and
%   imported into the subject library (UMITProjectStore.addImageReference).
%   Its content hash equals the folder's commitment.
%
%   Read-only towards the folder: DataParams.registration and every other
%   folder file stay byte-identical. Refused unless the plan state is
%   'canAdd' (Umitoolbox:ProjectManagement:cannotAdoptReference); a
%   reference already in the library is never added twice.
%
%   Name-Value options:
%       Name        - library display name. Default: the reference name.
%       Description - Default: the reference description.
%       SetActive   - also make it the subject's active reference.
%                     Default: false.
%
%   See also planReferenceAdoption, importImageReferenceFile.

arguments
    saveFolder (1,1) string
    opts.Name (1,1) string = ""
    opts.Description (1,1) string = ""
    opts.SetActive (1,1) logical = false
end

saveFolder = char(saveFolder);
plan = planReferenceAdoption(saveFolder);
if ~strcmp(plan.state, 'canAdd')
    error('Umitoolbox:ProjectManagement:cannotAdoptReference', ...
        'The folder reference cannot be added to the library: %s', plan.message);
end

[context, store] = UMITProjectStore.resolveProjectBinding(saveFolder);
reference = loadDataParams(saveFolder).registration.reference;
name = char(opts.Name);
if isempty(name)
    name = char(string(reference.name));
end
if isempty(name)
    name = sprintf('Reference of %s/%s', context.subjectID, context.sessionID);
end
description = char(opts.Description);
if isempty(description)
    description = char(string(reference.description));
end

ImageReference = genImageReferenceStruct(reference.image, ...
    'Name', name, ...
    'Description', description, ...
    'ProjectUUID', context.projectUUID, ...
    'ProjectName', context.projectName, ...
    'SubjectUUID', context.subjectUUID, ...
    'SubjectID', context.subjectID, ...
    'SessionUUID', context.sessionUUID, ...
    'SessionID', context.sessionID, ...
    'SourceFolder', saveFolder, ...
    'SourceType', 'adoptedFromFolder', ...
    'AxisConvention', reference.axisConvention, ...
    'CreatedBy', 'adoptFolderReference');
if ~strcmp(ImageReference.contentHash, plan.committedHash)
    error('Umitoolbox:ProjectManagement:cannotAdoptReference', ...
        ['The folder''s reference image does not reproduce its committed hash ' ...
         '(%s); nothing was added.'], UMITLedgerCore.shortHash(plan.committedHash));
end

stagedPath = [tempname '.mat'];
cleanup = onCleanup(@() iDelete(stagedPath));
save(stagedPath, 'ImageReference', '-mat');
resourceUUID = store.addImageReference(context.subjectID, stagedPath, ...
    struct('displayName', name, 'description', description));
if opts.SetActive
    store.setActiveImageReference(context.subjectID, resourceUUID);
end
end

function iDelete(path)
if isfile(path)
    delete(path);
end
end
