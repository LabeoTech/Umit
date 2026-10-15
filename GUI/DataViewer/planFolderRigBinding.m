function plan = planFolderRigBinding(saveFolder, rigUUID)
%PLANFOLDERRIGBINDING Preview binding a folder to a library rig (read-only).
%
%   plan = planFolderRigBinding(saveFolder, rigUUID)
%
%   Builds the binding DataViewer's Bind/Rebind action would write and
%   reports, BEFORE anything is written, which committed parts it would
%   change. Never writes. The binding takes the OPTICS of the chosen library
%   rig and keeps the folder's current camera-coregistration binding: the
%   transform applied to the data can only change through a re-import
%   (rig-ledger redesign, brief Phase 4).
%
%   PLAN fields:
%       action        - 'bind' (folder Unknown/NoRig), 'rebind', or 'none'
%                       (already bound to this rig with these optics)
%       rigUUID, rigLabel - the chosen rig
%       binding       - binding for UMITRigLedger.bind/rebind
%       conflict      - UMITRigLedger.previewRebind result (rebind only)
%       confirmTitle, confirmText, confirmButton - text for the confirmation
%
%   Errors: an invalid ledger (Umitoolbox:RigLedger:invalidLedger), an
%   unknown rig (Umitoolbox:UMITRigStore:openFailed) or an archived rig
%   (Umitoolbox:DataViewer:archivedRig).

saveFolder = char(string(saveFolder));
info = UMITRigLedger.read(saveFolder);
store = UMITRigStore.open(char(string(rigUUID)));
RigInfo = store.getRigInfo();
if strcmpi(RigInfo.status, 'archived')
    error('Umitoolbox:DataViewer:archivedRig', ...
        'Rig "%s" (%s) is archived and cannot be bound to a folder. Restore it in the Rig Manager first.', ...
        RigInfo.displayName, RigInfo.rigID);
end

coregistration = [];
if ~isempty(info.binding) && ~isempty(info.binding.coregistrationHash)
    coregistration = UMITRigLedger.findSnapshot(info, 'coregistration', ...
        info.binding.coregistrationHash);
end
plan = struct();
plan.rigUUID = RigInfo.uuid;
plan.rigLabel = iRigLabel(RigInfo.displayName, RigInfo.rigID);
plan.binding = UMITRigLedger.makeBinding(UMITRigLedger.sourceRigFromRigInfo(RigInfo), ...
    UMITRigLedger.buildOpticsSnapshot(RigInfo), coregistration);
plan.conflict = struct('hasConflict', false, 'parts', {{}}, 'message', '');

if isempty(info.binding)
    plan.action = 'bind';
    plan.confirmTitle = 'Bind Rig';
    plan.confirmButton = 'Bind';
    plan.confirmText = sprintf('Bind rig %s to this folder?\n\nFolder: %s', ...
        plan.rigLabel, saveFolder);
    return
end

current = info.binding.sourceRig;
currentLabel = iRigLabel(current.displayName, current.rigID);
if strcmpi(current.uuid, RigInfo.uuid) && ...
        strcmp(info.binding.opticsHash, plan.binding.optics.hash)
    plan.action = 'none';
    plan.confirmTitle = 'Rebind Rig';
    plan.confirmButton = '';
    plan.confirmText = sprintf('The folder is already bound to rig %s with these optics.', ...
        plan.rigLabel);
    return
end

plan.action = 'rebind';
plan.conflict = UMITRigLedger.previewRebind(saveFolder, plan.binding);
plan.confirmTitle = 'Rebind Rig';
plan.confirmButton = 'Rebind';
if plan.conflict.hasConflict
    plan.confirmText = sprintf(['Rebinding this folder from rig %s to rig %s CHANGES ' ...
        'COMMITTED PARTS:\n\n%s\n\nThe folder will become CONFLICTED: rig-dependent ' ...
        'processing (e.g. HemoCompute) is blocked until you revert the binding, use a ' ...
        'fresh save folder, or acknowledge the conflict. The conflict is recorded ' ...
        'permanently.\n\nRebind anyway?'], currentLabel, plan.rigLabel, plan.conflict.message);
else
    plan.confirmText = sprintf(['Rebind this folder from rig %s to rig %s?\n\n' ...
        'No committed part changes, so no conflict is created.'], currentLabel, plan.rigLabel);
end
end

function label = iRigLabel(displayName, rigID)
if isempty(displayName) || strcmp(displayName, rigID)
    label = sprintf('"%s"', rigID);
else
    label = sprintf('"%s" (%s)', displayName, rigID);
end
end
