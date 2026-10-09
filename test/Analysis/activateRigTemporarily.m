function fixture = activateRigTemporarily(rigUUID)
%ACTIVATERIGTEMPORARILY Make an existing Rig the sole Active/Default Rig.
%
%   fixture = activateRigTemporarily(rigUUID)
%
%   Snapshots whichever Rig is Active (or the raw default-pointer file, if
%   none is) before activating rigUUID via UMITRigStore.setDefaultRig, so
%   deactivateRigTemporarily(fixture) can restore the prior ambient
%   UMITRigStore state exactly -- including the previously-Active Rig's own
%   status, not just the default-pointer file (DFR-20260819-010: restoring
%   only the pointer file leaves the original Rig permanently demoted to
%   "available", corrupting the real ambient store for later work).

schema = getUMITRigSchema();
defaultFile = fullfile(UMITRigStore.getRigsRoot(), ...
    schema.store.internalFolder, schema.store.defaultFile);

fixture = struct();
fixture.PreviousActiveRigUUID = '';
fixture.HadDefaultFile = isfile(defaultFile);
fixture.DefaultFileBackup = '';
if fixture.HadDefaultFile
    fixture.DefaultFileBackup = [tempname '.mat'];
    copyfile(defaultFile, fixture.DefaultFileBackup);
end

lifecycle = UMITRigStore.validateLifecycle();
if lifecycle.activeCount == 1
    fixture.PreviousActiveRigUUID = lifecycle.activeRigUUID;
end

UMITRigStore.setDefaultRig(rigUUID);
end
