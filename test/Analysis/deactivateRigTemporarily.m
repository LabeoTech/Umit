function deactivateRigTemporarily(fixture)
%DEACTIVATERIGTEMPORARILY Undo activateRigTemporarily.
%
%   deactivateRigTemporarily(fixture)
%
%   Restores whichever Rig was Active before activateRigTemporarily ran (if
%   any), via UMITRigStore.activateRig() so both the Rig's own status and
%   the store-level default pointer are restored together. If no Rig was
%   Active before, the raw default-pointer file is restored byte-for-byte
%   instead (including "no file existed"), leaving pre-existing ambient
%   state exactly as found rather than introducing new corruption.

if ~isstruct(fixture)
    return
end

if isfield(fixture, 'PreviousActiveRigUUID') && ...
        ~isempty(fixture.PreviousActiveRigUUID)
    try %#ok<TRYNC>
        UMITRigStore.open(fixture.PreviousActiveRigUUID).activateRig();
    end
else
    schema = getUMITRigSchema();
    defaultFile = fullfile(UMITRigStore.getRigsRoot(), ...
        schema.store.internalFolder, schema.store.defaultFile);
    if isfile(defaultFile)
        delete(defaultFile);
    end
    if isfield(fixture, 'HadDefaultFile') && fixture.HadDefaultFile && ...
            isfile(fixture.DefaultFileBackup)
        copyfile(fixture.DefaultFileBackup, defaultFile);
    end
end

if isfield(fixture, 'DefaultFileBackup') && ...
        ~isempty(fixture.DefaultFileBackup) && isfile(fixture.DefaultFileBackup)
    delete(fixture.DefaultFileBackup);
end
end
