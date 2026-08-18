function displayName = resolveRigDisplayName(rigUUID, rigID)
%RESOLVERIGDISPLAYNAME Resolve the current human-readable name for one Rig.
%
%   displayName = resolveRigDisplayName(rigUUID, rigID)
%
%   A session/dataset stores a Rig's stable identity (rigUUID/rigID), while
%   UMITRigStore owns its editable, current human-readable display name.
%   This looks the display name up live so callers never show a stale
%   cached label; it falls back to rigID when the independent Rig store is
%   unavailable (moved rig root, transient error, etc.).

displayName = char(string(rigID));
if isempty(strtrim(displayName))
    displayName = '<none>';
    return
end

try
    if ~isempty(strtrim(char(string(rigUUID))))
        rigStore = UMITRigStore.open(rigUUID);
    else
        rigStore = UMITRigStore.openByRigID(rigID);
    end

    rigInfo = rigStore.getRigInfo();
    if isfield(rigInfo, 'displayName')
        candidate = char(string(rigInfo.displayName));
        if ~isempty(strtrim(candidate))
            displayName = candidate;
        end
    end
catch
    % Keep the stable Rig ID when the Rig store cannot be read.
end
end
