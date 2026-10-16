function lines = formatRigLedgerHistory(saveFolder)
%FORMATRIGLEDGERHISTORY Human-readable event history of a folder rig ledger.
%
%   lines = formatRigLedgerHistory(saveFolder)
%
%   Read-only. Returns a cellstr: a header (state, bound rig, hadConflict)
%   followed by one line per ledger event, oldest first, with timestamps in
%   UTC as stored. Snapshot hashes are not shown: a rebind says whether the
%   optics and coregistration changed. Rigs are named as recorded in the ledger, so the history
%   reads the same on a machine whose rig library is empty. An invalid
%   ledger throws Umitoolbox:RigLedger:invalidLedger; a folder without a
%   ledger returns a single explanatory line.

info = UMITRigLedger.read(saveFolder);
if ~info.exists
    lines = {sprintf('The folder has no rig ledger (%s): its rig is unknown.', ...
        UMITRigLedger.FileName)};
    return
end

lines = {sprintf('Folder: %s', info.folder), sprintf('State: %s', info.state)};
if ~isempty(info.binding)
    lines{end+1} = sprintf('Bound rig: %s', iRig(info.binding.sourceRig));
end
if info.hadConflict
    lines{end+1} = 'This folder had a rig conflict (permanent mark).';
end
lines{end+1} = sprintf('Ledger %s, created %s', info.ledgerUUID, info.createdOn);
lines{end+1} = '';

for iEvent = 1:numel(info.events)
    event = info.events{iEvent};
    switch event.type
        case 'noRig'
            detail = 'recorded "no rig"';
        case 'bind'
            if isempty(event.coregistrationHash)
                detail = sprintf('bound %s (optics only, no camera coregistration)', ...
                    iRig(event.sourceRig));
            else
                detail = sprintf('bound %s (optics and camera coregistration)', ...
                    iRig(event.sourceRig));
            end
        case 'rebind'
            detail = sprintf('rebound to %s: optics %s, camera coregistration %s', ...
                iRig(event.sourceRig), ...
                iChange(event.previousOpticsHash, event.opticsHash), ...
                iChange(event.previousCoregistrationHash, event.coregistrationHash));
        case 'commit'
            part = event.part;
            if strcmp(part, 'filterSet')
                part = sprintf('filter set "%s"', event.filterSetID);
            elseif strcmp(part, 'coregistration')
                part = 'camera coregistration';
            end
            detail = sprintf('committed %s (first used by %s)', part, event.firstUsedBy);
        case 'freshImport'
            detail = 'fresh import of the raw data: all rig commitments reset';
        case 'conflict'
            detail = sprintf('CONFLICT on %s', strjoin(cellstr(event.parts(:)'), ', '));
        case 'acknowledgeConflict'
            detail = sprintf('conflict acknowledged (%s)', strjoin(cellstr(event.parts(:)'), ', '));
            if ~isempty(event.comment)
                detail = sprintf('%s: "%s"', detail, event.comment);
            end
        otherwise
            detail = event.type;
    end
    lines{end+1} = sprintf('%s  %-20s %s  [%s]', event.timestamp, event.type, ...
        detail, event.origin); %#ok<AGROW>
end
end

function text = iChange(previousHash, newHash)
% Snapshot hashes are compared, never shown.
if isempty(previousHash) && isempty(newHash)
    text = 'none';
elseif isempty(previousHash)
    text = 'added';
elseif isempty(newHash)
    text = 'removed';
elseif strcmp(previousHash, newHash)
    text = 'unchanged';
else
    text = 'changed';
end
end

function text = iRig(sourceRig)
if isempty(sourceRig.displayName) || strcmp(sourceRig.displayName, sourceRig.rigID)
    text = sprintf('rig "%s"', sourceRig.rigID);
else
    text = sprintf('rig "%s" (%s)', sourceRig.displayName, sourceRig.rigID);
end
end
