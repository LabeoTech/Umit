function summary = describeFolderRigLedger(saveFolder)
%DESCRIBEFOLDERRIGLEDGER DataViewer summary of a folder's rig ledger.
%
%   summary = describeFolderRigLedger(saveFolder)
%
%   Read-only (never writes the ledger, AcqInfos.mat or the rig library;
%   rig-ledger redesign, brief 1.8). Never throws for an invalid ledger: it
%   is reported as state 'Invalid'. SUMMARY fields:
%       folder, exists, isValid, issues   - ledger file and integrity
%       state            - 'Unknown' (no ledger), 'NoRig', 'BoundUncommitted',
%                          'BoundCommitted', 'Conflicted' or 'Invalid'
%       rigUUID, rigID, rigDisplayName    - the bound rig ('' when unbound)
%       hadConflict      - permanent mark: a conflict occurred at some point
%       conflictParts    - open conflict parts (cellstr)
%       committedLines   - one line per committed part (cellstr)
%       displayText      - short text for the folder context label
%       tooltipLines     - detailed lines for a tooltip (cellstr)
%       canBind          - Bind/Rebind is possible (ledger readable)
%       canAcknowledge   - there is an open conflict
%       canViewHistory   - the ledger exists and is readable
%       canRepair        - the ledger is invalid and can be quarantined
%                          (UMITRigLedger.quarantineInvalidLedger)

saveFolder = char(string(saveFolder));
summary = struct('folder', saveFolder, 'exists', false, 'isValid', true, ...
    'issues', {{}}, 'state', 'Unknown', 'rigUUID', '', 'rigID', '', ...
    'rigDisplayName', '', 'hadConflict', false, 'conflictParts', {{}}, ...
    'committedLines', {{}}, 'displayText', '', 'tooltipLines', {{}}, ...
    'canBind', false, 'canAcknowledge', false, 'canViewHistory', false, ...
    'canRepair', false);

if isempty(saveFolder) || ~isfolder(saveFolder)
    summary.displayText = '<none>';
    summary.tooltipLines = {'Rig: no accessible save folder.'};
    return
end

report = UMITRigLedger.validate(saveFolder);
summary.exists = report.exists;
summary.isValid = report.isValid;
summary.issues = report.issues;
if ~report.isValid
    summary.state = 'Invalid';
    summary.canRepair = true;
    summary.displayText = 'invalid rig ledger';
    summary.tooltipLines = [{sprintf('Rig ledger: INVALID (%s)', UMITRigLedger.FileName)}, ...
        cellfun(@(issue) ['  ' issue], report.issues, 'UniformOutput', false), ...
        {['Rig-dependent processing is blocked until the file is restored from a backup ' ...
        'or repaired (Management > Rig > Folder Rig > Repair Rig Ledger...).']}];
    return
end

info = report.info;
summary.state = info.state;
summary.hadConflict = info.hadConflict;
summary.conflictParts = info.conflictParts;
summary.canBind = true;
summary.canAcknowledge = strcmp(info.state, 'Conflicted');
summary.canViewHistory = info.exists;
if ~isempty(info.binding)
    sourceRig = info.binding.sourceRig;
    summary.rigUUID = sourceRig.uuid;
    summary.rigID = sourceRig.rigID;
    summary.rigDisplayName = sourceRig.displayName;
    if isempty(summary.rigDisplayName)
        summary.rigDisplayName = sourceRig.rigID;
    end
end

if ~isempty(info.commitments.optics)
    summary.committedLines{end+1} = 'optics';
end
if ~isempty(info.commitments.coregistration)
    summary.committedLines{end+1} = 'camera coregistration';
end
for iSet = 1:numel(info.commitments.filterSets)
    entry = info.commitments.filterSets(iSet);
    summary.committedLines{end+1} = sprintf('filter set "%s"', entry.id);
end

switch info.state
    case 'Unknown'
        summary.displayText = 'unknown (no rig ledger)';
        stateText = 'Unknown: the folder has no rig ledger.';
    case 'NoRig'
        summary.displayText = 'none (no rig)';
        stateText = 'No rig: imported without a rig binding.';
    case 'BoundUncommitted'
        summary.displayText = summary.rigDisplayName;
        stateText = 'Bound; nothing has used the binding yet (rebinding is free).';
    case 'BoundCommitted'
        summary.displayText = [summary.rigDisplayName ' (committed)'];
        stateText = 'Bound and committed.';
    case 'Conflicted'
        summary.displayText = [summary.rigDisplayName ' - CONFLICT'];
        stateText = sprintf(['CONFLICTED (%s): rig-dependent processing is blocked until ' ...
            'the binding is reverted or the conflict is acknowledged.'], ...
            strjoin(info.conflictParts, ', '));
end
if summary.hadConflict && ~strcmp(info.state, 'Conflicted')
    summary.displayText = [summary.displayText ' [had conflict]'];
end

summary.tooltipLines = {['Rig ledger: ' stateText]};
if ~isempty(info.binding)
    summary.tooltipLines{end+1} = sprintf('Bound rig: %s (%s)', ...
        summary.rigDisplayName, summary.rigID);
end
if ~isempty(summary.committedLines)
    summary.tooltipLines{end+1} = ['Committed: ' strjoin(summary.committedLines, '; ')];
end
if summary.hadConflict
    summary.tooltipLines{end+1} = ['This folder had a rig conflict (permanent mark): ' ...
        'see Management > Rig > Folder Rig > Rig Ledger History.'];
end
end
