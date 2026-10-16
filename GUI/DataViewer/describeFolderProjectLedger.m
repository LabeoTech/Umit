function summary = describeFolderProjectLedger(folder)
%DESCRIBEFOLDERPROJECTLEDGER DataViewer summary of a folder's project ledger.
%
%   summary = describeFolderProjectLedger(folder)
%
%   Read-only (never writes the ledger or the store). Never throws for an
%   invalid ledger: it is reported as state 'Invalid'. This describes the
%   ledger alone; whether the bound project exists on this machine
%   (orphan, mismatch, conflict) is decided by
%   UMITProjectStore.resolveProjectBinding. SUMMARY fields:
%       folder, exists, isValid, issues  - ledger file and integrity
%       state          'Unknown' (no ledger), 'Bound', 'Unbound' (binding
%                      removed, history kept) or 'Invalid'
%       binding        the bound identity (UUIDs, folderRole, boundOn,
%                      modifiedOn); [] when not bound
%       recordedNames  project name, subject ID and session ID as recorded
%                      in the latest bind/rebind snapshot, for display only
%                      (frozen at that time; never validated or refreshed)
%       displayText    short text for the folder context label
%       tooltipLines   detailed lines (cellstr)
%       canViewHistory the ledger exists and is readable
%
%   See also formatProjectLedgerHistory, UMITProjectLedger,
%   describeFolderRigLedger.

folder = char(string(folder));
summary = struct('folder', folder, 'exists', false, 'isValid', true, ...
    'issues', {{}}, 'state', 'Unknown', 'binding', [], ...
    'recordedNames', struct('projectName', '', 'subjectID', '', 'sessionID', ''), ...
    'displayText', '', 'tooltipLines', {{}}, 'canViewHistory', false);

if isempty(folder) || ~isfolder(folder)
    summary.displayText = '<none>';
    summary.tooltipLines = {'Project: no accessible folder.'};
    return
end

report = UMITProjectLedger.validate(folder);
summary.exists = report.exists;
summary.isValid = report.isValid;
summary.issues = report.issues;
if ~report.isValid
    summary.state = 'Invalid';
    summary.displayText = 'invalid project ledger';
    summary.tooltipLines = [{sprintf('Project ledger: INVALID (%s)', UMITProjectLedger.FileName)}, ...
        cellfun(@(issue) ['  ' issue], report.issues, 'UniformOutput', false), ...
        {'Project operations on this folder are blocked until the file is restored from a backup.'}];
    return
end

info = report.info;
summary.state = info.state;
summary.binding = info.binding;
summary.canViewHistory = info.exists;
summary.recordedNames = iRecordedNames(info.snapshot);

switch info.state
    case 'Unknown'
        summary.displayText = 'none (not bound to a project)';
        lines = {'Project ledger: none. The folder was never bound to a project.'};
    case 'Unbound'
        summary.displayText = 'none (binding removed)';
        lines = {'Project ledger: the folder''s project binding was removed; its history is kept.'};
    case 'Bound'
        names = summary.recordedNames;
        summary.displayText = iBoundText(names);
        lines = {sprintf('Project ledger: bound (%s since %s).', ...
            info.binding.folderRole, info.binding.boundOn)};
        lines{end+1} = sprintf('Project %s, subject %s, session %s.', ...
            info.binding.projectUUID, info.binding.subjectUUID, info.binding.sessionUUID);
        if ~isempty(names.projectName)
            lines{end+1} = sprintf('Recorded at binding: project "%s", subject "%s", session "%s".', ...
                names.projectName, names.subjectID, names.sessionID);
        end
end
summary.tooltipLines = lines;
end

function names = iRecordedNames(snapshot)
% Display only: the snapshot may hold anything (it is never validated).
names = struct('projectName', '', 'subjectID', '', 'sessionID', '');
if ~isstruct(snapshot) || ~isscalar(snapshot)
    return
end
names.subjectID = iText(snapshot, 'subjectID');
names.sessionID = iText(snapshot, 'sessionID');
if isfield(snapshot, 'project') && isstruct(snapshot.project) && isscalar(snapshot.project)
    names.projectName = iText(snapshot.project, 'projectName');
end
end

function text = iBoundText(names)
if isempty(names.projectName)
    text = 'bound';
else
    text = sprintf('%s / %s / %s', names.projectName, names.subjectID, names.sessionID);
end
end

function text = iText(value, name)
text = '';
if isfield(value, name)
    try
        text = char(string(value.(name)));
    catch
        text = '';
    end
    if size(text, 1) > 1
        text = '';
    end
end
end
