function lines = formatProjectLedgerHistory(folder)
%FORMATPROJECTLEDGERHISTORY Human-readable event history of a project ledger.
%
%   lines = formatProjectLedgerHistory(folder)
%
%   Read-only. Returns a cellstr: a header (state, bound identity) followed
%   by one line per ledger event, oldest first, with timestamps in UTC as
%   stored. Project, subject and session names come from each event's
%   snapshot, as recorded at that time (display only; they may be out of
%   date and are never validated). An invalid ledger throws
%   Umitoolbox:ProjectLedger:invalidLedger; a folder without a ledger
%   returns a single explanatory line.
%
%   See also describeFolderProjectLedger, formatRigLedgerHistory.

info = UMITProjectLedger.read(folder);
if ~info.exists
    lines = {sprintf('The folder has no project ledger (%s): it was never bound to a project.', ...
        UMITProjectLedger.FileName)};
    return
end

lines = {sprintf('Folder: %s', info.folder), sprintf('State: %s', info.state)};
if ~isempty(info.binding)
    lines{end+1} = sprintf('Bound as %s: project %s, subject %s, session %s', ...
        info.binding.folderRole, info.binding.projectUUID, ...
        info.binding.subjectUUID, info.binding.sessionUUID);
end
lines{end+1} = sprintf('Ledger %s, created %s', info.ledgerUUID, info.createdOn);
lines{end+1} = '';

for iEvent = 1:numel(info.events)
    event = info.events{iEvent};
    switch event.type
        case 'bind'
            detail = sprintf('bound to %s', iNames(event));
        case 'rebind'
            detail = sprintf('rebound to %s on %s', iNames(event), iMachine(event));
        case 'unbind'
            detail = sprintf('binding removed (%s)', iReason(event));
        case 'orphanRemoval'
            detail = sprintf(['binding removed on %s: its project was not found ' ...
                'on that machine'], iMachine(event));
        case 'freshImport'
            detail = 'fresh import of the raw data; the binding is unchanged';
        otherwise
            detail = event.type;
    end
    lines{end+1} = sprintf('%s  %-14s %s  [%s]', event.timestamp, event.type, ...
        detail, event.origin); %#ok<AGROW>
end
end

function text = iNames(event)
names = struct('project', '', 'subject', '', 'session', '');
if isfield(event, 'snapshot') && isstruct(event.snapshot) && isscalar(event.snapshot)
    snapshot = event.snapshot;
    names.subject = iText(snapshot, 'subjectID');
    names.session = iText(snapshot, 'sessionID');
    if isfield(snapshot, 'project') && isstruct(snapshot.project) && isscalar(snapshot.project)
        names.project = iText(snapshot.project, 'projectName');
    end
end
if isempty(names.project)
    text = sprintf('project %s', event.binding.projectUUID);
else
    text = sprintf('project "%s", subject "%s", session "%s"', ...
        names.project, names.subject, names.session);
end
text = sprintf('%s (%s)', text, event.binding.folderRole);
end

function text = iMachine(event)
text = 'an unknown machine';
if isfield(event, 'machine')
    if isstruct(event.machine) && isfield(event.machine, 'hostName')
        value = iText(event.machine, 'hostName');
    else
        value = iText(event, 'machine');
    end
    if ~isempty(value)
        text = value;
    end
end
end

function text = iReason(event)
switch iText(event, 'reason')
    case 'sessionRemoved'
        text = 'the session was removed from its project';
    case 'removedFromProject'
        text = 'the folder was removed from its project';
    case 'relocated'
        text = 'the session was moved to another folder';
    case 'folderUnbound'
        text = 'the folder was unbound';
    case 'staleBinding'
        text = 'stale link cleared: no session of the project bound this folder any more';
    case ''
        text = 'no reason recorded';
    otherwise
        text = iText(event, 'reason');
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
