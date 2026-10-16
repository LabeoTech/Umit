function info = describeProjectBinding(dataFolder)
%DESCRIBEPROJECTBINDING Classify a folder's project binding. Never writes.
%
%   info = describeProjectBinding(dataFolder)
%
%   Combines the folder's project ledger (UMITProjectLedger) with the
%   project store on this machine. The state is decided from the ledger,
%   the store's metadata and error identifiers, never from message text
%   (audit item 14). INFO fields:
%       folder      the folder, as given
%       state       one of:
%           'unbound'           no ledger, or its binding was removed
%           'invalid'           the ledger fails integrity validation
%           'orphan'            the bound project is not on this machine
%           'conflict'          the project UUID is registered more than
%                               once below the projects root
%           'projectUnreadable' the project exists but cannot be opened
%           'stale'             the project is here, but none of its
%                               sessions binds this folder any more (the
%                               session was removed, moved to another
%                               folder or rebound while the folder was
%                               offline; audit item 9)
%           'pathMismatch'      the session binds this binding, but its
%                               registered path differs (moved or copied
%                               folder)
%           'bound'             consistent; CONTEXT is filled
%       binding     the ledger binding ([] when unbound or invalid)
%       context     UMITProjectStore.resolveProjectBinding context ('bound')
%       identifier  error identifier behind the state ('' if none)
%       message     human explanation ('' when bound or unbound)
%       actions     the repair actions that fit the state (cellstr):
%                   'restoreLedger', 'removeOrphanBinding',
%                   'removeStaleBinding', 'relocate', 'repair',
%                   'bindToSession'
%
%   See also UMITProjectStore.resolveProjectBinding,
%   UMITProjectStore.removeOrphanProjectBinding,
%   UMITProjectStore.removeStaleProjectBinding, UMITProjectLedger.

folder = char(string(dataFolder));
info = struct('folder', folder, 'state', 'unbound', 'binding', [], ...
    'context', [], 'identifier', '', 'message', '', 'actions', {{}});

switch UMITProjectStore.getProjectBindingState(folder)
    case 'unbound'
        info.actions = {'bindToSession'};
        return
    case 'invalid'
        report = UMITProjectLedger.validate(folder);
        info.state = 'invalid';
        info.identifier = 'Umitoolbox:ProjectLedger:invalidLedger';
        info.message = strjoin(report.issues, ' | ');
        info.actions = {'restoreLedger'};
        return
end

binding = UMITProjectStore.readProjectBinding(folder);
info.binding = binding;

projects = UMITProjectStore.listProjects();
matches = sum(strcmpi(cellstr(projects.ProjectUUID), binding.projectUUID));
if matches == 0
    info.state = 'orphan';
    info.message = sprintf('Project %s is not on this machine.', binding.projectUUID);
    info.actions = {'removeOrphanBinding', 'bindToSession'};
    return
elseif matches > 1
    info.state = 'conflict';
    info.message = sprintf('Project %s is registered more than once on this machine.', ...
        binding.projectUUID);
    return
end

try
    store = UMITProjectStore.open(binding.projectUUID);
catch ME
    info.state = 'projectUnreadable';
    info.identifier = ME.identifier;
    info.message = ME.message;
    return
end

[pathField, bindingField] = iRoleFields(binding.folderRole);
try
    [SessionInfo, SubjectInfo] = store.getSessionInfoByUUID(binding.sessionUUID);
catch ME
    if ~any(strcmp(ME.identifier, {'Umitoolbox:UMITProjectStore:sessionNotFound', ...
            'Umitoolbox:UMITProjectStore:subjectNotFound'}))
        rethrow(ME)
    end
    info = iStale(info, ME.identifier, sprintf( ...
        'Session %s is no longer in project %s.', binding.sessionUUID, binding.projectUUID));
    return
end
if ~strcmpi(SubjectInfo.uuid, binding.subjectUUID)
    info = iStale(info, '', sprintf('Session %s now belongs to another subject.', ...
        binding.sessionUUID));
    return
end
if ~strcmpi(SessionInfo.(bindingField), binding.bindingUUID)
    info = iStale(info, '', sprintf(['Session %s is bound to another folder (or to ' ...
        'none) as its %s.'], SessionInfo.sessionID, pathField));
    return
end
if ~strcmp(iComparable(SessionInfo.(pathField)), iComparable(folder))
    info.state = 'pathMismatch';
    info.identifier = 'Umitoolbox:UMITProjectStore:bindingPathMismatch';
    info.message = sprintf(['Session %s registers this binding at "%s"; the folder ' ...
        'was moved or copied.'], SessionInfo.sessionID, SessionInfo.(pathField));
    info.actions = {'relocate', 'repair'};
    return
end

info.state = 'bound';
info.context = UMITProjectStore.resolveProjectBinding(folder);
end

function info = iStale(info, identifier, message)
info.state = 'stale';
info.identifier = identifier;
info.message = message;
info.actions = {'removeStaleBinding', 'bindToSession'};
end

function [pathField, bindingField] = iRoleFields(folderRole)
if strcmp(folderRole, 'rawDataFolder')
    pathField = 'rawDataFolder';
    bindingField = 'rawDataBindingUUID';
else
    pathField = 'processedDataFolder';
    bindingField = 'processedDataBindingUUID';
end
end

function text = iComparable(folder)
% Same comparison rule as the store (iNormalizeComparisonPath): canonical
% absolute path, case-insensitive on Windows and macOS.
text = char(string(folder));
if isempty(text)
    return
end
text = char(java.io.File(text).getCanonicalPath());
if ispc || ismac
    text = lower(text);
end
end
