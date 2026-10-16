function [wasCleared, message, isLockError] = offerStaleLockClearing(parentFigure, ME, opts)
%OFFERSTALELOCKCLEARING After a "locked" error, offer to clear a stale lock.
%
%   wasCleared = offerStaleLockClearing(parentFigure, ME, 'Folder', folder)
%   wasCleared = offerStaleLockClearing(parentFigure, ME, 'Store', store)
%   [wasCleared, message, isLockError] = offerStaleLockClearing(..., 'Interactive', false)
%
%   Every write to a project, or to a folder's project or rig ledger, takes
%   a short-lived lock that is removed when the write ends. A crash, a
%   killed MATLAB or a lost network drive can leave the lock behind; every
%   later write then fails with a "locked" error. Call this from the error
%   handler of a write: when ME (or one of its causes) is such an error, it
%   explains the situation, asks whether to clear the lock, and clears it
%   only if it is old enough to be abandoned (project: 60 minutes; folder
%   ledgers: 10 minutes). A recent lock is kept: it may belong to a write
%   that is still running. No data or metadata are changed.
%
%   Any other error: returns false at once (no dialog), so callers can
%   call it unconditionally before showing their own error alert.
%
%   Name-Value options:
%       Folder      - the data folder whose ledger may be locked. A ledger
%                     "locked" error names its folder; that folder is used
%                     when it exists, so callers that do not know which
%                     folder was written (the Project Manager) can omit it.
%       Store       - the open UMITProjectStore whose project may be locked
%       Interactive - false: no dialogs (tests); the lock is cleared if
%                     stale. Default true.
%
%   Outputs: WASCLEARED (a lock was removed: the caller may retry), the
%   MESSAGE shown to the user, and ISLOCKERROR (ME was a lock error and the
%   user has been informed: the caller should not show its own alert).
%
%   Recognized errors: Umitoolbox:UMITProjectStore:projectLocked,
%   Umitoolbox:ProjectLedger:ledgerLocked, Umitoolbox:RigLedger:ledgerLocked.
%
%   See also UMITProjectStore.clearStaleLock, UMITProjectLedger.clearStaleLock,
%   UMITRigLedger.clearStaleLock.

arguments
    parentFigure
    ME (1,1) MException
    opts.Folder = ''
    opts.Store = []
    opts.Interactive (1,1) logical = true
end

wasCleared = false;
message = '';
[kind, lockError] = iLockKind(ME);
isLockError = ~isempty(kind);
if ~isLockError
    return
end

folder = iLockedFolder(lockError, opts.Folder);
switch kind
    case 'project'
        target = 'the project';
        holder = iProjectHolder(opts.Store);
    case 'projectLedger'
        target = sprintf('the project ledger of "%s"', folder);
        holder = iLedgerHolder(fullfile(folder, UMITProjectLedger.LockFolderName));
    otherwise
        target = sprintf('the rig ledger of "%s"', folder);
        holder = iLedgerHolder(fullfile(folder, UMITRigLedger.LockFolderName));
end
Target = [upper(target(1)) target(2:end)];

if opts.Interactive
    question = sprintf(['%s is locked%s.\n\nA lock is normally removed when the ' ...
        'operation that took it ends. A lock left behind by a crash or an ' ...
        'interrupted network connection blocks every change until it is cleared.\n\n' ...
        'Clear it only if no other UMIT window or MATLAB session is working on it. ' ...
        'A recent lock is kept.'], Target, holder);
    choice = uiconfirm(parentFigure, question, 'Locked', ...
        'Options', {'Clear Stale Lock', 'Cancel'}, ...
        'DefaultOption', 2, 'CancelOption', 2, 'Icon', 'warning');
    if ~strcmp(choice, 'Clear Stale Lock')
        return
    end
end

icon = 'success';
try
    switch kind
        case 'project'
            if isempty(opts.Store)
                error('Umitoolbox:DataViewer:noProjectStore', ...
                    'The project is not open here; clear its lock from the Project Manager.');
            end
            hadLock = ~isempty(fieldnames(opts.Store.getLockInfo()));
            opts.Store.clearStaleLock();
        case 'projectLedger'
            hadLock = UMITProjectLedger.clearStaleLock(folder);
        otherwise
            hadLock = UMITRigLedger.clearStaleLock(folder);
    end
    wasCleared = hadLock;
    if hadLock
        message = sprintf('The lock on %s was cleared. Try the operation again.', target);
    else
        message = sprintf('%s is no longer locked. Try the operation again.', Target);
    end
catch clearError
    icon = 'warning';
    if endsWith(clearError.identifier, ':lockNotStale')
        message = sprintf(['The lock on %s is recent, so it was kept: another UMIT ' ...
            'window or MATLAB session may still be writing. Try again later.\n\n%s'], ...
            target, clearError.message);
    else
        message = sprintf('The lock on %s could not be cleared: %s', target, ...
            clearError.message);
    end
end

if opts.Interactive
    uialert(parentFigure, message, 'Locked', 'Icon', icon);
end
end

function [kind, lockError] = iLockKind(ME)
% The lock error in ME or, depth-first, in its causes.
kind = '';
lockError = ME;
switch ME.identifier
    case 'Umitoolbox:UMITProjectStore:projectLocked'
        kind = 'project';
    case 'Umitoolbox:ProjectLedger:ledgerLocked'
        kind = 'projectLedger';
    case 'Umitoolbox:RigLedger:ledgerLocked'
        kind = 'rigLedger';
end
for iCause = 1:numel(ME.cause)
    if ~isempty(kind)
        return
    end
    [kind, lockError] = iLockKind(ME.cause{iCause});
end
end

function folder = iLockedFolder(lockError, fallback)
% A ledger "locked" error names its folder ('The <ledger> of "<folder>" is
% locked'); it is the folder that is locked, whatever the caller passed.
folder = char(string(fallback));
token = regexp(lockError.message, 'of "([^"]+)" is locked', 'tokens', 'once');
if ~isempty(token) && isfolder(token{1})
    folder = token{1};
end
end

function text = iProjectHolder(store)
text = '';
try
    info = store.getLockInfo();
    text = sprintf(' by "%s" since %s (user %s on %s)', info.operation, ...
        char(info.createdOn), info.userName, info.hostName);
catch
end
end

function text = iLedgerHolder(lockFolder)
text = '';
if isfolder(lockFolder)
    text = sprintf(' (lock age %s)', char(UMITLedgerCore.lockAge(lockFolder), 'hh:mm:ss'));
end
end
