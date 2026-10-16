function choice = confirmFreshImport(plan, opts)
%CONFIRMFRESHIMPORT Ask the user to confirm a fresh import (lists the deletions).
%
%   choice = confirmFreshImport(plan)
%   choice = confirmFreshImport(plan, 'AskBackup', false)
%
%   PLAN comes from planFreshImport. The dialog lists every file that will
%   be deleted and every unknown file that will be kept.
%       AskBackup true (default): the dialog also decides the backup and
%           returns 'GENBACKUP' (zip, then delete), 'ERASE' (delete without
%           a zip) or '' (cancelled: nothing may be touched). Used by the
%           importers when they are called directly with backupOpts = '' on
%           a desktop.
%       AskBackup false: the dialog only informs and asks to continue; it
%           returns 'CONTINUE' or ''. Used by DataViewer, where the import
%           step's backupOpts parameter is the only backup decision.
%
%   See also planFreshImport, applyFreshImport, prepareFreshImport.

arguments
    plan (1,1) struct
    opts.AskBackup (1,1) logical = true
end

maxListed = 25;
lines = {sprintf(['"%s" already holds an import. Re-importing is a FRESH IMPORT: ' ...
    'everything derived from the previous data is deleted (alignment, ROIs, ' ...
    'events, analysis outputs, history) and the rig commitments are reset.'], ...
    plan.saveFolder), '', sprintf('%d file(s) will be deleted:', numel(plan.toDelete))};
lines = [lines, iList(plan.toDelete, maxListed)];
if ~isempty(plan.unknown)
    lines = [lines, {'', sprintf('%d unknown item(s) are kept:', numel(plan.unknown))}, ...
        iList(plan.unknown, maxListed)];
end
lines = [lines, {'', 'Raw acquisition files and the rig/project ledgers are kept.'}];

if ~opts.AskBackup
    lines = [lines, {'', ['Whether the deleted files are zipped first is set by ' ...
        'the import parameter backupOpts (GENBACKUP or ERASE).']}];
    answer = questdlg(strjoin(lines, newline), 'Fresh import', 'Continue', 'Cancel', 'Cancel');
    choice = '';
    if strcmp(answer, 'Continue')
        choice = 'CONTINUE';
    end
    return
end

answer = questdlg(strjoin(lines, newline), 'Fresh import', ...
    'Zip, delete and import', 'Delete and import', 'Cancel', 'Cancel');
switch answer
    case 'Zip, delete and import'
        choice = 'GENBACKUP';
    case 'Delete and import'
        choice = 'ERASE';
    otherwise
        choice = '';
end
end

function lines = iList(names, maxListed)
lines = cellfun(@(n) ['  ' n], names(1:min(end, maxListed)), 'UniformOutput', false);
if numel(names) > maxListed
    lines{end+1} = sprintf('  ... and %d more', numel(names) - maxListed);
end
end
