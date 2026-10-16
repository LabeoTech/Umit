function report = applyFreshImport(SaveFolder, opts)
%APPLYFRESHIMPORT Reset a SaveFolder before a re-import (fresh import).
%
%   report = applyFreshImport(SaveFolder, 'Backup', 'GENBACKUP', 'Origin', name)
%
%   Deletes every artefact planFreshImport lists for SaveFolder, so the
%   folder is equivalent to a first import (brief section 1.6), and logs the
%   fresh import in both ledgers. The caller has already obtained the user's
%   confirmation (GUI) or the explicit pipeline option; this function never
%   prompts.
%
%   Order:
%       1) GENBACKUP (or a name): zip every file to be deleted into the
%          SaveFolder (bkp_<yyyyMMdd_HHmmss>.zip or <name>.zip, made
%          unique); the zip is kept. 'ERASE' makes no zip.
%       2) Ledgers: a freshImport event in UMITRigLedger.json (all rig
%          commitments are reset) when it exists, and in
%          UMITProjectLedger.json when the folder is bound to a project.
%          Unbound folders get no project ledger.
%       3) Delete the files.
%   Raw acquisition files, ledgers, earlier backup zips, subfolders and
%   unknown files are never touched (planFreshImport).
%
%   Name-Value options:
%       Backup             - 'GENBACKUP' (default), 'ERASE', or a zip name.
%       Origin             - Calling tool, recorded in the ledger events.
%                            Required.
%       RecordLedgerEvents - false skips step 2 (calibration imports into
%                            temporary folders). Default: true.
%
%   Output REPORT: saveFolder, deleted (cellstr), backupPath, rigEvent,
%   projectEvent (logical), kept (planFreshImport lists).
%
%   See also planFreshImport, confirmFreshImport, UMITRigLedger,
%   UMITProjectLedger.

arguments
    SaveFolder (1,1) string
    opts.Backup (1,1) string = "GENBACKUP"
    opts.Origin (1,1) string
    opts.RecordLedgerEvents (1,1) logical = true
end

if ~isfield(opts, 'Origin') || strlength(opts.Origin) == 0
    error('Umitoolbox:DataImport:invalidInput', 'applyFreshImport needs a nonempty ''Origin''.');
end
SaveFolder = char(SaveFolder);
plan = planFreshImport(SaveFolder);
report = struct('saveFolder', SaveFolder, 'deleted', {{}}, 'backupPath', '', ...
    'rigEvent', false, 'projectEvent', false, 'kept', plan);

backup = char(opts.Backup);
if ~strcmpi(backup, 'ERASE') && ~isempty(plan.toDelete)
    if strcmpi(backup, 'GENBACKUP') || isempty(backup)
        backup = ['bkp_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
    end
    backup = regexprep(backup, '\.zip$', '', 'ignorecase');
    zipPath = fullfile(SaveFolder, [backup '.zip']);
    suffix = 0;
    while isfile(zipPath)
        suffix = suffix + 1;
        zipPath = fullfile(SaveFolder, sprintf('%s_%03d.zip', backup, suffix));
    end
    try
        zip(zipPath, plan.toDelete, SaveFolder);
    catch ME
        error('Umitoolbox:DataImport:freshImportBackupFailed', ...
            'The backup of "%s" failed; nothing was deleted. %s', SaveFolder, ME.message);
    end
    report.backupPath = zipPath;
end

origin = char(opts.Origin);
if opts.RecordLedgerEvents
    if plan.rigLedgerExists
        UMITRigLedger.recordFreshImport(SaveFolder, 'Origin', origin);
        report.rigEvent = true;
    end
    report.projectEvent = UMITProjectLedger.recordFreshImport(SaveFolder, 'Origin', origin);
end

for iFile = 1:numel(plan.toDelete)
    path = fullfile(SaveFolder, plan.toDelete{iFile});
    if isfile(path)
        delete(path);
    end
    report.deleted{end+1} = plan.toDelete{iFile};
end
end
