function [report, plan] = prepareFreshImport(SaveFolder, opts)
%PREPAREFRESHIMPORT Importer gate: decide the backup, then reset the SaveFolder.
%
%   [report, plan] = prepareFreshImport(SaveFolder, 'backupOpts', b, 'Origin', name)
%
%   Called by the importers (ImagesClassification, importFromTif) before
%   they write anything. Re-importing into a SaveFolder that already holds
%   an import (planFreshImport.hasPreviousImport) is a FRESH IMPORT (brief
%   section 1.6): every artefact derived from the previous data is deleted.
%   BACKUPOPTS is the single decision (Bruno 2026-10-16):
%       'GENBACKUP'  - zip the files to delete, then delete them (the
%                      pipeline default);
%       'ERASE'      - delete them without a zip;
%       other name   - zip under that name, then delete;
%       ''           - a direct call on a MATLAB desktop shows the
%                      confirmation listing what will be deleted
%                      (confirmFreshImport: zip / delete / cancel); Cancel
%                      throws Umitoolbox:DataImport:freshImportCancelled. In
%                      batch MATLAB '' means GENBACKUP.
%   Pipelines are not refused (a deviation from brief 1.6, decided
%   2026-10-16): the default GENBACKUP keeps the deleted files in a zip.
%   A SaveFolder without a previous import passes through unchanged
%   (REPORT is []), whatever BACKUPOPTS says.
%
%   Name-Value options:
%       backupOpts         - '' | 'GENBACKUP' | 'ERASE' | zip name.
%       Origin             - Importer name (ledger events).
%       RecordLedgerEvents - See applyFreshImport. Default true.
%
%   See also planFreshImport, applyFreshImport, confirmFreshImport.

arguments
    SaveFolder (1,1) string
    opts.backupOpts (1,1) string = ""
    opts.Origin (1,1) string
    opts.RecordLedgerEvents (1,1) logical = true
end

report = [];
plan = planFreshImport(SaveFolder);
if ~plan.hasPreviousImport
    return
end

backup = strtrim(char(opts.backupOpts));
if isempty(backup)
    if usejava('desktop') && ~batchStartupOptionUsed
        backup = confirmFreshImport(plan);
        if isempty(backup)
            error('Umitoolbox:DataImport:freshImportCancelled', ...
                'The fresh import into "%s" was cancelled; nothing was changed.', ...
                plan.saveFolder);
        end
    else
        backup = 'GENBACKUP';
    end
end

report = applyFreshImport(SaveFolder, 'Backup', backup, 'Origin', opts.Origin, ...
    'RecordLedgerEvents', opts.RecordLedgerEvents);
end
