function result = promptRepairRigLedger(parentFigure, saveFolder, opts)
%PROMPTREPAIRRIGLEDGER DataViewer Repair Rig Ledger workflow.
%
%   result = promptRepairRigLedger(parentFigure, saveFolder)
%   result = promptRepairRigLedger(..., 'Interactive', false, 'Choice', choice)
%
%   For a folder whose rig ledger is invalid, shows why, what repairing does
%   and what is lost, and asks for confirmation. On "Repair", moves the
%   invalid file aside with UMITRigLedger.quarantineInvalidLedger (Origin
%   'DataViewer'): nothing is edited or deleted, and the folder is left
%   without a rig so that Bind / Rebind Rig can bind one again. RESULT has
%   fields changed, message and severity (see runFolderRigBindingDialog).
%
%   Name-Value options:
%       Interactive - false: no dialogs (tests). Default true.
%       Choice      - with Interactive false: 'repair' or 'cancel'
%                     (default 'cancel').

arguments
    parentFigure
    saveFolder (1,:) char
    opts.Interactive (1,1) logical = true
    opts.Choice (1,1) string {mustBeMember(opts.Choice, ["repair", "cancel"])} = "cancel"
end

result = struct('changed', false, 'message', '', 'severity', 'info');
try
    summary = describeFolderRigLedger(saveFolder);
    if ~summary.canRepair
        result.message = sprintf('The rig ledger of this folder is not invalid (rig: %s); there is nothing to repair.', ...
            summary.displayText);
        if opts.Interactive
            uialert(parentFigure, result.message, 'Repair Rig Ledger', 'Icon', 'info');
        end
        return
    end

    issueLines = cellfun(@(issue) ['  - ' issue], summary.issues, 'UniformOutput', false);
    explanation = sprintf([ ...
        'The rig ledger of this folder (%s) is invalid, so rig-dependent processing ' ...
        '(for example HemoCompute) is blocked:\n\n%s\n\n' ...
        'Repairing moves the invalid file aside, in the same folder, as ' ...
        'UMITRigLedger.invalid_<date>T<time>.json. Nothing is edited or deleted. ' ...
        'The folder then has no rig: bind one again with Bind / Rebind Rig.\n\n' ...
        'WHAT IS LOST: the folder''s rig commitments (optics, camera coregistration and ' ...
        'filter-set snapshots) and its conflict history. Outputs computed before are ' ...
        'not recomputed, and a rig change made while the ledger was invalid cannot be ' ...
        'detected for them.\n\n' ...
        'Other options: restore UMITRigLedger.json from a backup copy of the folder, ' ...
        'or re-import the raw data into a fresh save folder.'], ...
        UMITRigLedger.FileName, strjoin(issueLines, newline));

    accepted = false;
    if opts.Interactive
        dialog = uifigure('Name', 'Repair Rig Ledger', 'WindowStyle', 'modal', ...
            'Position', [100 100 600 430]);
        closeDialog = onCleanup(@() iDeleteIfValid(dialog));
        movegui(dialog, 'center');
        try
            initializeUmitTheme(dialog);
        catch
            % Theme is cosmetic.
        end
        grid = uigridlayout(dialog, [2 2], 'RowHeight', {'1x', 'fit'}, ...
            'ColumnWidth', {'1x', '1x'});
        textArea = uitextarea(grid, 'Value', splitlines(explanation), 'Editable', 'off');
        textArea.Layout.Column = [1 2];
        uibutton(grid, 'Text', 'Repair', 'ButtonPushedFcn', @(~, ~) iClose(true));
        uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) iClose(false));
        dialog.CloseRequestFcn = @(~, ~) iClose(false);
        uiwait(dialog);
        % Delete now: the nested-function callbacks keep this workspace (and so
        % closeDialog) alive after return, and a lingering modal dialog would
        % hide the alerts that follow.
        iDeleteIfValid(dialog);
    else
        accepted = opts.Choice == "repair";
    end

    if ~accepted
        result.message = 'Rig ledger repair cancelled; the folder was left unchanged.';
        return
    end
    quarantine = UMITRigLedger.quarantineInvalidLedger(saveFolder, 'Origin', 'DataViewer');
    [~, name, ext] = fileparts(quarantine.quarantinedFile);
    result.changed = true;
    result.message = sprintf(['The invalid rig ledger was moved aside as %s. The folder has no ' ...
        'rig: bind one with Management > Rig > Folder Rig > Bind / Rebind Rig.'], [name ext]);
    result.severity = 'warning';
catch ME
    result.message = sprintf('Could not repair the rig ledger: %s', ME.message);
    result.severity = 'error';
    isLockError = false;
    lockMessage = '';
    if opts.Interactive
        [~, lockMessage, isLockError] = offerStaleLockClearing(parentFigure, ME, 'Folder', saveFolder);
    end
    if isLockError
        result.message = lockMessage;
        result.severity = 'warning';
    else
        if opts.Interactive
            uialert(parentFigure, ME.message, 'Repair Rig Ledger failed', 'Icon', 'error');
        end
    end
end

    function iClose(isAccepted)
        accepted = isAccepted;
        uiresume(dialog);
    end
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
