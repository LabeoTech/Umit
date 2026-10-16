function result = promptAcknowledgeRigConflict(parentFigure, saveFolder)
%PROMPTACKNOWLEDGERIGCONFLICT DataViewer Acknowledge Conflict workflow.
%
%   result = promptAcknowledgeRigConflict(parentFigure, saveFolder)
%
%   Shows the open conflict of the folder, what acknowledging means, and an
%   optional comment field. On "Acknowledge", writes
%   UMITRigLedger.acknowledgeConflict with Origin 'DataViewer'. RESULT has
%   fields changed, message and severity (see runFolderRigBindingDialog).

result = struct('changed', false, 'message', '', 'severity', 'info');
saveFolder = char(string(saveFolder));
try
    summary = describeFolderRigLedger(saveFolder);
    if ~summary.canAcknowledge
        result.message = sprintf('The folder has no open rig conflict (rig: %s).', ...
            summary.displayText);
        uialert(parentFigure, result.message, 'Acknowledge Rig Conflict', 'Icon', 'info');
        return
    end

    explanation = sprintf([ ...
        'This folder''s rig binding conflicts with its commitments (%s).\n\n%s\n\n' ...
        'Acknowledging accepts the CURRENT binding (%s) for this folder. Outputs ' ...
        'computed before the rebind are not changed; step skipping will recompute ' ...
        'rig-dependent results. The folder stays permanently marked "had conflict".\n\n' ...
        'Other fixes: revert the binding (Bind / Rebind back to the committed rig), or ' ...
        'start a fresh save folder.'], strjoin(summary.conflictParts, ', '), ...
        strjoin(summary.tooltipLines, newline), summary.rigDisplayName);

    dialog = uifigure('Name', 'Acknowledge Rig Conflict', 'WindowStyle', 'modal', ...
        'Position', [100 100 560 420]);
    closeDialog = onCleanup(@() iDeleteIfValid(dialog));
    movegui(dialog, 'center');
    try
        initializeUmitTheme(dialog);
    catch
        % Theme is cosmetic.
    end
    grid = uigridlayout(dialog, [4 2], 'RowHeight', {'1x', 'fit', 70, 'fit'}, ...
        'ColumnWidth', {'1x', '1x'});
    textArea = uitextarea(grid, 'Value', splitlines(explanation), 'Editable', 'off');
    textArea.Layout.Column = [1 2];
    commentLabel = uilabel(grid, 'Text', 'Comment (optional, stored in the ledger):');
    commentLabel.Layout.Column = [1 2];
    comment = uitextarea(grid);
    comment.Layout.Column = [1 2];
    accepted = false;
    commentText = {};
    uibutton(grid, 'Text', 'Acknowledge', 'ButtonPushedFcn', @(~, ~) iClose(true));
    uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) iClose(false));
    dialog.CloseRequestFcn = @(~, ~) iClose(false);
    uiwait(dialog);
    % Delete now: the nested-function callbacks keep this workspace (and so
    % closeDialog) alive after return, and a lingering modal dialog would
    % hide the alerts that follow.
    iDeleteIfValid(dialog);

    if ~accepted
        result.message = 'Conflict acknowledgement cancelled.';
        return
    end
    UMITRigLedger.acknowledgeConflict(saveFolder, 'Origin', 'DataViewer', ...
        'Comment', strjoin(commentText, newline));
    result.changed = true;
    result.message = sprintf(['Rig conflict acknowledged: the folder now uses rig %s. ' ...
        'It stays marked "had conflict".'], summary.rigDisplayName);
    result.severity = 'warning';
catch ME
    result.message = sprintf('Could not acknowledge the rig conflict: %s', ME.message);
    result.severity = 'error';
    [~, lockMessage, isLockError] = offerStaleLockClearing(parentFigure, ME, 'Folder', saveFolder);
    if isLockError
        result.message = lockMessage;
        result.severity = 'warning';
    else
        uialert(parentFigure, ME.message, 'Acknowledge Rig Conflict failed', 'Icon', 'error');
    end
end

    function iClose(isAccepted)
        accepted = isAccepted;
        commentText = strtrim(cellstr(comment.Value));
        commentText = commentText(~cellfun(@isempty, commentText));
        uiresume(dialog);
    end
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
