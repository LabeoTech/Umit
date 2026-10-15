function result = runFolderRigBindingDialog(parentFigure, saveFolder)
%RUNFOLDERRIGBINDINGDIALOG DataViewer Bind/Rebind Rig workflow.
%
%   result = runFolderRigBindingDialog(parentFigure, saveFolder)
%
%   1) The user picks a library rig (promptSelectLibraryRig, Active rig
%      pre-selected; the rig library's only role here is that
%      pre-selection).
%   2) planFolderRigBinding previews the change and, BEFORE confirmation,
%      which committed parts would conflict.
%   3) On confirmation, applyFolderRigBinding writes the ledger.
%   RESULT has fields changed (logical), message and severity ('info',
%   'warning' or 'error') for DataViewer's status bar. Errors are shown in
%   a uialert and reported in RESULT, never rethrown.

result = struct('changed', false, 'message', '', 'severity', 'info');
saveFolder = char(string(saveFolder));
try
    summary = describeFolderRigLedger(saveFolder);
    if ~summary.canBind
        result.message = sprintf('The rig of this folder cannot be changed: %s', ...
            strjoin(summary.tooltipLines, ' '));
        result.severity = 'warning';
        uialert(parentFigure, result.message, 'Bind / Rebind Rig', 'Icon', 'warning');
        return
    end

    rigUUID = promptSelectLibraryRig(parentFigure, ...
        'Title', 'Bind / Rebind Rig', ...
        'Message', sprintf('Choose the library rig for this folder. Current rig: %s.', ...
        summary.displayText));
    if isempty(rigUUID)
        result.message = 'Rig binding cancelled.';
        return
    end

    plan = planFolderRigBinding(saveFolder, rigUUID);
    if strcmp(plan.action, 'none')
        result.message = plan.confirmText;
        return
    end
    icon = 'question';
    if plan.conflict.hasConflict
        icon = 'warning';
    end
    choice = uiconfirm(parentFigure, plan.confirmText, plan.confirmTitle, ...
        'Options', {plan.confirmButton, 'Cancel'}, ...
        'DefaultOption', 2, 'CancelOption', 2, 'Icon', icon);
    if ~strcmp(choice, plan.confirmButton)
        result.message = 'Rig binding cancelled.';
        return
    end

    applied = applyFolderRigBinding(saveFolder, rigUUID);
    result.changed = true;
    result.message = applied.message;
    if applied.conflict.hasConflict
        result.severity = 'warning';
    end
catch ME
    result.message = sprintf('Could not change the folder rig: %s', ME.message);
    result.severity = 'error';
    uialert(parentFigure, ME.message, 'Bind / Rebind Rig failed', 'Icon', 'error');
end
end
