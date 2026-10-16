function showProjectLedgerHistoryDialog(parentFigure, folder)
%SHOWPROJECTLEDGERHISTORYDIALOG Show a folder's project ledger history (read-only).
%
%   showProjectLedgerHistoryDialog(parentFigure, folder)
%
%   Displays formatProjectLedgerHistory in a modal window. Never writes. An
%   invalid ledger is reported in a uialert.
%
%   See also formatProjectLedgerHistory, showRigLedgerHistoryDialog.

try
    lines = formatProjectLedgerHistory(folder);
catch ME
    uialert(parentFigure, ME.message, 'Project Ledger History', 'Icon', 'error');
    return
end

dialog = uifigure('Name', 'Project Ledger History', 'WindowStyle', 'modal', ...
    'Position', [100 100 900 460]);
closeDialog = onCleanup(@() iDeleteIfValid(dialog));
movegui(dialog, 'center');
try
    initializeUmitTheme(dialog);
catch
    % Theme is cosmetic.
end
grid = uigridlayout(dialog, [2 1], 'RowHeight', {'1x', 'fit'});
uitextarea(grid, 'Value', lines, 'Editable', 'off', 'FontName', 'Consolas');
uibutton(grid, 'Text', 'Close', 'ButtonPushedFcn', @(~, ~) uiresume(dialog));
dialog.CloseRequestFcn = @(~, ~) uiresume(dialog);
uiwait(dialog);
iDeleteIfValid(dialog);
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
