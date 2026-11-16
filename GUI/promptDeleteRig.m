function result = promptDeleteRig(parentFigure, rigUUID, opts)
%PROMPTDELETERIG Rig Manager Delete Rig workflow (archived rigs only).
%
%   result = promptDeleteRig(parentFigure, rigUUID)
%   result = promptDeleteRig(parentFigure, rigUUID, 'Interactive', false, ...
%       'ConfirmRigID', typedID)
%
%   Permanently deletes an ARCHIVED rig with UMITRigStore.deleteRig after a
%   confirmation dialog that says what is lost, offers "Export Rig first..."
%   (promptExportRig, so a .umitrig copy can be imported again) and enables
%   "Delete permanently" only when the rig ID has been typed exactly. The
%   shared spectra and filter sets, and the rig ledgers of folders that were
%   bound to the rig, are not touched.
%
%   RESULT has fields changed (true when the rig was deleted), message,
%   severity ('info', 'success', 'warning' or 'error') and rigID. Errors are
%   shown in a uialert and reported in RESULT, never rethrown.
%
%   Name-Value options:
%       Interactive  - false: no dialog (tests). Default true.
%       ConfirmRigID - with Interactive false, the typed confirmation; the
%                      rig is deleted only when it equals the rig ID.

arguments
    parentFigure
    rigUUID (1,:) char
    opts.Interactive (1,1) logical = true
    opts.ConfirmRigID (1,:) char = ''
end

result = struct('changed', false, 'message', '', 'severity', 'info', 'rigID', '');
try
    store = UMITRigStore.open(rigUUID);
    rigInfo = store.getRigInfo();
    result.rigID = rigInfo.rigID;
    if ~strcmp(rigInfo.status, 'archived')
        result.message = sprintf(['Only an archived Rig can be deleted. Archive ' ...
            'Rig "%s" first.'], rigInfo.rigID);
        result.severity = 'warning';
        if opts.Interactive
            uialert(parentFigure, result.message, 'Delete Rig', 'Icon', 'warning');
        end
        return
    end

    if opts.Interactive
        confirmed = iAskConfirmation(parentFigure, rigInfo);
    else
        confirmed = strcmp(opts.ConfirmRigID, rigInfo.rigID);
    end
    if ~confirmed
        result.message = 'Rig deletion cancelled; the rig was left unchanged.';
        return
    end

    store.deleteRig();
    result.changed = true;
    result.severity = 'success';
    result.message = sprintf('Deleted rig "%s" (%s).', rigInfo.displayName, rigInfo.rigID);
catch ME
    result.message = sprintf('Could not delete the rig: %s', ME.message);
    result.severity = 'error';
    if opts.Interactive
        uialert(parentFigure, ME.message, 'Delete Rig failed', 'Icon', 'error');
    end
end
end

function confirmed = iAskConfirmation(parentFigure, rigInfo)
%IASKCONFIRMATION Modal dialog; true only for "Delete permanently" with the ID typed.
confirmed = false;
explanation = sprintf([ ...
    'Delete rig "%s" (%s) permanently?\n\n' ...
    'Deleted: the rig record and its camera-coregistration files and history.\n' ...
    'Not deleted: the shared spectra and filter sets, and the rig ledger of every ' ...
    'folder that used this rig (analysis results stay valid).\n\n' ...
    'Consequence: the dual-camera coregistration tool cannot save a calibration for a ' ...
    'folder bound to this rig until the rig is imported again.\n\n' ...
    'This cannot be undone. Use "Export Rig first..." to keep a .umitrig copy that ' ...
    'Import Rig can install again.'], rigInfo.displayName, rigInfo.rigID);

dialog = uifigure('Name', 'Delete Rig', 'WindowStyle', 'modal', ...
    'Position', [100 100 560 400], 'Resize', 'off');
closeDialog = onCleanup(@() iDeleteIfValid(dialog));
if ~isempty(parentFigure) && isvalid(parentFigure) && ...
        exist('placeAppInsideCaller', 'file') == 2
    try
        placeAppInsideCaller(parentFigure, dialog, 'center');
    catch
        movegui(dialog, 'center');
    end
else
    movegui(dialog, 'center');
end
try
    initializeUmitTheme(dialog);
catch
    % Theme is cosmetic.
end

grid = uigridlayout(dialog, [4 3], 'RowHeight', {'1x', 'fit', 30, 30}, ...
    'ColumnWidth', {'1x', '1x', '1x'});
textArea = uitextarea(grid, 'Value', splitlines(explanation), 'Editable', 'off');
textArea.Layout.Column = [1 3];
prompt = uilabel(grid, 'Text', sprintf('Type the Rig ID (%s) to confirm:', rigInfo.rigID));
prompt.Layout.Column = [1 3];
typed = uieditfield(grid, 'text', 'ValueChangingFcn', @(~, event) iTyping(event.Value));
typed.Layout.Column = [1 3];
exportButton = uibutton(grid, 'Text', 'Export Rig first...', ...
    'ButtonPushedFcn', @(~, ~) promptExportRig(dialog, 'RigUUID', rigInfo.uuid));
deleteButton = uibutton(grid, 'Text', 'Delete permanently', 'Enable', 'off', ...
    'FontColor', [0.7 0 0], 'ButtonPushedFcn', @(~, ~) iClose(true));
cancelButton = uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) iClose(false));
dialog.CloseRequestFcn = @(~, ~) iClose(false);
exportButton.Layout.Row = 4;
deleteButton.Layout.Row = 4;
cancelButton.Layout.Row = 4;

uiwait(dialog);
% Delete now: the nested-function callbacks keep this workspace (and so
% closeDialog) alive after return, and a lingering modal dialog would hide
% the alerts that follow.
iDeleteIfValid(dialog);

    function iTyping(text)
        deleteButton.Enable = matlab.lang.OnOffSwitchState(strcmp(text, rigInfo.rigID));
    end

    function iClose(isConfirmed)
        confirmed = isConfirmed;
        uiresume(dialog);
    end
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
