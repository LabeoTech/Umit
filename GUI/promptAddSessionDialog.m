function sessionInfo = promptAddSessionDialog(callerFigure, subjectInfo, saveFolder)
%PROMPTADDSESSIONDIALOG Prompt for new-session metadata against one SaveFolder.
%
%   sessionInfo = promptAddSessionDialog(callerFigure, subjectInfo, saveFolder)
%
%   Shared "Add Session" dialog (Session ID, display name, description,
%   acquisition date, read-only Rig preview) used by ProjectManagerTool and
%   ProjectBindingManager so both present the identical form. Returns []
%   when the user cancels; otherwise a struct with fields sessionID,
%   displayName, description, acquisitionDate (processedDataFolder is not
%   set here -- callers add it before calling UMITProjectStore.addSession).

sessionInfo = [];
inferredDate = getSessionAcquisitionDatePrefill(saveFolder);
[previewRigUUID, previewRigID] = ...
    UMITRigStore.readDatasetRigAssociation(saveFolder);
if isempty(previewRigID)
    try
        previewRigInfo = UMITRigStore.getActiveRig().getRigInfo();
        previewRigUUID = previewRigInfo.uuid;
        previewRigID = previewRigInfo.rigID;
    catch
        previewRigUUID = '';
        previewRigID = '';
    end
end
rigPreviewText = resolveRigDisplayName(previewRigUUID, previewRigID);

dialogFigure = uifigure( ...
    'Name', sprintf('Add Session to %s', subjectInfo.subjectID), ...
    'Visible', 'off', ...
    'WindowStyle', 'modal', ...
    'Position', [100, 100, 620, 440]);
cleanupFigure = onCleanup(@() iDeleteIfValid(dialogFigure));
dialogFigure.UserData = [];

grid = uigridlayout(dialogFigure, [7, 2]);
grid.RowHeight = {28, 28, 90, 28, 28, '1x', 32};
grid.ColumnWidth = {180, '1x'};
grid.Padding = [14, 14, 14, 14];
grid.RowSpacing = 8;
grid.ColumnSpacing = 8;

uilabel(grid, 'Text', 'Session ID:', 'HorizontalAlignment', 'right');
sessionIDField = uieditfield(grid, 'text');

displayLabel = uilabel(grid, 'Text', 'Display name:', ...
    'HorizontalAlignment', 'right');
displayLabel.Layout.Row = 2;
displayLabel.Layout.Column = 1;
displayNameField = uieditfield(grid, 'text');
displayNameField.Layout.Row = 2;
displayNameField.Layout.Column = 2;

descriptionLabel = uilabel(grid, 'Text', 'Description:', ...
    'HorizontalAlignment', 'right');
descriptionLabel.Layout.Row = 3;
descriptionLabel.Layout.Column = 1;
descriptionField = uitextarea(grid);
descriptionField.Layout.Row = 3;
descriptionField.Layout.Column = 2;

dateLabel = uilabel(grid, 'Text', 'Acquisition date:', ...
    'HorizontalAlignment', 'right');
dateLabel.Layout.Row = 4;
dateLabel.Layout.Column = 1;
dateField = uidatepicker(grid, ...
    'DisplayFormat', 'yyyy-MM-dd', ...
    'Value', inferredDate);
dateField.Layout.Row = 4;
dateField.Layout.Column = 2;

rigLabel = uilabel(grid, 'Text', 'Rig:', 'HorizontalAlignment', 'right');
rigLabel.Layout.Row = 5;
rigLabel.Layout.Column = 1;
rigField = uieditfield(grid, 'text', ...
    'Value', rigPreviewText, 'Editable', 'off', ...
    'Tooltip', ['Rig assignment is resolved automatically when the ' ...
    'session is created. Change it afterward via DataViewer''s ' ...
    'Assign Rig menu or RigManagerTool.']);
rigField.Layout.Row = 5;
rigField.Layout.Column = 2;

folderLabel = uilabel(grid, ...
    'Text', sprintf('SaveFolder: %s', saveFolder), ...
    'Tooltip', saveFolder);
folderLabel.Layout.Row = 6;
folderLabel.Layout.Column = [1, 2];

controls = struct( ...
    'sessionID', sessionIDField, ...
    'displayName', displayNameField, ...
    'description', descriptionField, ...
    'acquisitionDate', dateField);

addButton = uibutton(grid, 'push', ...
    'Text', 'Add Session', ...
    'ButtonPushedFcn', @(~, ~) iCompleteDialog(dialogFigure, controls));
addButton.Layout.Row = 7;
addButton.Layout.Column = 1;

cancelButton = uibutton(grid, 'push', ...
    'Text', 'Cancel', ...
    'ButtonPushedFcn', @(~, ~) iCompleteDialog(dialogFigure, []));
cancelButton.Layout.Row = 7;
cancelButton.Layout.Column = 2;

dialogFigure.CloseRequestFcn = @(~, ~) iCompleteDialog(dialogFigure, []);
placeAppInsideCaller(callerFigure, dialogFigure, 'center');
dialogFigure.Visible = 'on';
drawnow
figure(dialogFigure);
focus(sessionIDField);
uiwait(dialogFigure);

if isvalid(dialogFigure)
    sessionInfo = dialogFigure.UserData;
end
clear cleanupFigure
end

function iCompleteDialog(dialogFigure, controls)
if isempty(dialogFigure) || ~isvalid(dialogFigure)
    return
end

if isempty(controls)
    dialogFigure.UserData = [];
    uiresume(dialogFigure);
    return
end

sessionID = strtrim(char(string(controls.sessionID.Value)));
if isempty(sessionID)
    uialert(dialogFigure, 'Session ID cannot be empty.', ...
        'Invalid Session ID');
    focus(controls.sessionID);
    return
end

description = strjoin(string(controls.description.Value(:)), newline);
dialogFigure.UserData = struct( ...
    'sessionID', sessionID, ...
    'displayName', char(string(controls.displayName.Value)), ...
    'description', char(description), ...
    'acquisitionDate', controls.acquisitionDate.Value);
uiresume(dialogFigure);
end

function iDeleteIfValid(fig)
if ~isempty(fig) && isvalid(fig)
    delete(fig);
end
end
