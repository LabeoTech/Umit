function rigUUID = promptSelectLibraryRig(parentFigure, varargin)
%PROMPTSELECTLIBRARYRIG Modal dialog to pick one Active or Available library rig.
%
%   rigUUID = promptSelectLibraryRig(parentFigure)
%   rigUUID = promptSelectLibraryRig(parentFigure, 'Title', title, ...
%       'Message', message, 'PreselectUUID', rigUUID)
%
%   Lists the readable Active and Available rigs of the rig library (see
%   buildLibraryRigChoices). PreselectUUID is selected when offered,
%   otherwise the Active rig. Returns the chosen rig UUID, or '' when the
%   user cancels or the library has no rig to offer (the caller decides
%   what to report). PARENTFIGURE positions the dialog and may be [] or
%   still hidden; the empty-library alert is shown only on a visible one.

p = inputParser;
p.FunctionName = 'promptSelectLibraryRig';
addParameter(p, 'Title', 'Choose a rig', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Message', 'Choose a library rig.', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'PreselectUUID', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, varargin{:});

rigUUID = '';
[labels, rigUUIDs, defaultIndex] = buildLibraryRigChoices( ...
    UMITRigStore.listRigs(), p.Results.PreselectUUID);
if isempty(rigUUIDs)
    % uialert errors on a hidden figure (e.g. an app still in its
    % startupFcn); the caller reports the empty selection then.
    if ~isempty(parentFigure) && isvalid(parentFigure) && ...
            strcmp(parentFigure.Visible, 'on')
        uialert(parentFigure, ['The rig library has no Active or Available rig. ' ...
            'Create or import one in the Rig Manager.'], char(p.Results.Title), ...
            'Icon', 'warning');
    end
    return
end

dialog = uifigure('Name', char(p.Results.Title), 'WindowStyle', 'modal', ...
    'Position', [100 100 460 190], 'Resize', 'off');
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

grid = uigridlayout(dialog, [3 2], 'RowHeight', {'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'});
message = uilabel(grid, 'Text', char(p.Results.Message), 'WordWrap', 'on');
message.Layout.Column = [1 2];
dropdown = uidropdown(grid, 'Items', labels, 'ItemsData', 1:numel(labels), ...
    'Value', defaultIndex);
dropdown.Layout.Column = [1 2];
uibutton(grid, 'Text', 'OK', 'ButtonPushedFcn', @(~, ~) iChoose(true));
uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) iChoose(false));
dialog.CloseRequestFcn = @(~, ~) iChoose(false);

uiwait(dialog);
% Delete now: the nested-function callbacks keep this workspace (and so
% closeDialog) alive after return, and a lingering modal dialog would hide
% the caller's next uiconfirm.
iDeleteIfValid(dialog);

    function iChoose(accepted)
        if accepted
            rigUUID = rigUUIDs{dropdown.Value};
        end
        uiresume(dialog);
    end
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
