function [filterSetIDs, accepted] = promptSelectFilterSets(parentFigure, filterSets)
%PROMPTSELECTFILTERSETS Modal dialog to choose which filter sets an exported rig carries.
%
%   [filterSetIDs, accepted] = promptSelectFilterSets(parentFigure, filterSets)
%
%   FILTERSETS is a table with variables ID and DisplayName (the user-defined
%   rows of UMITRigStore.listFilterSets). Every set is pre-selected; the user
%   may deselect any, or all, and then export (ACCEPTED true) or cancel
%   (ACCEPTED false). FILTERSETIDS is a cell array of the chosen IDs, empty
%   when none is chosen or the dialog is cancelled. PARENTFIGURE positions the
%   dialog and may be [].

filterSetIDs = {};
accepted = false;
if isempty(filterSets) || height(filterSets) == 0
    accepted = true;
    return
end

ids = cellstr(string(filterSets.ID));
labels = cellstr(strcat(string(filterSets.DisplayName), " (", string(filterSets.ID), ")"));

dialog = uifigure('Name', 'Export Rig: filter sets', 'WindowStyle', 'modal', ...
    'Position', [100 100 480 320], 'Resize', 'off');
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

grid = uigridlayout(dialog, [3 2], 'RowHeight', {'fit', '1x', 30}, ...
    'ColumnWidth', {'1x', '1x'});
message = uilabel(grid, 'Text', ['Filter sets are not part of a rig. Choose the ' ...
    'filter sets (with the filter spectra they use) to include in the ' ...
    'exported file, so that they can be used on the other computer. ' ...
    'Built-in filter sets are not listed: every installation has them.'], ...
    'WordWrap', 'on');
message.Layout.Column = [1 2];
listbox = uilistbox(grid, 'Items', labels, 'ItemsData', ids, ...
    'Multiselect', 'on', 'Value', ids);
listbox.Layout.Column = [1 2];
uibutton(grid, 'Text', 'Export', 'ButtonPushedFcn', @(~, ~) iClose(true));
uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) iClose(false));
dialog.CloseRequestFcn = @(~, ~) iClose(false);

uiwait(dialog);
% Delete now: the nested-function callbacks keep this workspace (and so
% closeDialog) alive after return, and a lingering modal dialog would hide
% the caller's next dialog.
iDeleteIfValid(dialog);

    function iClose(isAccepted)
        accepted = isAccepted;
        if isAccepted
            chosen = listbox.Value;
            if ischar(chosen) || isstring(chosen)
                chosen = cellstr(chosen);
            end
            filterSetIDs = chosen(:)';
        end
        uiresume(dialog);
    end
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
