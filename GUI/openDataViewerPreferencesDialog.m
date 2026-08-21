function dialogFigure = openDataViewerPreferencesDialog( ...
    parentFigure, colormapNames, colormapKeys, applyFunction, varargin)
%OPENDATAVIEWERPREFERENCESDIALOG Open the compact DataViewer settings UI.
%
%   openDataViewerPreferencesDialog(parent, names, keys, applyFunction)
%   loads persisted user-level preferences. OK saves the edited values and
%   invokes applyFunction(preferences). Cancel and the window close button
%   discard edits.

if nargin < 4
    applyFunction = [];
end

parser = inputParser;
parser.FunctionName = 'openDataViewerPreferencesDialog';
addRequired(parser, 'parentFigure', @iIsFigure);
addRequired(parser, 'colormapNames', @iIsTextCollection);
addRequired(parser, 'colormapKeys', @iIsTextCollection);
addRequired(parser, 'applyFunction', ...
    @(value) isempty(value) || isa(value, 'function_handle'));
addParameter(parser, 'PreferencesFolder', '', @iIsTextScalar);
parse(parser, parentFigure, colormapNames, colormapKeys, applyFunction, ...
    varargin{:});

colormapNames = cellstr(string(colormapNames(:)));
colormapKeys = cellstr(string(colormapKeys(:)));
if isempty(colormapKeys) || numel(colormapNames) ~= numel(colormapKeys)
    error('Umitoolbox:DataViewerPreferences:InvalidColormaps', ...
        'Colormap names and keys must be nonempty collections of equal size.');
end

preferencesFolder = char(string(parser.Results.PreferencesFolder));
preferences = DataViewerPreferences.load( ...
    'AvailableColormaps', colormapKeys, ...
    'PreferencesFolder', preferencesFolder);

existingDialog = findall(groot, 'Type', 'figure', ...
    'Tag', 'DataViewerPreferencesDialog');
if ~isempty(existingDialog)
    figure(existingDialog(1));
    dialogFigure = existingDialog(1);
    return
end

dialogFigure = uifigure( ...
    'Name', 'DataViewer Preferences', ...
    'Tag', 'DataViewerPreferencesDialog', ...
    'Resize', 'off', ...
    'WindowStyle', 'modal', ...
    'Position', iCenteredPosition(parentFigure, [520 281]));

layout = uigridlayout(dialogFigure, [7 3]);
layout.RowHeight = {28, 28, 28, 28, 28, 32, 32};
layout.ColumnWidth = {145, '1x', 85};
layout.Padding = [14 14 14 12];
layout.RowSpacing = 8;
layout.ColumnSpacing = 8;

themeLabel = uilabel(layout, 'Text', 'Theme');
themeLabel.Layout.Row = 1;
themeLabel.Layout.Column = 1;

themeDropDown = uidropdown(layout, ...
    'Items', {'Light', 'Dark'}, ...
    'ItemsData', {'light', 'dark'}, ...
    'Value', preferences.theme, ...
    'Tag', 'ThemeDropDown');
themeDropDown.Layout.Row = 1;
themeDropDown.Layout.Column = [2 3];

colormapLabel = uilabel(layout, 'Text', 'Default colormap');
colormapLabel.Layout.Row = 2;
colormapLabel.Layout.Column = 1;

colormapDropDown = uidropdown(layout, ...
    'Items', colormapNames, ...
    'ItemsData', colormapKeys, ...
    'Value', preferences.defaultColormap, ...
    'Tag', 'DefaultColormapDropDown');
colormapDropDown.Layout.Row = 2;
colormapDropDown.Layout.Column = [2 3];

reopenCheckBox = uicheckbox(layout, ...
    'Text', 'Reopen the last opened file on startup', ...
    'Value', preferences.reopenLastFile, ...
    'Tag', 'ReopenLastFileCheckBox');
reopenCheckBox.Layout.Row = 3;
reopenCheckBox.Layout.Column = [1 3];

rememberCheckBox = uicheckbox(layout, ...
    'Text', 'Remember the last SaveFolder for file dialogs', ...
    'Value', preferences.rememberLastSaveFolder, ...
    'Tag', 'RememberLastSaveFolderCheckBox');
rememberCheckBox.Layout.Row = 4;
rememberCheckBox.Layout.Column = [1 3];

folderLabel = uilabel(layout, 'Text', 'Default folder');
folderLabel.Layout.Row = 5;
folderLabel.Layout.Column = 1;

folderField = uieditfield(layout, 'text', ...
    'Value', preferences.defaultFolder, ...
    'Editable', 'off', ...
    'Tag', 'DefaultFolderField');
folderField.Layout.Row = 5;
folderField.Layout.Column = [2 3];

browseButton = uibutton(layout, 'push', ...
    'Text', 'Browse...', ...
    'Tag', 'BrowseDefaultFolderButton', ...
    'ButtonPushedFcn', @browseForDefaultFolder);
browseButton.Layout.Row = 6;
browseButton.Layout.Column = 3;

resetButton = uibutton(layout, 'push', ...
    'Text', 'Reset', ...
    'Tag', 'ResetDataViewerPreferencesButton', ...
    'ButtonPushedFcn', @resetControls);
resetButton.Layout.Row = 7;
resetButton.Layout.Column = 1;

cancelButton = uibutton(layout, 'push', ...
    'Text', 'Cancel', ...
    'Tag', 'CancelDataViewerPreferencesButton', ...
    'ButtonPushedFcn', @(~, ~) delete(dialogFigure));
cancelButton.Layout.Row = 7;
cancelButton.Layout.Column = 2;

okButton = uibutton(layout, 'push', ...
    'Text', 'OK', ...
    'Tag', 'SaveDataViewerPreferencesButton', ...
    'ButtonPushedFcn', @saveAndClose);
okButton.Layout.Row = 7;
okButton.Layout.Column = 3;

initializeUmitTheme(dialogFigure);

    function browseForDefaultFolder(~, ~)
        startFolder = folderField.Value;
        if ~isfolder(startFolder)
            startFolder = DataViewerPreferences.resolveDialogFolder( ...
                pwd, preferences);
        end

        selectedFolder = uigetdir(startFolder, 'Select default folder');
        if ~isequal(selectedFolder, 0)
            folderField.Value = selectedFolder;
        end
    end

    function resetControls(~, ~)
        defaults = DataViewerPreferences.defaults();
        if ~ismember(defaults.defaultColormap, colormapKeys)
            defaults.defaultColormap = colormapKeys{1};
        end
        themeDropDown.Value = defaults.theme;
        colormapDropDown.Value = defaults.defaultColormap;
        reopenCheckBox.Value = defaults.reopenLastFile;
        rememberCheckBox.Value = defaults.rememberLastSaveFolder;
        folderField.Value = defaults.defaultFolder;
    end

    function saveAndClose(~, ~)
        edited = preferences;
        edited.theme = themeDropDown.Value;
        edited.defaultColormap = colormapDropDown.Value;
        edited.reopenLastFile = reopenCheckBox.Value;
        edited.rememberLastSaveFolder = rememberCheckBox.Value;
        edited.defaultFolder = folderField.Value;

        try
            [~, saved] = DataViewerPreferences.save(edited, ...
                'AvailableColormaps', colormapKeys, ...
                'PreferencesFolder', preferencesFolder);
        catch ME
            uialert(dialogFigure, ME.message, ...
                'Unable to Save Preferences', 'Icon', 'error');
            return
        end

        if ~isempty(applyFunction)
            try
                applyFunction(saved);
            catch ME
                warning('Umitoolbox:DataViewerPreferences:ApplyFailed', ...
                    'Preferences were saved but could not be applied: %s', ...
                    ME.message);
            end
        end
        delete(dialogFigure);
    end
end

function position = iCenteredPosition(parentFigure, dialogSize)
parentPosition = parentFigure.Position;
position = [ ...
    parentPosition(1) + (parentPosition(3) - dialogSize(1)) / 2, ...
    parentPosition(2) + (parentPosition(4) - dialogSize(2)) / 2, ...
    dialogSize];
end

function tf = iIsFigure(value)
tf = isscalar(value) && isgraphics(value) && ...
    isa(value, 'matlab.ui.Figure');
end

function tf = iIsTextCollection(value)
tf = ischar(value) || isstring(value) || iscellstr(value);
end

function tf = iIsTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
