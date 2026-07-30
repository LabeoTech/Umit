function mode = initializeUmitTheme(figureHandle, varargin)
%INITIALIZEUMITTHEME Initialize a uifigure from the saved umIT theme.
%
%   mode = initializeUmitTheme(figureHandle)
%   mode = initializeUmitTheme(figureHandle, 'ShowMenu', true)
%   mode = initializeUmitTheme(figureHandle, ...
%       'ThemeChangedFcn', @(mode) updateAppColors(mode))
%
%   Light mode leaves the App Designer appearance unchanged at startup.
%   Dark mode is applied immediately. When ShowMenu is true, a Theme menu
%   is added to the supplied figure. Menu selections affect this figure
%   only; other open apps adopt the preference when they next start. An
%   optional ThemeChangedFcn can update app-specific graphics that the
%   generic theme engine deliberately leaves untouched.

p = inputParser;
p.FunctionName = 'initializeUmitTheme';
addRequired(p, 'figureHandle', @iIsFigure);
addParameter(p, 'ShowMenu', false, ...
    @(value) islogical(value) && isscalar(value));
addParameter(p, 'PreferencesFolder', '', @iIsTextScalar);
addParameter(p, 'ThemeChangedFcn', [], ...
    @(value) isempty(value) || isa(value, 'function_handle'));
parse(p, figureHandle, varargin{:});

preferencesFolder = p.Results.PreferencesFolder;
if ~isempty(p.Results.ThemeChangedFcn)
    setappdata(figureHandle, 'UmitThemeChangedFcn', ...
        p.Results.ThemeChangedFcn);
end

mode = getUmitTheme('PreferencesFolder', preferencesFolder);
setappdata(figureHandle, 'UmitTheme', mode);

if strcmp(mode, 'dark')
    setGUIcolorScheme(figureHandle, mode);
end

iInvokeThemeChangedFcn(figureHandle, mode);

if p.Results.ShowMenu
    iEnsureThemeMenu(figureHandle, preferencesFolder, mode);
end

end

function iEnsureThemeMenu(figureHandle, preferencesFolder, mode)
themeMenu = findall(figureHandle, 'Type', 'uimenu', ...
    'Tag', 'UmitThemeMenu');

if isempty(themeMenu)
    themeMenu = uimenu(figureHandle, ...
        'Text', 'Theme', ...
        'Tag', 'UmitThemeMenu');
else
    themeMenu = themeMenu(1);
end

lightMenu = findall(themeMenu, 'Type', 'uimenu', ...
    'Tag', 'UmitLightThemeMenu');
if isempty(lightMenu)
    lightMenu = uimenu(themeMenu, ...
        'Text', 'Light', ...
        'Tag', 'UmitLightThemeMenu');
else
    lightMenu = lightMenu(1);
end

darkMenu = findall(themeMenu, 'Type', 'uimenu', ...
    'Tag', 'UmitDarkThemeMenu');
if isempty(darkMenu)
    darkMenu = uimenu(themeMenu, ...
        'Text', 'Dark', ...
        'Tag', 'UmitDarkThemeMenu');
else
    darkMenu = darkMenu(1);
end

lightMenu.MenuSelectedFcn = @(~, ~) iSelectTheme( ...
    figureHandle, 'light', preferencesFolder);
darkMenu.MenuSelectedFcn = @(~, ~) iSelectTheme( ...
    figureHandle, 'dark', preferencesFolder);

iUpdateMenuChecks(lightMenu, darkMenu, mode);
end

function iSelectTheme(figureHandle, mode, preferencesFolder)
if ~isvalid(figureHandle)
    return
end

try
    setUmitTheme(mode, 'PreferencesFolder', preferencesFolder);
catch ME
    warning('Umitoolbox:Preferences:WriteFailed', ...
        'Could not save the umIT theme: %s', ME.message);
    return
end

setGUIcolorScheme(figureHandle, mode);
iInvokeThemeChangedFcn(figureHandle, mode);

lightMenu = findall(figureHandle, 'Type', 'uimenu', ...
    'Tag', 'UmitLightThemeMenu');
darkMenu = findall(figureHandle, 'Type', 'uimenu', ...
    'Tag', 'UmitDarkThemeMenu');

if ~isempty(lightMenu) && ~isempty(darkMenu)
    iUpdateMenuChecks(lightMenu(1), darkMenu(1), mode);
end
end

function iInvokeThemeChangedFcn(figureHandle, mode)
if ~isappdata(figureHandle, 'UmitThemeChangedFcn')
    return
end

themeChangedFcn = getappdata(figureHandle, 'UmitThemeChangedFcn');
if ~isa(themeChangedFcn, 'function_handle')
    return
end

try
    themeChangedFcn(mode);
catch ME
    warning('Umitoolbox:Theme:CallbackFailed', ...
        'The app-specific theme callback failed: %s', ME.message);
end
end

function iUpdateMenuChecks(lightMenu, darkMenu, mode)
if strcmp(mode, 'dark')
    lightMenu.Checked = 'off';
    darkMenu.Checked = 'on';
else
    lightMenu.Checked = 'on';
    darkMenu.Checked = 'off';
end
end

function tf = iIsFigure(value)
tf = isscalar(value) && isgraphics(value) && ...
    isa(value, 'matlab.ui.Figure');
end

function tf = iIsTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
