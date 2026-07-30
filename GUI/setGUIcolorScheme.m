function setGUIcolorScheme(parentHandle, mode)
%SETGUICOLORSCHEME Apply an umIT light or dark palette to a GUI subtree.
%
%   setGUIcolorScheme(parentHandle, mode)
%
%   The function updates UI component colors without changing plotted data
%   colors. Icon and image paths whose filename stem ends in "Black" or
%   "White" are switched to the matching theme when the paired file exists.

if ~(isscalar(parentHandle) && isgraphics(parentHandle))
    error('Umitoolbox:setGUIcolorScheme:InvalidHandle', ...
        'parentHandle must be a valid scalar graphics handle.');
end

mode = iNormalizeMode(mode);
palette = iLoadPalette(mode);
allHandles = findall(parentHandle);

for idx = 1:numel(allHandles)
    graphicHandle = allHandles(idx);

    if ~isvalid(graphicHandle)
        continue
    end

    className = class(graphicHandle);
    if iIsExcludedElement(className)
        continue
    end

    iSwitchTaggedAsset(graphicHandle, 'Icon', mode);
    iSwitchTaggedAsset(graphicHandle, 'ImageSource', mode);

    if ~startsWith(className, 'matlab.ui.')
        continue
    end

    iApplyComponentColors(graphicHandle, className, palette);
end

figureHandle = iFindParentFigure(parentHandle);
if ~isempty(figureHandle)
    setappdata(figureHandle, 'UmitTheme', mode);
end

end

function mode = iNormalizeMode(mode)
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('Umitoolbox:setGUIcolorScheme:InvalidTheme', ...
        'Theme must be "light" or "dark".');
end

mode = lower(strtrim(char(mode)));
if ~ismember(mode, {'light', 'dark'})
    error('Umitoolbox:setGUIcolorScheme:InvalidTheme', ...
        'Theme must be "light" or "dark".');
end
end

function palette = iLoadPalette(mode)
paletteFile = fullfile(fileparts(mfilename('fullpath')), ...
    'GUIcolorScheme.json');

try
    colorScheme = jsondecode(fileread(paletteFile));
catch ME
    error('Umitoolbox:setGUIcolorScheme:InvalidPalette', ...
        'Could not read the GUI color palette: %s', ME.message);
end

if strcmp(mode, 'dark')
    palette = colorScheme.DarkMode;
else
    palette = colorScheme.LightMode;
end
end

function tf = iIsExcludedElement(className)
excludedClasses = {
    'matlab.ui.container.Menu'
    'matlab.ui.control.Table'
    };

tf = ismember(className, excludedClasses);
end

function iApplyComponentColors(graphicHandle, className, palette)
if isa(graphicHandle, 'matlab.ui.Figure')
    iSafeSet(graphicHandle, 'Color', palette.FigureBackground);
    return
end

if isa(graphicHandle, 'matlab.ui.control.UIAxes')
    iApplyAxesColors(graphicHandle, palette);
    return
end

if isa(graphicHandle, 'matlab.ui.container.Tab') || ...
        isa(graphicHandle, 'matlab.ui.container.TabGroup')
    % Tab and tab-group foreground/font colors are app-defined and must
    % remain untouched. Only apply a background when MATLAB exposes one.
    iSafeSet(graphicHandle, 'BackgroundColor', ...
        palette.ContainerBackground);
    return
end

if startsWith(className, 'matlab.ui.container.')
    background = palette.ContainerBackground;
elseif isa(graphicHandle, 'matlab.ui.control.Image')
    background = palette.DataBackground;
elseif iUsesTransparentBackground(className)
    background = 'none';
else
    background = palette.ControlBackground;
end

iSafeSet(graphicHandle, 'BackgroundColor', background);
iSafeSet(graphicHandle, 'ForegroundColor', palette.Foreground);
iSafeSet(graphicHandle, 'FontColor', palette.Foreground);
iSafeSet(graphicHandle, 'BorderColor', palette.Border);
iSafeSet(graphicHandle, 'HighlightColor', palette.Border);
iSafeSet(graphicHandle, 'ShadowColor', palette.Border);
end

function tf = iUsesTransparentBackground(className)
transparentClasses = {
    'matlab.ui.control.Label'
    'matlab.ui.control.CheckBox'
    'matlab.ui.control.RadioButton'
    'matlab.ui.control.Slider'
    'matlab.ui.control.Hyperlink'
    };

tf = ismember(className, transparentClasses);
end

function iApplyAxesColors(axesHandle, palette)
iSafeSet(axesHandle, 'Color', palette.DataBackground);
iSafeSet(axesHandle, 'XColor', palette.Foreground);
iSafeSet(axesHandle, 'YColor', palette.Foreground);
iSafeSet(axesHandle, 'ZColor', palette.Foreground);
iSafeSet(axesHandle, 'GridColor', palette.Grid);
iSafeSet(axesHandle, 'MinorGridColor', palette.MinorGrid);

iSafeSetTextColor(axesHandle.Title, palette.Foreground);
iSafeSetTextColor(axesHandle.XLabel, palette.Foreground);
iSafeSetTextColor(axesHandle.YLabel, palette.Foreground);
iSafeSetTextColor(axesHandle.ZLabel, palette.Foreground);
end

function iSafeSetTextColor(textHandle, color)
if ~isempty(textHandle) && isvalid(textHandle)
    iSafeSet(textHandle, 'Color', color);
end
end

function iSafeSet(graphicHandle, propertyName, value)
if ~isprop(graphicHandle, propertyName)
    return
end

try
    graphicHandle.(propertyName) = value;
catch
    % Some R2024b component properties accept a narrower set of values.
    % A single unsupported property must not prevent the remaining GUI
    % elements from being themed.
end
end

function iSwitchTaggedAsset(graphicHandle, propertyName, mode)
if ~isprop(graphicHandle, propertyName)
    return
end

try
    currentAsset = graphicHandle.(propertyName);
catch
    return
end

if isstring(currentAsset)
    if ~isscalar(currentAsset)
        return
    end
    currentAsset = char(currentAsset);
elseif ~(ischar(currentAsset) && isrow(currentAsset))
    return
end

if isempty(currentAsset)
    return
end

[assetFolder, assetStem, assetExtension] = fileparts(currentAsset);
currentSuffix = regexpi(assetStem, '(black|white)$', 'match', 'once');
if isempty(currentSuffix)
    return
end

if strcmp(mode, 'dark')
    desiredSuffix = 'White';
else
    desiredSuffix = 'Black';
end

if strcmpi(currentSuffix, desiredSuffix)
    return
end

newStem = [assetStem(1:end-numel(currentSuffix)) desiredSuffix];
newAsset = iResolvePairedAsset(assetFolder, newStem, assetExtension);
if isempty(newAsset)
    return
end

try
    graphicHandle.(propertyName) = newAsset;
catch
    % Leave the current asset unchanged when the component rejects a path.
end
end

function assetPath = iResolvePairedAsset(assetFolder, assetStem, assetExtension)
assetName = [assetStem assetExtension];
assetPath = fullfile(assetFolder, assetName);

if isfile(assetPath)
    return
end

if isempty(assetFolder)
    pathMatch = which(assetName);
    if ~isempty(pathMatch)
        assetPath = pathMatch;
        return
    end
elseif isfolder(assetFolder)
    folderEntries = dir(assetFolder);
    entryNames = {folderEntries(~[folderEntries.isdir]).name};
    matchIdx = find(strcmpi(entryNames, assetName), 1);
    if ~isempty(matchIdx)
        assetPath = fullfile(assetFolder, entryNames{matchIdx});
        return
    end
end

assetPath = '';
end

function figureHandle = iFindParentFigure(graphicHandle)
if isa(graphicHandle, 'matlab.ui.Figure')
    figureHandle = graphicHandle;
else
    figureHandle = ancestor(graphicHandle, 'figure');
end
end
