function preferencesFile = setUmitTheme(mode, varargin)
%SETUMITTHEME Persist the per-user umIT GUI theme.
%
%   preferencesFile = setUmitTheme(mode)
%   preferencesFile = setUmitTheme(mode, 'PreferencesFolder', folder)

p = inputParser;
p.FunctionName = 'setUmitTheme';
addRequired(p, 'mode', @iIsTextScalar);
addParameter(p, 'PreferencesFolder', '', @iIsTextScalar);
parse(p, mode, varargin{:});

mode = lower(strtrim(char(p.Results.mode)));
if ~ismember(mode, {'light', 'dark'})
    error('Umitoolbox:setUmitTheme:InvalidTheme', ...
        'Theme must be "light" or "dark".');
end

preferences = readUmitPreferences(p.Results.PreferencesFolder, true);
preferences.schemaVersion = 1;
preferences.theme = mode;
preferencesFile = writeUmitPreferences( ...
    preferences, p.Results.PreferencesFolder);

end

function tf = iIsTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
