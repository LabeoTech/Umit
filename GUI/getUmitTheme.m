function mode = getUmitTheme(varargin)
%GETUMITTHEME Return the persisted per-user umIT GUI theme.
%
%   mode = getUmitTheme()
%   mode = getUmitTheme('PreferencesFolder', folder)
%
%   The default preference is "light". The optional PreferencesFolder
%   input is intended for isolated tests and specialized deployments.

p = inputParser;
p.FunctionName = 'getUmitTheme';
addParameter(p, 'PreferencesFolder', '', @iIsTextScalar);
parse(p, varargin{:});

mode = 'light';

try
    preferences = readUmitPreferences(p.Results.PreferencesFolder, true);
catch ME
    warning('Umitoolbox:Preferences:ReadFailed', ...
        'Could not read umIT preferences; using light mode. %s', ...
        ME.message);
    return
end

if ~isfield(preferences, 'theme') || ~iIsTextScalar(preferences.theme)
    warning('Umitoolbox:Preferences:InvalidTheme', ...
        'The saved umIT theme is invalid; using light mode.');
    return
end

savedMode = lower(strtrim(char(preferences.theme)));
if ismember(savedMode, {'light', 'dark'})
    mode = savedMode;
else
    warning('Umitoolbox:Preferences:InvalidTheme', ...
        'The saved umIT theme "%s" is invalid; using light mode.', ...
        savedMode);
end

end

function tf = iIsTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
