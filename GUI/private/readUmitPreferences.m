function preferences = readUmitPreferences(preferencesFolder, createIfMissing)
%READUMITPREFERENCES Read preferences, falling back to safe defaults.

if nargin < 1
    preferencesFolder = '';
end
if nargin < 2
    createIfMissing = true;
end

preferences = defaultUmitPreferences();
[preferencesFile, ~] = resolveUmitPreferencesFile( ...
    preferencesFolder, createIfMissing);

if isempty(preferencesFile)
    return
end

if ~isfile(preferencesFile)
    if createIfMissing
        try
            writeUmitPreferences(preferences, preferencesFolder);
        catch ME
            warning('Umitoolbox:Preferences:WriteDefaultsFailed', ...
                'Could not create the default preferences file: %s', ...
                ME.message);
        end
    end
    return
end

try
    decoded = jsondecode(fileread(preferencesFile));
    if ~(isstruct(decoded) && isscalar(decoded))
        error('Preferences JSON must contain one object.');
    end
    preferences = decoded;
catch ME
    warning('Umitoolbox:Preferences:InvalidFile', ...
        'Could not parse "%s"; using defaults. %s', ...
        preferencesFile, ME.message);
end

end
