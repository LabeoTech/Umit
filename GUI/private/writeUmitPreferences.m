function preferencesFile = writeUmitPreferences( ...
    preferences, preferencesFolder)
%WRITEUMITPREFERENCES Atomically write the umIT preference JSON file.

if nargin < 2
    preferencesFolder = '';
end

if ~(isstruct(preferences) && isscalar(preferences))
    error('Umitoolbox:Preferences:InvalidData', ...
        'Preferences must be a scalar structure.');
end

[preferencesFile, resolvedFolder] = resolveUmitPreferencesFile( ...
    preferencesFolder, true);
temporaryFile = tempname(resolvedFolder);
cleanupTemporary = onCleanup(@() iDeleteIfPresent(temporaryFile));

try
    jsonText = jsonencode(preferences, 'PrettyPrint', true);
catch ME
    error('Umitoolbox:Preferences:EncodeFailed', ...
        'Could not encode preferences: %s', ME.message);
end

fileId = fopen(temporaryFile, 'w', 'n', 'UTF-8');
if fileId < 0
    error('Umitoolbox:Preferences:OpenFailed', ...
        'Could not open a temporary preferences file for writing.');
end

cleanupFile = onCleanup(@() iCloseIfOpen(fileId));
writeCount = fwrite(fileId, [jsonText newline], 'char');
if writeCount ~= numel([jsonText newline])
    error('Umitoolbox:Preferences:WriteFailed', ...
        'Could not write the complete preferences file.');
end

fclose(fileId);
clear cleanupFile

[moved, message] = movefile(temporaryFile, preferencesFile, 'f');
if ~moved
    error('Umitoolbox:Preferences:ReplaceFailed', ...
        'Could not replace "%s": %s', preferencesFile, message);
end

clear cleanupTemporary

end

function iCloseIfOpen(fileId)
if fileId > 0
    try
        fclose(fileId);
    catch
    end
end
end

function iDeleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end
