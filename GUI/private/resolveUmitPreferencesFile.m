function [preferencesFile, preferencesFolder] = ...
    resolveUmitPreferencesFile(preferencesFolder, createFolder)
%RESOLVEUMITPREFERENCESFILE Resolve the per-user preference JSON path.

if nargin < 1 || isempty(preferencesFolder)
    preferencesFolder = getUmitFolder('preferences', ...
        'create', createFolder);
else
    preferencesFolder = char(preferencesFolder);

    if ~isfolder(preferencesFolder)
        if createFolder
            [created, message] = mkdir(preferencesFolder);
            if ~created
                error('Umitoolbox:Preferences:CreateFolderFailed', ...
                    'Could not create preferences folder "%s": %s', ...
                    preferencesFolder, message);
            end
        else
            preferencesFolder = '';
        end
    end
end

if isempty(preferencesFolder)
    preferencesFile = '';
else
    preferencesFile = fullfile(preferencesFolder, ...
        'appPreferences.json');
end

end
