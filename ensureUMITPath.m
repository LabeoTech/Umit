function result = ensureUMITPath(options)
%ENSUREUMITPATH Ensure the production umIT folder tree is on the MATLAB path.
%
%   ensureUMITPath adds any missing production folders below the umIT
%   toolbox root, removes test folders from the MATLAB search path, and
%   saves the updated path. Test folders are excluded because they may
%   contain mocks that intentionally shadow MATLAB built-ins.
%
%   ensureUMITPath('Persist', false) adds missing production folders and
%   removes test folders for the current session without saving the updated
%   search path.
%
%   result = ensureUMITPath(...) returns a structure describing whether the
%   path changed, which folders were missing or removed, and the SAVEPATH
%   outcome.

arguments
    options.Persist (1, 1) logical = true
end

toolboxRoot = fileparts(mfilename('fullpath'));
allToolboxFolders = splitPathEntries(genpath(toolboxRoot));
testRoot = fullfile(toolboxRoot, 'test');
isTestFolder = cellfun(@(folder) isSameOrDescendant(folder, testRoot), ...
    allToolboxFolders);
requiredFolders = allToolboxFolders(~isTestFolder);
currentFolders = splitPathEntries(path);

testFoldersOnPath = currentFolders(cellfun( ...
    @(folder) isSameOrDescendant(folder, testRoot), currentFolders));

if ispc
    isPresent = cellfun( ...
        @(folder) any(strcmpi(folder, currentFolders)), requiredFolders);
else
    isPresent = ismember(requiredFolders, currentFolders);
end

missingFolders = requiredFolders(~isPresent);
result = struct( ...
    'Changed', ~isempty(missingFolders) || ~isempty(testFoldersOnPath), ...
    'MissingFolders', {missingFolders}, ...
    'RemovedFolders', {testFoldersOnPath}, ...
    'SaveAttempted', false, ...
    'SaveStatus', NaN);

if ~result.Changed
    return
end

if ~isempty(testFoldersOnPath)
    rmpath(testFoldersOnPath{:});
end

if ~isempty(missingFolders)
    addpath(missingFolders{:});
end

if ~options.Persist
    return
end

result.SaveAttempted = true;
result.SaveStatus = savepath;
if result.SaveStatus ~= 0
    warning('umIT:Path:SaveFailed', ...
        ['The umIT folders were added for this MATLAB session, but ' ...
         'MATLAB could not save the updated path permanently.']);
end
end

function entries = splitPathEntries(pathValue)
entries = strsplit(pathValue, pathsep);
entries = entries(~cellfun('isempty', entries));
end

function tf = isSameOrDescendant(folder, rootFolder)
folder = stripTrailingSeparators(folder);
rootFolder = stripTrailingSeparators(rootFolder);
rootPrefix = [rootFolder filesep];

if ispc
    tf = strcmpi(folder, rootFolder) || ...
        strncmpi(folder, rootPrefix, length(rootPrefix));
else
    tf = strcmp(folder, rootFolder) || ...
        strncmp(folder, rootPrefix, length(rootPrefix));
end
end

function folder = stripTrailingSeparators(folder)
while length(folder) > 1 && any(folder(end) == ['/' '\'])
    folder(end) = [];
end
end
