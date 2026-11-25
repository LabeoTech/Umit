function result = ensureUMITPath(options)
%ENSUREUMITPATH Ensure only the production umIT folders are on the MATLAB path.
%
%   ensureUMITPath adds any missing production folders below the umIT
%   toolbox root, removes every excluded folder of the checkout from the
%   MATLAB search path, cleans excluded or vanished checkout entries out of
%   the saved path file (pathdef.m), and saves the path.
%
%   ensureUMITPath('Persist', false) fixes the current session only and
%   never saves the search path; it still reports what the saved path holds.
%
%   result = ensureUMITPath(...) returns a structure describing whether the
%   session path changed, which folders were missing or removed, which
%   excluded or vanished checkout entries the saved path file holds
%   (StaleSavedEntries), and the SAVEPATH outcome.
%
%   Excluded from the path (never added, removed when found):
%       - any folder whose name starts with "." (.git, .github, .claude, ...);
%       - docs, claude-context and test (test holds mocks that shadow MATLAB
%         built-ins such as memory);
%       - local underscore folders (_tasks, _logs, _apps_exported,
%         _branchComp, ...), except the two development tool folders _tools
%         and _mlapp_tools (validateTaskChanges, exportMatlabAppToDelta, ...).
%
%   Everything else below the toolbox root that GENPATH returns stays on
%   the path. Do not use addpath(genpath(<repo>)): it walks .git and docs.
%
%   One-time cleanup of a saved path that still names excluded folders or
%   folders that no longer exist (MATLAB prints "Name is nonexistent or not
%   a directory" at every startup): run ensureUMITPath once from the
%   repository root with the default Persist = true. Vanished folders are
%   already dropped from the session path at startup, so only the saved
%   file keeps them; saving the path rewrites it without them.

arguments
    options.Persist (1, 1) logical = true
    options.PathDefFile (1, 1) string = ""   % saved-path file to inspect (default: pathdef.m)
end

toolboxRoot = fileparts(mfilename('fullpath'));
allToolboxFolders = splitPathEntries(genpath(toolboxRoot));
isExcluded = cellfun(@(folder) isExcludedFolder(folder, toolboxRoot), ...
    allToolboxFolders);
requiredFolders = allToolboxFolders(~isExcluded);
currentFolders = splitPathEntries(path);

isInToolbox = cellfun(@(folder) isSameOrDescendant(folder, toolboxRoot), ...
    currentFolders);
excludedFoldersOnPath = currentFolders(isInToolbox & cellfun( ...
    @(folder) isExcludedFolder(folder, toolboxRoot), currentFolders));

if ispc
    isPresent = cellfun( ...
        @(folder) any(strcmpi(folder, currentFolders)), requiredFolders);
else
    isPresent = ismember(requiredFolders, currentFolders);
end

missingFolders = requiredFolders(~isPresent);
staleSavedEntries = findStaleSavedEntries(options.PathDefFile, toolboxRoot);

result = struct( ...
    'Changed', ~isempty(missingFolders) || ~isempty(excludedFoldersOnPath), ...
    'MissingFolders', {missingFolders}, ...
    'RemovedFolders', {excludedFoldersOnPath}, ...
    'StaleSavedEntries', {staleSavedEntries}, ...
    'SaveAttempted', false, ...
    'SaveStatus', NaN);

needsSave = result.Changed || ~isempty(staleSavedEntries);
if ~needsSave
    return
end

if ~isempty(excludedFoldersOnPath)
    rmpath(excludedFoldersOnPath{:});
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
        ['The umIT folders were fixed for this MATLAB session, but ' ...
         'MATLAB could not save the updated path permanently.']);
end
end

function stale = findStaleSavedEntries(pathDefFile, toolboxRoot)
%FINDSTALESAVEDENTRIES Checkout entries of the saved path that are excluded or vanished.
%
% Reads the saved path file as text (pathdef.m by default) and returns the entries
% below the toolbox root that name an excluded folder or a folder that no longer
% exists. MATLAB drops vanished folders from the session path at startup, so the
% saved file is the only place that still lists them.

stale = cell(1, 0);

pathDefFile = char(pathDefFile);
if isempty(pathDefFile)
    pathDefFile = which('pathdef');
end
if isempty(pathDefFile) || ~isfile(pathDefFile)
    return
end

text = fileread(pathDefFile);
tokens = regexp(text, '''([^'';]+);''', 'tokens');
for iToken = 1:numel(tokens)
    entry = tokens{iToken}{1};
    if ~isSameOrDescendant(entry, toolboxRoot)
        continue
    end
    if isExcludedFolder(entry, toolboxRoot) || ~isfolder(entry)
        stale{end + 1} = entry; %#ok<AGROW>
    end
end
end

function tf = isExcludedFolder(folder, toolboxRoot)
%ISEXCLUDEDFOLDER True for a folder of the checkout that must not be on the path.

relativeFolder = stripTrailingSeparators(folder);
rootFolder = stripTrailingSeparators(toolboxRoot);
if length(relativeFolder) <= length(rootFolder)
    tf = false;
    return
end
relativeFolder = relativeFolder(length(rootFolder) + 2:end);
segments = strsplit(relativeFolder, {'/', '\'});

developmentToolFolders = {'_tools', '_mlapp_tools'};
excludedTopLevelFolders = {'docs', 'claude-context', 'test'};

isDotFolder = any(startsWith(segments, '.'));
topLevel = segments{1};
isExcludedTopLevel = any(strcmpi(topLevel, excludedTopLevelFolders));
isLocalUnderscoreFolder = startsWith(topLevel, '_') && ...
    ~any(strcmpi(topLevel, developmentToolFolders));

tf = isDotFolder || isExcludedTopLevel || isLocalUnderscoreFolder;
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
