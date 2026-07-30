function result = ensureUMITPath(options)
%ENSUREUMITPATH Ensure the complete umIT folder tree is on the MATLAB path.
%
%   ensureUMITPath adds any missing folders below the umIT toolbox root and
%   saves the updated MATLAB search path. If the path is already complete,
%   it makes no changes and does not call SAVEPATH.
%
%   ensureUMITPath('Persist', false) adds missing folders for the current
%   MATLAB session without saving the updated search path.
%
%   result = ensureUMITPath(...) returns a structure describing whether the
%   path changed, which folders were missing, and the SAVEPATH outcome.

arguments
    options.Persist (1, 1) logical = true
end

toolboxRoot = fileparts(mfilename('fullpath'));
requiredFolders = splitPathEntries(genpath(toolboxRoot));
currentFolders = splitPathEntries(path);

if ispc
    isPresent = cellfun( ...
        @(folder) any(strcmpi(folder, currentFolders)), requiredFolders);
else
    isPresent = ismember(requiredFolders, currentFolders);
end

missingFolders = requiredFolders(~isPresent);
result = struct( ...
    'Changed', ~isempty(missingFolders), ...
    'MissingFolders', {missingFolders}, ...
    'SaveAttempted', false, ...
    'SaveStatus', NaN);

if isempty(missingFolders)
    return
end

addpath(missingFolders{:});

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
