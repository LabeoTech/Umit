function saveMatAtomic(filePath, variableName, value)
%SAVEMATATOMIC Save one named variable to a MAT file atomically.
%
%   saveMatAtomic(filePath, variableName, value)
%
%   Writes the requested variable to a temporary file in the destination
%   folder, verifies that the temporary MAT file can be loaded, and then
%   replaces the destination. If replacement fails, the previous file is
%   restored whenever possible.
%
%   Inputs:
%       filePath     - Destination MAT-file path.
%       variableName - Valid MATLAB variable name stored in the file.
%       value        - Value assigned to the stored variable.
%
%   Notes:
%       - The temporary file is created beside the destination so the final
%         rename remains on the same filesystem.
%       - Existing files are moved to a temporary backup until the new file
%         has been installed successfully.

errID = 'Umitoolbox:saveMatAtomic:invalidInput';

if ~(ischar(filePath) || (isstring(filePath) && isscalar(filePath)))
    error(errID, '"filePath" must be a character vector or string scalar.');
end

if ~(ischar(variableName) || (isstring(variableName) && isscalar(variableName)))
    error(errID, '"variableName" must be a character vector or string scalar.');
end

filePath = char(string(filePath));
variableName = char(string(variableName));

if ~isvarname(variableName)
    error(errID, '"variableName" must be a valid MATLAB variable name.');
end

[parentFolder, ~, extension] = fileparts(filePath);
if isempty(parentFolder)
    parentFolder = pwd;
    filePath = fullfile(parentFolder, filePath);
end

if ~strcmpi(extension, '.mat')
    error(errID, 'Destination file must use the .mat extension.');
end

if ~isfolder(parentFolder)
    error(errID, 'Destination folder does not exist: %s', parentFolder);
end

tempPath = [tempname(parentFolder), '.mat'];
backupPath = [tempname(parentFolder), '.bak.mat'];
backupCreated = false;

cleanupObj = onCleanup(@() iDeleteIfPresent(tempPath));

payload = struct();
payload.(variableName) = value;
save(tempPath, '-mat', '-struct', 'payload', variableName);

loaded = load(tempPath, variableName, '-mat');
if ~isfield(loaded, variableName)
    error('Umitoolbox:saveMatAtomic:verificationFailed', ...
        'Temporary MAT file did not contain variable "%s".', variableName);
end

if isfile(filePath)
    [ok, message] = movefile(filePath, backupPath, 'f');
    if ~ok
        error('Umitoolbox:saveMatAtomic:backupFailed', ...
            'Could not back up existing file "%s": %s', filePath, message);
    end
    backupCreated = true;
end

[ok, message] = movefile(tempPath, filePath, 'f');
if ~ok
    if backupCreated && isfile(backupPath)
        movefile(backupPath, filePath, 'f');
    end
    error('Umitoolbox:saveMatAtomic:replaceFailed', ...
        'Could not install new file "%s": %s', filePath, message);
end

if backupCreated && isfile(backupPath)
    try
        delete(backupPath);
    catch ME
        warning('Umitoolbox:saveMatAtomic:backupCleanupFailed', ...
            'New file was saved, but temporary backup cleanup failed: %s', ...
            ME.message);
    end
end

clear cleanupObj

end

function iDeleteIfPresent(filePath)
%IDELETEIFPRESENT Delete a temporary file when it still exists.

if isfile(filePath)
    delete(filePath);
end

end
