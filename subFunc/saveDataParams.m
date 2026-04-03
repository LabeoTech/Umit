function saveDataParams(folderPath, DataParams)
%SAVEDATAPARAMS Validate and save DataParams.mat to a folder.
%
%   saveDataParams(folderPath, DataParams)
%
%   Validates the supplied DataParams structure, updates the
%   "lastModified" timestamp, and saves it as DataParams.mat in the
%   specified folder.
%
%   Inputs:
%       folderPath  - Path to the folder.
%       DataParams  - Folder-global DataParams structure.
%
%   Notes:
%       - The folder must already exist.
%       - The variable is saved with the name "DataParams".

if ~(ischar(folderPath) || isstring(folderPath))
    error('saveDataParams:InvalidFolderPath', ...
        'folderPath must be a char vector or string scalar.');
end

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('saveDataParams:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

validateDataParams(DataParams);
DataParams.lastModified = datetime('now');

filePath = fullfile(folderPath, 'DataParams.mat');
save(filePath, 'DataParams', '-mat');
end