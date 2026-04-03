function DataParams = loadDataParams(folderPath)
%LOADDATAPARAMS Load DataParams.mat from a folder.
%
%   DataParams = loadDataParams(folderPath)
%
%   Loads the folder-global DataParams structure from DataParams.mat.
%
%   Input:
%       folderPath - Path to the save folder.
%
%   Output:
%       DataParams - Loaded DataParams structure.
%
%   Error:
%       Throws an error if the folder or DataParams.mat file does not exist,
%       or if the file does not contain the variable "DataParams".

if ~(ischar(folderPath) || isstring(folderPath))
    error('loadDataParams:InvalidFolderPath', ...
        'folderPath must be a char vector or string scalar.');
end

folderPath = char(string(folderPath));

if ~isfolder(folderPath)
    error('loadDataParams:InvalidFolder', ...
        'Folder does not exist: %s', folderPath);
end

filePath = fullfile(folderPath, 'DataParams.mat');

if ~exist(filePath, 'file')
    error('loadDataParams:MissingFile', ...
        'DataParams.mat was not found in folder: %s', folderPath);
end

S = load(filePath, 'DataParams');

if ~isfield(S, 'DataParams')
    error('loadDataParams:MissingVariable', ...
        'File "%s" does not contain variable "DataParams".', filePath);
end

DataParams = S.DataParams;
validateDataParams(DataParams);
end