function out = getFileList(folder, extension)
% GETFILELIST retrieves a list of relevant files in a folder based on extension.
%
% Inputs:
%   folder (char): The path to the folder containing the files.
%   extension (char): The desired file extension or 'all' to list all files.
%
% Output:
%   out (cell): List of valid filenames with the specified extension.
%
% Example usage:
%   fileList = getFileList('/path/to/folder', '.mat');
%
% Note: This function does not change the algorithm for file retrieval.



% Get a list of .mat files
matFiles = dir(fullfile(folder, '*.mat'));

% Check if "datSize" exists inside each .mat file
warning('off'); % Disable warning to suppress messages from .mat files that do not have "datSize" variable.
validMatFiles = false(size(matFiles));
for ii = 1:length(matFiles)
    data = load(fullfile(folder, matFiles(ii).name), 'dim_names');
    validMatFiles(ii) = ~isempty(fieldnames(data));
end
warning('on');

fileList = {};
% if any(strcmpi(extension, {'.mat', '.dat'})) && ~any(validMatFiles)
%     return
% end

if ischar(extension)
    extension = {extension};
end

for ii = 1:length(extension)
    ext = extension{ii};
    matFileNames = {matFiles(validMatFiles).name}';
    [~, matNames, ~] = cellfun(@(x) fileparts(x), matFileNames, 'UniformOutput', false);
    
    % Get a list of .dat files
    datFiles = dir(fullfile(folder, '*.dat'));
    [~, datNames, ~] = arrayfun(@(x) fileparts(x.name), datFiles, 'UniformOutput', false);
    
    switch ext
        case '.mat'
            fileList{ii} = cellfun(@(x) [x, '.mat'], setdiff(matNames, datNames), 'UniformOutput', false);
        case '.dat'
            fileList{ii} = cellfun(@(x) [x, '.dat'], intersect(matNames, datNames), 'UniformOutput', false);
        case 'all'
            matList = cellfun(@(x) [x, '.mat'], setdiff(matNames, datNames), 'UniformOutput', false);
            datList = cellfun(@(x) [x, '.dat'], intersect(matNames, datNames), 'UniformOutput', false);
            fileList{ii} = vertcat(datList, matList);
        otherwise
            thisExt = dir(fullfile(folder, ['*' ext]));
            fileList{ii} = {thisExt.name};
    end
end

% Unpack filelist:
out = {};
for ii = 1:length(fileList)
    % Unpack list:
    list = [fileList{ii}];
    if size(list, 1) < size(list, 2)
        list = list';
    end
    out = [out; list];
end
end
