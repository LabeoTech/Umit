function out = getFileList(folder, extension)
% GETFILELIST retrieves a list of relevant files in a folder based on extension (.dat or .mat).
%
% Inputs:
%   folder (char): The path to the folder containing the files.
%   extension (char): The desired file extension or 'all' to list all files.
%
% Output:
%   out (cell): List of valid filenames with the specified extension. Empty
%   cell array if no files are found.
%
% Example usage:
%   fileList = getFileList('/path/to/folder', '.mat');
%
% Note: This function does not change the algorithm for file retrieval.



% Get a list of .dat files
datFiles = dir(fullfile(folder, '*.dat'));
% Get a list of .mat files
matFiles = dir(fullfile(folder,'*.mat'));
% Filter list of mat files to get only those with the proper toolbox format

% Here, we check if the file contains the variable "datSize" to check if it
% is valid or not
warning('off'); % Disable warning to suppress messages from .mat files that do not have "datSize" variable.
validMatFiles = false(size(matFiles));
for ii = 1:length(matFiles)
    data = load(fullfile(folder, matFiles(ii).name), 'datSize');
    validMatFiles(ii) = ~isempty(fieldnames(data));
end
matFiles = matFiles(validMatFiles);
warning('on');

switch lower(erase(extension,'.'))
    case 'dat'
        out = datFiles;
    case 'mat'
        out = matFiles;
    case 'all'
        out = [datFiles;matFiles];
    otherwise
        error('Unknown file extension "%s"!',extension)
end

if isempty(out)
    out = cell.empty(0,1);
else
    out = {out.name};
end

end