function out = getFileList(folder, varargin)
% GETFILELIST retrieves a list of relevant files in a folder based on extension (.dat or .datstat).
%
% Inputs:
%   folder (char|cell): The path to the folder containing the files. If a
%   cell array of folders is provided, the list of unique file names is
%   retrieved.
%   Optional:
%   extension (char|default = 'all'): The desired file extension. If not
%   set, the function returns all .dat and .datstat files.
%
% Output:
%   out (cell): List of valid filenames with the specified extension. Empty
%   cell array if no files are found.
%
% Example usage:
%   fileList = getFileList('/path/to/folder', '.dat'); % Gets all ".dat"
%   files from the given folder.
%
p = inputParser;
addRequired(p,'folder',@(x) all(isfolder(x)))
addOptional(p,'extension','all',@(x) ismember(lower(x),{'.dat','.datstat','all'}))
parse(p,folder,varargin{:});
extension = lower(p.Results.extension);

if ischar(folder);folder = {folder};end
out = {};

for ii = 1:length(folder)
    out = [out;getList(folder{ii})];%#ok
end
% Get unique file names:
out = unique(out, 'stable');
    function fileList = getList(folder)
        % Get a list of .dat files
        datFiles = dir(fullfile(folder, '*.dat'));
        % Get a list of .datstat files
        statFiles = dir(fullfile(folder,'*.datstat'));
        % Filter list of mat files to get only those with the proper toolbox format
        
        switch extension
            case '.dat'
                fileList = {datFiles.name}';
            case '.datstat'
                fileList = {statFiles.name}';
            case 'all'
                allFiles = vertcat(datFiles, statFiles);
                fileList = {allFiles.name}';
            otherwise
                error('File extension "%s" not supported by this function!',extension)
        end
    end
end