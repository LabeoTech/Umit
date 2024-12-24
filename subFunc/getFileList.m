function out = getFileList(folder, varargin)
% GETFILELIST retrieves a list of relevant files in a folder based on extension (.dat or .datstat).
%
% Inputs:
%   folder (char): The path to the folder containing the files.
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
addRequired(p,'folder',@isfolder)
addOptional(p,'extension','all',@(x) ismember(lower(x),{'.dat','.datstat','all'}))
parse(p,folder,varargin{:});
extension = lower(p.Results.extension);
% Get a list of .dat files
datFiles = dir(fullfile(folder, '*.dat'));
% Get a list of .datstat files
statFiles = dir(fullfile(folder,'*.datstat'));
% Filter list of mat files to get only those with the proper toolbox format

switch extension
    case '.dat'
        out = {datFiles.name}';
    case '.datstat'
        out = {statFiles.name}';
    case 'all'
        allFiles = vertcat(datFiles, statFiles);
        out = {allFiles.name}';
    otherwise
        error('File extension "%s" not supported by this function!',extension)
end


end