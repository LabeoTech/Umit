function fileList = getFileList(folder, extension)
% GETFILELIST creates a list of .mat/.dat files created using umIT.
% The files considered here are the ones with the .mat files containing the
% variable "dataHistory". This is the criterion used to identify the
% pertinent files.
% Inputs:
%   folder (char): path to the folder containing the files:
%   extension (char): extension of the files to list: {'.mat','.dat','all'}.
%       Set to 'all', to retrive a list of .mat AND .dat files.
% Output:
%   fileList (cell): list of valid filenames with file extension. Empty, if
%       no valid files are found.

matFileList = dir(fullfile(folder, '*.mat'));
% Check if "dataHistory" exists inside each file
idxValid = false(size(matFileList));
warning('off')
for ii = 1:length(matFileList)
    a = load(fullfile(folder, matFileList(ii).name), 'dataHistory');
    idxValid(ii) = ~isempty(fieldnames(a));
end
clear a
warning('on')
if ~any(idxValid)
    fileList = {};
    return
end
% matFilesMap = arrayfun(@(x) matfile(fullfile(folder,x.name)), matFileList, 'UniformOutput',false);
% idxValid = cellfun(@(x) isprop(x, 'dataHistory'), matFilesMap);
matList = {matFileList(idxValid).name}';
[~,matNames,~] = cellfun(@(x) fileparts(x), matList, 'UniformOutput', false);

datFileList = dir(fullfile(folder, '*.dat'));
[~,datNames,~] = arrayfun(@(x) fileparts(x.name), datFileList, 'UniformOutput', false);

switch extension
    case '.mat'
        fileList = cellfun(@(x) [x, '.mat'],setdiff(matNames,datNames), 'UniformOutput',false);
    case '.dat'
        fileList = cellfun(@(x) [x, '.dat'],intersect(matNames,datNames), 'UniformOutput',false);
    otherwise
        matList = cellfun(@(x) [x, '.mat'],setdiff(matNames,datNames), 'UniformOutput',false);
        datList = cellfun(@(x) [x, '.dat'],intersect(matNames,datNames), 'UniformOutput',false);
        fileList = vertcat(datList, matList);
end
end