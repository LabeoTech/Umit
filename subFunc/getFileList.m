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
datList = {};
matFileList = dir(fullfile(folder, '*.mat'));
% Check if "dataHistory" exists inside each file
matFilesMap = arrayfun(@(x) matfile(fullfile(folder,x.name)), matFileList, 'UniformOutput',false);
idxValid = cellfun(@(x) isprop(x, 'dataHistory'), matFilesMap);
matList = {matFileList(idxValid).name}';

if ~strcmpi(extension,'.mat')
    datFileList = dir(fullfile(folder, '*.dat'));
    [~,datNames,~] = arrayfun(@(x) fileparts(x.name), datFileList, 'UniformOutput', false);
    [~,matNames,~] = cellfun(@(x) fileparts(x), matList, 'UniformOutput', false);
    datList = {datFileList(ismember(datNames,matNames)).name}';
end
switch extension
    case '.mat'
        fileList = matList;
    case '.dat'
        fileList = datList;
    otherwise
        fileList = vertcat(datList,matList);
end

end


