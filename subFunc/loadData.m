function [outFile, Info] = loadData(DatFileName)
%LOADDATA opens binary data stored as .DAT files and their metadata OR
% .mat files stored as ".datstat".
% Inputs:
%   datFileName: full path of .DAT file
% Output:
%   outFile (numerical array) : data matrix containing the mapped
%       data.
%   Info (struct): file meta data. For the new data format, this is
%       the same as the structure from the "AcqInfos.mat" file. For the older
%       data format, it corresponds to the content of the ".mat" file
%       associated with the .dat file. If .datstat, Info will be empty.

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) isfile(x) & (endsWith(x,'.dat') | endsWith(x,'.datstat')));
parse(p, DatFileName);
%%%%%%
%Initialize variables:
DatFileName = p.Results.DatFileName;
clear p
%%%%
Info = []; 
if isempty(fileparts(DatFileName))
    % Force full path:
    DatFileName = fullfile(pwd,DatFileName);
end

fprintf('Opening file "%s" ...\n',DatFileName);
% Load .dat files:
if endsWith(DatFileName,'.dat')
    [outFile,Info] = loadDat(DatFileName);
else
    % Loads the .mat structure in the "datstat" file.
    outFile = load(DatFileName, '-mat');            
end
disp('Done.');
end

% Local functions
function [data, AcqInfoStream] = loadDat(filename)
% Opens .dat files. Manages versions and retrocompatibility.
% Load AcqInfos:
load(fullfile(fileparts(filename),'AcqInfos.mat'));%#ok
datLen = [];
% For retrocompatibility:
fName = {strrep(filename,'.dat','.mat'), strrep(filename, '.dat', '_info.mat')};
fName = fName(cellfun(@isfile,fName));
if ~isempty(fName)
    % Load using info from meta data ".mat" file:
    matInfo = load(fName{1});
    datLen = matInfo.datLength;
    % Update information from "AcqInfos"
    AcqInfoStream.Heigth = matInfo.datSize(1);
    AcqInfoStream.Width = matInfo.datSize(2);
    AcqInfoStream.FrameRateHz = matInfo.Freq;
end
fid = fopen(filename);
data = fread(fid, inf, '*single');
data = reshape(data, AcqInfoStream.Height, AcqInfoStream.Width, datLen);
fclose(fid);
end