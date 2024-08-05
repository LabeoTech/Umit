function [outFile, Info] = loadDat(DatFileName)
%LOADDAT opens binary data stored as .DAT files and their metadata. 
% Inputs:
%   datFileName: full path of .DAT file
% Output:
%   outFile (numerical array) : data matrix containing the mapped
%       data.
%   Info (struct): file meta data. For the new data% format, this is 
%       the same as the structure from the "AcqInfos.mat" file. For the older 
%       data format, it corresponds to the content of the ".mat" file 
%       associated with the .dat file.


% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) isfile(x) & endsWith(x,'.dat'));
parse(p, DatFileName);
%%%%%%
%Initialize variables:
DatFileName = p.Results.DatFileName;
clear p
%%%%
fprintf('Opening file "%s" ...\n',DatFileName);
% Load AcqInfos:
load(fullfile(fileparts(DatFileName),'AcqInfos.mat'));%#ok
Info = AcqInfoStream;
datLen = [];
% For retrocompatibility:
fName = {strrep(DatFileName,'.dat','.mat'), strrep(DatFileName, '.dat', '_info.mat')};
fName = fName(cellfun(@isfile,fName));
if ~isempty(fName)
    % Load using info from meta data ".mat" file:
    matInfo = load(fName{1});
    datLen = matInfo.datLength;
    % Update information from "AcqInfos"
    Info.Heigth = matInfo.datSize(1);
    Info.Width = matInfo.datSize(2);
    Info.FrameRateHz = matInfo.Freq;
end
fid = fopen(DatFileName);
outFile = fread(fid, inf, '*single');
outFile = reshape(outFile, Info.Height, Info.Width, datLen);
fclose(fid);
disp('Done.');
end