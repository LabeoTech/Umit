function outFile = loadDat(DatFileName)
%LOADDAT opens binary data stored as .DAT files and their metadata. 
% Inputs:
%   datFileName: full path of .DAT file
% Output:
%   outFile (numerical array) : data matrix containing the mapped
%       data.

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
% For retrocompatibility:
fName = {strrep(DatFileName,'.dat','.mat'), strrep(DatFileName, '.dat', '_info.mat')};
fName = fName(cellfun(@isfile,fName));
if isempty(fName)
    % Load info from "AcqInfos.mat" file:
    Info = loadAcqInfo(fileparts(DatFileName));    
    datSize = [Info.AcqInfoStream.Width, Info.AcqInfoStream.Height];
    datLen = [];
else
    % Load using info from meta data ".mat" file:
    Info = load(fName{1});
    datSize = Info.datSize;
    datLen = Info.datLength;    
end
fid = fopen(DatFileName);
outFile = fread(fid, inf, '*single');
outFile = reshape(outFile, datSize(1), datSize(2),datLen);
fclose(fid);
disp('Done.');
end