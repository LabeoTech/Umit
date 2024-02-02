function [mmFile, mmFileInfo]= mapDat(DatFileName)
%MAPDAT creates a memmory map file of binary data stored as .DAT
%files.
% DATFILENAME: full path of .DAT file

% MMFILE: memmapfile containing the mapped data 
% MMFILEINFO (optional): structure with file meta data. For the new data
% format, this is the same as the structure from the "AcqInfos.mat" file.
% For the older data format, it corresponds to the content of the ".mat"
% file associated with the .dat file.

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) isfile(x) & endsWith(x,'.dat'));
parse(p, DatFileName);
%Initialize variables:
DatFileName = p.Results.DatFileName;
clear p
%%%%
% Get file format:
[format, mmFileInfo]= getFileFormat(DatFileName);
% Create memmapfile:
mmFile = memmapfile(DatFileName, 'Format', format);
end

% Local functions
function [format, Info] = getFileFormat(DatFileName)
% Gathers information about the data file size to create the "Format"
% property for the memmap file:
load(fullfile(fileparts(DatFileName),'AcqInfos.mat'));%#ok
Info = AcqInfoStream;
% For retrocompatibility:
fName = {strrep(DatFileName,'.dat','.mat'), strrep(DatFileName, '.dat', '_info.mat')};
fName = fName(cellfun(@isfile,fName));
if isempty(fName)    
    ds = [Info.Height, Info.Width];
    % Calculate data length:
    fileInfo = dir(DatFileName);   
    dl = fileInfo.bytes/prod([ds,4]);    
else
    % Load using info from meta data ".mat" file:
    matInfo = load(fName{1});
    ds = matInfo.datSize;
    dl = matInfo.datLength;
    % Update Info:
    Info.Height = matInfo.datSize(1);
    Info.Width = matInfo.datSize(2);
    Info.Length = matInfo.datLength(1);
    Info.FrameRateHz = matInfo.Freq;
end

dt = 'single'; % Data type as "single"
dn = 'data'; % Data name
% Create "Format" cell array:
format = {dt double([ds dl]) dn};

end

