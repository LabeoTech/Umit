function mmFile = mapDat(DatFileName)
%MAPDAT creates a memmory map file of binary data stored as .DAT
%files.
% DATFILENAME: full path of .DAT file

% MMFILE: memmapfile containing the mapped data 

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) isfile(x) & endsWith(x,'.dat'));
parse(p, DatFileName);
%Initialize variables:
DatFileName = p.Results.DatFileName;
clear p
%%%%
% Get file format:
format = getFileFormat(DatFileName);
% Create memmapfile:
mmFile = memmapfile(DatFileName, 'Format', format);
end

% Local functions
function format = getFileFormat(DatFileName)
% Gathers information about the data file size to create the "Format"
% property for the memmap file:

% For retrocompatibility:
fName = {strrep(DatFileName,'.dat','.mat'), strrep(DatFileName, '.dat', '_info.mat')};
fName = fName(cellfun(@isfile,fName));
if isempty(fName)
    % Load info from "AcqInfos.mat" file:
    Info = loadAcqInfo(fileparts(DatFileName));    
    ds = [Info.AcqInfoStream.Width, Info.AcqInfoStream.Height];
    % Calculate data length:
    fileInfo = dir(DatFileName);   
    dl = fileInfo.bytes/prod([ds,4]);
else
    % Load using info from meta data ".mat" file:
    Info = load(fName{1});
    ds = Info.datSize;
    dl = Info.datLength;
end

dt = 'single'; % Data type as "single"
dn = 'data'; % Data name
% Create "Format" cell array:
format = {dt double([ds dl]) dn};

end

