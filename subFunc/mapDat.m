function [mmFile, mmFileInfo]= mapDat(DatFileName)
%MAPDAT creates a memmory map file of binary data stored as .DAT
%files.
% DATFILENAME: full path of .DAT file

% MMFILE: memmapfile containing the mapped data in field "data", with the
% shape Info.dimSizes and the class Info.dataClass.
% MMFILEINFO (optional): .dat metadata from loadMetaData. The data offset,
% class, and shape of every readable .dat kind (headered, legacy sidecar,
% AcqInfos-bound) come from there.

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) isfile(x) & endsWith(x,'.dat'));
parse(p, DatFileName);
%Initialize variables:
DatFileName = p.Results.DatFileName;
clear p
%%%%
mmFileInfo = loadMetaData(DatFileName);
% Create memmapfile:
mmFile = memmapfile(DatFileName, 'Offset', mmFileInfo.dataOffset, ...
    'Format', {mmFileInfo.dataClass, double(mmFileInfo.dimSizes), 'data'});
end
