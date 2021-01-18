function mmFile = mapDatFile(datFileName, metaDatFileName)
%MAPDATFILE creates a memmory map file of binary data stored as .DAT
%files.
% DATFILENAME: full path of .DAT file
% METADATFILENAME : fullpath of .MAT file containing the metadata of
% DATFILE.
% MMFILE: memmapfile containing the mapped data with dimensions specified
%  as [METADATFILE.DATSIZE  METADATFILE.DATLENGTH]. The dataset found in
%  MMFILE.DATA.data.

arguments
    datFileName {mustBeDatFile(datFileName)}
    metaDatFileName {mustBeMATfileOrEmpty(metaDatFileName)} = ''
end
if isempty(metaDatFileName)
    metaDatFileName = strrep(datFileName, '.dat', '_info.mat');
end
format = getFileFormat(metaDatFileName);
mmFile = memmapfile(datFileName, 'Format', format);
end

% Local functions
function format = getFileFormat(metaDatFileName)
md = matfile(metaDatFileName);
dt = md.Datatype;
ds = md.datSize;
dl = md.datLength;
dn = md.datName;
if iscell(dt)
    format = cell(numel(dt), 3);
    for i = 1:length(md.Datatype)
        format(i,:) = {dt{i}, double([ds{i} dl{i}]) dn{i}};
    end
else
    format = {dt double([ds dl]) dn};
end
end

% Validation function
function mustBeDatFile(datFileName)
if ~isfile(datFileName) && ~strcmp(datFileName(end-3:end), '.dat')
    errID = 'IsaToolbox:InvalidFile';
    msg = [datFileName ' is not a .DAT file!'];
    throwAsCaller(MException(errID,msg))
end
end

function mustBeMATfileOrEmpty(metaDatFile)
if ~isfile(metaDatFile) && ~isempty(metaDatFile) || isfile(metaDatFile) && ~strcmp(metaDatFile(end-3:end), '.mat')
    errID = 'IsaToolbox:InvalidFile';
    msg = 'Invalid File!';
    throwAsCaller(MException(errID,msg))
end
end
