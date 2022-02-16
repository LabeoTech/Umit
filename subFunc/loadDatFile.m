function [outFile, metaData] = loadDatFile(DatFileName, varargin)
%LOADDATFILE opens binary data stored as .DAT files and their metadata. 
% Inputs:
%   datFileName: full path of .DAT file
%   metaDatFileName(Optional) : fullpath of .MAT file containing the metadata of
%       "datFile". Enter an empty char('') to automatically look for the .mat
%       metaData file.
% Outputs:
%   outFile (numerical array) : data matrix containing the mapped
%       data.
%   metaData (struct): metaData associated with "outFile".

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @mustBeDatFile);
addOptional(p, 'metaDatFileName', '', @(x) ischar(x) || isstring(x));
parse(p, DatFileName, varargin{:});
%%%%%%
%Initialize variables:
DatFileName = p.Results.DatFileName;
metaDatFileName = p.Results.metaDatFileName;
%%%%
% Get metaData:
if isempty(metaDatFileName)
    metaDatFileName = strrep(DatFileName, '.dat', '.mat');
    if ~isfile(metaDatFileName)% For retro-compatibility.
        metaDatFileName = strrep(DatFileName, '.dat', '_info.mat'); 
    end
end
% Get metaData file:
metaData = load(metaDatFileName);

% Read binary file:
fid = fopen(DatFileName);
outFile = fread(fid, inf, '*single');
outFile = reshape(outFile, [metaData.datSize metaData.datLength]);
end

% Validation function
function mustBeDatFile(datFileName)
if ~isfile(datFileName) && ~strcmp(datFileName(end-3:end), '.dat')
    errID = 'umIToolbox:loadDatFile:InvalidInput';
    msg = [datFileName ' is not a .DAT file!'];
    throwAsCaller(MException(errID,msg))
end
end