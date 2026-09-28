function hdr = readDatHeader(filename)
%READDATHEADER Read and decode the self-describing header of a .dat file.
%
%   hdr = readDatHeader(filename)
%
%   Reads the first 512 bytes of FILENAME (fewer if the file is shorter)
%   and returns decodeDatHeader(bytes). See decodeDatHeader for the fields
%   of HDR and the errors raised for corrupt or invalid headers.
%
%   Input:
%       filename - Path of an existing .dat file (char or string).
%
%   Output:
%       hdr      - Decoded and validated header description.
%
%   Errors:
%       Umitoolbox:readDatHeader:fileNotFound - FILENAME does not exist.
%       Umitoolbox:readDatHeader:noHeader     - The file does not start with
%           the header magic (legacy headerless or unrelated file). Choosing
%           between header and legacy reading belongs to the file-opening
%           entry point, which uses isDatWithHeader.
%
%   This function does not warn when hdr.writeComplete is false.
%
%   See also: decodeDatHeader, isDatWithHeader, datHeaderSchema

filename = char(string(filename));
if ~isfile(filename)
    error('Umitoolbox:readDatHeader:fileNotFound', 'File not found: %s', filename);
end

[~, codes] = datHeaderSchema(1);
constants = codes.constants;

fid = fopen(filename, 'r');
if fid < 0
    error('Umitoolbox:readDatHeader:openFailed', ...
        'Cannot open file for reading: %s', filename);
end
cleanupObj = onCleanup(@() fclose(fid));
bytes = fread(fid, constants.headerLength, '*uint8');
clear cleanupObj

bytes = bytes(:).';
nMagic = numel(constants.magic);
if numel(bytes) < nMagic || ~isequal(bytes(1:nMagic), constants.magic)
    error('Umitoolbox:readDatHeader:noHeader', ...
        'File does not start with a .dat header (legacy or unrelated file): %s', filename);
end

hdr = decodeDatHeader(bytes);
end
