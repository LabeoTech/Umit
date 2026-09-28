function tf = isDatWithHeader(filename)
%ISDATWITHHEADER True if a .dat file starts with the self-describing header magic.
%
%   tf = isDatWithHeader(filename)
%
%   Reads only the first 6 bytes of the file and compares them with the
%   header magic (0x89 'U' 'M' 'D' CR LF). Returns false for legacy
%   headerless files and for files shorter than 6 bytes. A true result
%   only means the file claims to have a header; use readDatHeader to
%   decode and validate it.
%
%   Input:
%       filename - Path of an existing file (char or string).
%
%   Output:
%       tf       - Logical scalar.
%
%   See also: readDatHeader, datHeaderSchema

filename = char(string(filename));
if ~isfile(filename)
    error('Umitoolbox:isDatWithHeader:fileNotFound', ...
        'File not found: %s', filename);
end

[~, codes] = datHeaderSchema(1);
magic = codes.constants.magic;

fid = fopen(filename, 'r');
if fid < 0
    error('Umitoolbox:isDatWithHeader:openFailed', ...
        'Cannot open file for reading: %s', filename);
end
cleanupObj = onCleanup(@() fclose(fid));

firstBytes = fread(fid, numel(magic), '*uint8');
tf = numel(firstBytes) == numel(magic) && isequal(firstBytes(:).', magic);

clear cleanupObj
end
