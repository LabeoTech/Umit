function hdr = decodeDatHeader(bytes)
%DECODEDATHEADER Parse and validate the bytes of a self-describing .dat header.
%
%   hdr = decodeDatHeader(bytes)
%
%   Input:
%       bytes - uint8 vector holding the start of a .dat file. At least
%               the 512 header bytes are needed; extra bytes are ignored.
%
%   Output:
%       hdr - Scalar struct with fields:
%               headerVersion     - Header version (1)
%               writeComplete     - Logical; true once the writer finished
%               headerLength      - Header size in bytes (512)
%               dataClass         - MATLAB class of the data, e.g. 'single'
%               bytesPerValue     - Bytes per data value
%               frameRateHz       - Frame rate; NaN without a T axis
%               exposureMsec      - Exposure in ms; NaN when unknown
%               channelName       - Channel label (char)
%               dimNames          - Present axes in slot order, e.g.
%                                   {'Y','X','T'} or {'Y','X','E'}
%               dimSizes          - Sizes of the present axes (double row)
%               dataOffset        - Byte offset of the data (= headerLength)
%               expectedDataBytes - prod(dimSizes) * bytesPerValue
%
%   The byte-order mark is read before any other multi-byte field, so
%   headers written in either byte order decode correctly.
%
%   Errors:
%       Umitoolbox:decodeDatHeader:noHeader      - BYTES do not start with
%                                                  the header magic.
%       Umitoolbox:decodeDatHeader:corruptHeader - The magic matches but the
%           header is unusable: fewer than 512 bytes, invalid byte-order
%           mark, unknown version, header length other than 512, nonzero
%           reserved flag bits, a nonzero reserved axis slot, or a channel
%           name without a null terminator. Never falls back to legacy
%           reading when the magic matches.
%       Umitoolbox:validateDatHeader:*           - The decoded description
%                                                  breaks a validation rule.
%
%   This function does not warn when writeComplete is false; warning
%   about possibly incomplete files is the job of the code that opens a
%   file for data access.
%
%   See also: readDatHeader, encodeDatHeader, validateDatHeader, datHeaderSchema

[layout, codes] = datHeaderSchema(1);
constants = codes.constants;

if ~isa(bytes, 'uint8') || ~(isvector(bytes) || isempty(bytes))
    error('Umitoolbox:decodeDatHeader:invalidInput', ...
        'Header bytes must be a uint8 vector.');
end
bytes = bytes(:).';

nMagic = numel(constants.magic);
if numel(bytes) < nMagic || ~isequal(bytes(1:nMagic), constants.magic)
    error('Umitoolbox:decodeDatHeader:noHeader', ...
        'The bytes do not start with the .dat header magic.');
end

if numel(bytes) < constants.headerLength
    iCorrupt(sprintf('only %d of the %d header bytes are present.', ...
        numel(bytes), constants.headerLength));
end
bytes = bytes(1:constants.headerLength);

% --- byte order ---------------------------------------------------------
bomField = layout(strcmp({layout.name}, 'byteOrderMark'));
bomBytes = bytes(bomField.offset + (1:bomField.nBytes));
expectedLE = typecast(iToLittleEndian(constants.byteOrderMark), 'uint8');
if isequal(bomBytes, expectedLE)
    fileIsLittleEndian = true;
elseif isequal(bomBytes, fliplr(expectedLE))
    fileIsLittleEndian = false;
else
    iCorrupt(sprintf('invalid byte-order mark [%s].', num2str(double(bomBytes))));
end

readField = @(name) iGet(bytes, layout, name, fileIsLittleEndian);

% --- preamble -----------------------------------------------------------
headerVersion = double(readField('headerVersion'));
if headerVersion ~= constants.headerVersion
    iCorrupt(sprintf('unknown header version %d.', headerVersion));
end

headerLength = double(readField('headerLength'));
if headerLength ~= constants.headerLength
    iCorrupt(sprintf('header length is %d instead of %d.', headerLength, constants.headerLength));
end

flags = readField('flags');
reservedFlagBits = bitset(flags, constants.writeCompleteBit + 1, 0);
if reservedFlagBits ~= 0
    iCorrupt(sprintf('reserved flag bits are set (flags = 0x%04X).', double(flags)));
end
writeComplete = logical(bitget(flags, constants.writeCompleteBit + 1));

% --- dimension table ----------------------------------------------------
axisSizes = double(readField('axisSizes'));
reservedSlots = [codes.axis.reserved];
if any(axisSizes(reservedSlots) ~= 0)
    iCorrupt('a reserved axis slot (5 to 7) is nonzero.');
end
present = ~reservedSlots & axisSizes > 0;
dimNames = {codes.axis(present).name};
dimSizes = axisSizes(present);

% --- channel metadata ---------------------------------------------------
dtypeCode = double(readField('dtype'));
classIdx = find([codes.dataClass.code] == dtypeCode, 1);
if isempty(classIdx)
    dataClass = sprintf('unknown (code %d)', dtypeCode);
    bytesPerValue = NaN;
else
    dataClass = codes.dataClass(classIdx).className;
    bytesPerValue = codes.dataClass(classIdx).bytesPerValue;
end

nameBytes = readField('channelName');
nullIdx = find(nameBytes == 0, 1);
if isempty(nullIdx)
    iCorrupt('channelName has no null terminator.');
end
channelName = char(nameBytes(1:nullIdx - 1));
if isempty(channelName)
    channelName = '';
end

hdr = struct();
hdr.headerVersion = headerVersion;
hdr.writeComplete = writeComplete;
hdr.headerLength = headerLength;
hdr.dataClass = dataClass;
hdr.bytesPerValue = bytesPerValue;
hdr.frameRateHz = double(readField('frameRateHz'));
hdr.exposureMsec = double(readField('exposureMsec'));
hdr.channelName = channelName;
hdr.dimNames = dimNames;
hdr.dimSizes = dimSizes;
hdr.dataOffset = headerLength;
hdr.expectedDataBytes = prod(dimSizes) * bytesPerValue;

validateDatHeader(hdr, 'ReservedBytes', readField('reservedChannel'));

end

% =========================================================================
function iCorrupt(reason)
error('Umitoolbox:decodeDatHeader:corruptHeader', 'Corrupt .dat header: %s', reason);
end

function value = iGet(bytes, layout, fieldName, fileIsLittleEndian)
% Read FIELDNAME from BYTES as its LAYOUT class, honoring the file byte order.
f = layout(strcmp({layout.name}, fieldName));
value = typecast(bytes(f.offset + (1:f.nBytes)), f.class);
[~, ~, endian] = computer;
hostIsLittleEndian = endian == 'L';
if fileIsLittleEndian ~= hostIsLittleEndian
    value = swapbytes(value);
end
end

function value = iToLittleEndian(value)
[~, ~, endian] = computer;
if endian == 'B'
    value = swapbytes(value);
end
end
