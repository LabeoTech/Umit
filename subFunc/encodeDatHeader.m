function bytes = encodeDatHeader(hdr)
%ENCODEDATHEADER Build the 512-byte version-1 .dat header from a description.
%
%   bytes = encodeDatHeader(hdr)
%
%   Returns the header as a 1x512 uint8 row vector, little-endian. The
%   function performs no file I/O: data-level writers place these bytes at
%   the start of a .dat file together with the data they describe.
%
%   Input:
%       hdr - Scalar struct using the field names of decodeDatHeader.
%             Required fields:
%               dataClass    - MATLAB class of the data, e.g. 'single'
%               frameRateHz  - Frame rate (finite > 0 with a T axis,
%                              NaN without)
%               exposureMsec - Exposure in ms, NaN when unknown
%               channelName  - ASCII text; longer than 31 characters is
%                              truncated with a warning; '' is allowed
%               dimNames     - Present axes in slot order, e.g. {'Y','X','T'}
%               dimSizes     - Sizes of those axes
%             Optional field:
%               writeComplete - Logical, default false. Sets bit 0 of flags.
%
%   Output:
%       bytes - 1x512 uint8 header.
%
%   Each axis is encoded in its dedicated dimension-table slot. Slots of
%   absent axes, reserved slots, and all reserved bytes are zero. HDR is
%   validated with validateDatHeader before any byte is produced.
%
%   See also: decodeDatHeader, validateDatHeader, datHeaderSchema

[layout, codes] = datHeaderSchema(1);
constants = codes.constants;

% --- channel name (normalized before validation) -----------------------
if isstruct(hdr) && isscalar(hdr) && isfield(hdr, 'channelName')
    hdr.channelName = iNormalizeChannelName(hdr.channelName, constants.channelNameMaxChars);
end

% --- write-complete flag ------------------------------------------------
writeComplete = false;
if isstruct(hdr) && isscalar(hdr) && isfield(hdr, 'writeComplete') && ~isempty(hdr.writeComplete)
    if ~(islogical(hdr.writeComplete) || isnumeric(hdr.writeComplete)) || ...
            ~isscalar(hdr.writeComplete)
        error('Umitoolbox:validateDatHeader:invalidHeader', ...
            'Invalid .dat header: writeComplete must be a logical scalar.');
    end
    writeComplete = logical(hdr.writeComplete);
end

validateDatHeader(hdr);

% --- values -------------------------------------------------------------
dataClass = char(string(hdr.dataClass));
dtypeCode = codes.dataClass(strcmp({codes.dataClass.className}, dataClass)).code;

dimNames = cellstr(hdr.dimNames);
dimSizes = double(hdr.dimSizes(:).');
axisSizes = zeros(1, constants.nAxisSlots);
for k = 1:numel(dimNames)
    slot = codes.axis(strcmp({codes.axis.name}, dimNames{k})).slot;
    axisSizes(slot + 1) = dimSizes(k);
end

flags = uint16(0);
if writeComplete
    flags = bitset(flags, constants.writeCompleteBit + 1);
end

nameBytes = zeros(1, 32, 'uint8');
nameBytes(1:numel(hdr.channelName)) = uint8(hdr.channelName);

% --- assemble ------------------------------------------------------------
bytes = zeros(1, constants.headerLength, 'uint8');
bytes = iPut(bytes, layout, 'magic', constants.magic);
bytes = iPut(bytes, layout, 'byteOrderMark', constants.byteOrderMark);
bytes = iPut(bytes, layout, 'headerVersion', constants.headerVersion);
bytes = iPut(bytes, layout, 'flags', flags);
bytes = iPut(bytes, layout, 'headerLength', constants.headerLength);
bytes = iPut(bytes, layout, 'dtype', dtypeCode);
bytes = iPut(bytes, layout, 'frameRateHz', hdr.frameRateHz);
bytes = iPut(bytes, layout, 'exposureMsec', hdr.exposureMsec);
bytes = iPut(bytes, layout, 'channelName', nameBytes);
bytes = iPut(bytes, layout, 'axisSizes', axisSizes);
% reservedChannel stays zero.

end

% =========================================================================
function name = iNormalizeChannelName(name, maxChars)
if isstring(name) && isscalar(name)
    name = char(name);
end
if isempty(name) && (ischar(name) || isnumeric(name) || isstring(name))
    name = '';
    return
end
if ~ischar(name) || ~isrow(name)
    error('Umitoolbox:encodeDatHeader:invalidChannelName', ...
        'channelName must be a character row vector.');
end
codesInName = double(name);
if any(codesInName > 127) || any(codesInName == 0)
    error('Umitoolbox:encodeDatHeader:invalidChannelName', ...
        'channelName must contain only ASCII characters (no null characters): "%s".', name);
end
if numel(name) > maxChars
    warning('Umitoolbox:encodeDatHeader:channelNameTruncated', ...
        'channelName "%s" is longer than %d characters and was truncated to "%s".', ...
        name, maxChars, name(1:maxChars));
    name = name(1:maxChars);
end
end

function bytes = iPut(bytes, layout, fieldName, value)
% Write VALUE into the bytes of FIELDNAME as little-endian LAYOUT class.
f = layout(strcmp({layout.name}, fieldName));
values = cast(value(:).', f.class);
[~, ~, endian] = computer;
if endian == 'B'
    values = swapbytes(values);
end
raw = typecast(values, 'uint8');
assert(numel(raw) == f.nBytes, 'Umitoolbox:encodeDatHeader:internal', ...
    'Field %s encodes to %d bytes instead of %d.', fieldName, numel(raw), f.nBytes);
bytes(f.offset + (1:f.nBytes)) = raw;
end
