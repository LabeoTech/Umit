function [layout, codes] = datHeaderSchema(version)
%DATHEADERSCHEMA Return the byte layout and code tables of a .dat header.
%
%   layout = datHeaderSchema(version)
%   [layout, codes] = datHeaderSchema(version)
%
%   Single source of truth for the self-describing .dat header defined in
%   docs/dev/dat-header-spec.md. Every other header function reads byte
%   offsets, classes, codes, and fixed values from here.
%
%   Input:
%       version - Header version (numeric scalar). Only version 1 exists.
%
%   Outputs:
%       layout  - Struct array, one element per header field, with fields:
%                   name    - Field name (char)
%                   offset  - Zero-based byte offset in the header
%                   class   - MATLAB class of the stored values
%                   nBytes  - Number of bytes the field occupies
%                 The fields are contiguous and cover the whole header.
%       codes   - Structure with:
%                   codes.dataClass - Struct array with fields code,
%                                     className, bytesPerValue: the data
%                                     class codes stored in "dtype".
%                   codes.axis      - Struct array with fields slot,
%                                     name, reserved: the axis assigned to
%                                     each dimension-table slot (slot is
%                                     zero-based; reserved slots have an
%                                     empty name).
%                   codes.constants - Fixed values of the version: magic,
%                                     byteOrderMark, headerVersion,
%                                     headerLength, writeCompleteBit,
%                                     channelNameMaxChars, nAxisSlots.
%
%   Notes:
%       - Multi-byte fields are written little-endian; readers use the
%         byte-order mark to accept either order.
%       - The dimension table ("axisSizes") holds one uint32 size per
%         slot. Size 0 means the slot's axis is absent.
%
%   See also: encodeDatHeader, decodeDatHeader, readDatHeader,
%             validateDatHeader, isDatWithHeader

if ~isnumeric(version) || ~isscalar(version)
    error('Umitoolbox:datHeaderSchema:unknownHeaderVersion', ...
        'Header version must be a numeric scalar.');
end

switch version
    case 1
        layout = struct( ...
            'name',   {'magic', 'byteOrderMark', 'headerVersion', 'flags', ...
                       'headerLength', 'dtype', 'frameRateHz', 'exposureMsec', ...
                       'channelName', 'reservedChannel', 'axisSizes'}, ...
            'offset', {0, 6, 8, 10, 12, 16, 17, 21, 25, 57, 480}, ...
            'class',  {'uint8', 'uint16', 'uint16', 'uint16', 'uint32', ...
                       'uint8', 'single', 'single', 'uint8', 'uint8', 'uint32'}, ...
            'nBytes', {6, 2, 2, 2, 4, 1, 4, 4, 32, 423, 32});

        codes = struct();
        codes.dataClass = struct( ...
            'code',          {1, 2, 3, 4, 5, 6, 7}, ...
            'className',     {'single', 'double', 'uint8', 'uint16', 'int16', 'uint32', 'int32'}, ...
            'bytesPerValue', {4, 8, 1, 2, 2, 4, 4});
        codes.axis = struct( ...
            'slot',     {0, 1, 2, 3, 4, 5, 6, 7}, ...
            'name',     {'Y', 'X', 'T', 'E', 'F', '', '', ''}, ...
            'reserved', {false, false, false, false, false, true, true, true});

        constants = struct();
        constants.magic = uint8([137 85 77 68 13 10]); % 0x89 'U' 'M' 'D' CR LF
        constants.byteOrderMark = uint16(258);         % 0x0102
        constants.headerVersion = 1;
        constants.headerLength = 512;
        constants.writeCompleteBit = 0;                % bit 0 of flags
        constants.channelNameMaxChars = 31;            % 32 bytes incl. null
        constants.nAxisSlots = 8;
        codes.constants = constants;

    otherwise
        error('Umitoolbox:datHeaderSchema:unknownHeaderVersion', ...
            'Unknown .dat header version: %g. Supported version: 1.', version);
end

end
