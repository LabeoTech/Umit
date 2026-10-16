function [hash, image] = computeImageContentHash(image)
%COMPUTEIMAGECONTENTHASH Content hash of a 2-D reference image.
%
%   hash = computeImageContentHash(image)
%   [hash, image] = computeImageContentHash(image)
%
%   Identifies a reference image by its content, not by any file: the same
%   image gives the same hash whatever its source (project library, .dat
%   extraction, exported file) and however often it is re-saved.
%
%   Canonical form (normative; also in IMPLEMENTATION_NOTES_project.md):
%       1) Normalize exactly as genImageReferenceStruct does: convert to
%          single and set non-finite values to 0. Negative zero is also
%          written as +0, so the bytes do not depend on the sign of zero.
%       2) Hash input, in this order:
%            - the class name 'single' as UTF-8 bytes, then one 0x00 byte;
%            - the number of dimensions as uint32, little-endian;
%            - each dimension size as uint64, little-endian (size(image));
%            - the image values as single, little-endian, column-major.
%       3) hash = 'sha256:' + lowercase hex SHA-256 of that byte sequence
%          (UMITLedgerCore.hashBytes).
%
%   Input:
%       image - Non-empty numeric or logical 2-D matrix.
%
%   Outputs:
%       hash  - 'sha256:<64 hex>'.
%       image - The normalized single image that was hashed.
%
%   See also genImageReferenceStruct, validateImageReferenceStruct,
%   UMITLedgerCore.

if ~(isnumeric(image) || islogical(image)) || isempty(image) || ~ismatrix(image)
    error('Umitoolbox:computeImageContentHash:invalidInput', ...
        'The image must be a non-empty numeric or logical 2-D matrix.');
end

image = single(image);
image(~isfinite(image)) = 0;
image(image == 0) = 0;              % -0 -> +0

sizeVector = size(image);
values = image(:);
if iIsBigEndian()
    values = swapbytes(values);
    dimCount = swapbytes(uint32(numel(sizeVector)));
    dims = swapbytes(uint64(sizeVector(:)));
else
    dimCount = uint32(numel(sizeVector));
    dims = uint64(sizeVector(:));
end

bytes = [uint8('single'), uint8(0), ...
    typecast(dimCount, 'uint8'), ...
    typecast(dims.', 'uint8'), ...
    typecast(values.', 'uint8')];
hash = UMITLedgerCore.hashBytes(bytes);
end

function tf = iIsBigEndian()
[~, ~, endian] = computer();
tf = strcmpi(endian, 'B');
end
