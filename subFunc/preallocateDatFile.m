function preallocateDatFile(filename, dims, precision)
%PREALLOCATEDATFILE Safely preallocate a raw binary file on disk.
%
%   preallocateDatFile(filename, dims)
%   preallocateDatFile(filename, dims, precision)
%
%   This function writes zeros to disk to guarantee correct allocation and
%   later random-access writes. The chunk size is dynamically computed to
%   respect available RAM using calculateMaxChunkSize.
%
%   Inputs:
%       filename  - Output raw binary file path.
%       dims      - Numeric vector of array dimensions stored on disk.
%       precision - Numeric class stored in the file.
%                   Default: 'single'
%
%   Notes:
%       - This function only preallocates the binary file.
%       - It no longer writes any metadata sidecar file.
%       - Raw .dat files in this project are typically stored as YXT single.

if nargin < 3 || isempty(precision)
    precision = 'single';
end

validateattributes(filename, {'char','string'}, {'nonempty'}, ...
    mfilename, 'filename');
filename = char(string(filename));

validateattributes(dims, {'numeric'}, ...
    {'vector','nonempty','real','finite','positive','integer'}, ...
    mfilename, 'dims');

if ~(ischar(precision) || (isstring(precision) && isscalar(precision)))
    error('preallocateDatFile:InvalidPrecision', ...
        'precision must be a character vector or string scalar.');
end
precision = char(string(precision));

bytesPerElement = getByteSize(precision);
nElements = prod(double(dims));
totalBytes = nElements * bytesPerElement;

nChunks = calculateMaxChunkSize(totalBytes, 1, 0.2);
chunkElements = ceil(nElements / nChunks);

fid = fopen(filename, 'w');
cid = onCleanup(@() safeFclose(fid)); %#ok<NASGU>
if fid < 0
    error('preallocateDatFile:OpenFailed', ...
        'Cannot open file for writing: %s', filename);
end

disp('Preallocating .dat file...')

for k = 1:nChunks
    startIdx = (k-1) * chunkElements + 1;
    endIdx   = min(k * chunkElements, nElements);
    currentChunkSize = endIdx - startIdx + 1;

    fwrite(fid, zeros(currentChunkSize, 1, precision), precision);
end

fclose(fid);
disp('Empty .dat file created.')
end