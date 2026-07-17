function checksum = computeFileChecksum(filePath, varargin)
%COMPUTEFILECHECKSUM Compute a hexadecimal checksum for a file.
%
%   checksum = computeFileChecksum(filePath)
%   checksum = computeFileChecksum(filePath, 'Algorithm', algorithm)
%
%   Inputs:
%       filePath  - Existing file to hash.
%
%   Name-Value options:
%       Algorithm - Java MessageDigest algorithm. Default: 'SHA-256'.
%
%   Output:
%       checksum  - Lowercase hexadecimal digest.

errID = 'Umitoolbox:computeFileChecksum:invalidInput';

p = inputParser;
p.FunctionName = 'computeFileChecksum';
addRequired(p, 'filePath', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Algorithm', 'SHA-256', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, filePath, varargin{:});

filePath = char(string(p.Results.filePath));
algorithm = char(string(p.Results.Algorithm));

if ~isfile(filePath)
    error(errID, 'File does not exist: %s', filePath);
end

try
    digest = java.security.MessageDigest.getInstance(algorithm);
catch ME
    error(errID, 'Unsupported checksum algorithm "%s": %s', ...
        algorithm, ME.message);
end

fid = fopen(filePath, 'rb');
if fid < 0
    error(errID, 'Could not open file for reading: %s', filePath);
end

cleanupObj = onCleanup(@() fclose(fid));
chunkSize = 1024 * 1024;

while true
    bytes = fread(fid, chunkSize, '*uint8');
    if isempty(bytes)
        break
    end
    digest.update(typecast(bytes(:), 'int8'));
end

digestBytes = typecast(digest.digest(), 'uint8');
checksum = lower(reshape(dec2hex(digestBytes, 2).', 1, []));

clear cleanupObj

end
