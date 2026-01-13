function nChunks = calculateMaxChunkSize(numElements, overheadFactor, precision)
% CALCULATEMAXCHUNKSIZE Estimate number of chunks to avoid out-of-memory.
%
%   nChunks = CALCULATEMAXCHUNKSIZE(numElements, overheadFactor)
%   nChunks = CALCULATEMAXCHUNKSIZE(numElements, overheadFactor, precision)
%
%   INPUTS:
%       numElements    - Total number of data elements to process (scalar)
%       overheadFactor - Multiplicative factor for processing memory
%       precision      - Optional: numeric type as string (default: 'single')
%                        Supported types: 'double','single','logical',
%                        'int8','uint8','int16','uint16','int32','uint32',
%                        'int64','uint64'
%
%   OUTPUT:
%       nChunks        - Recommended number of chunks to safely process data
%
%   NOTES:
%       - Memory usage is estimated as:
%         bytesPerElement * numElements * overheadFactor
%       - Uses available physical RAM to compute the safe number of chunks.
%       - If RAM detection fails, nChunks defaults to 1.

% --- Default precision ---
if nargin < 3
    precision = 'single';
end

% --- Determine bytes per element based on numeric type ---
switch lower(precision)
    case {'double','int64','uint64'}
        bytesPerElement = 8;
    case {'single','int32','uint32'}
        bytesPerElement = 4;
    case {'int16','uint16'}
        bytesPerElement = 2;
    case {'int8','uint8','logical'}
        bytesPerElement = 1;
    otherwise
        error('Unsupported numeric type: %s', precision);
end

% --- Estimate memory needed for processing ---
estimatedBytes = numElements * bytesPerElement * overheadFactor;

% --- Detect available RAM ---
try
    if ispc
        [~, sys] = memory;
        bytesAvailable = sys.PhysicalMemory.Available;
    elseif isunix && ~ismac
        [status, result] = system('free -b | awk ''/Mem:/ {print $7}''');
        if status
            error('Failed to retrieve RAM info on Linux');
        end
        bytesAvailable = str2double(result);
    elseif ismac
        [status, result] = system('vm_stat | awk ''/Pages free:/ {print $3}''');
        if status
            error('Failed to retrieve RAM info on macOS');
        end
        page_size = str2double(system('getconf PAGESIZE'));
        bytesAvailable = str2double(result) * page_size;
    else
        error('Unsupported platform');
    end

    % --- Compute number of chunks ---
    nChunks = max(1, ceil(estimatedBytes / bytesAvailable));

catch ME
    % --- Fallback if RAM detection fails ---
    ex = MException('calculateMaxChunkSize:NoOSAccess', ...
                    ['Failed to assess available RAM for OS: ', computer]);
    ex = addCause(ex, ME);
    warning(getReport(ex,'extended','hyperlinks','off'));
    nChunks = 1;
end

end
