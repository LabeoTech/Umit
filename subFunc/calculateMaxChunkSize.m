function nChunks = calculateMaxChunkSize(data, sizeFactor, RAMoverhead)
% CALCULATEMAXCHUNKSIZE Estimate the number of chunks required to safely
% process data without exceeding system RAM.
%
% This function computes a conservative chunk count based on an estimate
% of peak memory usage, available system RAM, and a user-defined RAM
% overhead. The input data can be provided either as a numeric array or
% directly as a scalar specifying the required memory footprint in bytes.
%
% Inputs:
%   data :
%       Either:
%         - Numeric array to be processed. Its in-memory size is determined
%           using WHOS and multiplied by sizeFactor to estimate peak usage.
%         - Scalar specifying the estimated required memory in BYTES
%           (already including any sizeFactor or algorithm-specific scaling).
%
%   sizeFactor (scalar > 0):
%       Multiplicative factor used to estimate peak memory usage relative
%       to the input data size. This factor is only applied when DATA is
%       provided as a numeric array. It is ignored when DATA is a scalar.
%
%   RAMoverhead (scalar in [0,1)):
%       Fraction of total physical RAM to reserve for the operating system,
%       MATLAB internals, memory fragmentation, and temporary allocations.
%
% Output:
%   nChunks (scalar integer):
%       Number of chunks required to process the data while staying within
%       the computed RAM budget. Guaranteed to be at least 1.
%
% Notes:
%   - The function uses a conservative RAM budget based on available RAM
%     clamped to a fixed fraction of total system RAM.
%   - This approach avoids instability caused by transient memory peaks
%     and allocator cleanup during execution.


if ~exist('RAMoverhead','var')
    RAMoverhead = .2;
end


% -------------------------------------------------------------------------
% Data memory footprint
% -------------------------------------------------------------------------
if ~isscalar(data)
    s = whos('data');
    requiredBytes = sizeFactor * s.bytes;
else
    requiredBytes = data * sizeFactor;
end

% -------------------------------------------------------------------------
% Query system memory
% -------------------------------------------------------------------------
try
    switch computer
        case 'PCWIN64'
            [~, sys] = memory;
            totalRAM     = sys.PhysicalMemory.Total;
            availableRAM = sys.PhysicalMemory.Available;

        case 'GLNXA64'
            [~, t] = system("free -b | grep Mem | awk '{print $2}'");
            [~, a] = system("free -b | grep Mem | awk '{print $7}'");
            totalRAM     = str2double(t);
            availableRAM = str2double(a);

        case 'MACI64'
            [~, t] = system('sysctl -n hw.memsize');
            totalRAM = str2double(t);

            [~, p] = system("vm_stat | grep 'Pages free' | awk '{print $3}'");
            pageSize = str2double(system('getconf PAGESIZE'));
            availableRAM = str2double(p) * pageSize;

        otherwise
            error('Unsupported platform');
    end
catch ME
    warning('calculateMaxChunkSize:NoOSAccess', ...
            'RAM query failed, falling back to 1 chunk.\n%s', ME.message);
    nChunks = 1;
    return
end

% -------------------------------------------------------------------------
% Conservative RAM budgeting
% -------------------------------------------------------------------------
reservedBytes = RAMoverhead * totalRAM;

% Do not allow MATLAB to assume more than a safe fraction of total RAM
maxUsableFraction = 0.7;   % tunable, but stable
usableRAM = min(availableRAM, totalRAM * maxUsableFraction) - reservedBytes;

if usableRAM <= 0   
    nChunks = [];
    warning('calculateMaxChunkSize:NoUsableRAM', ...
            'No usable RAM after overhead reservation.');   
        return
end

% -------------------------------------------------------------------------
% Chunk estimation
% -------------------------------------------------------------------------
nChunks = max(1, ceil(requiredBytes / usableRAM));

end
