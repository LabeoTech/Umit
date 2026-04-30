function nChunks = calculateMaxChunkSize(data, sizeFactor, RAMoverhead)
%CALCULATEMAXCHUNKSIZE Estimate the number of chunks required for RAM-safe processing.
%
%   nChunks = calculateMaxChunkSize(data, sizeFactor, RAMoverhead)
%
%   This function estimates how many chunks are needed to process data while
%   staying within a conservative RAM budget. The estimate is based on the
%   required memory footprint, a user-defined peak-memory scaling factor, and
%   a reserved RAM overhead.
%
%   Inputs:
%       data :
%           Either:
%               - Numeric array:
%                   The in-memory size is measured with WHOS and multiplied
%                   by sizeFactor.
%
%               - Scalar numeric value:
%                   Interpreted as the estimated required memory footprint
%                   in bytes before applying sizeFactor. The final estimated
%                   requirement is:
%
%                       requiredBytes = data * sizeFactor
%
%       sizeFactor :
%           Positive scalar multiplier used to estimate peak memory use.
%           This should account for temporary arrays, type conversion,
%           intermediate buffers, and output allocation.
%
%       RAMoverhead :
%           Scalar in [0,1). Fraction of total physical RAM reserved for the
%           operating system, MATLAB internals, memory fragmentation, and
%           unrelated temporary allocations. Default is 0.2.
%
%   Output:
%       nChunks :
%           Number of chunks required to keep the estimated memory use within
%           the computed RAM budget. Guaranteed to be at least 1 when usable
%           RAM can be estimated. Returns [] if no usable RAM remains after
%           overhead reservation.
%
%   Notes:
%       - The RAM budget is conservative and clamps usable memory to at most
%         70% of total physical RAM.
%       - For scalar byte estimates, pass sizeFactor = 1 if the estimate
%         already includes all algorithm-specific temporary-memory scaling.


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
