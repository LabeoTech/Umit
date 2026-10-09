function [userview, systemview] = memory
%MEMORY Test-only shadow of the PCWIN64 built-in `memory` function.
%
%   Used exclusively by TestCalculateMaxChunkSize (and any caller-level test
%   simulating RAM exhaustion, e.g. TestHemoCompute) to force specific RAM
%   query outcomes without touching real system state. Only takes effect
%   while this folder is prepended to the path via
%   matlab.unittest.fixtures.PathFixture; outside tests the real built-in is
%   used as normal.
%
%   Controlled via environment variables set by the calling test:
%       UMIT_TEST_MOCK_MEMORY_MODE      'fail' or 'fixed'
%       UMIT_TEST_MOCK_MEMORY_TOTAL     total physical RAM in bytes (mode 'fixed')
%       UMIT_TEST_MOCK_MEMORY_AVAILABLE available RAM in bytes (mode 'fixed')

mode = getenv('UMIT_TEST_MOCK_MEMORY_MODE');

switch mode
    case 'fail'
        error('UmitTest:MockMemory:SimulatedFailure', ...
            'Simulated RAM query failure for testing.');
    case 'fixed'
        totalRAM = str2double(getenv('UMIT_TEST_MOCK_MEMORY_TOTAL'));
        availableRAM = str2double(getenv('UMIT_TEST_MOCK_MEMORY_AVAILABLE'));
    otherwise
        % Mode not set: this folder is on the path outside a test (e.g.
        % addpath(genpath(pwd))). Delegate to the real built-in.
        % BUILTIN bypasses this shadow. (Locating the real function with
        % WHICH and CD-ing to its folder fails: WHICH reports it as
        % 'built-in (<path>)', which is not a folder.)
        if nargout == 0
            builtin('memory');
        else
            [userview, systemview] = builtin('memory');
        end
        return
end

userview = struct('MemUsedMATLAB', 0, 'MemAvailableAllArrays', availableRAM);
systemview = struct('PhysicalMemory', struct('Total', totalRAM, 'Available', availableRAM));

end
