function [availBytes, totalBytes] = queryRAM()
%QUERYRAM Return available and total physical RAM in bytes.
%
%   [availBytes, totalBytes] = queryRAM()
%
% Outputs:
%   availBytes - Estimated currently available physical memory in bytes.
%   totalBytes - Total physical memory in bytes.
%
% Notes:
%   - On Windows, this uses MEMORY.
%   - On Linux, this parses /proc/meminfo.
%   - On macOS, this uses vm_stat and sysctl.
%   - If RAM cannot be queried, outputs are [].

    availBytes = [];
    totalBytes = [];

    try
        if ispc
            [~, s] = memory();

            if isstruct(s)
                if isfield(s, 'PhysicalMemory') && isstruct(s.PhysicalMemory)

                    if isfield(s.PhysicalMemory, 'Available')
                        availBytes = double(s.PhysicalMemory.Available);
                    end

                    if isfield(s.PhysicalMemory, 'Total')
                        totalBytes = double(s.PhysicalMemory.Total);
                    end
                end
            end
            return
        end

        if isunix && ~ismac
            meminfoPath = '/proc/meminfo';
            if isfile(meminfoPath)
                txt = fileread(meminfoPath);

                tokAvail = regexp(txt, 'MemAvailable:\s+(\d+)\s+kB', 'tokens', 'once');
                tokTotal = regexp(txt, 'MemTotal:\s+(\d+)\s+kB', 'tokens', 'once');

                if ~isempty(tokAvail)
                    availBytes = str2double(tokAvail{1}) * 1024;
                end
                if ~isempty(tokTotal)
                    totalBytes = str2double(tokTotal{1}) * 1024;
                end
            end
            return
        end

        if ismac
            [status1, pageSizeTxt] = system('sysctl -n hw.pagesize');
            [status2, totalMemTxt] = system('sysctl -n hw.memsize');
            [status3, vmStatTxt]   = system('vm_stat');

            if status1 == 0
                pageSize = str2double(strtrim(pageSizeTxt));
            else
                pageSize = 4096;
            end

            if status2 == 0
                totalBytes = str2double(strtrim(totalMemTxt));
            end

            if status3 == 0
                freeTok   = regexp(vmStatTxt, 'Pages free:\s+(\d+)\.', 'tokens', 'once');
                specTok   = regexp(vmStatTxt, 'Pages speculative:\s+(\d+)\.', 'tokens', 'once');
                inactiveTok = regexp(vmStatTxt, 'Pages inactive:\s+(\d+)\.', 'tokens', 'once');

                freePages = 0;
                specPages = 0;
                inactivePages = 0;

                if ~isempty(freeTok),     freePages = str2double(freeTok{1}); end
                if ~isempty(specTok),     specPages = str2double(specTok{1}); end
                if ~isempty(inactiveTok), inactivePages = str2double(inactiveTok{1}); end

                availBytes = double(freePages + specPages + inactivePages) * double(pageSize);
            end
            return
        end

    catch
        availBytes = [];
        totalBytes = [];
    end
end