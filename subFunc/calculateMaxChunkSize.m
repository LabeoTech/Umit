function nChunks = calculateMaxChunkSize(data,sizeFactor)
% CALCULATEMAXCHUNKSIZE reads the available memory and calculates the
% number of chunks of the data to be processed by a given function.
% Inputs:
%   data(numerical array): data to be processed.
%   sizeFactor(scalar): multiplication factor of the size of data used to
%       estimate the amount of memory to be used by the data processing.
% Outputs:
%   nChunks(scalar): number of chunks of data to be processed in order to
%   avoid out-of-memory errors.


% Calculate the data byte size
s = whos('data');
% clear data;
dataByteSize = s.bytes;
maxBytes = sizeFactor*dataByteSize;
% Get available RAM:
switch computer
    case 'PCWIN64'
        [~,sys]= memory;
        bytesAvailable = sys.PhysicalMemory.Available;
    case 'GLNXA64'
        [status, result] = system('free -b | grep Mem | awk ''{print $7}''');
        if status
            error('Failed to retrieve RAM information in Linux System!')
        end
        bytesAvailable = str2double(result);
    case 'MACI64'
        [status, result] = system('vm_stat | grep "Pages free" | awk ''{print $3}''');
        if status
            error('Failed to retrieve RAM information in MacOS!')
        end
        page_size = str2double(system('getconf PAGESIZE'));
        bytesAvailable = str2double(result)*page_size;
    otherwise
        error('Unsupported Platform!');
end

% Calculate the number of chunks to be used given the available memory:
nChunks = ceil(maxBytes/bytesAvailable);
end
