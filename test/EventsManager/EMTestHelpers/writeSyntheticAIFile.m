function filePath = writeSyntheticAIFile(folderPath, data, varargin)
%WRITESYNTHETICAIFILE Write a synthetic ai_00000.bin file.
%
% The file layout matches EventsManager/setAnalogIN expectations:
%   - 5 x int32 header values (20 bytes)
%   - double payload
%   - payload length compatible with blocks of 1e4 x nChannels x nBlocks
%
% Syntax:
%   filePath = writeSyntheticAIFile(folderPath, data)
%   filePath = writeSyntheticAIFile(..., 'FileName', 'ai_00000.bin')
%
% Inputs:
%   folderPath - Destination folder.
%   data       - Numeric matrix of size Nsamples x Nchannels.
%
% Output:
%   filePath   - Full path to the written file.

    p = inputParser;
    addRequired(p, 'folderPath', @(x) ischar(x) || isStringScalar(x));
    addRequired(p, 'data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addParameter(p, 'FileName', 'ai_00000.bin', @(x) ischar(x) || isStringScalar(x));
    parse(p, folderPath, data, varargin{:});

    folderPath = convertStringsToChars(folderPath);
    fileName = convertStringsToChars(p.Results.FileName);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    data = double(p.Results.data);
    [nSamples, nChan] = size(data);
    blockLen = 1e4;
    nPad = mod(-nSamples, blockLen);
    if nPad > 0
        data = [data; zeros(nPad, nChan)]; %#ok<AGROW>
    end

    filePath = fullfile(folderPath, fileName);
    fid = fopen(filePath, 'w');
    assert(fid ~= -1, 'Failed to open %s for writing.', filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fwrite(fid, zeros(5,1,'int32'), 'int32');

    tmp = reshape(data, blockLen, [], nChan);
    tmp = permute(tmp, [1 3 2]);
    fwrite(fid, tmp(:), 'double');
end
