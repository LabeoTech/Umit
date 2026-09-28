function [outFile, Info] = loadData(fileName)
%LOADDATA Load raw .dat data or derived .umt data.
%
%   [outFile, Info] = loadData(fileName)
%
%   Inputs:
%       fileName - Full path to a .dat or .umt file.
%
%   Outputs:
%       outFile  - For .dat files, numeric array loaded from disk.
%                  For .umt files, struct loaded from the MAT-backed file.
%       Info     - Unified metadata structure from loadMetaData.
%
%   Notes:
%       - .dat metadata are resolved through loadMetaData.
%       - Headered .dat files are read from their data offset in the class
%         stored in the header and reshaped to Info.dimSizes.
%       - Headerless .dat files (legacy sidecar or AcqInfos-bound) are
%         read as single precision from byte 0 and reshaped to
%         Y x X x frames.
%       - .umt files are MAT-files with a custom extension.

p = inputParser;
p.FunctionName = 'loadData';
addRequired(p, 'fileName', @(x) ischar(x) || isstring(x));
parse(p, fileName);

fileName = char(string(p.Results.fileName));

if isempty(fileparts(fileName))
    fileName = fullfile(pwd, fileName);
end

if ~isfile(fileName)
    error('Umitoolbox:loadData:fileNotFound', ...
        'File not found: "%s".', fileName);
end

[~, ~, ext] = fileparts(fileName);
ext = lower(ext);

if ~ismember(ext, {'.dat', '.umt'})
    error('Umitoolbox:loadData:invalidExtension', ...
        'Supported extensions are ".dat" and ".umt".');
end

fprintf('Opening file "%s" ...\n', fileName);

Info = loadMetaData(fileName);

switch ext
    case '.dat'
        outFile = iLoadDat(fileName, Info);

    case '.umt'
        outFile = load(fileName, '-mat');
        validateUMTStruct(outFile);
end

disp('Done.');
end

% =========================================================================
% Local helpers
% =========================================================================
function data = iLoadDat(fileName, Info)
%ILOADDAT Load .dat file using unified metadata from loadMetaData.

fid = fopen(fileName, 'r');
if fid == -1
    error('Umitoolbox:loadData:fileOpenFailed', ...
        'Could not open file "%s".', fileName);
end
cleanupObj = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

if strcmp(Info.format, 'header')
    % Self-describing file: read exactly the described array, in its stored
    % class, starting after the header.
    nValues = prod(Info.dimSizes);
    fseek(fid, Info.dataOffset, 'bof');
    data = fread(fid, nValues, ['*' Info.dataClass], 0, 'ieee-le');
    if numel(data) ~= nValues
        error('Umitoolbox:loadData:invalidFileLength', ...
            'File "%s" holds %d of the %d values its header describes.', ...
            fileName, numel(data), nValues);
    end
    data = reshape(data, Info.dimSizes);
    return
end

% Headerless files (legacy sidecar or AcqInfos-bound): single precision
% from byte 0, reshaped to Y x X x frames.
data = fread(fid, inf, '*single');

frameSize = datAxisSize(Info, 'Y') * datAxisSize(Info, 'X');
if frameSize <= 0 || mod(frameSize, 1) ~= 0
    error('Umitoolbox:loadData:invalidFrameSize', ...
        'Invalid frame dimensions in metadata.');
end

if mod(numel(data), frameSize) ~= 0
    error('Umitoolbox:loadData:invalidFileLength', ...
        'File size is incompatible with metadata frame dimensions.');
end

nFrames = numel(data) / frameSize;

% Use the actual on-disk temporal length. loadMetaData already derives
% datLength from file size, but recompute here to keep loading strict.
data = reshape(data, datAxisSize(Info, 'Y'), datAxisSize(Info, 'X'), nFrames);
end
