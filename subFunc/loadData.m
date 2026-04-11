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
%       - .dat files are assumed to be single precision.
%       - .dat metadata are resolved through loadMetaData.
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

data = fread(fid, inf, '*single');

frameSize = double(Info.Height) * double(Info.Width);
if frameSize <= 0 || mod(frameSize, 1) ~= 0
    error('Umitoolbox:loadData:invalidFrameSize', ...
        'Invalid frame dimensions in metadata.');
end

if mod(numel(data), frameSize) ~= 0
    error('Umitoolbox:loadData:invalidFileLength', ...
        'File size is incompatible with metadata frame dimensions.');
end

data = reshape(data, Info.Height, Info.Width, Info.datLength);
end