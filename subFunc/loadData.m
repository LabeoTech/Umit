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
%       Info     - For .dat files, AcqInfoStream metadata.
%                  For .umt files, [].
%
%   Notes:
%       - .dat files are reconstructed using folder-global AcqInfos.mat.
%       - .umt files are MAT-files with a custom extension.
%       - Legacy sidecar metadata (.mat or _info.mat) is supported for
%         backward compatibility when present.

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

Info = [];

switch ext
    case '.dat'
        [outFile, Info] = loadDat(fileName);

    case '.umt'
        outFile = load(fileName, '-mat');
        validateUMTStruct(outFile);
end

disp('Done.');
end

% -------------------------------------------------------------------------
function [data, AcqInfoStream] = loadDat(fileName)
%LOADDAT Load .dat file using folder-global AcqInfos.mat.

folderPath = fileparts(fileName);
acqInfoFile = fullfile(folderPath, 'AcqInfos.mat');

if ~isfile(acqInfoFile)
    error('Umitoolbox:loadData:missingAcqInfos', ...
        'Missing "AcqInfos.mat" in folder "%s".', folderPath);
end

S = load(acqInfoFile, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream')
    error('Umitoolbox:loadData:invalidAcqInfos', ...
        '"AcqInfos.mat" does not contain variable "AcqInfoStream".');
end

AcqInfoStream = S.AcqInfoStream;

if ~isstruct(AcqInfoStream) || ~isscalar(AcqInfoStream) || ...
        ~isfield(AcqInfoStream, 'Height') || ~isfield(AcqInfoStream, 'Width')
    error('Umitoolbox:loadData:invalidAcqInfos', ...
        '"AcqInfoStream" must be a scalar struct containing "Height" and "Width".');
end

% Retrocompatibility with old sidecar metadata
[folderPath, baseName, ~] = fileparts(fileName);
legacyFiles = { ...
    fullfile(folderPath, [baseName, '.mat']), ...
    fullfile(folderPath, [baseName, '_info.mat'])};

legacyFiles = legacyFiles(cellfun(@isfile, legacyFiles));

if ~isempty(legacyFiles)
    matInfo = load(legacyFiles{1});

    if isfield(matInfo, 'datSize') && numel(matInfo.datSize) >= 2
        AcqInfoStream.Height = matInfo.datSize(1);
        AcqInfoStream.Width  = matInfo.datSize(2);
    end

    if isfield(matInfo, 'Freq')
        AcqInfoStream.FrameRateHz = matInfo.Freq;
    end
end

fid = fopen(fileName, 'r');
if fid == -1
    error('Umitoolbox:loadData:fileOpenFailed', ...
        'Could not open file "%s".', fileName);
end
cleanupObj = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

data = fread(fid, inf, '*single');

frameSize = double(AcqInfoStream.Height) * double(AcqInfoStream.Width);
if frameSize <= 0 || mod(frameSize, 1) ~= 0
    error('Umitoolbox:loadData:invalidFrameSize', ...
        'Invalid frame dimensions in AcqInfoStream.');
end

if mod(numel(data), frameSize) ~= 0
    error('Umitoolbox:loadData:invalidFileLength', ...
        ['File size is incompatible with frame dimensions from ' ...
         '"AcqInfos.mat".']);
end

nFrames = numel(data) / frameSize;
data = reshape(data, AcqInfoStream.Height, AcqInfoStream.Width, nFrames);
end