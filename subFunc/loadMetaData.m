function Info = loadMetaData(fileName)
%LOADMETADATA Build a unified metadata structure for .dat or .umt files.
%
%   Info = loadMetaData(fileName)
%
%   Inputs:
%       fileName - Full path or relative path to a .dat or .umt file.
%
%   Output:
%       Info     - Metadata structure.
%
%   Supported sources:
%       1) .dat files with folder-global AcqInfos.mat metadata
%       2) .dat files with legacy sidecar metadata (.mat or _info.mat)
%       3) .umt files with embedded metadata and core metadata from
%          AcqInfos.mat when available
%
%   Notes:
%       - All .dat files are assumed to be single precision.
%       - For .dat files, the number of frames is inferred from file size.
%       - For .umt files, embedded metadata is optional.

p = inputParser;
p.FunctionName = 'loadMetaData';
addRequired(p, 'fileName', @(x) ischar(x) || isstring(x));
parse(p, fileName);

fileName = char(string(p.Results.fileName));

if isempty(fileparts(fileName))
    fileName = fullfile(pwd, fileName);
end

if ~isfile(fileName)
    error('Umitoolbox:loadMetaData:fileNotFound', ...
        'File not found: "%s".', fileName);
end

[folderPath, baseName, ext] = fileparts(fileName);
ext = lower(ext);

if ~ismember(ext, {'.dat', '.umt'})
    error('Umitoolbox:loadMetaData:invalidExtension', ...
        'Supported extensions are ".dat" and ".umt".');
end

switch ext
    case '.dat'
        Info = iLoadDatMetaData(fileName, folderPath, baseName);

    case '.umt'
        Info = iLoadUMTMetaData(fileName, folderPath);
end

end

% =========================================================================
% Local helpers
% =========================================================================
function Info = iLoadDatMetaData(fileName, folderPath, baseName)
%ILOADDATMETADATA Build metadata for a .dat file.

AcqInfoStream = iLoadAcqInfo(folderPath);
legacyInfo = iLoadLegacySidecar(folderPath, baseName);

Info = struct();

if ~isempty(AcqInfoStream)
    Info = AcqInfoStream;
    Info.AcqInfoStream = AcqInfoStream;
end

if isfield(legacyInfo, 'datSize') && numel(legacyInfo.datSize) >= 2
    Info.Height = legacyInfo.datSize(1);
    Info.Width  = legacyInfo.datSize(2);
end

if isfield(legacyInfo, 'Freq')
    Info.FrameRateHz = legacyInfo.Freq;
end

if ~isfield(Info, 'Height') || ~isfield(Info, 'Width')
    error('Umitoolbox:loadMetaData:missingFrameSize', ...
        ['Failed to build metadata for "%s". Height and Width were not ' ...
         'found in AcqInfos.mat or legacy sidecar metadata.'], ...
        fileName);
end

Height = double(Info.Height);
Width  = double(Info.Width);

validateattributes(Height, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
    'loadMetaData', 'Height');
validateattributes(Width, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
    'loadMetaData', 'Width');

fileInfo = dir(fileName);
frameSize = Height * Width;
expectedFrameBytes = 4 * frameSize; % single precision only

if mod(fileInfo.bytes, expectedFrameBytes) ~= 0
    error('Umitoolbox:loadMetaData:invalidFileLength', ...
        ['File size is incompatible with frame dimensions for single-precision ' ...
         'data in "%s".'], ...
        fileName);
end

nFrames = fileInfo.bytes / expectedFrameBytes;

Info.fileName = fileName;
Info.folderPath = folderPath;
Info.FileType = '.dat';
Info.datFile = fileName;
Info.datSize = [Height, Width];
Info.datLength = nFrames;
Info.Datatype = 'single';

Info.Height = Height;
Info.Width = Width;
Info.Length = nFrames;

if isfield(Info, 'FrameRateHz')
    Info.Freq = Info.FrameRateHz;
else
    Info.Freq = [];
end

if ~isempty(fieldnames(legacyInfo))
    Info.legacyMetaData = legacyInfo;
end
end

function Info = iLoadUMTMetaData(fileName, folderPath)
%ILOADUMTMETADATA Build metadata for a .umt file.

umt = load(fileName, '-mat');
validateUMTStruct(umt);

Info = struct();

if isfield(umt, 'metaData') && isstruct(umt.metaData) && isscalar(umt.metaData)
    Info = umt.metaData;
elseif isfield(umt, 'metadata') && isstruct(umt.metadata) && isscalar(umt.metadata)
    Info = umt.metadata;
elseif isfield(umt, 'Info') && isstruct(umt.Info) && isscalar(umt.Info)
    Info = umt.Info;
end

AcqInfoStream = iLoadAcqInfo(folderPath);
if ~isempty(AcqInfoStream)
    coreFields = {'Height', 'Width', 'Length', 'FrameRateHz', 'ExposureMsec'};

    for ii = 1:numel(coreFields)
        if isfield(AcqInfoStream, coreFields{ii})
            Info.(coreFields{ii}) = AcqInfoStream.(coreFields{ii});
        end
    end

    Info.AcqInfoStream = AcqInfoStream;

    if ~isfield(Info, 'datSize') && isfield(AcqInfoStream, 'Height') && isfield(AcqInfoStream, 'Width')
        Info.datSize = [AcqInfoStream.Height, AcqInfoStream.Width];
    end

    if ~isfield(Info, 'datLength') && isfield(AcqInfoStream, 'Length')
        Info.datLength = AcqInfoStream.Length;
    end

    if ~isfield(Info, 'Freq') && isfield(AcqInfoStream, 'FrameRateHz')
        Info.Freq = AcqInfoStream.FrameRateHz;
    end
end

Info.fileName = fileName;
Info.folderPath = folderPath;
Info.FileType = '.umt';
end

function AcqInfoStream = iLoadAcqInfo(folderPath)
%ILOADACQINFO Load AcqInfoStream from AcqInfos.mat when available.

AcqInfoStream = struct();

acqInfoFile = fullfile(folderPath, 'AcqInfos.mat');
if ~isfile(acqInfoFile)
    return
end

S = load(acqInfoFile, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream')
    error('Umitoolbox:loadMetaData:invalidAcqInfos', ...
        '"AcqInfos.mat" does not contain variable "AcqInfoStream".');
end

AcqInfoStream = S.AcqInfoStream;

if ~isstruct(AcqInfoStream) || ~isscalar(AcqInfoStream)
    error('Umitoolbox:loadMetaData:invalidAcqInfos', ...
        '"AcqInfoStream" must be a scalar struct.');
end
end

function legacyInfo = iLoadLegacySidecar(folderPath, baseName)
%ILOADLEGACYSIDECAR Load first available legacy metadata sidecar.

legacyInfo = struct();

legacyFiles = { ...
    fullfile(folderPath, [baseName, '.mat']), ...
    fullfile(folderPath, [baseName, '_info.mat'])};

legacyFiles = legacyFiles(cellfun(@isfile, legacyFiles));

if ~isempty(legacyFiles)
    legacyInfo = load(legacyFiles{1});
end
end