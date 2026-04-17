function Info = loadMetaData(fileName)
%LOADMETADATA Build a flat compatibility metadata structure for .dat or .umt files.
%
%   Info = loadMetaData(fileName)
%
%   Inputs:
%       fileName - Full path or relative path to a .dat or .umt file.
%
%   Output:
%       Info     - Flat metadata structure.
%
%   Supported sources:
%       1) .dat files with legacy sidecar metadata (.mat or _info.mat)
%       2) .dat files with folder-global AcqInfos.mat metadata
%       3) .umt files with embedded metadata and/or AcqInfos.mat
%
%   Compatibility rules for .dat files:
%       - If a legacy sidecar metadata file is found and recognized, its
%         fields are copied first and take precedence.
%       - Fields from AcqInfoStream are then appended only when they are
%         missing from the legacy metadata.
%       - If no valid legacy metadata is found, the output is built from
%         AcqInfoStream and legacy-compatible fields are derived from it.
%       - All .dat files are assumed to be single precision unless a
%         legacy metadata file already defines Datatype.
%
%   Notes:
%       - The output is always a flat structure. No nested metadata
%         structures are returned.
%       - For raw .dat files, legacy-compatible fields such as dim_names,
%         datSize, datLength, Freq, and Datatype are always provided.
%       - For .umt files, embedded metadata is optional. When missing, core
%         metadata are derived from the first data entry when possible.

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
%ILOADDATMETADATA Build flat metadata for a .dat file.

acqInfo = iLoadAcqInfo(folderPath);
legacyInfo = iLoadLegacySidecar(folderPath, baseName);

if ~isempty(fieldnames(legacyInfo))
    % Legacy metadata takes precedence. Append only missing fields from
    % AcqInfoStream to preserve branch-specific semantics.
    Info = legacyInfo;
    Info = iAppendMissingFields(Info, acqInfo);
else
    % Newer branch: build from AcqInfoStream and derive legacy-compatible fields.
    Info = acqInfo;
end

% -------------------------------------------------------------------------
% Derive/complete legacy-compatible fields
% -------------------------------------------------------------------------
if ~isfield(Info, 'Height') && isfield(Info, 'datSize') && numel(Info.datSize) >= 1
    Info.Height = Info.datSize(1);
end

if ~isfield(Info, 'Width') && isfield(Info, 'datSize') && numel(Info.datSize) >= 2
    Info.Width = Info.datSize(2);
end

if ~isfield(Info, 'FrameRateHz') && isfield(Info, 'Freq')
    Info.FrameRateHz = Info.Freq;
end

if ~isfield(Info, 'Freq') && isfield(Info, 'FrameRateHz')
    Info.Freq = Info.FrameRateHz;
end

if ~isfield(Info, 'Height') || ~isfield(Info, 'Width')
    error('Umitoolbox:loadMetaData:missingFrameSize', ...
        ['Failed to build metadata for "%s". Height and Width were not ' ...
         'found in legacy metadata or AcqInfos.mat.'], ...
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

if ~isfield(Info, 'datSize') || isempty(Info.datSize)
    Info.datSize = [Height, Width];
end

if ~isfield(Info, 'dim_names') || isempty(Info.dim_names)
    Info.dim_names = {'Y', 'X', 'T'};
end

if ~isfield(Info, 'Datatype') || isempty(Info.Datatype)
    Info.Datatype = 'single';
end

% Always validate/infer datLength from actual file size if missing.
if ~isfield(Info, 'datLength') || isempty(Info.datLength)
    fileInfo = dir(fileName);
    frameSize = Height * Width;
    bytesPerElement = getByteSize('single');

    if mod(fileInfo.bytes, frameSize * bytesPerElement) ~= 0
        error('Umitoolbox:loadMetaData:invalidFileLength', ...
            ['File size is incompatible with frame dimensions for single-precision ' ...
             'data in "%s".'], ...
            fileName);
    end

    Info.datLength = fileInfo.bytes / (frameSize * bytesPerElement);
end

if ~isfield(Info, 'Length') || isempty(Info.Length)
    Info.Length = Info.datLength;
end

% Flat convenience fields
Info.fileName = fileName;
Info.folderPath = folderPath;
Info.FileType = '.dat';
Info.datFile = fileName;

end

function Info = iLoadUMTMetaData(fileName, folderPath)
%ILOADUMTMETADATA Build flat metadata for a .umt file.

umt = iLoadUMTFromFile(fileName);
validateUMTStruct(umt, 'requireEventInfo', false);

Info = struct();

% Embedded metadata, if any, comes first.
embedded = iExtractEmbeddedMetadata(umt);
Info = iAppendMissingFields(Info, embedded);

% Then append missing fields from AcqInfos.mat.
acqInfo = iLoadAcqInfo(folderPath);
Info = iAppendMissingFields(Info, acqInfo);

% Derive compatibility fields from the first entry when possible.
entryNames = fieldnames(umt.data);
if ~isempty(entryNames)
    firstEntry = umt.data.(entryNames{1});
    dimNames = cellstr(string(firstEntry.dimNames));
    dimSizes = iGetDeclaredDimensionSizes(firstEntry.value, dimNames);

    if ~isfield(Info, 'dim_names') || isempty(Info.dim_names)
        Info.dim_names = dimNames;
    end

    if ~isfield(Info, 'datSize') || isempty(Info.datSize)
        idxKeep = ~strcmp(dimNames, 'T');
        if any(strcmp(dimNames, 'E'))
            Info.datSize = dimSizes(idxKeep);
        else
            Info.datSize = dimSizes(idxKeep);
        end
    end

    idxT = find(strcmp(dimNames, 'T'), 1, 'first');
    if ~isempty(idxT) && (~isfield(Info, 'datLength') || isempty(Info.datLength))
        Info.datLength = dimSizes(idxT);
    end
end

if ~isfield(Info, 'Height') && isfield(Info, 'datSize') && numel(Info.datSize) >= 1
    Info.Height = Info.datSize(1);
end

if ~isfield(Info, 'Width') && isfield(Info, 'datSize') && numel(Info.datSize) >= 2
    Info.Width = Info.datSize(2);
end

if ~isfield(Info, 'Length') && isfield(Info, 'datLength')
    Info.Length = Info.datLength;
end

if ~isfield(Info, 'FrameRateHz') && isfield(Info, 'Freq')
    Info.FrameRateHz = Info.Freq;
end

if ~isfield(Info, 'Freq') && isfield(Info, 'FrameRateHz')
    Info.Freq = Info.FrameRateHz;
end

if ~isfield(Info, 'Datatype') || isempty(Info.Datatype)
    Info.Datatype = 'single';
end

Info.fileName = fileName;
Info.folderPath = folderPath;
Info.FileType = '.umt';

end

function acqInfo = iLoadAcqInfo(folderPath)
%ILOADACQINFO Load and flatten AcqInfoStream from AcqInfos.mat when available.

acqInfo = struct();

acqInfoFile = fullfile(folderPath, 'AcqInfos.mat');
if ~isfile(acqInfoFile)
    return
end

S = load(acqInfoFile, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream')
    error('Umitoolbox:loadMetaData:invalidAcqInfos', ...
        '"AcqInfos.mat" does not contain variable "AcqInfoStream".');
end

acqInfo = S.AcqInfoStream;

if ~isstruct(acqInfo) || ~isscalar(acqInfo)
    error('Umitoolbox:loadMetaData:invalidAcqInfos', ...
        '"AcqInfoStream" must be a scalar struct.');
end

end

function legacyInfo = iLoadLegacySidecar(folderPath, baseName)
%ILOADLEGACYSIDECAR Load the first valid legacy metadata sidecar.

legacyInfo = struct();

candidateFiles = { ...
    fullfile(folderPath, [baseName, '.mat']), ...
    fullfile(folderPath, [baseName, '_info.mat'])};

for iFile = 1:numel(candidateFiles)
    if ~isfile(candidateFiles{iFile})
        continue
    end

    raw = load(candidateFiles{iFile});
    candidate = iExtractLegacyMetadata(raw);

    if ~isempty(fieldnames(candidate))
        legacyInfo = candidate;
        return
    end
end

end

function out = iExtractLegacyMetadata(raw)
%IEXTRACTLEGACYMETADATA Extract a flat legacy metadata struct when recognized.

out = struct();

if iLooksLikeLegacyMetadata(raw)
    out = raw;
    return
end

fn = fieldnames(raw);
for iField = 1:numel(fn)
    candidate = raw.(fn{iField});
    if isstruct(candidate) && isscalar(candidate) && iLooksLikeLegacyMetadata(candidate)
        out = candidate;
        return
    end
end

end

function tf = iLooksLikeLegacyMetadata(S)
%ILOOKSLIKELEGACYMETADATA Heuristic test for legacy metadata payloads.

if ~isstruct(S) || ~isscalar(S)
    tf = false;
    return
end

anchorFields = {'dim_names','datSize','datLength','Freq','Datatype'};
tf = sum(isfield(S, anchorFields)) >= 2;

end

function out = iAppendMissingFields(out, src)
%IAPPENDMISSINGFIELDS Append fields from src into out without overwriting.

if isempty(src) || ~isstruct(src) || ~isscalar(src)
    return
end

srcFields = fieldnames(src);
for iField = 1:numel(srcFields)
    if ~isfield(out, srcFields{iField})
        out.(srcFields{iField}) = src.(srcFields{iField});
    end
end

end

function embedded = iExtractEmbeddedMetadata(umt)
%IEXTRACTEMBEDDEDMETADATA Extract flat embedded metadata from a UMT struct.

embedded = struct();

candidateFields = {'metaData','metadata','Info'};
for iField = 1:numel(candidateFields)
    if isfield(umt, candidateFields{iField}) && ...
            isstruct(umt.(candidateFields{iField})) && ...
            isscalar(umt.(candidateFields{iField}))
        embedded = umt.(candidateFields{iField});
        return
    end
end

end

function umt = iLoadUMTFromFile(fileName)
%ILOADUMTFROMFILE Load the first scalar UMT struct found in a .umt file.

S = load(fileName, '-mat');

if isstruct(S) && isscalar(S) && all(ismember({'version','kind','data'}, fieldnames(S)))
    umt = S;
    return
end

fn = fieldnames(S);
umt = [];
for iField = 1:numel(fn)
    candidate = S.(fn{iField});
    if isstruct(candidate) && isscalar(candidate) && ...
            all(ismember({'version','kind','data'}, fieldnames(candidate)))
        umt = candidate;
        break
    end
end

if isempty(umt)
    error('Umitoolbox:loadMetaData:invalidUMT', ...
        'No scalar UMT struct was found in "%s".', fileName);
end

end

function dimSizes = iGetDeclaredDimensionSizes(value, dimNames)
%IGETDECLAREDDIMENSIONSIZES Return sizes compatible with declared dimNames.

nDimsExpected = numel(dimNames);
sz = size(value);

if numel(sz) < nDimsExpected
    sz(end+1:nDimsExpected) = 1;
elseif numel(sz) > nDimsExpected
    sz = sz(1:nDimsExpected);
end

dimSizes = sz;

end