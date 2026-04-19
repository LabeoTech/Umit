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
%       - For raw and derived .dat files, legacy-compatible fields such as
%         dim_names, datSize, datLength, Freq, and Datatype are always
%         provided.
%       - For .dat files, Height and Width are strict metadata, but the
%         temporal length is inferred from actual file size.
%       - Valid legacy event-split .dat metadata are preserved and are not
%         collapsed to continuous YXT.
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
if ~isfield(Info, 'FrameRateHz') && isfield(Info, 'Freq')
    Info.FrameRateHz = Info.Freq;
end

if ~isfield(Info, 'Freq') && isfield(Info, 'FrameRateHz')
    Info.Freq = Info.FrameRateHz;
end

if ~isfield(Info, 'dim_names') || isempty(Info.dim_names)
    Info.dim_names = {'Y', 'X', 'T'};
else
    Info.dim_names = cellstr(string(Info.dim_names));
end

if ~isfield(Info, 'Datatype') || isempty(Info.Datatype)
    Info.Datatype = 'single';
end

assert(strcmpi(char(string(Info.Datatype)), 'single'), ...
    'Umitoolbox:loadMetaData:unsupportedDatatype', ...
    'Only single-precision .dat files are currently supported.');

% Preserve valid legacy datSize dimensionality. Do not collapse event-split
% metadata to YX. Only synthesize [Height Width] when datSize is absent.
if isfield(Info, 'datSize') && ~isempty(Info.datSize)
    Info.datSize = double(Info.datSize(:).');
end

if (~isfield(Info, 'datSize') || isempty(Info.datSize)) && ...
        isfield(Info, 'Height') && isfield(Info, 'Width')
    Info.datSize = [double(Info.Height), double(Info.Width)];
end

% Resolve Height and Width from dim_names + datSize when possible. This
% supports both continuous YXT and legacy event-split metadata.
if isfield(Info, 'datSize') && ~isempty(Info.datSize)
    idxY = find(strcmp(Info.dim_names, 'Y'), 1, 'first');
    idxX = find(strcmp(Info.dim_names, 'X'), 1, 'first');

    if numel(Info.datSize) == numel(Info.dim_names)
        if ~isfield(Info, 'Height') && ~isempty(idxY)
            Info.Height = Info.datSize(idxY);
        end
        if ~isfield(Info, 'Width') && ~isempty(idxX)
            Info.Width = Info.datSize(idxX);
        end
    elseif numel(Info.datSize) == numel(Info.dim_names) - 1
        % Common case where datSize stores only non-T dimensions.
        if ~isfield(Info, 'Height') && numel(Info.datSize) >= 1
            Info.Height = Info.datSize(1);
        end
        if ~isfield(Info, 'Width') && numel(Info.datSize) >= 2
            Info.Width = Info.datSize(2);
        end
    end
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

% Always infer datLength from the actual file size. This keeps loading
% strict on file integrity but flexible on temporal length for derived
% files such as Y X (T-1) outputs, while preserving legacy event-split
% dimensionality when valid metadata are present.
fileInfo = dir(fileName);
bytesPerElement = getByteSize('single');

idxT = find(strcmp(Info.dim_names, 'T'), 1, 'first');
if isempty(idxT)
    error('Umitoolbox:loadMetaData:missingTimeDimension', ...
        'dim_names must contain a T dimension for .dat files.');
end

if ~isfield(Info, 'datSize') || isempty(Info.datSize)
    nonTProd = Height * Width;
elseif numel(Info.datSize) == numel(Info.dim_names)
    nonTProd = prod(Info.datSize(setdiff(1:numel(Info.datSize), idxT)));
elseif numel(Info.datSize) == numel(Info.dim_names) - 1
    nonTProd = prod(Info.datSize);
else
    error('Umitoolbox:loadMetaData:invalidDatSize', ...
        ['datSize for "%s" is incompatible with dim_names. Expected either ' ...
         'all dimensions or all non-T dimensions.'], ...
        fileName);
end

if mod(fileInfo.bytes, nonTProd * bytesPerElement) ~= 0
    error('Umitoolbox:loadMetaData:invalidFileLength', ...
        ['File size is incompatible with declared non-T dimensions for ' ...
         'single-precision data in "%s".'], ...
        fileName);
end

actualLength = fileInfo.bytes / (nonTProd * bytesPerElement);
legacyLength = [];
if isfield(Info, 'datLength') && ~isempty(Info.datLength)
    legacyLength = double(Info.datLength);
end

Info.datLength = actualLength;

% Preserve any original declared length when it differs from the actual
% on-disk temporal length. This is useful for derived files.
if ~isempty(legacyLength) && isfinite(legacyLength) && legacyLength ~= actualLength
    Info.OriginalLength = legacyLength;
end

if isfield(Info, 'Length') && ~isempty(Info.Length) && double(Info.Length) ~= actualLength
    Info.OriginalLength = double(Info.Length);
end
Info.Length = actualLength;

% If datSize explicitly includes T, keep the multi-dimensional layout but
% update the T slot to the actual on-disk length.
if isfield(Info, 'datSize') && ~isempty(Info.datSize) && numel(Info.datSize) == numel(Info.dim_names)
    Info.datSize(idxT) = actualLength;
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
        Info.datSize = dimSizes(idxKeep);
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
