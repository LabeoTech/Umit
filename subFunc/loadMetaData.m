function Info = loadMetaData(fileName)
%LOADMETADATA Build concise file-facing metadata for .dat or .umt files.
%
%   Info = loadMetaData(fileName)
%
%   Inputs:
%       fileName - Full path or relative path to a .dat or .umt file.
%
%   Output:
%       Info     - File-facing metadata structure.
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
%       - The output is a compact flat structure with fields that describe
%         the associated file directly. Acquisition-session fields such as
%         analog input, stimulation, and illumination definitions are not
%         copied to Info.
%       - For .dat files, legacy-compatible fields such as dim_names,
%         datSize, datLength, Freq, and Datatype are always provided.
%       - For .dat files, Height and Width are strict metadata. The
%         temporal length is inferred from actual file size and, for
%         AcqInfos-driven data without legacy sidecars, must match one
%         imported/base timeline described by AcqInfos.mat.
%       - Valid legacy event-split .dat metadata are preserved and are not
%         collapsed to continuous YXT.
%       - For .umt files, embedded metadata is optional. When missing, core
%         metadata are derived from the first data entry when possible.
%       - When legacy sidecar metadata are used for .dat files, the fields
%         FrameRateHz, Width, Height, and Length are refreshed from the
%         legacy fields Freq, datSize, datLength, and dim_names for forward
%         compatibility.

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
hasLegacySidecar = ~isempty(fieldnames(legacyInfo));

if hasLegacySidecar
    % Legacy metadata takes precedence. Append only missing fields from
    % AcqInfoStream to preserve source-specific semantics.
    Info = legacyInfo;
    Info = iAppendMissingFields(Info, acqInfo);

    % Refresh forward-compatible core fields from the legacy payload.
    [Info, updatedFields] = iUpdateLegacyCoreFieldsForForwardCompatibility(Info);

    if ~isempty(updatedFields)
        fprintf(['loadMetaData: Updated legacy metadata field(s) for ' ...
            'forward compatibility in "%s": %s\n'], ...
            fileName, strjoin(updatedFields, ', '));
    end
else
    % Current metadata model: build from AcqInfoStream and derive legacy-compatible fields.
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
% strict on file integrity. For AcqInfos-driven .dat files without
% legacy sidecars, the inferred temporal length must also match one known
% imported/base timeline. Legacy sidecar metadata are preserved for
% backwards compatibility, including legacy event-split dimensionality.
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
Info.datLength = actualLength;
Info.Length = actualLength;

% Current data-format rule: AcqInfos-driven .dat files are valid only when their
% inferred T matches one known imported/base timeline. Legacy sidecar files
% are kept backwards-compatible and are not constrained by ImportedChannels.
if ~hasLegacySidecar
    timelineInfo = resolveDatTimeline(actualLength, acqInfo, 'DatFile', fileName);
    Info.FrameRateHz = timelineInfo.FrameRateHz;
    Info.Freq = timelineInfo.Freq;
    Info.TimelineSource = timelineInfo.SourceType;
    Info.TimelineSourceIndex = timelineInfo.SourceIndex;
end

% If datSize explicitly includes T, keep the multi-dimensional layout but
% update the T slot to the actual on-disk length.
if isfield(Info, 'datSize') && ~isempty(Info.datSize) && numel(Info.datSize) == numel(Info.dim_names)
    Info.datSize(idxT) = actualLength;
end

% Keep only metadata that directly describes this .dat file.
importedEntry = iFindImportedChannelForInfo(acqInfo, fileName, actualLength);
Info = iFinalizeDatInfo(Info, acqInfo, fileName, folderPath, actualLength, ...
    importedEntry, hasLegacySidecar);

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


function importedEntry = iFindImportedChannelForInfo(acqInfo, fileName, actualLength)
%IFINDIMPORTEDCHANNELFORINFO Return the imported-channel entry for Info.
%
% The exact DatFile match is authoritative. If there is no exact match, a
% unique length match is used only when it identifies one imported channel.

importedEntry = struct();

if isempty(acqInfo) || ~isstruct(acqInfo) || ~isscalar(acqInfo) || ...
        ~isfield(acqInfo, 'ImportedChannels') || isempty(acqInfo.ImportedChannels)
    return
end

raw = acqInfo.ImportedChannels(:).';
[~, datBase, datExt] = fileparts(fileName);
if isempty(datExt)
    datExt = '.dat';
end
datName = [datBase, datExt];

if isfield(raw, 'DatFile')
    datFiles = cell(1, numel(raw));
    for iEntry = 1:numel(raw)
        [~, thisBase, thisExt] = fileparts(char(string(raw(iEntry).DatFile)));
        if isempty(thisExt)
            thisExt = '.dat';
        end
        datFiles{iEntry} = [thisBase, thisExt];
    end

    idxFile = find(strcmpi(datFiles, datName));
    if numel(idxFile) == 1
        importedEntry = raw(idxFile);
        return
    end
end

if isfield(raw, 'Length')
    lenList = nan(1, numel(raw));
    for iEntry = 1:numel(raw)
        if ~isempty(raw(iEntry).Length)
            lenList(iEntry) = double(raw(iEntry).Length);
        end
    end

    idxLength = find(lenList == double(actualLength));
    if numel(idxLength) == 1
        importedEntry = raw(idxLength);
    end
end

end

function Info = iFinalizeDatInfo(rawInfo, acqInfo, fileName, folderPath, actualLength, importedEntry, hasLegacySidecar)
%IFINALIZEDATINFO Keep only file-facing metadata fields for .dat files.

Info = struct();

Info.datFile = fileName;
Info.folderPath = folderPath;
Info.FileType = '.dat';

Info.Height = double(rawInfo.Height);
Info.Width = double(rawInfo.Width);
Info.Length = double(actualLength);

if isfield(rawInfo, 'FrameRateHz') && ~isempty(rawInfo.FrameRateHz)
    Info.FrameRateHz = double(rawInfo.FrameRateHz);
elseif isfield(rawInfo, 'Freq') && ~isempty(rawInfo.Freq)
    Info.FrameRateHz = double(rawInfo.Freq);
end

if ~isfield(Info, 'FrameRateHz') || isempty(Info.FrameRateHz)
    error('Umitoolbox:loadMetaData:missingFrameRate', ...
        'Failed to resolve FrameRateHz for "%s".', fileName);
end

Info.Datatype = char(string(rawInfo.Datatype));
Info.dim_names = cellstr(string(rawInfo.dim_names));

if isfield(rawInfo, 'datSize') && ~isempty(rawInfo.datSize)
    Info.datSize = double(rawInfo.datSize(:).');
else
    Info.datSize = [Info.Height, Info.Width];
end

Info.datLength = Info.Length;
Info.Freq = Info.FrameRateHz;

if isfield(rawInfo, 'datName') && ~isempty(rawInfo.datName)
    Info.datName = rawInfo.datName;
else
    Info.datName = 'data';
end

if isfield(rawInfo, 'FirstDim') && ~isempty(rawInfo.FirstDim)
    Info.FirstDim = rawInfo.FirstDim;
else
    Info.FirstDim = 'y';
end

if ~isempty(fieldnames(importedEntry)) && isfield(importedEntry, 'ExposureMsec') && ...
        ~isempty(importedEntry.ExposureMsec)
    Info.ExposureMsec = importedEntry.ExposureMsec;
elseif isfield(rawInfo, 'ExposureMsec') && ~isempty(rawInfo.ExposureMsec)
    Info.ExposureMsec = rawInfo.ExposureMsec;
end

if ~isempty(fieldnames(importedEntry)) && isfield(importedEntry, 'CamIdx') && ...
        ~isempty(importedEntry.CamIdx)
    Info.CamIdx = importedEntry.CamIdx;
elseif isfield(rawInfo, 'CamIdx') && ~isempty(rawInfo.CamIdx)
    Info.CamIdx = rawInfo.CamIdx;
end

% Acquisition-wide dual-camera flag. Session-level fields are otherwise
% deliberately excluded from this flat, file-facing Info, but MultiCam is
% needed by callers that decide whether dual-camera coregistration applies
% to the current file (e.g. DataViewer's currentDatSourceIsMultiCam).
if isfield(rawInfo, 'MultiCam') && ~isempty(rawInfo.MultiCam)
    Info.MultiCam = rawInfo.MultiCam;
end

if hasLegacySidecar
    Info.MetadataSource = 'legacy_sidecar';
else
    Info.MetadataSource = 'AcqInfos.mat';
end

if isfield(rawInfo, 'TimelineSource') && ~isempty(rawInfo.TimelineSource)
    Info.TimelineSource = rawInfo.TimelineSource;
end

if isfield(rawInfo, 'TimelineSourceIndex') && ~isempty(rawInfo.TimelineSourceIndex)
    Info.TimelineSourceIndex = rawInfo.TimelineSourceIndex;
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

function [Info, updatedFields] = iUpdateLegacyCoreFieldsForForwardCompatibility(Info)
%IUPDATELEGACYCOREFIELDSFORFORWARDCOMPATIBILITY
% Refresh core fields from legacy metadata using:
%   - FrameRateHz = Freq
%   - Width, Height, Length from datCat = [datSize datLength] indexed by dim_names

updatedFields = {};

if ~isfield(Info, 'dim_names') || isempty(Info.dim_names)
    return
end

if ~isfield(Info, 'datSize') || isempty(Info.datSize) || ...
        ~isfield(Info, 'datLength') || isempty(Info.datLength)
    return
end

dimNames = cellstr(string(Info.dim_names));
datSize = double(Info.datSize(:).');
datLength = double(Info.datLength);

datCat = [datSize, datLength];

if numel(datCat) ~= numel(dimNames)
    error('Umitoolbox:loadMetaData:invalidLegacyMetadata', ...
        ['Legacy metadata are inconsistent: [datSize datLength] must have ' ...
         'the same number of elements as dim_names.']);
end

% -------------------------------------------------------------------------
% FrameRateHz <- Freq
% -------------------------------------------------------------------------
if isfield(Info, 'Freq') && ~isempty(Info.Freq)
    newFrameRateHz = double(Info.Freq);

    if ~isfield(Info, 'FrameRateHz') || isempty(Info.FrameRateHz) || ...
            ~isequal(double(Info.FrameRateHz), newFrameRateHz)
        Info.FrameRateHz = newFrameRateHz;
        updatedFields{end+1} = 'FrameRateHz'; %#ok<AGROW>
    end
end

if isfield(Info, 'FrameRateHz') && ~isempty(Info.FrameRateHz)
    Info.Freq = double(Info.FrameRateHz);
end

% -------------------------------------------------------------------------
% Height <- Y, Width <- X, Length <- T
% -------------------------------------------------------------------------
idxY = find(strcmp(dimNames, 'Y'), 1, 'first');
idxX = find(strcmp(dimNames, 'X'), 1, 'first');
idxT = find(strcmp(dimNames, 'T'), 1, 'first');

if ~isempty(idxY)
    newHeight = datCat(idxY);
    if ~isfield(Info, 'Height') || isempty(Info.Height) || ...
            ~isequal(double(Info.Height), newHeight)
        Info.Height = newHeight;
        updatedFields{end+1} = 'Height'; %#ok<AGROW>
    end
end

if ~isempty(idxX)
    newWidth = datCat(idxX);
    if ~isfield(Info, 'Width') || isempty(Info.Width) || ...
            ~isequal(double(Info.Width), newWidth)
        Info.Width = newWidth;
        updatedFields{end+1} = 'Width'; %#ok<AGROW>
    end
end

if ~isempty(idxT)
    newLength = datCat(idxT);
    if ~isfield(Info, 'Length') || isempty(Info.Length) || ...
            ~isequal(double(Info.Length), newLength)
        Info.Length = newLength;
        updatedFields{end+1} = 'Length'; %#ok<AGROW>
    end
end

end
