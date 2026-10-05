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
%   Supported sources, in detection order for .dat files:
%       1) .dat files with a self-describing header (see readDatHeader).
%          Described from the header alone; AcqInfos.mat and sidecars in
%          the folder are ignored.
%       2) Legacy .dat files with sidecar metadata (.mat or _info.mat).
%       3) .umt files: embedded metadata, then the first entry (dimNames,
%          sizes, meta.FrameRateHz); AcqInfos.mat is not used
%
%   Headerless .dat files without a legacy sidecar (formerly described
%   only by the folder-global AcqInfos.mat) are no longer supported and
%   raise Umitoolbox:loadMetaData:acqInfosBoundUnsupported: re-import the
%   raw data to obtain headered files.
%
%   .dat Info schema (all .dat kinds):
%       filePath      - Full path of the .dat file
%       format        - 'header' or 'legacySidecar'
%       dataOffset    - Byte offset of the data (512 or 0)
%       dataClass     - MATLAB class of the stored values
%       dimNames      - Axis names in memory order, e.g. {'Y','X','T'}
%       dimSizes      - Full shape in dimNames order
%       frameRateHz   - Frame rate; NaN without a T axis
%       exposureMsec  - Exposure in ms; NaN when unknown
%       channelName   - Channel label ('' for headerless files)
%       writeComplete - Header write-complete flag (true for headerless)
%   Use datAxisSize(Info, axisName) for single axis sizes.
%
%   For .dat files, Info holds exactly these schema fields. The former
%   names (Height, Width, Length, datLength, datSize, dim_names, Datatype,
%   Freq, FrameRateHz, ExposureMsec, ...) were removed in .dat header
%   Phase 8a; see the field mapping in the custom-function update guide.
%
%   For headered files whose write-complete flag is not set, a warning
%   (Umitoolbox:loadMetaData:incompleteFile) is issued once. A headered
%   file marked complete must have exactly the size its header describes.
%
%   Legacy sidecar .dat files:
%       - The sidecar metadata (.mat or _info.mat) are read with their old
%         field names and converted to the schema. Fields missing from the
%         sidecar are completed from AcqInfos.mat, as before.
%       - Height and Width come from the sidecar; the temporal length is
%         inferred from the actual file size.
%       - Valid legacy event-split metadata are preserved and are not
%         collapsed to continuous YXT.
%       - Files are single precision unless the sidecar defines Datatype.
%
%   Notes:
%       - For .umt files, Info keeps the UMT metadata names (for example
%         FrameRateHz, Freq, dim_names). Embedded metadata is optional;
%         when missing, core metadata are derived from the first entry.

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
        if isDatWithHeader(fileName)
            Info = iLoadHeaderedDatMetaData(fileName);
        else
            [resolved, ownExposureMsec] = iLoadDatMetaData(fileName, folderPath, baseName);
            Info = iDatSchemaFromResolved(resolved, ownExposureMsec);
        end

    case '.umt'
        Info = iLoadUMTMetaData(fileName, folderPath);
end

end

% =========================================================================
% .dat Info schema
% =========================================================================
function Info = iLoadHeaderedDatMetaData(fileName)
%ILOADHEADEREDDATMETADATA Describe a headered .dat file from its header alone.

hdr = readDatHeader(fileName);
fileInfo = dir(fileName);
expectedBytes = hdr.dataOffset + hdr.expectedDataBytes;

if hdr.writeComplete
    validateDatHeader(hdr, 'FileBytes', fileInfo.bytes);
elseif fileInfo.bytes ~= expectedBytes
    warning('Umitoolbox:loadMetaData:incompleteFile', ...
        ['"%s" is not marked write-complete and may be incomplete: the file has ' ...
         '%d bytes, its header describes %d bytes.'], fileName, fileInfo.bytes, expectedBytes);
else
    warning('Umitoolbox:loadMetaData:incompleteFile', ...
        '"%s" is not marked write-complete and may be incomplete.', fileName);
end

Info = struct();
Info.filePath = fileName;
Info.format = 'header';
Info.dataOffset = hdr.dataOffset;
Info.dataClass = hdr.dataClass;
Info.dimNames = hdr.dimNames;
Info.dimSizes = hdr.dimSizes;
Info.frameRateHz = hdr.frameRateHz;
Info.exposureMsec = hdr.exposureMsec;
Info.channelName = hdr.channelName;
Info.writeComplete = hdr.writeComplete;

end

function Info = iDatSchemaFromResolved(resolved, ownExposureMsec)
%IDATSCHEMAFROMRESOLVED Build the .dat Info schema from resolved legacy sidecar metadata.
%
% ownExposureMsec is this file's own exposure, resolved by iLoadDatMetaData
% (NaN when unknown).

Info = struct();
Info.filePath = resolved.datFile;
Info.format = 'legacySidecar';
Info.dataOffset = 0;
Info.dataClass = resolved.Datatype;
Info.dimNames = resolved.dim_names;

datSize = double(resolved.datSize(:).');
idxT = find(strcmp(resolved.dim_names, 'T'), 1, 'first');
if numel(datSize) == numel(resolved.dim_names)
    Info.dimSizes = datSize;
else
    % datSize holds only the non-T axes: insert the resolved length at T.
    Info.dimSizes = [datSize(1:idxT-1), double(resolved.Length), datSize(idxT:end)];
end

Info.frameRateHz = double(resolved.FrameRateHz);
Info.exposureMsec = double(ownExposureMsec);
Info.channelName = '';
Info.writeComplete = true;

end

% =========================================================================
% Local helpers
% =========================================================================
function [Info, ownExposureMsec] = iLoadDatMetaData(fileName, folderPath, baseName)
%ILOADDATMETADATA Build flat metadata for a headerless legacy sidecar .dat file.
%
% ownExposureMsec is the exposure of this file itself: the speckle exposure
% when there is file-specific evidence that the file holds speckle data,
% otherwise the general exposure; NaN when unknown.

legacyInfo = iLoadLegacySidecar(folderPath, baseName);
if isempty(fieldnames(legacyInfo))
    % Headerless without a sidecar: described only by AcqInfos.mat
    % (dev-era imports made before the headered importers). Unsupported.
    error('Umitoolbox:loadMetaData:acqInfosBoundUnsupported', ...
        ['"%s" has no header and no legacy metadata .mat file, so it can no ' ...
         'longer be read. Files of this kind were described only by ' ...
         'AcqInfos.mat and are not supported anymore. Re-import the raw data ' ...
         'with the umIToolbox/DataViewer importer, which writes headered .dat files.'], ...
        fileName);
end
acqInfo = iLoadAcqInfo(folderPath);

% Legacy metadata takes precedence. Append only missing fields from
% AcqInfoStream to preserve source-specific semantics.
Info = legacyInfo;

% A sidecar without dim_names describes a Y-X-T file. Default it before
% merging AcqInfos.mat so the refresh below takes Height, Width,
% Length, and FrameRateHz from the sidecar, which has precedence.
if (~isfield(Info, 'dim_names') || isempty(Info.dim_names)) && ...
        isfield(Info, 'datSize') && isfield(Info, 'datLength') && ...
        numel(Info.datSize) + numel(Info.datLength) == 3
    Info.dim_names = {'Y', 'X', 'T'};
end

Info = iAppendMissingFields(Info, acqInfo);

% Refresh forward-compatible core fields from the legacy payload.
[Info, updatedFields] = iUpdateLegacyCoreFieldsForForwardCompatibility(Info);

if ~isempty(updatedFields)
    fprintf(['loadMetaData: Updated legacy metadata field(s) for ' ...
        'forward compatibility in "%s": %s\n'], ...
        fileName, strjoin(updatedFields, ', '));
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
% strict on file integrity. Legacy sidecar metadata are preserved for
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
elseif iIsLegacyEventSplitDatSize(Info)
    error('Umitoolbox:loadMetaData:legacyEventSplitUnsupported', ...
        ['This file appears to be a legacy transformed/event-split dataset ' ...
         '(non-YXT, non-3D-image-time-series) created by an older version of ' ...
         'UMIT. This format is not currently supported for opening in this ' ...
         'version. There is currently no automated way to convert it to the ' ...
         'current format \x2014 reprocessing from the original raw/continuous ' ...
         'data is recommended if this analysis is needed. File: "%s".'], ...
        fileName);
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

% If datSize explicitly includes T, keep the multi-dimensional layout but
% update the T slot to the actual on-disk length.
if isfield(Info, 'datSize') && ~isempty(Info.datSize) && numel(Info.datSize) == numel(Info.dim_names)
    Info.datSize(idxT) = actualLength;
end

% Keep only metadata that directly describes this .dat file.
importedEntry = iFindImportedChannelForInfo(acqInfo, fileName, actualLength);
isOwnSpeckle = iIsOwnSpeckleFile(fileName, importedEntry, legacyInfo);
Info = iFinalizeDatInfo(Info, acqInfo, fileName, folderPath, actualLength, ...
    importedEntry);

if isOwnSpeckle && isfield(Info, 'ExposureSpeckleMsec') && ~isempty(Info.ExposureSpeckleMsec)
    ownExposureMsec = double(Info.ExposureSpeckleMsec);
elseif isfield(Info, 'ExposureMsec') && ~isempty(Info.ExposureMsec)
    ownExposureMsec = double(Info.ExposureMsec);
else
    ownExposureMsec = NaN;
end

end

function tf = iIsOwnSpeckleFile(fileName, importedEntry, legacyInfo)
%IISOWNSPECKLEFILE File-specific evidence that this .dat holds speckle data.
%
% Unlike iIsSpeckleDataFile (which also accepts folder-level AcqInfos.mat
% fields and so flags every file of an acquisition that has a speckle
% channel), only the file's name, its own ImportedChannels entry, and its
% own legacy sidecar count here.

[~, fileBase] = fileparts(fileName);
tf = iContainsSpeckleText(fileBase) || ...
    iContainsSpeckleText(importedEntry) || ...
    iHasNonEmptyField(importedEntry, 'ExposureSpeckleMsec') || ...
    iHasNonEmptyField(legacyInfo, 'ExposureSpeckleMsec');
end

function tf = iIsLegacyEventSplitDatSize(Info)
%IISLEGACYEVENTSPLITDATSIZE Detect the pre-Astrocyte-fix event-split layout.
%
% This legacy convention stores 4 dim_names (a non-YXT dimension, e.g. 'E',
% plus Y/X/T) with datSize and datLength each holding 2 elements -- a 2-and-2
% split that matches neither currently-supported datSize/dim_names
% convention. Detection is metadata-shape-only so it never over-matches
% currently-loadable 3D YXT or full/({numel(dim_names)-1})-length datSize
% files, which are already handled by the branches above this check.

tf = false;

if ~isfield(Info, 'dim_names') || numel(Info.dim_names) ~= 4
    return
end

if ~any(~ismember(Info.dim_names, {'Y', 'X', 'T'}))
    return
end

if ~isfield(Info, 'datSize') || ~isfield(Info, 'datLength') || ...
        isempty(Info.datSize) || isempty(Info.datLength)
    return
end

tf = numel(Info.datSize) == 2 && numel(Info.datLength) == 2;

end

function Info = iLoadUMTMetaData(fileName, folderPath)
%ILOADUMTMETADATA Build flat metadata for a .umt file.

umt = iLoadUMTFromFile(fileName);
validateUMTStruct(umt, 'requireEventInfo', false);

Info = struct();

% Embedded metadata, if any, comes first. AcqInfos.mat is not merged: it
% describes the raw acquisition, not this file's data (.dat header Phase
% 7a); the file's own entries describe the rest.
embedded = iExtractEmbeddedMetadata(umt);
Info = iAppendMissingFields(Info, embedded);

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
    if ~isempty(idxT)
        if ~isfield(Info, 'datLength') || isempty(Info.datLength)
            Info.datLength = dimSizes(idxT);
        end
        % Length mirrors this entry's own temporal extent (datLength).
        Info.Length = Info.datLength;
    end

    % Frame rate: the entry's own meta.FrameRateHz, when present.
    if ~isfield(Info, 'FrameRateHz') && ~isfield(Info, 'Freq') && ...
            isfield(firstEntry, 'meta') && isstruct(firstEntry.meta) && ...
            isfield(firstEntry.meta, 'FrameRateHz') && ~isempty(firstEntry.meta.FrameRateHz)
        Info.FrameRateHz = double(firstEntry.meta.FrameRateHz);
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

function Info = iFinalizeDatInfo(rawInfo, acqInfo, fileName, folderPath, actualLength, importedEntry)
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

exposureMsec = [];
if ~isempty(fieldnames(importedEntry)) && isfield(importedEntry, 'ExposureMsec') && ...
        ~isempty(importedEntry.ExposureMsec)
    exposureMsec = importedEntry.ExposureMsec;
elseif isfield(rawInfo, 'ExposureMsec') && ~isempty(rawInfo.ExposureMsec)
    exposureMsec = rawInfo.ExposureMsec;
end

if ~isempty(exposureMsec)
    Info.ExposureMsec = exposureMsec;
end

% Keep the speckle-specific exposure name available to legacy analysis code.
% ExposureMsec remains the canonical normalized field; this is only a
% file-facing compatibility alias for data identified as speckle.
if iIsSpeckleDataFile(fileName, importedEntry, rawInfo, acqInfo)
    speckleExposure = iFirstNonEmptyField(importedEntry, {'ExposureSpeckleMsec'});
    if isempty(speckleExposure)
        speckleExposure = iFirstNonEmptyField(rawInfo, {'ExposureSpeckleMsec'});
    end
    if isempty(speckleExposure)
        speckleExposure = exposureMsec;
    end

    if ~isempty(speckleExposure)
        Info.ExposureSpeckleMsec = speckleExposure;
    end
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

Info.MetadataSource = 'legacy_sidecar';

end

function value = iFirstNonEmptyField(S, fieldNames)
%IFIRSTNONEMPTYFIELD Return the first non-empty field found in a struct.

value = [];
if ~isstruct(S) || ~isscalar(S)
    return
end

for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    if isfield(S, fieldName) && ~isempty(S.(fieldName))
        value = S.(fieldName);
        return
    end
end

end

function tf = iIsSpeckleDataFile(fileName, importedEntry, rawInfo, acqInfo)
%IISSPECKLEDATAFILE Identify a file using speckle channel metadata.

[~, baseName] = fileparts(fileName);
tf = iContainsSpeckleText(baseName) || ...
    iContainsSpeckleText(importedEntry) || ...
    iHasNonEmptyField(importedEntry, 'ExposureSpeckleMsec') || ...
    iHasNonEmptyField(rawInfo, 'ExposureSpeckleMsec') || ...
    iHasNonEmptyField(acqInfo, 'ExposureSpeckleMsec');
if tf || ~isstruct(acqInfo) || ~isscalar(acqInfo)
    return
end

tf = iContainsSpeckleText(acqInfo);
if tf
    return
end

acqFields = fieldnames(acqInfo);
for iField = 1:numel(acqFields)
    candidate = acqInfo.(acqFields{iField});
    if isstruct(candidate) && isscalar(candidate) && iContainsSpeckleText(candidate)
        tf = true;
        return
    end
end

end

function tf = iHasNonEmptyField(S, fieldName)
%IHASNONEMPTYFIELD Check whether a scalar struct has a populated field.

tf = isstruct(S) && isscalar(S) && isfield(S, fieldName) && ...
    ~isempty(S.(fieldName));

end

function tf = iContainsSpeckleText(value)
%ICONTAINSSPECKLETEXT Check common channel labels for "speckle".

tf = false;
if ischar(value) || (isstring(value) && isscalar(value))
    tf = contains(char(string(value)), 'speckle', 'IgnoreCase', true);
    return
end

if ~isstruct(value) || ~isscalar(value)
    return
end

labelFields = {'Tag', 'Color', 'Name', 'datName'};
for iField = 1:numel(labelFields)
    fieldName = labelFields{iField};
    if isfield(value, fieldName) && ...
            (ischar(value.(fieldName)) || (isstring(value.(fieldName)) && isscalar(value.(fieldName)))) && ...
            contains(char(string(value.(fieldName))), 'speckle', 'IgnoreCase', true)
        tf = true;
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
