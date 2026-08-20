function [resolvedAcqInfo, report] = resolveImportedChannelFallback( ...
        AcqInfoStream, saveFolder, varargin)
%RESOLVEIMPORTEDCHANNELFALLBACK Infer legacy imported-channel metadata.
%
%   [acqInfo, report] = resolveImportedChannelFallback(AcqInfoStream, saveFolder)
%   [...] = resolveImportedChannelFallback(..., 'RequiredChannels', {'red','green'})
%
%   Returns an in-memory AcqInfoStream with ImportedChannels only when that
%   field is absent. Explicit ImportedChannels metadata is authoritative and
%   is returned without normalization or replacement. For older classified
%   acquisitions, the fallback resolves the red, green, and yellow channel
%   conventions from, in precedence order: per-channel legacy .mat timing,
%   legacy IlluminationN declarations, and canonical channel filenames.
%   The .dat file size remains authoritative for each inferred Length.
%
%   This helper never writes AcqInfos.mat or any sidecar metadata. When
%   RequiredChannels is supplied, incomplete or conflicting evidence raises
%   a diagnostic suitable for PipelineManager preflight.

p = inputParser;
p.FunctionName = 'resolveImportedChannelFallback';
addRequired(p, 'AcqInfoStream', @(x) isstruct(x) && isscalar(x));
addRequired(p, 'saveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'RequiredChannels', {}, @(x) iscell(x) || isstring(x) || ischar(x));
parse(p, AcqInfoStream, saveFolder, varargin{:});

resolvedAcqInfo = AcqInfoStream;
report = struct('Source', 'none', 'InferredChannels', {{}}, ...
    'Messages', {{}});

% An explicitly present field, including an explicitly empty struct array,
% is the acquisition owner's decision and must not be replaced.
if isfield(AcqInfoStream, 'ImportedChannels')
    report.Source = 'explicit';
    return
end

saveFolder = char(string(saveFolder));
required = iNormalizeRequiredChannels(p.Results.RequiredChannels);
legacy = iReadLegacyIlluminations(AcqInfoStream);

if isempty(required)
    required = unique({legacy.Name}, 'stable');
    if isempty(required)
        required = iCanonicalFilesInFolder(saveFolder);
    end
end

inferred = struct('DatFile', {}, 'Tag', {}, 'Color', {}, 'Length', {}, ...
    'FrameRateHz', {}, 'ExposureMsec', {}, 'CamIdx', {});
missing = {};

for iChannel = 1:numel(required)
    channelName = required{iChannel};
    legacyEntries = legacy(strcmp({legacy.Name}, channelName));
    entry = iResolveOneChannel(AcqInfoStream, saveFolder, channelName, legacyEntries);
    if isempty(entry)
        missing{end+1} = channelName; %#ok<AGROW>
    else
        inferred(end+1) = entry; %#ok<AGROW>
    end
end

if isempty(inferred)
    return
end

if ~isempty(missing)
    error('Umitoolbox:resolveImportedChannelFallback:incompleteEvidence', ...
        ['Cannot infer ImportedChannels for %s. Missing usable legacy ' ...
         'IlluminationN metadata or canonical .dat files in "%s".'], ...
        strjoin(missing, ', '), saveFolder);
end

resolvedAcqInfo.ImportedChannels = inferred;
report.Source = 'fallback';
report.InferredChannels = {inferred.DatFile};
report.Messages = {'ImportedChannels inferred in memory; AcqInfos.mat was not modified.'};

end

% =========================================================================
% Local helpers
% =========================================================================
function required = iNormalizeRequiredChannels(value)

if ischar(value) || (isstring(value) && isscalar(value))
    value = {char(string(value))};
elseif isstring(value)
    value = cellstr(value(:).');
end

required = cell(1, numel(value));
for iValue = 1:numel(value)
    required{iValue} = iNormalizeChannelName(value{iValue});
end
required = unique(required, 'stable');

end

function names = iCanonicalFilesInFolder(folderPath)

names = {};
for name = {'red', 'green', 'yellow'}
    if isfile(fullfile(folderPath, [name{1} '.dat']))
        names{end+1} = name{1}; %#ok<AGROW>
    end
end

end

function legacy = iReadLegacyIlluminations(acqInfo)

legacy = struct('Name', {}, 'DatFile', {}, 'FrameRateHz', {}, ...
    'ExposureMsec', {}, 'CamIdx', {});
fields = fieldnames(acqInfo);
matches = regexp(fields, '^Illumination(\d+)$', 'tokens', 'once');
idx = find(~cellfun(@isempty, matches));
if isempty(idx)
    return
end

sequence = cellfun(@(x) str2double(x{1}), matches(idx));
[~, order] = sort(sequence);
idx = idx(order);
legacy = repmat(struct('Name', '', 'DatFile', '', 'FrameRateHz', nan, ...
    'ExposureMsec', nan, 'CamIdx', nan), 1, numel(idx));
count = 0;

for iField = 1:numel(idx)
    raw = acqInfo.(fields{idx(iField)});
    if ~isstruct(raw) || ~isscalar(raw)
        continue
    end
    sourceName = iFirstTextField(raw, {'Color', 'Tag', 'Name'});
    if isempty(sourceName)
        continue
    end
    try
        name = iNormalizeChannelName(sourceName);
    catch
        continue
    end

    datFile = iFirstTextField(raw, {'DatFile', 'FileName', 'Filename'});
    if isempty(datFile)
        datFile = [name '.dat'];
    else
        [~, base, ext] = fileparts(datFile);
        if isempty(ext)
            ext = '.dat';
        end
        datFile = [base ext];
    end

    count = count + 1;
    legacy(count) = struct( ...
        'Name', name, ...
        'DatFile', datFile, ...
        'FrameRateHz', iFirstNumericField(raw, {'FrameRateHz', 'Freq'}), ...
        'ExposureMsec', iFirstNumericField(raw, {'ExposureMsec', 'Exposure'}), ...
        'CamIdx', iFirstNumericField(raw, {'CamIdx', 'CameraIndex'}));
end
legacy = legacy(1:count);

end

function entry = iResolveOneChannel(acqInfo, folderPath, channelName, legacy)

entry = struct([]);
[datFile, camIdx, color, legacyFreq, exposure] = iResolveLegacyIdentity( ...
    acqInfo, channelName, legacy);
datPath = fullfile(folderPath, datFile);
if ~isfile(datPath)
    return
end

[length, sidecarFreq] = iResolveLegacyDatMetadata(datPath, acqInfo);
if ~isempty(sidecarFreq)
    frameRateHz = sidecarFreq;
elseif ~isempty(legacyFreq)
    frameRateHz = legacyFreq;
elseif isfield(acqInfo, 'FrameRateHz') && ~isempty(acqInfo.FrameRateHz)
    frameRateHz = double(acqInfo.FrameRateHz);
else
    error('Umitoolbox:resolveImportedChannelFallback:missingFrameRate', ...
        ['Cannot infer FrameRateHz for "%s". Provide a legacy channel ' ...
         'sidecar, IlluminationN.FrameRateHz, or AcqInfoStream.FrameRateHz.'], datFile);
end

validateattributes(frameRateHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    'resolveImportedChannelFallback', 'FrameRateHz');

entry = struct( ...
    'DatFile', datFile, ...
    'Tag', channelName, ...
    'Color', color, ...
    'Length', length, ...
    'FrameRateHz', double(frameRateHz), ...
    'ExposureMsec', exposure, ...
    'CamIdx', camIdx);

end

function [datFile, camIdx, color, frameRateHz, exposure] = iResolveLegacyIdentity( ...
        acqInfo, channelName, legacy)

datFile = [channelName '.dat'];
color = channelName;
frameRateHz = [];
exposure = [];

if isempty(legacy)
    camIdx = iDefaultCamIdx(acqInfo, channelName);
    return
end

datFiles = unique({legacy.DatFile}, 'stable');
if numel(datFiles) ~= 1
    error('Umitoolbox:resolveImportedChannelFallback:conflictingEvidence', ...
        'Legacy illumination records for "%s" map to multiple .dat files.', channelName);
end
datFile = datFiles{1};

camValues = [legacy.CamIdx];
camValues = camValues(~isnan(camValues));
if isempty(camValues)
    camIdx = iDefaultCamIdx(acqInfo, channelName);
elseif numel(unique(camValues)) ~= 1
    error('Umitoolbox:resolveImportedChannelFallback:conflictingEvidence', ...
        'Legacy illumination records for "%s" specify conflicting CamIdx values.', channelName);
else
    camIdx = double(camValues(1));
end

freqValues = [legacy.FrameRateHz];
freqValues = freqValues(~isnan(freqValues));
if ~isempty(freqValues)
    if numel(unique(freqValues)) ~= 1
        error('Umitoolbox:resolveImportedChannelFallback:conflictingEvidence', ...
            'Legacy illumination records for "%s" specify conflicting FrameRateHz values.', channelName);
    end
    frameRateHz = double(freqValues(1));
elseif isfield(acqInfo, 'FrameRateHz') && ~isempty(acqInfo.FrameRateHz)
    % Repeated legacy illumination records are combined into one canonical
    % .dat file by ImagesClassification, so its native rate scales likewise.
    frameRateHz = double(acqInfo.FrameRateHz) * numel(legacy);
end

exposureValues = [legacy.ExposureMsec];
exposureValues = exposureValues(~isnan(exposureValues));
if ~isempty(exposureValues)
    if numel(unique(exposureValues)) ~= 1
        error('Umitoolbox:resolveImportedChannelFallback:conflictingEvidence', ...
            'Legacy illumination records for "%s" specify conflicting exposure values.', channelName);
    end
    exposure = double(exposureValues(1));
elseif isfield(acqInfo, 'ExposureMsec') && ~isempty(acqInfo.ExposureMsec)
    exposure = double(acqInfo.ExposureMsec);
end

end

function camIdx = iDefaultCamIdx(acqInfo, channelName)

if isfield(acqInfo, 'MultiCam') && logical(acqInfo.MultiCam)
    error('Umitoolbox:resolveImportedChannelFallback:missingCameraIndex', ...
        ['Cannot infer CamIdx for "%s" in a multicamera acquisition. ' ...
         'Provide legacy IlluminationN.CamIdx or ImportedChannels metadata.'], channelName);
end
camIdx = 1;

end

function length = iInferDatLength(datPath, acqInfo)

if ~isfield(acqInfo, 'Height') || ~isfield(acqInfo, 'Width')
    error('Umitoolbox:resolveImportedChannelFallback:missingFrameSize', ...
        'AcqInfoStream must contain Height and Width to infer .dat channel lengths.');
end
height = double(acqInfo.Height);
width = double(acqInfo.Width);
validateattributes(height, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'});
validateattributes(width, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'});

datatype = 'single';
if isfield(acqInfo, 'Datatype') && ~isempty(acqInfo.Datatype)
    datatype = char(string(acqInfo.Datatype));
end
bytes = iBytesPerElement(datatype);
fileInfo = dir(datPath);
frameBytes = height * width * bytes;
if mod(fileInfo.bytes, frameBytes) ~= 0
    error('Umitoolbox:resolveImportedChannelFallback:invalidFileLength', ...
        'File size for "%s" is incompatible with AcqInfoStream frame dimensions.', datPath);
end
length = fileInfo.bytes / frameBytes;
if length < 1
    error('Umitoolbox:resolveImportedChannelFallback:emptyChannelFile', ...
        'Channel file "%s" does not contain any complete frames.', datPath);
end

if ~isfield(acqInfo, 'Length') || isempty(acqInfo.Length)
    error('Umitoolbox:resolveImportedChannelFallback:missingExpectedLength', ...
        ['Cannot safely infer the frame geometry for legacy channel "%s" without ' ...
         'a per-channel metadata sidecar or AcqInfoStream.Length.'], datPath);
end

expectedLength = double(acqInfo.Length);
validateattributes(expectedLength, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', 'integer'});
if length ~= expectedLength
    error('Umitoolbox:resolveImportedChannelFallback:ambiguousLegacyGeometry', ...
        ['The frame count inferred for legacy channel "%s" from AcqInfoStream ' ...
         'Height/Width is %g, but AcqInfoStream.Length is %g. The acquisition ' ...
         'may have been spatially binned or cropped. Restore its per-channel ' ...
         'metadata sidecar or reprocess it from raw data.'], ...
        datPath, length, expectedLength);
end

end

function [length, frameRateHz] = iResolveLegacyDatMetadata(datPath, acqInfo)

[folderPath, baseName] = fileparts(datPath);
hasSidecar = isfile(fullfile(folderPath, [baseName '.mat'])) || ...
    isfile(fullfile(folderPath, [baseName '_info.mat']));

if hasSidecar
    info = loadMetaData(datPath);
    if strcmp(info.MetadataSource, 'legacy_sidecar')
        length = double(info.Length);
        frameRateHz = double(info.FrameRateHz);
        return
    end
end

length = iInferDatLength(datPath, acqInfo);
frameRateHz = [];

end

function value = iFirstTextField(S, names)

value = '';
for iName = 1:numel(names)
    if isfield(S, names{iName}) && ~isempty(S.(names{iName}))
        value = char(string(S.(names{iName})));
        return
    end
end

end

function value = iFirstNumericField(S, names)

value = nan;
for iName = 1:numel(names)
    if isfield(S, names{iName}) && ~isempty(S.(names{iName}))
        raw = S.(names{iName});
        if isnumeric(raw) && isscalar(raw) && isfinite(raw)
            value = double(raw);
            return
        end
    end
end

end

function name = iNormalizeChannelName(value)

name = lower(strtrim(char(string(value))));
if contains(name, 'red')
    name = 'red';
elseif contains(name, 'green')
    name = 'green';
elseif contains(name, 'yellow') || contains(name, 'amber')
    name = 'yellow';
else
    error('Umitoolbox:resolveImportedChannelFallback:unsupportedChannel', ...
        'Only red, green, yellow, and amber legacy channel names are supported.');
end

end

function nBytes = iBytesPerElement(datatype)

switch lower(char(string(datatype)))
    case {'single', 'float32'}
        nBytes = 4;
    case {'double', 'float64'}
        nBytes = 8;
    case {'uint8', 'int8', 'logical'}
        nBytes = 1;
    case {'uint16', 'int16'}
        nBytes = 2;
    case {'uint32', 'int32'}
        nBytes = 4;
    case {'uint64', 'int64'}
        nBytes = 8;
    otherwise
        error('Umitoolbox:resolveImportedChannelFallback:unsupportedDatatype', ...
            'Unsupported AcqInfoStream.Datatype: "%s".', char(string(datatype)));
end

end
