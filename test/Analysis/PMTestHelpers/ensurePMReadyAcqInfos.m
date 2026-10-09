function AcqInfoStream = ensurePMReadyAcqInfos(saveFolder, datFiles, varargin)
%ENSUREPMREADYACQINFOS Register .dat files so PipelineManager accepts a fixture.
%
%   AcqInfoStream = ensurePMReadyAcqInfos(saveFolder)
%   AcqInfoStream = ensurePMReadyAcqInfos(saveFolder, datFiles)
%   AcqInfoStream = ensurePMReadyAcqInfos(..., 'FrameRateHz', 10)
%
%   PipelineManager.executePipeline refuses any SaveFolder that
%   ISLEGACYSCHEMAFOLDER reports as legacy, i.e. whose AcqInfos.mat has no
%   AcqInfoStream.ImportedChannels registry and whose channels cannot be
%   inferred by resolveImportedChannelFallback. Many per-function Analysis
%   fixtures were written for direct function calls, where that registry is
%   never consulted, so they are legacy-schema as far as PipelineManager is
%   concerned and cannot be driven through the pipeline at all.
%
%   ENSUREPMREADYACQINFOS upgrades such a fixture in place by appending an
%   ImportedChannels entry for each .dat file, deriving each channel's Length
%   from the file size and the fixture's own frame geometry. It is a
%   test-fixture utility only: it makes an existing fixture representative of
%   a current-schema dataset. It does not invent acquisition metadata beyond
%   the per-channel registry entries.
%
%   DATFILES is a char, string array or cellstr of .dat file names inside
%   SAVEFOLDER. When omitted, every .dat file in SAVEFOLDER is registered.
%
%   Name-Value options:
%       'FrameRateHz' - Frame rate recorded for each channel. Default: the
%                       fixture's AcqInfoStream.FrameRateHz.
%       'Datatype'    - Element type used to derive Length from file bytes.
%                       Default: the fixture's AcqInfoStream.Datatype, or
%                       'single'.
%
%   See also buildPMForScenario, appendImportedChannelInfo, isLegacySchemaFolder.

p = inputParser;
addRequired(p, 'saveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addOptional(p, 'datFiles', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'FrameRateHz', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'Datatype', '', @(x) ischar(x) || isstring(x));
parse(p, saveFolder, datFiles, varargin{:});

saveFolder = char(string(p.Results.saveFolder));
acqInfoPath = fullfile(saveFolder, 'AcqInfos.mat');

assert(isfile(acqInfoPath), ...
    'Umitoolbox:ensurePMReadyAcqInfos:missingAcqInfos', ...
    'AcqInfos.mat not found in "%s".', saveFolder);

loaded = load(acqInfoPath, 'AcqInfoStream');
assert(isfield(loaded, 'AcqInfoStream'), ...
    'Umitoolbox:ensurePMReadyAcqInfos:invalidAcqInfos', ...
    'File "%s" does not contain AcqInfoStream.', acqInfoPath);
AcqInfoStream = loaded.AcqInfoStream;

% -------------------------------------------------------------------------
% Resolve the .dat files to register
% -------------------------------------------------------------------------
requested = p.Results.datFiles;
if isempty(requested)
    found = dir(fullfile(saveFolder, '*.dat'));
    requested = {found.name};
else
    requested = cellstr(string(requested));
end

assert(~isempty(requested), ...
    'Umitoolbox:ensurePMReadyAcqInfos:noDatFiles', ...
    'No .dat files to register in "%s".', saveFolder);

% -------------------------------------------------------------------------
% Resolve geometry used to derive each channel's Length
% -------------------------------------------------------------------------
frameRate = p.Results.FrameRateHz;
if isempty(frameRate)
    assert(isfield(AcqInfoStream, 'FrameRateHz'), ...
        'Umitoolbox:ensurePMReadyAcqInfos:missingFrameRate', ...
        'AcqInfoStream has no FrameRateHz; pass ''FrameRateHz'' explicitly.');
    frameRate = AcqInfoStream.FrameRateHz;
end

datatype = char(string(p.Results.Datatype));
if isempty(datatype)
    if isfield(AcqInfoStream, 'Datatype') && ~isempty(AcqInfoStream.Datatype)
        datatype = char(string(AcqInfoStream.Datatype));
    else
        datatype = 'single';
    end
end

bytesPerElem = iBytesPerElement(datatype);

assert(isfield(AcqInfoStream, 'Height') && isfield(AcqInfoStream, 'Width'), ...
    'Umitoolbox:ensurePMReadyAcqInfos:missingGeometry', ...
    'AcqInfoStream must define Height and Width.');
pixelsPerFrame = double(AcqInfoStream.Height) * double(AcqInfoStream.Width);

% -------------------------------------------------------------------------
% Register each channel
% -------------------------------------------------------------------------
for iFile = 1:numel(requested)
    datName = char(string(requested{iFile}));
    datPath = fullfile(saveFolder, datName);

    fileInfo = dir(datPath);
    assert(~isempty(fileInfo), ...
        'Umitoolbox:ensurePMReadyAcqInfos:datFileNotFound', ...
        'Data file "%s" not found.', datPath);

    if isDatWithHeader(datPath)
        % Headered file (.dat header): its length is in the header.
        hdr = readDatHeader(datPath);
        nFrames = hdr.dimSizes(strcmp(hdr.dimNames, 'T'));
    else
        nFrames = fileInfo.bytes / (pixelsPerFrame * bytesPerElem);
    end
    assert(nFrames > 0 && nFrames == round(nFrames), ...
        'Umitoolbox:ensurePMReadyAcqInfos:frameCountMismatch', ...
        ['Could not derive an integer frame count for "%s" ' ...
        '(bytes=%d, pixels/frame=%d, bytes/element=%d).'], ...
        datName, fileInfo.bytes, pixelsPerFrame, bytesPerElem);

    AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, struct( ...
        'DatFile', datName, ...
        'Length', nFrames, ...
        'FrameRateHz', frameRate), ...
        'Overwrite', true);
end

save(acqInfoPath, 'AcqInfoStream');
end

% =========================================================================
% Local helpers
% =========================================================================
function nBytes = iBytesPerElement(datatype)
%IBYTESPERELEMENT Bytes per element for a .dat element type.

switch lower(strtrim(datatype))
    case {'uint8', 'int8'}
        nBytes = 1;
    case {'uint16', 'int16'}
        nBytes = 2;
    case {'single', 'uint32', 'int32'}
        nBytes = 4;
    case {'double', 'uint64', 'int64'}
        nBytes = 8;
    otherwise
        error('Umitoolbox:ensurePMReadyAcqInfos:unsupportedDatatype', ...
            'Unsupported .dat element type "%s".', datatype);
end
end
