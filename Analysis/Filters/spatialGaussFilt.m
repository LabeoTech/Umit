function outData = spatialGaussFilt(data, SaveFolder, varargin)
%SPATIALGAUSSFILT Apply a spatial Gaussian filter to image data.
%
%   outData = spatialGaussFilt(data, SaveFolder)
%   outData = spatialGaussFilt(data, SaveFolder, 'Sigma', sigma)
%
%   This function applies a spatial Gaussian filter using IMGAUSSFILT.
%
%   Supported execution modes:
%       1) STANDARD MODE (in-memory)
%          - Triggered when "data" is a numeric array or a UMT struct
%       2) LOW-RAM MODE (file-backed)
%          - Triggered when "data" is a .dat filename
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X or Y x X x T
%       2) Filename to a .dat file storing Y x X x T data
%       3) UMT struct with image entries using dimensions:
%              {'Y','X'}
%              {'Y','X','T'}
%              {'Y','X','E'}
%              {'Y','X','T','E'}
%       4) Filename to a .umt file containing one UMT struct
%
%   Input/output behavior:
%       - If the input is a numeric array, the output is a numeric array.
%       - If the input is a .dat filename, the output is a .dat filename.
%       - If the input is a UMT struct, the output is a UMT struct.
%       - If the input is a .umt filename, the file is loaded in RAM and
%         the output is a UMT struct.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder used for file resolution and AcqInfos.mat lookup.
%
%   Name-Value parameters:
%       Sigma      - Positive scalar Gaussian sigma. Default: 1
%
%   Output:
%       outData    - Filtered output with the same representation type as
%                    the input.
%
%   Notes:
%       - Raw .dat files are assumed to store continuous YXT data in single
%         precision.
%       - UMT entries with an E dimension preserve shared top-level
%         eventInfo unchanged.
%       - NaN values are replaced by zero before filtering and restored
%         afterward, preserving the original algorithm behavior.
%       - NaN masks are tracked per frame: a pixel that is NaN in one frame
%         does not cause valid data in other frames to be discarded.

default_Output = 'spatialGaussFilt.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Sigma', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
Sigma = double(p.Results.Sigma);

if ~isfolder(SaveFolder)
    error('spatialGaussFilt:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Case 1: Numeric array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty'}, mfilename, 'data');

    if ~(ismatrix(data) || ndims(data) == 3)
        error('spatialGaussFilt:InvalidArrayInput', ...
            'Numeric input must be a YX image or a YXT image time series.');
    end

    outData = iSpatialGaussBlock(data, Sigma);
    return
end

% -------------------------------------------------------------------------
% Case 2: File input
% -------------------------------------------------------------------------
if ischar(data) || (isstring(data) && isscalar(data))

    dataFile = char(string(data));

    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('spatialGaussFilt:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);

    switch ext
        case '.dat'
            outData = iSpatialGaussDatFile(dataFile, SaveFolder, Sigma, default_Output);
            return

        case '.umt'
            warning('spatialGaussFilt:UMTFileLoadsInRAM', ...
                ['RAM-safe mode is not available for data stored in this format. ' ...
                 'Loading the UMT content into RAM.']);
            data = loadData(dataFile);

        otherwise
            error('spatialGaussFilt:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

% -------------------------------------------------------------------------
% Case 3: UMT struct
% -------------------------------------------------------------------------
if ~isstruct(data)
    error('spatialGaussFilt:UnsupportedInputType', ...
        ['Input "data" must be a YX/YXT array, a .dat filename, ' ...
         'a UMT struct, or a .umt filename containing a UMT struct.']);
end

[entryNames, entryData, entryDims, labels, sourceEventInfo, hasE, entryMetas] = ...
    iExtractValidUMTData(data);

outStruct = data;
outStruct.data = struct();

for iEntry = 1:numel(entryNames)

    value = entryData{iEntry};
    dimNames = entryDims{iEntry};

    switch strjoin(dimNames, '')
        case 'YX'
            filtData = iSpatialGaussBlock(value, Sigma);

        case 'YXT'
            filtData = iSpatialGaussBlock(value, Sigma);

        case 'YXE'
            filtData = zeros(size(value), 'like', value);
            for iEvent = 1:size(value, 3)
                filtData(:,:,iEvent) = iSpatialGaussBlock(value(:,:,iEvent), Sigma);
            end

        case 'YXTE'
            filtData = zeros(size(value), 'like', value);
            for iEvent = 1:size(value, 4)
                filtData(:,:,:,iEvent) = iSpatialGaussBlock(value(:,:,:,iEvent), Sigma);
            end

        otherwise
            error('spatialGaussFilt:InvalidUMTEntryDims', ...
                'Unsupported dimNames in entry "%s".', entryNames{iEntry});
    end

    outStruct = genUMTStruct( ...
        outStruct, ...
        'value', filtData, ...
        'entryName', entryNames{iEntry}, ...
        'dimNames', dimNames, ...
        'meta', entryMetas{iEntry}, ...
        'overwrite', true);
end

if ~isempty(fieldnames(labels))
    outStruct.labels = labels;
elseif isfield(outStruct, 'labels')
    outStruct = rmfield(outStruct, 'labels');
end

if any(hasE)
    outStruct = appendUMTEventInfo(outStruct, ...
        'eventID', sourceEventInfo.eventID, ...
        'repetitionIndex', sourceEventInfo.repetitionIndex, ...
        'eventName', sourceEventInfo.eventName, ...
        'eventAxisMode', sourceEventInfo.eventAxisMode, ...
        'overwrite', true);
else
    if isfield(outStruct, 'eventInfo')
        outStruct = rmfield(outStruct, 'eventInfo');
    end
    validateUMTStruct(outStruct, 'requireEventInfo', true);
end

outData = outStruct;

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Apply a spatial Gaussian filter to image data.');

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'Image','ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YX/YXT array, .dat filename, ' ...
             'UMT struct, or .umt file containing one UMT struct.'], ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder used for file resolution and AcqInfos.mat lookup.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'Sigma', ...
            'parameter', ...
            'Positive scalar Gaussian sigma.', ...
            'kind', 'parameter', ...
            'default', 1, ...
            'allowed', [eps, Inf], ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'Image','ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Spatially filtered output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Local helper: spatial filtering with NaN restore
% =========================================================================
function outBlock = iSpatialGaussBlock(inBlock, sigma)
%ISPATIALGAUSSBLOCK Apply imgaussfilt while preserving spatial NaN regions.

if ismatrix(inBlock)
    spatialMask = isnan(inBlock);

    if any(spatialMask(:))
        work = inBlock;
        work(spatialMask) = 0;
        outBlock = imgaussfilt(work, sigma, 'FilterDomain', 'spatial');
        outBlock(spatialMask) = NaN;
    else
        outBlock = imgaussfilt(inBlock, sigma, 'FilterDomain', 'spatial');
    end

elseif ndims(inBlock) == 3
    spatialMask = isnan(inBlock);

    if any(spatialMask(:))
        work = inBlock;
        work(spatialMask) = 0;
        outBlock = imgaussfilt(work, sigma, 'FilterDomain', 'spatial');
        outBlock(spatialMask) = NaN;
    else
        outBlock = imgaussfilt(inBlock, sigma, 'FilterDomain', 'spatial');
    end

else
    error('spatialGaussFilt:InvalidBlockDims', ...
        'iSpatialGaussBlock expects a 2-D or 3-D block.');
end
end

% =========================================================================
% Local helper: low-RAM .dat execution
% =========================================================================
function outFile = iSpatialGaussDatFile(inFile, SaveFolder, sigma, defaultOutput)
%ISPATIALGAUSSDATFILE Apply spatial filtering to a raw YXT .dat file.

[Ny, Nx, Nt] = iGetRawDatInfo(SaveFolder, inFile);

% Write through a scratch file so the declared pipeline output only appears
% once the run has completed, and so the input can safely be the file the
% declared output would overwrite.
outFile = fullfile(SaveFolder, defaultOutput);
[~, outStem, outExt] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [outStem, '_writing', outExt]);
preallocateDatFile(tmpFile, [Ny, Nx, Nt], 'single');

frameBytes = Ny * Nx * getByteSize('single');
totalBytes = frameBytes * Nt;
nChunks = calculateMaxChunkSize(totalBytes, 2, 0.1);
chunkFrames = ceil(Nt / nChunks);
nChunks = ceil(Nt / chunkFrames);

fidIn  = fopen(inFile, 'r');
assert(fidIn ~= -1, 'spatialGaussFilt:OpenInputFailed', ...
    'Failed to open input file "%s".', inFile);
cIn = onCleanup(@() safeFclose(fidIn)); %#ok<NASGU>

fidOut = fopen(tmpFile, 'r+');
assert(fidOut ~= -1, 'spatialGaussFilt:OpenOutputFailed', ...
    'Failed to open output file "%s".', tmpFile);
cOut = onCleanup(@() safeFclose(fidOut)); %#ok<NASGU>

for c = 1:nChunks
    tStart = (c-1) * chunkFrames + 1;
    tEnd   = min(tStart + chunkFrames - 1, Nt);
    nThisChunk = tEnd - tStart + 1;

    fseek(fidIn, (tStart-1) * frameBytes, 'bof');
    slab = fread(fidIn, [Nx*Ny, nThisChunk], '*single');
    slab = reshape(slab, Ny, Nx, nThisChunk);

    % NaN masking is per chunk (per frame), not collapsed across the whole
    % file, so the result does not depend on chunk boundaries.
    spatialMask = isnan(slab);
    if any(spatialMask(:))
        slab = iApplyMaskedGauss(slab, sigma, spatialMask);
    else
        slab = imgaussfilt(slab, sigma, 'FilterDomain', 'spatial');
    end

    fseek(fidOut, (tStart-1) * frameBytes, 'bof');
    fwrite(fidOut, slab, 'single');
end

clear cIn cOut; % close fidIn/fidOut via safeFclose before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'spatialGaussFilt:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end

% =========================================================================
% Local helper: zero-fill / filter / restore with a supplied spatial mask
% =========================================================================
function outSlab = iApplyMaskedGauss(inSlab, sigma, spatialMask)
%IAPPLYMASKEDGAUSS Filter a YXT slab, keeping masked pixels out of the kernel.
%
% spatialMask must be the same size as inSlab (a per-frame NaN mask, not
% collapsed across T), so a pixel that is NaN in one frame does not cause
% valid data in other frames to be discarded.

outSlab = inSlab;
outSlab(spatialMask) = 0;
outSlab = imgaussfilt(outSlab, sigma, 'FilterDomain', 'spatial');
outSlab(spatialMask) = NaN;
end

% =========================================================================
% Local helper: raw .dat dimensions from AcqInfos.mat / file size
% =========================================================================
function [Ny, Nx, Nt] = iGetRawDatInfo(SaveFolder, inFile)
%IGETRAWDATINFO Return YXT dimensions for a raw .dat file.
%
% Prefer loadMetaData(...) so Nt follows the actual file size rather than a
% potentially stale AcqInfoStream.Length value.

if ~isfolder(SaveFolder)
    error('spatialGaussFilt:MissingSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

meta = loadMetaData(inFile);

if ~isfield(meta, 'Height') || ~isfield(meta, 'Width') || ~isfield(meta, 'datLength')
    error('spatialGaussFilt:InvalidMetaData', ...
        'loadMetaData did not return Height, Width, and datLength for "%s".', inFile);
end

Ny = double(meta.Height);
Nx = double(meta.Width);
Nt = double(meta.datLength);
end

% =========================================================================
% Local helper: validate/extract UMT data
% =========================================================================
function [entryNames, entryData, entryDims, labels, eventInfo, hasE, entryMetas] = iExtractValidUMTData(umt)
%IEXTRACTVALIDUMTDATA Validate and extract image-backed UMT entries.

validateUMTStruct(umt, 'requireEventInfo', false);

if ~strcmpi(umt.kind, 'image')
    error('spatialGaussFilt:InvalidUMTKind', ...
        ['Operation aborted. UMT input must have kind = "image". ' ...
         'This function does not support non-image UMT structures.']);
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error('spatialGaussFilt:EmptyUMTData', ...
        'Operation aborted. UMT data is empty.');
end

entryData = cell(size(entryNames));
entryDims = cell(size(entryNames));
entryMetas = cell(size(entryNames));
hasE = false(size(entryNames));

allowed = { ...
    {'Y','X'}, ...
    {'Y','X','T'}, ...
    {'Y','X','E'}, ...
    {'Y','X','T','E'}};

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    isAllowed = any(cellfun(@(x) isequal(thisDims, x), allowed));
    if ~isAllowed
        error('spatialGaussFilt:InvalidUMTEntry', ...
            ['Operation aborted. All UMT entries must use dimNames ' ...
             '{''Y'',''X''}, {''Y'',''X'',''T''}, {''Y'',''X'',''E''}, or ' ...
             '{''Y'',''X'',''T'',''E''}. Invalid entry: "%s".'], ...
            entryNames{iEntry});
    end

    entryData{iEntry} = thisEntry.value;
    entryDims{iEntry} = thisDims;
    hasE(iEntry) = any(strcmp(thisDims, 'E'));

    if isfield(thisEntry, 'meta') && isstruct(thisEntry.meta) && isscalar(thisEntry.meta)
        entryMetas{iEntry} = thisEntry.meta;
    else
        entryMetas{iEntry} = struct();
    end
end

if any(hasE)
    if ~isfield(umt, 'eventInfo')
        error('spatialGaussFilt:MissingEventInfo', ...
            ['Operation aborted. The input UMT contains entries with an E ' ...
             'dimension but has no shared top-level eventInfo.']);
    end
    eventInfo = umt.eventInfo;
else
    eventInfo = struct();
end

if isfield(umt, 'labels')
    labels = umt.labels;
else
    labels = struct();
end
end