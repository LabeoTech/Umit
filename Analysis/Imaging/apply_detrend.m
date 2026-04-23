function outData = apply_detrend(data, SaveFolder, varargin)
%APPLY_DETREND Apply linear detrending along the time dimension.
%
%   outData = apply_detrend(data, SaveFolder)
%
%   This function applies the existing linear detrend algorithm along the
%   time dimension T. The algorithm is unchanged from the legacy version:
%       1) Estimate a baseline from the first N frames
%       2) Estimate a terminal level from the last N frames
%       3) Build a linear trend
%       4) Subtract the trend while preserving the initial baseline offset
%
%   Supported execution modes:
%       1) STANDARD MODE (in-memory)
%          - Triggered when "data" is a numeric array or a UMT struct
%       2) LOW-RAM MODE (file-backed)
%          - Triggered when "data" is a .dat filename
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X x T
%       2) Numeric array with dimensions Y x X x T x E
%       3) Filename to a .dat file storing continuous Y x X x T data
%       4) UMT struct with image entries using dimensions:
%              {'Y','X','T'}
%              {'Y','X','T','E'}
%       5) Filename to a .umt file containing one UMT struct
%
%   Input/output behavior:
%       - Numeric array in  -> numeric array out
%       - .dat filename in  -> .dat filename out
%       - UMT struct in     -> UMT struct out
%       - .umt filename in  -> UMT struct out
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder used for file resolution and metadata lookup.
%
%   Output:
%       outData    - Detrended output with the same representation type as
%                    the input.
%
%   Notes:
%       - Raw .dat files are assumed to store continuous YXT data in single
%         precision.
%       - Low-RAM mode is implemented only for raw .dat input.
%       - For event-split data, the current convention is YXTE.
%       - If events.mat exists and contains a baseline period, that value is
%         converted to a frame count and used to determine the detrend
%         baseline window. Otherwise a default of 7 frames is used.

default_Output = 'data_detrended.dat'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));

if ~isfolder(SaveFolder)
    error('apply_detrend:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Case 1: Numeric array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)

    validateattributes(data, {'numeric','logical'}, {'nonempty'}, mfilename, 'data');

    if ~(ndims(data) == 3 || ndims(data) == 4)
        error('apply_detrend:InvalidArrayInput', ...
            'Numeric input must be YXT or YXTE.');
    end

    frames = iGetDetrendFrameCount(SaveFolder, size(data, 3));

    if ndims(data) == 3
        outData = iApplyDetrendToYXT(data, frames);
    else
        outData = iApplyDetrendToYXTE(data, frames);
    end
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
            error('apply_detrend:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);

    switch ext
        case '.dat'
            outData = iApplyDetrendDatFile(dataFile, SaveFolder, default_Output);
            return

        case '.umt'
            warning('apply_detrend:UMTFileLoadsInRAM', ...
                ['RAM-safe mode is not available for data stored in this format. ' ...
                 'Loading the UMT content into RAM.']);
            data = loadData(dataFile);

        otherwise
            error('apply_detrend:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

% -------------------------------------------------------------------------
% Case 3: UMT struct
% -------------------------------------------------------------------------
if ~isstruct(data)
    error('apply_detrend:UnsupportedInputType', ...
        ['Input "data" must be a YXT/YXTE array, a .dat filename, ' ...
         'a UMT struct, or a .umt filename containing a UMT struct.']);
end

[entryNames, entryData, entryDims, labels, sourceEventInfo, hasE] = ...
    iExtractValidUMTData(data);

outStruct = data;
outStruct.data = struct();

for iEntry = 1:numel(entryNames)

    value = entryData{iEntry};
    dimNames = entryDims{iEntry};
    frames = iGetDetrendFrameCount(SaveFolder, size(value, 3));

    switch strjoin(dimNames, '')
        case 'YXT'
            detrended = iApplyDetrendToYXT(value, frames);

        case 'YXTE'
            detrended = iApplyDetrendToYXTE(value, frames);

        otherwise
            error('apply_detrend:InvalidUMTEntryDims', ...
                'Unsupported dimNames in entry "%s".', entryNames{iEntry});
    end

    outStruct.data.(entryNames{iEntry}) = struct( ...
        'value', detrended, ...
        'dimNames', {dimNames});
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
disp('Finished detrend!');

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Apply linear detrending along the time dimension.');

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT array, YXTE array, .dat filename, ' ...
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
            'Folder used for file resolution and metadata lookup.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Detrended output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: apply legacy detrend algorithm to one YXT block
% =========================================================================
function outBlock = iApplyDetrendToYXT(inBlock, frames)
%IAPPLYDETRENDTOYXT Apply the existing detrend algorithm to YXT data.

origSz = size(inBlock);
slab = reshape(inBlock, [], origSz(3));
slab = iApplyDetrend2D(slab, frames);
outBlock = reshape(slab, origSz);
end

% =========================================================================
% Helper: apply legacy detrend algorithm to YXTE data
% =========================================================================
function outBlock = iApplyDetrendToYXTE(inBlock, frames)
%IAPPLYDETRENDTOYXTE Apply the existing detrend algorithm to YXTE data.

outBlock = zeros(size(inBlock), 'like', inBlock);

for iEvent = 1:size(inBlock, 4)
    outBlock(:,:,:,iEvent) = iApplyDetrendToYXT(inBlock(:,:,:,iEvent), frames);
end
end

% =========================================================================
% Helper: core 2-D detrend logic [pixels x time]
% =========================================================================
function out2D = iApplyDetrend2D(in2D, frames)
%IAPPLYDETREND2D Apply the unchanged detrend logic to [Npix x Nt] data.

Nt = size(in2D, 2);

frames = min(frames, Nt);
frames = max(frames, 3);

if mod(frames, 2) == 0
    frames = frames - 1;
end

if frames >= Nt
    frames = max(3, Nt - 1);
    if mod(frames, 2) == 0
        frames = frames - 1;
    end
end

if frames < 3 || Nt < 3
    out2D = in2D;
    return
end

delta_y = median(in2D(:, end-frames+1:end), 2, 'omitnan') - ...
          median(in2D(:, 1:frames), 2, 'omitnan');

delta_x = Nt - frames;
M = delta_y ./ delta_x;
b = median(in2D(:, 1:frames), 2, 'omitnan');

trend = M .* linspace(-2, Nt-3, Nt) + b;
out2D = in2D - trend + b;
end

% =========================================================================
% Helper: determine baseline frame count
% =========================================================================
function frames = iGetDetrendFrameCount(SaveFolder, Nt)
%IGETDETRENDFRAMECOUNT Determine the detrend baseline window in frames.
%
% Prefer loadMetaData(...) over direct AcqInfos.mat reading so the reported
% temporal metadata stays consistent with the actual on-disk .dat size.

frames = 7;

freqHz = [];
baselineSec = [];

candidateFiles = [dir(fullfile(SaveFolder, '*.dat')); dir(fullfile(SaveFolder, '*.umt'))];
if ~isempty(candidateFiles)
    try
        meta = loadMetaData(fullfile(SaveFolder, candidateFiles(1).name));
        if isfield(meta, 'Freq') && ~isempty(meta.Freq)
            freqHz = double(meta.Freq);
        end
    catch
    end
end

eventsFile = fullfile(SaveFolder, 'events.mat');
if isfile(eventsFile)
    try
        evObj = EventsManager(SaveFolder, '', 'csv');
        if ~isempty(evObj.baselinePeriod)
            baselineSec = double(evObj.baselinePeriod);
        end
    catch
    end
end

if ~isempty(freqHz) && ~isempty(baselineSec)
    frames = round(baselineSec * freqHz);
end

frames = max(frames, 3);

if mod(frames, 2) == 0
    frames = frames + 1;
end

if nargin > 1 && ~isempty(Nt)
    frames = min(frames, Nt);
    if mod(frames, 2) == 0 && frames > 3
        frames = frames - 1;
    end
    frames = max(min(frames, Nt), 3);
end
end

% =========================================================================
% Helper: low-RAM .dat execution for continuous YXT data
% =========================================================================
function outFile = iApplyDetrendDatFile(inFile, SaveFolder, defaultOutput)
%IAPPLYDETRENDDATFILE Apply detrending to a raw continuous YXT .dat file.

[Ny, Nx, Nt] = iGetRawDatInfo(SaveFolder, inFile);
frames = iGetDetrendFrameCount(SaveFolder, Nt);
[~,defOutFilename,ext] = fileparts(defaultOutput);
outFile = fullfile(fileparts(inFile), [defOutFilename, '_PREALLOCATED_FILE', ext]);
preallocateDatFile(outFile, [Ny, Nx, Nt], 'single');

fidIn = fopen(inFile, 'r');
assert(fidIn ~= -1, 'apply_detrend:OpenInputFailed', ...
    'Failed to open input file "%s".', inFile);
cIn = onCleanup(@() safeFclose(fidIn)); %#ok<NASGU>

fidOut = fopen(outFile, 'r+');
assert(fidOut ~= -1, 'apply_detrend:OpenOutputFailed', ...
    'Failed to open output file "%s".', outFile);
cOut = onCleanup(@() safeFclose(fidOut)); %#ok<NASGU>

nChunks = calculateMaxChunkSize(Nx * Ny * Nt * 4, 2, 0.3);
chunkX = ceil(Nx / nChunks);

for c = 1:nChunks
    xStart = (c-1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    fprintf('Chunk %i/%i [Reading file ...]\n', c, nChunks)
    slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, 'single');

    fprintf('Chunk %i/%i [Detrending data ...]\n', c, nChunks)
    slab = reshape(slab, Ny * numel(xIdx), Nt);
    slab = iApplyDetrend2D(slab, frames);
    slab = reshape(slab, Ny, numel(xIdx), Nt);

    fprintf('Chunk %i/%i [Writing to file ...]\n', c, nChunks)
    spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, 'single', slab);
    fprintf('Chunk %i/%i [Completed]\n', c, nChunks)
end

fclose(fidIn);
fclose(fidOut);
disp('Finished detrend!');
end

% =========================================================================
% Helper: raw .dat dimensions from AcqInfos.mat / file size
% =========================================================================
function [Ny, Nx, Nt] = iGetRawDatInfo(SaveFolder, inFile)
%IGETRAWDATINFO Return YXT dimensions for a raw .dat file.
%
% Prefer loadMetaData(...) so Nt follows the real file size rather than a
% potentially stale AcqInfoStream.Length value.

if ~isfolder(SaveFolder)
    error('apply_detrend:MissingSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

meta = loadMetaData(inFile);

if ~isfield(meta, 'Height') || ~isfield(meta, 'Width') || ~isfield(meta, 'datLength')
    error('apply_detrend:InvalidMetaData', ...
        'loadMetaData did not return Height, Width, and datLength for "%s".', inFile);
end

Ny = double(meta.Height);
Nx = double(meta.Width);
Nt = double(meta.datLength);
end

% =========================================================================
% Helper: validate/extract UMT data
% =========================================================================
function [entryNames, entryData, entryDims, labels, eventInfo, hasE] = iExtractValidUMTData(umt)
%IEXTRACTVALIDUMTDATA Validate and extract image-backed UMT entries.

validateUMTStruct(umt, 'requireEventInfo', false);

if ~strcmpi(umt.kind, 'image')
    error('apply_detrend:InvalidUMTKind', ...
        ['Operation aborted. UMT input must have kind = "image". ' ...
         'This function does not support non-image UMT structures.']);
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error('apply_detrend:EmptyUMTData', ...
        'Operation aborted. UMT data is empty.');
end

entryData = cell(size(entryNames));
entryDims = cell(size(entryNames));
hasE = false(size(entryNames));

allowed = { ...
    {'Y','X','T'}, ...
    {'Y','X','T','E'}};

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    isAllowed = any(cellfun(@(x) isequal(thisDims, x), allowed));
    if ~isAllowed
        error('apply_detrend:InvalidUMTEntry', ...
            ['Operation aborted. All UMT entries must use dimNames ' ...
             '{''Y'',''X'',''T''} or {''Y'',''X'',''T'',''E''}. ' ...
             'Invalid entry: "%s".'], ...
            entryNames{iEntry});
    end

    entryData{iEntry} = thisEntry.value;
    entryDims{iEntry} = thisDims;
    hasE(iEntry) = any(strcmp(thisDims, 'E'));
end

if any(hasE)
    if ~isfield(umt, 'eventInfo')
        error('apply_detrend:MissingEventInfo', ...
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