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

default_Output = 'data_detrended.dat';

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
    disp('Finished detrend!');
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

[entryNames, entryData, entryDims, labels, sourceEventInfo, hasE, entryMetas] = ...
    iExtractValidUMTData(data);

outStruct = data;
outStruct.data = struct();
if isfield(outStruct, 'eventInfo')
    % Strip before rebuilding entries one at a time below: with mixed
    % YXT/YXTE entries, a non-E entry can be appended before any E entry
    % exists yet, which would make a carried-forward eventInfo temporarily
    % inconsistent with the entries seen so far. eventInfo is re-attached
    % once, after every entry is in place, via appendUMTEventInfo below.
    outStruct = rmfield(outStruct, 'eventInfo');
end

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

    outStruct = genUMTStruct( ...
        outStruct, ...
        'value', detrended, ...
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

% Unchanged from the legacy algorithm (see function header): the slope M is
% estimated over (Nt - frames) samples but the axis below spans (Nt - 1)
% samples starting at -2, so trend(1) = -2*M + b rather than 0. This is a
% preserved quirk of the original implementation, not a derived formula --
% do not "fix" the axis without re-validating against legacy output.
trend = M .* linspace(-2, Nt-3, Nt) + b;
out2D = in2D - trend + b;
end

% =========================================================================
% Helper: determine baseline frame count
% =========================================================================
function frames = iGetDetrendFrameCount(SaveFolder, Nt, dataFile)
%IGETDETRENDFRAMECOUNT Determine the detrend baseline window in frames.
%
% When dataFile is given, prefer loadMetaData(dataFile) so the reported
% frame rate stays consistent with that specific file's on-disk size.
% Never select an arbitrary, unrelated file from SaveFolder: when there is
% no file being processed (a raw numeric array or an in-RAM UMT struct),
% fall back to the single authoritative AcqInfos.mat instead.

frames = 7;

freqHz = [];
baselineSec = [];

if nargin > 2 && ~isempty(dataFile)
    try
        meta = loadMetaData(dataFile);
        if isfield(meta, 'Freq') && ~isempty(meta.Freq)
            freqHz = double(meta.Freq);
        end
    catch ME
        warning('apply_detrend:FrameRateFromFileFailed', ...
            'Could not resolve frame rate from "%s": %s', dataFile, ME.message);
    end
else
    try
        acqInfoFile = fullfile(SaveFolder, 'AcqInfos.mat');
        if isfile(acqInfoFile)
            S = load(acqInfoFile, 'AcqInfoStream');
            if isfield(S, 'AcqInfoStream') && isfield(S.AcqInfoStream, 'FrameRateHz') && ...
                    ~isempty(S.AcqInfoStream.FrameRateHz)
                freqHz = double(S.AcqInfoStream.FrameRateHz);
            end
        end
    catch ME
        warning('apply_detrend:FrameRateFromAcqInfosFailed', ...
            'Could not resolve frame rate from "%s": %s', acqInfoFile, ME.message);
    end
end

eventsFile = fullfile(SaveFolder, 'events.mat');
if isfile(eventsFile)
    try
        evObj = EventsManager(SaveFolder);
        if ~isempty(evObj.baselinePeriod)
            baselineSec = double(evObj.baselinePeriod);
        end
    catch ME
        warning('apply_detrend:BaselinePeriodResolutionFailed', ...
            'Could not resolve baselinePeriod from "%s": %s', eventsFile, ME.message);
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

[Ny, Nx, Nt] = getRawDatInfo(SaveFolder, inFile);
frames = iGetDetrendFrameCount(SaveFolder, Nt, inFile);

% Write through a scratch file so the declared pipeline output only appears
% once the run has completed, and so the input can safely be the file that
% the declared output would overwrite (a pipeline re-run).
outFile = fullfile(SaveFolder, defaultOutput);
[~, defOutFilename, ext] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [defOutFilename, '_writing', ext]);
preallocateDatFile(tmpFile, [Ny, Nx, Nt], 'single');

fidIn = fopen(inFile, 'r');
assert(fidIn ~= -1, 'apply_detrend:OpenInputFailed', ...
    'Failed to open input file "%s".', inFile);
cIn = onCleanup(@() safeFclose(fidIn));

fidOut = fopen(tmpFile, 'r+');
assert(fidOut ~= -1, 'apply_detrend:OpenOutputFailed', ...
    'Failed to open output file "%s".', tmpFile);
cOut = onCleanup(@() safeFclose(fidOut));

nChunks = calculateMaxChunkSize(Nx * Ny * Nt * 4, 2, 0.3);
chunkX = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkX);

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

clear cIn cOut; % close fidIn/fidOut via safeFclose before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'apply_detrend:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);

disp('Finished detrend!');
end

% =========================================================================
% Helper: validate/extract UMT data
% =========================================================================
function [entryNames, entryData, entryDims, labels, eventInfo, hasE, entryMetas] = iExtractValidUMTData(umt)
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
entryMetas = cell(size(entryNames));
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

    if isfield(thisEntry, 'meta') && isstruct(thisEntry.meta) && isscalar(thisEntry.meta)
        entryMetas{iEntry} = thisEntry.meta;
    else
        entryMetas{iEntry} = struct();
    end
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