function outData = apply_detrend(data, SaveFolder, varargin)
%APPLY_DETREND Apply linear detrending along the time dimension.
%
%   outData = apply_detrend(data, SaveFolder)
%   outData = apply_detrend(data, SaveFolder, 'FrameRateHz', rate)
%
%   FrameRateHz is the frame rate of DATA, used to convert the events.mat
%   baseline period into frames. PipelineManager injects it from the data;
%   a .dat input's header provides it otherwise. In-RAM arrays need it
%   explicitly when events.mat defines a baseline period; AcqInfos.mat is
%   not used.
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
%          - Triggered when "data" is a numeric array
%       2) LOW-RAM MODE (file-backed)
%          - Triggered when "data" is a .dat filename
%          - The file is read, detrended, and written in X slabs, so the
%            recording is never loaded whole
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X x T
%       2) Numeric array with dimensions Y x X x T x E
%       3) Filename to a .dat file with axes Y-X-T or Y-X-T-E
%       UMT structs and .umt files are not supported.
%
%   Input/output behavior:
%       - Numeric array in  -> numeric array out, same size
%       - .dat filename in  -> .dat filename out ("data_detrended.dat" in
%         SaveFolder), same axes and sizes
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder used for file resolution, events.mat lookup,
%                    and the .dat output.
%
%   Output:
%       outData    - Detrended output with the same representation type
%                    and dimensions as the input.
%
%   Notes:
%       - Raw .dat files are assumed to store single-precision data.
%       - For event-split data (E axis) every trial is detrended on its own
%         along T, ignored event instances included, so the E axis of a .dat
%         output still matches events.mat.
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
addParameter(p, 'FrameRateHz', []);
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
explicitRate = p.Results.FrameRateHz;

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

    frames = iGetDetrendFrameCount(SaveFolder, size(data, 3), explicitRate, data);

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
            outData = iApplyDetrendDatFile(dataFile, SaveFolder, default_Output, explicitRate);
            return

        otherwise
            error('apply_detrend:UnsupportedInputFile', ...
                'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    end
end

error('apply_detrend:UnsupportedInputType', ...
    ['Input "data" must be a YXT/YXTE array or a .dat filename. ' ...
     'UMT structs and .umt files are not supported.']);

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
            ['Input data. Accepted forms: YXT or YXTE array, or a .dat ' ...
             'filename with axes Y-X-T or Y-X-T-E.'], ...
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
            'Folder used for file resolution, events.mat lookup, and the .dat output.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

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
function frames = iGetDetrendFrameCount(SaveFolder, Nt, explicitRate, rateData)
%IGETDETRENDFRAMECOUNT Determine the detrend baseline window in frames.
%
% With an events.mat baseline period, the window is that period in frames
% of the data's own frame rate: the explicit FrameRateHz (injected by
% PipelineManager), else the .dat header of RATEDATA; without one,
% resolveDataInfoValue raises an error. AcqInfos.mat is not used. Without a
% baseline period the default window of 7 frames is used. NT is the length
% of one trace (the trial length for event-split data).

frames = 7;
baselineSec = [];

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

if ~isempty(baselineSec)
    freqHz = resolveDataInfoValue('frameRateHz', explicitRate, rateData, 'apply_detrend');
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
% Helper: low-RAM .dat execution for YXT and YXTE data
% =========================================================================
function outFile = iApplyDetrendDatFile(inFile, SaveFolder, defaultOutput, explicitRate)
%IAPPLYDETRENDDATFILE Apply detrending to a Y-X-T or Y-X-T-E .dat file.
%
% The file is processed in X slabs (all T and E of a slab at once), so the
% recording is never loaded whole. Every trial of an event-split file is
% detrended on its own, exactly as the in-memory path does.

slabIn = spatialSlabIO('open', inFile);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));
assertDatLayout(slabIn.Info, {{'Y','X','T'}, {'Y','X','T','E'}}, 'apply_detrend');
Ny = slabIn.Ny;
Nx = slabIn.Nx;
Nt = datAxisSize(slabIn.Info, 'T');
Ne = 1;
if any(strcmp(cellstr(string(slabIn.Info.dimNames)), 'E'))
    Ne = datAxisSize(slabIn.Info, 'E');
end
frames = iGetDetrendFrameCount(SaveFolder, Nt, explicitRate, inFile);

% Write through a scratch file so the declared pipeline output only appears
% once the run has completed, and so the input can safely be the file that
% the declared output would overwrite (a pipeline re-run).
outFile = fullfile(SaveFolder, defaultOutput);
[~, defOutFilename, ext] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [defOutFilename, '_writing', ext]);
slabOut = spatialSlabIO('create', tmpFile, ...
    datHeaderFromInfo(slabIn.Info, defOutFilename, 'dataClass', 'single'));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

% Peak memory of a slab: the slab, the detrended copy, and the trend and
% difference temporaries of the 2-D algorithm.
nChunks = calculateMaxChunkSize(Nx * Ny * Nt * Ne * 4, 4, 0.3);
chunkX = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkX);

for c = 1:nChunks
    xStart = (c-1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    fprintf('Chunk %i/%i [Reading file ...]\n', c, nChunks)
    slab = single(spatialSlabIO('read', slabIn, xIdx));
    slabSize = size(slab);

    fprintf('Chunk %i/%i [Detrending data ...]\n', c, nChunks)
    slab = reshape(slab, Ny * numel(xIdx), Nt, Ne);
    for iEvent = 1:Ne
        slab(:, :, iEvent) = iApplyDetrend2D(slab(:, :, iEvent), frames);
    end
    slab = reshape(slab, slabSize);

    fprintf('Chunk %i/%i [Writing to file ...]\n', c, nChunks)
    spatialSlabIO('write', slabOut, xIdx, slab);
    fprintf('Chunk %i/%i [Completed]\n', c, nChunks)
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close the input reader before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'apply_detrend:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);

disp('Finished detrend!');
end

