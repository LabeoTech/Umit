function outData = GSR(data, SaveFolder, varargin)
% GSR Perform Global Signal Regression on image time series.
%
% This function removes global fluctuations from imaging data by regressing
% out the global mean signal.
%
% The function supports two execution modes:
%
%   1) STANDARD MODE (in-memory)
%      - Triggered when "data" is a single array
%      - Entire dataset is processed in RAM
%
%   2) LOW-RAM MODE (streaming)
%      - Triggered when "data" is a .dat filename
%      - Data are processed in spatial chunks directly from disk
%
% Syntax:
%   outData = GSR(data, SaveFolder)
%   outData = GSR(data, SaveFolder, 'Name', Value, ...)
%   info    = GSR('pipelineInfo')
%
% Inputs:
%   data :
%       Either:
%         - single array [Y, X, T] or [Y, X, T, E]   -> STANDARD MODE
%         - Character/string .dat filename with axes
%           Y-X-T or Y-X-T-E                         -> LOW-RAM MODE
%       UMT structs and .umt files are not supported.
%
%   SaveFolder :
%       Folder containing AcqInfos.mat and, optionally, DataParams.mat.
%
% Name-Value parameters:
%   UseMask :
%       Logical scalar. If true, GSR is computed only inside the logical
%       mask stored in DataParams.mat. If DataParams.mat is missing, or if
%       DataParams.mask.logical is empty or all true, a warning is raised
%       and the function falls back to using the full frame.
%
% Output:
%   outData :
%       - STANDARD MODE: corrected data array, same size as the input
%       - LOW-RAM MODE : full path to corrected .dat file ("GSR.dat" in
%                        SaveFolder) with the same axes as the input
%
% Notes:
%   - All data inputs are assumed to be single precision.
%   - Invalid traces are identified from the first frame only. A NaN in
%     frame 1 indicates that the whole pixel trace is invalid across time.
%   - For event-split data (E axis), every E slice is regressed as an
%     independent Y-X-T recording: its own dataset mean, global signal and
%     invalid traces (first frame of that slice). The mask is shared by all
%     E slices. Ignored event instances are processed like the rest, so the
%     E axis of a .dat output still matches events.mat.

% -------------------------------------------------------------------------
% pipelineInfo query
% -------------------------------------------------------------------------
if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localGetPipelineInfo();
    return
end

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'data', @(x) ...
    (isa(x, 'single') && (ndims(x) == 3 || ndims(x) == 4)) || ...
    ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'SaveFolder', @isfolder);
addParameter(p, 'UseMask', false, @(x) islogical(x) && isscalar(x));
parse(p, data, SaveFolder, varargin{:});

dataIn = p.Results.data;
SaveFolder = p.Results.SaveFolder;
UseMask = p.Results.UseMask;

clear p

% -------------------------------------------------------------------------
% Resolve execution mode and metadata
% -------------------------------------------------------------------------
bLowRAM = ischar(dataIn) || (isstring(dataIn) && isscalar(dataIn));

if bLowRAM
    dataFile = localResolveDataFile(char(string(dataIn)), SaveFolder);
    metaData = loadMetaData(dataFile);
    assertDatLayout(metaData, {{'Y','X','T'}, {'Y','X','T','E'}}, 'GSR');

    assert(strcmpi(metaData.dataClass, 'single'), ...
        'Umitoolbox:GSR:InvalidInput', ...
        'GSR currently supports only single-precision .dat inputs.');

    frameSize = [datAxisSize(metaData, 'Y'), datAxisSize(metaData, 'X')];
else
    dataFile = '';
    metaData = struct();
    frameSize = size(dataIn, [1 2]);
end

% -------------------------------------------------------------------------
% Load or build logical mask
% -------------------------------------------------------------------------
logical_mask = localGetLogicalMask(SaveFolder, frameSize, UseMask);

% -------------------------------------------------------------------------
% Dispatch execution mode
% -------------------------------------------------------------------------
if bLowRAM
    outData = GSR_lowRAMmode(dataFile, SaveFolder, metaData, logical_mask);
else
    outData = GSR_standardMode(dataIn, logical_mask);
end

disp('Finished GSR.');

% -------------------------------------------------------------------------
% Local pipelineInfo factory
% -------------------------------------------------------------------------
    function info = localGetPipelineInfo()

        allowedUseMask = [true, false];

        info = PipelineManager.createPipelineInfo( ...
            mfilename, ...
            'Perform Global Signal Regression on image time series.');

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'UnknownDataType','ImageTimeSeries','ProcessedData'}, ...
            ['Input image time series (YXT or YXTE). Accepts in-memory data ' ...
             'or a file-backed .dat input.'], ...
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
            'Folder containing AcqInfos.mat and optional DataParams.mat.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'UseMask', ...
            'parameter', ...
            'If true, use DataParams.mask.logical when available.', ...
            'kind', 'parameter', ...
            'default', false, ...
            'allowed', allowedUseMask, ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'GSR-corrected data.', ...
            'GSR.dat', ...
            1, ...
            'isData', true, ...
            'saveFileName', 'GSR.dat');
    end

% -------------------------------------------------------------------------
% Local helper: resolve logical mask from DataParams
% -------------------------------------------------------------------------
    function logical_mask = localGetLogicalMask(saveFolder, frameSizeLocal, useMaskLocal)

        logical_mask = true(frameSizeLocal);

        if ~useMaskLocal
            return
        end

        dataParamsFile = fullfile(saveFolder, 'DataParams.mat');
        if ~isfile(dataParamsFile)
            warning('Umitoolbox:GSR:MissingDataParams', ...
                'DataParams.mat not found in "%s". Falling back to UseMask = false.', saveFolder);
            return
        end

        S = load(dataParamsFile, 'DataParams');
        assert(isfield(S, 'DataParams'), ...
            'Umitoolbox:GSR:InvalidInput', ...
            'File "%s" does not contain variable "DataParams".', dataParamsFile);
        if ~(isfield(S, 'DataParams') && isfield(S.DataParams, 'mask') && isfield(S.DataParams.mask, 'logical'))
            warning('Umitoolbox:GSR:MissingMask', ...
                'Logical mask was not set. Falling back to UseMask = false.');
            return
        end

        DataParams = S.DataParams;
        logical_mask = DataParams.mask.logical;

        % Check for empty masks
        if isempty(logical_mask)
            warning('Umitoolbox:GSR:TrivialMask', ...
                'DataParams.mask.logical is empty. Falling back to UseMask = false.');
            logical_mask = true(frameSizeLocal);
            return
        end
        % Test size
        assert(isequal(size(logical_mask), frameSizeLocal), ...
            'Umitoolbox:GSR:InvalidInput', ...
            'Logical mask size does not match the data frame size.');
        % Warn if mask is all True
        if all(logical_mask(:))
            warning('Umitoolbox:GSR:TrivialMask', ...
                'DataParams.mask.logical is all TRUE. Falling back to UseMask = false.');
            logical_mask = true(frameSizeLocal);
            return
        end

        logical_mask = logical(logical_mask);
    end

% -------------------------------------------------------------------------
% Local helper: resolve data file path
% -------------------------------------------------------------------------
    function dataFileOut = localResolveDataFile(dataFileIn, saveFolder)

        if isfile(dataFileIn)
            dataFileOut = dataFileIn;
            return
        end

        candidate = fullfile(saveFolder, dataFileIn);
        if isfile(candidate)
            dataFileOut = candidate;
            return
        end

        error('Umitoolbox:GSR:MissingInput', ...
            'Input data file "%s" was not found.', dataFileIn);
    end

end

%--------------------------------------------------------------------------
% Local functions
%--------------------------------------------------------------------------
function outData = GSR_standardMode(data, logical_mask)
% GSR_STANDARDMODE In-memory Global Signal Regression.
%
% Note: This version mirrors the low-RAM mode by using the same
% sum/count approach for dataset mean and global-signal estimation (shared
% helpers). Every E slice of a Y x X x T x E array is regressed as an
% independent Y x X x T recording.

disp('Calculating Global Signal Regression...');
stats = iAccumulateStats(data, logical_mask);
[mData, Sig] = iFinalizeStats(stats);

outData = iRegressSlab(data, Sig, mData);
end


function outFileName = GSR_lowRAMmode(dataFile, SaveFolder, metaData, logical_mask)
% GSR_LOWRAMMODE Disk-streamed Global Signal Regression.
%
% This function performs GSR using a low-RAM, two-pass strategy:
%   1) First pass computes the dataset mean and global signal over time
%      (per E slice for event-split data)
%   2) Second pass regresses it out chunk-by-chunk
%
% The numerical intent matches GSR_standardMode.

% -------------------------------------------------------------------------
% Estimate data size and chunking
% -------------------------------------------------------------------------
Ny = datAxisSize(metaData, 'Y');
Nx = datAxisSize(metaData, 'X');
Nt = datAxisSize(metaData, 'T');
Ne = 1;
if any(strcmp(cellstr(string(metaData.dimNames)), 'E'))
    Ne = datAxisSize(metaData, 'E');
end

dataBytes = Ny * Nx * Nt * Ne * 4;
nChunks = calculateMaxChunkSize(dataBytes, 3, .1);

chunkSizePixels = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkSizePixels);

% -------------------------------------------------------------------------
% File handles
% -------------------------------------------------------------------------
slabIn = spatialSlabIO('open', dataFile, 'Info', metaData);
c_in = onCleanup(@() spatialSlabIO('close', slabIn));

outFileName = fullfile(SaveFolder, 'GSR.dat');
slabOut = spatialSlabIO('create', outFileName, datHeaderFromInfo(metaData, 'GSR'));
c_out = onCleanup(@() spatialSlabIO('close', slabOut));

h = waitbar(0, 'GSR: computing global signal (pass 1)...');
h.Name = 'GSR (pass 1/2)';
cleanupWaitbar = onCleanup(@() iCloseWaitbarSafely(h));

% =====================================================================
% PASS 1 - Compute dataset mean and global signal
% =====================================================================
stats = iZeroStats(Nt, Ne);

for ii = 1:nChunks
    waitbar(ii / nChunks, h, 'GSR: computing global signal...');

    pxStart = (ii-1) * chunkSizePixels + 1;
    pxEnd   = min(ii * chunkSizePixels, Nx);
    idxX    = pxStart:pxEnd;

    slab = spatialSlabIO('read', slabIn, idxX);
    stats = iAddStats(stats, iAccumulateStats(slab, logical_mask(:, idxX)));

    clear slab
end

[mData, Sig] = iFinalizeStats(stats);

% =====================================================================
% PASS 2 - Regress global signal chunk-by-chunk
% =====================================================================
h.Name = 'GSR (pass 2/2)';

for ii = 1:nChunks
    waitbar(ii / nChunks, h, 'GSR: regressing signal...');

    pxStart = (ii-1) * chunkSizePixels + 1;
    pxEnd   = min(ii * chunkSizePixels, Nx);
    idxX    = pxStart:pxEnd;

    slab = spatialSlabIO('read', slabIn, idxX);
    slab = iRegressSlab(slab, Sig, mData);
    spatialSlabIO('write', slabOut, idxX, slab);

    clear slab
end

spatialSlabIO('finalize', slabOut);
end

% =========================================================================
% Helpers shared by standard and low-RAM mode
% =========================================================================
function stats = iZeroStats(Nt, Ne)
%IZEROSTATS Empty accumulators: one column per E slice.

stats = struct( ...
    'dataSum',     zeros(1, Ne), ...
    'dataCount',   zeros(1, Ne), ...
    'globalSum',   zeros(Nt, Ne), ...
    'globalCount', zeros(1, Ne));
end

function stats = iAccumulateStats(slab, maskSlab)
%IACCUMULATESTATS Sums behind the dataset mean and the global signal.
%
% SLAB is Y x X x T (x E); MASKSLAB is the logical Y x X mask of the same
% columns. Every E slice is accumulated on its own. Invalid traces are
% identified from the first frame of the slice, and excluded from the
% global signal.

Nt = size(slab, 3);
Ne = size(slab, 4);
stats = iZeroStats(Nt, Ne);

for iE = 1:Ne
    trial = slab(:,:,:,iE);

    stats.dataSum(iE)   = sum(double(trial(:)), 'omitnan');
    stats.dataCount(iE) = sum(~isnan(trial(:)));

    idx_invalid_trace = isnan(trial(:,:,1));
    maskIdx = maskSlab(:) & ~idx_invalid_trace(:);

    if any(maskIdx)
        trial2D = reshape(trial, [], Nt);
        stats.globalSum(:, iE) = sum(double(trial2D(maskIdx, :)), 1).';
        stats.globalCount(iE)  = sum(maskIdx);
    end
end
end

function stats = iAddStats(stats, chunkStats)
%IADDSTATS Add the accumulators of one chunk.

stats.dataSum     = stats.dataSum     + chunkStats.dataSum;
stats.dataCount   = stats.dataCount   + chunkStats.dataCount;
stats.globalSum   = stats.globalSum   + chunkStats.globalSum;
stats.globalCount = stats.globalCount + chunkStats.globalCount;
end

function [mData, Sig] = iFinalizeStats(stats)
%IFINALIZESTATS Dataset mean (1 x E) and normalized global signal (T x E).

Ne = numel(stats.dataSum);
Nt = size(stats.globalSum, 1);
mData = zeros(1, Ne, 'single');
Sig = zeros(Nt, Ne, 'single');

for iE = 1:Ne
    where = '';
    if Ne > 1
        where = sprintf(' (E slice %d)', iE);
    end

    if stats.dataCount(iE) == 0
        error('Umitoolbox:GSR:InvalidInput', ...
            'Input data contain only NaN values%s.', where);
    end

    if stats.globalCount(iE) == 0
        error('Umitoolbox:GSR:InvalidInput', ...
            'Logical mask does not contain any valid pixels for GSR%s.', where);
    end

    mData(iE) = single(stats.dataSum(iE) / stats.dataCount(iE));
    thisSig   = single(stats.globalSum(:, iE) ./ stats.globalCount(iE));

    sigMean = mean(thisSig);
    if ~isfinite(sigMean) || sigMean == 0
        error('Umitoolbox:GSR:InvalidInput', ...
            'Global signal mean is zero or invalid. Cannot normalize signal%s.', where);
    end

    Sig(:, iE) = thisSig ./ sigMean;
end
end

function slab = iRegressSlab(slab, Sig, mData)
%IREGRESSSLAB Regress the global signal out of a Y x X x T (x E) array.
%
% Each E slice uses its own signal (column of SIG) and dataset mean (entry
% of MDATA). Invalid traces (NaN in the first frame) are zeroed for the
% regression and restored to NaN afterwards.

slabSz = size(slab);
Nt = slabSz(3);
constantTerm = ones(Nt, 1, 'single');

for iE = 1:size(slab, 4)
    trial2D = reshape(slab(:,:,:,iE), [], Nt);
    idx_invalid_trace = isnan(trial2D(:, 1));

    trial2D(idx_invalid_trace, :) = 0;

    X = [constantTerm, Sig(:, iE)];
    A = X * (X \ trial2D');
    trial2D = trial2D - A';

    trial2D = trial2D + mData(iE);
    trial2D(idx_invalid_trace, :) = NaN;

    slab(:,:,:,iE) = reshape(trial2D, slabSz(1), slabSz(2), Nt);
end
end

function iCloseWaitbarSafely(h)
%ICLOSEWAITBARSAFELY Close a waitbar handle if it still exists.

if isgraphics(h)
    close(h);
end
end
