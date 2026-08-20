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
%         - single array [Y, X, T]              -> STANDARD MODE
%         - Character/string .dat filename      -> LOW-RAM MODE
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
%       - STANDARD MODE: corrected data array [Y, X, T]
%       - LOW-RAM MODE : full path to corrected .dat file
%
% Notes:
%   - All data inputs are assumed to be single precision.
%   - Invalid traces are identified from the first frame only. A NaN in
%     frame 1 indicates that the whole pixel trace is invalid across time.

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
    (isa(x, 'single') && ndims(x) == 3) || ...
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

    assert(strcmpi(metaData.Datatype, 'single'), ...
        'Umitoolbox:GSR:InvalidInput', ...
        'GSR currently supports only single-precision .dat inputs.');

    frameSize = metaData.datSize;
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
            {'UnknownDataType','ImageTimeSeries'}, ...
            'Input image time series. Accepts in-memory data or a file-backed .dat input.', ...
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
            'ImageTimeSeries', ...
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
% Note: This version mirrors the low-RAM mode more closely by using the same
% sum/count approach for dataset mean and global-signal estimation.

szData = size(data);

% Identify invalid traces from the first frame only
idx_invalid_trace = isnan(data(:,:,1));

% Compute global mean from original data using sum/count
dataSum = sum(double(data(:)), 'omitnan');
dataCount = sum(~isnan(data(:)));

if dataCount == 0
    error('Umitoolbox:GSR:InvalidInput', ...
        'Input data contain only NaN values.');
end

mData = single(dataSum / dataCount);

% Reshape to [pixels x time]
data = reshape(data, [], szData(3));
idx_invalid_trace = idx_invalid_trace(:);

% Build valid mask for global-signal estimation
maskIdx = logical_mask(:) & ~idx_invalid_trace;
assert(any(maskIdx), ...
    'Umitoolbox:GSR:InvalidInput', ...
    'Logical mask does not contain any valid pixels for GSR.');

% Replace invalid traces with zeros before regression
data(idx_invalid_trace, :) = 0;

% Compute global signal using sum/count
disp('Calculating Global Signal Regression...');
globalSum = sum(double(data(maskIdx, :)), 1);
globalCount = sum(maskIdx);

Sig = single(globalSum ./ globalCount);

sigMean = mean(Sig);
assert(isfinite(sigMean) && sigMean ~= 0, ...
    'Umitoolbox:GSR:InvalidInput', ...
    'Global signal mean is zero or invalid. Cannot normalize signal.');
Sig = Sig ./ sigMean;

% Regression
X = [ones(szData(3),1,'single'), Sig(:)];
A = X * (X \ data');
data = data - A';

% Restore mean
data = data + mData;

% Restore invalid traces
data(idx_invalid_trace, :) = NaN;

% Reshape back
outData = reshape(data, szData);
end


function outFileName = GSR_lowRAMmode(dataFile, SaveFolder, metaData, logical_mask)
% GSR_LOWRAMMODE Disk-streamed Global Signal Regression.
%
% This function performs GSR using a low-RAM, two-pass strategy:
%   1) First pass computes the dataset mean and global signal over time
%   2) Second pass regresses it out chunk-by-chunk
%
% The numerical intent matches GSR_standardMode.

% -------------------------------------------------------------------------
% Estimate data size and chunking
% -------------------------------------------------------------------------
dataBytes = prod([metaData.datSize, metaData.datLength, 4]);
nChunks = calculateMaxChunkSize(dataBytes, 3, .1);

Ny = metaData.datSize(1);
Nx = metaData.datSize(2);
Nt = metaData.datLength;

chunkSizePixels = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkSizePixels);

% -------------------------------------------------------------------------
% File handles
% -------------------------------------------------------------------------
fid_in = fopen(dataFile, 'r');
assert(fid_in ~= -1, 'Umitoolbox:GSR:FileOpenError', ...
    'Could not open input file "%s".', dataFile);
c_in = onCleanup(@() safeFclose(fid_in));

outFileName = fullfile(SaveFolder, 'GSR.dat');
preallocateDatFile(outFileName, [Ny, Nx, Nt], metaData.Datatype);

fid_out = fopen(outFileName, 'r+');
assert(fid_out ~= -1, 'Umitoolbox:GSR:FileOpenError', ...
    'Could not open output file "%s".', outFileName);
c_out = onCleanup(@() safeFclose(fid_out));

% -------------------------------------------------------------------------
% Prepare accumulators
% -------------------------------------------------------------------------
globalSum   = zeros(1, Nt, 'double');
globalCount = zeros(1, Nt, 'double');
dataSum     = 0;
dataCount   = 0;

h = waitbar(0, 'GSR: computing global signal (pass 1)...');
h.Name = 'GSR (pass 1/2)';
cleanupWaitbar = onCleanup(@() iCloseWaitbarSafely(h));

% =====================================================================
% PASS 1 - Compute dataset mean and global signal
% =====================================================================
for ii = 1:nChunks
    waitbar(ii / nChunks, h, 'GSR: computing global signal...');

    pxStart = (ii-1) * chunkSizePixels + 1;
    pxEnd   = min(ii * chunkSizePixels, Nx);
    idxX    = pxStart:pxEnd;

    % Read slab
    slab = spatialSlabIO('read', fid_in, Ny, Nx, Nt, idxX, metaData.Datatype);

    % Accumulate global mean from original data
    dataSum   = dataSum + sum(double(slab(:)), 'omitnan');
    dataCount = dataCount + sum(~isnan(slab(:)));

    % Identify invalid traces from first frame only
    idx_invalid_trace = isnan(slab(:,:,1));

    % Reshape slab to [pixels x time]
    slab = reshape(slab, [], Nt);
    idx_invalid_trace = idx_invalid_trace(:);

    % Build valid mask for global-signal estimation
    maskSlab = logical_mask(:, idxX);
    maskIdx  = maskSlab(:) & ~idx_invalid_trace;

    if any(maskIdx)
        validData = double(slab(maskIdx, :));
        globalSum   = globalSum + sum(validData, 1);
        globalCount = globalCount + sum(maskIdx);
        clear validData
    end

    clear slab idx_invalid_trace maskSlab maskIdx
end

if dataCount == 0
    if isgraphics(h), close(h); end
    error('Umitoolbox:GSR:InvalidInput', ...
        'Input data contain only NaN values.');
end

if ~any(globalCount > 0)
    if isgraphics(h), close(h); end
    error('Umitoolbox:GSR:InvalidInput', ...
        'Logical mask does not contain any valid pixels for GSR.');
end

mData = single(dataSum / dataCount);
Sig   = single(globalSum ./ globalCount);

sigMean = mean(Sig);
if ~isfinite(sigMean) || sigMean == 0
    if isgraphics(h), close(h); end
    error('Umitoolbox:GSR:InvalidInput', ...
        'Global signal mean is zero or invalid. Cannot normalize signal.');
end
Sig = Sig ./ sigMean;

% -------------------------------------------------------------------------
% Regression design matrix (constant + global signal)
% -------------------------------------------------------------------------
X = [ones(Nt,1,'single'), Sig(:)];
clear Sig

% =====================================================================
% PASS 2 - Regress global signal chunk-by-chunk
% =====================================================================
h.Name = 'GSR (pass 2/2)';

for ii = 1:nChunks
    waitbar(ii / nChunks, h, 'GSR: regressing signal...');

    pxStart = (ii-1) * chunkSizePixels + 1;
    pxEnd   = min(ii * chunkSizePixels, Nx);
    idxX    = pxStart:pxEnd;

    % Read slab
    slab = spatialSlabIO('read', fid_in, Ny, Nx, Nt, idxX, metaData.Datatype);

    % Identify invalid traces from first frame only
    idx_invalid_trace = isnan(slab(:,:,1));
    slabSz = size(slab);

    % Reshape slab to [pixels x time]
    slab = reshape(slab, [], Nt);
    idx_invalid_trace = idx_invalid_trace(:);

    % Replace invalid traces before regression
    slab(idx_invalid_trace, :) = 0;

    % Regression
    A = X * (X \ slab');
    slab = slab - A';

    % Restore mean
    slab = slab + mData;

    % Restore invalid traces
    slab(idx_invalid_trace, :) = NaN;

    % Reshape back and write corrected slab
    slab = reshape(slab, slabSz);
    spatialSlabIO('write', fid_out, Ny, Nx, Nt, idxX, metaData.Datatype, slab);

    clear slab A idx_invalid_trace
end

if isgraphics(h)
    close(h);
end
end

function iCloseWaitbarSafely(h)
%ICLOSEWAITBARSAFELY Close a waitbar handle if it still exists.

if isgraphics(h)
    close(h);
end
end
