function outData = GSR(data, SaveFolder, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GSR  Perform Global Signal Regression on image time series.
%
% This function removes global fluctuations from fluorescence or imaging
% data by regressing out the global mean signal.
%
% The function supports TWO execution modes:
%
%   1) STANDARD MODE (in-memory)
%      - Triggered when "data" is a numeric array
%      - Entire dataset is processed in RAM
%
%   2) LOW-RAM MODE (streaming)
%      - Triggered when "data" is a filename (.dat)
%      - Data are processed in spatial chunks directly from disk
%
% -------------------------------------------------------------------------
% Inputs:
%
%   data :
%       Either:
%         - Numeric array [Y, X, T]              -> STANDARD MODE
%         - Char array pointing to a .dat file   -> LOW-RAM MODE
%
%   metaData :
%       Metadata associated with the data.
%       Must be a struct or matlab.io.MatFile.
%
%   opts (optional) :
%       Structure with optional parameters:
%         - b_UseMask (logical):
%               If true, GSR is computed only inside a logical mask.
%         - MaskFile (char):
%               .mat file containing variable "logical_mask".
%
%   object (optional) :
%       Acquisition or Modality object (used by PipelineManager only).
%
% -------------------------------------------------------------------------
% Output:
%
%   outData :
%       - STANDARD MODE: corrected data array [Y, X, T]
%       - LOW-RAM MODE : filename of corrected .dat file
%
% -------------------------------------------------------------------------
% Notes:
%
%   Invalid traces are identified from the FIRST FRAME only. A NaN in
%   frame 1 indicates that the whole pixel trace is invalid across time.
%   This avoids allocating a full logical array of size(data), which can
%   be too large for big datasets.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults (PipelineManager compatibility)
default_Output = 'GSR.dat'; %#ok
default_opts   = struct('b_UseMask', false);
opts_values    = struct('b_UseMask',[true,false]); %#ok

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x)==3) || ischar(x));
addRequired(p,'SaveFolder',@isfolder);
addOptional(p,'opts',default_opts,@(x) isstruct(x));

parse(p,data,SaveFolder,varargin{:});

opts = p.Results.opts;

% -------------------------------------------------------------------------
% Load or build logical mask
% -------------------------------------------------------------------------
if opts.b_UseMask    
    a = load(SaveFolder,'logical_mask');
    if isempty(fieldnames(a))
        error('umIToolbox:GSR:MissingInput',...
              'Variable "logical_mask" not found in mask file.');
    end

    logical_mask = a.logical_mask;
    assert(isequal(size(logical_mask),metaData.datSize),...
        'umIToolbox:GSR:InvalidInput',...
        'Logical mask size does not match data frame size.');
else    
    logical_mask = true(size(data,[1 2])); % TODO: This will fail for data as filename. TO BE CHANGED.    
end

% -------------------------------------------------------------------------
% Dispatch execution mode
% -------------------------------------------------------------------------
if ischar(data)
    % =========================
    % LOW-RAM MODE
    % =========================
    outData = GSR_lowRAMmode(data, metaData, logical_mask);
else
    % =========================
    % STANDARD MODE
    % =========================
    outData = GSR_standardMode(data, logical_mask);
end

disp('Finished GSR.');
end


%--------------------------------------------------------------------------
% Local functions
%--------------------------------------------------------------------------
function outData = GSR_standardMode(data, logical_mask)
% GSR_STANDARDMODE  In-memory Global Signal Regression.
%
% Performs GSR assuming the entire dataset fits in RAM.

szData = size(data);

% Identify invalid traces from the first frame only
idx_invalid_trace = isnan(data(:,:,1));

% Compute global mean from original data before zero-filling invalid traces
mData = mean(data,'all','omitnan');

% Reshape to [pixels x time]
data = reshape(data,[],szData(3));
idx_invalid_trace = idx_invalid_trace(:);

% Build valid mask for global-signal estimation
maskIdx = logical_mask(:) & ~idx_invalid_trace;
assert(any(maskIdx), ...
    'umIToolbox:GSR:InvalidInput', ...
    'Logical mask does not contain any valid pixels for GSR.');

% Replace invalid traces with zeros before regression
data(idx_invalid_trace,:) = 0;

% Compute global signal
disp('Calculating Global Signal Regression...');
Sig = mean(data(maskIdx,:),1);

sigMean = mean(Sig);
assert(isfinite(sigMean) && sigMean ~= 0, ...
    'umIToolbox:GSR:InvalidInput', ...
    'Global signal mean is zero or invalid. Cannot normalize signal.');
Sig = Sig ./ sigMean;

% Regression
X = [ones(szData(3),1), Sig(:)];
A = X * (X \ data');
data = data - A';

% Restore mean
data = data + mData;

% Restore invalid traces
data(idx_invalid_trace,:) = NaN;

% Reshape back
outData = reshape(data,szData);
end


function outFileName = GSR_lowRAMmode(fluoFile, metaData, logical_mask)
% GSR_LOWRAMMODE  Disk-streamed Global Signal Regression.
%
% This function performs GSR using a low-RAM, two-pass strategy:
%   1) First pass computes the dataset mean and global signal over time
%   2) Second pass regresses it out chunk-by-chunk
%
% The numerical intent matches GSR_standardMode.

% -------------------------------------------------------------------------
% Estimate data size and chunking
% -------------------------------------------------------------------------
dataBytes = prod([metaData.datSize, metaData.datLength,...
                  getByteSize(metaData.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes,3,.1);

Ny = metaData.datSize(1);
Nx = metaData.datSize(2);
Nt = metaData.datLength;

chunkSizePixels = ceil(Nx / nChunks);

% -------------------------------------------------------------------------
% File handles
% -------------------------------------------------------------------------
fid_in = fopen(fluoFile,'r');
assert(fid_in ~= -1, 'umIToolbox:GSR:FileOpenError', ...
    'Could not open input file "%s".', fluoFile);
c_in = onCleanup(@() safeFclose(fid_in)); %#ok<NASGU>

outFileName = 'GSRCORRECTED_DATA.dat';
preallocateDatFile(outFileName,metaData);

fid_out = fopen(outFileName,'r+');
assert(fid_out ~= -1, 'umIToolbox:GSR:FileOpenError', ...
    'Could not open output file "%s".', outFileName);
c_out = onCleanup(@() safeFclose(fid_out)); %#ok<NASGU>

% -------------------------------------------------------------------------
% Prepare accumulators
% -------------------------------------------------------------------------
globalSum   = zeros(1,Nt,'double');
globalCount = zeros(1,Nt,'double');
dataSum     = 0;
dataCount   = 0;

h = waitbar(0,'GSR: computing global signal (pass 1)...');
h.Name = 'GSR (pass 1/2)';

% =====================================================================
% PASS 1 - Compute dataset mean and global signal
% =====================================================================
for ii = 1:nChunks
    waitbar(ii/nChunks,h,'GSR: computing global signal...');

    pxStart = (ii-1)*chunkSizePixels + 1;
    pxEnd   = min(ii*chunkSizePixels,Nx);
    idxX    = pxStart:pxEnd;

    % Read slab
    slab = spatialSlabIO('read',fid_in,Ny,Nx,Nt,idxX,metaData.Datatype);

    % Accumulate global mean from original data
    dataSum   = dataSum + sum(double(slab(:)),'omitnan');
    dataCount = dataCount + sum(~isnan(slab(:)));

    % Identify invalid traces from first frame only
    idx_invalid_trace = isnan(slab(:,:,1));

    % Reshape slab to [pixels x time]
    slab = reshape(slab,[],Nt);
    idx_invalid_trace = idx_invalid_trace(:);

    % Build valid mask for global-signal estimation
    maskSlab = logical_mask(:,idxX);
    maskIdx  = maskSlab(:) & ~idx_invalid_trace;

    if any(maskIdx)
        validData = double(slab(maskIdx,:));
        globalSum   = globalSum + sum(validData,1);
        globalCount = globalCount + sum(maskIdx);
        clear validData
    end

    clear slab idx_invalid_trace maskSlab maskIdx
end

if dataCount == 0
    if isgraphics(h), close(h); end
    error('umIToolbox:GSR:InvalidInput', ...
          'Input data contain only NaN values.');
end

if ~any(globalCount > 0)
    if isgraphics(h), close(h); end
    error('umIToolbox:GSR:InvalidInput', ...
          'Logical mask does not contain any valid pixels for GSR.');
end

mData = dataSum / dataCount;
Sig   = globalSum ./ globalCount;

sigMean = mean(Sig);
if ~isfinite(sigMean) || sigMean == 0
    if isgraphics(h), close(h); end
    error('umIToolbox:GSR:InvalidInput', ...
          'Global signal mean is zero or invalid. Cannot normalize signal.');
end
Sig = Sig ./ sigMean;

% -------------------------------------------------------------------------
% Regression design matrix (constant + global signal)
% -------------------------------------------------------------------------
X = [ones(Nt,1), Sig(:)];
clear Sig

% =====================================================================
% PASS 2 - Regress global signal chunk-by-chunk
% =====================================================================
h.Name = 'GSR (pass 2/2)';

for ii = 1:nChunks
    waitbar(ii/nChunks,h,'GSR: regressing signal...');

    pxStart = (ii-1)*chunkSizePixels + 1;
    pxEnd   = min(ii*chunkSizePixels,Nx);
    idxX    = pxStart:pxEnd;

    % Read slab
    slab = spatialSlabIO('read',fid_in,Ny,Nx,Nt,idxX,metaData.Datatype);

    % Identify invalid traces from first frame only
    idx_invalid_trace = isnan(slab(:,:,1));
    slabSz = size(slab);

    % Reshape slab to [pixels x time]
    slab = reshape(slab,[],Nt);
    idx_invalid_trace = idx_invalid_trace(:);

    % Replace invalid traces before regression
    slab(idx_invalid_trace,:) = 0;

    % Regression
    A = X * (X \ slab');
    slab = slab - A';

    % Restore mean
    slab = slab + mData;

    % Restore invalid traces
    slab(idx_invalid_trace,:) = NaN;

    % Reshape back and write corrected slab
    slab = reshape(slab,slabSz);
    spatialSlabIO('write',fid_out,Ny,Nx,Nt,idxX,metaData.Datatype,slab);

    clear slab A idx_invalid_trace
end

if isgraphics(h)
    close(h);
end
end