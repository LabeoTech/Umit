function outData = GSR(data, metaData, varargin)
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
%         - Numeric array [Y, X, T]              ? STANDARD MODE
%         - Char array pointing to a .dat file   ? LOW-RAM MODE
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults (PipelineManager compatibility)
default_Output = 'GSR.dat'; %#ok
default_opts   = struct('b_UseMask', false, 'MaskFile','ImagingReferenceFrame.mat');
opts_values    = struct('b_UseMask',[true,false],'MaskFile',{{'ImagingReferenceFrame.mat'}}); %#ok
default_object = '';

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x)==3) || ischar(x));
addRequired(p,'metaData',@(x) isstruct(x) || isa(x,'matlab.io.MatFile'));
addOptional(p,'opts',default_opts,@(x) isstruct(x));
addOptional(p,'object',default_object,...
    @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality'));

parse(p,data,metaData,varargin{:});

data     = p.Results.data;
metaData = p.Results.metaData;
opts     = p.Results.opts;
object   = p.Results.object;
clear p

% -------------------------------------------------------------------------
% Validate metadata
% -------------------------------------------------------------------------
errID  = 'umIToolbox:GSR:InvalidInput';
errMsg = 'Data must be an image time series with dimensions {Y,X,T}.';
assert(all(ismember(metaData.dim_names,{'Y','X','T'})),errID,errMsg);

% -------------------------------------------------------------------------
% Load or build logical mask
% -------------------------------------------------------------------------
if opts.b_UseMask
    opts.MaskFile = findMyROIfile(opts.MaskFile,object);
    a = load(opts.MaskFile,'logical_mask');
    if isempty(fieldnames(a))
        error('umIToolbox:GSR:MissingInput',...
              'Variable "logical_mask" not found in mask file.');
    end
    logical_mask = a.logical_mask;
    assert(isequal(size(logical_mask),metaData.datSize),...
        'umIToolbox:GSR:InvalidInput',...
        'Logical mask size does not match data frame size.');
else
    logical_mask = true(metaData.datSize);
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

% Replace NaNs
idx_nan = isnan(data);
data(idx_nan) = 0;

% Reshape
szData = size(data);
mData  = mean(data,'all');
data   = reshape(data,[],szData(3));

% Compute global signal
disp('Calculating Global Signal Regression...');
Sig = mean(data(logical_mask(:),:),1);
Sig = Sig ./ mean(Sig);

% Regression
X = [ones(szData(3),1), Sig'];
A = X * (X \ data');
data = data - A';

% Reshape back
outData = reshape(data,szData);
outData = outData + mData;

% Restore NaNs
outData(idx_nan) = NaN;


end


function outFileName = GSR_lowRAMmode(fluoFile, metaData, logical_mask)
% GSR_LOWRAMMODE  Disk-streamed Global Signal Regression.
%
% This function performs GSR using a low-RAM, two-pass strategy:
%   1) First pass computes the global signal over time
%   2) Second pass regresses it out chunk-by-chunk
%
% The numerical result is identical to GSR_standardMode.

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
c_in = onCleanup(@() safeFclose(fid_in));
outFileName = 'GSRCORRECTED_DATA.dat';
preallocateDatFile(outFileName,metaData);
fid_out = fopen(outFileName,'r+');
c_out = onCleanup(@() safeFclose(fid_out));

% -------------------------------------------------------------------------
% Prepare global signal accumulator
% -------------------------------------------------------------------------
globalSum   = zeros(1,Nt,'double');
globalCount = 0;

h = waitbar(0,'GSR: computing global signal (pass 1)...');
h.Name = 'GSR (pass 1/2)';
% =====================================================================
% PASS 1 — Compute global signal
% =====================================================================
for ii = 1:nChunks
    waitbar(ii/nChunks,h,'GSR: computing global signal')
    pxStart = (ii-1)*chunkSizePixels + 1;
    pxEnd   = min(ii*chunkSizePixels,Nx);
    idxX    = pxStart:pxEnd;

    % Read slab
    waitbar(ii/nChunks,h, 'GSR: reading file...',ii,nChunks);
    slab = spatialSlabIO('read',fid_in,Ny,Nx,Nt,idxX,metaData.Datatype);

    % Reshape to [pixels x time]
    slab = reshape(slab,[],Nt);

    % Mask pixels
    maskSlab = logical_mask(:,idxX);
    maskIdx  = maskSlab(:);
    waitbar(ii/nChunks,h,'GSR: computing global signal')
    if any(maskIdx)
        globalSum   = globalSum + sum(double(slab(maskIdx,:)),1);
        globalCount = globalCount + sum(maskIdx);
    end
    clear slab
end

% Normalize global signal
Sig = globalSum ./ globalCount;
mData = mean(Sig);
Sig = Sig ./ mean(Sig);

% -------------------------------------------------------------------------
% Regression design matrix (constant + global signal)
% -------------------------------------------------------------------------
X = [ones(Nt,1), Sig(:)];
clear Sig
% =====================================================================
% PASS 2 — Regress global signal chunk-by-chunk
% =====================================================================
h.Name = 'GSR (pass 2/2)';

for ii = 1:nChunks
    waitbar(ii/nChunks,h, 'GSR: regressing signal',ii,nChunks);
    
    pxStart = (ii-1)*chunkSizePixels + 1;
    pxEnd   = min(ii*chunkSizePixels,Nx);
    idxX    = pxStart:pxEnd;

    % Read slab
    waitbar(ii/nChunks,h, 'GSR: reading file...',ii,nChunks);
    slab = spatialSlabIO('read',fid_in,Ny,Nx,Nt,idxX,metaData.Datatype);

    slabSz = size(slab);
    slab   = reshape(slab,[],Nt);

    % Replace NaNs
    idx_nan = isnan(slab);
    slab(idx_nan) = 0;
    waitbar(ii/nChunks,h, 'GSR: regressing signal...',ii,nChunks);    
    % Regression
    A = X * (X \ slab');
    slab = slab - A';
    clear A
    % Restore mean
    slab = slab + mData;

    % Restore NaNs
    slab(idx_nan) = NaN;
    clear idx_nan
    % Reshape back
    slab = reshape(slab,slabSz);
    waitbar(ii/nChunks,h, 'GSR: writing to file...',ii,nChunks);
    % Write corrected slab
    spatialSlabIO('write',fid_out,Ny,Nx,Nt,idxX,metaData.Datatype,slab);
    clear slab
    waitbar(ii/nChunks,h)
end

close(h);

% -------------------------------------------------------------------------
% Cleanup
% -------------------------------------------------------------------------
fclose(fid_in);
fclose(fid_out);

end



