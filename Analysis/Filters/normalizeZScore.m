function outData = normalizeZScore(data, metaData)
% NORMALIZEZSCORE normalizes an image time series with zero mean and unit
% standard deviation (z-score).
%
% Supports two execution modes:
%   1) STANDARD MODE (in-memory)
%      - Triggered when "data" is a numeric array
%   2) LOW-RAM MODE (file-backed)
%      - Triggered when "data" is a .dat filename
%
% Limitations:
%   The data must be an Image time series with dimensions {Y,X,T}.
%   Event-split data are not supported.
%
% Inputs:
%   data     : numeric array (Y,X,T) or a .dat filename
%   metaData : structure or matlab.io.MatFile with associated metadata
%
% Output:
%   outData  : normalized array (Y,X,T) or output filename in low-RAM mode

% Defaults:
default_Output = 'normZ.dat'; %#ok This line is here just for Pipeline management.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) (isnumeric(x) && ndims(x) == 3) || ischar(x));
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') || isstruct(x));
parse(p,data, metaData);

data     = p.Results.data;
metaData = p.Results.metaData;
clear p

% Validate dimensions
errID  = 'umIToolbox:normalizeZScore:InvalidInput';
errMsg = 'Input data must be an image time series with dimensions "Y","X","T".';
assert(all(ismember({'Y','X','T'}, metaData.dim_names)), errID, errMsg);
assert(~any(strcmpi(metaData.dim_names,'E')), errID, ...
    'normalizeZScore does not support event-split data.');

disp('Calculating z-score normalization...');

% -------------------------------------------------------------------------
% Low-RAM / file-backed mode
% -------------------------------------------------------------------------
if ischar(data)

    inFile  = data;
    outData = fullfile(fileparts(inFile), 'NORMALIZE_ZSCORE.dat');

    Ny = metaData.datSize(1);
    Nx = metaData.datSize(2);
    Nt = metaData.datLength;

    % Preallocate output file
    preallocateDatFile(outData, metaData);

    fidIn  = fopen(inFile,  'r');
    cIn = onCleanup(@() safeFclose(fidIn)); %#ok<NASGU>
    fidOut = fopen(outData, 'r+');
    cOut = onCleanup(@() safeFclose(fidOut)); %#ok<NASGU>

    % Chunking along X
    dataBytes = prod(metaData.datSize) * getByteSize('single');
    nChunks   = calculateMaxChunkSize(dataBytes, 2);
    chunkX    = ceil(Nx / nChunks);

    for c = 1:nChunks
        xStart = (c-1)*chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, Nx);
        xIdx   = xStart:xEnd;

        % Read slab: [Y, Xc, T]
        slab = spatialSlabIO( ...
            'read', fidIn, Ny, Nx, Nt, xIdx, 'single');

        % Reshape to pixels x time
        slab2D = reshape(slab, [], Nt);

        % Z-score across time
        mu  = mean(slab2D, 2, 'omitnan');
        sig = std(slab2D, 0, 2, 'omitnan');
        sig(sig == 0) = 1; % Avoid division by zero for constant traces

        slab2D = (slab2D - mu) ./ sig;
        slab   = reshape(slab2D, size(slab));

        % Write back to output file
        spatialSlabIO( ...
            'write', fidOut, Ny, Nx, Nt, xIdx, 'single', slab);
    end

    fclose(fidIn);
    fclose(fidOut);

% -------------------------------------------------------------------------
% In-memory mode
% -------------------------------------------------------------------------
else
    orig_sz = size(data);
    data2D = reshape(data, [], orig_sz(3));

    mu  = mean(data2D, 2, 'omitnan');
    sig = std(data2D, 0, 2, 'omitnan');
    sig(sig == 0) = 1; % Avoid division by zero for constant traces

    data2D = (data2D - mu) ./ sig;
    outData = reshape(data2D, orig_sz);
end

disp('Finished normalizeZScore.');
end