function outData = spatialGaussFilt(data, metaData, varargin)
% SPATIALGAUSSFILT performs a spatial Gaussian filter on an image time series.
%
% Limitations:
%   The data must be an image time series with dimensions {Y,X,T}.
%
% Inputs:
%   data: numerical matrix or filename of a .dat file containing image time series (Y×X×T)
%   metaData: .mat file or struct with metadata
%   opts (optional): structure with extra parameters
%       - Sigma: Gaussian sigma (positive scalar)
%
% Outputs:
%   outData: filtered data (if input is array). For file input, a new .dat file is created.

% Defaults:
default_Output = 'spatFilt.dat'; %#ok Pipeline management
default_opts = struct('Sigma',1);
opts_values = struct('Sigma',[0,Inf]); %#ok Pipeline management

% Input parsing
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) || ischar(x));
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') || isstruct(x));
addOptional(p,'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p, data, metaData, varargin{:});

outData = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p

% Validate Sigma
errID = 'umIToolbox:spatialGaussFilt:InvalidInput';
assert(opts.Sigma>0, errID,'Sigma must be a positive value!');

%% 
if ischar(outData)
    % --- File mode ---
    inFile = outData;
    outFile = fullfile(fileparts(inFile), 'DATASPATIALFILTERED.dat');

    fprintf('Spatial Gaussian filtering on file: %s\n', inFile);

    % Open file
    fidIn = fopen(inFile,'r');
    cIn = onCleanup(@() safeFclose(fidIn));
    if fidIn == -1
        error('Cannot open input file: %s', inFile);
    end

    % Metadata
    Ny = metaData.datSize(1);
    Nx = metaData.datSize(2);
    Nt = metaData.datLength;

    % Preallocate output .dat file
    preallocateDatFile(outFile, metaData);
    fidOut = fopen(outFile,'r+');
    cOut = onCleanup(@() safeFclose(fidOut));
    if fidOut == -1
        fclose(fidIn);
        error('Cannot create output file: %s', outFile);
    end

    % Calculate optimal chunk size along time
    frameBytes = Ny*Nx*getByteSize(metaData.Datatype); % single
    totalBytes = frameBytes*Nt;
    nChunks = calculateMaxChunkSize(totalBytes,2,.1);
    chunkFrames = ceil(Nt/nChunks);

    fprintf('Filtering %d frames in %d chunk(s)...\n', Nt, nChunks);
    

    for c = 1:nChunks
        
        tStart = (c-1)*chunkFrames + 1;
        tEnd = min(tStart + chunkFrames - 1, Nt);
        nThisChunk = tEnd - tStart + 1;

        % Move to start of chunk
        fseek(fidIn, (tStart-1)*frameBytes, 'bof');
        fprintf('Chunk #%i [Reading file ...]\n',c)
        % Read chunk and reshape
        slab = fread(fidIn, [Nx*Ny, nThisChunk], '*single');
        slab = reshape(slab, Ny, Nx, nThisChunk);
        
        fprintf('Chunk #%i [Filtering...]\n',c)
        % Apply Gaussian filter per frame
        slab = imgaussfilt(slab, opts.Sigma, 'FilterDomain','spatial');
           
        % Write chunk back to output
        fseek(fidOut, (tStart-1)*frameBytes, 'bof');
        
        fprintf('Chunk #%i [Writing to file...]\n',c)
        fwrite(fidOut, slab, 'single'); % transpose back
        fprintf('Chunk #%i [Completed]\n',c)
        clear slab
        
    end

    fclose(fidIn);
    fclose(fidOut);
    fprintf('Finished spatial Gaussian filtering.\n');

    outData = outFile;

else
    % --- Array mode ---
    errMsg = 'Data must be a 3D matrix with dimensions {Y,X,T}.';
    assert(ndims(outData)==3, errID, errMsg);

    idx_nan = isnan(outData);
    outData(idx_nan) = 0;

    fprintf('Filtering %d frames in memory...\n', size(outData,3));
    outData = imgaussfilt(outData,opts.Sigma, 'FilterDomain','spatial');
%     w = waitbar(0,'Filtering frames...');
%     for t = 1:size(outData,3)
%         waitbar(t/size(outData,3), w);
%         outData(:,:,t) = imgaussfilt(outData(:,:,t), opts.Sigma, 'FilterDomain','spatial');
%     end
%     close(w);

    % Restore NaNs
    outData(idx_nan) = NaN;
    fprintf('Finished spatial Gaussian filtering.\n');
end

end
