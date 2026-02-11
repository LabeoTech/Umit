function outData = apply_detrend(data, metaData)
% APPLY_DETREND applies a linear detrend along the time dimension of image
% time series or event-split image time series.
%
% Inputs:
%   data: 3D or 4D numerical matrix (Y,X,T) or (E,Y,X,T), or filename of .dat file
%   metaData: struct or matfile containing metadata
%
% Outputs:
%   outData: detrended data (or .dat filename if input was a file)

default_Output = 'data_detrended.dat'; %#ok

% Input validation
errID = 'umIToolbox:apply_detrend:WrongInput';
if isnumeric(data)
    dims = ndims(data);
    assert(dims==3 || dims==4, errID, 'Data must be 3D or 4D.');
    outData = data;
    bIsFile = false;
else
    assert(ischar(data), errID, 'Data must be numeric or a filename.');
    inFile = data;
    bIsFile = true;
end

% Metadata and dimensions
dataSize = [metaData.datSize, metaData.datLength];
bHasEvents = any(strcmpi('E', metaData.dim_names));
if bHasEvents
    Ne = dataSize(1); Ny = dataSize(2); Nx = dataSize(3); Nt = dataSize(4);
else
    Ne = 1; Ny = dataSize(1); Nx = dataSize(2); Nt = dataSize(3);
end

% Determine baseline frames
if isfield(metaData, 'preEventTime_sec')
    frames = round(metaData.preEventTime_sec*metaData.Freq);
    frames = max(frames,3);
    if mod(frames,2)==0, frames = frames+1; end
else
    frames = 7;
end

if bIsFile
    outFile = fullfile(fileparts(inFile), 'DATADETRENDED.dat');
    preallocateDatFile(outFile, metaData);

    fidIn  = fopen(inFile,'r');
    cIn = onCleanup(@() safeFclose(fidIn));
    fidOut = fopen(outFile,'r+');
    cOut = onCleanup(@() safeFclose(fidOut));
    
    if bHasEvents
        
        % --- Event-split data: fread trial by trial ---
        
        for e = 1:Ne            
            % Read trial
            fprintf('Chunk %i/%i [Reading file ...]\n',e,Ne)
            trialData = readTrial(fidIn,e,[Ne,Ny,Nx,Nt],'single');
            fprintf('Chunk %i/%i [Detrending data ...]\n',e,Ne)
            trialData = reshape(trialData, Ny*Nx, Nt);
            % Detrending
            delta_y = median(trialData(:,end-frames+1:end),2,'omitnan') - ...
                      median(trialData(:,1:frames),2,'omitnan');
            delta_x = Nt - frames;
            M = delta_y ./ delta_x;
            b = median(trialData(:,1:frames),2,'omitnan');
            
            trend = M .* linspace(-2,Nt-3,Nt) + b;
            
            trialData = trialData - trend + b;            
            clear delta_x delta_y M b trend
            
            trialData = reshape(trialData, Ny, Nx, Nt);
            
            % Write to outFile
            fprintf('Chunk %i/%i [Writing to file ...]\n',e,Ne)
            writeTrial_YXTE(fidOut,e,trialData,[Ny,Nx,Nt,Ne],'single');
            clear trialData
        end
        fclose(fidOut);
        permuteDat_YXTE_to_EYXT_inplace(outFile,[Ny,Nx,Nt,Ne],'single')
        fprintf('Chunk %i/%i [Completed]\n',e,Ne)
        
    else
        % --- Non-event data: chunk along X dimension ---
        nChunks = calculateMaxChunkSize(Nx*Ny*Nt*4,2,.3);
        chunkX  = ceil(Nx / nChunks);
        
        for c = 1:nChunks
            xStart = (c-1)*chunkX + 1;
            xEnd   = min(xStart + chunkX - 1, Nx);
            xIdx   = xStart:xEnd;
            fprintf('Chunk %i/%i [Reading file ...]\n',c,nChunks)
            slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, metaData.Datatype);
            fprintf('Chunk %i/%i [Detrending data ...]\n',c,nChunks)
            slab = reshape(slab, Ny*numel(xIdx), Nt);
            % Detrending
            delta_y = median(slab(:,end-frames+1:end),2,'omitnan') - ...
                      median(slab(:,1:frames),2,'omitnan');
            delta_x = Nt - frames;
            M = delta_y ./ delta_x;
            b = median(slab(:,1:frames),2,'omitnan');
            
            trend = M .* linspace(-2,Nt-3,Nt) + b;
            
            slab = slab - trend + b;
            clear delta_x delta_y M b trend
            slab = reshape(slab, Ny, numel(xIdx), Nt);     
            
            % Write to file
            fprintf('Chunk %i/%i [Writing to file ...]\n',c,nChunks)
            spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, metaData.Datatype, slab);
            fprintf('Chunk %i/%i [Completed]\n',c,nChunks)
            clear slab
        end
        fclose(fidOut);
    end

    fclose(fidIn);    
    outData = outFile;

else
    % --- Array mode ---
    orig_sz = size(outData);

    if bHasEvents
        for e = 1:orig_sz(1)
            slab = squeeze(outData(e,:,:,:));
            slab = reshape(slab, Ny*Nx, Nt);
            delta_y = median(slab(:,end-frames+1:end),2,'omitnan') - ...
                      median(slab(:,1:frames),2,'omitnan');
            delta_x = Nt - frames;
            M = delta_y ./ delta_x;
            b = median(slab(:,1:frames),2,'omitnan');
            trend = M .* linspace(-2,Nt-3,Nt) + b;

            slab = slab - trend + b;
            outData(e,:,:,:) = reshape(slab, Ny, Nx, Nt);
        end
    else
        slab = reshape(outData, Ny*Nx, Nt);
        delta_y = median(slab(:,end-frames+1:end),2,'omitnan') - ...
                  median(slab(:,1:frames),2,'omitnan');
        delta_x = Nt - frames;
        M = delta_y ./ delta_x;
        b = median(slab(:,1:frames),2,'omitnan');
        trend = M .* linspace(-2,Nt-3,Nt) + b;

        slab = slab - trend + b;
        outData = reshape(slab, orig_sz);
    end
end

disp('Finished detrend!');
end
