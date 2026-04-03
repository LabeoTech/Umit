function preallocateDatFile(filename, metaData)
% PREALLOCATEDATFILE  Safely preallocate a raw binary file on disk.
%
% This function writes zeros to disk to guarantee correct allocation
% and later random-access writes. The chunk size is dynamically
% computed to respect system RAM using calculateMaxChunkSize.

    % Determine file dimensions and precision
    dims = [metaData.datSize metaData.datLength];
    if isfield(metaData,'Datatype')    
        precision = metaData.Datatype;
    else
        precision = 'single';  % default
    end

    bytesPerElement = getByteSize(precision);
    nElements = prod(dims);

    % Estimate total required memory in bytes
    totalBytes = nElements * bytesPerElement;

    % Compute optimal number of chunks (RAM safe)
    nChunks = calculateMaxChunkSize(totalBytes, 1, 0.2);  % 20% RAM reserved
    chunkElements = ceil(nElements / nChunks);

    fid = fopen(filename, 'w');
    cid = onCleanup(@() safeFclose(fid));
    if fid < 0
        error('Cannot open file for writing: %s', filename);
    end

    disp('Preallocating .dat file...')

    % Preallocate chunks
    for k = 1:nChunks
        % Determine actual chunk size (last chunk may be smaller)
        startIdx = (k-1)*chunkElements + 1;
        endIdx   = min(k*chunkElements, nElements);
        currentChunkSize = endIdx - startIdx + 1;

        fwrite(fid, zeros(currentChunkSize, 1, precision), precision);
        
    end

    fclose(fid);
    disp('Empty .dat file created.')

    % Save meta data file
    [filepath, name, ext] = fileparts(filename);
    if isa(metaData,'matlab.io.MatFile')
        metaData.Properties.Writable = true;
    end
    metaData.datFile = [name ext];
    save(fullfile(filepath, [name, '.mat']), '-struct','metaData')
end
