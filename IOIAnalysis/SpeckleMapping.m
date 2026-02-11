function DatOut = SpeckleMapping(folderPath, sType, channel, bSaveMap, bLogScale, bRAMSafeMode)
% SPECKLEMAPPING computes speckle contrast maps from a 3D acquisition.
% Supports low-RAM processing with bRAMSafeMode flag.
%
% INPUTS:
%   folderPath: folder containing speckle.dat and .mat metadata
%   sType: 'spatial' or 'temporal'
%   channel: optional, channel to process ('speckle' default)
%   bSaveMap: boolean, save a .tiff map
%   bLogScale: boolean, apply -log10 scaling
%   bRAMSafeMode: boolean, enable low RAM hybrid processing
%
% OUTPUT:
%   DatOut: speckle map (or filename if bRAMSafeMode is true)

if nargin < 3, channel='speckle'; end
if nargin < 4, bSaveMap=1; end
if nargin < 5, bLogScale=1; end
if nargin < 6, bRAMSafeMode=0; end

if ~strcmp(folderPath(end), filesep)
    folderPath = [folderPath filesep];
end

channel = lower(channel);
datFile = [folderPath channel '.dat'];
matFile = [folderPath channel '.mat'];

if ~exist(datFile,'file')
    error('%s file is missing.', datFile);
end

Infos = load(matFile);
Ny = Infos.datSize(1);
Nx = Infos.datSize(2);
Nt = Infos.datLength;

md = Infos;
md.datLength = 1;
md.Freq = 0;
md.dim_names = {'Y','X'};

if bRAMSafeMode
    % --- Low RAM mode: hybrid chunked approach ---
    outFilename = [folderPath 'SPECKLEMAP.dat'];
    frameOut = zeros(Ny,Nx,Infos.Datatype);
    mData = zeros(Ny,Nx,Infos.Datatype);
    fidIn = fopen(datFile,'r');
    cIn = onCleanup(@() safeFclose(fidIn));
    disp('Pass 1 - Calculating temporal mean...')
    for t = 1:Nt
        fseek(fidIn, (t-1)*Ny*Nx*getByteSize('single'),'bof');
        frame = fread(fidIn, Ny*Nx, ['*' Infos.Datatype]);
        frame = reshape(frame, Ny, Nx);
        mData = mData + frame;        
    end
    mData = mData/Nt;
    disp('Pass 2 - Calculating Speckle Contrast ...')
    switch lower(sType)
        case 'spatial'
            % Iterate over time using fread + fseek                       
            Kernel = single(fspecial('disk',2)>0);                
            
            for t = 1:Nt
                fseek(fidIn, (t-1)*Ny*Nx*getByteSize('single'),'bof');
                frame = fread(fidIn, Ny*Nx, ['*' Infos.Datatype]);
                frame = reshape(frame, Ny, Nx);
                frame = frame./mData;
                               
                frame = stdfilt(frame, Kernel);
                frame = remOutlier(frame);
                
                
                frameOut = frameOut + frame;
                
            end
            frameOut = frameOut./Nt;            
            
        case 'temporal'
            % Iterate over X chunks using spatialSlabIO                        
            Kernel = ones(1,1,5,'single');
            
            nChunks = calculateMaxChunkSize(Nx*Ny*Nt*getByteSize(Infos.Datatype), 2);
            chunkX = ceil(Nx / nChunks);
            
            for c = 1:nChunks
                xStart = (c-1)*chunkX + 1;
                xEnd   = min(xStart + chunkX - 1, Nx);
                xIdx   = xStart:xEnd;
                
                slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, Infos.Datatype);
                slab = slab./mean(slab,3);
                slab = stdfilt(slab, Kernel);
                slab = remOutlier(slab);               
                frameOut = frameOut + mean(slab,3);                
            end
            
        otherwise
            error('Invalid sType. Use "spatial" or "temporal".');
    end
    if bLogScale
        frameOut = -log10(frame);
    end
    % Close file handle
    fclose(fidIn);
    % Save Speckle map to file
    md = genMetaData(frameOut,{'Y','X'},Infos);
    md.datFile = outFilename;
    save2Dat(outFilename,single(frameOut),md);          
    
    DatOut = outFilename; % return filename for low RAM
    
else
    % --- Standard execution: full RAM loading (original approach) ---
    fid = fopen(datFile,'r');
    dat = fread(fid, inf, '*single');
    fclose(fid);
    dat = reshape(dat, Ny, Nx, Nt);
    dat = dat ./ mean(dat,3);

    disp('Mapping computation...');
    switch lower(sType)
        case 'spatial'
            Kernel = single(fspecial('disk',2)>0);
        case 'temporal'
            Kernel = ones(1,1,5,'single');
        otherwise
            error('Invalid sType. Use "spatial" or "temporal".');
    end

    DatOut = stdfilt(dat, Kernel);
    DatOut = remOutlier(DatOut);
    % Calculate time average of Speckle Contrast values
    DatOut = mean(DatOut, 3);
    if bLogScale
        DatOut = -log10(DatOut);
    end
    DatOut =  single(DatOut);
end

% Save TIFF map if requested
if bSaveMap
    obj = Tiff(fullfile(folderPath,'std_speckle.tiff'),'w');
    setTag(obj,'ImageWidth',size(DatOut,2));
    setTag(obj,'ImageLength',size(DatOut,1));
    setTag(obj,'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(obj,'SampleFormat',Tiff.SampleFormat.IEEEFP);
    setTag(obj,'BitsPerSample',32);
    setTag(obj,'SamplesPerPixel',1);
    setTag(obj,'Compression',Tiff.Compression.None);
    setTag(obj,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    write(obj, DatOut);
    close(obj);
end

disp('Done');
end

function data = remOutlier(data)
% Replace values above 99th percentile with 99th percentile
dataSort = data(~isnan(data));
dataSort = sort(dataSort(:));
idx99 = floor(0.99*length(dataSort));
val99th = dataSort(idx99);
data(data>val99th) = val99th;
end
