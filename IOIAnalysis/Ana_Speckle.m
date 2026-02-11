function varargout = Ana_Speckle(Folder, bNormalize, varargin)
% ANa_SPECKLE calculates blood flow from laser speckle data.
%
% Usage:
%   Ana_Speckle(Folder, bNormalize)
%   data = Ana_Speckle(Folder, bNormalize)
%   [data, metaData] = Ana_Speckle(Folder, bNormalize, 'filename', bLowRAM)
%
% Inputs:
%   Folder      : Path to the folder containing the dataset.
%   bNormalize  : Logical flag to normalize the output by the temporal mean.
%   filename    : (Optional) Name of the .dat file to process (default: 'speckle').
%   bLowRAM     : (Optional) Logical flag to enable low RAM hybrid processing (default: false).
%
% Outputs:
%   data        : 3D matrix containing blood flow data (Y × X × T-1).
%   metaData    : Structure containing metadata of the processed dataset.
%
% Behavior:
%   - If no output is requested, the function saves "Flow.dat" and "Flow.mat" in the Folder.
%   - Low RAM mode uses a two-pass approach to compute temporal mean and flow.
%   - Standard mode loads the full dataset into memory for faster computation.
%   - The function applies temporal median filtering and converts contrast to flow.


disp('Running Ana speckle...');

%% Parse optional inputs
p = inputParser;
addRequired(p,'Folder', @ischar);
addRequired(p,'bNormalize', @islogical);
addOptional(p,'filename','speckle', @ischar);
addParameter(p,'bRAMsafe', false, @islogical);

parse(p, Folder, bNormalize, varargin{:});
Folder     = p.Results.Folder;
bNormalize = p.Results.bNormalize;
filename   = erase(p.Results.filename, '.dat'); % remove extension if given
bRAMsafe    = p.Results.bRAMsafe;

%% Load metadata
FileList = dir(fullfile(Folder, [filename '.dat']));
if isempty(FileList)
    disp(['No speckle data files found in ' Folder]);
    disp('Speckle Analysis will not run');
    return;
end

Iptr = matfile(fullfile(Folder, [filename '.mat']));
ny = Iptr.datSize(1,1);
nx = Iptr.datSize(1,2);
nt = Iptr.datLength;
tFreq = Iptr.Freq;
speckle_int_time = Iptr.tExposure/1000;

datFile = fullfile(Folder, [filename '.dat']);
OPTIONS.GPU = 0; OPTIONS.Power2Flag = 0; OPTIONS.Brep = 0;
%% Low RAM / standard execution
if bRAMsafe
    % --- Low RAM mode: two-pass processing ---    
    outFile = fullfile(Folder, 'FLOWDATA.dat');
    
    % Determine output size: (T-1) × Y × X   
    outMeta = struct();
    outMeta.datFile = 'FLOWDATA.dat';
    outMeta.datSize = [ny, nx];
    outMeta.datLength = nt-1;
    outMeta.Freq = tFreq;
    outMeta.Datatype = 'single';
    outMeta.dim_names = {'Y','X','T'};
    
    % Preallocate .dat file
    preallocateDatFile(outFile, outMeta);
    
    fidIn  = fopen([Folder filesep filename '.dat'],'r');
    cIn = onCleanup(@() safeFclose(fidIn));
    fidOut = fopen(outFile,'r+');
    cOut = onCleanup(@() safeFclose(fidOut));
    
    % Pass 1: Calculate temporal mean
    fprintf("PASS 1: Calculating temporal mean\n")
    MeanMap = zeros(ny, nx, 'single');
    lastPct = -1;
    for t = 1:nt
        fseek(fidIn, (t-1)*ny*nx*getByteSize('single'),'bof');
        frame = fread(fidIn, ny*nx, '*single');
        frame = reshape(frame, ny, nx);
        MeanMap = MeanMap + frame;
        % Print progress
        pct = floor(100 * t / nt);
        if pct ~= lastPct && mod(pct,10) == 0
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end
    MeanMap = MeanMap / nt;
    fprintf('\nPASS 2: Computing flow...\n')
    % Pass 2: Compute flow and write each frame directly
    lastPct = -1;

    for t = 1:nt-1
        fseek(fidIn, (t-1)*ny*nx*getByteSize('single'),'bof'); % frame t+1
        frameNext = fread(fidIn, ny*nx, '*single');
        frameNext = reshape(frameNext, ny, nx);
        if bNormalize
            % Normalize by temporal mean
            frameNext = frameNext ./ MeanMap;
        end
        % Compute std and smoothed mean for contrast
        speckle_window = fspecial('disk',2) > 0;
        std_laser  = imgaussfilt(stdfilt(frameNext, speckle_window),1);
        mean_laser = imgaussfilt(convnfft(frameNext,speckle_window,'same',1:2,OPTIONS)/sum(speckle_window(:)),1);
        contrast   = std_laser ./ mean_laser;
        flow       = single(private_flow_from_contrast(contrast, speckle_int_time));
        
        % Write frame to output file at correct offset
        fseek(fidOut, (t-1)*ny*nx*getByteSize('single'),'bof');
        fwrite(fidOut, flow, 'single');
        % Print progress
        pct = floor(100 * t / nt);
        if pct ~= lastPct && mod(pct,10) == 0
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end
    % Pass 3: Temporal median filter (chunked over X)
    fW = ceil(0.5 * tFreq);
    
    % Determine chunking along X
    nChunks = calculateMaxChunkSize(ny * nx * (nt-1) * getByteSize('single'), 2);
    chunkX  = ceil(nx / nChunks);
    fprintf('\nPASS 3: Applying Temporal median filter...')
    
    for c = 1:nChunks
        xStart = (c-1)*chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, nx);
        xIdx   = xStart:xEnd;
        
        % Read slab: Y × Xchunk × (T-1)
        slab = spatialSlabIO( ...
            'read', fidOut, ny, nx, nt-1, xIdx, 'single');
        
        % Temporal median filter along T
        slab = medfilt1(slab, fW, [], 3, 'truncate');
        
        % Write filtered slab back to disk
        spatialSlabIO( ...
            'write', fidOut, ny, nx, nt-1, xIdx, 'single',slab);
        % Print progress
        pct = floor(100 * c / nChunks);
        fprintf('%d%% ', pct);        
    end
    
    fprintf("\nFinished")
    fclose(fidIn);
    fclose(fidOut);
    
    % Save meta data
    save(fullfile(Folder, 'Flow.mat'), '-struct', 'outMeta');
    
    % Return filename if requested
    if nargout > 0
        varargout{1} = outFile;
        varargout{2} = outMeta;
    end
    return
% 
else
    % --- Standard: full RAM ---
    fid = fopen(datFile,'r');
    dat = fread(fid, inf, '*single');
    fclose(fid);
    dat = reshape(dat, ny, nx, nt);
    MeanMap = mean(dat,3);

    datOut = zeros(nt-1, ny, nx, 'single');
    speckle_window = fspecial('disk',2)>0;


    for t = 1:nt-1
        tmp_laser = dat(:,:,t);
        if bNormalize
            tmp_laser = tmp_laser ./ MeanMap;
        end
        std_laser = imgaussfilt(stdfilt(tmp_laser,speckle_window),1);
        mean_laser = imgaussfilt(convnfft(tmp_laser,speckle_window,'same',1:2,OPTIONS)/sum(speckle_window(:)),1);
        contrast = std_laser ./ mean_laser;
        datOut(t,:,:) = single(private_flow_from_contrast(contrast,speckle_int_time));
    end

    % Temporal median filter
    fW = ceil(0.5*tFreq);
    datOut = medfilt1(datOut, fW, [], 1, 'truncate');
end

%% Finalize output
datOut = permute(datOut, [2 3 1]);

% Meta data
metaData = struct();
metaData.datName = 'data';
metaData.Stim = Iptr.Stim;
fn = fieldnames(Iptr);
fn = fn(startsWith(fn,'stim_','IgnoreCase',true));
for i = 1:length(fn), metaData.(fn{i}) = Iptr.(fn{i}); end
metaData.datLength = nt-1;
metaData.datSize = Iptr.datSize;
metaData.Freq = Iptr.Freq;
metaData.Datatype = class(datOut);
metaData.dim_names = Iptr.dim_names;
metaData.datFile = fullfile(Folder,'Flow.dat');

if nargout > 0
    varargout{1} = datOut;
    if nargout > 1
        varargout{2} = metaData;
    end
else
    disp('Saving data to file : "Flow.dat"...')
    fFlow = fopen(metaData.datFile,'w');
    fwrite(fFlow, datOut, 'single');
    fclose(fFlow);
    save(fullfile(Folder,'Flow.mat'), '-struct', 'metaData');
end

fprintf('Done!\n');
end

%% Local function
function speed = private_flow_from_contrast(contrast,T)
contrast(isnan(contrast)|contrast<0)=0;
contrast2 = contrast(3:end-2,3:end-2);
mmean = mean(contrast2(:));
sstd  = std(contrast2(:));
tau=(logspace(-15,0,60).^.5);
K  = ((tau/(2*T)).*(1-exp(-2*T*ones(size(tau))./tau))).^(1/2);
[~, index1] = find(K>(mmean-3*sstd),1);
[~, index2] = find(K>(mmean+3*sstd),1);
if isempty(index1), index1=1; end
if isempty(index2)||index2==index1, index2=60; end
Tau2=(logspace(log10(tau(index1)),log10(tau(index2)),40));
K  = ((Tau2/(2*T)).*(1-exp(-2*T*ones(size(Tau2))./Tau2))).^(1/2);
Tau2=[Tau2(1) Tau2 Tau2(end)];
K=[0 K 1e30];
speed=1./interp1(K,Tau2,contrast);
end
