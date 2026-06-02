function varargout = HemoCorrection(Folder, FileData, fMetaData, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEMOCORRECTION  Remove hemodynamic fluctuations from fluorescence signals.
%
% This function performs pixel-wise hemodynamic regression on fluorescence
% image time series using one or more intrinsic/hemodynamic channels.
%
% The function supports TWO execution modes:
%
%   1) STANDARD MODE (in-memory)
%      - Triggered when FileData is a numeric array
%      - All data are loaded into RAM
%      - Faster execution but higher memory usage
%
%   2) LOW-RAM MODE (streaming / hybrid)
%      - Triggered when FileData is a filename (.dat)
%      - Data are streamed directly from binary files in spatial chunks
%      - Intermediate results are written back to disk
%      - Minimal RAM footprint, suitable for very large datasets
%
% Hemodynamic channels are temporally resampled to match the fluorescence
% channel before regression. The fluorescence channel is always the reference
% timeline. If a hemodynamic channel has a higher sampling rate than the
% fluorescence channel, an anti-aliasing low-pass filter is applied before
% temporal downsampling.
%
% -------------------------------------------------------------------------
% Inputs:
%
%   Folder :
%       Path to the folder containing the dataset.
%
%   FileData :
%       Either:
%         - 3D numeric array [Y, X, T] containing fluorescence data
%           STANDARD MODE
%         - String or char array pointing to a .dat fluorescence file
%           LOW-RAM MODE
%
%   fMetaData :
%       Metadata associated with the fluorescence data. Can be:
%         - Struct loaded from .mat
%         - matlab.io.MatFile object
%
%   varargin :
%       Optional inputs:
%
%       cList :
%           Cell array of strings specifying which hemodynamic channels
%           to use (e.g. {'red','green','amber'}).
%           If empty or omitted, a dialog box prompts the user.
%
%       sFreq :
%           Low-pass cutoff frequency (Hz) applied to hemodynamic channels.
%           If 0 or omitted, no temporal filtering is applied.
%
% -------------------------------------------------------------------------
% Outputs:
%
%   If an output argument is requested:
%       - Returns the corrected fluorescence data (array or filename,
%         depending on execution mode).
%
%   If no output argument is requested:
%       - The original fluorescence data file is overwritten in-place.
%
% -------------------------------------------------------------------------
% Examples:
%
%   % Interactive channel selection (standard mode)
%   fCorr = HemoCorrection(pwd, fData, fMetaData);
%
%   % Explicit channel list (standard mode)
%   fCorr = HemoCorrection(pwd, fData, fMetaData, {'Green'});
%
%   % Low-RAM streaming mode (file-based)
%   HemoCorrection(pwd, 'fluo_475.dat', fMetaData, {'Red','Green'});
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'Folder', @isfolder);
addRequired(p,'FileData',@(x) (isnumeric(x) & ndims(x) == 3) | ischar(x));
addRequired(p,'fMetaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x));
addOptional(p, 'cList', {}, @(x) isempty(x) || (iscell(x) && all(cellfun(@ischar, x))));
addOptional(p, 'sFreq', 0,@(x) isnumeric(x) & isscalar(x));
addParameter(p,'outFilename','HEMOCORRECTED_DATA.dat',@ischar)
% Parse inputs:
parse(p,Folder, FileData, fMetaData, varargin{:});
cList = p.Results.cList;
sFreq = p.Results.sFreq;
outFilename = p.Results.outFilename;

if( ~strcmp(Folder(end),filesep) )
    Folder = strcat(Folder, filesep);
end

if isempty(cList)
    cList = dir([Folder '*.dat']);
    fn = {};
    for ind = 1:size(cList,1)
        if( ~strcmp(cList(ind).name(1),'f') )
            fn{end+1} = cList(ind).name; %#ok<AGROW>
        end
    end
    
    [idx, tf] = listdlg('PromptString',{'Select channels to be used to',...
        'compute hemodynamic correction.',''},...
        'ListString',fn);
    
    if( tf == 0 )
        return;
    end
    
    fn = fn(idx);
    clear cList idx ind tf;
else
    fn = {};
    
    for ind = 1:size(cList,2)
        tag = lower(cList{ind});
        switch tag
            case 'red'
                if( exist([Folder 'rChan.dat'], 'file') )
                    fn{end+1} = 'rChan.dat'; %#ok<AGROW>
                else
                    fn{end+1} = 'red.dat'; %#ok<AGROW>
                end
            case {'amber', 'yellow'}
                if( exist([Folder 'yChan.dat'], 'file') )
                    fn{end+1} = 'yChan.dat'; %#ok<AGROW>
                else
                    fn{end+1} = 'yellow.dat'; %#ok<AGROW>
                end
            case 'green'
                if( exist([Folder 'gChan.dat'], 'file') )
                    fn{end+1} = 'gChan.dat'; %#ok<AGROW>
                else
                    fn{end+1} = 'green.dat'; %#ok<AGROW>
                end
        end
    end
end

if sFreq
    freq = fMetaData.Freq;
    if ( sFreq >= freq/2 )
        sFreq = 0; % No Temporal filter applied
    end
end

fn = fullfile(Folder,fn);
if ischar(FileData)
    % Execute in Low RAM mode
    
    % Set all file names as full path
    FileData = fullfile(Folder,FileData);
    outFilename = fullfile(Folder,outFilename);
    
    % Execute LowRAM mode
    outFileData = HemoCorrection_lowRAMmode(outFilename, FileData, fMetaData, fn,sFreq);
    % Output file name
    [~,outFileData,ext] = fileparts(outFileData);
    varargout{1} = [outFileData,ext];
else
    % Execute in Standard mode
    FileData = HemoCorrection_standardMode(FileData,fMetaData, fn,sFreq);
    if nargout
        % Output fluo file
        varargout{1} = FileData;
    else
        % Overwrite fluo file
        [~,filename,ext] = fileparts(fMetaData.datFile);
        fid = fopen([Folder filename ext], 'w');
        fwrite(fid, FileData, 'single');
        fclose(fid);
    end
end

end

%% ========================================================================
% Local functions
% =========================================================================
function fData = HemoCorrection_standardMode(fData,fMetaData, colorList,LPcutoffFreq)
% HEMOCORRECTION_STANDARDMODE  In-memory hemodynamic correction.
%
% This function performs hemodynamic correction assuming all fluorescence
% and hemodynamic channels can be loaded into RAM. Hemodynamic channels are
% temporally resampled to match the fluorescence channel when needed.

%--------------------------------------------------------------------------
Ny = fMetaData.datSize(1,1);
Nx = fMetaData.datSize(1,2);
Nt = fMetaData.datLength;
Np = Nx * Ny;
%--------------------------------------------------------------------------
% Resolve hemodynamic metadata before RAM estimation.
%--------------------------------------------------------------------------
numChannels = numel(colorList);
cMetaData = cell(1,numChannels);
maxNt = Nt;
for kk = 1:numChannels
    [cFolder,cName] = fileparts(colorList{kk});
    cMatFile = fullfile(cFolder, [cName '.mat']);
    if( ~exist(cMatFile, 'file') )
        error('Missing metadata file for hemodynamic channel: %s', cMatFile);
    end
    cMetaData{kk} = load(cMatFile);
    if( any(cMetaData{kk}.datSize ~= fMetaData.datSize) )
        error('Hemodynamic channel %s has incompatible spatial dimensions.', colorList{kk});
    end
    maxNt = max(maxNt, cMetaData{kk}.datLength);
end
%--------------------------------------------------------------------------
% Validate that all selected channels span the same recording duration.
% This prevents interpolation from hiding cropped or truncated channels.
%--------------------------------------------------------------------------
durationToleranceSec = 1e-3;
channelLabelList = [{'fluorescence'}, colorList(:)'];
lengthList = zeros(1, numChannels + 1);
freqList = zeros(1, numChannels + 1);

lengthList(1) = size(fData, 3);
freqList(1) = double(fMetaData.Freq);

assert(lengthList(1) == double(fMetaData.datLength), ...
    ['Fluorescence data length (%d frames) does not match fluorescence ' ...
     'metadata (%d frames).'], ...
    lengthList(1), double(fMetaData.datLength));

for kk = 1:numChannels
    fileInfo = dir(colorList{kk});
    assert(~isempty(fileInfo), ...
        'Hemodynamic channel file was not found: %s', colorList{kk});

    bytesPerFrame = prod(double(cMetaData{kk}.datSize)) * ...
        getByteSize(cMetaData{kk}.Datatype);
    actualNt = fileInfo.bytes / bytesPerFrame;

    assert(isfinite(actualNt) && actualNt > 0 && abs(actualNt - round(actualNt)) < 1e-9, ...
        ['Hemodynamic channel %s has a file size that is incompatible with ' ...
         'its spatial metadata and datatype.'], ...
        colorList{kk});

    actualNt = round(actualNt);

    assert(actualNt == double(cMetaData{kk}.datLength), ...
        ['Hemodynamic channel %s contains %d frames on disk, but its ' ...
         'metadata declares %d frames.'], ...
        colorList{kk}, actualNt, double(cMetaData{kk}.datLength));

    lengthList(kk + 1) = actualNt;
    freqList(kk + 1) = double(cMetaData{kk}.Freq);
end

durationListSec = lengthList ./ freqList;
durationDeltaSec = abs(durationListSec - durationListSec(1));

if any(durationDeltaSec > durationToleranceSec)
    durationMsg = '';
    for kk = 1:numel(channelLabelList)
        durationMsg = sprintf('%s%s: %d frames at %.10g Hz = %.10g s\n', ...
            durationMsg, channelLabelList{kk}, lengthList(kk), freqList(kk), durationListSec(kk));
    end

    assert(false, ...
        ['Cannot resample hemodynamic channels because the selected ' ...
         'channels do not span the same recording duration within %.4g s. ' ...
         'This usually means that one file was cropped, truncated, or no ' ...
         'longer matches its metadata.\n%s'], ...
        durationToleranceSec, durationMsg);
end

%--------------------------------------------------------------------------
% Normalize fluorescence
%--------------------------------------------------------------------------
fData = reshape(fData, prod(fMetaData.datSize), []);
m_fData = mean(fData, 2);
fData = (fData - m_fData) ./ m_fData;
%--------------------------------------------------------------------------

% Calculate available RAM. Use the largest temporal length because repeated
% illuminations can create hemodynamic channels longer than fluorescence.
dataBytes = max(numel(fData) * getByteSize(class(fData)), ...
    prod(fMetaData.datSize) * maxNt * getByteSize('single'));
nChunks = calculateMaxChunkSize(dataBytes, 2 + numChannels, 0.15);
chunkSizePixels = ceil(Nx / nChunks);
if nChunks > 1
    fid = {};
    for ii = 1:numChannels
        fid{ii} = fopen(colorList{ii},'r'); %#ok<AGROW>
    end
end

% Spatial filter settings
spatSigma = 1;
pad = ceil(3 * spatSigma);
h = waitbar(0, 'Fitting Hemodynamics...');
h_out = onCleanup(@() delete(h));
for ii = 1:nChunks
    
    if nChunks == 1
        h.Name = 'Hemodynamic Correction';drawnow()
        % There is enough RAM to load all the channels (No chunking needed)
        HemoData = zeros(numChannels, Np, Nt, 'single');
        padStart = 0;
        padStop = 0;
        indList = 1:Np;
    else
        % Some chunking needed
        h.Name = ['Hemodynamic Corr. (chunk ' num2str(ii) '/' num2str(nChunks) ')'];drawnow()
        % ----- X-range for this chunk
        pxStart = (ii - 1) * chunkSizePixels + 1;
        pxEnd   = min(ii * chunkSizePixels, Nx);
        idxPixels = pxStart:pxEnd;
        % ----- Padding
        padStart = min(pad, pxStart - 1);
        padStop  = min(pad, Nx - pxEnd);
        idxPixels_with_pad = (pxStart - padStart):(pxEnd + padStop);
        % ----- Build pixel index list (matches reshape order)
        [COL, ROW] = meshgrid(idxPixels, 1:Ny);
        indList = sub2ind([Ny, Nx], ROW(:), COL(:));
        % ----- Preallocate hemodynamic block
        HemoData = zeros(numChannels, numel(indList), Nt, 'single');
    end
    for kk = 1:numChannels
        [~,colorName,ext] = fileparts(colorList{kk});
        if nChunks == 1
            tmp = loadDatFile(colorList{kk});
        else
            tmp = spatialSlabIO('read',fid{kk}, Ny, Nx, cMetaData{kk}.datLength, ...
                idxPixels_with_pad, cMetaData{kk}.Datatype);
        end
        
        tmp = iResampleHemoToFluoTimeline(tmp, cMetaData{kk}, fMetaData, LPcutoffFreq);
        tmp_sz = size(tmp);
        
        tmp = imgaussfilt(tmp, spatSigma, 'Padding', 'symmetric');
        
        % Crop padding
        tmp = tmp(:, padStart+1:end-padStop, :);
        
        % Back to [pixels vs time]
        tmp = reshape(tmp, [], tmp_sz(3));
        
        % Normalize each pixel
        m = mean(tmp, 2);
        tmp = (tmp - m) ./ m;
        
        % Store
        HemoData(kk, :, :) = tmp;
        clear tmp m tmp_sz
        waitbar(.99, h, ['Loaded hemodynamic channel [' ,colorName, ext ']']);drawnow()
    end
    %% ---------------------------------------------------------------
    %  Hemodynamic regression
    % ---------------------------------------------------------------
    warning('off', 'MATLAB:rankDeficientMatrix');
    waitbar(0, h, 'Performing Hemodynamic correction...');drawnow()
    nPixelsThisChunk = numel(indList);
    for indP = 1:nPixelsThisChunk
        if size(HemoData,1) == 1
            X = [ ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :))' ];
        else
            X = [ ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :)) ];
        end
        B = X' \ fData(indList(indP), :)';
        fData(indList(indP), :) = fData(indList(indP), :) - (X' * B)';
        % Update waitbar
        if mod(indP, 500) == 0
            waitbar(indP / nPixelsThisChunk, h);
        end
    end
    
    warning('on', 'MATLAB:rankDeficientMatrix');
end

close(h);

% Close Hemo Data files
if exist('fid','var')
    for kk = 1:length(fid)
        fclose(fid{kk});
    end
end
% Undo normalization
fData = fData .* m_fData + m_fData;
fData = reshape(fData, fMetaData.datSize(1), fMetaData.datSize(2), []);

end


function outFilename = HemoCorrection_lowRAMmode(outFilename, fluoFile, fMetaData, colorList,LPcutoffFreq)
% HEMOCORRECTION_LOWRAMMODE  Disk-streamed hemodynamic correction.
%
% This function performs hemodynamic correction using a low-RAM streaming
% strategy. Hemodynamic channels are read with their own metadata and
% temporally resampled to the fluorescence timebase inside each spatial slab.

% Get Fluo file handle
f_fid = fopen(fluoFile,'r');
c_f = onCleanup(@() safeFclose(f_fid));
% Get file handles and metadata of hemodynamic channels
numChannels = length(colorList);
h_fid = cell(1,numChannels);
c_r = cell(1,length(colorList));
cMetaData = cell(1,numChannels);
maxNt = fMetaData.datLength;
for k = 1:length(colorList)
    h_fid{k} = fopen(colorList{k},'r');
    c_r{k} = onCleanup(@() safeFclose(h_fid{k}));
    [cFolder,cName] = fileparts(colorList{k});
    cMatFile = fullfile(cFolder, [cName '.mat']);
    if( ~exist(cMatFile, 'file') )
        error('Missing metadata file for hemodynamic channel: %s', cMatFile);
    end
    cMetaData{k} = load(cMatFile);
    if( any(cMetaData{k}.datSize ~= fMetaData.datSize) )
        error('Hemodynamic channel %s has incompatible spatial dimensions.', colorList{k});
    end
    maxNt = max(maxNt, cMetaData{k}.datLength);
end

% Validate that all selected channels span the same recording duration.
% This prevents interpolation from hiding cropped or truncated channels.
durationToleranceSec = 1e-3;
channelLabelList = [{'fluorescence'}, colorList(:)'];
lengthList = zeros(1, numChannels + 1);
freqList = zeros(1, numChannels + 1);

fluoInfo = dir(fluoFile);
assert(~isempty(fluoInfo), ...
    'Fluorescence file was not found: %s', fluoFile);

fluoBytesPerFrame = prod(double(fMetaData.datSize)) * getByteSize(fMetaData.Datatype);
actualFluoNt = fluoInfo.bytes / fluoBytesPerFrame;

assert(isfinite(actualFluoNt) && actualFluoNt > 0 && abs(actualFluoNt - round(actualFluoNt)) < 1e-9, ...
    ['Fluorescence file size is incompatible with its spatial metadata ' ...
     'and datatype.']);

actualFluoNt = round(actualFluoNt);

assert(actualFluoNt == double(fMetaData.datLength), ...
    ['Fluorescence file contains %d frames on disk, but its metadata ' ...
     'declares %d frames.'], ...
    actualFluoNt, double(fMetaData.datLength));

lengthList(1) = actualFluoNt;
freqList(1) = double(fMetaData.Freq);

for kk = 1:numChannels
    fileInfo = dir(colorList{kk});
    assert(~isempty(fileInfo), ...
        'Hemodynamic channel file was not found: %s', colorList{kk});

    bytesPerFrame = prod(double(cMetaData{kk}.datSize)) * ...
        getByteSize(cMetaData{kk}.Datatype);
    actualNt = fileInfo.bytes / bytesPerFrame;

    assert(isfinite(actualNt) && actualNt > 0 && abs(actualNt - round(actualNt)) < 1e-9, ...
        ['Hemodynamic channel %s has a file size that is incompatible with ' ...
         'its spatial metadata and datatype.'], ...
        colorList{kk});

    actualNt = round(actualNt);

    assert(actualNt == double(cMetaData{kk}.datLength), ...
        ['Hemodynamic channel %s contains %d frames on disk, but its ' ...
         'metadata declares %d frames.'], ...
        colorList{kk}, actualNt, double(cMetaData{kk}.datLength));

    lengthList(kk + 1) = actualNt;
    freqList(kk + 1) = double(cMetaData{kk}.Freq);
end

durationListSec = lengthList ./ freqList;
durationDeltaSec = abs(durationListSec - durationListSec(1));

if any(durationDeltaSec > durationToleranceSec)
    durationMsg = '';
    for kk = 1:numel(channelLabelList)
        durationMsg = sprintf('%s%s: %d frames at %.10g Hz = %.10g s\n', ...
            durationMsg, channelLabelList{kk}, lengthList(kk), freqList(kk), durationListSec(kk));
    end

    assert(false, ...
        ['Cannot resample hemodynamic channels because the selected ' ...
         'channels do not span the same recording duration within %.4g s. ' ...
         'This usually means that one file was cropped, truncated, or no ' ...
         'longer matches its metadata.\n%s'], ...
        durationToleranceSec, durationMsg);
end

% RAM management
Ny = fMetaData.datSize(1,1);
Nx = fMetaData.datSize(1,2);
Nt = fMetaData.datLength;

% Estimate data volume in bytes from the largest temporal input.
dataBytes = prod([fMetaData.datSize, maxNt, getByteSize(fMetaData.Datatype)]);
% Calculate number of chunks for data processing
nChunks = calculateMaxChunkSize(dataBytes,2+numel(colorList),.1);
chunkSizePixels = ceil(Nx / nChunks);

% Spatial filter settings
spatSigma = 1;
pad = ceil(3 * spatSigma);   % 3σ Gaussian support

% ---------------------------------------------------------------------
% Preallocate output binary file and set to write
% ---------------------------------------------------------------------
preallocateDatFile(outFilename,fMetaData);
fid_out = fopen(outFilename, 'r+');
c_out = onCleanup(@() safeFclose(fid_out));

% ========================================================================
%                          Main chunk loop
% ========================================================================
h = waitbar(0, 'Fitting Hemodynamics...');
h_out = onCleanup(@() delete(h));
for ii = 1:nChunks
    h.Name = ['Hemodynamic Corr. (chunk ' num2str(ii) '/' num2str(nChunks) ')'];drawnow()
    % ----- X-range for this chunk
    pxStart = (ii - 1) * chunkSizePixels + 1;
    pxEnd   = min(ii * chunkSizePixels, Nx);
    idxPixels = pxStart:pxEnd;
    
    % ----- Padding
    padStart = min(pad, pxStart - 1);
    padStop  = min(pad, Nx - pxEnd);
    idxPixels_with_pad = (pxStart - padStart):(pxEnd + padStop);
            
    % ----- Preallocate hemodynamic block
    Np = numel(idxPixels) * Ny;
    HemoData = zeros(numChannels, Np, Nt, 'single');
    
    % ---------------------------------------------------------------------
    % Normalize chunk of fluo channel
    % ---------------------------------------------------------------------
    waitbar(.99, h, 'Reading fluo channel...');drawnow()
    fData = spatialSlabIO('read', f_fid, Ny, Nx, Nt, idxPixels, fMetaData.Datatype);
    waitbar(.99, h, 'Normalizing fluo channel...');drawnow()
    f_slabSz = size(fData);
    fData = reshape(fData,[],Nt);
    m_fData = mean(fData,2);
    fData = (fData - m_fData)./m_fData;
    
    %% ---------------------------------------------------------------
    %  Load, filter, normalize hemodynamic channels
    % ---------------------------------------------------------------    
    for kk = 1:numChannels
        [~,colorName,ext] = fileparts(colorList{kk});
        waitbar(.99, h, ['Reading file [', colorName, ext, ']']);drawnow()
        
        % Read padded spatial slab using this channel's own length/type.
        tmp = spatialSlabIO('read',h_fid{kk}, Ny, Nx, cMetaData{kk}.datLength, ...
            idxPixels_with_pad, cMetaData{kk}.Datatype);
        
        waitbar(.99, h, ['Resampling hemodynamic file [', colorName, ext, ']']);drawnow()
        tmp = iResampleHemoToFluoTimeline(tmp, cMetaData{kk}, fMetaData, LPcutoffFreq);
        tmp_sz = size(tmp);
        
        waitbar(.99, h, 'Applying spatial filter to hemodynamic data...');drawnow()
        % Apply spatial filter        
        tmp = imgaussfilt(tmp, spatSigma, 'Padding', 'symmetric');
        
        % Crop padding
        tmp = tmp(:, padStart+1:end-padStop, :);
        
        % Back to [pixels by time]
        tmp = reshape(tmp, [], tmp_sz(3));
        waitbar(.99, h, 'Normalizing hemodynamic data...');drawnow()
        
        % Normalize each pixel
        m = mean(tmp, 2);
        tmp = (tmp - m) ./ m;
        
        % Store
        HemoData(kk, :, :) = tmp;
                
        clear tmp m tmp_sz
    end
    
    %% ---------------------------------------------------------------
    %  Hemodynamic regression
    % ---------------------------------------------------------------    
    waitbar(0, h, 'Performing Hemodynamic correction...');drawnow()
    warning('off', 'MATLAB:rankDeficientMatrix');
    for indP = 1:Np
        if size(HemoData,1) == 1
            X = [ ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :))' ];
        else
            X = [ ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :)) ];
        end
        
        B = X' \ fData(indP, :)';
        fData(indP, :) = fData(indP, :) - (X' * B)';
        % Update waitbar
        if mod(indP, 500) == 0
            waitbar(indP / Np, h);
        end
    end
    clear B X HemoData
    warning('on', 'MATLAB:rankDeficientMatrix');
    
    % ---------------------------------------------------------------
    %  Undo normalization in fluo channel and reshape for saving
    % ---------------------------------------------------------------
    fData = fData .* m_fData + m_fData;
    fData = reshape(fData,f_slabSz);
   
    % ---------------------------------------------------------------
    %  Write to .dat file
    % ---------------------------------------------------------------
    waitbar(0.99, h, 'Writing corrected fluo to file...');drawnow()
    spatialSlabIO('write', fid_out, Ny, Nx, Nt, idxPixels, fMetaData.Datatype, fData);   
    clear fData 
    
end
close(h);
% ---------------------------------------------------------------
%  Close .dat files
% ---------------------------------------------------------------

% Close fluo data file
fclose(f_fid);
% Close output data file
fclose(fid_out);
% Close Hemo Data files
for kk = 1:length(h_fid)
    fclose(h_fid{kk});
end

end


function tmp = iResampleHemoToFluoTimeline(tmp, cMetaData, fMetaData, LPcutoffFreq)
%IRESAMPLEHEMOTOFLUOTIMELINE Match one hemodynamic slab to fluorescence T.
%
%   tmp = iResampleHemoToFluoTimeline(tmp, cMetaData, fMetaData, LPcutoffFreq)
%   applies anti-aliasing/user low-pass filtering when needed and resamples
%   tmp along the third dimension to fMetaData.datLength. The fluorescence
%   channel is the reference timeline.

NtHemo = cMetaData.datLength;
NtFluo = fMetaData.datLength;
freqHemo = cMetaData.Freq;
freqFluo = fMetaData.Freq;

if( size(tmp,3) ~= NtHemo )
    error('Input hemodynamic slab length does not match its metadata.');
end

% Apply user-requested low-pass when valid. If no user filter is requested
% and the hemodynamic channel is being downsampled, apply anti-aliasing.
cutoffFreq = 0;
if( LPcutoffFreq > 0 )
    cutoffFreq = LPcutoffFreq;
elseif( freqHemo > freqFluo && NtHemo > NtFluo )
    cutoffFreq = 0.45 * freqFluo;
end

if( cutoffFreq > 0 && cutoffFreq < freqHemo/2 )
    sz = size(tmp);
    tmp = reshape(tmp, [], sz(3));
    f = fdesign.lowpass('N,F3dB', 4, cutoffFreq, freqHemo);
    lpass = design(f, 'butter');
    tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
    tmp = reshape(tmp, sz);
end

if( NtHemo ~= NtFluo )
    sz = size(tmp);
    xHemo = linspace(0, 1, NtHemo);
    xFluo = linspace(0, 1, NtFluo);
    tmp = reshape(tmp, [], NtHemo);
    tmp = interp1(xHemo, single(tmp)', xFluo, 'linear', 'extrap')';
    tmp = reshape(single(tmp), sz(1), sz(2), NtFluo);
else
    tmp = single(tmp);
end

end
