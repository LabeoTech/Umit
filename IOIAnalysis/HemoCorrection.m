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
% -------------------------------------------------------------------------
% Inputs:
%
%   Folder :
%       Path to the folder containing the dataset.
%
%   FileData :
%       Either:
%         - 3D numeric array [Y, X, T] containing fluorescence data
%           ? STANDARD MODE
%         - String or char array pointing to a .dat fluorescence file
%           ? LOW-RAM MODE
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
%   HemoCorrection(pwd, 'fluo_475.dat', 'fluo_475.mat', {'Red','Green'});
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
            fn{end+1} = cList(ind).name;
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
                    fn{end+1} = 'rChan.dat';
                else
                    fn{end+1} = 'red.dat';
                end
            case {'amber', 'yellow'}
                if( exist([Folder 'yChan.dat'], 'file') )
                    fn{end+1} = 'yChan.dat';
                else
                    fn{end+1} = 'yellow.dat';
                end
            case 'green'
                if( exist([Folder 'gChan.dat'], 'file') )
                    fn{end+1} = 'gChan.dat';
                else
                    fn{end+1} = 'green.dat';
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
        eval(['fid = fopen(''' Folder filename ext ''', ''w'');']);
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
% and hemodynamic channels can be loaded into RAM.
%
% Features:
%   - Automatic chunking along X dimension if RAM is limited
%   - Optional temporal low-pass filtering of hemodynamic channels
%   - Spatial Gaussian filtering with symmetric padding
%
% Inputs:
%   fData :
%       Fluorescence data array [Y, X, T]
%
%   fMetaData :
%       Metadata structure associated with fData
%
%   colorList :
%       Cell array of filenames for hemodynamic channels (.dat)
%
%   LPcutoffFreq :
%       Low-pass cutoff frequency in Hz (0 disables filtering)
%
% Output:
%   fData :
%       Hemodynamically corrected fluorescence data [Y, X, T]
%
% Notes:
%   - This mode is faster than low-RAM mode but requires sufficient memory.
%   - Data normalization and denormalization are performed internally.

%--------------------------------------------------------------------------
Ny = fMetaData.datSize(1);
Nx = fMetaData.datSize(2);
Nt = fMetaData.datLength;
Np = Nx *Ny;
%--------------------------------------------------------------------------
% Normalize fluorescence
%--------------------------------------------------------------------------
fData = reshape(fData, prod(fMetaData.datSize(1:2)), []);
m_fData = mean(fData, 2);
fData = (fData- m_fData) ./ m_fData;
%--------------------------------------------------------------------------

% Design temporal filter
if LPcutoffFreq
    f = fdesign.lowpass('N,F3dB', 4, LPcutoffFreq, fMetaData.Freq);
    lpass = design(f, 'butter');
end

% Calculate available RAM
numChannels = numel(colorList);
nChunks = calculateMaxChunkSize(fData, numChannels,0.15);
chunkSizePixels = ceil(Nx / nChunks);
if nChunks > 1
    fid = {};
    for ii = 1:numChannels
        fid{ii} = fopen(colorList{ii},'r');%#ok
    end
end

% Spatial filter settings
spatSigma = 1;

h = waitbar(0, 'Fitting Hemodynamics...');
for ii = 1:nChunks
    
    if nChunks == 1
        h.Name = ['Hemodynamic Correction'];drawnow()
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
            % slower
            tmp = spatialSlabIO('read',fid{kk}, Ny, Nx, Nt, idxPixels_with_pad, 'single');
        end
        
        tmp_sz = size(tmp);
        
        % Reshape for temporal filtering
        tmp = reshape(tmp, [], tmp_sz(3));
        
        % Temporal filtering (optional)
        if LPcutoffFreq
            waitbar(.99, h, ['Applying temporal filter [' ,colorName, ext ']']);drawnow()
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
        end
        
        % Restore spatial shape for Gaussian blur
        tmp = reshape(tmp, tmp_sz);
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
    end
    clear tmp
    %% ---------------------------------------------------------------
    %  Hemodynamic regression
    % ---------------------------------------------------------------
    warning('off', 'MATLAB:rankDeficientMatrix');
    waitbar(0, h, 'Performing Hemodynamic correction...');drawnow()
    for indF = 1:Np
        X = [ ones(1, Nt); linspace(0, 1, Nt);squeeze(HemoData(:, indF, :)) ];
        B = X' \ fData(indList(indF), :)';
        fData(indList(indF), :) = fData(indList(indF), :) - (X' * B)';
        % Update waitbar
        if mod(indF, 500) == 0
            waitbar(indF / Np, h);
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
% strategy. All data are read directly from binary files in spatial chunks,
% processed, and written back to a preallocated output file.
%
% This mode is intended for very large datasets that do not fit in memory.
%
% Inputs:
%   outFilename :
%       Output .dat file where corrected fluorescence will be written
%
%   fluoFile :
%       Input fluorescence .dat file
%
%   fMetaData :
%       Metadata .mat file associated with the fluorescence data
%
%   colorList :
%       Cell array of hemodynamic channel .dat filenames
%
%   LPcutoffFreq :
%       Low-pass cutoff frequency in Hz (0 disables filtering)
%
% Output:
%   outFilename :
%       Filename of the corrected fluorescence data (.dat)
%
% Implementation details:
%   - Spatial chunking along X dimension
%   - Temporal processing performed per chunk
%   - Uses spatialSlabIO for safe random-access disk I/O
%   - Output file is preallocated before processing
%
% Notes:
%   - This mode minimizes peak RAM usage at the cost of execution time.
%   - Suitable for long recordings and high-resolution imaging data.


% Get Fluo file handle
f_fid = fopen(fluoFile,'r');
c_f = onCleanup(@() safeFclose(f_fid));
% Get file handles of hemodynamic channels
numChannels = length(colorList);
h_fid = cell(1,numChannels);
c_r = cell(1,length(colorList));
for k = 1:length(colorList)
    h_fid{k} = fopen(colorList{k},'r');
    c_r{k} = onCleanup(@() safeFclose(h_fid{k}));
end

% RAM management
Ny = fMetaData.datSize(1);
Nx = fMetaData.datSize(2);
Nt = fMetaData.datLength;

% Estimate data volume in bytes from fluo file
dataBytes = prod([fMetaData.datSize, fMetaData.datLength, getByteSize(fMetaData.Datatype)]);
% Calculate number of chunks for data processing
nChunks = calculateMaxChunkSize(dataBytes,2+numel(colorList),.1);
chunkSizePixels = ceil(Nx / nChunks);

% Spatial filter settings
spatSigma = 1;
pad = ceil(3 * spatSigma);   % 3σ Gaussian support

% Design temporal filter
if LPcutoffFreq
    f = fdesign.lowpass('N,F3dB', 4, LPcutoffFreq, fMetaData.Freq);
    lpass = design(f, 'butter');
end
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
    Np = numel(idxPixels)*Ny;
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
        
        % Read padded spatial slab
        tmp = spatialSlabIO('read',h_fid{kk}, Ny, Nx, Nt, idxPixels_with_pad, fMetaData.Datatype);
        tmp_sz = size(tmp);
                        
        % Temporal filtering (optional)
        if LPcutoffFreq
            
            waitbar(.99, h, ['Applying temporal filter [', colorName, ext,  ']']);drawnow()
            % Reshape to pixel by time for temporal filtering
            tmp = reshape(tmp, [], tmp_sz(3));
            % Apply temporal filter
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
            % Restore spatial shape for Gaussian blur
            tmp = reshape(tmp, tmp_sz);
        end
        
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
    for indF = 1:Np
        
        X = [ ones(1, Nt); linspace(0, 1, Nt);squeeze(HemoData(:, indF, :)) ];
        B = X' \ fData(indF, :)';
        fData(indF, :) = fData(indF, :) - (X' * B)';
        % Update waitbar
        if mod(indF, 500) == 0
            waitbar(indF / Np, h);
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