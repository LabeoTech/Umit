function varargout = HemoCorrection(Folder, fData, fMetaData, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% General Infos:
%
% This function is used to remove the hemodynamic fluctuations from any
% fluorescence signal.
%
% Inputs:
% 1. Folder: Folder contaning the dataset to work with.
% 2. fData (3D num. array): Image time series (with dimensions "Y","X","T") 
%       with the fluorescence channel to be corrected.
% 3. fMetaData(struct): meta data structure (from .mat) file associated
%       with the "fluoData".
% 4. Varargin -> if empty: a dialog box will be prompt to ask user which
%                   channels to use to do the correction.
%             -> cell array of string: to specify which channels to use.%                   
% 3. Optional: lowpass filter
%           -> Value of the cutoff frequency to use on a lowpass filter
%           applied to intrinsic signals. THIS PARAMETER IS OPTIONAL
%           To keep the data as is, just use the function with 3-4 parameters
% Ouput:
% - If an output is set, the result of the correction will be given back
% through this output. All the data in the folder will remain unchanged.
% - If no output is specified, the file with the "fData" will be overwritten
% with the corrected data.
%
% Exemples:
%
% 1- HemoCorrection(pwd, fData, fMetaData);
% The fluorescence data in "fData" will be corrected and a dialog
% box will be used in order to select which channels must be used to
% compute the correction.
% 2- NewFluo = HemoCorrection(pwd, fData, fMetaData, {'Green'});
% Only the green channel will be used to compute the correction. The
% corrected data will be stored in the NewFluo variable.
% 3- HemoCorrection(pwd, fData, fMetaData, {'Red, 'Green', 'Amber'});
% All three hemodynamic channels will be used to compute the correction.

if( ~strcmp(Folder(end),filesep) )
    Folder = strcat(Folder, filesep);
end

if( nargin <= 3 )
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
    tmp = varargin{1};
    for ind = 1:size(tmp,2)
        tag = lower(tmp{ind});
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
%--------------------------------------------------------------------------
% NbFrames = fMetaData.datLength;
% HemoData = zeros(size(fn,2), prod(fMetaData.datSize), NbFrames, 'single');
%--------------------------------------------------------------------------

if( nargin <= 4 )
    bFilt = false;
    sFreq = fMetaData.Freq/2;
else
    sFreq = varargin{2};
    bFilt = true;
    if( sFreq >= fMetaData.Freq/2 )
        bFilt = false;
        sFreq = fMetaData.Freq/2;
    end
end


%--------------------------------------------------------------------------
% Get file handles of hemodynamic channels
fid = cell(1,length(fn));
for k = 1:length(fn)
    fid{k} = fopen(fullfile(Folder,fn{k}),'r');
end

% RAM management
Ny = size(fData, 1);
Nx = size(fData, 2);
Nt = size(fData, 3);
numChannels = length(fn);

if bFilt    
    overheadFactor = 8;     
else    
    overheadFactor = 3;
end
nChunks = calculateMaxChunkSize(fData, overheadFactor);
chunkSizePixels = ceil(Nx / nChunks);

% Spatial filter settings
spatSigma = 1;
pad = ceil(3 * spatSigma);   % 3σ Gaussian support

% Design temporal filter 
if bFilt
    f = fdesign.lowpass('N,F3dB', 4, sFreq, fMetaData.Freq);
    lpass = design(f, 'butter');
end

% Normalize fluorescence
fData = reshape(fData, prod(fMetaData.datSize(1:2)), []);
m_fData = mean(fData, 2);
fData = (fData - m_fData) ./ m_fData;

%% ========================================================================
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
    % ----- Build pixel index list (matches reshape order)
    [COL, ROW] = meshgrid(idxPixels, 1:Ny);
    indList = sub2ind([Ny, Nx], ROW(:), COL(:));
    Np = numel(indList);
    % ----- Preallocate hemodynamic block
    HemoData = zeros(numChannels, Np, Nt, 'single');

    %% ---------------------------------------------------------------
    %  Load, filter, normalize hemodynamic channels
    % ---------------------------------------------------------------
    for kk = 1:numChannels
        waitbar(.99, h, ['Reading file [' fn{kk} ']']);drawnow()
        % Read padded spatial slab
        tmp = readSpatialSlab(fid{kk}, Ny, Nx, Nt, idxPixels_with_pad, 'single');
        tmp_sz = size(tmp);

        % Reshape to [pixels × time] for temporal filtering
        tmp = reshape(tmp, [], tmp_sz(3));

        % Temporal filtering (optional)
        if bFilt
            waitbar(.99, h, ['Applying temporal filter [' fn{kk} ']']);drawnow()
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
        end

        % Restore spatial shape for Gaussian blur
        tmp = reshape(tmp, tmp_sz);
        tmp = imgaussfilt(tmp, spatSigma, 'Padding', 'symmetric');

        % Crop padding
        tmp = tmp(:, padStart+1:end-padStop, :);

        % Back to [pixels × time]
        tmp = reshape(tmp, [], tmp_sz(3));

        % Normalize each pixel
        m = mean(tmp, 2);
        tmp = (tmp - m) ./ m;

        % Store
        HemoData(kk, :, :) = tmp;
    end
    
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
for kk = 1:length(fid)
    fclose(fid{kk});
end
% Undo normalization
fData = fData .* m_fData + m_fData;
fData = reshape(fData, fMetaData.datSize(1), fMetaData.datSize(2), []);



if( nargout == 0 )
    [~,filename,ext] = fileparts(fMetaData.datFile);
    eval(['fid = fopen(''' Folder filename ext ''', ''w'');']);
    fwrite(fid, fData, 'single');
    fclose(fid);
else
    varargout{:} = fData;
end

fprintf('Finished Hemodynamic Correction.\n');
end