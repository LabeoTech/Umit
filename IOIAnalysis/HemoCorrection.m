function varargout = HemoCorrection(Folder, FileData, fMetaData, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEMOCORRECTION Remove hemodynamic fluctuations from fluorescence signals.
%
% This function performs pixel-wise hemodynamic regression on fluorescence
% image time series using one or more intrinsic/hemodynamic channels.
%
% The function supports two execution modes:
%
%   1) STANDARD MODE (in-memory)
%      - Triggered when FileData is a numeric array
%      - All data are loaded into RAM
%
%   2) LOW-RAM MODE (streaming / hybrid)
%      - Triggered when FileData is a filename (.dat)
%      - Data are streamed in spatial chunks and written to disk
%
% Inputs
% ------
% Folder :
%   Path to the folder containing the dataset.
%
% FileData :
%   Either:
%       - 3D numeric array [Y, X, T]
%       - .dat fluorescence filename
%
% fMetaData :
%   Metadata associated with the fluorescence signal.
%
% varargin :
%   Optional inputs:
%       cList :
%           Cell array specifying hemodynamic channels to use. Supports:
%               {'red','green','amber'}
%           and custom filenames such as:
%               {'myReference.dat'}
%
%       sFreq :
%           Low-pass cutoff frequency in Hz. Set to 0 to disable filtering.
%
% Output
% ------
% If an output is requested:
%   - corrected fluorescence array or output filename
%
% If no output is requested:
%   - standard mode overwrites original fluorescence file

p = inputParser;
addRequired(p, 'Folder', @isfolder);
addRequired(p, 'FileData', @(x) (isnumeric(x) && ndims(x) == 3) || ischar(x));
addRequired(p, 'fMetaData', @(x) isa(x,'matlab.io.MatFile') || isstruct(x));
addOptional(p, 'cList', {}, @(x) isempty(x) || (iscell(x) && all(cellfun(@ischar, x))));
addOptional(p, 'sFreq', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'outFilename', 'HEMOCORRECTED_DATA.dat', @ischar);
parse(p, Folder, FileData, fMetaData, varargin{:});

cList = p.Results.cList;
sFreq = p.Results.sFreq;
outFilename = p.Results.outFilename;

if ~strcmp(Folder(end), filesep)
    Folder = [Folder filesep];
end

% -------------------------------------------------------------------------
% Resolve hemodynamic channel files
% -------------------------------------------------------------------------
if isempty(cList)
    cList = dir([Folder '*.dat']);
    fn = {};
    for ind = 1:numel(cList)
        if ~strcmp(cList(ind).name(1), 'f')
            fn{end+1} = cList(ind).name; %#ok<AGROW>
        end
    end

    [idx, tf] = listdlg('PromptString', {'Select channels to be used to', ...
        'compute hemodynamic correction.', ''}, ...
        'ListString', fn);

    if tf == 0
        return;
    end

    fn = fn(idx);
else
    fn = {};
    for ind = 1:numel(cList)
        tag = lower(cList{ind});

        switch tag
            case 'red'
                if exist([Folder 'rChan.dat'], 'file')
                    fn{end+1} = 'rChan.dat'; %#ok<AGROW>
                else
                    fn{end+1} = 'red.dat'; %#ok<AGROW>
                end

            case {'amber', 'yellow'}
                if exist([Folder 'yChan.dat'], 'file')
                    fn{end+1} = 'yChan.dat'; %#ok<AGROW>
                else
                    fn{end+1} = 'yellow.dat'; %#ok<AGROW>
                end

            case 'green'
                if exist([Folder 'gChan.dat'], 'file')
                    fn{end+1} = 'gChan.dat'; %#ok<AGROW>
                else
                    fn{end+1} = 'green.dat'; %#ok<AGROW>
                end

            otherwise
                [~, name, ext] = fileparts(cList{ind});
                if isempty(ext)
                    ext = '.dat';
                end
                fn{end+1} = [name, ext]; %#ok<AGROW>
        end
    end
end

% Expand to full paths
fn = fullfile(Folder, fn);

% -------------------------------------------------------------------------
% Validate optional temporal filter
% -------------------------------------------------------------------------
if sFreq
    freq = fMetaData.Freq;
    if sFreq >= freq/2
        sFreq = 0;
    end
end

% -------------------------------------------------------------------------
% Dispatch by input type
% -------------------------------------------------------------------------
if ischar(FileData)
    FileData = fullfile(Folder, FileData);
    outFilename = fullfile(Folder, outFilename);

    outFileData = HemoCorrection_lowRAMmode(outFilename, FileData, fMetaData, fn, sFreq);
    [~, outFileData, ext] = fileparts(outFileData);
    varargout{1} = [outFileData, ext];
else
    FileData = HemoCorrection_standardMode(FileData, fMetaData, fn, sFreq);

    if nargout
        varargout{1} = FileData;
    else
        [~, filename, ext] = fileparts(fMetaData.datFile);
        fid = fopen([Folder filename ext], 'w');
        fwrite(fid, FileData, 'single');
        fclose(fid);
    end
end
end


%% ========================================================================
% Local functions
% =========================================================================
function fData = HemoCorrection_standardMode(fData, fMetaData, colorList, LPcutoffFreq)
% HEMOCORRECTION_STANDARDMODE In-memory hemodynamic correction.

Ny = fMetaData.datSize(1);
Nx = fMetaData.datSize(2);
Nt = fMetaData.datLength;
Np = Nx * Ny;

% Normalize fluorescence
fData = reshape(fData, prod(fMetaData.datSize(1:2)), []);
m_fData = mean(fData, 2);
fData = (fData - m_fData) ./ m_fData;

% Design temporal filter
if LPcutoffFreq
    f = fdesign.lowpass('N,F3dB', 4, LPcutoffFreq, fMetaData.Freq);
    lpass = design(f, 'butter');
end

% Estimate chunking
numChannels = numel(colorList);
nChunks = calculateMaxChunkSize(fData, numChannels, 0.15);
chunkSizePixels = ceil(Nx / nChunks);

if nChunks > 1
    fid = {};
    for ii = 1:numChannels
        fid{ii} = fopen(colorList{ii}, 'r'); %#ok<AGROW>
    end
end

spatSigma = 1;
pad = ceil(3 * spatSigma);

h = waitbar(0, 'Fitting Hemodynamics...');
h_out = onCleanup(@() delete(h)); %#ok<NASGU>

for ii = 1:nChunks
    if nChunks == 1
        h.Name = 'Hemodynamic Correction'; drawnow()
        HemoData = zeros(numChannels, Np, Nt, 'single');
        padStart = 0;
        padStop = 0;
        indList = 1:Np;
    else
        h.Name = ['Hemodynamic Corr. (chunk ' num2str(ii) '/' num2str(nChunks) ')']; drawnow()

        pxStart = (ii - 1) * chunkSizePixels + 1;
        pxEnd   = min(ii * chunkSizePixels, Nx);
        idxPixels = pxStart:pxEnd;

        padStart = min(pad, pxStart - 1);
        padStop  = min(pad, Nx - pxEnd);
        idxPixels_with_pad = (pxStart - padStart):(pxEnd + padStop);

        [COL, ROW] = meshgrid(idxPixels, 1:Ny);
        indList = sub2ind([Ny, Nx], ROW(:), COL(:));
        HemoData = zeros(numChannels, numel(indList), Nt, 'single');
    end

    for kk = 1:numChannels
        [~, colorName, ext] = fileparts(colorList{kk});

        if nChunks == 1
            tmp = loadDatFile(colorList{kk});
        else
            tmp = spatialSlabIO('read', fid{kk}, Ny, Nx, Nt, idxPixels_with_pad, 'single');
        end

        tmp_sz = size(tmp);

        if LPcutoffFreq
            tmp = reshape(tmp, [], tmp_sz(3));
            waitbar(.99, h, ['Applying temporal filter [' colorName ext ']']); drawnow()
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
            tmp = reshape(tmp, tmp_sz);
        end

        tmp = imgaussfilt(tmp, spatSigma, 'Padding', 'symmetric');

        tmp = tmp(:, padStart+1:end-padStop, :);
        tmp = reshape(tmp, [], tmp_sz(3));

        m = mean(tmp, 2);
        tmp = (tmp - m) ./ m;

        HemoData(kk, :, :) = tmp;
    end
    clear tmp

    warning('off', 'MATLAB:rankDeficientMatrix');
    waitbar(0, h, 'Performing Hemodynamic correction...'); drawnow()

    for indP = 1:numel(indList)
        if size(HemoData,1) == 1
            X = [ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :))'];
        else
            X = [ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :))];
        end

        B = X' \ fData(indList(indP), :)';
        fData(indList(indP), :) = fData(indList(indP), :) - (X' * B)';

        if mod(indP, 500) == 0
            waitbar(indP / numel(indList), h);
        end
    end

    warning('on', 'MATLAB:rankDeficientMatrix');
end

close(h);

if exist('fid', 'var')
    for kk = 1:numel(fid)
        fclose(fid{kk});
    end
end

fData = fData .* m_fData + m_fData;
fData = reshape(fData, fMetaData.datSize(1), fMetaData.datSize(2), []);
end


function outFilename = HemoCorrection_lowRAMmode(outFilename, fluoFile, fMetaData, colorList, LPcutoffFreq)
% HEMOCORRECTION_LOWRAMMODE Disk-streamed hemodynamic correction.

f_fid = fopen(fluoFile, 'r');
c_f = onCleanup(@() safeFclose(f_fid)); %#ok<NASGU>

numChannels = numel(colorList);
h_fid = cell(1, numChannels);
c_r = cell(1, numChannels);
for k = 1:numChannels
    h_fid{k} = fopen(colorList{k}, 'r');
    c_r{k} = onCleanup(@() safeFclose(h_fid{k})); %#ok<NASGU>
end

Ny = fMetaData.datSize(1);
Nx = fMetaData.datSize(2);
Nt = fMetaData.datLength;

dataBytes = prod([fMetaData.datSize, fMetaData.datLength, getByteSize(fMetaData.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes, 2 + numel(colorList), .1);
chunkSizePixels = ceil(Nx / nChunks);

spatSigma = 1;
pad = ceil(3 * spatSigma);

if LPcutoffFreq
    f = fdesign.lowpass('N,F3dB', 4, LPcutoffFreq, fMetaData.Freq);
    lpass = design(f, 'butter');
end

preallocateDatFile(outFilename, fMetaData);
fid_out = fopen(outFilename, 'r+');
c_out = onCleanup(@() safeFclose(fid_out)); %#ok<NASGU>

h = waitbar(0, 'Fitting Hemodynamics...');
h_out = onCleanup(@() delete(h)); %#ok<NASGU>

for ii = 1:nChunks
    h.Name = ['Hemodynamic Corr. (chunk ' num2str(ii) '/' num2str(nChunks) ')']; drawnow()

    pxStart = (ii - 1) * chunkSizePixels + 1;
    pxEnd   = min(ii * chunkSizePixels, Nx);
    idxPixels = pxStart:pxEnd;

    padStart = min(pad, pxStart - 1);
    padStop  = min(pad, Nx - pxEnd);
    idxPixels_with_pad = (pxStart - padStart):(pxEnd + padStop);

    Np = numel(idxPixels) * Ny;
    HemoData = zeros(numChannels, Np, Nt, 'single');

    waitbar(.99, h, 'Reading fluo channel...'); drawnow()
    fData = spatialSlabIO('read', f_fid, Ny, Nx, Nt, idxPixels, fMetaData.Datatype);

    waitbar(.99, h, 'Normalizing fluo channel...'); drawnow()
    f_slabSz = size(fData);
    fData = reshape(fData, [], Nt);
    m_fData = mean(fData, 2);
    fData = (fData - m_fData) ./ m_fData;

    for kk = 1:numChannels
        [~, colorName, ext] = fileparts(colorList{kk});
        waitbar(.99, h, ['Reading file [' colorName ext ']']); drawnow()

        tmp = spatialSlabIO('read', h_fid{kk}, Ny, Nx, Nt, idxPixels_with_pad, fMetaData.Datatype);
        tmp_sz = size(tmp);

        if LPcutoffFreq
            waitbar(.99, h, ['Applying temporal filter [' colorName ext ']']); drawnow()
            tmp = reshape(tmp, [], tmp_sz(3));
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
            tmp = reshape(tmp, tmp_sz);
        end

        waitbar(.99, h, 'Applying spatial filter to hemodynamic data...'); drawnow()
        tmp = imgaussfilt(tmp, spatSigma, 'Padding', 'symmetric');

        tmp = tmp(:, padStart+1:end-padStop, :);
        tmp = reshape(tmp, [], tmp_sz(3));

        waitbar(.99, h, 'Normalizing hemodynamic data...'); drawnow()
        m = mean(tmp, 2);
        tmp = (tmp - m) ./ m;

        HemoData(kk, :, :) = tmp;
        clear tmp m tmp_sz
    end

    waitbar(0, h, 'Performing Hemodynamic correction...'); drawnow()
    warning('off', 'MATLAB:rankDeficientMatrix');

    for indP = 1:Np
        if size(HemoData,1) == 1
            X = [ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :))'];
        else
            X = [ones(1, Nt); linspace(0, 1, Nt); squeeze(HemoData(:, indP, :))];
        end

        B = X' \ fData(indP, :)';
        fData(indP, :) = fData(indP, :) - (X' * B)';

        if mod(indP, 500) == 0
            waitbar(indP / Np, h);
        end
    end
    clear B X HemoData
    warning('on', 'MATLAB:rankDeficientMatrix');

    fData = fData .* m_fData + m_fData;
    fData = reshape(fData, f_slabSz);

    waitbar(0.99, h, 'Writing corrected fluo to file...'); drawnow()
    spatialSlabIO('write', fid_out, Ny, Nx, Nt, idxPixels, fMetaData.Datatype, fData);
    clear fData
end

close(h);

fclose(f_fid);
fclose(fid_out);
for kk = 1:numel(h_fid)
    fclose(h_fid{kk});
end
end