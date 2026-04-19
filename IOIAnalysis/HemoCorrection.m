function varargout = HemoCorrection(data, SaveFolder, varargin)
%HEMOCORRECTION Remove hemodynamic fluctuations from fluorescence signals.
%
%   out = HemoCorrection(data, SaveFolder)
%   out = HemoCorrection(data, SaveFolder, 'ChannelList', {'red'})
%   out = HemoCorrection(data, SaveFolder, 'LowPassFreq', 2)
%
%   This function performs pixel-wise hemodynamic regression on
%   fluorescence image time series using one or more intrinsic/hemodynamic
%   reference channels located in the same SaveFolder.
%
%   Supported inputs:
%       - Numeric fluorescence data with dimensions Y X T
%       - Raw .dat fluorescence filename
%
%   Inputs:
%       data       - 3D numeric YXT array or raw .dat filename.
%       SaveFolder - Folder containing AcqInfos.mat and the hemodynamic
%                    reference channel files.
%
%   Name-Value parameters:
%       'ChannelList' - Cell array of channel tags or filenames. Supported
%                       tags: 'red', 'green', 'yellow', 'amber'.
%                       If empty, a list dialog is shown.
%       'LowPassFreq' - Low-pass cutoff frequency in Hz. Set to 0 to
%                       disable filtering.
%
%   Output:
%       - Numeric input  -> corrected numeric YXT data
%       - Raw .dat input -> corrected output filename
%
%   Notes:
%       - Only numeric arrays and raw .dat files are supported.
%       - Metadata are resolved directly from AcqInfos.mat.
%       - Reference channels must be stored in the same SaveFolder.

% Default output for pipeline management:
default_Output = 'fluoHemoCorr.dat'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    varargout{1} = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'HemoCorrection';
addRequired(p, 'data', @(x) (isnumeric(x) && ndims(x) == 3) || ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addParameter(p, 'ChannelList', {}, @(x) isempty(x) || (iscell(x) && all(cellfun(@(c) ischar(c) || (isstring(c) && isscalar(c)), x))));
addParameter(p, 'LowPassFreq', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, data, SaveFolder, varargin{:});

channelList = p.Results.ChannelList;
lowPassFreq = double(p.Results.LowPassFreq);

acqFile = fullfile(SaveFolder, 'AcqInfos.mat');
assert(isfile(acqFile), ...
    'Umitoolbox:HemoCorrection:missingAcqInfos', ...
    'AcqInfos.mat was not found in "%s".', SaveFolder);

md = load(acqFile, 'AcqInfoStream');
assert(isfield(md, 'AcqInfoStream') && isstruct(md.AcqInfoStream), ...
    'Umitoolbox:HemoCorrection:invalidAcqInfos', ...
    'AcqInfos.mat does not contain a valid AcqInfoStream structure.');
acq = md.AcqInfoStream;

assert(isfield(acq, 'Height') && isfield(acq, 'Width') && isfield(acq, 'Length') && isfield(acq, 'FrameRateHz'), ...
    'Umitoolbox:HemoCorrection:missingAcqFields', ...
    'AcqInfoStream must contain Height, Width, Length, and FrameRateHz.');

fMetaData = struct();
fMetaData.datSize = [double(acq.Height), double(acq.Width)];
fMetaData.datLength = double(acq.Length);
fMetaData.Freq = double(acq.FrameRateHz);
fMetaData.Datatype = 'single';

if ischar(data) || (isstring(data) && isscalar(data))
    [~,~,ext] = fileparts(char(string(data)));
    assert(strcmpi(ext, '.dat'), ...
        'Umitoolbox:HemoCorrection:unsupportedInputFile', ...
        'Only raw .dat fluorescence files are supported.');
    fileData = fullfile(SaveFolder, char(string(data)));
    assert(isfile(fileData), ...
        'Umitoolbox:HemoCorrection:missingInputFile', ...
        'Input fluorescence file was not found: %s', fileData);
else
    fileData = data;
end

assert(lowPassFreq < fMetaData.Freq/2 || lowPassFreq == 0, ...
    'Umitoolbox:HemoCorrection:invalidLowPassFreq', ...
    'LowPassFreq must be smaller than the Nyquist frequency.');

% Resolve hemodynamic reference channels. Keep the interactive list dialog
% for standalone use when the caller does not provide ChannelList.
if isempty(channelList)
    datFiles = dir(fullfile(SaveFolder, '*.dat'));
    available = {};
    for iFile = 1:numel(datFiles)
        thisName = datFiles(iFile).name;
        if ~(ischar(data) || (isstring(data) && isscalar(data))) || ~strcmpi(thisName, char(string(data)))
            available{end+1} = thisName; %#ok<AGROW>
        end
    end

    assert(~isempty(available), ...
        'Umitoolbox:HemoCorrection:noReferenceChannelsFound', ...
        'No hemodynamic reference channel files were found in "%s".', SaveFolder);

    [idx, tf] = listdlg('PromptString', {'Select channels to be used to', ...
        'compute hemodynamic correction.', ''}, ...
        'ListString', available);

    if tf == 0
        error('Umitoolbox:HemoCorrection:selectionCancelled', ...
            'Hemodynamic channel selection was cancelled by the user.');
    end

    resolvedFiles = available(idx);
else
    resolvedFiles = cell(1, numel(channelList));
    for iChan = 1:numel(channelList)
        tag = lower(char(string(channelList{iChan})));
        switch tag
            case 'red'
                resolvedFiles{iChan} = 'red.dat';
            case {'yellow', 'amber'}
                resolvedFiles{iChan} = 'yellow.dat';
            case 'green'
                resolvedFiles{iChan} = 'green.dat';
            otherwise
                [~, name, ext] = fileparts(char(string(channelList{iChan})));
                if isempty(ext)
                    ext = '.dat';
                end
                resolvedFiles{iChan} = [name ext];
        end
    end
end

resolvedFiles = fullfile(SaveFolder, resolvedFiles);
for iFile = 1:numel(resolvedFiles)
    assert(isfile(resolvedFiles{iFile}), ...
        'Umitoolbox:HemoCorrection:missingReferenceChannel', ...
        'Hemodynamic reference channel was not found in SaveFolder: %s', resolvedFiles{iFile});
end

if ischar(data) || (isstring(data) && isscalar(data))
    outPath = fullfile(SaveFolder, default_Output);
    if isfile(outPath)
        [folderPath, baseName, ext] = fileparts(outPath);
        outPath = fullfile(folderPath, [baseName '_preallocData' ext]);
    end

    outFileData = HemoCorrection_lowRAMmode(outPath, fileData, fMetaData, resolvedFiles, lowPassFreq);
    [~, baseName, ext] = fileparts(outFileData);
    varargout{1} = [baseName ext];
else
    correctedData = HemoCorrection_standardMode(fileData, fMetaData, resolvedFiles, lowPassFreq);
    if nargout > 0
        varargout{1} = correctedData;
    else
        error('Umitoolbox:HemoCorrection:noOutputForNumericInput', ...
            ['Numeric input mode requires an output argument in the refactored version. ' ...
             'Use a raw .dat input to write corrected data to disk.']);
    end
end

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'HemoCorrection', ...
            'Remove hemodynamic fluctuations from fluorescence image time series.');

        info = PipelineManager.addInput(info, ...
            'data', ...
            {'ImageTimeSeries'}, ...
            'Fluorescence image time series input.', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            {'parameter'}, ...
            'Folder containing AcqInfos.mat and hemodynamic channels.', ...
            'kind', 'parameter', ...
            'position', 2, ...
            'callType', 'positional', ...
            'default', '', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'ChannelList', ...
            'parameter', ...
            'List of hemodynamic channels or filenames.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', {{}}, ...
            'dataType', 'cell');

        info = PipelineManager.addInput(info, ...
            'LowPassFreq', ...
            'parameter', ...
            'Low-pass cutoff frequency in Hz.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', 0, ...
            'dataType', 'numeric');

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            {'ImageTimeSeries'}, ...
            'data', ...
            'Corrected fluorescence image time series.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

%% ========================================================================
% Local functions
% =========================================================================
function fData = HemoCorrection_standardMode(fData, fMetaData, colorList, LPcutoffFreq)
%HEMOCORRECTION_STANDARDMODE In-memory hemodynamic correction.

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
            tmp = loadData(colorList{kk});
        else
            tmp = spatialSlabIO('read', fid{kk}, Ny, Nx, Nt, idxPixels_with_pad, 'single');
        end

        tmp_sz = size(tmp);

        if LPcutoffFreq
            tmp = reshape(tmp, [], tmp_sz(3));
            waitbar(.99, h, ['Applying temporal filter [' colorName ext ']']); drawnow()
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))' ;
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
%HEMOCORRECTION_LOWRAMMODE Disk-streamed hemodynamic correction.

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

preallocateDatFile(outFilename, [Ny, Nx, Nt], fMetaData.Datatype);
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
            tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))' ;
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
