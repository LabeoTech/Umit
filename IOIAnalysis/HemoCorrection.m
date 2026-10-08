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
%       SaveFolder - Folder containing the hemodynamic reference channel
%                    files (AcqInfos.mat is not read).
%
%   Name-Value parameters:
%       'ChannelList' - Cell array of channel tags or filenames. Supported
%                       tags: 'red', 'green', 'yellow', 'amber'.
%                       If empty, a list dialog is shown.
%       'LowPassFreq' - Low-pass cutoff frequency in Hz. Set to 0 to
%                       disable filtering.
%       'FrameRateHz' - Frame rate of numeric fluorescence input (Hz),
%                       required for numeric input (the data's own rate;
%                       run_HemoCorrection passes it).
%
%   Output:
%       - Numeric input  -> corrected numeric YXT data
%       - Raw .dat input -> corrected output filename
%
%   Notes:
%       - Metadata are resolved through loadMetaData for file inputs.
%       - Numeric inputs take their frame rate from 'FrameRateHz'
%         (AcqInfos.mat is not used for it).
%       - Hemodynamic reference channels are temporally resampled to match
%         the fluorescence timeline before regression.

% Default output for pipeline management:
default_Output = 'fluoHemoCorr.dat'; 

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
% Frame rate of numeric fluorescence input (.dat header Phase 6b-2): the
% data's own rate, never AcqInfos.mat.
addParameter(p, 'FrameRateHz', []);
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
if ~strcmp(SaveFolder(end), filesep)
    SaveFolder = [SaveFolder filesep];
end

channelList = p.Results.ChannelList;
lowPassFreq = double(p.Results.LowPassFreq);

if ischar(data) || (isstring(data) && isscalar(data))
    [fileData, inputFileName] = iResolveFileInSaveFolder(SaveFolder, data);
    [~,~,ext] = fileparts(fileData);
    assert(strcmpi(ext, '.dat'), ...
        'Umitoolbox:HemoCorrection:unsupportedInputFile', ...
        'Only raw .dat fluorescence files are supported.');
    assert(isfile(fileData), ...
        'Umitoolbox:HemoCorrection:missingInputFile', ...
        'Input fluorescence file was not found: %s', fileData);
    fMetaData = iNormalizeDatMeta(loadMetaData(fileData));
else
    fileData = data;
    inputFileName = '';
    fMetaData = iResolveNumericYXTMetadata(data, ...
        resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, data, 'HemoCorrection'));
end

assert(lowPassFreq < fMetaData.Freq/2 || lowPassFreq == 0, ...
    'Umitoolbox:HemoCorrection:invalidLowPassFreq', ...
    'LowPassFreq must be smaller than the fluorescence Nyquist frequency.');

% Resolve hemodynamic reference channels. Keep the interactive list dialog
% for standalone use when the caller does not provide ChannelList.
if isempty(channelList)
    datFiles = dir(fullfile(SaveFolder, '*.dat'));
    available = {};
    for iFile = 1:numel(datFiles)
        thisName = datFiles(iFile).name;
        if isempty(inputFileName) || ~strcmpi(thisName, inputFileName)
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
    % Write through a fixed-name scratch file, then move it onto the
    % declared output. Renaming the output when it already exists would
    % make every pipeline re-run write to a different file and leave the
    % stale original in place.
    outPath = fullfile(SaveFolder, default_Output);
    [~, baseName, ext] = fileparts(default_Output);
    tmpPath = fullfile(SaveFolder, [baseName '_writing' ext]);

    HemoCorrection_lowRAMmode(tmpPath, baseName, fileData, fMetaData, resolvedFiles, lowPassFreq);

    [moveOk, moveMsg] = movefile(tmpPath, outPath, 'f');
    assert(moveOk, 'Umitoolbox:HemoCorrection:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', tmpPath, outPath, moveMsg);

    varargout{1} = default_Output;
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
            'Folder containing the hemodynamic channels.', ...
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

numChannels = numel(colorList);
cMetaData = cell(1, numChannels);
maxNt = Nt;
for kk = 1:numChannels
    cMetaData{kk} = iNormalizeDatMeta(loadMetaData(colorList{kk}));
    iValidateSpatialMatch(cMetaData{kk}, fMetaData, colorList{kk});
    maxNt = max(maxNt, cMetaData{kk}.datLength);
end

% Assert that every reference channel spans the same recording duration as
% the fluorescence channel before any temporal resampling. Interpolation is
% only valid for channels covering the same acquisition interval.
durationTolSec = 1e-3;
fluoDurationSec = double(fMetaData.datLength) / double(fMetaData.Freq);
assert(isfinite(fluoDurationSec) && fluoDurationSec > 0, ...
    'Umitoolbox:HemoCorrection:InvalidChannelDuration', ...
    'Fluorescence channel has invalid duration metadata: Length=%g, FrameRateHz=%g.', ...
    double(fMetaData.datLength), double(fMetaData.Freq));
for kk = 1:numChannels
    refDurationSec = double(cMetaData{kk}.datLength) / double(cMetaData{kk}.Freq);
    assert(isfinite(refDurationSec) && refDurationSec > 0, ...
        'Umitoolbox:HemoCorrection:InvalidChannelDuration', ...
        'Reference channel "%s" has invalid duration metadata: Length=%g, FrameRateHz=%g.', ...
        colorList{kk}, double(cMetaData{kk}.datLength), double(cMetaData{kk}.Freq));

    assert(abs(refDurationSec - fluoDurationSec) <= durationTolSec, ...
        'Umitoolbox:HemoCorrection:DurationMismatch', ...
        ['Reference channel "%s" does not span the same recording duration as ' ...
         'the fluorescence channel. Reference: Length=%g, FrameRateHz=%g, ' ...
         'Duration=%0.6f s. Fluorescence duration=%0.6f s.'], ...
        colorList{kk}, double(cMetaData{kk}.datLength), double(cMetaData{kk}.Freq), ...
        refDurationSec, fluoDurationSec);
end

% Normalize fluorescence
fData = reshape(fData, prod(fMetaData.datSize(1:2)), []);
m_fData = mean(fData, 2);
fData = (fData - m_fData) ./ m_fData;

% Estimate chunking. Use the largest temporal input because reference
% channels can have a different timeline from the fluorescence channel.
dataBytes = max(numel(fData) * getByteSize(class(fData)), ...
    prod(fMetaData.datSize) * maxNt * getByteSize('single'));
nChunks = calculateMaxChunkSize(dataBytes, 2 + numChannels, 0.15);
chunkSizePixels = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkSizePixels);

if nChunks > 1
    refSlab = cell(1, numChannels);
    for ii = 1:numChannels
        refSlab{ii} = spatialSlabIO('open', colorList{ii}, 'Info', cMetaData{ii});
    end
    c_ref = onCleanup(@() cellfun(@(r) spatialSlabIO('close', r), refSlab));
end

spatSigma = 1;
pad = ceil(3 * spatSigma);

h = waitbar(0, 'Fitting Hemodynamics...');
h_out = onCleanup(@() delete(h)); 

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
            tmp = spatialSlabIO('read', refSlab{kk}, idxPixels_with_pad);
        end

        tmp = iResampleHemoToFluoTimeline(tmp, cMetaData{kk}, fMetaData, LPcutoffFreq);
        tmp_sz = size(tmp);

        tmp = imgaussfilt(tmp, spatSigma, 'Padding', 'symmetric');

        tmp = tmp(:, padStart+1:end-padStop, :);
        tmp = reshape(tmp, [], tmp_sz(3));

        m = mean(tmp, 2);
        tmp = (tmp - m) ./ m;

        HemoData(kk, :, :) = tmp;
        clear tmp m tmp_sz
        waitbar(.99, h, ['Loaded hemodynamic channel [' colorName ext ']']); drawnow()
    end

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

if exist('refSlab', 'var')
    for kk = 1:numel(refSlab)
        spatialSlabIO('close', refSlab{kk});
    end
end

fData = fData .* m_fData + m_fData;
fData = reshape(fData, fMetaData.datSize(1), fMetaData.datSize(2), []);
end


function outFilename = HemoCorrection_lowRAMmode(outFilename, outBaseName, fluoFile, fMetaData, colorList, LPcutoffFreq)
%HEMOCORRECTION_LOWRAMMODE Disk-streamed hemodynamic correction.
%
% OUTFILENAME is the scratch file written here; OUTBASENAME is the name of
% the declared output it is moved onto (the header channelName).

fluoSlab = spatialSlabIO('open', fluoFile, 'Info', fMetaData);
c_f = onCleanup(@() spatialSlabIO('close', fluoSlab)); 

numChannels = numel(colorList);
refSlab = cell(1, numChannels);
c_r = cell(1, numChannels);
cMetaData = cell(1, numChannels);
maxNt = fMetaData.datLength;
for k = 1:numChannels
    cMetaData{k} = iNormalizeDatMeta(loadMetaData(colorList{k}));
    refSlab{k} = spatialSlabIO('open', colorList{k}, 'Info', cMetaData{k});
    c_r{k} = onCleanup(@() spatialSlabIO('close', refSlab{k})); 
    iValidateSpatialMatch(cMetaData{k}, fMetaData, colorList{k});
    maxNt = max(maxNt, cMetaData{k}.datLength);
end

% Assert that every reference channel spans the same recording duration as
% the fluorescence channel before any temporal resampling. Interpolation is
% only valid for channels covering the same acquisition interval.
durationTolSec = 1e-3;
fluoDurationSec = double(fMetaData.datLength) / double(fMetaData.Freq);
assert(isfinite(fluoDurationSec) && fluoDurationSec > 0, ...
    'Umitoolbox:HemoCorrection:InvalidChannelDuration', ...
    'Fluorescence channel has invalid duration metadata: Length=%g, FrameRateHz=%g.', ...
    double(fMetaData.datLength), double(fMetaData.Freq));
for k = 1:numChannels
    refDurationSec = double(cMetaData{k}.datLength) / double(cMetaData{k}.Freq);
    assert(isfinite(refDurationSec) && refDurationSec > 0, ...
        'Umitoolbox:HemoCorrection:InvalidChannelDuration', ...
        'Reference channel "%s" has invalid duration metadata: Length=%g, FrameRateHz=%g.', ...
        colorList{k}, double(cMetaData{k}.datLength), double(cMetaData{k}.Freq));

    assert(abs(refDurationSec - fluoDurationSec) <= durationTolSec, ...
        'Umitoolbox:HemoCorrection:DurationMismatch', ...
        ['Reference channel "%s" does not span the same recording duration as ' ...
         'the fluorescence channel. Reference: Length=%g, FrameRateHz=%g, ' ...
         'Duration=%0.6f s. Fluorescence duration=%0.6f s.'], ...
        colorList{k}, double(cMetaData{k}.datLength), double(cMetaData{k}.Freq), ...
        refDurationSec, fluoDurationSec);
end

Ny = fMetaData.datSize(1);
Nx = fMetaData.datSize(2);
Nt = fMetaData.datLength;

dataBytes = prod([fMetaData.datSize, maxNt, getByteSize(fMetaData.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes, 2 + numel(colorList), .1);
chunkSizePixels = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkSizePixels);

spatSigma = 1;
pad = ceil(3 * spatSigma);

slabOut = spatialSlabIO('create', outFilename, datHeaderFromInfo(fMetaData, outBaseName));
c_out = onCleanup(@() spatialSlabIO('close', slabOut)); 

h = waitbar(0, 'Fitting Hemodynamics...');
h_out = onCleanup(@() delete(h)); 

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
    fData = spatialSlabIO('read', fluoSlab, idxPixels);

    waitbar(.99, h, 'Normalizing fluo channel...'); drawnow()
    f_slabSz = size(fData);
    fData = reshape(fData, [], Nt);
    m_fData = mean(fData, 2);
    fData = (fData - m_fData) ./ m_fData;

    for kk = 1:numChannels
        [~, colorName, ext] = fileparts(colorList{kk});
        waitbar(.99, h, ['Reading file [' colorName ext ']']); drawnow()

        tmp = spatialSlabIO('read', refSlab{kk}, idxPixels_with_pad);

        waitbar(.99, h, ['Resampling hemodynamic file [' colorName ext ']']); drawnow()
        tmp = iResampleHemoToFluoTimeline(tmp, cMetaData{kk}, fMetaData, LPcutoffFreq);
        tmp_sz = size(tmp);

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
    spatialSlabIO('write', slabOut, idxPixels, fData);
    clear fData
end

close(h);

spatialSlabIO('close', fluoSlab);
spatialSlabIO('finalize', slabOut);
for kk = 1:numel(refSlab)
    spatialSlabIO('close', refSlab{kk});
end
end


function tmp = iResampleHemoToFluoTimeline(tmp, cMetaData, fMetaData, LPcutoffFreq)
%IRESAMPLEHEMOTOFLUOTIMELINE Match one hemodynamic slab to fluorescence T.

NtHemo = cMetaData.datLength;
NtFluo = fMetaData.datLength;
freqHemo = cMetaData.Freq;
freqFluo = fMetaData.Freq;

if size(tmp,3) ~= NtHemo
    error('Umitoolbox:HemoCorrection:InvalidSlabLength', ...
        'Input hemodynamic slab length does not match its metadata.');
end

cutoffFreq = 0;
if LPcutoffFreq > 0
    cutoffFreq = LPcutoffFreq;
elseif freqHemo > freqFluo && NtHemo > NtFluo
    cutoffFreq = 0.45 * freqFluo;
end

if cutoffFreq > 0 && cutoffFreq < freqHemo/2
    sz = size(tmp);
    tmp = reshape(tmp, [], sz(3));
    f = fdesign.lowpass('N,F3dB', 4, cutoffFreq, freqHemo);
    lpass = design(f, 'butter');
    tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))' ;
    tmp = reshape(tmp, sz);
end

if NtHemo ~= NtFluo
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


function fMetaData = iResolveNumericYXTMetadata(data, frameRateHz)
%IRESOLVENUMERICYXTMETADATA Build file-like metadata for numeric YXT input.
%
% FRAMERATEHZ is the data's own frame rate (the FrameRateHz Name-Value).

% Y and X come from the array itself: AcqInfos.mat Height/Width is the raw
% acquisition size and no longer matches aligned data. The reference
% channels are checked against this size (iValidateSpatialMatch).
height = double(size(data, 1));
width = double(size(data, 2));

nT = double(size(data, 3));

fMetaData = struct();
fMetaData.datSize = [height, width];
fMetaData.Height = height;
fMetaData.Width = width;
fMetaData.datLength = nT;
fMetaData.Length = nT;
fMetaData.Freq = double(frameRateHz);
fMetaData.FrameRateHz = double(frameRateHz);
fMetaData.Datatype = 'single';
fMetaData.dim_names = {'Y','X','T'};
end


function meta = iNormalizeDatMeta(meta)
%INORMALIZEDATMETA Add this function's internal size fields to loadMetaData output.
%
% The internal fields (datSize = [Y X], datLength, Freq, Datatype, Height,
% Width) are derived from the .dat Info schema only. The schema fields are
% kept so the struct can be passed to spatialSlabIO('open', ..., 'Info', meta).

ny = datAxisSize(meta, 'Y');
nx = datAxisSize(meta, 'X');
meta.datSize = [ny, nx];
meta.datLength = datAxisSize(meta, 'T');
meta.Freq = double(meta.frameRateHz);
meta.Datatype = char(string(meta.dataClass));
meta.Height = ny;
meta.Width = nx;
end


function iValidateSpatialMatch(cMetaData, fMetaData, fileName)
%IVALIDATESPATIALMATCH Validate reference/fluorescence spatial dimensions.

assert(isequal(double(cMetaData.datSize(1:2)), double(fMetaData.datSize(1:2))), ...
    'Umitoolbox:HemoCorrection:SpatialMismatch', ...
    'Hemodynamic reference "%s" has incompatible spatial dimensions.', fileName);
end


function [filePath, fileName] = iResolveFileInSaveFolder(SaveFolder, fileInput)
%IRESOLVEFILEINSAVEFOLDER Resolve a filename or full path to a .dat file.

fileInput = char(string(fileInput));
if isfile(fileInput)
    filePath = fileInput;
else
    filePath = fullfile(SaveFolder, fileInput);
end

[~, baseName, ext] = fileparts(filePath);
if isempty(ext)
    ext = '.dat';
    filePath = [filePath ext];
end
fileName = [baseName ext];
end
