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
%       - Metadata are resolved through loadMetaData for file inputs.
%       - Numeric inputs are validated against the timelines described by
%         AcqInfos.mat.
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
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
if ~strcmp(SaveFolder(end), filesep)
    SaveFolder = [SaveFolder filesep];
end

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
    fMetaData = iResolveNumericYXTMetadata(data, acq);
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

if nChunks > 1
    fid = {};
    for ii = 1:numChannels
        fid{ii} = fopen(colorList{ii}, 'r'); %#ok<AGROW>
        assert(fid{ii} ~= -1, ...
            'Umitoolbox:HemoCorrection:FileOpenFailed', ...
            'Could not open hemodynamic reference file "%s".', colorList{ii});
    end
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
            tmp = spatialSlabIO('read', fid{kk}, Ny, Nx, cMetaData{kk}.datLength, ...
                idxPixels_with_pad, cMetaData{kk}.Datatype);
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
assert(f_fid ~= -1, ...
    'Umitoolbox:HemoCorrection:FileOpenFailed', ...
    'Could not open fluorescence file "%s".', fluoFile);
c_f = onCleanup(@() safeFclose(f_fid)); 

numChannels = numel(colorList);
h_fid = cell(1, numChannels);
c_r = cell(1, numChannels);
cMetaData = cell(1, numChannels);
maxNt = fMetaData.datLength;
for k = 1:numChannels
    h_fid{k} = fopen(colorList{k}, 'r');
    assert(h_fid{k} ~= -1, ...
        'Umitoolbox:HemoCorrection:FileOpenFailed', ...
        'Could not open hemodynamic reference file "%s".', colorList{k});
    c_r{k} = onCleanup(@() safeFclose(h_fid{k})); 
    cMetaData{k} = iNormalizeDatMeta(loadMetaData(colorList{k}));
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

spatSigma = 1;
pad = ceil(3 * spatSigma);

preallocateDatFile(outFilename, [Ny, Nx, Nt], fMetaData.Datatype);
fid_out = fopen(outFilename, 'r+');
assert(fid_out ~= -1, ...
    'Umitoolbox:HemoCorrection:FileOpenFailed', ...
    'Could not create output file "%s".', outFilename);
c_out = onCleanup(@() safeFclose(fid_out)); 

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
    fData = spatialSlabIO('read', f_fid, Ny, Nx, Nt, idxPixels, fMetaData.Datatype);

    waitbar(.99, h, 'Normalizing fluo channel...'); drawnow()
    f_slabSz = size(fData);
    fData = reshape(fData, [], Nt);
    m_fData = mean(fData, 2);
    fData = (fData - m_fData) ./ m_fData;

    for kk = 1:numChannels
        [~, colorName, ext] = fileparts(colorList{kk});
        waitbar(.99, h, ['Reading file [' colorName ext ']']); drawnow()

        tmp = spatialSlabIO('read', h_fid{kk}, Ny, Nx, cMetaData{kk}.datLength, ...
            idxPixels_with_pad, cMetaData{kk}.Datatype);

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


function fMetaData = iResolveNumericYXTMetadata(data, AcqInfoStream)
%IRESOLVENUMERICYXTMETADATA Build file-like metadata for numeric YXT input.

assert(isfield(AcqInfoStream, 'Height') && isfield(AcqInfoStream, 'Width'), ...
    'Umitoolbox:HemoCorrection:InvalidAcqInfos', ...
    'AcqInfoStream must contain Height and Width.');

height = double(AcqInfoStream.Height);
width = double(AcqInfoStream.Width);
assert(isequal([size(data,1), size(data,2)], [height, width]), ...
    'Umitoolbox:HemoCorrection:InvalidNumericInput', ...
    'Numeric fluorescence input does not match AcqInfos.mat Height/Width.');

timelineInfo = resolveDatTimeline(size(data,3), AcqInfoStream);

fMetaData = struct();
fMetaData.datSize = [height, width];
fMetaData.Height = height;
fMetaData.Width = width;
fMetaData.datLength = double(timelineInfo.Length);
fMetaData.Length = double(timelineInfo.Length);
fMetaData.Freq = double(timelineInfo.FrameRateHz);
fMetaData.FrameRateHz = double(timelineInfo.FrameRateHz);
fMetaData.Datatype = 'single';
fMetaData.dim_names = {'Y','X','T'};
end


function meta = iNormalizeDatMeta(meta)
%INORMALIZEDATMETA Ensure loadMetaData output has legacy-compatible fields.

if ~isfield(meta, 'datSize') || isempty(meta.datSize)
    meta.datSize = [double(meta.Height), double(meta.Width)];
else
    meta.datSize = double(meta.datSize(:).');
end

if ~isfield(meta, 'datLength') || isempty(meta.datLength)
    meta.datLength = double(meta.Length);
else
    meta.datLength = double(meta.datLength);
end

if ~isfield(meta, 'Freq') || isempty(meta.Freq)
    meta.Freq = double(meta.FrameRateHz);
else
    meta.Freq = double(meta.Freq);
end

if ~isfield(meta, 'Datatype') || isempty(meta.Datatype)
    meta.Datatype = 'single';
else
    meta.Datatype = char(string(meta.Datatype));
end

if ~isfield(meta, 'Height') || isempty(meta.Height)
    meta.Height = meta.datSize(1);
end

if ~isfield(meta, 'Width') || isempty(meta.Width)
    meta.Width = meta.datSize(2);
end
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
