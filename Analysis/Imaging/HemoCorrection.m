function outData = HemoCorrection(data, SaveFolder, varargin)
%HEMOCORRECTION Apply hemodynamic correction to a fluorescence channel.
%
%   outData = HemoCorrection(data, SaveFolder)
%   outData = HemoCorrection(data, SaveFolder, 'Algorithm', 'Ratiometric', ...)
%
%   This function supports two correction approaches:
%
%   1) LinearRegression
%      Pixel-wise linear regression of fluorescence onto one or more
%      reflectance channels.
%
%   2) Ratiometric
%      Single-reference correction based on normalized fluorescence and
%      normalized reference reflectance.
%
%   Supported input modes:
%       - numeric fluorescence array [Y, X, T] or [Y, X, T, E]
%       - fluorescence .dat filename with axes Y-X-T or Y-X-T-E
%       UMT structs and .umt files are not supported.
%
%   Inputs:
%       data       - Fluorescence array or .dat filename (see above).
%       SaveFolder - Folder containing the reference files (continuous
%                    Y-X-T channels such as red.dat), and events.mat for
%                    event-split (E axis) fluorescence. AcqInfos.mat is not
%                    read.
%
%   Name-Value parameters:
%       'Algorithm' - 'LinearRegression' or 'Ratiometric'
%       'Red'       - logical scalar
%       'Green'     - logical scalar
%       'Amber'     - logical scalar
%       'Other'     - custom channel filename
%       'LowPassFreq' - low-pass cutoff (Hz) applied to the reference
%                     channels before the regression. LinearRegression on
%                     Y-X-T fluorescence only; it must be 0 for Ratiometric
%                     and for event-split (Y-X-T-E) fluorescence. It must be
%                     below the fluorescence Nyquist frequency. 0 (default)
%                     applies only the automatic anti-alias filter of a
%                     reference faster than the fluorescence. A reference
%                     whose Nyquist frequency is at or below the cutoff
%                     (a slower reference) is not filtered, with a warning.
%       'FrameRateHz' - frame rate of the FLUORESCENCE data (Hz), never of a
%                     reference channel. PipelineManager injects it from the
%                     data input; a .dat input's header provides it
%                     otherwise. Numeric input needs it explicitly;
%                     AcqInfos.mat is not used. Every reference channel uses
%                     its own .dat header rate and is resampled to this one.
%
%   Output:
%       - numeric input: corrected fluorescence array (same size)
%       - .dat input   : hemoCorr_fluo.dat (same axes as the input)
%
%   Event-split fluorescence (Y-X-T-E):
%       Each E slice is a trial that is corrected independently: its own
%       fluorescence and reference means, and its own regression over the
%       trial's frames. The E axis must hold one slice per event instance of
%       events.mat (as saved by split_data_by_event). The continuous
%       reference channels are cut by the same events: every trial's frame
%       times, taken from events.mat at the fluorescence frame rate, are
%       sampled in each reference channel by linear interpolation, so a
%       reference with a different frame rate is up- or down-sampled onto the
%       fluorescence timeline of every trial. A reference faster than the
%       fluorescence is low-pass filtered (0.45 x the fluorescence rate)
%       before it is sampled. Continuous (Y-X-T) fluorescence resamples the
%       whole reference to the fluorescence length the same way.

% Default output for pipeline management:
default_Output = 'hemoCorr_fluo.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'HemoCorrection';
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addParameter(p, 'Algorithm', 'LinearRegression', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Red', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Green', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Amber', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Other', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'LowPassFreq', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'FrameRateHz', []);
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
if ~strcmp(SaveFolder(end), filesep)
    SaveFolder = [SaveFolder filesep];
end

algorithm = char(string(p.Results.Algorithm));
useRed = p.Results.Red;
useGreen = p.Results.Green;
useAmber = p.Results.Amber;
otherChan = char(string(p.Results.Other));
lowPassFreq = double(p.Results.LowPassFreq);

assert(ismember(lower(algorithm), {'linearregression', 'ratiometric'}), ...
    'Umitoolbox:HemoCorrection:InvalidInput', ...
    'Unknown correction algorithm "%s".', algorithm);

% Y-X-T runs through the continuous paths below; Y-X-T-E (event-split)
% fluorescence is corrected trial by trial (localCorrectEventSplit).
hasE = iClassifyInput(data, SaveFolder);

assert(lowPassFreq == 0 || (~hasE && strcmpi(algorithm, 'linearregression')), ...
    'Umitoolbox:HemoCorrection:UnsupportedLowPassFreq', ...
    ['LowPassFreq is only supported by the LinearRegression algorithm on ' ...
     'Y-X-T fluorescence. Set LowPassFreq to 0 for Ratiometric correction ' ...
     'or event-split (Y-X-T-E) fluorescence.']);

% Build selected channel list.
channelList = {};
if useRed
    channelList{end+1} = 'red';
end
if useGreen
    channelList{end+1} = 'green';
end
if useAmber
    % Normalize to the canonical name also used by HemoCompute
    % (UMITRigStore.normalizeIlluminationName maps 'amber' -> 'yellow').
    channelList{end+1} = 'yellow';
end
if ~isempty(otherChan)
    channelList{end+1} = otherChan;
end

fprintf('Performing hemodynamic correction in fluo channel using %s algorithm...\n', ...
    algorithm);

% Frame rate of the fluorescence data itself (explicit or injected, else
% the .dat header); AcqInfos.mat is not used. Resolved where needed, after
% the input checks.
resolveRate = @() resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, ...
    iDataFileOrArray(SaveFolder, data), mfilename);

if hasE
    outData = localCorrectEventSplit(data, SaveFolder, algorithm, channelList, ...
        p.Results.FrameRateHz, default_Output);
    fprintf('Finished hemodynamic correction.\n');
    return
end

switch lower(algorithm)
    case 'linearregression'
        if isnumeric(data)
            outData = localLinearRegression(data, SaveFolder, channelList, ...
                lowPassFreq, resolveRate(), default_Output);
        else
            % File input: the frame rate comes from the .dat header.
            outData = localLinearRegression(data, SaveFolder, channelList, ...
                lowPassFreq, [], default_Output);
        end

    case 'ratiometric'
        localAssertSingleReference(channelList);

        refFile = localResolveReferenceFile(SaveFolder, channelList{1});
        fprintf('Using channel "%s" in hemodynamic correction...\n', refFile);

        if ischar(data) || (isstring(data) && isscalar(data))
            [fluoPath, fluoName] = localResolveFileInSaveFolder(SaveFolder, data);
            fluoMeta = localNormalizeDatMeta(loadMetaData(fluoPath));
            frameRateHz = resolveRate();
            fluoMeta.Freq = frameRateHz;
            fluoMeta.FrameRateHz = frameRateHz;
            outData = localRatiometricLowRAM(SaveFolder, fluoName, fluoMeta, refFile, default_Output);
        else
            fluoMeta = localResolveNumericYXTMetadata(data, resolveRate());
            outData = localRatiometricStandard(SaveFolder, data, fluoMeta, refFile);
        end

    otherwise
        error('Umitoolbox:HemoCorrection:InvalidInput', ...
            'Unknown correction algorithm "%s".', algorithm);
end

fprintf('Finished hemodynamic correction.\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'HemoCorrection', ...
            'Apply hemodynamic correction to a fluorescence channel.');

        info = PipelineManager.addInput(info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'Fluorescence image time series input (YXT or event-split YXTE).', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing the reference channels (and events.mat for event-split data).', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'Algorithm', ...
            'parameter', ...
            'Correction algorithm.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', 'LinearRegression', ...
            'allowed', {'LinearRegression','Ratiometric'}, ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'Red', ...
            'parameter', ...
            'Use red.dat as reference.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', true, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Green', ...
            'parameter', ...
            'Use green.dat as reference.', ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'namevalue', ...
            'default', true, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Amber', ...
            'parameter', ...
            'Use yellow.dat as reference.', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', true, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Other', ...
            'parameter', ...
            'Custom reference filename.', ...
            'kind', 'parameter', ...
            'position', 7, ...
            'callType', 'namevalue', ...
            'default', '', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'LowPassFreq', ...
            'parameter', ...
            ['Low-pass cutoff (Hz) applied to the reference channels before ' ...
             'the regression. LinearRegression algorithm on Y-X-T ' ...
             'fluorescence only: must be 0 with Ratiometric or with ' ...
             'event-split (Y-X-T-E) data. Must be below half the ' ...
             'fluorescence frame rate. 0 = no extra filtering.'], ...
            'kind', 'parameter', ...
            'position', 8, ...
            'callType', 'namevalue', ...
            'default', 0, ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            ['Frame rate of the fluorescence data (Hz), injected from the ' ...
             '"data" input (never from a reference channel).'], ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'sourceInput', 'data');

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Hemodynamically corrected fluorescence output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

function refFile = localReferenceFileName(channelTag)
%LOCALREFERENCEFILENAME Map a channel tag or filename to a .dat file name.

tag = lower(char(string(channelTag)));

switch tag
    case 'red'
        refFile = 'red.dat';
    case {'amber', 'yellow'}
        refFile = 'yellow.dat';
    case 'green'
        refFile = 'green.dat';
    otherwise
        [~, name, ext] = fileparts(char(string(channelTag)));
        if isempty(ext)
            ext = '.dat';
        end
        refFile = [name, ext];
end
end

function refFile = localResolveReferenceFile(SaveFolder, channelTag)
%LOCALRESOLVEREFERENCEFILE Resolve channel tag or filename in SaveFolder.

refFile = localReferenceFileName(channelTag);

assert(isfile(fullfile(SaveFolder, refFile)), ...
    'Umitoolbox:HemoCorrection:FileNotFound', ...
    'Reference channel file "%s" not found.', refFile);
end

% =========================================================================
% LinearRegression (Y-X-T fluorescence)
% =========================================================================
function outData = localLinearRegression(data, SaveFolder, channelList, ...
        lowPassFreq, frameRateHz, defaultOutput)
%LOCALLINEARREGRESSION Pixel-wise regression onto the reference channels.
%
% Numeric DATA is corrected in RAM (FRAMERATEHZ is its frame rate). A .dat
% DATA is streamed and written to DEFAULTOUTPUT in SaveFolder; its frame
% rate comes from the .dat header.

isFile = ischar(data) || (isstring(data) && isscalar(data));
if isFile
    fileData = localResolveFileInSaveFolder(SaveFolder, data);
    fMetaData = localNormalizeDatMeta(loadMetaData(fileData));
else
    fileData = data;
    fMetaData = localResolveNumericYXTMetadata(data, frameRateHz);
end

assert(lowPassFreq < fMetaData.Freq/2 || lowPassFreq == 0, ...
    'Umitoolbox:HemoCorrection:invalidLowPassFreq', ...
    'LowPassFreq must be smaller than the fluorescence Nyquist frequency.');

assert(~isempty(channelList), ...
    'Umitoolbox:HemoCorrection:InvalidInput', ...
    'No reference channel is selected. Enable Red, Green, or Amber, or set Other.');

resolvedFiles = cell(1, numel(channelList));
for iChan = 1:numel(channelList)
    resolvedFiles{iChan} = localReferenceFileName(channelList{iChan});
end

resolvedFiles = fullfile(SaveFolder, resolvedFiles);
for iFile = 1:numel(resolvedFiles)
    assert(isfile(resolvedFiles{iFile}), ...
        'Umitoolbox:HemoCorrection:missingReferenceChannel', ...
        'Hemodynamic reference channel was not found in SaveFolder: %s', resolvedFiles{iFile});
end

if isFile
    % Write through a fixed-name scratch file, then move it onto the
    % declared output. Renaming the output when it already exists would
    % make every pipeline re-run write to a different file and leave the
    % stale original in place.
    outPath = fullfile(SaveFolder, defaultOutput);
    [~, baseName, ext] = fileparts(defaultOutput);
    tmpPath = fullfile(SaveFolder, [baseName '_writing' ext]);

    localLinearRegressionLowRAM(tmpPath, baseName, fileData, fMetaData, resolvedFiles, lowPassFreq);

    [moveOk, moveMsg] = movefile(tmpPath, outPath, 'f');
    assert(moveOk, 'Umitoolbox:HemoCorrection:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', tmpPath, outPath, moveMsg);

    outData = defaultOutput;
else
    outData = localLinearRegressionStandard(fileData, fMetaData, resolvedFiles, lowPassFreq);
end
end

function fData = localLinearRegressionStandard(fData, fMetaData, colorList, LPcutoffFreq)
%LOCALLINEARREGRESSIONSTANDARD In-memory hemodynamic correction.

Ny = fMetaData.datSize(1);
Nx = fMetaData.datSize(2);
Nt = fMetaData.datLength;
Np = Nx * Ny;

numChannels = numel(colorList);
cMetaData = cell(1, numChannels);
maxNt = Nt;
for kk = 1:numChannels
    cMetaData{kk} = localNormalizeDatMeta(loadMetaData(colorList{kk}));
    localValidateSpatialMatch(cMetaData{kk}, fMetaData, colorList{kk});
    maxNt = max(maxNt, cMetaData{kk}.datLength);
end

localAssertReferenceDurations(fMetaData, cMetaData, colorList);
localWarnSkippedLowPass(cMetaData, colorList, LPcutoffFreq);

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

% Restore the caller's warning state, also on error.
warnState = warning('off', 'MATLAB:rankDeficientMatrix');
c_warn = onCleanup(@() warning(warnState));

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

        tmp = localResampleReferenceToFluoTimeline(tmp, cMetaData{kk}, fMetaData, LPcutoffFreq);
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

function outFilename = localLinearRegressionLowRAM(outFilename, outBaseName, fluoFile, fMetaData, colorList, LPcutoffFreq)
%LOCALLINEARREGRESSIONLOWRAM Disk-streamed hemodynamic correction.
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
    cMetaData{k} = localNormalizeDatMeta(loadMetaData(colorList{k}));
    refSlab{k} = spatialSlabIO('open', colorList{k}, 'Info', cMetaData{k});
    c_r{k} = onCleanup(@() spatialSlabIO('close', refSlab{k}));
    localValidateSpatialMatch(cMetaData{k}, fMetaData, colorList{k});
    maxNt = max(maxNt, cMetaData{k}.datLength);
end

localAssertReferenceDurations(fMetaData, cMetaData, colorList);
localWarnSkippedLowPass(cMetaData, colorList, LPcutoffFreq);

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

% Restore the caller's warning state, also on error.
warnState = warning('off', 'MATLAB:rankDeficientMatrix');
c_warn = onCleanup(@() warning(warnState));

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
        tmp = localResampleReferenceToFluoTimeline(tmp, cMetaData{kk}, fMetaData, LPcutoffFreq);
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

function localAssertReferenceDurations(fMetaData, cMetaData, colorList)
%LOCALASSERTREFERENCEDURATIONS Every reference spans the fluorescence duration.
%
% Interpolation is only valid for channels covering the same acquisition
% interval, so this runs before any temporal resampling.

durationTolSec = 1e-3;
fluoDurationSec = double(fMetaData.datLength) / double(fMetaData.Freq);
assert(isfinite(fluoDurationSec) && fluoDurationSec > 0, ...
    'Umitoolbox:HemoCorrection:InvalidChannelDuration', ...
    'Fluorescence channel has invalid duration metadata: Length=%g, FrameRateHz=%g.', ...
    double(fMetaData.datLength), double(fMetaData.Freq));
for kk = 1:numel(cMetaData)
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
end

function localWarnSkippedLowPass(cMetaData, colorList, LPcutoffFreq)
%LOCALWARNSKIPPEDLOWPASS Warn when LowPassFreq cannot filter a reference.
%
% A cutoff at or above a reference's Nyquist frequency cannot be applied to
% it (localResampleReferenceToFluoTimeline skips it), which only happens
% for a reference slower than the fluorescence.

if LPcutoffFreq <= 0
    return
end
for kk = 1:numel(cMetaData)
    refNyquist = double(cMetaData{kk}.Freq) / 2;
    if LPcutoffFreq >= refNyquist
        [~, refName, refExt] = fileparts(colorList{kk});
        warning('Umitoolbox:HemoCorrection:LowPassFreqSkipped', ...
            ['Low-pass filtering of reference channel "%s" is skipped: ' ...
             'LowPassFreq (%g Hz) is not below its Nyquist frequency ' ...
             '(%g Hz, half its %g Hz frame rate), so the channel cannot ' ...
             'contain frequencies above the cutoff. Use a LowPassFreq ' ...
             'below %g Hz to filter it.'], ...
            [refName refExt], LPcutoffFreq, refNyquist, ...
            double(cMetaData{kk}.Freq), refNyquist);
    end
end
end

% =========================================================================
% Ratiometric (Y-X-T fluorescence)
% =========================================================================
function outData = localRatiometricStandard(SaveFolder, data, fluoMeta, refFile)
%LOCALRATIOMETRICSTANDARD Ratiometric correction in standard mode.

refPath = fullfile(SaveFolder, refFile);
refMeta = localNormalizeDatMeta(loadMetaData(refPath));
refData = loadData(refPath);

localValidateSpatialMatch(refMeta, fluoMeta, refFile);
refData = localResampleReferenceToFluoTimeline(refData, refMeta, fluoMeta, 0);

assert(isequal(size(data), size(refData)), ...
    'Umitoolbox:HemoCorrection:InvalidInput', ...
    'Reference and fluorescence channels must have the same size after temporal alignment.');

datSize = size(data);

% Fluorescence normalization
fData = reshape(data, [], datSize(end));
mData = mean(fData, 2, 'omitnan');
fData = (fData - mData) ./ mData;

% Reference normalization
ref2D = reshape(refData, [], datSize(end));
mRef = mean(ref2D, 2, 'omitnan');
ref2D = (ref2D - mRef) ./ mRef;

% Ratiometric correction
fData = fData - ref2D;

% Restore fluorescence mean
fData = (fData .* mData) + mData;

outData = reshape(fData, datSize);
end

function outFile = localRatiometricLowRAM(SaveFolder, ~, fluoMeta, refFile, defaultOutput)
%LOCALRATIOMETRICLOWRAM Ratiometric correction in low-RAM mode.

% Write through a scratch file, then move it onto the declared output.
% Renaming the output when it already exists would make every pipeline
% re-run write to a different file and leave the stale original in place.
outFile = fullfile(SaveFolder, defaultOutput);
[~, outBaseName, outExt] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [outBaseName '_writing' outExt]);

refPath = fullfile(SaveFolder, refFile);
refMeta = localNormalizeDatMeta(loadMetaData(refPath));
localValidateSpatialMatch(refMeta, fluoMeta, refFile);

Nx = fluoMeta.datSize(2);
Nt = fluoMeta.datLength;

% Read each input through the file its metadata were resolved from.
slabFluo = spatialSlabIO('open', fluoMeta.filePath, 'Info', fluoMeta);
cFluo = onCleanup(@() spatialSlabIO('close', slabFluo));

slabRef = spatialSlabIO('open', refMeta.filePath, 'Info', refMeta);
cRef = onCleanup(@() spatialSlabIO('close', slabRef));

slabOut = spatialSlabIO('create', tmpFile, datHeaderFromInfo(fluoMeta, outBaseName));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

% Chunk along X to control RAM.
dataBytes = prod([fluoMeta.datSize, max(fluoMeta.datLength, refMeta.datLength), getByteSize(fluoMeta.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes, 3, .1);
chunkX = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkX);

for ii = 1:nChunks
    xStart = (ii - 1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    % Read slabs using each file's own temporal length.
    fSlab = spatialSlabIO('read', slabFluo, xIdx);
    rSlab = spatialSlabIO('read', slabRef, xIdx);
    rSlab = localResampleReferenceToFluoTimeline(rSlab, refMeta, fluoMeta, 0);

    assert(all(size(fSlab) == size(rSlab)), ...
        'Umitoolbox:HemoCorrection:InvalidInput', ...
        'Reference and fluorescence slabs must have the same size after temporal alignment.');

    slabSz = size(fSlab);

    % Reshape to [pixels x time].
    fSlab = reshape(fSlab, [], Nt);
    rSlab = reshape(rSlab, [], Nt);

    % Normalize each pixel trace.
    mFluo = mean(fSlab, 2, 'omitnan');
    fSlab = (fSlab - mFluo) ./ mFluo;

    mRef = mean(rSlab, 2, 'omitnan');
    rSlab = (rSlab - mRef) ./ mRef;

    % Ratiometric correction.
    fSlab = fSlab - rSlab;

    % Restore fluorescence mean.
    fSlab = (fSlab .* mFluo) + mFluo;

    % Write corrected slab.
    fSlab = reshape(fSlab, slabSz);
    spatialSlabIO('write', slabOut, xIdx, fSlab);
end

% Close every handle before the move: on Windows an open handle blocks it,
% and the input may be the file the declared output overwrites.
spatialSlabIO('finalize', slabOut);
spatialSlabIO('close', slabFluo);
spatialSlabIO('close', slabRef);

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:HemoCorrection:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);

outFile = defaultOutput;
end

% =========================================================================
% Shared helpers
% =========================================================================
function refData = localResampleReferenceToFluoTimeline(refData, refMeta, fluoMeta, LPcutoffFreq)
%LOCALRESAMPLEREFERENCETOFLUOTIMELINE Match reference data to fluorescence T.
%
% LPCUTOFFFREQ > 0 low-pass filters the reference at that cutoff (Hz).
% With 0, a reference faster (and longer) than the fluorescence is
% low-pass filtered at 0.45 x the fluorescence rate. A cutoff at or above
% the reference Nyquist frequency is not applied.

NtRef = refMeta.datLength;
NtFluo = fluoMeta.datLength;
freqRef = refMeta.Freq;
freqFluo = fluoMeta.Freq;

if size(refData,3) ~= NtRef
    error('Umitoolbox:HemoCorrection:InvalidReferenceLength', ...
        'Reference data length does not match its metadata.');
end

% Refuse to resample channels that do not cover the same recording span.
% Interpolation is valid for different sampling rates, not for cropped or
% truncated channels.
refDurationSec = double(NtRef) / double(freqRef);
fluoDurationSec = double(NtFluo) / double(freqFluo);
durationTolSec = 1e-3;
assert(abs(refDurationSec - fluoDurationSec) <= durationTolSec, ...
    'Umitoolbox:HemoCorrection:DurationMismatch', ...
    ['Reference channel does not span the same recording duration as the ' ...
     'fluorescence channel. Reference: Length=%d, FrameRateHz=%g, ' ...
     'Duration=%0.6f s. Fluorescence: Length=%d, FrameRateHz=%g, ' ...
     'Duration=%0.6f s.'], ...
    NtRef, freqRef, refDurationSec, NtFluo, freqFluo, fluoDurationSec);

cutoffFreq = 0;
if LPcutoffFreq > 0
    cutoffFreq = LPcutoffFreq;
elseif freqRef > freqFluo && NtRef > NtFluo
    cutoffFreq = 0.45 * freqFluo;
end

if cutoffFreq > 0 && cutoffFreq < freqRef/2
    sz = size(refData);
    refData = reshape(refData, [], sz(3));
    f = fdesign.lowpass('N,F3dB', 4, cutoffFreq, freqRef);
    lpass = design(f, 'butter');
    refData = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(refData')))' ;
    refData = reshape(refData, sz);
end

if NtRef ~= NtFluo
    sz = size(refData);
    xRef = linspace(0, 1, NtRef);
    xFluo = linspace(0, 1, NtFluo);
    refData = reshape(refData, [], NtRef);
    refData = interp1(xRef, single(refData)', xFluo, 'linear', 'extrap')';
    refData = reshape(single(refData), sz(1), sz(2), NtFluo);
else
    refData = single(refData);
end
end

function meta = localResolveNumericYXTMetadata(data, frameRateHz)
%LOCALRESOLVENUMERICYXTMETADATA Build file-like metadata for numeric YXT input.
%
% FRAMERATEHZ is the data's own frame rate (resolveDataInfoValue).

% Y and X come from the array itself: AcqInfos.mat Height/Width is the raw
% acquisition size and no longer matches aligned data. The reference
% channel is checked against this size (localValidateSpatialMatch).
height = double(size(data, 1));
width = double(size(data, 2));

nT = double(size(data, 3));

meta = struct();
meta.datSize = [height, width];
meta.Height = height;
meta.Width = width;
meta.datLength = nT;
meta.Length = nT;
meta.Freq = double(frameRateHz);
meta.FrameRateHz = double(frameRateHz);
meta.Datatype = 'single';
meta.dim_names = {'Y','X','T'};
end

function meta = localNormalizeDatMeta(meta)
%LOCALNORMALIZEDATMETA Add this function's internal size fields to loadMetaData output.
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

function localValidateSpatialMatch(refMeta, fluoMeta, refFile)
%LOCALVALIDATESPATIALMATCH Validate reference/fluorescence spatial dimensions.

assert(isequal(double(refMeta.datSize(1:2)), double(fluoMeta.datSize(1:2))), ...
    'Umitoolbox:HemoCorrection:SpatialMismatch', ...
    'Reference channel "%s" has incompatible spatial dimensions.', refFile);
end

function [filePath, fileName] = localResolveFileInSaveFolder(SaveFolder, fileInput)
%LOCALRESOLVEFILEINSAVEFOLDER Resolve a filename or full path to a .dat file.

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

function dataOut = iDataFileOrArray(SaveFolder, data)
%IDATAFILEORARRAY The fluorescence .dat path for file input, else the array.
dataOut = data;
if ischar(data) || (isstring(data) && isscalar(data))
    dataOut = localResolveFileInSaveFolder(SaveFolder, data);
end
end

function hasE = iClassifyInput(data, SaveFolder)
%ICLASSIFYINPUT Validate the fluorescence input; true when it has an E axis.
%
% Accepted: a Y-X-T or Y-X-T-E numeric array, or a .dat file with those axes.

if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty'}, ...
        'HemoCorrection', 'data');
    assert(ndims(data) == 3 || ndims(data) == 4, ...
        'Umitoolbox:HemoCorrection:InvalidInput', ...
        'Numeric input must be a Y x X x T or Y x X x T x E array.');
    hasE = ndims(data) == 4;
elseif ischar(data) || (isstring(data) && isscalar(data))
    filePath = localResolveFileInSaveFolder(SaveFolder, data);
    [~, ~, ext] = fileparts(filePath);
    assert(strcmpi(ext, '.dat'), ...
        'Umitoolbox:HemoCorrection:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    assert(isfile(filePath), ...
        'Umitoolbox:HemoCorrection:FileNotFound', ...
        'Fluorescence file "%s" was not found.', filePath);
    info = loadMetaData(filePath);
    assertDatLayout(info, {{'Y','X','T'}, {'Y','X','T','E'}}, 'HemoCorrection');
    hasE = any(strcmp(cellstr(string(info.dimNames)), 'E'));
else
    error('Umitoolbox:HemoCorrection:UnsupportedInputType', ...
        ['Input "data" must be a YXT or YXTE array or a .dat filename. ' ...
         'UMT structs and .umt files are not supported.']);
end
end

function localAssertSingleReference(channelList)
%LOCALASSERTSINGLEREFERENCE Ratiometric correction uses exactly one reference.

assert(isscalar(channelList), ...
    'Umitoolbox:HemoCorrection:InvalidInput', ...
    ['Ratiometric correction requires exactly one reference channel, ' ...
     'but %d are currently selected (%s). Red, Green, and Amber each ' ...
     'default to true; set all but one to false before using ' ...
     'Algorithm=''Ratiometric''.'], ...
    numel(channelList), strjoin(channelList, ', '));
end

% =========================================================================
% Event-split (Y-X-T-E) fluorescence
% =========================================================================
function outData = localCorrectEventSplit(data, SaveFolder, algorithm, channelList, frameRateArg, defaultOutput)
%LOCALCORRECTEVENTSPLIT Hemodynamic correction of event-split fluorescence.
%
% Every E slice is a trial corrected on its own. The reference channels are
% continuous recordings: each trial's frame times (events.mat, at the
% FLUORESCENCE frame rate) are sampled in every reference by linear
% interpolation, so references of any frame rate land on the fluorescence
% timeline of the trial. A faster reference is first low-pass filtered, as in
% the continuous path. LinearRegression follows the algorithm of
% localLinearRegression, Ratiometric that of localRatiometricStandard, each
% applied to one trial at a time.

isRatiometric = strcmpi(algorithm, 'ratiometric');
assert(~isempty(channelList), ...
    'Umitoolbox:HemoCorrection:InvalidInput', ...
    'No reference channel is selected. Enable Red, Green, or Amber, or set Other.');
if isRatiometric
    localAssertSingleReference(channelList);
end

% ---- Fluorescence ------------------------------------------------------
isFile = ischar(data) || (isstring(data) && isscalar(data));
if isFile
    fluoPath = localResolveFileInSaveFolder(SaveFolder, data);
    fluoInfo = loadMetaData(fluoPath);
    fluoRate = resolveDataInfoValue('frameRateHz', frameRateArg, fluoPath, mfilename);
    Ny = datAxisSize(fluoInfo, 'Y');
    Nx = datAxisSize(fluoInfo, 'X');
    Nt = datAxisSize(fluoInfo, 'T');
    Ne = datAxisSize(fluoInfo, 'E');
    mapInfo = fluoInfo;
else
    fluoRate = resolveDataInfoValue('frameRateHz', frameRateArg, data, mfilename);
    [Ny, Nx, Nt, Ne] = size(data);
    mapInfo = struct('filePath', 'input data', 'dimNames', {{'Y','X','T','E'}}, ...
        'dimSizes', [Ny, Nx, Nt, Ne]);
end

mapping = resolveDatEventMapping(mapInfo, SaveFolder);
if ~strcmpi(mapping.status, 'matched')
    reason = mapping.message;
    if isempty(reason)
        reason = ['Its E axis holds one slice per condition (aggregated data), ' ...
                  'not one per event instance.'];
    end
    error('Umitoolbox:HemoCorrection:EventsNotMatched', ...
        ['Event-split fluorescence needs one E slice per event instance of the ' ...
         'events.mat in "%s". %s'], SaveFolder, reason);
end

% ---- Reference channels (continuous Y-X-T files) ---------------------------
nRef = numel(channelList);
refPaths = cell(1, nRef);
refMeta = cell(1, nRef);
refRate = zeros(1, nRef);
refNt = zeros(1, nRef);
fluoSizeMeta = struct('datSize', [Ny, Nx]);
for k = 1:nRef
    refFile = localResolveReferenceFile(SaveFolder, channelList{k});
    refPaths{k} = fullfile(SaveFolder, refFile);
    refMeta{k} = localNormalizeDatMeta(loadMetaData(refPaths{k}));
    assertDatLayout(refMeta{k}, {{'Y','X','T'}}, 'HemoCorrection');
    localValidateSpatialMatch(refMeta{k}, fluoSizeMeta, refFile);
    refRate(k) = refMeta{k}.Freq;
    refNt(k) = refMeta{k}.datLength;
    assert(isfinite(refRate(k)) && refRate(k) > 0, ...
        'Umitoolbox:HemoCorrection:InvalidReferenceRate', ...
        'Reference channel "%s" has no valid frame rate in its header.', refFile);
end

refDurationSec = refNt ./ refRate;
assert(max(refDurationSec) - min(refDurationSec) <= 1e-3, ...
    'Umitoolbox:HemoCorrection:DurationMismatch', ...
    ['Reference channels do not span the same recording duration (%s s). ' ...
     'Event trials are cut from channels that share one acquisition clock.'], ...
    mat2str(refDurationSec, 6));

% ---- Trial times on the fluorescence timeline ------------------------------
% The recording length is inferred from the references (two frames of slack
% so a last trial is not truncated); only the first Nt frames of each trial
% are used, as split_data_by_event cropped every trial to the shortest one.
datLen = ceil(refDurationSec(1) * fluoRate) + 2;
frMat = EventsManager(SaveFolder).getFrameMatrix(datLen, '', [], ...
    'FrameRateHz', fluoRate, 'IncludeIgnored', true);
assert(size(frMat, 1) == Ne && size(frMat, 2) >= Nt, ...
    'Umitoolbox:HemoCorrection:EventsNotMatched', ...
    ['The trials of events.mat (%d trials, up to %d frames at %g Hz) do not fit ' ...
     'the event-split fluorescence (%d trials of %d frames).'], ...
    size(frMat, 1), size(frMat, 2), fluoRate, Ne, Nt);
frMat = frMat(:, 1:Nt);
assert(~any(isnan(frMat(:))), ...
    'Umitoolbox:HemoCorrection:EventsNotMatched', ...
    'Some trial frames of events.mat fall outside the recording.');
tFluo = (frMat - 1) / fluoRate;

% Per reference: where each trial frame falls in the reference (linear
% interpolation between frames k0 and k1 with weight w), and the anti-alias
% filter when the reference is faster than the fluorescence.
tables = repmat(struct('k0', [], 'k1', [], 'w', [], 'lowpass', []), 1, nRef);
for k = 1:nRef
    pos = tFluo * refRate(k) + 1;
    assert(all(pos(:) >= 0) && all(pos(:) <= refNt(k) + 1), ...
        'Umitoolbox:HemoCorrection:EventsOutsideReference', ...
        ['Trial times of events.mat fall more than one frame outside reference ' ...
         'channel "%s".'], channelList{k});
    pos = min(max(pos, 1), refNt(k));
    k0 = floor(pos);
    tables(k).k0 = k0;
    tables(k).k1 = min(k0 + 1, refNt(k));
    tables(k).w = single(pos - k0);
    tables(k).lowpass = iAntiAliasFilter(refRate(k), fluoRate);
end

% ---- Slabs ---------------------------------------------------------------
spatSigma = 1;
pad = 0;
if ~isRatiometric
    pad = ceil(3 * spatSigma);
end

bytesPerX = Ny * 4 * (3 * max(refNt) + (4 + nRef) * Nt * Ne);
nChunks = calculateMaxChunkSize(bytesPerX * Nx, 1, .1);
chunkX = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkX);

outFile = fullfile(SaveFolder, defaultOutput);
if isFile
    [~, outStem, outExt] = fileparts(defaultOutput);
    tmpFile = fullfile(SaveFolder, [outStem '_writing' outExt]);
    slabFluo = spatialSlabIO('open', fluoPath, 'Info', fluoInfo);
    cFluo = onCleanup(@() spatialSlabIO('close', slabFluo));
    slabOut = spatialSlabIO('create', tmpFile, ...
        datHeaderFromInfo(fluoInfo, outStem, 'dataClass', 'single'));
    cOut = onCleanup(@() spatialSlabIO('close', slabOut));
else
    outData = zeros(Ny, Nx, Nt, Ne, 'single');
end
slabRef = cell(1, nRef);
for k = 1:nRef
    slabRef{k} = spatialSlabIO('open', refPaths{k}, 'Info', refMeta{k});
end
cRef = onCleanup(@() cellfun(@(h) spatialSlabIO('close', h), slabRef));

for c = 1:nChunks
    xStart = (c - 1) * chunkX + 1;
    xEnd = min(xStart + chunkX - 1, Nx);
    xIdx = xStart:xEnd;
    padStart = min(pad, xStart - 1);
    padStop = min(pad, Nx - xEnd);
    xIdxPad = (xStart - padStart):(xEnd + padStop);
    nX = numel(xIdx);
    Np = Ny * nX;
    fprintf('Hemodynamic correction (event-split): chunk %i/%i\n', c, nChunks);

    if isFile
        fSlab = single(spatialSlabIO('read', slabFluo, xIdx));
    else
        fSlab = single(data(:, xIdx, :, :));
    end
    fSlab = reshape(fSlab, Np, Nt, Ne);

    refData = zeros(nRef, Np, Nt, Ne, 'single');
    for k = 1:nRef
        raw = single(spatialSlabIO('read', slabRef{k}, xIdxPad));
        trials = iSampleReferenceTrials(raw, tables(k));      % Ny x nPad x Nt x Ne
        if ~isRatiometric
            nPad = size(trials, 2);
            trials = imgaussfilt(reshape(trials, Ny, nPad, Nt * Ne), spatSigma, ...
                'Padding', 'symmetric');
            trials = reshape(trials, Ny, nPad, Nt, Ne);
            trials = trials(:, padStart + 1:end - padStop, :, :);
        end
        trials = reshape(trials, Np, Nt, Ne);

        % Normalize each trial of each pixel: (x - mean) / mean.
        if isRatiometric
            m = mean(trials, 2, 'omitnan');
        else
            m = mean(trials, 2);
        end
        refData(k, :, :, :) = reshape((trials - m) ./ m, 1, Np, Nt, Ne);
    end

    if isRatiometric
        mFluo = mean(fSlab, 2, 'omitnan');
        fSlab = (((fSlab - mFluo) ./ mFluo) - reshape(refData, Np, Nt, Ne)) .* mFluo + mFluo;
    else
        mFluo = mean(fSlab, 2);
        fSlab = iRegressTrials((fSlab - mFluo) ./ mFluo, refData) .* mFluo + mFluo;
    end

    fSlab = reshape(fSlab, Ny, nX, Nt, Ne);
    if isFile
        spatialSlabIO('write', slabOut, xIdx, fSlab);
    else
        outData(:, xIdx, :, :) = fSlab;
    end
end

if isFile
    % Close every handle before the move: on Windows an open handle blocks it,
    % and the input may be the file the declared output overwrites.
    spatialSlabIO('finalize', slabOut);
    spatialSlabIO('close', slabFluo);
    cellfun(@(h) spatialSlabIO('close', h), slabRef);

    [moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
    assert(moveOk, 'Umitoolbox:HemoCorrection:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
    outData = defaultOutput;
end
end

function lp = iAntiAliasFilter(refRate, fluoRate)
%IANTIALIASFILTER Low-pass filter for a reference faster than the fluorescence.
%
% Same rule as the continuous path: 0.45 x the fluorescence rate. Empty when
% the reference is not faster (or the cutoff is not usable).

lp = [];
if refRate > fluoRate
    cutoff = 0.45 * fluoRate;
    if cutoff > 0 && cutoff < refRate / 2
        f = fdesign.lowpass('N,F3dB', 4, cutoff, refRate);
        d = design(f, 'butter');
        lp = struct('sos', d.sosMatrix, 'scale', d.ScaleValues);
    end
end
end

function trials = iSampleReferenceTrials(raw, tbl)
%ISAMPLEREFERENCETRIALS Sample a continuous reference slab at every trial's frame times.
%
% RAW is Y x X x T(reference); the result is Y x X x Nt x Ne, interpolating
% linearly between reference frames K0 and K1 with weight W (up- or
% down-sampling onto the fluorescence timeline of each trial).

if ~isempty(tbl.lowpass)
    sz = size(raw);
    flat = double(reshape(raw, [], sz(3)));
    raw = reshape(single(filtfilt(tbl.lowpass.sos, tbl.lowpass.scale, flat.')).', sz);
end

[Ne, Nt] = size(tbl.k0);
trials = zeros(size(raw, 1), size(raw, 2), Nt, Ne, 'single');
for e = 1:Ne
    a = raw(:, :, tbl.k0(e, :));
    b = raw(:, :, tbl.k1(e, :));
    trials(:, :, :, e) = a + (b - a) .* reshape(tbl.w(e, :), 1, 1, Nt);
end
end

function fNorm = iRegressTrials(fNorm, refData)
%IREGRESSTRIALS Per-trial, per-pixel regression of the normalized references.
%
% FNORM is Np x Nt x Ne (normalized fluorescence); REFDATA is nRef x Np x Nt x Ne.
% The design matrix of localLinearRegression: constant, linear drift, and the
% reference traces. The fit is subtracted from every trial.

[Np, Nt, Ne] = size(fNorm);
nRef = size(refData, 1);
baseTerms = single([ones(Nt, 1), linspace(0, 1, Nt).']);

warnState = warning('off', 'MATLAB:rankDeficientMatrix');
cleanupWarn = onCleanup(@() warning(warnState));

for e = 1:Ne
    for p = 1:Np
        X = [baseTerms, reshape(refData(:, p, :, e), nRef, Nt).'];
        y = reshape(fNorm(p, :, e), Nt, 1);
        fNorm(p, :, e) = (y - X * (X \ y)).';
    end
end
end
