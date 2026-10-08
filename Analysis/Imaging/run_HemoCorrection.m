function outData = run_HemoCorrection(data, SaveFolder, varargin)
%RUN_HEMOCORRECTION Apply hemodynamic correction to a fluorescence channel.
%
%   outData = run_HemoCorrection(data, SaveFolder)
%   outData = run_HemoCorrection(data, SaveFolder, 'Algorithm', 'Ratiometric', ...)
%
%   This wrapper supports two correction approaches:
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
p.FunctionName = 'run_HemoCorrection';
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addParameter(p, 'Algorithm', 'LinearRegression', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Red', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Green', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Amber', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Other', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
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

assert(ismember(lower(algorithm), {'linearregression', 'ratiometric'}), ...
    'Umitoolbox:run_HemoCorrection:InvalidInput', ...
    'Unknown correction algorithm "%s".', algorithm);

% Y-X-T runs through the legacy paths below; Y-X-T-E (event-split) fluorescence
% is corrected trial by trial (localCorrectEventSplit).
hasE = iClassifyInput(data, SaveFolder);

% Build selected channel list.
channelList = {};
if useRed
    channelList{end+1} = 'red'; 
end
if useGreen
    channelList{end+1} = 'green'; 
end
if useAmber
    % Normalize to the canonical name also used by run_HemoCompute
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
            outData = HemoCorrection(data, SaveFolder, 'ChannelList', channelList, ...
                'FrameRateHz', resolveRate());
        else
            outData = HemoCorrection(data, SaveFolder, 'ChannelList', channelList);
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
        error('Umitoolbox:run_HemoCorrection:InvalidInput', ...
            'Unknown correction algorithm "%s".', algorithm);
end

if ischar(outData) || (isstring(outData) && isscalar(outData))
    outData = localNormalizeFileOutput( ...
        SaveFolder, outData, default_Output);
end

fprintf('Finished hemodynamic correction.\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'run_HemoCorrection', ...
            'Apply wrapper-based hemodynamic correction to a fluorescence channel.');

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

function outFile = localNormalizeFileOutput(SaveFolder, producedFile, canonicalFile)
%LOCALNORMALIZEFILEOUTPUT Move a core-produced file to the wrapper contract.

producedFile = char(string(producedFile));
if isfile(producedFile)
    sourcePath = producedFile;
else
    sourcePath = fullfile(SaveFolder, producedFile);
end
destinationPath = fullfile(SaveFolder, canonicalFile);

assert(isfile(sourcePath), ...
    'Umitoolbox:run_HemoCorrection:MissingOutputFile', ...
    'Hemodynamic correction did not write the reported output file: %s', ...
    sourcePath);

if ~strcmpi(string(sourcePath), string(destinationPath))
    [moveOk, moveMsg] = movefile(sourcePath, destinationPath, 'f');
    assert(moveOk, 'Umitoolbox:run_HemoCorrection:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', ...
        sourcePath, destinationPath, moveMsg);
end

outFile = canonicalFile;
end

function refFile = localResolveReferenceFile(SaveFolder, channelTag)
%LOCALRESOLVEREFERENCEFILE Resolve channel tag or filename in SaveFolder.

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

assert(isfile(fullfile(SaveFolder, refFile)), ...
    'Umitoolbox:run_HemoCorrection:FileNotFound', ...
    'Reference channel file "%s" not found.', refFile);
end

function outData = localRatiometricStandard(SaveFolder, data, fluoMeta, refFile)
%LOCALRATIOMETRICSTANDARD Ratiometric correction in standard mode.

refPath = fullfile(SaveFolder, refFile);
refMeta = localNormalizeDatMeta(loadMetaData(refPath));
refData = loadData(refPath);

localValidateSpatialMatch(refMeta, fluoMeta, refFile);
refData = localResampleReferenceToFluoTimeline(refData, refMeta, fluoMeta);

assert(isequal(size(data), size(refData)), ...
    'Umitoolbox:run_HemoCorrection:InvalidInput', ...
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
    rSlab = localResampleReferenceToFluoTimeline(rSlab, refMeta, fluoMeta);

    assert(all(size(fSlab) == size(rSlab)), ...
        'Umitoolbox:run_HemoCorrection:InvalidInput', ...
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
assert(moveOk, 'Umitoolbox:run_HemoCorrection:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);

outFile = defaultOutput;
end

function refData = localResampleReferenceToFluoTimeline(refData, refMeta, fluoMeta)
%LOCALRESAMPLEREFERENCETOFLUOTIMELINE Match reference data to fluorescence T.

NtRef = refMeta.datLength;
NtFluo = fluoMeta.datLength;
freqRef = refMeta.Freq;
freqFluo = fluoMeta.Freq;

if size(refData,3) ~= NtRef
    error('Umitoolbox:run_HemoCorrection:InvalidReferenceLength', ...
        'Reference data length does not match its metadata.');
end

% Refuse to resample channels that do not cover the same recording span.
% Interpolation is valid for different sampling rates, not for cropped or
% truncated channels.
refDurationSec = double(NtRef) / double(freqRef);
fluoDurationSec = double(NtFluo) / double(freqFluo);
durationTolSec = 1e-3;
assert(abs(refDurationSec - fluoDurationSec) <= durationTolSec, ...
    'Umitoolbox:run_HemoCorrection:DurationMismatch', ...
    ['Reference channel does not span the same recording duration as the ' ...
     'fluorescence channel. Reference: Length=%d, FrameRateHz=%g, ' ...
     'Duration=%0.6f s. Fluorescence: Length=%d, FrameRateHz=%g, ' ...
     'Duration=%0.6f s.'], ...
    NtRef, freqRef, refDurationSec, NtFluo, freqFluo, fluoDurationSec);

if freqRef > freqFluo && NtRef > NtFluo
    cutoffFreq = 0.45 * freqFluo;
    if cutoffFreq > 0 && cutoffFreq < freqRef/2
        sz = size(refData);
        refData = reshape(refData, [], sz(3));
        f = fdesign.lowpass('N,F3dB', 4, cutoffFreq, freqRef);
        lpass = design(f, 'butter');
        refData = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(refData')))' ;
        refData = reshape(refData, sz);
    end
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
    'Umitoolbox:run_HemoCorrection:SpatialMismatch', ...
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
        'run_HemoCorrection', 'data');
    assert(ndims(data) == 3 || ndims(data) == 4, ...
        'Umitoolbox:run_HemoCorrection:InvalidInput', ...
        'Numeric input must be a Y x X x T or Y x X x T x E array.');
    hasE = ndims(data) == 4;
elseif ischar(data) || (isstring(data) && isscalar(data))
    filePath = localResolveFileInSaveFolder(SaveFolder, data);
    [~, ~, ext] = fileparts(filePath);
    assert(strcmpi(ext, '.dat'), ...
        'Umitoolbox:run_HemoCorrection:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    assert(isfile(filePath), ...
        'Umitoolbox:run_HemoCorrection:FileNotFound', ...
        'Fluorescence file "%s" was not found.', filePath);
    info = loadMetaData(filePath);
    assertDatLayout(info, {{'Y','X','T'}, {'Y','X','T','E'}}, 'run_HemoCorrection');
    hasE = any(strcmp(cellstr(string(info.dimNames)), 'E'));
else
    error('Umitoolbox:run_HemoCorrection:UnsupportedInputType', ...
        ['Input "data" must be a YXT or YXTE array or a .dat filename. ' ...
         'UMT structs and .umt files are not supported.']);
end
end

function localAssertSingleReference(channelList)
%LOCALASSERTSINGLEREFERENCE Ratiometric correction uses exactly one reference.

assert(isscalar(channelList), ...
    'Umitoolbox:run_HemoCorrection:InvalidInput', ...
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
% the continuous path. LinearRegression follows the algorithm of the
% IOIAnalysis HemoCorrection core, Ratiometric that of localRatiometricStandard,
% each applied to one trial at a time.

isRatiometric = strcmpi(algorithm, 'ratiometric');
assert(~isempty(channelList), ...
    'Umitoolbox:run_HemoCorrection:InvalidInput', ...
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
    error('Umitoolbox:run_HemoCorrection:EventsNotMatched', ...
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
    assertDatLayout(refMeta{k}, {{'Y','X','T'}}, 'run_HemoCorrection');
    localValidateSpatialMatch(refMeta{k}, fluoSizeMeta, refFile);
    refRate(k) = refMeta{k}.Freq;
    refNt(k) = refMeta{k}.datLength;
    assert(isfinite(refRate(k)) && refRate(k) > 0, ...
        'Umitoolbox:run_HemoCorrection:InvalidReferenceRate', ...
        'Reference channel "%s" has no valid frame rate in its header.', refFile);
end

refDurationSec = refNt ./ refRate;
assert(max(refDurationSec) - min(refDurationSec) <= 1e-3, ...
    'Umitoolbox:run_HemoCorrection:DurationMismatch', ...
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
    'Umitoolbox:run_HemoCorrection:EventsNotMatched', ...
    ['The trials of events.mat (%d trials, up to %d frames at %g Hz) do not fit ' ...
     'the event-split fluorescence (%d trials of %d frames).'], ...
    size(frMat, 1), size(frMat, 2), fluoRate, Ne, Nt);
frMat = frMat(:, 1:Nt);
assert(~any(isnan(frMat(:))), ...
    'Umitoolbox:run_HemoCorrection:EventsNotMatched', ...
    'Some trial frames of events.mat fall outside the recording.');
tFluo = (frMat - 1) / fluoRate;

% Per reference: where each trial frame falls in the reference (linear
% interpolation between frames k0 and k1 with weight w), and the anti-alias
% filter when the reference is faster than the fluorescence.
tables = repmat(struct('k0', [], 'k1', [], 'w', [], 'lowpass', []), 1, nRef);
for k = 1:nRef
    pos = tFluo * refRate(k) + 1;
    assert(all(pos(:) >= 0) && all(pos(:) <= refNt(k) + 1), ...
        'Umitoolbox:run_HemoCorrection:EventsOutsideReference', ...
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
    assert(moveOk, 'Umitoolbox:run_HemoCorrection:OutputMoveFailed', ...
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
% The design matrix of the IOIAnalysis core: constant, linear drift, and the
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
