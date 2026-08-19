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
%       - numeric fluorescence array [Y, X, T]
%       - fluorescence .dat filename
%
%   Inputs:
%       data       - Numeric YXT fluorescence array or .dat filename.
%       SaveFolder - Folder containing AcqInfos.mat and reference files.
%
%   Name-Value parameters:
%       'Algorithm' - 'LinearRegression' or 'Ratiometric'
%       'Red'       - logical scalar
%       'Green'     - logical scalar
%       'Amber'     - logical scalar
%       'Other'     - custom channel filename
%
%   Output:
%       - standard mode: corrected fluorescence array
%       - low-RAM mode : corrected fluorescence filename

% Default output for pipeline management:
default_Output = 'hemoCorr_fluo.dat'; 

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'run_HemoCorrection';
addRequired(p, 'data', @(x) (isnumeric(x) && ndims(x) == 3) || ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addParameter(p, 'Algorithm', 'LinearRegression', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Red', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Green', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Amber', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Other', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
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

acqFile = fullfile(SaveFolder, 'AcqInfos.mat');
assert(isfile(acqFile), ...
    'Umitoolbox:run_HemoCorrection:MissingAcqInfos', ...
    'AcqInfos.mat was not found in "%s".', SaveFolder);

md = load(acqFile, 'AcqInfoStream');
assert(isfield(md, 'AcqInfoStream') && isstruct(md.AcqInfoStream), ...
    'Umitoolbox:run_HemoCorrection:InvalidAcqInfos', ...
    'AcqInfos.mat does not contain a valid AcqInfoStream structure.');
acq = md.AcqInfoStream;

% Build selected channel list.
channelList = {};
if useRed
    channelList{end+1} = 'red'; 
end
if useGreen
    channelList{end+1} = 'green'; 
end
if useAmber
    channelList{end+1} = 'amber'; 
end
if ~isempty(otherChan)
    channelList{end+1} = otherChan; %#ok<AGROW>
end

fprintf('Performing hemodynamic correction in fluo channel using %s algorithm...\n', ...
    algorithm);

switch lower(algorithm)
    case 'linearregression'
        outData = HemoCorrection(data, SaveFolder, 'ChannelList', channelList);

    case 'ratiometric'
        assert(numel(channelList) == 1, ...
            'Umitoolbox:run_HemoCorrection:InvalidInput', ...
            'Ratiometric correction requires exactly one reference channel.');

        refFile = localResolveReferenceFile(SaveFolder, channelList{1});
        fprintf('Using channel "%s" in hemodynamic correction...\n', refFile);

        if ischar(data) || (isstring(data) && isscalar(data))
            [fluoPath, fluoName] = localResolveFileInSaveFolder(SaveFolder, data);
            fluoMeta = localNormalizeDatMeta(loadMetaData(fluoPath));
            outData = localRatiometricLowRAM(SaveFolder, fluoName, fluoMeta, refFile, default_Output);
        else
            fluoMeta = localResolveNumericYXTMetadata(data, acq);
            outData = localRatiometricStandard(SaveFolder, data, fluoMeta, refFile);
        end

    otherwise
        error('Umitoolbox:run_HemoCorrection:InvalidInput', ...
            'Unknown correction algorithm "%s".', algorithm);
end

fprintf('Finished hemodynamic correction.\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'run_HemoCorrection', ...
            'Apply wrapper-based hemodynamic correction to a fluorescence channel.');

        info = PipelineManager.addInput(info, ...
            'data', ...
            'ImageTimeSeries', ...
            'Fluorescence image time series input.', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing AcqInfos.mat and reference channels.', ...
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
            'allowed', {false,true}, ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Green', ...
            'parameter', ...
            'Use green.dat as reference.', ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'namevalue', ...
            'default', true, ...
            'allowed', {false,true}, ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Amber', ...
            'parameter', ...
            'Use yellow.dat as reference.', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', true, ...
            'allowed', {false,true}, ...
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

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            'ImageTimeSeries', ...
            'data', ...
            'Hemodynamically corrected fluorescence output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
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

function outFile = localRatiometricLowRAM(SaveFolder, fluoFile, fluoMeta, refFile, defaultOutput)
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

Ny = fluoMeta.datSize(1);
Nx = fluoMeta.datSize(2);
Nt = fluoMeta.datLength;

fidFluo = fopen(fullfile(SaveFolder, fluoFile), 'r');
assert(fidFluo ~= -1, ...
    'Umitoolbox:run_HemoCorrection:FileOpenError', ...
    'Could not open fluorescence file "%s".', fluoFile);
cFluo = onCleanup(@() safeFclose(fidFluo)); 

fidRef = fopen(refPath, 'r');
assert(fidRef ~= -1, ...
    'Umitoolbox:run_HemoCorrection:FileOpenError', ...
    'Could not open reference file "%s".', refFile);
cRef = onCleanup(@() safeFclose(fidRef)); 

preallocateDatFile(tmpFile, [Ny, Nx, Nt], fluoMeta.Datatype);

fidOut = fopen(tmpFile, 'r+');
assert(fidOut ~= -1, ...
    'Umitoolbox:run_HemoCorrection:FileOpenError', ...
    'Could not create output file "%s".', tmpFile);
cOut = onCleanup(@() safeFclose(fidOut)); 

% Chunk along X to control RAM.
dataBytes = prod([fluoMeta.datSize, max(fluoMeta.datLength, refMeta.datLength), getByteSize(fluoMeta.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes, 3, .1);
chunkX = ceil(Nx / nChunks);

for ii = 1:nChunks
    xStart = (ii - 1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    % Read slabs using each file's own temporal length.
    fSlab = spatialSlabIO('read', fidFluo, Ny, Nx, fluoMeta.datLength, xIdx, fluoMeta.Datatype);
    rSlab = spatialSlabIO('read', fidRef,  Ny, Nx, refMeta.datLength, xIdx, refMeta.Datatype);
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
    spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, fluoMeta.Datatype, fSlab);
end

% Close every handle before the move: on Windows an open handle blocks it,
% and the input may be the file the declared output overwrites.
fclose(fidOut);
fclose(fidFluo);
fclose(fidRef);

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

function meta = localResolveNumericYXTMetadata(data, AcqInfoStream)
%LOCALRESOLVENUMERICYXTMETADATA Build file-like metadata for numeric YXT input.

assert(isfield(AcqInfoStream, 'Height') && isfield(AcqInfoStream, 'Width'), ...
    'Umitoolbox:run_HemoCorrection:InvalidAcqInfos', ...
    'AcqInfoStream must contain Height and Width.');

height = double(AcqInfoStream.Height);
width = double(AcqInfoStream.Width);
assert(isequal([size(data,1), size(data,2)], [height, width]), ...
    'Umitoolbox:run_HemoCorrection:InvalidNumericInput', ...
    'Numeric fluorescence input does not match AcqInfos.mat Height/Width.');

timelineInfo = resolveDatTimeline(size(data,3), AcqInfoStream);

meta = struct();
meta.datSize = [height, width];
meta.Height = height;
meta.Width = width;
meta.datLength = double(timelineInfo.Length);
meta.Length = double(timelineInfo.Length);
meta.Freq = double(timelineInfo.FrameRateHz);
meta.FrameRateHz = double(timelineInfo.FrameRateHz);
meta.Datatype = 'single';
meta.dim_names = {'Y','X','T'};
end

function meta = localNormalizeDatMeta(meta)
%LOCALNORMALIZEDATMETA Ensure loadMetaData output has expected fields.

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
