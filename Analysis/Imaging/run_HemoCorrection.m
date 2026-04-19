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
default_Output = 'hemoCorr_fluo.dat'; %#ok<NASGU>

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

assert(isfield(acq, 'Height') && isfield(acq, 'Width') && isfield(acq, 'Length') && isfield(acq, 'FrameRateHz'), ...
    'Umitoolbox:run_HemoCorrection:InvalidAcqInfos', ...
    'AcqInfoStream must contain Height, Width, Length, and FrameRateHz.');

metaData = struct();
metaData.datSize = [double(acq.Height), double(acq.Width)];
metaData.datLength = double(acq.Length);
metaData.Freq = double(acq.FrameRateHz);
metaData.Datatype = 'single';

% Build selected channel list.
channelList = {};
if useRed
    channelList{end+1} = 'red'; %#ok<AGROW>
end
if useGreen
    channelList{end+1} = 'green'; %#ok<AGROW>
end
if useAmber
    channelList{end+1} = 'amber'; %#ok<AGROW>
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
            outData = localRatiometricLowRAM(SaveFolder, char(string(data)), metaData, refFile, default_Output);
        else
            outData = localRatiometricStandard(SaveFolder, data, refFile);
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
            'Folder containing AcqInfos.mat and reference channels.', ...
            'kind', 'parameter', ...
            'position', 2, ...
            'callType', 'positional', ...
            'default', '', ...
            'dataType', 'char');

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
            {'ImageTimeSeries'}, ...
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

function outData = localRatiometricStandard(SaveFolder, data, refFile)
%LOCALRATIOMETRICSTANDARD Ratiometric correction in standard mode.

refData = loadData(fullfile(SaveFolder, refFile));

assert(all(size(data) == size(refData)), ...
    'Umitoolbox:run_HemoCorrection:InvalidInput', ...
    'Reference and fluorescence channels must have the same size.');

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

function outFile = localRatiometricLowRAM(SaveFolder, fluoFile, metaData, refFile, defaultOutput)
%LOCALRATIOMETRICLOWRAM Ratiometric correction in low-RAM mode.

outFile = fullfile(SaveFolder, defaultOutput);
if isfile(outFile)
    [folderPath, baseName, ext] = fileparts(outFile);
    outFile = fullfile(folderPath, [baseName '_preallocData' ext]);
end

Ny = metaData.datSize(1);
Nx = metaData.datSize(2);
Nt = metaData.datLength;

fidFluo = fopen(fullfile(SaveFolder, fluoFile), 'r');
assert(fidFluo ~= -1, ...
    'Umitoolbox:run_HemoCorrection:FileOpenError', ...
    'Could not open fluorescence file "%s".', fluoFile);
cFluo = onCleanup(@() safeFclose(fidFluo)); %#ok<NASGU>

fidRef = fopen(fullfile(SaveFolder, refFile), 'r');
assert(fidRef ~= -1, ...
    'Umitoolbox:run_HemoCorrection:FileOpenError', ...
    'Could not open reference file "%s".', refFile);
cRef = onCleanup(@() safeFclose(fidRef)); %#ok<NASGU>

preallocateDatFile(outFile, [Ny, Nx, Nt], metaData.Datatype);

fidOut = fopen(outFile, 'r+');
assert(fidOut ~= -1, ...
    'Umitoolbox:run_HemoCorrection:FileOpenError', ...
    'Could not create output file "%s".', outFile);
cOut = onCleanup(@() safeFclose(fidOut)); %#ok<NASGU>

% Chunk along X to control RAM.
dataBytes = prod([metaData.datSize, metaData.datLength, getByteSize(metaData.Datatype)]);
nChunks = calculateMaxChunkSize(dataBytes, 3, .1);
chunkX = ceil(Nx / nChunks);

for ii = 1:nChunks
    xStart = (ii - 1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    % Read slabs.
    fSlab = spatialSlabIO('read', fidFluo, Ny, Nx, Nt, xIdx, metaData.Datatype);
    rSlab = spatialSlabIO('read', fidRef,  Ny, Nx, Nt, xIdx, metaData.Datatype);

    assert(all(size(fSlab) == size(rSlab)), ...
        'Umitoolbox:run_HemoCorrection:InvalidInput', ...
        'Reference and fluorescence slabs must have the same size.');

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
    spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, metaData.Datatype, fSlab);
end

[~, outName, outExt] = fileparts(outFile);
outFile = [outName outExt];
end
