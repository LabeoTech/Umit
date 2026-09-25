function outData = run_BloodFlow(SaveFolder, data, varargin)
%RUN_BLOODFLOW Estimate time-resolved blood flow from raw speckle data.
%
%   outData = run_BloodFlow(SaveFolder, data)
%   outData = run_BloodFlow(SaveFolder, data, 'Name', Value, ...)
%   info    = run_BloodFlow('pipelineInfo')
%
%   This function computes a Laser Speckle Contrast Image (LSCI) with the
%   same contrast algorithm as SPECKLEMAPPING, but without averaging the
%   contrast over time. The input Y x X x T shape is therefore preserved.
%   Blood flow is then estimated from the speckle contrast K as:
%
%       Flow = 1 / (T * K^2)
%
%   where T is the camera exposure time in seconds. The output unit is 1/s.
%
%   Algorithm:
%       1) Each pixel is divided by its temporal mean (removal of static
%          structure, as in SPECKLEMAPPING).
%       2) The local standard deviation of the normalized intensity is
%          computed over either a spatial disk of diameter KernelSize
%          ('Spatial', per frame) or a KernelSize-frame temporal window
%          ('Temporal'), both with symmetric padding at the borders. This
%          is the speckle contrast K, and matches the kernels used by
%          SPECKLEMAPPING (5 x 5 disk / 5 frames at the default size).
%       3) K is converted to flow with 1/(T*K^2). Pixels with K <= 0 or a
%          non-finite K are set to NaN.
%       4) Optional (bNormalize): the finished flow is divided by its own
%          per-pixel temporal mean, as in ANA_SPECKLE. NaN samples are
%          ignored when computing that mean.
%
%   Numerical notes:
%       The normalized intensity is close to 1 while K can be as small as
%       1e-5, so a naive local variance loses most of its significant
%       digits to cancellation (STDFILT uses running sums). K is therefore
%       computed in double precision on the intensity centered on 0
%       (normalized - 1; the std is unchanged by the shift), and the
%       temporal std uses MOVSTD instead of STDFILT. The output is single.
%
%   Inputs:
%       SaveFolder - Folder containing the speckle file and metadata.
%       data       - Either:
%                    1) Numeric Y x X x T speckle array held in RAM
%                       (standard mode). The exposure time is resolved
%                       from "<SpeckleFileName>.dat" in SaveFolder.
%                    2) .dat filename (RAM-safe mode). The file is read
%                       from SaveFolder (or used as-is when it is a full
%                       path) and processed in chunks. The exposure time
%                       is resolved from that same file.
%
%   Name-Value parameters:
%       sType           - Speckle contrast mode.
%                         Allowed: 'Spatial', 'Temporal'
%                         Default: 'Temporal'
%       SpeckleFileName - Basename or filename of the speckle file used to
%                         resolve metadata in standard mode.
%                         Default: 'speckle'
%       ExposureMsec    - Exposure time in milliseconds. 'auto' reads
%                         ExposureSpeckleMsec (or ExposureMsec) from the
%                         file metadata. A positive number overrides it.
%                         Default: 'auto'
%       KernelSize      - Size of the contrast window: a disk inscribed in
%                         a KernelSize x KernelSize square, i.e.
%                         fspecial('disk', (KernelSize-1)/2) > 0
%                         ('Spatial'), or KernelSize frames ('Temporal').
%                         Odd integer >= 3 (STDFILT requires an odd window;
%                         a size of 1 gives K = 0).
%                         Default: 5
%       bNormalize      - Logical scalar. If true, normalize the finished
%                         flow by its own per-pixel temporal mean (output
%                         becomes relative flow, unitless).
%                         Default: false
%
%   Output:
%       outData    - Standard mode: single Y x X x T blood-flow array.
%                    RAM-safe mode: full path to the raw BloodFlow.dat
%                    output in SaveFolder.
%
%   See also SPECKLEMAPPING, RUN_SPECKLEMAPPING, RUN_ANA_SPECKLE, STDFILT.

% Default output for pipeline management.
default_Output = 'BloodFlow.dat';

allowedSType = {'Spatial', 'Temporal'};

% -------------------------------------------------------------------------
% pipelineInfo query
% -------------------------------------------------------------------------
if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) ...
        && strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addRequired(p, 'data', @(x) (isnumeric(x) && isreal(x) && ndims(x) == 3) || ...
    ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'sType', 'Temporal', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    any(strcmpi(char(string(x)), allowedSType)));
addParameter(p, 'SpeckleFileName', 'speckle', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ExposureMsec', 'auto', @localIsValidExposure);
addParameter(p, 'KernelSize', 5, @localIsValidKernelSize);
addParameter(p, 'bNormalize', false, @(x) islogical(x) && isscalar(x));
parse(p, SaveFolder, data, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
sType = lower(char(string(p.Results.sType)));
speckleFileName = char(string(p.Results.SpeckleFileName));
exposureOpt = p.Results.ExposureMsec;
kernelSize = double(p.Results.KernelSize);
bNormalize = p.Results.bNormalize;

% RAM-safe mode is inferred from the type of "data" (array vs. filename),
% as in the other speckle wrappers.
bRAMsafeMode = ischar(data) || isstring(data);

if bRAMsafeMode
    datFile = localResolveDatFile(SaveFolder, char(string(data)));
else
    datFile = localResolveDatFile(SaveFolder, speckleFileName);
end
mdIn = loadMetaData(datFile);

exposureSec = localResolveExposure(mdIn, exposureOpt, datFile) / 1000;

fprintf('Calculating blood flow (%s speckle contrast, kernel size %d)...\n', ...
    sType, kernelSize);

if bRAMsafeMode
    outData = localRunRAMSafe(datFile, mdIn, SaveFolder, default_Output, ...
        sType, kernelSize, exposureSec, bNormalize);
else
    data = double(data);
    K = localSpeckleContrast(data, mean(data, 3, 'omitnan'), sType, kernelSize);
    clear data
    outData = localFlowFromContrast(K, exposureSec);
    clear K
    if bNormalize
        outData = localNormalizeByTemporalMean(outData);
    end
end

fprintf('Finished Blood Flow.\n');

    % =====================================================================
    % Local pipelineInfo factory (nested, shares allowedSType and
    % default_Output with the parent scope instead of redefining them)
    % =====================================================================
    function info = localPipelineInfo()

        info = PipelineManager.createPipelineInfo( ...
            mfilename, ...
            'Estimate time-resolved blood flow (1/(T*K^2)) from speckle contrast.');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing the speckle input and metadata.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'data', ...
            'ImageTimeSeries', ...
            'Numeric YXT speckle array or raw .dat filename.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'sType', ...
            'parameter', ...
            'Speckle contrast mode (spatial disk or temporal window).', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', 'Temporal', ...
            'allowed', allowedSType, ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'SpeckleFileName', ...
            'parameter', ...
            'Basename or filename of the speckle file used for metadata in standard mode.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', 'speckle', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'ExposureMsec', ...
            'parameter', ...
            'Exposure time in ms. "auto" reads it from the file metadata.', ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'namevalue', ...
            'default', 'auto', ...
            'allowed', {'auto', [0 Inf]});

        info = PipelineManager.addInput(info, ...
            'KernelSize', ...
            'parameter', ...
            'Contrast window size (odd, >= 3): disk of diameter N (Spatial) or N frames (Temporal).', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', 5, ...
            'allowed', [3 Inf], ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'bNormalize', ...
            'parameter', ...
            'If true, normalize the flow by its own per-pixel temporal mean.', ...
            'kind', 'parameter', ...
            'position', 7, ...
            'callType', 'namevalue', ...
            'default', false, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            {'ImageTimeSeries', 'ProcessedData'}, ...
            'data', ...
            'Blood-flow output (1/s, or relative if bNormalize): a single Y x X x T image time series.', ...
            default_Output, ...
            1, ...
            'isData', true, ...
            'saveFileName', default_Output);
    end

end

% =========================================================================
% Local functions
% =========================================================================
function outFile = localRunRAMSafe(datFile, mdIn, SaveFolder, default_Output, ...
    sType, kernelSize, exposureSec, bNormalize)
%LOCALRUNRAMSAFE Chunked, file-backed blood-flow computation.

assert(isfield(mdIn, 'Height') && isfield(mdIn, 'Width') && isfield(mdIn, 'Length'), ...
    'Umitoolbox:run_BloodFlow:InvalidMetadata', ...
    'Could not resolve Height/Width/Length from "%s".', datFile);

Ny = double(mdIn.Height);
Nx = double(mdIn.Width);
Nt = double(mdIn.Length);

if isfield(mdIn, 'Datatype') && ~isempty(mdIn.Datatype)
    dataType = char(string(mdIn.Datatype));
else
    dataType = 'single';
end

% Compute through a fixed-name scratch file, then move it onto the declared
% output so re-runs overwrite the same file.
outFile = fullfile(SaveFolder, default_Output);
[~, baseName] = fileparts(default_Output);
scratchFile = fullfile(SaveFolder, [baseName '_compute.dat']);
preallocateDatFile(scratchFile, [Ny, Nx, Nt], 'single');

fidIn = fopen(datFile, 'r');
assert(fidIn ~= -1, 'Umitoolbox:run_BloodFlow:FileOpenFailed', ...
    'Could not open "%s" for reading.', datFile);
cIn = onCleanup(@() safeFclose(fidIn));

fidOut = fopen(scratchFile, 'r+');
assert(fidOut ~= -1, 'Umitoolbox:run_BloodFlow:FileOpenFailed', ...
    'Could not open "%s" for writing.', scratchFile);
cOut = onCleanup(@() safeFclose(fidOut));

inBytesPerFrame = Ny * Nx * getByteSize(dataType);
outBytesPerFrame = Ny * Nx * getByteSize('single');
totalBytes = Ny * Nx * Nt * getByteSize('single');

switch sType
    case 'spatial'
        % Pass 1: temporal mean.
        fprintf('Pass 1/%d - Calculating temporal mean...\n', 2 + bNormalize)
        nChunks = calculateMaxChunkSize(totalBytes, 4, .1);
        chunkT = ceil(Nt / nChunks);
        nChunks = ceil(Nt / chunkT);

        sumData = zeros(Ny, Nx, 'double');
        countData = zeros(Ny, Nx, 'double');
        lastPct = -1;
        fprintf('0%% ');
        for c = 1:nChunks
            [slab, ~] = localReadFrames(fidIn, c, chunkT, Nt, Ny, Nx, ...
                inBytesPerFrame, dataType);
            sumData = sumData + sum(slab, 3, 'omitnan');
            countData = countData + sum(~isnan(slab), 3);
            lastPct = localPrintProgress(c, nChunks, lastPct);
        end
        meanData = sumData ./ countData;
        clear sumData countData slab

        % Pass 2: per-frame spatial contrast and flow. The per-pixel flow
        % mean needed by bNormalize is accumulated from the written values.
        fprintf('\nPass 2/%d - Calculating speckle contrast and flow...\n', 2 + bNormalize)
        nChunks = calculateMaxChunkSize(totalBytes, 20, .1);
        chunkT = ceil(Nt / nChunks);
        nChunks = ceil(Nt / chunkT);

        sumFlow = zeros(Ny, Nx, 'double');
        countFlow = zeros(Ny, Nx, 'double');
        lastPct = -1;
        fprintf('0%% ');
        for c = 1:nChunks
            [slab, tStart] = localReadFrames(fidIn, c, chunkT, Nt, Ny, Nx, ...
                inBytesPerFrame, dataType);
            slab = localSpeckleContrast(slab, meanData, 'spatial', kernelSize);
            slab = localFlowFromContrast(slab, exposureSec);
            if bNormalize
                sumFlow = sumFlow + sum(double(slab), 3, 'omitnan');
                countFlow = countFlow + sum(~isnan(slab), 3);
            end

            fseek(fidOut, (tStart-1) * outBytesPerFrame, 'bof');
            fwrite(fidOut, slab, 'single');
            lastPct = localPrintProgress(c, nChunks, lastPct);
        end

        % Pass 3 (bNormalize only): divide by the per-pixel flow mean.
        if bNormalize
            fprintf('\nPass 3/3 - Normalizing flow by its temporal mean...\n')
            meanFlow = sumFlow ./ countFlow;
            clear sumFlow countFlow

            lastPct = -1;
            fprintf('0%% ');
            for c = 1:nChunks
                [slab, tStart] = localReadFrames(fidOut, c, chunkT, Nt, Ny, Nx, ...
                    outBytesPerFrame, 'single');
                slab = single(slab ./ meanFlow);

                fseek(fidOut, (tStart-1) * outBytesPerFrame, 'bof');
                fwrite(fidOut, slab, 'single');
                lastPct = localPrintProgress(c, nChunks, lastPct);
            end
        end

    case 'temporal'
        % Single pass over X slabs: each slab holds the full time course.
        disp('Calculating speckle contrast and flow...')
        nChunks = calculateMaxChunkSize(totalBytes, 24, .15);
        chunkX = ceil(Nx / nChunks);
        nChunks = ceil(Nx / chunkX);

        lastPct = -1;
        fprintf('0%% ');
        for c = 1:nChunks
            xStart = (c-1) * chunkX + 1;
            xEnd = min(xStart + chunkX - 1, Nx);
            xIdx = xStart:xEnd;

            slab = double(spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, dataType));
            slab = localSpeckleContrast(slab, mean(slab, 3, 'omitnan'), 'temporal', kernelSize);
            slab = localFlowFromContrast(slab, exposureSec);
            if bNormalize
                % The slab holds the full time course of its pixels.
                slab = localNormalizeByTemporalMean(slab);
            end
            spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, 'single', slab);
            lastPct = localPrintProgress(c, nChunks, lastPct);
        end
end
fprintf('\n');

clear cIn cOut % close both files before replacing the declared output

[moveOk, moveMsg] = movefile(scratchFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:run_BloodFlow:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', scratchFile, outFile, moveMsg);

end

function [slab, tStart] = localReadFrames(fid, c, chunkT, Nt, Ny, Nx, bytesPerFrame, dataType)
%LOCALREADFRAMES Read one temporal chunk of frames as double.

tStart = (c-1) * chunkT + 1;
tEnd = min(tStart + chunkT - 1, Nt);
nFrames = tEnd - tStart + 1;

fseek(fid, (tStart-1) * bytesPerFrame, 'bof');
slab = fread(fid, Ny * Nx * nFrames, ['*' dataType]);
slab = reshape(double(slab), Ny, Nx, nFrames);

end

function K = localSpeckleContrast(I, meanI, sType, kernelSize)
%LOCALSPECKLECONTRAST Speckle contrast of Y x X x T intensity I (double).
%
%   I is normalized by its per-pixel temporal mean meanI and centered on 0
%   before the local std is taken, to avoid cancellation (see help text).

I = I ./ meanI - 1;

switch sType
    case 'spatial'
        % Same disk as SPECKLEMAPPING, generalized to any odd diameter.
        K = stdfilt(I, fspecial('disk', (kernelSize - 1) / 2) > 0);
    case 'temporal'
        % Same symmetric-padded window as STDFILT(I, ones(1,1,kernelSize)).
        halfWidth = (kernelSize - 1) / 2;
        K = movstd(padarray(I, [0 0 halfWidth], 'symmetric'), kernelSize, 0, 3, ...
            'Endpoints', 'discard');
end

end

function flow = localFlowFromContrast(K, exposureSec)
%LOCALFLOWFROMCONTRAST Convert speckle contrast to flow as 1/(T*K^2).

K = single(K);
K(~isfinite(K) | K <= 0) = NaN;
flow = 1 ./ (single(exposureSec) .* K.^2);

end

function flow = localNormalizeByTemporalMean(flow)
%LOCALNORMALIZEBYTEMPORALMEAN Divide flow by its per-pixel temporal mean.

flow = double(flow);
flow = single(flow ./ mean(flow, 3, 'omitnan'));

end

function exposureMsec = localResolveExposure(mdIn, exposureOpt, datFile)
%LOCALRESOLVEEXPOSURE Return the exposure time in milliseconds.

exposureMsec = localNumericExposure(exposureOpt);
if ~isempty(exposureMsec)
    return
end
if isfield(mdIn, 'ExposureSpeckleMsec') && ~isempty(mdIn.ExposureSpeckleMsec)
    exposureMsec = double(mdIn.ExposureSpeckleMsec);
elseif isfield(mdIn, 'ExposureMsec') && ~isempty(mdIn.ExposureMsec)
    exposureMsec = double(mdIn.ExposureMsec);
end

assert(~isempty(exposureMsec) && isscalar(exposureMsec) && ...
    isfinite(exposureMsec) && exposureMsec > 0, ...
    'Umitoolbox:run_BloodFlow:MissingExposure', ...
    ['Could not resolve a positive exposure time from the metadata of "%s". ' ...
     'Set the "ExposureMsec" parameter explicitly.'], datFile);

end

function datFile = localResolveDatFile(SaveFolder, fileName)
%LOCALRESOLVEDATFILE Resolve a .dat basename/filename against SaveFolder.

if ~endsWith(fileName, '.dat', 'IgnoreCase', true)
    fileName = [fileName '.dat'];
end

datFile = fileName;
if ~isfile(datFile)
    datFile = fullfile(SaveFolder, fileName);
end

assert(isfile(datFile), ...
    'Umitoolbox:run_BloodFlow:FileNotFound', ...
    'Speckle .dat file not found: "%s".', fileName);

end

function tf = localIsValidKernelSize(x)
%LOCALISVALIDKERNELSIZE Accept an odd integer scalar >= 3.

tf = isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && ...
    x >= 3 && x == round(x) && mod(x, 2) == 1;

end

function tf = localIsValidExposure(x)
%LOCALISVALIDEXPOSURE Accept 'auto' or a finite positive scalar.
%
%   Numeric text (e.g. '10') is accepted too: the PipelineManager dialog
%   edits this mixed 'auto'/number parameter in a text field.

isText = ischar(x) || (isstring(x) && isscalar(x));
tf = (isText && strcmpi(strtrim(char(string(x))), 'auto')) || ...
    ~isempty(localNumericExposure(x));

end

function exposureMsec = localNumericExposure(x)
%LOCALNUMERICEXPOSURE Positive finite exposure from a number or numeric text.
%
%   Returns [] when X is 'auto' or not a valid positive number.

exposureMsec = [];
if ischar(x) || (isstring(x) && isscalar(x))
    x = str2double(char(string(x)));
end
if isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0
    exposureMsec = double(x);
end

end

function lastPct = localPrintProgress(c, nChunks, lastPct)
%LOCALPRINTPROGRESS Print integer percent progress when it changes.

pct = floor(100 * c / nChunks);
if pct ~= lastPct
    fprintf('%d%% ', pct);
    lastPct = pct;
end

end
