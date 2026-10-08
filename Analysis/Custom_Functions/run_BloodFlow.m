function outData = run_BloodFlow(SaveFolder, data, varargin)
%RUN_BLOODFLOW Estimate time-resolved blood flow from raw speckle data.
%
%   outData = run_BloodFlow(SaveFolder, data)
%   outData = run_BloodFlow(SaveFolder, data, 'Name', Value, ...)
%   info    = run_BloodFlow('pipelineInfo')
%
%   This function computes a Laser Speckle Contrast Image (LSCI) with the
%   same contrast algorithm as SPECKLEMAPPING, but without averaging the
%   contrast over time. The input Y-X-T shape is therefore preserved.
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
%       SaveFolder - Folder of the speckle file and of the .dat output.
%       data       - Name of a Y-X-T speckle .dat file (looked up in
%                    SaveFolder, or used as-is when it is a full path).
%                    The file is processed in chunks and the exposure time is
%                    resolved from its header. Arrays and other layouts are
%                    not supported.
%
%   Name-Value parameters:
%       sType           - Speckle contrast mode.
%                         Allowed: 'Spatial', 'Temporal'
%                         Default: 'Temporal'
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
%       outData    - Full path to the Y-X-T BloodFlow.dat output in
%                    SaveFolder (single precision, same axes and rate as the
%                    input).
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
addRequired(p, 'data');
addParameter(p, 'sType', 'Temporal', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    any(strcmpi(char(string(x)), allowedSType)));
addParameter(p, 'ExposureMsec', 'auto', @localIsValidExposure);
addParameter(p, 'KernelSize', 5, @localIsValidKernelSize);
addParameter(p, 'bNormalize', false, @(x) islogical(x) && isscalar(x));
parse(p, SaveFolder, data, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
sType = lower(char(string(p.Results.sType)));
exposureOpt = p.Results.ExposureMsec;
kernelSize = double(p.Results.KernelSize);
bNormalize = p.Results.bNormalize;

assert(ischar(data) || (isstring(data) && isscalar(data)), ...
    'Umitoolbox:run_BloodFlow:UnsupportedInputType', ...
    'Input "data" must be the name or path of a Y-X-T .dat file (arrays are not supported).');

datFile = localResolveDatFile(SaveFolder, char(string(data)));
mdIn = loadMetaData(datFile);
assertDatLayout(mdIn, {{'Y','X','T'}}, 'run_BloodFlow');

exposureSec = localResolveExposure(mdIn, exposureOpt, datFile) / 1000;

fprintf('Calculating blood flow (%s speckle contrast, kernel size %d)...\n', ...
    sType, kernelSize);

outData = localRunDatFile(datFile, mdIn, SaveFolder, default_Output, ...
    sType, kernelSize, exposureSec, bNormalize);

fprintf('Finished Blood Flow.\n');

    % =====================================================================
    % Local pipelineInfo factory (nested, shares allowedSType and
    % default_Output with the parent scope instead of redefining them)
    % =====================================================================
    function info = localPipelineInfo()

        info = PipelineManager.createPipelineInfo( ...
            mfilename, ...
            'Estimate time-resolved blood flow (1/(T*K^2)) from speckle contrast.');
        info.version = '2.0.0';

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
            'Y-X-T speckle .dat file.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'file');

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
            'ExposureMsec', ...
            'parameter', ...
            'Exposure time in ms. "auto" reads it from the file metadata.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', 'auto', ...
            'allowed', {'auto', [0 Inf]});

        info = PipelineManager.addInput(info, ...
            'KernelSize', ...
            'parameter', ...
            'Contrast window size (odd, >= 3): disk of diameter N (Spatial) or N frames (Temporal).', ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'namevalue', ...
            'default', 5, ...
            'allowed', [3 Inf], ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'bNormalize', ...
            'parameter', ...
            'If true, normalize the flow by its own per-pixel temporal mean.', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', false, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            {'ImageTimeSeries', 'ProcessedData'}, ...
            'data', ...
            'Blood-flow output (1/s, or relative if bNormalize): a single Y-X-T .dat image time series.', ...
            default_Output, ...
            1, ...
            'isData', true, ...
            'saveFileName', default_Output);
    end

end

% =========================================================================
% Local functions
% =========================================================================
function outFile = localRunDatFile(datFile, mdIn, SaveFolder, default_Output, ...
    sType, kernelSize, exposureSec, bNormalize)
%LOCALRUNDATFILE Chunked, file-backed blood-flow computation.

Ny = datAxisSize(mdIn, 'Y');
Nx = datAxisSize(mdIn, 'X');
Nt = datAxisSize(mdIn, 'T');

% Compute through a fixed-name scratch file, then move it onto the declared
% output so re-runs overwrite the same file.
outFile = fullfile(SaveFolder, default_Output);
[~, baseName] = fileparts(default_Output);
scratchFile = fullfile(SaveFolder, [baseName '_compute.dat']);

slabIn = spatialSlabIO('open', datFile, 'Info', mdIn);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

slabOut = spatialSlabIO('create', scratchFile, ...
    datHeaderFromInfo(mdIn, baseName, 'dataClass', 'single'));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

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
            [slab, ~] = localReadInputFrames(slabIn, c, chunkT, Nt);
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
            [slab, tStart] = localReadInputFrames(slabIn, c, chunkT, Nt);
            slab = localSpeckleContrast(slab, meanData, 'spatial', kernelSize);
            slab = localFlowFromContrast(slab, exposureSec);
            if bNormalize
                sumFlow = sumFlow + sum(double(slab), 3, 'omitnan');
                countFlow = countFlow + sum(~isnan(slab), 3);
            end

            spatialSlabIO('write', slabOut, 1:Nx, slab, tStart:tStart + size(slab, 3) - 1);
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
                [slab, tStart] = localReadInputFrames(slabOut, c, chunkT, Nt);
                slab = single(slab ./ meanFlow);

                spatialSlabIO('write', slabOut, 1:Nx, slab, tStart:tStart + size(slab, 3) - 1);
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

            slab = double(spatialSlabIO('read', slabIn, xIdx));
            slab = localSpeckleContrast(slab, mean(slab, 3, 'omitnan'), 'temporal', kernelSize);
            slab = localFlowFromContrast(slab, exposureSec);
            if bNormalize
                % The slab holds the full time course of its pixels.
                slab = localNormalizeByTemporalMean(slab);
            end
            spatialSlabIO('write', slabOut, xIdx, slab);
            lastPct = localPrintProgress(c, nChunks, lastPct);
        end
end
fprintf('\n');

spatialSlabIO('finalize', slabOut);
clear cIn cOut % close both files before replacing the declared output

[moveOk, moveMsg] = movefile(scratchFile, outFile, 'f');
assert(moveOk, 'Umitoolbox:run_BloodFlow:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', scratchFile, outFile, moveMsg);

end

function [slab, tStart] = localReadInputFrames(slabIn, c, chunkT, Nt)
%LOCALREADINPUTFRAMES Read one temporal chunk of frames as double.
%
% Used on the input and, through its create handle, on this function's
% own scratch output.

tStart = (c-1) * chunkT + 1;
tEnd = min(tStart + chunkT - 1, Nt);
slab = double(spatialSlabIO('read', slabIn, 1:slabIn.Nx, tStart:tEnd));
slab = reshape(slab, slabIn.Ny, slabIn.Nx, tEnd - tStart + 1);

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
% exposureMsec is the file's own exposure (the speckle exposure for
% speckle data); NaN when unknown.
if isfield(mdIn, 'exposureMsec') && ~isempty(mdIn.exposureMsec) && ~isnan(mdIn.exposureMsec)
    exposureMsec = double(mdIn.exposureMsec);
end

assert(~isempty(exposureMsec) && isscalar(exposureMsec) && ...
    isfinite(exposureMsec) && exposureMsec > 0, ...
    'Umitoolbox:run_BloodFlow:MissingExposure', ...
    ['Could not resolve a positive exposure time from the metadata of "%s". ' ...
     'Set the "ExposureMsec" parameter explicitly.'], datFile);

end

function datFile = localResolveDatFile(SaveFolder, fileName)
%LOCALRESOLVEDATFILE Resolve a .dat basename/filename against SaveFolder.

[~, ~, ext] = fileparts(fileName);
assert(isempty(ext) || strcmpi(ext, '.dat'), ...
    'Umitoolbox:run_BloodFlow:UnsupportedInputFile', ...
    'Unsupported input file extension "%s". Only .dat files are supported.', ext);
if isempty(ext)
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
