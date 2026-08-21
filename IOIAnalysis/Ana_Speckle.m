function varargout = Ana_Speckle(SaveFolder, bNormalize, varargin)
%ANA_SPECKLE Calculate blood-flow maps from speckle data.
%
%   Ana_Speckle(SaveFolder, bNormalize)
%   data = Ana_Speckle(SaveFolder, bNormalize)
%   [data, metaData] = Ana_Speckle(SaveFolder, bNormalize, ...)
%
%   This function computes blood-flow maps from laser speckle data using
%   the existing algorithm:
%       1) always-on removal of static structure: each frame is divided by
%          MeanMap (the per-pixel temporal mean of the raw input), which
%          otherwise contaminates the local spatial std/mean ratio below
%       2) local contrast estimation from each frame
%       3) conversion from contrast to flow using private_flow_from_contrast
%       4) temporal median filtering
%       5) optional output-level normalization of the finished flow map by
%          its own per-pixel temporal mean
%
%   Inputs:
%       SaveFolder   - Folder containing the speckle .dat file and
%                      metadata sources resolvable by loadMetaData(...).
%       bNormalize   - Logical scalar. If true, normalize the finished flow
%                      map (after temporal median filtering) by its own
%                      per-pixel temporal mean. Does not affect the
%                      always-on MeanMap correction in step 1, which runs
%                      regardless of this flag.
%
%   Name-Value parameters:
%       'Filename'    - Basename or filename of the speckle .dat input.
%                       Default: 'speckle'
%       'RAMSafeMode' - Logical scalar. If true, use low-RAM file-backed
%                       execution. Default: false
%
%   Outputs:
%       data / outFile - Standard mode returns a single Y x X x T
%                        blood-flow array. Low-RAM mode returns the full
%                        path to the raw Flow.dat output.
%       metaData       - Flat compatibility metadata describing the output.
%
%   Notes:
%       - The MeanMap (or per-frame equivalent in RAMSafe mode) correction
%         is always applied, independent of bNormalize; bNormalize affects
%         only the output-level normalization in step 5.
%       - One flow frame is calculated from each input frame, so the output
%         retains the input temporal length T and follows the normal raw
%         .dat timeline contract.
%       - Raw .dat length is resolved through loadMetaData, which infers
%         datLength from the actual file size.
%       - ExposureSpeckleMsec is read from the metadata returned by
%         loadMetaData(...).

% Default output for pipeline management.
default_Output = 'Flow.dat';

if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) && ...
        strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    varargout{1} = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'Ana_Speckle';
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addRequired(p, 'bNormalize', @(x) islogical(x) && isscalar(x));
addParameter(p, 'Filename', 'speckle', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'RAMSafeMode', false, @(x) islogical(x) && isscalar(x));
parse(p, SaveFolder, bNormalize, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
bNormalize = p.Results.bNormalize;
filename = erase(char(string(p.Results.Filename)), '.dat');
bRAMsafe = p.Results.RAMSafeMode;

fprintf('Running Ana_Speckle...\n');

% -------------------------------------------------------------------------
% Resolve input file and metadata
% -------------------------------------------------------------------------
datFile = fullfile(SaveFolder, [filename '.dat']);
if ~isfile(datFile)
    error('Ana_Speckle:InputFileNotFound', ...
        'No speckle data file was found: "%s".', datFile);
end

Iptr = loadMetaData(datFile);

requiredFields = {'Height','Width','datLength','Freq','ExposureSpeckleMsec'};
for iField = 1:numel(requiredFields)
    assert(isfield(Iptr, requiredFields{iField}) && ~isempty(Iptr.(requiredFields{iField})), ...
        'Ana_Speckle:MissingMetaData', ...
        'loadMetaData did not return required field "%s" for "%s".', ...
        requiredFields{iField}, datFile);
end

ny = double(Iptr.Height);
nx = double(Iptr.Width);
nt = double(Iptr.datLength);
tFreq = double(Iptr.Freq);
speckle_int_time = double(Iptr.ExposureSpeckleMsec) / 1000;

OPTIONS.GPU = 0;
OPTIONS.Power2Flag = 0;
OPTIONS.Brep = 0;

outMeta = struct();
outMeta.datFile = fullfile(SaveFolder, default_Output);
outMeta.datSize = [ny, nx];
outMeta.datLength = nt;
outMeta.Freq = tFreq;
outMeta.Datatype = 'single';
outMeta.dim_names = {'Y','X','T'};
outMeta.Height = ny;
outMeta.Width = nx;
outMeta.Length = nt;
outMeta.FrameRateHz = tFreq;
outMeta.ExposureSpeckleMsec = double(Iptr.ExposureSpeckleMsec);

assert(nt >= 2, 'Ana_Speckle:InvalidInputLength', ...
    'Speckle input must contain at least 2 frames.');

%% ------------------------------------------------------------------------
% Low-RAM mode
% -------------------------------------------------------------------------
if bRAMsafe
    % Compute through a fixed-name raw scratch file (bounded RAM via slab
    % I/O), then move the completed file onto the declared Flow.dat output.
    % Renaming the declared output when it already exists would make every
    % pipeline re-run write to a different file and leave the stale original
    % in place.
    outFile = fullfile(SaveFolder, default_Output);
    [~, baseName] = fileparts(default_Output);
    computeScratchFile = fullfile(SaveFolder, [baseName '_compute.dat']);
    outMeta.datFile = outFile;

    preallocateDatFile(computeScratchFile, [ny, nx, nt], 'single');

    fidIn  = fopen(datFile, 'r');
    assert(fidIn ~= -1, 'Ana_Speckle:OpenInputFailed', ...
        'Could not open input file "%s".', datFile);
    cIn = onCleanup(@() safeFclose(fidIn));

    fidOut = fopen(computeScratchFile, 'r+');
    assert(fidOut ~= -1, 'Ana_Speckle:OpenOutputFailed', ...
        'Could not open output file "%s".', computeScratchFile);
    cOut = onCleanup(@() safeFclose(fidOut));

    % Pass 1: temporal mean
    fprintf('PASS 1/3: Calculating temporal mean\n');
    MeanMap = zeros(ny, nx, 'single');
    lastPct = -1;
    frameBytes = ny * nx * getByteSize('single');

    for t = 1:nt
        fseek(fidIn, (t-1) * frameBytes, 'bof');
        frame = fread(fidIn, ny * nx, '*single');
        frame = reshape(frame, ny, nx);
        MeanMap = MeanMap + frame;

        pct = floor(100 * t / nt);
        if pct ~= lastPct && mod(pct, 10) == 0
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end
    MeanMap = MeanMap / nt;

    % Pass 2: compute flow frame-by-frame
    fprintf('\nPASS 2/3: Computing flow...\n');
    lastPct = -1;
    speckle_window = fspecial('disk', 2) > 0;

    for t = 1:nt
        fseek(fidIn, (t-1) * frameBytes, 'bof');
        frameNext = fread(fidIn, ny * nx, '*single');
        frameNext = reshape(frameNext, ny, nx);
        frameNext = frameNext ./ MeanMap;

        std_laser  = imgaussfilt(stdfilt(frameNext, speckle_window), 1);
        mean_laser = imgaussfilt(convnfft(frameNext, speckle_window, 'same', 1:2, OPTIONS) / sum(speckle_window(:)), 1);
        contrast   = std_laser ./ mean_laser;
        flow       = single(private_flow_from_contrast(contrast, speckle_int_time));

        fseek(fidOut, (t-1) * frameBytes, 'bof');
        fwrite(fidOut, flow, 'single');

        pct = floor(100 * t / nt);
        if pct ~= lastPct && mod(pct, 10) == 0
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end

    % Pass 3: temporal median filter in X chunks
    fW = ceil(0.5 * tFreq);
    nChunks = calculateMaxChunkSize(ny * nx * nt * getByteSize('single'), 2);
    chunkX = ceil(nx / nChunks);
    nChunks = ceil(nx / chunkX);

    fprintf('\nPASS 3/3: Applying temporal median filter...\n');
    lastPct = -1;

    for c = 1:nChunks
        xStart = (c-1) * chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, nx);
        xIdx   = xStart:xEnd;

        slab = spatialSlabIO('read', fidOut, ny, nx, nt, xIdx, 'single');
        slab = medfilt1(slab, fW, [], 3, 'truncate');
        spatialSlabIO('write', fidOut, ny, nx, nt, xIdx, 'single', slab);

        pct = floor(100 * c / nChunks);
        if pct ~= lastPct
            fprintf('%d%% ', pct);
            lastPct = pct;
        end
    end

    % Pass 4: output-level normalization (streaming), only when requested.
    % Reuses the same X-chunking as Pass 3 rather than a new chunk scheme.
    if bNormalize
        fprintf('\nNormalizing output by its own temporal mean...\n');
        lastPct = -1;

        for c = 1:nChunks
            xStart = (c-1) * chunkX + 1;
            xEnd   = min(xStart + chunkX - 1, nx);
            xIdx   = xStart:xEnd;

            slab = spatialSlabIO('read', fidOut, ny, nx, nt, xIdx, 'single');
            slab = slab ./ mean(slab, 3);
            spatialSlabIO('write', fidOut, ny, nx, nt, xIdx, 'single', slab);

            pct = floor(100 * c / nChunks);
            if pct ~= lastPct
                fprintf('%d%% ', pct);
                lastPct = pct;
            end
        end
    end

    clear cIn cOut; % close fidIn/fidOut before replacing the declared output

    [moveOk, moveMsg] = movefile(computeScratchFile, outFile, 'f');
    assert(moveOk, 'Ana_Speckle:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', computeScratchFile, outFile, moveMsg);

    if nargout > 0
        varargout{1} = outFile;
        if nargout > 1
            varargout{2} = outMeta;
        end
    end
    fprintf('\nDone!\n');
    return
end

%% ------------------------------------------------------------------------
% Standard mode
% -------------------------------------------------------------------------
fid = fopen(datFile, 'r');
assert(fid ~= -1, 'Ana_Speckle:OpenInputFailed', ...
    'Could not open input file "%s".', datFile);
cStd = onCleanup(@() safeFclose(fid));

dat = fread(fid, inf, '*single');
dat = reshape(dat, ny, nx, nt);
MeanMap = mean(dat, 3);

datOut = zeros(nt, ny, nx, 'single');
speckle_window = fspecial('disk', 2) > 0;

for t = 1:nt
    tmp_laser = dat(:,:,t) ./ MeanMap;

    std_laser = imgaussfilt(stdfilt(tmp_laser, speckle_window), 1);
    mean_laser = imgaussfilt(convnfft(tmp_laser, speckle_window, 'same', 1:2, OPTIONS) / sum(speckle_window(:)), 1);
    contrast = std_laser ./ mean_laser;
    datOut(t,:,:) = single(private_flow_from_contrast(contrast, speckle_int_time));
end

% Temporal median filter
tic
fW = ceil(0.5 * tFreq);
datOut = medfilt1(datOut, fW, [], 1, 'truncate');

% Finalize to Y X T
datOut = permute(datOut, [2 3 1]);

% Output-level normalization by the flow map's own per-pixel temporal
% mean. Independent of the always-on MeanMap correction applied above.
if bNormalize
    datOut = datOut ./ mean(datOut, 3);
end

outMeta.datFile = fullfile(SaveFolder, default_Output);

if nargout > 0
    varargout{1} = datOut;
    if nargout > 1
        varargout{2} = outMeta;
    end
else
    fprintf('Saving data to file: "%s"...\n', default_Output);
    saveData(outMeta.datFile, datOut);
end

fprintf('Done!\n');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'Ana_Speckle', ...
            'Calculate blood-flow maps from speckle data.');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            {'parameter'}, ...
            'Folder containing the speckle input and metadata.', ...
            'kind', 'parameter', ...
            'position', 1, ...
            'callType', 'positional', ...
            'default', '', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'bNormalize', ...
            'parameter', ...
            'If true, normalize each frame by the temporal mean.', ...
            'kind', 'parameter', ...
            'position', 2, ...
            'callType', 'positional', ...
            'default', false, ...
            'allowed', {false,true}, ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Filename', ...
            'parameter', ...
            'Basename or filename of the speckle input.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', 'speckle', ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'RAMSafeMode', ...
            'parameter', ...
            'If true, use low-RAM file-backed execution.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', false, ...
            'allowed', {false,true}, ...
            'dataType', 'logical');

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            {'ImageTimeSeries', 'ProcessedData'}, ...
            'data', ...
            'Blood-flow output: a single Y x X x T image time series.', ...
            default_Output, ...
            1, ...
            'isData', true, ...
            'saveFileName', default_Output);
    end
end

%% Local function
function speed = private_flow_from_contrast(contrast,T)
contrast(isnan(contrast)|contrast<0)=0;
contrast2 = contrast(3:end-2,3:end-2);
mmean = mean(contrast2(:));
sstd  = std(contrast2(:));
tau=(logspace(-15,0,60).^.5);
K  = ((tau/(2*T)).*(1-exp(-2*T*ones(size(tau))./tau))).^(1/2);
[~, index1] = find(K>(mmean-3*sstd),1);
[~, index2] = find(K>(mmean+3*sstd),1);
if isempty(index1), index1=1; end
if isempty(index2)||index2==index1, index2=60; end
Tau2=(logspace(log10(tau(index1)),log10(tau(index2)),40));
K  = ((Tau2/(2*T)).*(1-exp(-2*T*ones(size(Tau2))./Tau2))).^(1/2);
Tau2=[Tau2(1) Tau2 Tau2(end)];
K=[0 K 1e30];
speed=1./interp1(K,Tau2,contrast);
end
