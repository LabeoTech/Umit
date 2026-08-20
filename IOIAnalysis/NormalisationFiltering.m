function varargout = NormalisationFiltering(FolderData, FileData, lowFreq, ...
    highFreq, bDivide, bExpfit, varargin)
%NORMALISATIONFILTERING Data normalization by low-pass filtering.
%
% This function can be used to normalize channels (delta F/F or delta R/R),
% or to perform low-pass filtering.
%
% Inputs:
%
% Option A: data to be normalized are opened within this function
%   1- FolderData:  Folder containing the data to be opened
%   2- FileData:    Channel/file to open
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore
%   5- bDivide:     if true, the high-passed signal is normalized by the
%                   low-passed signal; if false, the low-passed signal is
%                   subtracted
%   6- bExpFit:     if true, apply double exponential decay correction
%
% Option B: data to be normalized are given directly
%   1- FolderData:  Folder containing the data
%   2- FileData:    Data as a 3D matrix (Y, X, Time)
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore
%   5- bDivide:     if true, the high-passed signal is normalized by the
%                   low-passed signal; if false, the low-passed signal is
%                   subtracted
%   6- bExpFit:     if true, apply double exponential decay correction
%   7- Freq:        (Optional) sample rate
%
% Optional inputs:
%   Freq         - Sample rate override
%   saveFilename - Output filename in file mode
%
% Notes:
%   - File mode uses loadMetaData(...) to retrieve compatibility metadata.
%   - If dim_names are missing:
%         * non-event data assume legacy order {'Y','X','T'}
%         * event data assume legacy order {'E','Y','X','T'}
%     and a warning is raised.
%   - Non-event file data are expected to be stored as Y,X,T on disk.
%   - The core filtering algorithm is unchanged.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'FolderData', @isfolder)
addRequired(p, 'FileData', @(x) ischar(x) || isstring(x) || isnumeric(x))
addRequired(p, 'lowFreq', @(x) isscalar(x) && isnumeric(x))
addRequired(p, 'highFreq', @(x) isscalar(x) && isnumeric(x))
addRequired(p, 'bDivide', @(x) isscalar(x) && (isnumeric(x) || islogical(x)))
addRequired(p, 'bExpFit', @(x) isscalar(x) && (isnumeric(x) || islogical(x)))
addOptional(p, 'Freq', [], @(x) (isscalar(x) && isnumeric(x)) || isempty(x))
addOptional(p, 'saveFilename', '', @(x) ischar(x) || (isstring(x) && isscalar(x)))

parse(p, FolderData, FileData, lowFreq, highFreq, bDivide, bExpfit, varargin{:});

lowFreq   = p.Results.lowFreq;
highFreq  = p.Results.highFreq;
bDivide   = p.Results.bDivide;
bExpFit   = p.Results.bExpFit;
freqArg   = p.Results.Freq;
saveFile  = char(string(p.Results.saveFilename));

if ischar(p.Results.FileData) || (isstring(p.Results.FileData) && isscalar(p.Results.FileData))
    if nargout
        OutData = NormFiltFromFile( ...
            p.Results.FolderData, ...
            char(string(p.Results.FileData)), ...
            lowFreq, highFreq, bDivide, bExpFit, ...
            true, saveFile, freqArg);
        varargout{1} = OutData;
    else
        NormFiltFromFile( ...
            p.Results.FolderData, ...
            char(string(p.Results.FileData)), ...
            lowFreq, highFreq, bDivide, bExpFit, ...
            false, saveFile, freqArg);
    end
else
    OutData = NormFiltDirect( ...
        p.Results.FolderData, ...
        p.Results.FileData, ...
        lowFreq, highFreq, bDivide, bExpFit, ...
        freqArg);
    varargout{1} = OutData;
end

end


function OutData = NormFiltDirect(FolderData, OutData, lowFreq, highFreq, bDivide, bExpFit, Freq)
%NORMFILTDIRECT Filter an in-memory YXT array.

ExpFun = @(P,x) abs(P(1)).*exp(-abs(P(2)).*x) + abs(P(3)).*exp(-abs(P(4)).*x) ...
    - abs(P(5))*x + P(6);
Opt = optimset(@fminsearch);
Opt.Display = 'off';

if isempty(Freq)
    datFiles = dir(fullfile(FolderData, '*.dat'));
    umtFiles = dir(fullfile(FolderData, '*.umt'));

    if ~isempty(datFiles)
        meta = loadMetaData(fullfile(FolderData, datFiles(1).name));
        Freq = meta.Freq;
    elseif ~isempty(umtFiles)
        meta = loadMetaData(fullfile(FolderData, umtFiles(1).name));
        Freq = meta.Freq;
    else
        error('NormalisationFiltering:MissingFrequency', ...
            ['Freq was not provided and no supported data file (*.dat or *.umt) ' ...
             'was found in "%s" to infer it.'], ...
            FolderData);
    end
end

% Temporal filtering
if lowFreq > 0
    UseLPFilt = 1;
    if (1/lowFreq) > (size(OutData,3)/Freq)
        lowFreq = 1 / (2 * (size(OutData,3)/Freq));
    end
    f = fdesign.lowpass('N,F3dB', 4, lowFreq, Freq);
    lpass = design(f, 'butter');
else
    UseLPFilt = 0;
end

if highFreq > 0
    UseHPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, highFreq, Freq);
    hpass = design(f, 'butter');
else
    UseHPFilt = 0;
end

dims = size(OutData);
if bExpFit
    rng('shuffle');
    S = mean(reshape(OutData, [], dims(3)), 1);
    initB = double(rand(1,6) .* [30 1 20 1 1 mean(double(S))]);
    B = fminsearch(@(P) sum((double(S) - ExpFun(P, (1:size(S,2)))).^2), ...
        initB, Opt);
    Approx = ExpFun([B(1:4) 0 0], 1:size(S,2));
    Pred = [ones(1, size(S,2)); linspace(0,1,size(S,2)); Approx]';
end

PrcLims = round(linspace(1, dims(1), 11));
fprintf('Progress: ');
for ind = 1:dims(1)
    Signal = double(squeeze(OutData(ind,:,:)));

    if bExpFit
        B = Pred \ Signal';
        Approx = (Pred * B)';
        Signal = Signal ./ Approx;
    end

    if UseLPFilt
        LP_lowCutOff = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
    else
        LP_lowCutOff = ones(size(Signal));
    end

    if UseHPFilt
        LP_highCutOff = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
    else
        LP_highCutOff = Signal;
    end

    if bDivide
        OutData(ind,:,:) = single(LP_highCutOff ./ LP_lowCutOff);
    else
        OutData(ind,:,:) = single(LP_highCutOff - LP_lowCutOff);
    end

    if any(ind == PrcLims)
        idx = find(ind == PrcLims);
        fprintf('%d%%..', 10*(idx-1));
    end
end
fprintf('\n');
end


function OutData = NormFiltFromFile(FolderData, FileName, lowFreq, highFreq, bDivide, bExpFit, bReturn, outFile, freqOverride)
%NORMFILTFROMFILE Filter data from disk using compatibility metadata.

if nargin < 9
    freqOverride = [];
end

% -------------------------------------------------------------------------
% Resolve input file path
% -------------------------------------------------------------------------
if isempty(fileparts(FileName))
    inFile = fullfile(FolderData, FileName);
else
    inFile = FileName;
end

if ~isfile(inFile)
    error('NormalisationFiltering:InputFileNotFound', ...
        'Input file not found: "%s".', inFile);
end

[inFolder, inBase, ~] = fileparts(inFile);

% -------------------------------------------------------------------------
% Load metadata through compatibility layer
% -------------------------------------------------------------------------
fMetaData = loadMetaData(inFile);

if ~isfield(fMetaData, 'datSize') || ~isfield(fMetaData, 'datLength') || ...
        ~isfield(fMetaData, 'Datatype')
    error('NormalisationFiltering:InvalidMetaData', ...
        ['loadMetaData("%s") did not return the required fields: ' ...
         'datSize, datLength, Datatype.'], ...
        inFile);
end

storedSize = [fMetaData.datSize fMetaData.datLength];
dataType = char(string(fMetaData.Datatype));

if isfield(fMetaData, 'dim_names') && ~isempty(fMetaData.dim_names)
    dimNames = cellstr(string(fMetaData.dim_names));
else
    if numel(storedSize) == 4
        dimNames = {'E','Y','X','T'};
    elseif numel(storedSize) == 3
        dimNames = {'Y','X','T'};
    else
        error('NormalisationFiltering:MissingDimNames', ...
            'Failed to infer dim_names for "%s".', inFile);
    end
    warning('NormalisationFiltering:MissingDimNames', ...
        ['Metadata field dim_names was missing for "%s". ' ...
         'Assuming legacy order {%s}.'], ...
        inFile, strjoin(dimNames, ','));
end

if numel(dimNames) ~= numel(storedSize)
    if numel(storedSize) == 4
        dimNames = {'E','Y','X','T'};
    elseif numel(storedSize) == 3
        dimNames = {'Y','X','T'};
    else
        error('NormalisationFiltering:InvalidDimNames', ...
            'dim_names length does not match stored data size for "%s".', inFile);
    end
    warning('NormalisationFiltering:InvalidDimNames', ...
        ['Metadata dim_names length did not match stored size for "%s". ' ...
         'Assuming legacy order {%s}.'], ...
        inFile, strjoin(dimNames, ','));
end

if isfield(fMetaData, 'Freq') && ~isempty(fMetaData.Freq)
    Freq = fMetaData.Freq;
elseif ~isempty(freqOverride)
    Freq = freqOverride;
else
    error('NormalisationFiltering:MissingFrequency', ...
        'Failed to determine Freq for "%s".', inFile);
end

hasEvents = any(strcmpi(dimNames, 'E'));
idxT = find(strcmpi(dimNames, 'T'), 1, 'first');
assert(~isempty(idxT), 'NormalisationFiltering:MissingTimeDimension', ...
    'dim_names must contain a T dimension.');

% -------------------------------------------------------------------------
% Resolve output filename
% -------------------------------------------------------------------------
if isempty(outFile)
    outDat = fullfile(inFolder, [inBase '_NormFilt.dat']);
else
    outFile = char(string(outFile));
    [outFolder, outBase, outExt] = fileparts(outFile);

    if isempty(outFolder)
        outFolder = inFolder;
    end
    if isempty(outExt)
        outExt = '.dat';
    end
    if ~isfolder(outFolder)
        mkdir(outFolder);
    end

    outDat = fullfile(outFolder, [outBase outExt]);
end

% -------------------------------------------------------------------------
% Preallocate output file
% -------------------------------------------------------------------------
preallocateDatFile(outDat, storedSize, dataType);

fidIn  = fopen(inFile, 'r');
assert(fidIn ~= -1, 'NormalisationFiltering:OpenInputFailed', ...
    'Failed to open input file "%s".', inFile);
cIn = onCleanup(@() safeFclose(fidIn));

fidOut = fopen(outDat, 'r+');
assert(fidOut ~= -1, 'NormalisationFiltering:OpenOutputFailed', ...
    'Failed to open output file "%s".', outDat);
cOut = onCleanup(@() safeFclose(fidOut));

% ---------------- Filter design ----------------
if lowFreq > 0
    nTimeSamples = storedSize(idxT);
    if (1/lowFreq) > (nTimeSamples / Freq)
        lowFreq = 1 / (2 * (nTimeSamples / Freq));
    end

    f = fdesign.lowpass('N,F3dB', 4, lowFreq, Freq);
    lpass = design(f, 'butter');
    UseLPFilt = true;
else
    UseLPFilt = false;
end

if highFreq > 0
    f = fdesign.lowpass('N,F3dB', 4, highFreq, Freq);
    hpass = design(f, 'butter');
    UseHPFilt = true;
else
    UseHPFilt = false;
end

% ---------------- Exponential fit helpers ----------------
ExpFun = @(P,x) abs(P(1)).*exp(-abs(P(2)).*x) + ...
                abs(P(3)).*exp(-abs(P(4)).*x) ...
                - abs(P(5))*x + P(6);
Opt = optimset('Display', 'off');

fprintf('NormalisationFiltering (file mode)\n');

%% =========================================================
% EVENT MODE
%% =========================================================
if hasEvents

    idxY = find(strcmpi(dimNames, 'Y'), 1, 'first');
    idxX = find(strcmpi(dimNames, 'X'), 1, 'first');
    idxE = find(strcmpi(dimNames, 'E'), 1, 'first');

    assert(~isempty(idxY) && ~isempty(idxX) && ~isempty(idxE), ...
        'NormalisationFiltering:InvalidDimNames', ...
        'Event data must define Y, X, T, and E dimensions.');

    nElem = prod(storedSize);
    raw = fread(fidIn, nElem, ['*' dataType]);
    assert(numel(raw) == nElem, 'NormalisationFiltering:UnexpectedEOF', ...
        'Unexpected end of file while reading "%s".', inFile);

    storedData = reshape(raw, storedSize);
    permToYXTE = [idxY idxX idxT idxE];
    dataYXTE = permute(storedData, permToYXTE);

    [Ny, ~, Nt, Ne] = size(dataYXTE);

    for e = 1:Ne
        fprintf('Trial %d / %d\n', e, Ne);
        slab = dataYXTE(:,:,:,e);

        if bExpFit
            S = mean(reshape(slab, [], Nt), 1);
            initB = double(rand(1,6) .* [30 1 20 1 1 mean(double(S))]);
            B = fminsearch(@(P) sum((double(S) - ExpFun(P, 1:Nt)).^2), ...
                initB, Opt);
            Approx = ExpFun([B(1:4) 0 0], 1:Nt);
            Pred = [ones(1,Nt); linspace(0,1,Nt); Approx]';
        end

        for y = 1:Ny
            Signal = double(squeeze(slab(y,:,:)));

            if bExpFit
                B = Pred \ Signal';
                Signal = Signal ./ (Pred * B)';
                clear B
            end

            if UseLPFilt
                LP_low = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
            else
                LP_low = ones(size(Signal));
            end

            if UseHPFilt
                LP_high = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
            else
                LP_high = Signal;
            end

            if bDivide
                slab(y,:,:) = cast(LP_high ./ LP_low, dataType);
            else
                slab(y,:,:) = cast(LP_high - LP_low, dataType);
            end
        end

        dataYXTE(:,:,:,e) = slab;
    end

    storedOut = ipermute(dataYXTE, permToYXTE);
    fwrite(fidOut, storedOut, dataType);

%% =========================================================
% NO EVENTS: Y,X,T only
%% =========================================================
else
    if ~isequal(dimNames, {'Y','X','T'})
        error('NormalisationFiltering:InvalidDimNames', ...
            ['Non-event file data must use dim_names {''Y'',''X'',''T''}. ' ...
             'Received {%s} for "%s".'], ...
            strjoin(dimNames, ','), inFile);
    end

    Ny = storedSize(1);
    Nx = storedSize(2);
    Nt = storedSize(3);

    if bExpFit
        nChunks = calculateMaxChunkSize(prod(storedSize) * getByteSize(dataType), 2, .3);
    else
        nChunks = calculateMaxChunkSize(prod(storedSize) * getByteSize(dataType), 1, .3);
    end

    chunkX = ceil(Nx / nChunks);
    nChunks = ceil(Nx / chunkX);
    fprintf('Filtering data (%i chunk(s))...\n', nChunks)

    % Global exponential-fit pre-pass (F-13/F-14). The detrend fit must be
    % one whole-image fit shared by every chunk, not recomputed per chunk --
    % otherwise the fit (and therefore the output) depends on nChunks, which
    % itself varies with available RAM via calculateMaxChunkSize. Stream the
    % whole-image mean trace across X-chunks via a running NaN-safe sum plus
    % a running per-time-point count of non-NaN pixels, the same technique as
    % the F-8 HemoCompute clamp-floor pre-pass and the F-1/F-2 SpeckleMapping
    % fix, so this stays RAM-safe and so masked pixels do not drag the fit.
    if bExpFit
        sumS = zeros(1, Nt);
        countS = zeros(1, Nt);
        for c = 1:nChunks
            xStart = (c-1) * chunkX + 1;
            xEnd   = min(xStart + chunkX - 1, Nx);
            xIdx   = xStart:xEnd;

            slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, dataType);
            slab = reshape(double(slab), [], Nt);
            isValidS = ~isnan(slab);
            slab(~isValidS) = 0;
            sumS = sumS + sum(slab, 1);
            countS = countS + sum(isValidS, 1);
        end
        clear slab isValidS

        S = sumS ./ countS;
        rng('shuffle');
        initB = double(rand(1,6) .* [30 1 20 1 1 mean(S)]);
        B = fminsearch(@(P) sum((S - ExpFun(P, 1:Nt)).^2), initB, Opt);
        Approx = ExpFun([B(1:4) 0 0], 1:Nt);
        Pred = [ones(1,Nt); linspace(0,1,Nt); Approx]';
    end

    for c = 1:nChunks
        xStart = (c-1) * chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, Nx);
        xIdx   = xStart:xEnd;

        fprintf('Chunk #%i [Reading data from file...]\n', c)
        slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, dataType);

        % F-14: mask NaNs before filtering and restore them afterward, the
        % same policy iFilterArray applies for in-RAM inputs
        % (Analysis/Filters/normalizeLPF.m). Per-chunk masking is exact here
        % since the mask is per-element and filtering runs per pixel along T.
        idxNaN = isnan(slab);
        if any(idxNaN(:))
            slab(idxNaN) = 0;
        end

        fprintf('Chunk #%i [Temporal filtering...]\n', c)
        fprintf('Progress:')
        lastPct = -1;

        for y = 1:Ny
            Signal = double(squeeze(slab(y,:,:)));

            if bExpFit
                B = Pred \ Signal';
                Signal = Signal ./ (Pred * B)';
                clear B
            end

            if UseLPFilt
                LP_low = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
            else
                LP_low = ones(size(Signal));
            end

            if UseHPFilt
                LP_high = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
            else
                LP_high = Signal;
            end

            if bDivide
                slab(y,:,:) = cast(LP_high ./ LP_low, dataType);
            else
                slab(y,:,:) = cast(LP_high - LP_low, dataType);
            end

            pct = floor(100 * y / Ny);
            if pct ~= lastPct && mod(pct, 10) == 0
                fprintf('%d%% ', pct);
                lastPct = pct;
            end
        end

        if any(idxNaN(:))
            slab(idxNaN) = NaN;
        end

        fprintf('\nChunk #%i [Writing to file...]\n', c)
        spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, dataType, slab);
        fprintf('Chunk #%i [Completed]\n', c)
        clear slab
    end
end

fclose(fidIn);
fclose(fidOut);

if bReturn
    fid = fopen(outDat, 'r');
    assert(fid ~= -1, 'NormalisationFiltering:OpenReturnFileFailed', ...
        'Failed to reopen output file "%s".', outDat);
    OutData = fread(fid, inf, ['*' dataType]);
    fclose(fid);
    OutData = reshape(OutData, storedSize);

    if exist(outDat, 'file')
        delete(outDat)
    end
else
    OutData = [];
end

end
