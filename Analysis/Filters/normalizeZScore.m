function outData = normalizeZScore(data, SaveFolder, varargin)
%NORMALIZEZSCORE Normalize image time-series data to z-scores along T.
%
%   outData = normalizeZScore(data, SaveFolder)
%
%   This function normalizes image time-series data along the time
%   dimension using z-score normalization:
%
%       z = (x - mean(x, T)) ./ std(x, T)
%
%   The statistics are computed independently for every pixel and, for
%   event-split data, for every event instance (E slice), so each trial is
%   normalized over its own frames.
%
%   Supported execution modes:
%       1) STANDARD MODE (in-memory)
%          - Triggered when "data" is a numeric array
%       2) LOW-RAM MODE (file-backed)
%          - Triggered when "data" is a .dat filename
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X x T
%       2) Numeric array with dimensions Y x X x T x E (event-split data)
%       3) Filename to a .dat file with axes Y-X-T or Y-X-T-E
%
%   Input/output behavior:
%       - If the input is a numeric array, the output is a numeric array
%         with the same size.
%       - If the input is a .dat filename, the output is a .dat filename
%         ("normZ.dat" in SaveFolder) with the same axes as the input.
%       - UMT structs and .umt files are not supported.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder used for file resolution and outputs.
%
%   Output:
%       outData    - Z-score normalized output with the same
%                    representation type as the input.
%
%   Notes:
%       - In low-RAM mode, the input file is processed chunk-by-chunk along
%         the X dimension; the chunk size accounts for the E axis.
%       - NaN values are ignored when computing the mean and standard
%         deviation.
%       - For constant time traces, standard deviation is forced to 1 to
%         avoid division by zero, preserving the original algorithm.
%       - Every E slice is normalized, ignored event instances included, so
%         the E axis of a .dat output still matches events.mat.
%       - Raw .dat files are assumed to store data in single precision.

default_Output = 'normZ.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));

if ~isfolder(SaveFolder)
    error('normalizeZScore:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Case 1: In-memory numeric array (Y x X x T or Y x X x T x E)
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty'}, ...
        mfilename, 'data');
    if ~(ndims(data) == 3 || ndims(data) == 4)
        error('normalizeZScore:InvalidArrayInput', ...
            'Numeric input must be YXT or YXTE.');
    end

    outData = iZscoreAlongT(data);
    return
end

% -------------------------------------------------------------------------
% Case 2: File input
% -------------------------------------------------------------------------
if ischar(data) || (isstring(data) && isscalar(data))

    dataFile = char(string(data));

    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('normalizeZScore:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);

    if ~strcmpi(ext, '.dat')
        error('normalizeZScore:UnsupportedInputFile', ...
            'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    end

    outData = iZscoreDatFile(dataFile, SaveFolder, default_Output);
    return
end

error('normalizeZScore:UnsupportedInputType', ...
    ['Input "data" must be a YXT or YXTE array or a .dat filename. ' ...
     'UMT structs and .umt files are not supported.']);

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            ['Normalize image time-series data to zero mean and unit ' ...
             'standard deviation along the time dimension, per event ' ...
             'instance for event-split data.']);

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            'Input data. Accepted forms: YXT or YXTE array, or .dat filename.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder used for file resolution and outputs.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Z-score normalized output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: z-score along the 3rd dimension (T), independently per pixel and E
% =========================================================================
function outArray = iZscoreAlongT(inArray)
%IZSCOREALONGT Z-score a Y x X x T (x E) array along T.

origSz = size(inArray);
workArray = reshape(inArray, origSz(1) * origSz(2), origSz(3), []);

mu  = mean(workArray, 2, 'omitnan');
sig = std(workArray, 0, 2, 'omitnan');
sig(sig == 0) = 1;

outArray = reshape((workArray - mu) ./ sig, origSz);
end

% =========================================================================
% Helper: low-RAM .dat execution
% =========================================================================
function outFile = iZscoreDatFile(inFile, SaveFolder, defaultOutput)
%IZSCOREDATFILE Apply z-score normalization to a Y-X-T or Y-X-T-E .dat file.

slabIn = spatialSlabIO('open', inFile);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));
assertDatLayout(slabIn.Info, {{'Y','X','T'}, {'Y','X','T','E'}}, 'normalizeZScore');
Ny = slabIn.Ny;
Nx = slabIn.Nx;
Nt = datAxisSize(slabIn.Info, 'T');
Ne = 1;
if any(strcmp(cellstr(string(slabIn.Info.dimNames)), 'E'))
    Ne = datAxisSize(slabIn.Info, 'E');
end

% Write through a scratch file so the declared pipeline output only
% appears once the run has completed, and so the input can safely be the
% file that the declared output would overwrite (a pipeline re-run).
outFile = fullfile(SaveFolder, defaultOutput);
[~, outStem, outExt] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [outStem '_writing' outExt]);
slabOut = spatialSlabIO('create', tmpFile, ...
    datHeaderFromInfo(slabIn.Info, outStem, 'dataClass', 'single'));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

totalBytes = Ny * Nx * Nt * Ne * getByteSize('single');
nChunks = calculateMaxChunkSize(totalBytes, 2);
chunkX = ceil(Nx / nChunks);
nChunks = ceil(Nx / chunkX);

for c = 1:nChunks
    xStart = (c-1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    slab = single(spatialSlabIO('read', slabIn, xIdx));
    slab = iZscoreAlongT(slab);

    spatialSlabIO('write', slabOut, xIdx, slab);
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close the input reader before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'normalizeZScore:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
end
