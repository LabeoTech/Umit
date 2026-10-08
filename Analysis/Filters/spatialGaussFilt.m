function outData = spatialGaussFilt(data, SaveFolder, varargin)
%SPATIALGAUSSFILT Apply a spatial Gaussian filter to image data.
%
%   outData = spatialGaussFilt(data, SaveFolder)
%   outData = spatialGaussFilt(data, SaveFolder, 'Sigma', sigma)
%
%   This function applies a spatial Gaussian filter using IMGAUSSFILT. Every
%   Y x X frame is filtered on its own; frames are never mixed along T or E.
%
%   Supported execution modes:
%       1) STANDARD MODE (in-memory)
%          - Triggered when "data" is a numeric array
%       2) LOW-RAM MODE (file-backed)
%          - Triggered when "data" is a .dat filename
%          - The file is read, filtered, and written in blocks of whole
%            frames, so the recording is never loaded whole; progress is
%            printed to the command window, one line per step and block
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X x T
%       2) Numeric array with dimensions Y x X x T x E (event-split data)
%       3) Filename to a .dat file with axes Y-X-T or Y-X-T-E
%       UMT structs and .umt files are not supported.
%
%   Input/output behavior:
%       - If the input is a numeric array, the output is a numeric array
%         with the same size (single or double as the input; integer and
%         logical arrays are converted to single).
%       - If the input is a .dat filename, the output is a .dat filename
%         ("spatialGaussFilt.dat" in SaveFolder) with the same axes and sizes.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder used for file resolution and outputs.
%
%   Name-Value parameters:
%       Sigma      - Positive scalar Gaussian sigma. Default: 1
%
%   Output:
%       outData    - Filtered output with the same representation type and
%                    dimensions as the input.
%
%   Notes:
%       - Raw .dat files are assumed to store single-precision data.
%       - NaN values are replaced by zero before filtering and restored
%         afterward, preserving the original algorithm behavior.
%       - NaN masks are tracked per frame: a pixel that is NaN in one frame
%         does not cause valid data in other frames to be discarded, so the
%         result does not depend on how the frames are blocked.

default_Output = 'spatialGaussFilt.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Sigma', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
Sigma = double(p.Results.Sigma);

if ~isfolder(SaveFolder)
    error('spatialGaussFilt:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Case 1: YXT or YXTE array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty'}, mfilename, 'data');
    if ~(ndims(data) == 3 || ndims(data) == 4)
        error('spatialGaussFilt:InvalidArrayInput', ...
            'Numeric input must be YXT or YXTE.');
    end

    if ~isfloat(data)
        % imgaussfilt returns the class of its input: an integer array would
        % be rounded back to integers.
        data = single(data);
    end

    % Trailing axes (T, E) are flattened into frames and restored afterwards,
    % so the whole array is filtered by one call.
    inSize = size(data);
    outData = iFilterFrames(reshape(data, inSize(1), inSize(2), []), Sigma);
    outData = reshape(outData, inSize);
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
            error('spatialGaussFilt:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);

    switch ext
        case '.dat'
            outData = iSpatialGaussDatFile(dataFile, SaveFolder, Sigma, default_Output);
            return

        otherwise
            error('spatialGaussFilt:UnsupportedInputFile', ...
                'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    end
end

error('spatialGaussFilt:UnsupportedInputType', ...
    ['Input "data" must be a YXT or YXTE array or a .dat filename. ' ...
     'UMT structs and .umt files are not supported.']);

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Apply a spatial Gaussian filter to image data.');

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'Image','ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT or YXTE array, or a .dat ' ...
             'filename with axes Y-X-T or Y-X-T-E.'], ...
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
            'Folder used for file resolution and the .dat output.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'Sigma', ...
            'parameter', ...
            'Positive scalar Gaussian sigma.', ...
            'kind', 'parameter', ...
            'default', 1, ...
            'allowed', [eps, Inf], ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'Image','ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Spatially filtered output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Local helper: filter a stack of frames, keeping NaN pixels out of the kernel
% =========================================================================
function outBlock = iFilterFrames(inBlock, sigma)
%IFILTERFRAMES Apply imgaussfilt to a Y x X x F stack, frame by frame.
%
% NaN pixels are zeroed for the filtering and restored afterwards. The mask
% is per element (per frame), so a pixel that is NaN in one frame does not
% affect any other frame, whatever the blocking of the frames.

spatialMask = isnan(inBlock);

if any(spatialMask(:))
    outBlock = inBlock;
    outBlock(spatialMask) = 0;
    outBlock = imgaussfilt(outBlock, sigma, 'FilterDomain', 'spatial');
    outBlock(spatialMask) = NaN;
else
    outBlock = imgaussfilt(inBlock, sigma, 'FilterDomain', 'spatial');
end
end

% =========================================================================
% Local helper: low-RAM .dat execution
% =========================================================================
function outFile = iSpatialGaussDatFile(inFile, SaveFolder, sigma, defaultOutput)
%ISPATIALGAUSSDATFILE Apply spatial filtering to a Y-X-T or Y-X-T-E .dat file.
%
% Every Y-X frame of the flattened trailing axes (T, E) is filtered, in
% blocks of whole frames, and the output keeps the input's axes.

slabIn = spatialSlabIO('open', inFile);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));
assertDatLayout(slabIn.Info, {{'Y','X','T'}, {'Y','X','T','E'}}, 'spatialGaussFilt');
Ny = slabIn.Ny;
Nx = slabIn.Nx;
Nt = slabIn.nFrames;   % frames across all trailing axes

% Write through a scratch file so the declared pipeline output only appears
% once the run has completed, and so the input can safely be the file the
% declared output would overwrite.
outFile = fullfile(SaveFolder, defaultOutput);
[~, outStem, outExt] = fileparts(defaultOutput);
tmpFile = fullfile(SaveFolder, [outStem, '_writing', outExt]);
slabOut = spatialSlabIO('create', tmpFile, ...
    datHeaderFromInfo(slabIn.Info, outStem, 'dataClass', 'single'));
cOut = onCleanup(@() spatialSlabIO('close', slabOut));

% Live arrays of a block: the frames, their NaN mask (a quarter of the
% frames), the zero-filled copy, and the filtered output.
frameBytes = Ny * Nx * getByteSize('single');
totalBytes = frameBytes * Nt;
nChunks = calculateMaxChunkSize(totalBytes, 4, 0.1);
chunkFrames = ceil(Nt / nChunks);
nChunks = ceil(Nt / chunkFrames);
fprintf('Spatial Gaussian filter: %d frame(s) in %d chunk(s)\n', Nt, nChunks);

for c = 1:nChunks
    tStart = (c-1) * chunkFrames + 1;
    tEnd   = min(tStart + chunkFrames - 1, Nt);

    fprintf('Chunk %i/%i [Reading frames %d-%d ...]\n', c, nChunks, tStart, tEnd)
    slab = single(spatialSlabIO('read', slabIn, 1:Nx, tStart:tEnd));

    fprintf('Chunk %i/%i [Filtering ...]\n', c, nChunks)
    slab = iFilterFrames(slab, sigma);

    fprintf('Chunk %i/%i [Writing to file ...]\n', c, nChunks)
    spatialSlabIO('write', slabOut, 1:Nx, slab, tStart:tEnd);
end

spatialSlabIO('finalize', slabOut);
clear cIn cOut; % close the input reader before the move below

[moveOk, moveMsg] = movefile(tmpFile, outFile, 'f');
assert(moveOk, 'spatialGaussFilt:OutputMoveFailed', ...
    'Failed to move "%s" onto "%s": %s', tmpFile, outFile, moveMsg);
fprintf('Finished spatial Gaussian filter.\n');
end
