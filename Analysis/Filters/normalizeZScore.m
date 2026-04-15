function outData = normalizeZScore(data, SaveFolder, varargin)
%NORMALIZEZSCORE Normalize continuous image time-series data to z-scores.
%
%   outData = normalizeZScore(data, SaveFolder)
%
%   This function normalizes image time-series data along the time
%   dimension using z-score normalization:
%
%       z = (x - mean(x, T)) ./ std(x, T)
%
%   Supported execution modes:
%       1) STANDARD MODE (in-memory)
%          - Triggered when "data" is a numeric array or a UMT struct
%       2) LOW-RAM MODE (file-backed)
%          - Triggered when "data" is a .dat filename
%
%   Accepted input forms:
%       1) Numeric array with dimensions Y x X x T
%       2) Filename to a .dat file storing Y x X x T data
%       3) UMT struct with image entries of dimensions Y x X x T
%       4) Filename to a .umt or .mat file containing one UMT struct
%
%   Input/output behavior:
%       - If the input is a numeric YXT array, the output is a numeric YXT
%         array with the same size.
%       - If the input is a .dat filename, the output is a .dat filename.
%       - If the input is a UMT struct, the output is a UMT struct.
%       - If the input is a .umt or .mat filename, the file is loaded in
%         RAM and the output is a UMT struct.
%
%   Limitations:
%       - Only continuous image time-series data are supported.
%       - Event-split data are not supported.
%       - UMT entries must use dimNames {'Y','X','T'} exactly.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder used for file resolution and AcqInfos.mat lookup.
%
%   Output:
%       outData    - Z-score normalized output with the same representation
%                    type as the input.
%
%   Notes:
%       - In low-RAM mode, the input file is processed chunk-by-chunk along
%         the X dimension.
%       - For constant time traces, standard deviation is forced to 1 to
%         avoid division by zero, preserving the original algorithm.
%       - Raw .dat files are assumed to store YXT data in single precision.
%       - For UMT/.umt input, RAM-safe mode is not available.

default_Output = 'normZ.dat'; %#ok<NASGU>

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
% Case 1: In-memory numeric array
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, ...
        mfilename, 'data');
    outData = iZscoreArray(data);
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
    ext = lower(ext);

    switch ext
        case '.dat'
            outData = iZscoreDatFile(dataFile, SaveFolder);
            return

        case '.umt'
            warning('normalizeZScore:UMTFileLoadsInRAM', ...
                ['RAM-safe mode is not available for data stored in this format. ' ...
                 'Loading the UMT content into RAM.']);
            data = loadData(dataFile);

        otherwise
            error('normalizeZScore:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

% -------------------------------------------------------------------------
% Case 3: UMT struct
% -------------------------------------------------------------------------
if ~isstruct(data)
    error('normalizeZScore:UnsupportedInputType', ...
        ['Input "data" must be a YXT array, a .dat filename, ' ...
         'a UMT struct, or a .umt/.mat filename containing a UMT struct.']);
end

[entryNames, entryData, entryDims, labels] = iExtractValidUMTData(data);

outStruct = data;
outStruct.data = struct();

for iEntry = 1:numel(entryNames)
    outStruct.data.(entryNames{iEntry}) = struct( ...
        'value', iZscoreArray(entryData{iEntry}), ...
        'dimNames', {entryDims{iEntry}});
end

if ~isempty(fieldnames(labels))
    outStruct.labels = labels;
elseif isfield(outStruct, 'labels')
    outStruct = rmfield(outStruct, 'labels');
end

if isfield(outStruct, 'eventInfo')
    outStruct = rmfield(outStruct, 'eventInfo');
end

validateUMTStruct(outStruct, 'requireEventInfo', true);
outData = outStruct;

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            ['Normalize continuous image time-series data to zero mean and ' ...
             'unit standard deviation along the time dimension.']);

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT array, .dat filename, ' ...
             'UMT struct, or .umt/.mat file containing one UMT struct.'], ...
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
            'Folder used for file resolution and AcqInfos.mat lookup.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            'data', ...
            'Z-score normalized output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: z-score one in-memory YXT array
% =========================================================================
function outArray = iZscoreArray(inArray)
%IZSCOREARRAY Apply z-score normalization along the 3rd dimension.

origSz = size(inArray);
workArray = reshape(inArray, [], origSz(3));

mu  = mean(workArray, 2, 'omitnan');
sig = std(workArray, 0, 2, 'omitnan');
sig(sig == 0) = 1;

outArray = reshape((workArray - mu) ./ sig, origSz);
end

% =========================================================================
% Helper: low-RAM .dat execution
% =========================================================================
function outFile = iZscoreDatFile(inFile, SaveFolder)
%IZSCOREDATFILE Apply z-score normalization to a raw YXT .dat file.

[Ny, Nx, Nt] = iGetRawDatInfo(SaveFolder, inFile);

outFile = fullfile(fileparts(inFile), 'normZ.dat');
preallocateDatFile(outFile, [Ny, Nx, Nt], 'single');

fidIn  = fopen(inFile,  'r');
assert(fidIn ~= -1, 'normalizeZScore:OpenInputFailed', ...
    'Failed to open input file "%s".', inFile);
cIn = onCleanup(@() safeFclose(fidIn)); %#ok<NASGU>

fidOut = fopen(outFile, 'r+');
assert(fidOut ~= -1, 'normalizeZScore:OpenOutputFailed', ...
    'Failed to open output file "%s".', outFile);
cOut = onCleanup(@() safeFclose(fidOut)); %#ok<NASGU>

totalBytes = Ny * Nx * Nt * getByteSize('single');
nChunks = calculateMaxChunkSize(totalBytes, 2);
chunkX = ceil(Nx / nChunks);

for c = 1:nChunks
    xStart = (c-1) * chunkX + 1;
    xEnd   = min(xStart + chunkX - 1, Nx);
    xIdx   = xStart:xEnd;

    slab = spatialSlabIO('read', fidIn, Ny, Nx, Nt, xIdx, 'single');

    slab2D = reshape(slab, [], Nt);

    mu  = mean(slab2D, 2, 'omitnan');
    sig = std(slab2D, 0, 2, 'omitnan');
    sig(sig == 0) = 1;

    slab2D = (slab2D - mu) ./ sig;
    slab = reshape(single(slab2D), size(slab));

    spatialSlabIO('write', fidOut, Ny, Nx, Nt, xIdx, 'single', slab);
end

fclose(fidIn);
fclose(fidOut);
end

% =========================================================================
% Helper: extract raw .dat dimensions from AcqInfos.mat / file size
% =========================================================================
function [Ny, Nx, Nt] = iGetRawDatInfo(SaveFolder, inFile)
%IGETRAWDATINFO Return YXT dimensions for a raw .dat file.

acqFile = fullfile(SaveFolder, 'AcqInfos.mat');
if ~isfile(acqFile)
    error('normalizeZScore:MissingAcqInfos', ...
        'AcqInfos.mat was not found in SaveFolder "%s".', SaveFolder);
end

S = load(acqFile);
if isfield(S, 'AcqInfoStream')
    acqInfo = S.AcqInfoStream;
else
    fn = fieldnames(S);
    acqInfo = S.(fn{1});
end

if ~isfield(acqInfo, 'Height') || ~isfield(acqInfo, 'Width')
    error('normalizeZScore:InvalidAcqInfos', ...
        'AcqInfoStream must contain Height and Width.');
end

Ny = double(acqInfo.Height);
Nx = double(acqInfo.Width);

if isfield(acqInfo, 'Length') && ~isempty(acqInfo.Length)
    Nt = double(acqInfo.Length);
else
    fileInfo = dir(inFile);
    nElem = fileInfo.bytes / getByteSize('single');
    if mod(nElem, Ny * Nx) ~= 0
        error('normalizeZScore:InvalidRawDatLength', ...
            'File size is incompatible with YXT dimensions for "%s".', inFile);
    end
    Nt = nElem / (Ny * Nx);
end
end

% =========================================================================
% Helper: validate/extract UMT data
% =========================================================================
function [entryNames, entryData, entryDims, labels] = iExtractValidUMTData(umt)
%IEXTRACTVALIDUMTDATA Validate and extract continuous YXT image entries.

validateUMTStruct(umt, 'requireEventInfo', false);

if ~strcmpi(umt.kind, 'image')
    error('normalizeZScore:InvalidUMTKind', ...
        ['Operation aborted. UMT input must have kind = "image". ' ...
         'This function does not support non-image UMT structures.']);
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error('normalizeZScore:EmptyUMTData', ...
        'Operation aborted. UMT data is empty.');
end

entryData = cell(size(entryNames));
entryDims = cell(size(entryNames));

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    if ~isequal(thisDims, {'Y','X','T'})
        error('normalizeZScore:InvalidUMTEntry', ...
            ['Operation aborted. All UMT entries must use dimNames ' ...
             '{''Y'',''X'',''T''}. Invalid entry: "%s".'], ...
            entryNames{iEntry});
    end

    entryData{iEntry} = thisEntry.value;
    entryDims{iEntry} = thisDims;
end

if isfield(umt, 'eventInfo')
    error('normalizeZScore:EventDataNotSupported', ...
        'normalizeZScore does not support event-split data.');
end

if isfield(umt, 'labels')
    labels = umt.labels;
else
    labels = struct();
end
end