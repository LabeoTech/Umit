function [outData, metaData] = SpeckleTimeLapse(RawFolder, varargin)
% SPECKLETIMELAPSE Create time-lapse Laser Speckle Contrast Images from raw binary files.
%
%   [outData, metaData] = SPECKLETIMELAPSE(RawFolder)
%   reads raw UMIT/LabeoTech binary image files from RawFolder, estimates
%   the total recording duration, samples temporal windows at regular time
%   intervals, and calculates one Laser Speckle Contrast Image per sampled
%   time point.
%
%   [outData, metaData] = SPECKLETIMELAPSE(RawFolder, opts)
%   uses the options specified in opts.
%
%   Inputs:
%       RawFolder:
%           Folder containing raw acquisition files named img_*.bin and the
%           associated Info.txt file readable by ReadInfoFile.
%
%       opts:
%           Optional structure with the following fields:
%
%           frameWindowSize:
%               Number of frames used for each LSCI map. The function uses
%               floor(frameWindowSize/2) frames before the reference frame
%               and frameWindowSize - floor(frameWindowSize/2) - 1 frames
%               after it.
%               Default: 40.
%
%           timeStepMin:
%               Time interval, in minutes, between sampled reference frames.
%               Default: 30.
%
%           sType:
%               Speckle contrast algorithm type. Valid values are:
%                   'Spatial'  - spatial standard deviation using a 5x5 area.
%                   'Temporal' - temporal standard deviation using 5 frames.
%               Default: 'Temporal'.
%
%           bLogScale:
%               Apply -log10 scaling to each LSCI map.
%               Default: false.
%
%   Outputs:
%       outData:
%           Single-precision Y-by-X-by-T array containing one LSCI map per
%           sampled time point.
%
%       metaData:
%           Metadata structure describing the output data. The output frame
%           rate is intentionally set to 2 Hz as a dummy visualization frame
%           rate.
%
%   Assumptions:
%       - Raw files are named img_00000.bin, img_00001.bin, etc.
%       - The img_*.bin index series starts at zero and is consecutive.
%       - The acquisition contains a single camera.
%       - The Info.txt file contains exactly one Illumination field.

% Defaults:
default_Output = 'speckle_timelapse.dat'; %#ok<NASGU> % PipelineManager reference.
default_opts = struct('frameWindowSize',40,'timeStepMin',30,'sType','Temporal','bLogScale',false);
opts_values = struct('frameWindowSize',[1,Inf],'timeStepMin',[eps,Inf],'sType',{{'Spatial','Temporal'}},'bLogScale',[true,false]); %#ok<NASGU> % This is here only as a reference for PIPELINEMANAGER.m.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'RawFolder',@isfolder);
addOptional(p,'opts',default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p,RawFolder,varargin{:});

RawFolder = p.Results.RawFolder;
opts = p.Results.opts;
clear p

opts = fillMissingOptions(opts, default_opts);
opts = validateOptions(opts);

if ~strcmp(RawFolder(end), filesep)
    RawFolder = strcat(RawFolder, filesep);
end

% Validate required UMIT helper functions.
assert(exist('ReadInfoFile', 'file') == 2, ...
    'umIToolbox:SpeckleTimeLapse:MissingFunction', ...
    'ReadInfoFile was not found on the MATLAB path.');

% Read acquisition information from Info.txt.
AcqInfoStream = ReadInfoFile(RawFolder);

assert(isfield(AcqInfoStream, 'FrameRateHz'), ...
    'umIToolbox:SpeckleTimeLapse:MissingInfoField', ...
    'FrameRateHz was not found in the acquisition information.');

if ~isfield(AcqInfoStream, 'Camera_Model')
    % Backward compatibility with older acquisitions.
    AcqInfoStream.Camera_Model = 'D1024';
end

% Validate the main assumptions for this function.
validateSingleIlluminationField(AcqInfoStream);

if isfield(AcqInfoStream, 'MultiCam') && AcqInfoStream.MultiCam
    error('umIToolbox:SpeckleTimeLapse:UnsupportedMultiCamera', ...
        ['SpeckleTimeLapse currently supports single-camera acquisitions only. ' ...
         'The Info.txt file indicates MultiCam = true.']);
end

% Raw binary image file layout. This matches ImagesClassification.
hWima = 5;
imgFilesList = dir([RawFolder 'img_*.bin']);

assert(~isempty(imgFilesList), ...
    'umIToolbox:SpeckleTimeLapse:MissingRawFiles', ...
    'No raw image files matching img_*.bin were found in %s.', RawFolder);

[imgFilesList, imgFileIdx] = sortAndValidateImageFileList(imgFilesList);

% Read image dimensions from the first raw file header.
header = memmapfile([RawFolder imgFilesList(1).name], ...
    'Offset', 0, ...
    'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, ...
    'Repeat', 1);

nx = double(header.Data.header(2));
ny = double(header.Data.header(3));

assert(nx > 0 && ny > 0, ...
    'umIToolbox:SpeckleTimeLapse:InvalidImageSize', ...
    'Invalid image size detected from raw file header.');

bytesPerRawFrame = nx * ny * 2 + 3 * 8;
bytesFileHeader = hWima * 4;

% Count raw frames in each img_*.bin file.
nFramesPerFile = zeros(numel(imgFilesList), 1);

for indFile = 1:numel(imgFilesList)
    currentFile = [RawFolder imgFilesList(indFile).name];
    fileInfo = dir(currentFile);

    dataBytes = fileInfo.bytes - bytesFileHeader;

    assert(dataBytes >= 0 && mod(dataBytes, bytesPerRawFrame) == 0, ...
        'umIToolbox:SpeckleTimeLapse:InvalidRawFileSize', ...
        'The file %s has an unexpected size and cannot be parsed as raw image data.', ...
        imgFilesList(indFile).name);

    nFramesPerFile(indFile) = dataBytes / bytesPerRawFrame;
end

totalRawFrames = sum(nFramesPerFile);
totalTimeMin = totalRawFrames / AcqInfoStream.FrameRateHz / 60;

assert(totalRawFrames > 0, ...
    'umIToolbox:SpeckleTimeLapse:EmptyRecording', ...
    'No raw frames were found in %s.', RawFolder);

% In this function, there is only one illumination channel. Therefore, raw
% frame indices and imaging frame indices are equivalent.
totalFrames = totalRawFrames;
frameRateHz = AcqInfoStream.FrameRateHz;
 
assert(opts.frameWindowSize<totalFrames,'umIToolbox:SpeckleTimeLapse:InvalidInput',...
    'The frameWindowSize should not be greater than %i',totalFrames);
% Determine sampled time points and reference frame indices.
timePointsMin = opts.timeStepMin:opts.timeStepMin:totalTimeMin;
centerFrames = round(timePointsMin * 60 * frameRateHz) + 1;

nBefore = floor(opts.frameWindowSize / 2);
nAfter = opts.frameWindowSize - nBefore - 1;

validIdx = centerFrames - nBefore >= 1 & centerFrames + nAfter <= totalFrames;

timePointsMin = timePointsMin(validIdx);
centerFrames = centerFrames(validIdx);

assert(~isempty(centerFrames), ...
    'umIToolbox:SpeckleTimeLapse:NoValidTimePoints', ...
    ['No valid time points were found. Reduce frameWindowSize or ' ...
     'timeStepMin, or use a longer recording.']);

% Preallocate output as Y-by-X-by-time.
outData = zeros(ny, nx, numel(centerFrames), 'single');

fprintf('SpeckleTimeLapse: found %d raw files from img_%05d.bin to img_%05d.bin.\n', ...
    numel(imgFilesList), imgFileIdx(1), imgFileIdx(end));
fprintf('SpeckleTimeLapse: total recording time = %.2f min.\n', totalTimeMin);
fprintf('SpeckleTimeLapse: processing %d time points using %d frames per map.\n', ...
    numel(centerFrames), opts.frameWindowSize);

cumFramesPerFile = cumsum(nFramesPerFile);

% Calculate one LSCI map per sampled time point.
for indTime = 1:numel(centerFrames)
    frameWindow = (centerFrames(indTime) - nBefore):(centerFrames(indTime) + nAfter);

    [frameBlock, usedBinFiles] = readRawFrameBlock( ...
        RawFolder, ...
        imgFilesList, ...
        frameWindow, ...
        cumFramesPerFile, ...
        nFramesPerFile, ...
        bytesFileHeader, ...
        bytesPerRawFrame, ...
        ny, ...
        nx);

    % Calculate the LSCI map directly from the raw frame block. This avoids
    % the SpeckleMapping dependency because this function does not operate
    % on pre-classified .dat/.mat pairs.
    outData(:, :, indTime) = calculateSpeckleMap(frameBlock, opts.sType, opts.bLogScale);

    fprintf('[%3.0f%%] Time point %.2f min processed\n', ...
        100 * indTime / numel(centerFrames), ...
        timePointsMin(indTime));
end

% Create metadata for PipelineManager/save2Dat compatibility.
metaData = struct();
metaData.Datatype = 'single';
metaData.datName = 'data';
metaData.datFile = 'speckle_timelapse.dat';
metaData.datSize = [ny, nx];
metaData.datLength = size(outData, 3);
metaData.dim_names = {'Y', 'X', 'T'};
metaData.FirstDim = 'y';

% Use a dummy visualization frame rate, as requested.
metaData.Freq = 2;
metaData.FrameRateHz = 2;

% Forward-compatible convenience fields.
metaData.Width = nx;
metaData.Height = ny;
metaData.Length = size(outData, 3);

% Store useful provenance fields.
metaData.RawFolder = RawFolder;
metaData.SourceFrameRateHz = AcqInfoStream.FrameRateHz;
metaData.TotalRecordingTimeMin = totalTimeMin;
metaData.TimePointsMin = timePointsMin;
metaData.CenterFrames = centerFrames;
metaData.FrameWindowSize = opts.frameWindowSize;
metaData.SpeckleAlgorithm = opts.sType;
metaData.bLogScale = opts.bLogScale;
metaData.SourceFileFirstIndex = imgFileIdx(1);
metaData.SourceFileLastIndex = imgFileIdx(end);
metaData.SourceFileCount = numel(imgFilesList);

% If genMetaData exists, use it to align with the current UMIT metadata
% format while preserving the fields above as extra metadata.
if exist('genMetaData', 'file') == 2
    try
        metaData = genMetaData(outData, {'Y', 'X', 'T'}, metaData);
        metaData.Freq = 2;
        metaData.FrameRateHz = 2;
        metaData.Width = nx;
        metaData.Height = ny;
        metaData.Length = size(outData, 3);
        metaData.TotalRecordingTimeMin = totalTimeMin;
        metaData.TimePointsMin = timePointsMin;
        metaData.CenterFrames = centerFrames;
        metaData.FrameWindowSize = opts.frameWindowSize;
        metaData.SpeckleAlgorithm = opts.sType;
        metaData.bLogScale = opts.bLogScale;
        metaData.SourceFileFirstIndex = imgFileIdx(1);
        metaData.SourceFileLastIndex = imgFileIdx(end);
        metaData.SourceFileCount = numel(imgFilesList);
    catch ME
        warning('umIToolbox:SpeckleTimeLapse:GenMetaDataFailed', ...
            ['genMetaData failed with message:\n%s\n' ...
             'Using locally generated metadata instead.'], ME.message);
    end
end

fprintf('SpeckleTimeLapse: finished.\n');

end

%% ========================================================================
%  Local helper functions
%  ========================================================================

function opts = fillMissingOptions(opts, default_opts)
% FILLMISSINGOPTIONS Add missing option fields using defaults.

defaultFields = fieldnames(default_opts);

for indField = 1:numel(defaultFields)
    currentField = defaultFields{indField};

    if ~isfield(opts, currentField) || isempty(opts.(currentField))
        opts.(currentField) = default_opts.(currentField);
    end
end

end

function opts = validateOptions(opts)
% VALIDATEOPTIONS Validate SpeckleTimeLapse option values.

validateattributes(opts.frameWindowSize, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1}, mfilename, 'opts.frameWindowSize');

validateattributes(opts.timeStepMin, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0}, mfilename, 'opts.timeStepMin');

assert(ischar(opts.sType) || isstring(opts.sType), ...
    'umIToolbox:SpeckleTimeLapse:InvalidSType', ...
    'opts.sType must be either ''Spatial'' or ''Temporal''.');

opts.sType = validatestring(char(opts.sType), {'Spatial', 'Temporal'}, ...
    mfilename, 'opts.sType');

assert(islogical(opts.bLogScale) || isnumeric(opts.bLogScale), ...
    'umIToolbox:SpeckleTimeLapse:InvalidOption', ...
    'opts.bLogScale must be true or false.');

opts.bLogScale = logical(opts.bLogScale);

end

function validateSingleIlluminationField(AcqInfoStream)
% VALIDATESINGLEILLUMINATIONFIELD Confirm that the acquisition has one color only.

Tags = fieldnames(AcqInfoStream);
idxIllumination = contains(Tags, 'Illumination');
nIlluminations = sum(idxIllumination);

if nIlluminations ~= 1
    error('umIToolbox:SpeckleTimeLapse:InvalidIlluminationCount', ...
        ['SpeckleTimeLapse currently supports single-color acquisitions only. ' ...
         'Info.txt contains %d Illumination fields.'], nIlluminations);
end

end

function [imgFilesList, imgFileIdx] = sortAndValidateImageFileList(imgFilesList)
% SORTANDVALIDATEIMAGEFILELIST Sort and validate raw image file sequence.
%
%   The expected file naming convention is:
%       img_00000.bin
%       img_00001.bin
%       img_00002.bin
%       ...
%
%   The index series must start at zero and must be consecutive.

fileNames = {imgFilesList.name};
imgFileIdx = nan(size(fileNames));

for indFile = 1:numel(fileNames)
    token = regexp(fileNames{indFile}, '^img_(\d+)\.bin$', 'tokens', 'once');

    if isempty(token)
        error('umIToolbox:SpeckleTimeLapse:InvalidRawFileName', ...
            ['Invalid raw image file name detected: %s. ' ...
             'Expected names such as img_00000.bin, img_00001.bin, etc.'], ...
            fileNames{indFile});
    end

    imgFileIdx(indFile) = str2double(token{1});
end

[imgFileIdx, sortIdx] = sort(imgFileIdx);
imgFilesList = imgFilesList(sortIdx);

if imgFileIdx(1) ~= 0
    error('umIToolbox:SpeckleTimeLapse:InvalidRawFileSequence', ...
        ['The img_*.bin index series must start at zero. ' ...
         'First detected index was %d.'], imgFileIdx(1));
end

expectedIdx = 0:(numel(imgFileIdx) - 1);

if ~isequal(imgFileIdx, expectedIdx)
    missingIdx = setdiff(expectedIdx, imgFileIdx);

    if isempty(missingIdx)
        error('umIToolbox:SpeckleTimeLapse:InvalidRawFileSequence', ...
            'The img_*.bin index series is not consecutive.');
    else
        error('umIToolbox:SpeckleTimeLapse:MissingRawFiles', ...
            'The img_*.bin index series is missing at least index %d.', ...
            missingIdx(1));
    end
end

end

function [frameBlock, usedBinFiles] = readRawFrameBlock( ...
    RawFolder, ...
    imgFilesList, ...
    frameWindow, ...
    cumFramesPerFile, ...
    nFramesPerFile, ...
    bytesFileHeader, ...
    bytesPerRawFrame, ...
    ny, ...
    nx)
% READRAWFRAMEBLOCK Read selected global raw-frame indices from img_*.bin files.
%
%   [frameBlock, usedBinFiles] = READRAWFRAMEBLOCK(...) reads a set of
%   global frame indices directly from the raw binary files using grouped
%   fseek/fread operations. The output is returned as a Y-by-X-by-T
%   single-precision array.
%
%   This implementation is optimized for sparse time-lapse extraction from
%   long recordings. Requested frames are first mapped to their source files.
%   Each required file is then opened once, and all requested local frames
%   from that file are read by direct byte offset.
%
%   Outputs:
%       frameBlock:
%           Y-by-X-by-T single-precision array containing the requested
%           frames in chronological order.
%
%       usedBinFiles:
%           Cell array containing the img_*.bin filenames used to build this
%           frame block.

frameBlock = zeros(ny, nx, numel(frameWindow), 'single');

fileIdxPerFrame = zeros(size(frameWindow));
localIdxPerFrame = zeros(size(frameWindow));

% Map global frame indices to source files and local frame indices.
for indFrame = 1:numel(frameWindow)
    globalFrameIdx = frameWindow(indFrame);

    fileIdx = find(globalFrameIdx <= cumFramesPerFile, 1, 'first');

    assert(~isempty(fileIdx), ...
        'umIToolbox:SpeckleTimeLapse:FrameIndexOutOfRange', ...
        'Requested raw frame index %d exceeds the available raw frames.', globalFrameIdx);

    if fileIdx == 1
        localFrameIdx = globalFrameIdx;
    else
        localFrameIdx = globalFrameIdx - cumFramesPerFile(fileIdx - 1);
    end

    assert(localFrameIdx >= 1 && localFrameIdx <= nFramesPerFile(fileIdx), ...
        'umIToolbox:SpeckleTimeLapse:FrameIndexOutOfRange', ...
        'Invalid local frame index while reading raw image data.');

    fileIdxPerFrame(indFrame) = fileIdx;
    localIdxPerFrame(indFrame) = localFrameIdx;
end

neededFiles = unique(fileIdxPerFrame, 'stable');
usedBinFiles = {imgFilesList(neededFiles).name};

% Read all requested frames, opening each source file only once.
for indNeededFile = 1:numel(neededFiles)
    fileIdx = neededFiles(indNeededFile);
    currentFile = [RawFolder imgFilesList(fileIdx).name];

    fid = fopen(currentFile, 'r');

    if fid < 0
        error('umIToolbox:SpeckleTimeLapse:FileOpenFailed', ...
            'Failed to open raw image file: %s.', currentFile);
    end

    cleanupObj = onCleanup(@() fclose(fid));

    framePositions = find(fileIdxPerFrame == fileIdx);

    for indPosition = 1:numel(framePositions)
        outFrameIdx = framePositions(indPosition);
        localFrameIdx = localIdxPerFrame(outFrameIdx);

        % Each frame record contains:
        %   3 uint64 header values, followed by nx*ny uint16 pixels.
        % Skip the file header and the frame-level header.
        imageOffsetBytes = bytesFileHeader + ...
            (localFrameIdx - 1) * bytesPerRawFrame + ...
            3 * 8;

        fseekStatus = fseek(fid, imageOffsetBytes, 'bof');

        if fseekStatus ~= 0
            error('umIToolbox:SpeckleTimeLapse:FseekFailed', ...
                ['Failed to seek to frame %d in file %s. ' ...
                 'Requested byte offset was %d.'], ...
                localFrameIdx, currentFile, imageOffsetBytes);
        end

        currentImage = fread(fid, nx * ny, 'uint16=>single');

        if numel(currentImage) ~= nx * ny
            error('umIToolbox:SpeckleTimeLapse:UnexpectedEndOfFile', ...
                ['Unexpected end of file while reading frame %d from %s. ' ...
                 'Expected %d pixels, read %d.'], ...
                localFrameIdx, currentFile, nx * ny, numel(currentImage));
        end

        % Raw image is stored X-by-Y. Convert to UMIT convention Y-by-X.
        frameBlock(:, :, outFrameIdx) = reshape(currentImage, nx, ny).';
    end

    clear cleanupObj
end

end

function speckleMap = calculateSpeckleMap(frameBlock, sType, bLogScale)
% CALCULATESPECKLEMAP Calculate one Laser Speckle Contrast Image.
%
%   speckleMap = CALCULATESPECKLEMAP(frameBlock, sType, bLogScale)
%   calculates a single LSCI map from a Y-by-X-by-T raw frame block.
%
%   The algorithm follows the standard non-chunked SpeckleMapping logic:
%       1) Normalize each pixel by its temporal mean.
%       2) Apply stdfilt using either a spatial or temporal kernel.
%       3) Average the filtered stack over time.
%       4) Optionally apply -log10 scaling.
%
%   Inputs:
%       frameBlock:
%           Y-by-X-by-T single-precision raw image block.
%
%       sType:
%           'Spatial' or 'Temporal'.
%
%       bLogScale:
%           Boolean flag to apply -log10 scaling.
%
%   Outputs:
%       speckleMap:
%           Y-by-X single-precision speckle contrast map.

frameBlock = single(frameBlock);

% Normalize each pixel by its temporal mean. This preserves the behavior of
% SpeckleMapping and reduces intensity-offset effects across the frame block.
temporalMean = mean(frameBlock, 3, 'omitnan');
frameBlock = frameBlock ./ temporalMean;

switch lower(sType)
    case 'spatial'
        kernel = zeros(5, 5, 1, 'double');
        kernel(:, :, 1) = double(fspecial('disk', 2) > 0);

    case 'temporal'
        kernel = ones(1, 1, 5, 'double');

    otherwise
        error('umIToolbox:SpeckleTimeLapse:InvalidSType', ...
            'sType must be either ''Spatial'' or ''Temporal''.');
end

% Standard non-chunked SpeckleMapping algorithm.
speckleMap = single(stdfilt(frameBlock, kernel));
speckleMap = mean(speckleMap, 3, 'omitnan');

if bLogScale
    speckleMap = -log10(speckleMap);
end

speckleMap = single(speckleMap);

end