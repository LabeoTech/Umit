function [referenceImage, sourceInfo] = readImageReferenceSource(sourceFile, opts)
%READIMAGEREFERENCESOURCE Read an external file as a 2D ImageReference source.
%
%   referenceImage = readImageReferenceSource(sourceFile)
%   [referenceImage, sourceInfo] = readImageReferenceSource(sourceFile)
%   [referenceImage, sourceInfo] = readImageReferenceSource(sourceFile, ...
%       'maxFramesToAverage', 100)
%
%   Reads a supported external file and returns a 2D single image suitable
%   for createImageReference.
%
%   Supported sources:
%       .tif/.tiff - Single image or multi-page TIFF stack.
%       .dat       - UMIT DAT image time series using DatImageSource.
%       .umt       - UMIT image source using UMTImageSource.
%       .mat       - MAT-file variable selected from 2D numeric arrays.
%
%   Behavior:
%       - TIF stack: average up to maxFramesToAverage evenly spaced pages.
%       - DAT/UMT: average up to maxFramesToAverage evenly spaced frames.
%       - MAT: list 2D numeric variables and load the selected one.
%       - RGB TIFF pages are converted to grayscale.
%       - Output is always a 2D single array.
%
%   Inputs:
%       sourceFile - Path to supported source file.
%
%   Name-Value options:
%       maxFramesToAverage - Maximum number of pages/frames to average.
%                            Default: 100.
%
%   Outputs:
%       referenceImage - 2D single reference image.
%       sourceInfo     - Struct with source metadata.
%
%   Notes:
%       This function does not normalize the image intensity range. It keeps
%       the source numeric scale after conversion to single.
%
%       MAT-file support intentionally lists only 2D numeric arrays. This
%       avoids silently using higher-dimensional data or non-image variables.

arguments
    sourceFile (1,1) string
    opts.maxFramesToAverage (1,1) double {mustBePositive, mustBeInteger} = 100
end

sourceFile = strtrim(sourceFile);

if sourceFile == "" || ~isfile(sourceFile)
    error('readImageReferenceSource:FileNotFound', ...
        'Reference source file not found: %s', sourceFile);
end

[~, ~, ext] = fileparts(sourceFile);
ext = lower(ext);

switch ext
    case {'.tif', '.tiff'}
        [referenceImage, sourceInfo] = iReadTiffSource(sourceFile, opts.maxFramesToAverage);

    case '.dat'
        [referenceImage, sourceInfo] = iReadDatSource(sourceFile, opts.maxFramesToAverage);

    case '.umt'
        [referenceImage, sourceInfo] = iReadUmtSource(sourceFile, opts.maxFramesToAverage);

    case '.mat'
        [referenceImage, sourceInfo] = iReadMatSource(sourceFile);

    otherwise
        error('readImageReferenceSource:UnsupportedFileType', ...
            'Unsupported reference source type: %s', ext);
end

referenceImage = iValidateReferenceImage(referenceImage);

end

% =========================================================================
% Source readers
% =========================================================================

function [referenceImage, sourceInfo] = iReadTiffSource(tiffFile, maxFramesToAverage)
%IREADTIFFSOURCE Read a TIF/TIFF file as one reference image.

info = imfinfo(tiffFile);
nPages = numel(info);

if nPages < 1
    error('readImageReferenceSource:EmptyTiff', ...
        'TIF file contains no readable pages: %s', tiffFile);
end

nRead = min(nPages, maxFramesToAverage);

if nPages > nRead
    frameIdx = unique(round(linspace(1, nPages, nRead)));
else
    frameIdx = 1:nPages;
end

referenceImage = iAverageIndexedFrames( ...
    frameIdx, ...
    @(idx) iReadTiffPage(tiffFile, idx, info));

[sourceFolder, sourceName, sourceExt] = fileparts(tiffFile);

sourceInfo = iBaseSourceInfo(tiffFile, 'external_tiff', size(referenceImage));
sourceInfo.nSourceFrames = nPages;
sourceInfo.usedFrameIndex = frameIdx(:).';
sourceInfo.selectedEntry = '';
sourceInfo.selectedVariable = '';
sourceInfo.sourceFrame = iScalarSourceFrame(frameIdx);
sourceInfo.processingNotes = sprintf( ...
    'Created from TIF file "%s%s". Averaged %d of %d page(s).', ...
    sourceName, sourceExt, numel(frameIdx), nPages);
sourceInfo.sourceFolder = sourceFolder;

end

function [referenceImage, sourceInfo] = iReadDatSource(datFile, maxFramesToAverage)
%IREADDATSOURCE Read a DAT image time series as one reference image.

if exist('DatImageSource', 'class') ~= 8
    error('readImageReferenceSource:MissingDatImageSource', ...
        'DatImageSource is not available on the MATLAB path.');
end

src = DatImageSource(char(datFile));
[referenceImage, usedFrameIndex, nSourceFrames] = iAverageImageSource(src, maxFramesToAverage);

sourceInfo = iBaseSourceInfo(datFile, 'external_dat', size(referenceImage));
sourceInfo.nSourceFrames = nSourceFrames;
sourceInfo.usedFrameIndex = usedFrameIndex;
sourceInfo.selectedEntry = '';
sourceInfo.selectedVariable = '';
sourceInfo.sourceFrame = iScalarSourceFrame(usedFrameIndex);
sourceInfo.processingNotes = sprintf( ...
    'Created from DAT file. Averaged %d of %d frame(s).', ...
    numel(usedFrameIndex), nSourceFrames);

end

function [referenceImage, sourceInfo] = iReadUmtSource(umtFile, maxFramesToAverage)
%IREADUMTSOURCE Read a UMT image source as one reference image.

if exist('UMTImageSource', 'class') ~= 8
    error('readImageReferenceSource:MissingUMTImageSource', ...
        'UMTImageSource is not available on the MATLAB path.');
end

entrySummary = UMTImageSource.inspectEntries(char(umtFile));

if isempty(entrySummary) || height(entrySummary) == 0
    error('readImageReferenceSource:EmptyUmt', ...
        'No readable UMT entries were found.');
end

if height(entrySummary) == 1
    selectedEntry = char(entrySummary.Name(1));
else
    selectedEntry = iSelectUMTEntry(entrySummary);
end

src = UMTImageSource(char(umtFile), selectedEntry);
[referenceImage, usedFrameIndex, nSourceFrames] = iAverageImageSource(src, maxFramesToAverage);

sourceInfo = iBaseSourceInfo(umtFile, 'external_umt', size(referenceImage));
sourceInfo.nSourceFrames = nSourceFrames;
sourceInfo.usedFrameIndex = usedFrameIndex;
sourceInfo.selectedEntry = selectedEntry;
sourceInfo.selectedVariable = '';
sourceInfo.sourceFrame = iScalarSourceFrame(usedFrameIndex);
sourceInfo.processingNotes = sprintf( ...
    'Created from UMT entry "%s". Averaged %d of %d frame(s).', ...
    selectedEntry, numel(usedFrameIndex), nSourceFrames);

end

function [referenceImage, sourceInfo] = iReadMatSource(matFile)
%IREADMATSOURCE Read selected 2D numeric variable from a MAT file.

vars = whos('-file', matFile);

if isempty(vars)
    error('readImageReferenceSource:EmptyMatFile', ...
        'MAT file contains no variables.');
end

candidateMask = arrayfun(@iIs2DNumericMatVariable, vars);
candidateVars = vars(candidateMask);

if isempty(candidateVars)
    error('readImageReferenceSource:No2DNumericVariables', ...
        'MAT file contains no 2D numeric array variables.');
end

if numel(candidateVars) == 1
    selectedVariable = candidateVars(1).name;
else
    selectedVariable = iSelectMatVariable(candidateVars);
end

S = load(matFile, selectedVariable);
referenceImage = S.(selectedVariable);

if issparse(referenceImage)
    referenceImage = full(referenceImage);
end

referenceImage = iValidateReferenceImage(referenceImage);

sourceInfo = iBaseSourceInfo(matFile, 'external_mat', size(referenceImage));
sourceInfo.nSourceFrames = 1;
sourceInfo.usedFrameIndex = 1;
sourceInfo.selectedEntry = '';
sourceInfo.selectedVariable = selectedVariable;
sourceInfo.sourceFrame = [];
sourceInfo.processingNotes = sprintf( ...
    'Created from MAT file variable "%s".', selectedVariable);

end

% =========================================================================
% Frame averaging helpers
% =========================================================================

function [referenceImage, usedFrameIndex, nSourceFrames] = iAverageImageSource(src, maxFramesToAverage)
%IAVERAGEIMAGESOURCE Average frames from DatImageSource/UMTImageSource.

sourceSize = iGetImageSourceSize(src);

if numel(sourceSize) < 2
    error('readImageReferenceSource:InvalidSourceSize', ...
        'Image source size must include at least Y and X dimensions.');
end

if numel(sourceSize) < 3
    referenceImage = iValidateReferenceImage(src.getFrame(1));
    usedFrameIndex = 1;
    nSourceFrames = 1;
    return
end

nT = sourceSize(3);

if numel(sourceSize) >= 4
    nE = prod(sourceSize(4:end));
else
    nE = 1;
end

nSourceFrames = nT * nE;
nRead = min(nSourceFrames, maxFramesToAverage);

if nSourceFrames > nRead
    usedFrameIndex = unique(round(linspace(1, nSourceFrames, nRead)));
else
    usedFrameIndex = 1:nSourceFrames;
end

referenceImage = iAverageIndexedFrames( ...
    usedFrameIndex, ...
    @(idx) iReadImageSourceFrame(src, idx, nT, nE));

end

function frame = iReadImageSourceFrame(src, linearIdx, nT, nE)
%IREADIMAGESOURCEFRAME Read one frame from a source object.

if nE > 1
    [tIdx, eIdx] = ind2sub([nT, nE], linearIdx);
    frame = src.getFrame(tIdx, eIdx);
else
    frame = src.getFrame(linearIdx);
end

frame = iConvertPageToSingleGray(frame);

end

function referenceImage = iAverageIndexedFrames(frameIdx, readFrameFcn)
%IAVERAGEINDEXEDFRAMES Numerically stable running average of frames/pages.

meanImage = [];
nValid = 0;

for iFrame = 1:numel(frameIdx)
    thisImage = readFrameFcn(frameIdx(iFrame));
    thisImage = iConvertPageToSingleGray(thisImage);

    if isempty(meanImage)
        meanImage = zeros(size(thisImage), 'single');
    elseif ~isequal(size(thisImage), size(meanImage))
        error('readImageReferenceSource:FrameSizeMismatch', ...
            'All source frames/pages must have the same spatial size.');
    end

    nValid = nValid + 1;
    meanImage = meanImage + (thisImage - meanImage) ./ nValid;
end

referenceImage = meanImage;

end

% =========================================================================
% File/page conversion helpers
% =========================================================================

function img = iReadTiffPage(tiffFile, pageIdx, info)
%IREADTIFFPAGE Read one TIFF page.

try
    img = imread(tiffFile, pageIdx, 'Info', info);
catch
    img = imread(tiffFile, pageIdx);
end

end

function img = iConvertPageToSingleGray(img)
%ICONVERTPAGETOSINGLEGRAY Convert source page/frame to 2D single grayscale.

if isempty(img)
    error('readImageReferenceSource:EmptyFrame', ...
        'Encountered an empty source frame/page.');
end

img = squeeze(img);

if issparse(img)
    img = full(img);
end

img = single(img);

if ndims(img) == 2
    return
end

if ndims(img) == 3 && size(img, 3) == 3
    img = 0.2989 .* img(:, :, 1) + ...
          0.5870 .* img(:, :, 2) + ...
          0.1140 .* img(:, :, 3);
    return
end

if ndims(img) == 3 && size(img, 3) == 1
    img = img(:, :, 1);
    return
end

error('readImageReferenceSource:UnsupportedFrameShape', ...
    'Source frames/pages must be 2D grayscale or RGB images.');

end

function img = iValidateReferenceImage(img)
%IVALIDATEREFERENCEIMAGE Validate and convert reference image to single.

if issparse(img)
    img = full(img);
end

if ~(isnumeric(img) || islogical(img)) || ndims(img) ~= 2 || isempty(img)
    error('readImageReferenceSource:InvalidReferenceImage', ...
        'Reference image must be a non-empty 2D numeric or logical array.');
end

img = single(img);
img(~isfinite(img)) = 0;

end

% =========================================================================
% Source-size and selection helpers
% =========================================================================

function sourceSize = iGetImageSourceSize(src)
%IGETIMAGESOURCESIZE Return [Y X T ...] from a source object.

if ismethod(src, 'getDataSize')
    sourceSize = src.getDataSize();
elseif ismethod(src, 'getSize')
    sourceSize = src.getSize();
elseif isprop(src, 'DataSize')
    sourceSize = src.DataSize;
elseif isprop(src, 'Size')
    sourceSize = src.Size;
elseif isprop(src, 'dataSize')
    sourceSize = src.dataSize;
else
    error('readImageReferenceSource:UnknownSourceSize', ...
        ['Could not determine image source size. Expected source method ' ...
         'getDataSize/getSize or property DataSize/Size.']);
end

sourceSize = double(sourceSize(:).');

if numel(sourceSize) < 2 || any(~isfinite(sourceSize)) || any(sourceSize < 1)
    error('readImageReferenceSource:InvalidSourceSize', ...
        'Image source returned an invalid size vector.');
end

sourceSize = round(sourceSize);

end

function selectedEntry = iSelectUMTEntry(entrySummary)
%ISELECTUMTENTRY Ask user to select one UMT entry.

entryNames = string(entrySummary.Name(:));

displayList = entryNames;

if ismember('Size', entrySummary.Properties.VariableNames)
    displayList = displayList + " | " + string(entrySummary.Size(:));
elseif ismember('DimNames', entrySummary.Properties.VariableNames)
    displayList = displayList + " | " + string(entrySummary.DimNames(:));
end

[idx, ok] = listdlg( ...
    'PromptString', 'Select the UMT entry to use as reference source:', ...
    'SelectionMode', 'single', ...
    'ListString', cellstr(displayList), ...
    'ListSize', [520 260], ...
    'Name', 'Select UMT entry');

if ok ~= 1 || isempty(idx)
    error('readImageReferenceSource:SelectionCancelled', ...
        'UMT entry selection cancelled.');
end

selectedEntry = char(entryNames(idx));

end

function selectedVariable = iSelectMatVariable(candidateVars)
%ISELECTMATVARIABLE Ask user to select one 2D numeric MAT-file variable.

displayList = strings(numel(candidateVars), 1);

for iVar = 1:numel(candidateVars)
    displayList(iVar) = sprintf('%s | %s | %s', ...
        candidateVars(iVar).name, ...
        candidateVars(iVar).class, ...
        mat2str(candidateVars(iVar).size));
end

[idx, ok] = listdlg( ...
    'PromptString', 'Select the 2D numeric variable to use as reference image:', ...
    'SelectionMode', 'single', ...
    'ListString', cellstr(displayList), ...
    'ListSize', [520 260], ...
    'Name', 'Select MAT variable');

if ok ~= 1 || isempty(idx)
    error('readImageReferenceSource:SelectionCancelled', ...
        'MAT variable selection cancelled.');
end

selectedVariable = candidateVars(idx).name;

end

function tf = iIs2DNumericMatVariable(v)
%IIS2DNUMERICMATVARIABLE True for non-scalar 2D numeric arrays.

numericClasses = { ...
    'double', ...
    'single', ...
    'uint8', ...
    'uint16', ...
    'uint32', ...
    'uint64', ...
    'int8', ...
    'int16', ...
    'int32', ...
    'int64'};

tf = ismember(v.class, numericClasses) && ...
    numel(v.size) == 2 && ...
    all(v.size > 1);

end

function sourceInfo = iBaseSourceInfo(sourceFile, sourceType, imageSizeYX)
%IBASESOURCEINFO Create common sourceInfo fields.

[sourceFolder, sourceName, sourceExt] = fileparts(sourceFile);

sourceInfo = struct();
sourceInfo.sourceFile = char(sourceFile);
sourceInfo.sourceFolder = sourceFolder;
sourceInfo.sourceName = sourceName;
sourceInfo.sourceExt = sourceExt;
sourceInfo.sourceType = char(sourceType);
sourceInfo.sourceFrame = [];
sourceInfo.nSourceFrames = [];
sourceInfo.usedFrameIndex = [];
sourceInfo.selectedEntry = '';
sourceInfo.selectedVariable = '';
sourceInfo.imageSizeYX = imageSizeYX;
sourceInfo.processingNotes = '';

end

function sourceFrame = iScalarSourceFrame(frameIdx)
%ISCALARSOURCEFRAME Return sourceFrame only when exactly one frame was used.

if numel(frameIdx) == 1
    sourceFrame = frameIdx;
else
    sourceFrame = [];
end

end
