function outData = getDataFromROI(data, SaveFolder, varargin)
%GETDATAFROMROI Extract ROI-organized data from image-backed inputs.
%
%   outData = getDataFromROI(data, SaveFolder)
%   outData = getDataFromROI(data, SaveFolder, 'ROImasks_filename', fileName, ...)
%
%   This function extracts data from ROI masks stored in an ROI file and
%   returns the result as a UMT structure of kind "roi".
%
%   Supported inputs (any data that has Y and X axes):
%       1) .dat filename with axes Y-X, Y-X-T, Y-X-E, Y-X-F, or Y-X-T-E
%       2) UMT struct of kind "image" whose entries use those same layouts
%       3) Filename to a .umt file containing one such image UMT struct
%   Numeric arrays are not supported. Event-split .dat files carry no event
%   labels: their eventInfo comes from SaveFolder's events.mat
%   (resolveDatEventMapping: one slice per instance or per condition;
%   otherwise one condition, with a warning).
%
%   Inputs:
%       data       - Image-backed input in one of the supported forms above.
%       SaveFolder - Folder used for relative file resolution.
%
%   Name-Value parameters:
%       ROImasks_filename - UMIT .roi file name or full path. A bare
%                           filename is resolved inside SaveFolder.
%                           Default: 'myROI.roi'
%       SpatialAggFcn     - Spatial aggregation across ROI pixels.
%                           Supported:
%                               'none'
%                               'mean'
%                               'max'
%                               'min'
%                               'median'
%                               'mode'
%                               'sum'
%                               'std'
%                           Default: 'mean'
%       FrameRateHz       - Frame rate of DATA (Hz), recorded in the output
%                           entry meta. PipelineManager injects it from the
%                           data; a .dat header provides it otherwise.
%
%   Output:
%       outData    - UMT structure of kind "roi".
%
%   Output dimension conventions:
%       - If SpatialAggFcn ~= 'none':
%             ROI x ...
%       - If SpatialAggFcn == 'none':
%             ROI x Pixel x ...
%
%       where "..." preserves the non-spatial dimensions of the input (T, E,
%       and/or F). 'none' is not available for inputs with an F axis, since
%       the UMT schema has no {ROI,Pixel,F} layout. Labels of the preserved
%       axes of a UMT input (for example F) are carried to the output.
%
%   Notes:
%       - ROI files are read through loadROIFile(...), which migrates and
%         validates the current UMIT .roi schema. Pre-.roi ROI files are
%         not supported.
%       - Top-level eventInfo and per-entry meta present on a UMT input are
%         carried through to the roi UMT output unchanged. Nothing is
%         invented: continuous inputs produce no eventInfo, and an
%         event-split input without eventInfo is still accepted.
%       - A .dat file is streamed: only the columns that contain ROI pixels
%         are read, in blocks of consecutive frames, and every block is
%         aggregated as soon as it is read (the aggregation runs across the
%         pixels of each frame, so it is exact for every SpatialAggFcn). The
%         recording is never loaded whole; memory is one block plus
%         nROI x frames. 'none' keeps the ROI pixels of every frame.
%       - With SpatialAggFcn='none', the Pixel dimension is sized to the
%         largest ROI's pixel count; every smaller ROI is NaN-padded to
%         that length.

default_Output = 'ROI_data.umt';
validAgg = {'none','mean','max','min','median','mode','sum','std'};
if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));


addParameter(p, 'ROImasks_filename', 'myROI.roi', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'SpatialAggFcn', 'mean', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), validAgg));
addParameter(p, 'FrameRateHz', []);

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
roiFile = char(string(p.Results.ROImasks_filename));
spatialAggFcn = lower(char(string(p.Results.SpatialAggFcn)));

if ~isfolder(SaveFolder)
    error('Umitoolbox:getDataFromROI:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Locate and load the ROI file
% -------------------------------------------------------------------------
roiSet = iLoadROISet(roiFile, SaveFolder, 'getDataFromROI');

% -------------------------------------------------------------------------
% Resolve the input and extract the ROI-organized entries
% -------------------------------------------------------------------------
[entries, srcEventInfo, srcLabels] = iResolveInput( ...
    data, SaveFolder, p.Results.FrameRateHz, roiSet, spatialAggFcn);

outData = struct();
for iEntry = 1:numel(entries)
    if iEntry == 1
        outData = genUMTStruct( ...
            entries(iEntry).roiValue, ...
            'kind', 'roi', ...
            'entryName', entries(iEntry).name, ...
            'dimNames', entries(iEntry).roiDims, ...
            'labels', iBuildROILabels(roiSet, entries(iEntry).roiValue, ...
                entries(iEntry).roiDims, spatialAggFcn, srcLabels), ...
            'meta', entries(iEntry).meta);
    else
        outData = genUMTStruct( ...
            outData, ...
            'value', entries(iEntry).roiValue, ...
            'entryName', entries(iEntry).name, ...
            'dimNames', entries(iEntry).roiDims, ...
            'meta', entries(iEntry).meta);
    end
end

% Carry the source event metadata through unchanged: a UMT input's own
% eventInfo, or for an event-split .dat the mapping onto events.mat.
% Without an E dimension the schema forbids eventInfo, and an event-split
% UMT that carried none must not have one invented for it.
if any(arrayfun(@(e) any(strcmp(e.roiDims, 'E')), entries)) && ...
        isstruct(srcEventInfo) && ~isempty(fieldnames(srcEventInfo))
    outData = appendUMTEventInfo(outData, ...
        'eventInfo', srcEventInfo, ...
        'overwrite', true);
end

validateUMTStruct(outData, 'requireEventInfo', false);

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Extract ROI-organized data from image-backed inputs and return a roi UMT.');

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Image-backed input with Y and X axes. Accepted forms: .dat ' ...
             'filename, image UMT struct, or .umt file containing one image ' ...
             'UMT struct.'], ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'file');

        info = PipelineManager.addInput( ...
            info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder used for relative path resolution.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'ROImasks_filename', ...
            'parameter', ...
            'UMIT .roi file name or full path, resolved inside SaveFolder.', ...
            'kind', 'parameter', ...
            'default', 'myROI.roi', ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'SpatialAggFcn', ...
            'parameter', ...
            'Spatial aggregation across ROI pixels.', ...
            'kind', 'parameter', ...
            'default', 'mean', ...
            'allowed', validAgg, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'ROI-organized UMT output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function layouts = iSupportedLayouts()
%ISUPPORTEDLAYOUTS Image layouts that have Y and X and a ROI-schema output.
layouts = {{'Y','X'}, {'Y','X','T'}, {'Y','X','E'}, {'Y','X','F'}, {'Y','X','T','E'}};
end

function [entries, eventInfo, srcLabels] = iResolveInput(data, SaveFolder, explicitRate, roiSet, spatialAggFcn)
%IRESOLVEINPUT Resolve supported input forms to ROI-organized entries.
%
% entries(k) has name, roiValue, roiDims, and meta (the source entry's meta
% struct; for a .dat file it holds the data's own frame rate, meta.FrameRateHz,
% when known: the explicit or injected FrameRateHz, else the .dat header;
% AcqInfos.mat is not used). eventInfo carries the source's shared top-level
% event metadata, or an empty struct when there is none. srcLabels holds the
% labels of a UMT input.

eventInfo = struct();
srcLabels = struct();

if isnumeric(data) || islogical(data)
    error('Umitoolbox:getDataFromROI:UnsupportedInputType', ...
        ['Numeric arrays are not supported. Input "data" must be a .dat ' ...
         'filename, an image UMT struct, or a .umt filename.']);
end

if ischar(data) || (isstring(data) && isscalar(data))
    dataFile = char(string(data));

    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('Umitoolbox:getDataFromROI:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);

    switch ext
        case '.dat'
            [entries, dims, sizes] = iEntryFromDat(dataFile, roiSet, spatialAggFcn, explicitRate);
            if any(strcmp(dims, 'E'))
                eventInfo = iEventInfoFromFolder(dims, sizes, SaveFolder);
            end
            return

        case '.umt'
            data = loadData(dataFile);

        otherwise
            error('Umitoolbox:getDataFromROI:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

assert(isstruct(data) && isscalar(data), ...
    'Umitoolbox:getDataFromROI:UnsupportedInputType', ...
    ['Input "data" must be a .dat filename, an image UMT struct, ' ...
     'or a .umt file containing one image UMT struct.']);

validateUMTStruct(data, 'requireEventInfo', false);

assert(strcmpi(char(string(data.kind)), 'image'), ...
    'Umitoolbox:getDataFromROI:InvalidUMTKind', ...
    'Input UMT must have kind = "image".');

entryNames = fieldnames(data.data);
assert(~isempty(entryNames), ...
    'Umitoolbox:getDataFromROI:EmptyUMTData', ...
    'Input UMT data is empty.');

layouts = iSupportedLayouts();
entries = repmat(struct('name', '', 'roiValue', [], 'roiDims', {{}}, 'meta', struct()), ...
    numel(entryNames), 1);

for iEntry = 1:numel(entryNames)
    thisEntry = data.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    isAllowed = any(cellfun(@(x) isequal(thisDims, x), layouts));
    assert(isAllowed, ...
        'Umitoolbox:getDataFromROI:InvalidUMTEntryDims', ...
        ['Entry "%s" has unsupported dimNames. Supported image layouts are ' ...
         '{Y,X}, {Y,X,T}, {Y,X,E}, {Y,X,F}, and {Y,X,T,E}.'], ...
        entryNames{iEntry});

    [roiPix, trailingDims, trailingSz] = iGatherROIPixels(single(thisEntry.value), ...
        thisDims, roiSet);
    [roiValue, roiDims] = iAssembleOutput(iAggregateAll(roiPix, spatialAggFcn), roiPix, ...
        trailingDims, trailingSz, spatialAggFcn);

    entries(iEntry).name = entryNames{iEntry};
    entries(iEntry).roiValue = roiValue;
    entries(iEntry).roiDims = roiDims;
    if isfield(thisEntry, 'meta') && isstruct(thisEntry.meta) && isscalar(thisEntry.meta)
        entries(iEntry).meta = thisEntry.meta;
    end
end

if isfield(data, 'eventInfo')
    eventInfo = data.eventInfo;
end
if isfield(data, 'labels') && isstruct(data.labels)
    srcLabels = data.labels;
end

end

function [entries, dims, sizes] = iEntryFromDat(dataFile, roiSet, spatialAggFcn, explicitRate)
%IENTRYFROMDAT ROI-organized entry of a .dat file, streamed in frame blocks.

info = loadMetaData(dataFile);
assertDatLayout(info, iSupportedLayouts(), 'getDataFromROI');

dims = cellstr(string(info.dimNames(:).'));
sizes = double(info.dimSizes(:).');
Ny = datAxisSize(info, 'Y');
Nx = datAxisSize(info, 'X');
assert(isequal([Ny, Nx], roiSet.imageSizeYX), ...
    'Umitoolbox:getDataFromROI:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

trailingDims = dims(3:end);
trailingSz = sizes(3:end);
nFrames = max(1, prod(trailingSz));

roiNames = roiSet.names;
nROI = numel(roiNames);
masks = roiSet.masks;
bNone = strcmp(spatialAggFcn, 'none');

% Only the columns that any ROI touches are read. The ROI masks are cut to
% those columns once, so a pixel selection indexes a block directly.
xNeeded = find(any(any(cat(3, masks{:}), 3), 1));
masksCols = cellfun(@(m) m(:, xNeeded), masks, 'UniformOutput', false);

roiAgg = nan(nROI, nFrames, 'single');
roiPix = cell(nROI, 1);
if bNone
    for iROI = 1:nROI
        roiPix{iROI} = nan(nnz(masks{iROI}), nFrames, 'single');
    end
end

slabIn = spatialSlabIO('open', dataFile, 'Info', info);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

% Consecutive frames per block: the block as read and its single-precision
% copy are alive together. Blocks of frames (not of columns) keep the file
% read in one pass, one contiguous read per frame, whatever the memory budget.
blockBytes = double(Ny) * numel(xNeeded) * nFrames * 4;
nBlocks = calculateMaxChunkSize(blockBytes, 2, 0.1);
framesPerBlock = max(1, ceil(nFrames / nBlocks));

for t0 = 1:framesPerBlock:nFrames
    tIdx = t0:min(t0 + framesPerBlock - 1, nFrames);
    block = reshape(single(spatialSlabIO('read', slabIn, xNeeded, tIdx)), ...
        Ny * numel(xNeeded), numel(tIdx));

    for iROI = 1:nROI
        pixVals = block(masksCols{iROI}(:), :);
        if bNone
            roiPix{iROI}(:, tIdx) = pixVals;
        else
            roiAgg(iROI, tIdx) = iApplyAggFcn(pixVals, spatialAggFcn);
        end
    end
end

[roiValue, roiDims] = iAssembleOutput(roiAgg, roiPix, trailingDims, trailingSz, spatialAggFcn);

entries = struct('name', 'main', 'roiValue', roiValue, 'roiDims', {roiDims}, ...
    'meta', iRateMeta(explicitRate, dataFile));
end

% =========================================================================
% Local helper: load and normalize a UMIT .roi file
% =========================================================================
function roiSet = iLoadROISet(roiFile, SaveFolder, callerName)
%ILOADROISET Resolve and load a UMIT .roi file into names, masks, and size.

roiFile = char(string(roiFile));
if ~isfile(roiFile)
    roiFile = fullfile(SaveFolder, roiFile);
end

if ~isfile(roiFile)
    error(sprintf('Umitoolbox:%s:MissingROIFile', callerName), ...
        'ROI file was not found: "%s".', roiFile);
end

[~, ~, ext] = fileparts(roiFile);
if ~strcmpi(ext, '.roi')
    error(sprintf('Umitoolbox:%s:UnsupportedROIFile', callerName), ...
        ['ROI files must use the UMIT ".roi" format. Pre-.roi ROI files ' ...
         'are not supported. Received: "%s".'], roiFile);
end

% loadROIFile migrates and validates the schema, so masks are guaranteed to
% be 2-D and to match imageInfo.imageSizeYX, and ROI names are unique.
ROIFile = loadROIFile(roiFile);

if isempty(ROIFile.ROIs)
    error(sprintf('Umitoolbox:%s:EmptyROIFile', callerName), ...
        'ROI file "%s" does not contain any ROI.', roiFile);
end

roiSet = struct();
roiSet.filePath = roiFile;
roiSet.imageSizeYX = double(ROIFile.imageInfo.imageSizeYX(:).');
roiSet.names = cellstr(string({ROIFile.ROIs.name}))';
roiSet.masks = arrayfun(@(r) logical(r.mask), ROIFile.ROIs, ...
    'UniformOutput', false);
roiSet.masks = roiSet.masks(:);

end

function [roiPix, trailingDims, trailingSz] = iGatherROIPixels(value, dimNames, roiSet)
%IGATHERROIPIXELS ROI pixel values of an in-RAM image entry.
%
% roiPix{k} is nPixels-by-nFrames, with the frames of the non-spatial axes
% flattened in their stored order. trailingDims/trailingSz describe those axes.

[~,yxLoc] = ismember({'Y','X'}, dimNames);

dataSz = size(value);
if numel(dataSz) < numel(dimNames)
    dataSz(end+1:numel(dimNames)) = 1;
end

assert(isequal(dataSz(yxLoc), roiSet.imageSizeYX), ...
    'Umitoolbox:getDataFromROI:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

origDim = 1:numel(dimNames);
newDim = [yxLoc, setdiff(origDim, yxLoc, 'stable')];

value = permute(value, newDim);
permDims = dimNames(newDim);
permSz = dataSz(newDim);
trailingDims = permDims(3:end);
trailingSz = permSz(3:end);

value2D = reshape(value, prod(permSz(1:2)), []);
roiPix = cell(numel(roiSet.masks), 1);
for iROI = 1:numel(roiPix)
    roiPix{iROI} = value2D(roiSet.masks{iROI}(:), :);
end
end

function roiAgg = iAggregateAll(roiPix, spatialAggFcn)
%IAGGREGATEALL Aggregate every ROI across its pixels (nROI-by-nFrames).
if strcmp(spatialAggFcn, 'none')
    roiAgg = [];
    return
end
nFrames = size(roiPix{1}, 2);
roiAgg = nan(numel(roiPix), nFrames, 'single');
for iROI = 1:numel(roiPix)
    roiAgg(iROI, :) = iApplyAggFcn(roiPix{iROI}, spatialAggFcn);
end
end

function [roiValue, roiDims] = iAssembleOutput(roiAgg, roiPix, trailingDims, trailingSz, spatialAggFcn)
%IASSEMBLEOUTPUT Shape ROI values as ROI x [Pixel] x trailing axes.
%
% roiAgg is nROI-by-nFrames (any aggregation); roiPix holds the per-ROI pixel
% values (nPixels-by-nFrames) and is only used by SpatialAggFcn 'none'.

nROI = numel(roiPix);

if strcmp(spatialAggFcn, 'none')
    assert(~any(strcmp(trailingDims, 'F')), ...
        'Umitoolbox:getDataFromROI:NoneNotSupportedForF', ...
        ['SpatialAggFcn "none" is not available for inputs with an F axis: ' ...
         'the UMT schema has no {ROI,Pixel,F} layout. Use an aggregation function.']);

    pixelCounts = cellfun(@(c) size(c, 1), roiPix);
    maxPixel = max(pixelCounts);

    if isempty(trailingSz)
        roiValue = nan(nROI, maxPixel, 'single');
        for iROI = 1:nROI
            roiValue(iROI, 1:pixelCounts(iROI)) = single(roiPix{iROI});
        end
    else
        roiValue = nan([nROI, maxPixel, trailingSz], 'single');
        for iROI = 1:nROI
            % The trailing-dimension assignment indexes each declared
            % trailing dimension with its own ':' range, which needs a
            % shape-conforming right-hand side.
            thisVal = reshape(single(roiPix{iROI}), [pixelCounts(iROI), trailingSz]);
            idx = repmat({':'}, 1, ndims(roiValue));
            idx{1} = iROI;
            idx{2} = 1:pixelCounts(iROI);
            roiValue(idx{:}) = thisVal;
        end
    end

    roiDims = [{'ROI','Pixel'}, trailingDims];
else
    if isempty(trailingSz)
        roiValue = reshape(single(roiAgg), [nROI, 1]);
    else
        roiValue = reshape(single(roiAgg), [nROI, trailingSz]);
    end

    roiDims = [{'ROI'}, trailingDims];
end

end

function labels = iBuildROILabels(roiSet, roiValue, roiDims, spatialAggFcn, srcLabels)
%IBUILDROILABELS Build shared display/reference labels for roi output.

roiNames = roiSet.names;

labels = struct();
labels.ROI = roiNames(:).';

if strcmp(spatialAggFcn, 'none')
    pixelLen = size(roiValue, 2);
    labels.Pixel = arrayfun(@num2str, 1:pixelLen, 'UniformOutput', false);
end

% Labels of the preserved axes (for example F = Amplitude/Phase) stay valid.
for fieldName = {'T', 'E', 'F'}
    if isfield(srcLabels, fieldName{1}) && ismember(fieldName{1}, roiDims)
        labels.(fieldName{1}) = srcLabels.(fieldName{1});
    end
end

labelFields = fieldnames(labels);
usedDims = roiDims;

for iField = numel(labelFields):-1:1
    if ~ismember(labelFields{iField}, usedDims)
        labels = rmfield(labels, labelFields{iField});
    end
end

end

function out = iApplyAggFcn(vals, fcnName)
%IAPPLYAGGFCN Apply one supported aggregation across the first dimension.

switch fcnName
    case 'mean'
        out = mean(vals, 1, 'omitnan');
    case 'median'
        out = median(vals, 1, 'omitnan');
    case 'mode'
        out = mode(vals, 1);
    case 'std'
        out = std(vals, 0, 1, 'omitnan');
    case 'max'
        out = max(vals, [], 1, 'omitnan');
    case 'min'
        out = min(vals, [], 1, 'omitnan');
    case 'sum'
        out = sum(vals, 1, 'omitnan');
    otherwise
        out = vals;
end

out = single(out);
end

function eventInfo = iEventInfoFromFolder(dims, sz, SaveFolder)
%IEVENTINFOFROMFOLDER eventInfo of an event-split input without labels, from
% the SaveFolder's events.mat (resolveDatEventMapping, .dat header Phase 8c).
sizes = ones(1, numel(dims));
sizes(1:min(numel(sz), numel(dims))) = sz(1:min(numel(sz), numel(dims)));
mapping = resolveDatEventMapping(struct('filePath', 'input data', ...
    'dimNames', {dims}, 'dimSizes', sizes), SaveFolder);
if ~any(strcmpi(mapping.status, {'matched', 'aggregated'}))
    warning('Umitoolbox:getDataFromROI:eventsNotMatched', '%s', mapping.message);
end
eventInfo = mapping.eventInfo;
end

function meta = iRateMeta(explicitRate, dataFile)
%IRATEMETA Entry meta with the data's own frame rate, when it is known.
meta = struct();
if isempty(explicitRate) && isempty(dataFile)
    return
end
try
    meta.FrameRateHz = resolveDataInfoValue('frameRateHz', explicitRate, dataFile, 'getDataFromROI');
catch ME
    if ~endsWith(ME.identifier, ':missingFrameRateHz')
        rethrow(ME)
    end
end
end
