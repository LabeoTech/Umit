function outData = getDataFromROI(data, SaveFolder, varargin)
%GETDATAFROMROI Extract ROI-organized data from image-backed inputs.
%
%   outData = getDataFromROI(data, SaveFolder)
%   outData = getDataFromROI(data, SaveFolder, 'ROImasks_filename', fileName, ...)
%
%   This function extracts data from ROI masks stored in an ROI file and
%   returns the result as a UMT structure of kind "roi".
%
%   Supported inputs:
%       1) Numeric image data with dimensions Y x X x T
%       2) Raw .dat filename storing continuous Y x X x T data
%       3) UMT struct of kind "image"
%       4) Filename to a .umt file containing one image UMT struct
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
%       where "..." preserves the non-spatial dimensions of the input.
%
%   Notes:
%       - ROI files are read through loadROIFile(...), which migrates and
%         validates the current UMIT .roi schema. Pre-.roi ROI files are
%         not supported.
%       - Top-level eventInfo and per-entry meta present on a UMT input are
%         carried through to the roi UMT output unchanged. Nothing is
%         invented: continuous inputs produce no eventInfo, and an
%         event-split input without eventInfo is still accepted.
%       - Raw .dat input is treated as continuous YXT data.
%       - Event-aware processing for raw .dat should be handled upstream by
%         converting data to a UMT image structure with E as the last
%         dimension when needed.


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
% Resolve input to one or more image entries
% -------------------------------------------------------------------------
[entryNames, entryValues, entryDims, entryMetas, srcEventInfo] = ...
    iResolveImageInput(data, SaveFolder);

outData = struct();
roiEntryNames = entryNames;
roiEntryValues = cell(size(entryValues));
roiEntryDims = cell(size(entryDims));

for iEntry = 1:numel(entryNames)
    thisValue = entryValues{iEntry};
    thisDims = entryDims{iEntry};

    [roiValue, roiDims] = iExtractROIsFromEntry(thisValue, thisDims, roiSet, spatialAggFcn);

    roiEntryValues{iEntry} = roiValue;
    roiEntryDims{iEntry} = roiDims;
end

for iEntry = 1:numel(roiEntryNames)
    if iEntry == 1
        outData = genUMTStruct( ...
            roiEntryValues{iEntry}, ...
            'kind', 'roi', ...
            'entryName', roiEntryNames{iEntry}, ...
            'dimNames', roiEntryDims{iEntry}, ...
            'labels', iBuildROILabels(roiSet, roiEntryValues{iEntry}, roiEntryDims{iEntry}, spatialAggFcn), ...
            'meta', entryMetas{iEntry}, ...
            'SaveFolder', SaveFolder);
    else
        outData = genUMTStruct( ...
            outData, ...
            'value', roiEntryValues{iEntry}, ...
            'entryName', roiEntryNames{iEntry}, ...
            'dimNames', roiEntryDims{iEntry}, ...
            'meta', entryMetas{iEntry}, ...
            'SaveFolder', SaveFolder);
    end
end

% Carry the source event metadata through unchanged. Both conditions matter:
% without an E dimension the schema forbids eventInfo, and an event-split
% input that carried none must not have one invented for it.
if any(cellfun(@(d) any(strcmp(d, 'E')), roiEntryDims)) && ...
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
            ['Image-backed input. Accepted forms: YXT array, .dat filename, ' ...
             'image UMT struct, or .umt file containing one image UMT struct.'], ...
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

function [entryNames, entryValues, entryDims, entryMetas, eventInfo] = ...
    iResolveImageInput(data, SaveFolder)
%IRESOLVEIMAGEINPUT Resolve supported input forms to image entries.
%
% entryMetas carries each source entry's meta struct (empty when the input
% is a raw array or .dat). eventInfo carries the source UMT's shared
% top-level event metadata, or an empty struct when there is none.

% entryMetas is assigned on every return path below. eventInfo needs a
% default because only the UMT path can supply one.
eventInfo = struct();

if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, ...
        mfilename, 'data');

    entryNames = {'main'};
    entryValues = {single(data)};
    entryDims = {{'Y','X','T'}};
    entryMetas = {struct()};
    return
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
            rawData = loadData(dataFile);
            entryNames = {'main'};
            entryValues = {single(rawData)};
            entryDims = {{'Y','X','T'}};
            entryMetas = {struct()};
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
    ['Input "data" must be a YXT array, a .dat filename, an image UMT struct, ' ...
     'or a .umt file containing one image UMT struct.']);

validateUMTStruct(data, 'requireEventInfo', false);

assert(strcmpi(char(string(data.kind)), 'image'), ...
    'Umitoolbox:getDataFromROI:InvalidUMTKind', ...
    'Input UMT must have kind = "image".');

entryNames = fieldnames(data.data);
assert(~isempty(entryNames), ...
    'Umitoolbox:getDataFromROI:EmptyUMTData', ...
    'Input UMT data is empty.');

entryValues = cell(size(entryNames));
entryDims = cell(size(entryNames));
entryMetas = cell(size(entryNames));

allowedDims = { ...
    {'Y','X'}, ...
    {'Y','X','T'}, ...
    {'Y','X','E'}, ...
    {'Y','X','T','E'}};

for iEntry = 1:numel(entryNames)
    thisEntry = data.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    isAllowed = any(cellfun(@(x) isequal(thisDims, x), allowedDims));
    assert(isAllowed, ...
        'Umitoolbox:getDataFromROI:InvalidUMTEntryDims', ...
        ['Entry "%s" has unsupported dimNames. Supported image layouts are ' ...
         '{Y,X}, {Y,X,T}, {Y,X,E}, and {Y,X,T,E}.'], ...
        entryNames{iEntry});

    entryValues{iEntry} = single(thisEntry.value);
    entryDims{iEntry} = thisDims;

    if isfield(thisEntry, 'meta') && isstruct(thisEntry.meta) && isscalar(thisEntry.meta)
        entryMetas{iEntry} = thisEntry.meta;
    else
        entryMetas{iEntry} = struct();
    end
end

if isfield(data, 'eventInfo')
    eventInfo = data.eventInfo;
end

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

function [roiValue, roiDims] = iExtractROIsFromEntry(value, dimNames, roiSet, spatialAggFcn)
%IEXTRACTROISFROMENTRY Extract ROI-organized values from one image entry.

[~,yxLoc] = ismember({'Y','X'}, dimNames);

dataSz = size(value);
if numel(dataSz) < numel(dimNames)
    dataSz(end+1:numel(dimNames)) = 1;
end

refFrameSz = roiSet.imageSizeYX;
assert(isequal(dataSz(yxLoc), refFrameSz), ...
    'Umitoolbox:getDataFromROI:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

origDim = 1:numel(dimNames);
newDim = [yxLoc, setdiff(origDim, yxLoc, 'stable')];

value = permute(value, newDim);
permDims = dimNames(newDim);
permSz = dataSz(newDim);

value2D = reshape(value, prod(permSz(1:2)), []);
roiNames = roiSet.names;

if strcmp(spatialAggFcn, 'none')
    pixelCounts = zeros(numel(roiNames), 1);
    roiPixVals = cell(size(roiNames));

    for iROI = 1:numel(roiPixVals)
        roiMask = roiSet.masks{iROI};
        pixVals = value2D(roiMask(:), :);
        pixelCounts(iROI) = size(pixVals, 1);

        if numel(permSz) > 2
            % The trailing-dimension assignment below indexes each declared
            % trailing dimension with its own ':' range, which only accepts
            % a flattened RHS when there is a single trailing dimension.
            pixVals = reshape(pixVals, [size(pixVals,1), permSz(3:end)]);
        end

        roiPixVals{iROI} = single(pixVals);
    end

    maxPixel = max(pixelCounts);
    trailingSz = permSz(3:end);

    if isempty(trailingSz)
        roiValue = nan(numel(roiNames), maxPixel, 'single');
        for iROI = 1:numel(roiNames)
            roiValue(iROI, 1:pixelCounts(iROI)) = roiPixVals{iROI};
        end
    else
        roiValue = nan([numel(roiNames), maxPixel, trailingSz], 'single');
        for iROI = 1:numel(roiNames)
            thisVal = roiPixVals{iROI};
            idx = repmat({':'}, 1, ndims(roiValue));
            idx{1} = iROI;
            idx{2} = 1:pixelCounts(iROI);
            roiValue(idx{:}) = thisVal;
        end
    end

    roiDims = [{'ROI','Pixel'}, permDims(3:end)];

else
    roiPixVals = cell(size(roiNames));

    for iROI = 1:numel(roiPixVals)
        roiMask = roiSet.masks{iROI};
        pixVals = value2D(roiMask(:), :);
        pixVals = iApplyAggFcn(pixVals, spatialAggFcn);

        if numel(permSz) > 2
            % See the matching comment in the 'none' branch above: the
            % trailing-dimension assignment below needs a shape-conforming
            % RHS whenever there is more than one trailing dimension.
            pixVals = reshape(pixVals, [size(pixVals,1), permSz(3:end)]);
        end

        roiPixVals{iROI} = single(pixVals);
    end

    trailingSz = permSz(3:end);

    if isempty(trailingSz)
        roiValue = nan(numel(roiNames), 1, 'single');
        for iROI = 1:numel(roiNames)
            roiValue(iROI, 1) = roiPixVals{iROI};
        end
        roiValue = reshape(roiValue, [numel(roiNames), 1]);
    else
        roiValue = nan([numel(roiNames), trailingSz], 'single');
        for iROI = 1:numel(roiNames)
            idx = repmat({':'}, 1, ndims(roiValue));
            idx{1} = iROI;
            roiValue(idx{:}) = roiPixVals{iROI};
        end
    end

    roiDims = [{'ROI'}, permDims(3:end)];
end

end

function labels = iBuildROILabels(roiSet, roiValue, roiDims, spatialAggFcn)
%IBUILDROILABELS Build shared display/reference labels for roi output.

roiNames = roiSet.names;

labels = struct();
labels.ROI = roiNames(:).';

if strcmp(spatialAggFcn, 'none')
    pixelLen = size(roiValue, 2);
    labels.Pixel = arrayfun(@num2str, 1:pixelLen, 'UniformOutput', false);
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
