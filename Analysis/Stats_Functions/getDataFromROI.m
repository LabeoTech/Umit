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
%       ROImasks_filename - ROI file name or full path.
%                           Default: 'ROImasks_data.mat'
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
%       object            - Optional legacy Acquisition / Modality object.
%                           This is kept only for compatibility with the
%                           current ROI-file lookup logic.
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
%       - The ROI file management path still relies on findMyROIfile(...).
%         This lookup logic should be revised once ROI file creation and
%         management are fully modernized.
%       - Raw .dat input is treated as continuous YXT data.
%       - Event-aware processing for raw .dat should be handled upstream by
%         converting data to a UMT image structure with E as the last
%         dimension when needed.


default_Output = 'ROI_data.umt'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));

validAgg = {'none','mean','max','min','median','mode','sum','std'};
addParameter(p, 'ROImasks_filename', 'ROImasks_data.mat', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'SpatialAggFcn', 'mean', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), validAgg));
addParameter(p, 'object', [], ...
    @(x) isempty(x) || isa(x, 'Acquisition') || isa(x, 'Modality'));

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
roiFile = char(string(p.Results.ROImasks_filename));
spatialAggFcn = lower(char(string(p.Results.SpatialAggFcn)));
object = p.Results.object;

if ~isfolder(SaveFolder)
    error('Umitoolbox:getDataFromROI:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% -------------------------------------------------------------------------
% Locate ROI file
% -------------------------------------------------------------------------
% NOTE:
% This legacy ROI lookup path is intentionally kept for now because ROI
% file creation/management has not yet been fully modernized. It should be
% revised in a later pass.
roiFile = findMyROIfile(roiFile, object);

if ~isfile(roiFile)
    error('Umitoolbox:getDataFromROI:MissingROIFile', ...
        'ROI file was not found: "%s".', roiFile);
end

roiData = load(roiFile);
assert(isfield(roiData, 'ROI_info') && isstruct(roiData.ROI_info), ...
    'Umitoolbox:getDataFromROI:InvalidROIFile', ...
    'ROI file must contain a struct variable named "ROI_info".');
assert(isfield(roiData, 'img_info') && isstruct(roiData.img_info), ...
    'Umitoolbox:getDataFromROI:InvalidROIFile', ...
    'ROI file must contain a struct variable named "img_info".');
assert(isfield(roiData.img_info, 'imageData') && ~isempty(roiData.img_info.imageData), ...
    'Umitoolbox:getDataFromROI:InvalidROIFile', ...
    'ROI file img_info must contain a non-empty field "imageData".');

% -------------------------------------------------------------------------
% Resolve input to one or more image entries
% -------------------------------------------------------------------------
[entryNames, entryValues, entryDims] = iResolveImageInput(data, SaveFolder);

outData = struct();
roiEntryNames = entryNames;
roiEntryValues = cell(size(entryValues));
roiEntryDims = cell(size(entryDims));

for iEntry = 1:numel(entryNames)
    thisValue = entryValues{iEntry};
    thisDims = entryDims{iEntry};

    [roiValue, roiDims] = iExtractROIsFromEntry(thisValue, thisDims, roiData, spatialAggFcn);

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
            'labels', iBuildROILabels(roiData, roiEntryValues{iEntry}, roiEntryDims{iEntry}, spatialAggFcn));
    else
        outData = genUMTStruct( ...
            outData, ...
            'value', roiEntryValues{iEntry}, ...
            'entryName', roiEntryNames{iEntry}, ...
            'dimNames', roiEntryDims{iEntry});
    end
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
            'ROI file name or full path.', ...
            'kind', 'parameter', ...
            'default', 'ROImasks_data.mat', ...
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
            'object', ...
            'parameter', ...
            'Optional legacy Acquisition/Modality object for ROI file lookup.', ...
            'kind', 'parameter', ...
            'default', [], ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ProcessedData'}, ...
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

function [entryNames, entryValues, entryDims] = iResolveImageInput(data, SaveFolder)
%IRESOLVEIMAGEINPUT Resolve supported input forms to image entries.

if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, ...
        mfilename, 'data');

    entryNames = {'main'};
    entryValues = {single(data)};
    entryDims = {{'Y','X','T'}};
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
end

end

function [roiValue, roiDims] = iExtractROIsFromEntry(value, dimNames, roiData, spatialAggFcn)
%IEXTRACTROISFROMENTRY Extract ROI-organized values from one image entry.

[~,yxLoc] = ismember({'Y','X'}, dimNames);

dataSz = size(value);
if numel(dataSz) < numel(dimNames)
    dataSz(end+1:numel(dimNames)) = 1;
end

refFrameSz = size(roiData.img_info.imageData);
assert(isequal(dataSz(yxLoc), refFrameSz), ...
    'Umitoolbox:getDataFromROI:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

origDim = 1:numel(dimNames);
newDim = [yxLoc, setdiff(origDim, yxLoc, 'stable')];

value = permute(value, newDim);
permDims = dimNames(newDim);
permSz = dataSz(newDim);

value2D = reshape(value, prod(permSz(1:2)), []);
roiNames = {roiData.ROI_info.Name}';

if strcmp(spatialAggFcn, 'none')
    pixelCounts = zeros(numel(roiNames), 1);
    roiPixVals = cell(size(roiNames));

    for iROI = 1:numel(roiPixVals)
        roiMask = roiData.ROI_info(iROI).Stats.ROI_binary_mask;
        pixVals = value2D(roiMask(:), :);
        pixelCounts(iROI) = size(pixVals, 1);

        if numel(permSz) > 2
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
        roiMask = roiData.ROI_info(iROI).Stats.ROI_binary_mask;
        pixVals = value2D(roiMask(:), :);
        pixVals = iApplyAggFcn(pixVals, spatialAggFcn);

        if numel(permSz) > 2
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

function labels = iBuildROILabels(roiData, roiValue, roiDims, spatialAggFcn)
%IBUILDROILABELS Build shared display/reference labels for roi output.

roiNames = {roiData.ROI_info.Name}';

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
