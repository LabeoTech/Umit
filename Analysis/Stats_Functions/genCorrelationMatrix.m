function outFile = genCorrelationMatrix(data, SaveFolder, varargin)
%GENCORRELATIONMATRIX Generate ROI correlation matrices from image-backed data.
%
%   outFile = genCorrelationMatrix(data, SaveFolder)
%   outFile = genCorrelationMatrix(data, SaveFolder, 'ROImasks_filename', fileName, ...)
%
%   Supported inputs:
%       1) Numeric Y x X x T array
%       2) Raw .dat filename storing continuous Y x X x T data
%       3) Image UMT struct
%       4) .umt filename containing one image UMT struct
%
%   Name-Value parameters:
%       ROImasks_filename    - ROI file name or full path.
%                              Default: 'myROI.roimsk'
%       CorrAlgorithm        - Correlation algorithm:
%                              'centroid_vs_centroid'
%                              'centroid_vs_agg'
%                              'avg_vs_avg'
%                              Default: 'centroid_vs_centroid'
%       SpatialAggFcn        - Spatial aggregation used by
%                              'centroid_vs_agg':
%                              'mean','max','min','median'
%                              Default: 'mean'
%       b_FisherZ_transform  - Apply truncated Fisher Z transform.
%                              Default: false
%       b_genSPCMaps         - Generate seed-pixel correlation maps and
%                              save them as an image UMT file in SaveFolder.
%                              Default: false
%       object               - Optional legacy Acquisition / Modality object.
%                              Kept only for ROI-file lookup compatibility.
%
%   Output:
%       outFile - File manifest cell array containing the generated UMT
%                 file name(s) saved in SaveFolder.
%
%   Notes:
%       - The ROI file lookup still relies on findMyROIfile(...). That path
%         is intentionally preserved for now and should be revised later.
%       - The main output is a roi UMT correlation matrix file.
%       - SPC maps are saved as a second image UMT file when requested.


default_Output = {'corrMatrix.umt', 'corrMatrix_SPCMaps.umt'};

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outFile = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROImasks_filename', 'myROI.roimsk', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'CorrAlgorithm', 'centroid_vs_centroid', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'centroid_vs_centroid','centroid_vs_agg','avg_vs_avg'}));
addParameter(p, 'SpatialAggFcn', 'mean', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'mean','max','min','median'}));
addParameter(p, 'b_FisherZ_transform', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'b_genSPCMaps', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'object', [], @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality'));

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
roiFile = char(string(p.Results.ROImasks_filename));
corrAlgorithm = lower(char(string(p.Results.CorrAlgorithm)));
spatialAggFcn = lower(char(string(p.Results.SpatialAggFcn)));
bFisherZ = p.Results.b_FisherZ_transform;
bGenSPCMaps = p.Results.b_genSPCMaps;
object = p.Results.object;

assert(isfolder(SaveFolder), 'Umitoolbox:genCorrelationMatrix:InvalidSaveFolder', ...
    'SaveFolder "%s" does not exist.', SaveFolder);

% NOTE:
% This legacy ROI lookup path is intentionally preserved for now. It should
% be revised once ROI file creation/management are fully modernized.
roiFile = findMyROIfile(roiFile, object);
assert(isfile(roiFile), 'Umitoolbox:genCorrelationMatrix:MissingROIFile', ...
    'ROI file was not found: "%s".', roiFile);

roiData = load(roiFile);
assert(isfield(roiData, 'ROI_info') && isstruct(roiData.ROI_info), ...
    'Umitoolbox:genCorrelationMatrix:InvalidROIFile', ...
    'ROI file must contain a struct variable named "ROI_info".');
assert(isfield(roiData, 'img_info') && isstruct(roiData.img_info) && ...
    isfield(roiData.img_info, 'imageData'), ...
    'Umitoolbox:genCorrelationMatrix:InvalidROIFile', ...
    'ROI file must contain img_info.imageData.');

[value, dimNames] = iResolveImageInput(data, SaveFolder);
assert(isequal(dimNames, {'Y','X','T'}), ...
    'Umitoolbox:genCorrelationMatrix:WrongFormat', ...
    'Input data must be an Image time series with dimensions {''Y'',''X'',''T''}.');
assert(isequal(size(value,1), size(roiData.img_info.imageData,1)) && ...
    isequal(size(value,2), size(roiData.img_info.imageData,2)), ...
    'Umitoolbox:genCorrelationMatrix:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

roiNames = {roiData.ROI_info.Name}';
[centroidList, roiMasks] = iExtractROIGeometry(roiData);

corrMatrix = iComputeCorrelationMatrix(value, roiMasks, centroidList, corrAlgorithm, spatialAggFcn);
if bFisherZ
    corrMatrix = iZFisherTruncated(corrMatrix);
end

labels = struct();
labels.ROI = roiNames(:).';
corrUMT = genUMTStruct(corrMatrix, ...
    'kind', 'roi', ...
    'entryName', 'CorrMatrix', ...
    'dimNames', {'ROI','ROI'}, ...
    'labels', labels);

corrFile = 'corrMatrix.umt';
saveData(fullfile(SaveFolder, corrFile), corrUMT);
outFile = {corrFile};

if bGenSPCMaps
    spcMaps = iComputeSPCMaps(value, centroidList);
    if bFisherZ
        for iMap = 1:numel(spcMaps)
            spcMaps{iMap} = iZFisherTruncated(spcMaps{iMap});
        end
    end

    spcUMT = [];
    for iMap = 1:numel(spcMaps)
        entryName = matlab.lang.makeValidName(roiNames{iMap});
        if isempty(entryName)
            entryName = sprintf('ROI_%d', iMap);
        end
        if iMap == 1
            spcUMT = genUMTStruct(spcMaps{iMap}, ...
                'kind', 'image', ...
                'entryName', entryName, ...
                'dimNames', {'Y','X'});
        else
            spcUMT = genUMTStruct(spcUMT, ...
                'value', spcMaps{iMap}, ...
                'entryName', entryName, ...
                'dimNames', {'Y','X'});
        end
    end

    spcFile = 'corrMatrix_SPCMaps.umt';
    saveData(fullfile(SaveFolder, spcFile), spcUMT);
    outFile = [outFile, {spcFile}]; %#ok<AGROW>
end

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Generate ROI correlation matrices from image-backed inputs.');
        info.version = '1.0.0';

        info = PipelineManager.addInput(info, 'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            'Image-backed input.', ...
            'kind', 'input', 'position', 1, 'callType', 'positional', ...
            'isData', true, 'supportsFile', true, 'dataMode', 'either');

        info = PipelineManager.addInput(info, 'SaveFolder', 'SaveFolder', ...
            'Folder used for relative path resolution and output saving.', ...
            'kind', 'input', 'position', 2, 'callType', 'positional', 'isData', false);

        info = PipelineManager.addInput(info, 'ROImasks_filename', 'parameter', ...
            'ROI file name or full path.', ...
            'kind', 'parameter', 'default', 'myROI.roimsk', 'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'CorrAlgorithm', 'parameter', ...
            'ROI correlation algorithm.', ...
            'kind', 'parameter', 'default', 'centroid_vs_centroid', ...
            'allowed', {'centroid_vs_centroid','centroid_vs_agg','avg_vs_avg'}, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'SpatialAggFcn', 'parameter', ...
            'Spatial aggregation for centroid_vs_agg.', ...
            'kind', 'parameter', 'default', 'mean', ...
            'allowed', {'mean','max','min','median'}, 'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'b_FisherZ_transform', 'parameter', ...
            'Apply Fisher Z transform to correlation values.', ...
            'kind', 'parameter', 'default', false, 'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'b_genSPCMaps', 'parameter', ...
            'Generate and save SPC maps as a second UMT output file.', ...
            'kind', 'parameter', 'default', false, 'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'object', 'parameter', ...
            'Optional legacy Acquisition/Modality object for ROI lookup.', ...
            'kind', 'parameter', 'default', [], 'callType', 'namevalue');

        info = PipelineManager.addOutput(info, 'outFile', 'ProcessedData', ...
            'file', 'Generated UMT file manifest saved in SaveFolder.', ...
            default_Output, 1, 'isData', true, 'saveFileName', '');
    end
end

function [value, dimNames] = iResolveImageInput(data, SaveFolder)
%IRESOLVEIMAGEINPUT Resolve supported image input forms.

if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, mfilename, 'data');
    value = single(data);
    dimNames = {'Y','X','T'};
    return
end

if ischar(data) || (isstring(data) && isscalar(data))
    dataFile = char(string(data));
    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('Umitoolbox:genCorrelationMatrix:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);
    switch ext
        case '.dat'
            value = single(loadData(dataFile));
            dimNames = {'Y','X','T'};
            return
        case '.umt'
            data = loadData(dataFile);
        otherwise
            error('Umitoolbox:genCorrelationMatrix:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

assert(isstruct(data) && isscalar(data), ...
    'Umitoolbox:genCorrelationMatrix:UnsupportedInputType', ...
    'Unsupported input type for genCorrelationMatrix.');
validateUMTStruct(data, 'requireEventInfo', false);
assert(strcmpi(char(string(data.kind)), 'image'), ...
    'Umitoolbox:genCorrelationMatrix:InvalidUMTKind', ...
    'Input UMT must have kind = "image".');

entryNames = fieldnames(data.data);
assert(~isempty(entryNames), 'Umitoolbox:genCorrelationMatrix:EmptyUMTData', ...
    'Input UMT contains no image entries.');
entry = data.data.(entryNames{1});
value = single(entry.value);
dimNames = cellstr(string(entry.dimNames));
end

function [centroidList, roiMasks] = iExtractROIGeometry(roiData)
%IEXTRACTROIGEOMETRY Extract ROI masks and centroid pixels.

nROI = numel(roiData.ROI_info);
roiMasks = cell(nROI,1);
centroidList = zeros(nROI,1);

for iROI = 1:nROI
    roiMasks{iROI} = roiData.ROI_info(iROI).Stats.ROI_binary_mask(:);
    cIdx = find(bwmorph(roiData.ROI_info(iROI).Stats.ROI_binary_mask, 'shrink', Inf));
    centroidList(iROI) = cIdx(1);
end
end

function B = iComputeCorrelationMatrix(data, roiMasks, centroidList, corrAlgorithm, spatialAggFcn)
%ICOMPUTECORRELATIONMATRIX Compute ROI correlation matrix.

[nY, nX, nT] = size(data);
data2D = reshape(single(data), nY*nX, nT);
nROI = numel(roiMasks);

switch corrAlgorithm
    case 'centroid_vs_centroid'
        roiVals = data2D(centroidList, :);
        B = corrcoef(roiVals');

    case 'avg_vs_avg'
        roiVals = zeros(nROI, nT, 'single');
        for iROI = 1:nROI
            roiVals(iROI,:) = mean(data2D(roiMasks{iROI}, :), 1, 'omitnan');
        end
        B = corrcoef(roiVals');

    case 'centroid_vs_agg'
        B = zeros(nROI, nROI, 'single');
        for iROI = 1:nROI
            source = data2D(centroidList(iROI), :);
            for jROI = 1:nROI
                target = data2D(roiMasks{jROI}, :);
                rhoVals = zeros(1, size(target,1), 'single');
                for k = 1:size(target,1)
                    rhoTmp = corrcoef(source, target(k,:));
                    rhoVals(k) = rhoTmp(1,2);
                end
                switch spatialAggFcn
                    case 'mean'
                        B(iROI,jROI) = mean(rhoVals, 'omitnan');
                    case 'median'
                        B(iROI,jROI) = median(rhoVals, 'omitnan');
                    case 'min'
                        B(iROI,jROI) = min(rhoVals, [], 'omitnan');
                    case 'max'
                        B(iROI,jROI) = max(rhoVals, [], 'omitnan');
                end
            end
        end
end

B = single(B);
end

function SPCMaps = iComputeSPCMaps(data, centroidList)
%ICOMPUTESPCMAPS Compute seed-pixel correlation maps.

[nY, nX, ~] = size(data);
data2D = reshape(single(data), nY*nX, []);
nROI = numel(centroidList);
SPCMaps = cell(nROI,1);

for iROI = 1:nROI
    seed = data2D(centroidList(iROI), :);
    tmpOut = zeros(1, size(data2D,1), 'single');
    for j = 1:size(data2D,1)
        rho = corrcoef(seed, data2D(j,:));
        tmpOut(j) = rho(1,2);
    end
    SPCMaps{iROI} = reshape(tmpOut, nY, nX);
end
end

function out = iZFisherTruncated(data)
%IZFISHERTRUNCATED Truncate rho before Fisher Z transform.

lim = 0.998;
data(data < -lim) = -lim;
data(data > lim) = lim;
out = atanh(data);
end
