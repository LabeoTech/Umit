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
%       ROImasks_filename    - UMIT .roi file name or full path. A bare
%                              filename is resolved inside SaveFolder.
%                              Default: 'myROI.roi'
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
%
%   Output:
%       outFile - File manifest cell array containing the generated UMT
%                 file name(s) saved in SaveFolder.
%
%   Notes:
%       - ROI files are read through loadROIFile(...), which migrates and
%         validates the current UMIT .roi schema. Pre-.roi ROI files are
%         not supported.
%       - The main output is a roi UMT correlation matrix file.
%       - SPC maps are saved as a second image UMT file when requested.
%       - Traces that are entirely NaN (masked pixels) are reported as NaN
%         and do not propagate into the coefficients of the other ROIs.
%         Partially masked traces use the pairwise-complete estimator, which
%         requires the Statistics and Machine Learning Toolbox.


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
addParameter(p, 'ROImasks_filename', 'myROI.roi', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'CorrAlgorithm', 'centroid_vs_centroid', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'centroid_vs_centroid','centroid_vs_agg','avg_vs_avg'}));
addParameter(p, 'SpatialAggFcn', 'mean', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'mean','max','min','median'}));
addParameter(p, 'b_FisherZ_transform', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'b_genSPCMaps', false, @(x) islogical(x) && isscalar(x));

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
roiFile = char(string(p.Results.ROImasks_filename));
corrAlgorithm = lower(char(string(p.Results.CorrAlgorithm)));
spatialAggFcn = lower(char(string(p.Results.SpatialAggFcn)));
bFisherZ = p.Results.b_FisherZ_transform;
bGenSPCMaps = p.Results.b_genSPCMaps;

assert(isfolder(SaveFolder), 'Umitoolbox:genCorrelationMatrix:InvalidSaveFolder', ...
    'SaveFolder "%s" does not exist.', SaveFolder);

roiSet = iLoadROISet(roiFile, SaveFolder);

[value, dimNames] = iResolveImageInput(data, SaveFolder);
assert(isequal(dimNames, {'Y','X','T'}), ...
    'Umitoolbox:genCorrelationMatrix:WrongFormat', ...
    'Input data must be an Image time series with dimensions {''Y'',''X'',''T''}.');
assert(isequal([size(value,1), size(value,2)], roiSet.imageSizeYX), ...
    'Umitoolbox:genCorrelationMatrix:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

roiNames = roiSet.names;
[centroidList, roiMasks] = iExtractROIGeometry(roiSet);

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
    outFile = [outFile, {spcFile}];
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
            'UMIT .roi file name or full path, resolved inside SaveFolder.', ...
            'kind', 'parameter', 'default', 'myROI.roi', 'callType', 'namevalue');

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

function [centroidList, roiMasks] = iExtractROIGeometry(roiSet)
%IEXTRACTROIGEOMETRY Extract ROI masks and centroid pixels.

nROI = numel(roiSet.masks);
roiMasks = cell(nROI,1);
centroidList = zeros(nROI,1);

for iROI = 1:nROI
    thisMask = roiSet.masks{iROI};
    roiMasks{iROI} = thisMask(:);

    cIdx = find(bwmorph(thisMask, 'shrink', Inf));
    if isempty(cIdx)
        % bwmorph can shrink a thin or empty ROI to nothing. Fall back to
        % the first mask pixel so a degenerate ROI cannot abort the run
        % with a bare index error.
        cIdx = find(thisMask);
    end
    if isempty(cIdx)
        error('Umitoolbox:genCorrelationMatrix:EmptyROIMask', ...
            'ROI "%s" has an empty mask.', roiSet.names{iROI});
    end
    centroidList(iROI) = cIdx(1);
end
end

% =========================================================================
% Local helper: load and normalize a UMIT .roi file
% =========================================================================
function roiSet = iLoadROISet(roiFile, SaveFolder)
%ILOADROISET Resolve and load a UMIT .roi file into names, masks, and size.

roiFile = char(string(roiFile));
if ~isfile(roiFile)
    roiFile = fullfile(SaveFolder, roiFile);
end

if ~isfile(roiFile)
    error('Umitoolbox:genCorrelationMatrix:MissingROIFile', ...
        'ROI file was not found: "%s".', roiFile);
end

[~, ~, ext] = fileparts(roiFile);
if ~strcmpi(ext, '.roi')
    error('Umitoolbox:genCorrelationMatrix:UnsupportedROIFile', ...
        ['ROI files must use the UMIT ".roi" format. Pre-.roi ROI files ' ...
         'are not supported. Received: "%s".'], roiFile);
end

% loadROIFile migrates and validates the schema, so masks are guaranteed to
% be 2-D and to match imageInfo.imageSizeYX, and ROI names are unique.
ROIFile = loadROIFile(roiFile);

if isempty(ROIFile.ROIs)
    error('Umitoolbox:genCorrelationMatrix:EmptyROIFile', ...
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

function B = iComputeCorrelationMatrix(data, roiMasks, centroidList, corrAlgorithm, spatialAggFcn)
%ICOMPUTECORRELATIONMATRIX Compute ROI correlation matrix.

[nY, nX, nT] = size(data);
data2D = reshape(single(data), nY*nX, nT);
nROI = numel(roiMasks);

switch corrAlgorithm
    case 'centroid_vs_centroid'
        roiVals = data2D(centroidList, :);
        B = iCorrRows(roiVals, roiVals);

    case 'avg_vs_avg'
        roiVals = zeros(nROI, nT, 'single');
        for iROI = 1:nROI
            roiVals(iROI,:) = mean(data2D(roiMasks{iROI}, :), 1, 'omitnan');
        end
        B = iCorrRows(roiVals, roiVals);

    case 'centroid_vs_agg'
        % One centroid-vs-all-pixels correlation per target ROI, rather than
        % one corrcoef call per (seed, pixel) pair.
        sources = data2D(centroidList, :);
        B = zeros(nROI, nROI, 'single');
        for jROI = 1:nROI
            rhoVals = iCorrRows(sources, data2D(roiMasks{jROI}, :));
            B(:,jROI) = iAggregateRho(rhoVals, spatialAggFcn);
        end
end

B = single(B);
end

function agg = iAggregateRho(rhoVals, spatialAggFcn)
%IAGGREGATERHO Reduce per-pixel correlations along the pixel dimension.

switch spatialAggFcn
    case 'mean'
        agg = mean(rhoVals, 2, 'omitnan');
    case 'median'
        agg = median(rhoVals, 2, 'omitnan');
    case 'min'
        agg = min(rhoVals, [], 2, 'omitnan');
    case 'max'
        agg = max(rhoVals, [], 2, 'omitnan');
    otherwise
        error('Umitoolbox:genCorrelationMatrix:InvalidAggFcn', ...
            'Unknown spatial aggregation function "%s".', spatialAggFcn);
end
end

function R = iCorrRows(A, B)
%ICORRROWS Correlate every row of A against every row of B over time.
%
%   A is m-by-T, B is n-by-T and R is m-by-n. Traces that are entirely NaN
%   -- the normal state of masked pixels after GSR or normalization -- are
%   excluded and reported as NaN instead of poisoning the whole matrix, and
%   partially masked traces fall back to the pairwise-complete estimator.

X = A';
Y = B';

keepX = ~all(isnan(X), 1);
keepY = ~all(isnan(Y), 1);

R = nan(size(X,2), size(Y,2), 'single');
if ~any(keepX) || ~any(keepY)
    return
end

Xk = X(:, keepX);
Yk = Y(:, keepY);

if any(isnan(Xk), 'all') || any(isnan(Yk), 'all')
    % Only some samples of these traces are masked. CORR drops NaN samples
    % per pair, which is slower but is the only correct reduction here.
    Rk = corr(Xk, Yk, 'rows', 'pairwise');
else
    Rk = iFastCorr(Xk, Yk);
end

R(keepX, keepY) = single(Rk);
end

function R = iFastCorr(X, Y)
%IFASTCORR Correlate the columns of two NaN-free matrices as one product.

Xc = X - mean(X, 1);
Yc = Y - mean(Y, 1);

% Zero-variance traces divide by zero and yield NaN, matching CORRCOEF.
Xc = Xc ./ vecnorm(Xc, 2, 1);
Yc = Yc ./ vecnorm(Yc, 2, 1);

R = Xc' * Yc;

% Guard against rounding pushing a coefficient just outside [-1, 1].
R = max(min(R, 1), -1);
end

function SPCMaps = iComputeSPCMaps(data, centroidList)
%ICOMPUTESPCMAPS Compute seed-pixel correlation maps.

[nY, nX, ~] = size(data);
data2D = reshape(single(data), nY*nX, []);
nROI = numel(centroidList);
SPCMaps = cell(nROI,1);

% One seeds-by-all-pixels correlation instead of one corrcoef call per
% (seed, pixel) pair.
rho = iCorrRows(data2D(centroidList, :), data2D);

for iROI = 1:nROI
    SPCMaps{iROI} = reshape(rho(iROI,:), nY, nX);
end
end

function out = iZFisherTruncated(data)
%IZFISHERTRUNCATED Truncate rho before Fisher Z transform.

lim = 0.998;
data(data < -lim) = -lim;
data(data > lim) = lim;
out = atanh(data);
end
