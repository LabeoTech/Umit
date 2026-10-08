function [outData, spcFile] = genCorrelationMatrix(data, SaveFolder, varargin)
%GENCORRELATIONMATRIX Generate ROI correlation matrices from a .dat image time series.
%
%   outData = genCorrelationMatrix(data, SaveFolder)
%   [outData, spcFile] = genCorrelationMatrix(data, SaveFolder, 'Name', Value, ...)
%
%   Supported input:
%       Raw .dat filename storing continuous Y-X-T data. Arrays, UMT
%       structs, .umt files, and event-split (Y-X-T-E) files are not
%       supported.
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
%       outData - UMT struct of kind "roi": entry "CorrMatrix" with
%                 dimensions {'ROI','ROI'} and the ROI names as labels.
%       spcFile - SaveFolder-relative name of the image UMT file holding the
%                 SPC maps ('corrMatrix_SPCMaps.umt', one Y-X entry per
%                 ROI). Empty when b_genSPCMaps is false.
%
%   Notes:
%       - ROI files are read through loadROIFile(...), which migrates and
%         validates the current UMIT .roi schema. Pre-.roi ROI files are
%         not supported.
%       - The .dat file is streamed in X slabs sized for a fixed 128 MB
%         budget; the recording is never loaded whole. 'centroid_vs_centroid'
%         reads only the seed pixels. Without SPC maps, only the columns that
%         contain ROI pixels are read. The SPC maps (4*Y*X*nROI bytes) stay
%         resident in RAM until saved.
%       - Traces that are entirely NaN (masked pixels) are reported as NaN
%         and do not propagate into the coefficients of the other ROIs.
%         Partially masked traces use the pairwise-complete estimator, which
%         requires the Statistics and Machine Learning Toolbox.

default_Output = 'corrMatrix.umt';
spcFileName = 'corrMatrix_SPCMaps.umt';
spcFile = '';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
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

dataFile = iResolveDatFile(data, SaveFolder);
datInfo = loadMetaData(dataFile);
assertDatLayout(datInfo, {{'Y','X','T'}}, 'genCorrelationMatrix');
assert(isequal([datAxisSize(datInfo, 'Y'), datAxisSize(datInfo, 'X')], roiSet.imageSizeYX), ...
    'Umitoolbox:genCorrelationMatrix:IncompatibleSizes', ...
    'Input frame size is different from the frame size in the ROI file.');

roiNames = roiSet.names;
[centroidList, roiMasks] = iExtractROIGeometry(roiSet);

slabIn = spatialSlabIO('open', dataFile, 'Info', datInfo);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

[corrMatrix, spcMaps] = iStreamCorrelations(slabIn, roiMasks, centroidList, ...
    corrAlgorithm, spatialAggFcn, bGenSPCMaps);
if bFisherZ
    corrMatrix = iZFisherTruncated(corrMatrix);
end

labels = struct();
labels.ROI = roiNames(:).';
outData = genUMTStruct(corrMatrix, ...
    'kind', 'roi', ...
    'entryName', 'CorrMatrix', ...
    'dimNames', {'ROI','ROI'}, ...
    'labels', labels);

if bGenSPCMaps
    if bFisherZ
        spcMaps = iZFisherTruncated(spcMaps);
    end

    spcUMT = [];
    for iMap = 1:numel(roiNames)
        entryName = matlab.lang.makeValidName(roiNames{iMap});
        if isempty(entryName)
            entryName = sprintf('ROI_%d', iMap);
        end
        if iMap == 1
            spcUMT = genUMTStruct(spcMaps(:,:,iMap), ...
                'kind', 'image', ...
                'entryName', entryName, ...
                'dimNames', {'Y','X'});
        else
            spcUMT = genUMTStruct(spcUMT, ...
                'value', spcMaps(:,:,iMap), ...
                'entryName', entryName, ...
                'dimNames', {'Y','X'});
        end
    end

    saveData(fullfile(SaveFolder, spcFileName), spcUMT);
    spcFile = spcFileName;
end

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Generate ROI correlation matrices from a .dat image time series.');
        info.version = '1.0.0';

        info = PipelineManager.addInput(info, 'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            'Continuous Y-X-T .dat image time series (file input only).', ...
            'kind', 'input', 'position', 1, 'callType', 'positional', ...
            'isData', true, 'supportsFile', true, 'dataMode', 'file');

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

        info = PipelineManager.addOutput(info, 'outData', 'ProcessedData', ...
            'data', 'ROI correlation matrix UMT (kind roi).', ...
            default_Output, 1, 'isData', true);

        info = PipelineManager.addOutput(info, 'spcFile', 'ProcessedData', ...
            'file', ['Seed-pixel correlation maps (image UMT), saved in ' ...
            'SaveFolder when b_genSPCMaps is true.'], ...
            spcFileName, 2, 'isData', false, 'isRequired', false);
    end
end

function dataFile = iResolveDatFile(data, SaveFolder)
%IRESOLVEDATFILE Resolve the .dat filename input (the only supported form).

if ~(ischar(data) || (isstring(data) && isscalar(data)))
    error('Umitoolbox:genCorrelationMatrix:UnsupportedInputType', ...
        ['Input "data" must be a .dat filename. Arrays, UMT structs, and ' ...
         '.umt files are not supported.']);
end

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
if ~strcmpi(ext, '.dat')
    error('Umitoolbox:genCorrelationMatrix:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', ext);
end
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

% =========================================================================
% Streaming computation
% =========================================================================
function [B, spcMaps] = iStreamCorrelations(slabIn, roiMasks, centroidList, corrAlgorithm, spatialAggFcn, bSPC)
%ISTREAMCORRELATIONS ROI correlation matrix (and SPC maps) from X slabs.
%
%   Reads the seed traces first. 'centroid_vs_centroid' needs nothing else;
%   the other algorithms make one pass over the file, restricted to the
%   columns holding ROI pixels unless SPC maps (all pixels) are requested.
%   B is nROI-by-nROI; spcMaps is Y-by-X-by-nROI (empty without SPC).

Ny = slabIn.Ny;
Nx = slabIn.Nx;
Nt = datAxisSize(slabIn.Info, 'T');
nROI = numel(roiMasks);

[cy, cx] = ind2sub([Ny, Nx], centroidList);
seeds = iReadSeedTraces(slabIn, cy, cx, Nt);

if bSPC
    spcMaps = nan(Ny, Nx, nROI, 'single');
else
    spcMaps = zeros(Ny, Nx, 0, 'single');
end
needPass = bSPC || any(strcmp(corrAlgorithm, {'avg_vs_avg', 'centroid_vs_agg'}));

if needPass
    roiMask2D = cellfun(@(m) reshape(m, Ny, Nx), roiMasks, 'UniformOutput', false);

    if bSPC
        xList = 1:Nx;
    else
        xList = find(any(cat(3, roiMask2D{:}), [1 3]));
    end

    % Slab width for a fixed byte budget. The factor covers the slab, the
    % extracted pixel traces, and the temporary copies made by iCorrRows.
    slabBudgetBytes = 128 * 1024 * 1024;
    bytesPerX = Ny * Nt * 4 * 6;
    xPerSlab = max(1, floor(slabBudgetBytes / bytesPerX));

    traceSum = zeros(nROI, Nt);
    traceCount = zeros(nROI, Nt);
    rhoByROI = cell(nROI, 1);
    rhoFilled = zeros(nROI, 1);
    if strcmp(corrAlgorithm, 'centroid_vs_agg')
        for iROI = 1:nROI
            rhoByROI{iROI} = nan(nROI, nnz(roiMasks{iROI}), 'single');
        end
    end

    for k = 1:xPerSlab:numel(xList)
        xIdx = xList(k:min(k + xPerSlab - 1, numel(xList)));
        slab2D = reshape(single(spatialSlabIO('read', slabIn, xIdx)), ...
            Ny * numel(xIdx), Nt);

        if bSPC
            rho = iCorrRows(seeds, slab2D);
            spcMaps(:, xIdx, :) = permute(reshape(rho, nROI, Ny, numel(xIdx)), [2 3 1]);
        end

        for iROI = 1:nROI
            inSlab = find(reshape(roiMask2D{iROI}(:, xIdx), [], 1));
            if isempty(inSlab)
                continue
            end
            traces = slab2D(inSlab, :);

            switch corrAlgorithm
                case 'avg_vs_avg'
                    traceSum(iROI, :) = traceSum(iROI, :) + sum(double(traces), 1, 'omitnan');
                    traceCount(iROI, :) = traceCount(iROI, :) + sum(~isnan(traces), 1);

                case 'centroid_vs_agg'
                    cols = rhoFilled(iROI) + (1:numel(inSlab));
                    rhoByROI{iROI}(:, cols) = iCorrRows(seeds, traces);
                    rhoFilled(iROI) = rhoFilled(iROI) + numel(inSlab);
            end
        end
    end
end

switch corrAlgorithm
    case 'centroid_vs_centroid'
        B = iCorrRows(seeds, seeds);

    case 'avg_vs_avg'
        roiVals = single(traceSum ./ traceCount);
        B = iCorrRows(roiVals, roiVals);

    case 'centroid_vs_agg'
        % One centroid-vs-all-pixels correlation per target ROI, rather than
        % one corrcoef call per (seed, pixel) pair.
        B = zeros(nROI, nROI, 'single');
        for jROI = 1:nROI
            B(:, jROI) = iAggregateRho(rhoByROI{jROI}, spatialAggFcn);
        end
end

B = single(B);
end

function seeds = iReadSeedTraces(slabIn, cy, cx, Nt)
%IREADSEEDTRACES Time traces of the centroid pixels (nROI-by-T).

xUnique = unique(cx(:)).';
block = single(spatialSlabIO('read', slabIn, xUnique));

seeds = zeros(numel(cy), Nt, 'single');
for iROI = 1:numel(cy)
    seeds(iROI, :) = reshape(block(cy(iROI), xUnique == cx(iROI), :), 1, Nt);
end
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

function out = iZFisherTruncated(data)
%IZFISHERTRUNCATED Truncate rho before Fisher Z transform.

% atanh(+-1) is +-Inf, so rho is clamped just inside the unit interval before
% the transform; 0.998 is an empirical margin, not a derived bound.
lim = 0.998;
data(data < -lim) = -lim;
data(data > lim) = lim;
out = atanh(data);
end
