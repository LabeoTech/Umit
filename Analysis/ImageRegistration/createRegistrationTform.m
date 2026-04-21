function outFile = createRegistrationTform(SaveFolder, varargin)
%CREATEREGISTRATIONTFORM Estimate and store folder registration transform.
%
%   outFile = createRegistrationTform(SaveFolder)
%   outFile = createRegistrationTform(SaveFolder, 'Name', Value, ...)
%
%   This function estimates a 2D geometric transform that aligns a moving
%   image from the folder to the folder reference image stored in
%   DataParams.registration.referenceImage. The estimated transform and its
%   associated quality-control metadata are stored in
%   DataParams.registration.
%
%   Inputs:
%       SaveFolder - Folder containing DataParams.mat and the .dat files
%                    associated with the acquisition.
%
%   Name-Value parameters:
%       UseFile      - 'auto' or a .dat filename. Default: 'auto'
%                      When 'auto' is used, the function tries to use
%                      DataParams.registration.referenceFile as the moving
%                      file source.
%       DataParamsFile - Name of the DataParams MAT-file stored in the
%                        folder. Default: 'DataParams.mat'
%       ShowFigure   - Logical scalar. If true, the QC figure remains
%                      visible after creation. Default: true
%
%   Output:
%       outFile      - File manifest containing the updated DataParams MAT
%                      file and the generated QC .fig and .png files.
%
%   Notes:
%       - This function does not modify any .dat file.
%       - It only estimates the transform and updates DataParams.
%       - The resulting transform should be manually reviewed before a
%         later destructive application to the folder data.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'SaveFolder', @isfolder);
addParameter(p, 'UseFile', 'auto', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'DataParamsFile', 'DataParams.mat', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ShowFigure', true, @(x) islogical(x) && isscalar(x));
parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
useFile = char(string(p.Results.UseFile));
dataParamsFile = char(string(p.Results.DataParamsFile));
showFigure = p.Results.ShowFigure;

if ~endsWith(dataParamsFile, '.mat', 'IgnoreCase', true)
    dataParamsFile = [dataParamsFile '.mat'];
end

dataParamsPath = fullfile(SaveFolder, dataParamsFile);
assert(isfile(dataParamsPath), ...
    'Umitoolbox:createRegistrationTform:MissingDataParams', ...
    'DataParams file not found: "%s".', dataParamsPath);

S = load(dataParamsPath, 'DataParams');
assert(isfield(S, 'DataParams') && isstruct(S.DataParams) && isscalar(S.DataParams), ...
    'Umitoolbox:createRegistrationTform:InvalidDataParams', ...
    'DataParams file does not contain a valid scalar DataParams struct.');
DataParams = S.DataParams;
validateDataParams(DataParams);

assert(isfield(DataParams, 'registration') && isstruct(DataParams.registration), ...
    'Umitoolbox:createRegistrationTform:MissingRegistration', ...
    'DataParams.registration is missing or invalid.');
assert(isfield(DataParams.registration, 'referenceImage') && ~isempty(DataParams.registration.referenceImage), ...
    'Umitoolbox:createRegistrationTform:MissingReferenceImage', ...
    'DataParams.registration.referenceImage must be populated before estimating registration.');

refFr = single(DataParams.registration.referenceImage);
assert(ndims(refFr) == 2, ...
    'Umitoolbox:createRegistrationTform:InvalidReferenceImage', ...
    'DataParams.registration.referenceImage must be a 2D image.');

[targetFr, movingFileUsed] = iResolveMovingFrame(SaveFolder, useFile, DataParams);

radius = 0.05 * max(size(refFr));
refFrMask = imgaussfilt(refFr, 0.5) - imgaussfilt(refFr, radius);
targetFrMask = imgaussfilt(targetFr, 0.5) - imgaussfilt(targetFr, radius);

refFrMask(isnan(refFrMask)) = 0;
targetFrMask(isnan(targetFrMask)) = 0;

Rfixed = imref2d(size(refFr));
tformInit = [];

try
    tformInit = imregcorr(targetFrMask, refFrMask, 'similarity', 'Window', true);
catch
end

MIBefore = iComputeMI(refFrMask, targetFrMask);
if ~isempty(tformInit)
    tmpFr = imwarp(targetFrMask, tformInit, 'nearest', 'OutputView', Rfixed);
    MIInit = iComputeMI(refFrMask, tmpFr);
    if MIInit <= MIBefore
        tformInit = [];
    end
end

GF = [1.10, 1.05, 1.02, 1.01];
Eps = [1e-10, 1e-15, 1e-20, 1e-25];
IR = [6.25e-3, 6.25e-5, 6.25e-8, 6.25e-10];
MaxIter = 1000;
bestMI = -Inf;
bestIdx = 1;
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = MaxIter;

for i = 1:4
    optimizer.GrowthFactor = GF(i);
    optimizer.Epsilon = Eps(i);
    optimizer.InitialRadius = IR(i);

    if ~isempty(tformInit)
        tmpFr = imregister(targetFrMask, imref2d(size(targetFrMask)), refFrMask, ...
            imref2d(size(refFrMask)), 'similarity', optimizer, metric, ...
            'DisplayOptimization', false, 'InitialTransformation', tformInit);
    else
        tmpFr = imregister(targetFrMask, imref2d(size(targetFrMask)), refFrMask, ...
            imref2d(size(refFrMask)), 'similarity', optimizer, metric, ...
            'DisplayOptimization', false);
    end

    tmpMI = iComputeMI(refFrMask, tmpFr, metric.NumberOfHistogramBins);
    if tmpMI > bestMI
        bestMI = tmpMI;
        bestIdx = i;
    else
        break
    end
end

optimizer.GrowthFactor = GF(bestIdx);
optimizer.Epsilon = Eps(bestIdx);
optimizer.InitialRadius = IR(bestIdx);

if ~isempty(tformInit)
    tform = imregtform(targetFrMask, imref2d(size(targetFrMask)), refFrMask, ...
        imref2d(size(refFrMask)), 'similarity', optimizer, metric, ...
        'InitialTransformation', tformInit);
else
    tform = imregtform(targetFrMask, imref2d(size(targetFrMask)), refFrMask, ...
        imref2d(size(refFrMask)), 'similarity', optimizer, metric);
end

registeredMask = imwarp(targetFrMask, tform, 'nearest', 'OutputView', Rfixed);
MIAfter = iComputeMI(refFrMask, registeredMask, metric.NumberOfHistogramBins);
MIDelta = MIAfter - MIBefore;

[translationXY, rotationDeg, scaleXY, detVal] = iExtractTransformMetrics(tform);
[qcStatus, qcWarning] = iSummarizeQC(MIDelta, translationXY, rotationDeg, scaleXY, size(refFr));

ts = datestr(now, 'yyyy-mm-dd_HHMMSS');
figName = sprintf('registrationQC_%s.fig', ts);
pngName = sprintf('registrationQC_%s.png', ts);
figPath = fullfile(SaveFolder, figName);
pngPath = fullfile(SaveFolder, pngName);

fig = figure('Name', ['Registration QC - ' ts], 'Visible', 'off', 'Tag', 'registrationQC');
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile; imagesc(refFrMask); axis image off; colormap gray; title('Reference');
nexttile; imagesc(targetFrMask); axis image off; colormap gray; title('Moving');
nexttile([1 2]); imagesc(registeredMask); axis image off; colormap gray; title('Registered');
nexttile([1 2]); imshowpair(refFrMask, registeredMask); title('Overlay: reference vs registered');
sgtitle(fig, sprintf('QC status: %s', qcStatus), 'Interpreter', 'none');
savefig(fig, figPath);
exportgraphics(fig, pngPath, 'Resolution', 150);
if showFigure
    fig.Visible = 'on';
else
    close(fig);
end

DataParams.registration.isRegistered = false;
DataParams.registration.isReviewed = false;
DataParams.registration.tform = tform;
DataParams.registration.transformType = 'similarity';
DataParams.registration.method = 'imregcorr_init_plus_multimodal_imregtform';
DataParams.registration.referenceDescription = 'Folder reference image from DataParams.registration.referenceImage';
DataParams.registration.referenceFile = movingFileUsed;
DataParams.registration.referenceImage = refFr;
DataParams.registration.createdOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');
DataParams.registration.source = mfilename;
DataParams.registration.qcStatus = qcStatus;
DataParams.registration.qcWarning = qcWarning;
DataParams.registration.qcMetrics.MIBefore = MIBefore;
DataParams.registration.qcMetrics.MIAfter = MIAfter;
DataParams.registration.qcMetrics.MIDelta = MIDelta;
DataParams.registration.qcMetrics.translationXY_px = translationXY;
DataParams.registration.qcMetrics.rotationDeg = rotationDeg;
DataParams.registration.qcMetrics.scaleXY = scaleXY;
DataParams.registration.qcMetrics.determinant = detVal;
DataParams.registration.qcFigureFile = figName;
DataParams.registration.qcPreviewImageFile = pngName;
DataParams.registration.appliedOn = '';
DataParams.registration.appliedBy = '';
DataParams.registration.confirmationMode = '';
validateDataParams(DataParams);
save(dataParamsPath, 'DataParams');

outFile = {dataParamsFile, figName, pngName};
end

function [targetFr, movingFileUsed] = iResolveMovingFrame(SaveFolder, useFile, DataParams)
%IRESOLVEMOVINGFRAME Resolve the moving image used to estimate the tform.

movingFileUsed = '';

switch lower(useFile)
    case 'auto'
        assert(isfield(DataParams.registration, 'referenceFile') && ~isempty(DataParams.registration.referenceFile), ...
            'Umitoolbox:createRegistrationTform:MissingReferenceFile', ...
            ['UseFile="auto" requires DataParams.registration.referenceFile ' ...
             'to be populated.']);
        movingFileUsed = char(string(DataParams.registration.referenceFile));
    otherwise
        movingFileUsed = char(string(useFile));
end

if ~endsWith(movingFileUsed, '.dat', 'IgnoreCase', true)
    movingFileUsed = [movingFileUsed '.dat'];
end

movingPath = fullfile(SaveFolder, movingFileUsed);
assert(isfile(movingPath), ...
    'Umitoolbox:createRegistrationTform:MissingMovingFile', ...
    'Moving .dat file not found: "%s".', movingPath);

md = loadMetaData(movingPath);
assert(isfield(md, 'Height') && isfield(md, 'Width'), ...
    'Umitoolbox:createRegistrationTform:InvalidMovingMetadata', ...
    'Could not resolve Height/Width from moving file metadata.');

fid = fopen(movingPath, 'r');
assert(fid ~= -1, ...
    'Umitoolbox:createRegistrationTform:FileOpenFailed', ...
    'Could not open moving file "%s".', movingPath);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

firstFrame = fread(fid, double(md.Height) * double(md.Width), '*single');
assert(numel(firstFrame) == double(md.Height) * double(md.Width), ...
    'Umitoolbox:createRegistrationTform:InvalidMovingFile', ...
    'Could not read the first frame from "%s".', movingPath);

targetFr = reshape(firstFrame, md.Height, md.Width);
end

function MI = iComputeMI(A, B, nBins)
%ICOMPUTEMI Estimate mutual information from 2D histograms.

if nargin < 3 || isempty(nBins)
    nBins = 50;
end

if any(size(A) ~= size(B))
    B = imresize(B, size(A));
end

counts = histcounts2(A(:), B(:), nBins);
Pxy = counts / sum(counts(:));
Px = sum(Pxy, 2);
Py = sum(Pxy, 1);
PxPy = Px * Py;
idx = Pxy > 0 & PxPy > 0;
MI = sum(Pxy(idx) .* log2(Pxy(idx) ./ PxPy(idx)), 'all');
end

function [translationXY, rotationDeg, scaleXY, detVal] = iExtractTransformMetrics(tform)
%IEXTRACTTRANSFORMMETRICS Extract simple QC metrics from a 2D transform.

if isprop(tform, 'A')
    T = tform.A;
elseif isprop(tform, 'T')
    T = tform.T;
else
    translationXY = [NaN NaN];
    rotationDeg = NaN;
    scaleXY = [NaN NaN];
    detVal = NaN;
    return
end

A = double(T(1:2,1:2));
translationXY = double(T(3,1:2));
rotationDeg = atan2d(A(2,1), A(1,1));
scaleXY = [norm(A(:,1)) norm(A(:,2))];
detVal = det(A);
end

function [qcStatus, qcWarning] = iSummarizeQC(MIDelta, translationXY, rotationDeg, scaleXY, imageSize)
%ISUMMARIZEQC Generate a simple QC status and warning summary.

warnings = {};
transMag = norm(translationXY);
maxDim = max(imageSize);

if MIDelta <= 0
    warnings{end+1} = 'Mutual information did not improve.'; %#ok<AGROW>
end
if transMag > 0.25 * maxDim
    warnings{end+1} = 'Large translation relative to frame size.'; %#ok<AGROW>
end
if abs(rotationDeg) > 15
    warnings{end+1} = 'Large rotation estimate.'; %#ok<AGROW>
end
if any(abs(scaleXY - 1) > 0.05)
    warnings{end+1} = 'Scale deviates substantially from unity.'; %#ok<AGROW>
end

if isempty(warnings)
    if MIDelta > 0
        qcStatus = 'estimated_not_reviewed';
        qcWarning = '';
    else
        qcStatus = 'review_recommended';
        qcWarning = 'Registration metrics are inconclusive. Visual review is strongly recommended.';
    end
else
    qcStatus = 'warning_large_transform';
    qcWarning = strjoin(warnings, ' ');
end
end
