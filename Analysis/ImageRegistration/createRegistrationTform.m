function outFile = createRegistrationTform(SaveFolder, varargin)
%CREATEREGISTRATIONTFORM Estimate and store folder registration transform.
%
%   outFile = createRegistrationTform(SaveFolder)
%   outFile = createRegistrationTform(SaveFolder, 'Name', Value, ...)
%   info    = createRegistrationTform('pipelineInfo')
%
%   This function estimates a 2D geometric transform that aligns one image
%   from the current folder to the reference image stored in
%   DataParams.registration.referenceImage. The estimated transform and its
%   associated QC metadata are stored back into DataParams.mat.
%
%   This function does NOT modify any .dat files. It only estimates and
%   stores the transform, generates QC artifacts, and updates
%   DataParams.registration so the user can inspect the result before
%   deciding whether to apply the registration to the whole folder.
%
%   Inputs:
%       SaveFolder - Folder containing DataParams.mat, AcqInfos.mat, and
%                    one or more .dat files.
%
%   Name-Value parameters:
%       UseFile    - File used to build the moving image. Use:
%                      'auto'  : first eligible .dat file found in folder
%                      <name>  : specific .dat file name
%                    Default: 'auto'
%
%       ShowFigure - Logical scalar. If true, show the QC figure after
%                    creation. If false, create and save the figure
%                    invisibly.
%                    Default: true
%
%       RefStatistic - Statistic used to collapse a YXT file into a 2D
%                    moving image when needed.
%                    Allowed: 'mean', 'median', 'first'
%                    Default: 'mean'
%
%   Output:
%       outFile    - File manifest of generated/updated files:
%                    {DataParams.mat, QC .fig, QC .png}
%
%   Notes:
%       - The reference image must already be stored in
%         DataParams.registration.referenceImage.
%       - QC metrics are advisory only and do not replace visual review.
%       - The QC artifacts are timestamped and saved in SaveFolder.
%
%   See also:
%       applyRegistrationTformOnFolder
%
%   Examples:
%       createRegistrationTform('D:\Data\Mouse01\Session01');
%
%       createRegistrationTform( ...
%           'D:\Data\Mouse01\Session01', ...
%           'UseFile', 'green.dat', ...
%           'ShowFigure', false);

% -------------------------------------------------------------------------
% Shared declarations
% -------------------------------------------------------------------------
defaultUseFile = 'auto';
defaultShowFigure = true;
defaultRefStatistic = 'mean';

% -------------------------------------------------------------------------
% pipelineInfo query
% -------------------------------------------------------------------------
if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) ...
        && strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    outFile = localPipelineInfo();
    return
end

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'SaveFolder', @isfolder);
addParameter(p, 'UseFile', defaultUseFile, ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ShowFigure', defaultShowFigure, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'RefStatistic', defaultRefStatistic, ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'mean','median','first'}));

parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
UseFile = char(string(p.Results.UseFile));
ShowFigure = p.Results.ShowFigure;
RefStatistic = lower(char(string(p.Results.RefStatistic)));

clear p

% -------------------------------------------------------------------------
% Load and validate DataParams
% -------------------------------------------------------------------------
dataParamsFile = fullfile(SaveFolder, 'DataParams.mat');
assert(isfile(dataParamsFile), ...
    'Umitoolbox:createRegistrationTform:MissingDataParams', ...
    'DataParams.mat was not found in "%s".', SaveFolder);

tmp = load(dataParamsFile, 'DataParams');
assert(isfield(tmp, 'DataParams'), ...
    'Umitoolbox:createRegistrationTform:MissingDataParamsVariable', ...
    'DataParams.mat does not contain a variable named "DataParams".');

DataParams = tmp.DataParams;
validateDataParams(DataParams);

assert(~isempty(DataParams.registration.referenceImage), ...
    'Umitoolbox:createRegistrationTform:MissingReferenceImage', ...
    ['DataParams.registration.referenceImage is empty. A folder ' ...
     'reference image must be defined before estimating registration.']);

refFr = single(DataParams.registration.referenceImage);
assert(ndims(refFr) == 2, ...
    'Umitoolbox:createRegistrationTform:InvalidReferenceImage', ...
    'DataParams.registration.referenceImage must be a 2D image.');

% -------------------------------------------------------------------------
% Resolve moving-image source file
% -------------------------------------------------------------------------
targetFile = iResolveTargetFile(SaveFolder, UseFile);

% -------------------------------------------------------------------------
% Load moving image from file
% -------------------------------------------------------------------------
targetFr = iBuildMovingImageFromDat(targetFile, RefStatistic);

assert(~isempty(targetFr) && isnumeric(targetFr) && ismatrix(targetFr), ...
    'Umitoolbox:createRegistrationTform:InvalidMovingImage', ...
    'Failed to build a valid 2D moving image from "%s".', targetFile);

targetFr = single(targetFr);

% Resize moving image if needed so registration steps are well-defined.
if any(size(targetFr) ~= size(refFr))
    targetFr = imresize(targetFr, size(refFr));
end

% -------------------------------------------------------------------------
% Build filtered images for registration
% -------------------------------------------------------------------------
radius = 0.05 * max(size(refFr));
radius = max(radius, 1);

refFrMask = imgaussfilt(refFr, 0.5) - imgaussfilt(refFr, radius);
targetFrMask = imgaussfilt(targetFr, 0.5) - imgaussfilt(targetFr, radius);

refFrMask(~isfinite(refFrMask)) = 0;
targetFrMask(~isfinite(targetFrMask)) = 0;

Rfixed = imref2d(size(refFrMask));

% -------------------------------------------------------------------------
% Initial phase-correlation transform
% -------------------------------------------------------------------------
tformInit = imregcorr(targetFrMask, refFrMask, 'similarity', 'Window', true);

counts = histcounts2(refFrMask(:), targetFrMask(:), 50);
MIBefore = iMutualInformation(counts);

tmpRegistered = imwarp(targetFrMask, tformInit, 'nearest', 'OutputView', Rfixed);
counts = histcounts2(refFrMask(:), tmpRegistered(:), 50);
MIAfterInit = iMutualInformation(counts);

if MIAfterInit < MIBefore
    tformInit = [];
end

% -------------------------------------------------------------------------
% Optimize registration parameters
% -------------------------------------------------------------------------
GF = [1.10, 1.05, 1.02, 1.01];
Eps = [1e-10, 1e-15, 1e-20, 1e-25];
IR = [6.25e-3, 6.25e-5, 6.25e-8, 6.25e-10];
MaxIter = 1000;

[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = MaxIter;

bestMI = -inf;
bestIdx = 1;

for i = 1:numel(GF)
    optimizer.GrowthFactor = GF(i);
    optimizer.Epsilon = Eps(i);
    optimizer.InitialRadius = IR(i);

    if ~isempty(tformInit)
        tmpFr = imregister( ...
            targetFrMask, imref2d(size(targetFrMask)), ...
            refFrMask, imref2d(size(refFrMask)), ...
            'similarity', optimizer, metric, ...
            'DisplayOptimization', false, ...
            'InitialTransformation', tformInit);
    else
        tmpFr = imregister( ...
            targetFrMask, imref2d(size(targetFrMask)), ...
            refFrMask, imref2d(size(refFrMask)), ...
            'similarity', optimizer, metric, ...
            'DisplayOptimization', false);
    end

    counts = histcounts2(refFrMask(:), tmpFr(:), metric.NumberOfHistogramBins);
    tmpMI = iMutualInformation(counts);

    if tmpMI > bestMI
        bestMI = tmpMI;
        bestIdx = i;
    end
end

optimizer.GrowthFactor = GF(bestIdx);
optimizer.Epsilon = Eps(bestIdx);
optimizer.InitialRadius = IR(bestIdx);

if ~isempty(tformInit)
    tform = imregtform( ...
        targetFrMask, imref2d(size(targetFrMask)), ...
        refFrMask, imref2d(size(refFrMask)), ...
        'similarity', optimizer, metric, ...
        'InitialTransformation', tformInit);
else
    tform = imregtform( ...
        targetFrMask, imref2d(size(targetFrMask)), ...
        refFrMask, imref2d(size(refFrMask)), ...
        'similarity', optimizer, metric);
end

registeredMask = imwarp(targetFrMask, tform, 'nearest', 'OutputView', Rfixed);
counts = histcounts2(refFrMask(:), registeredMask(:), metric.NumberOfHistogramBins);
MIAfter = iMutualInformation(counts);
MIDelta = MIAfter - MIBefore;

% -------------------------------------------------------------------------
% Derive QC metrics
% -------------------------------------------------------------------------
[translationXY, rotationDeg, scaleXY, determinantVal] = iExtractTformMetrics(tform);

[qcStatus, qcWarning] = iBuildQCStatus( ...
    MIDelta, translationXY, rotationDeg, scaleXY, size(refFrMask));

% -------------------------------------------------------------------------
% Save QC artifacts
% -------------------------------------------------------------------------
ts = datestr(now, 'yyyy-mm-dd_HHMMSS');
qcFigFile = fullfile(SaveFolder, ['registrationQC_' ts '.fig']);
qcPngFile = fullfile(SaveFolder, ['registrationQC_' ts '.png']);

fig = figure('Name', ['Registration QC - ' ts], 'Visible', 'off', 'Tag', 'registrationQC');
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(refFrMask);
axis image off;
colormap gray;
title('Reference');

nexttile;
imagesc(targetFrMask);
axis image off;
colormap gray;
title('Moving');

nexttile([1 2]);
imagesc(registeredMask);
axis image off;
colormap gray;
title('Registered');

nexttile([1 2]);
imshowpair(refFrMask, registeredMask);
title('Overlay: reference vs registered');

savefig(fig, qcFigFile);
exportgraphics(fig, qcPngFile, 'Resolution', 150);

if ShowFigure
    fig.Visible = 'on';
else
    close(fig);
end

% -------------------------------------------------------------------------
% Update DataParams.registration
% -------------------------------------------------------------------------
DataParams.registration.isRegistered = false;
DataParams.registration.isReviewed = false;
DataParams.registration.tform = tform;
DataParams.registration.transformType = 'similarity';
DataParams.registration.method = ...
    'imregcorr initialization + multimodal imregtform optimization';
DataParams.registration.referenceDescription = ...
    'Folder reference image stored in DataParams.registration.referenceImage';
[~, usedName, usedExt] = fileparts(targetFile);
DataParams.registration.referenceFile = [usedName usedExt];
DataParams.registration.referenceImage = refFr;
DataParams.registration.createdOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');
DataParams.registration.source = mfilename;
DataParams.registration.qcStatus = qcStatus;
DataParams.registration.qcWarning = qcWarning;
DataParams.registration.qcMetrics = struct( ...
    'MIBefore', MIBefore, ...
    'MIAfter', MIAfter, ...
    'MIDelta', MIDelta, ...
    'translationXY_px', translationXY, ...
    'rotationDeg', rotationDeg, ...
    'scaleXY', scaleXY, ...
    'determinant', determinantVal);
DataParams.registration.qcFigureFile = qcFigFile;
DataParams.registration.qcPreviewImageFile = qcPngFile;
DataParams.registration.appliedOn = '';
DataParams.registration.appliedBy = '';
DataParams.registration.confirmationMode = '';
DataParams.registration.notes = '';

DataParams.lastModified = datestr(now, 'yyyy-mm-dd HH:MM:SS');

validateDataParams(DataParams);
save(dataParamsFile, 'DataParams');

outFile = {dataParamsFile, qcFigFile, qcPngFile};

function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for this function.

info = PipelineManager.createPipelineInfo( ...
    mfilename, ...
    ['Estimate and store a folder registration transform using the ' ...
     'reference image in DataParams.']);

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    ['Folder containing DataParams.mat, AcqInfos.mat, and the source ' ...
     '.dat files.'], ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'UseFile', ...
    'parameter', ...
    'Specific .dat file used to build the moving image, or ''auto''.', ...
    'kind', 'parameter', ...
    'default', defaultUseFile, ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'ShowFigure', ...
    'parameter', ...
    'If true, show the QC figure after creation.', ...
    'kind', 'parameter', ...
    'default', defaultShowFigure, ...
    'allowed', [true false], ...
    'callType', 'namevalue');

info = PipelineManager.addInput( ...
    info, ...
    'RefStatistic', ...
    'parameter', ...
    'Statistic used to build a 2D moving image from a YXT file.', ...
    'kind', 'parameter', ...
    'default', defaultRefStatistic, ...
    'allowed', {'mean','median','first'}, ...
    'callType', 'namevalue');

info = PipelineManager.addOutput( ...
    info, ...
    'outFile', ...
    'Unknown', ...
    'file', ...
    ['Generated registration QC artifacts and updated DataParams.mat ' ...
     'saved in SaveFolder.'], ...
    {'DataParams.mat','registrationQC_<timestamp>.fig','registrationQC_<timestamp>.png'}, ...
    1, ...
    'isData', false, ...
    'saveFileName', '');

end

end

% =========================================================================
% Local helpers
% =========================================================================

function targetFile = iResolveTargetFile(saveFolder, useFile)
%IRESOLVETARGETFILE Resolve the .dat file used to estimate registration.

if strcmpi(useFile, 'auto')
    fileList = dir(fullfile(saveFolder, '*.dat'));
    assert(~isempty(fileList), ...
        'Umitoolbox:createRegistrationTform:MissingDatFiles', ...
        'No .dat files were found in "%s".', saveFolder);

    [~, idx] = sort({fileList.name});
    fileList = fileList(idx);
    targetFile = fullfile(saveFolder, fileList(1).name);
    return
end

useFile = char(string(useFile));
if ~endsWith(useFile, '.dat', 'IgnoreCase', true)
    useFile = [useFile '.dat'];
end

targetFile = fullfile(saveFolder, useFile);
assert(isfile(targetFile), ...
    'Umitoolbox:createRegistrationTform:UseFileNotFound', ...
    'Requested file "%s" was not found in "%s".', useFile, saveFolder);

end

function img2D = iBuildMovingImageFromDat(datFile, refStatistic)
%IBUILDMOVINGIMAGEFROMDAT Build one 2D moving image from a YXT .dat file.

data = loadData(datFile);

assert(isnumeric(data) && ndims(data) == 3, ...
    'Umitoolbox:createRegistrationTform:InvalidDatDims', ...
    'The moving file must resolve to a 3D YXT array.');

switch refStatistic
    case 'mean'
        img2D = mean(data, 3, 'omitnan');
    case 'median'
        img2D = median(data, 3, 'omitnan');
    case 'first'
        img2D = data(:,:,1);
    otherwise
        error('Umitoolbox:createRegistrationTform:InvalidRefStatistic', ...
            'Unsupported RefStatistic "%s".', refStatistic);
end

img2D = single(img2D);

end

function MI = iMutualInformation(counts)
%IMUTUALINFORMATION Compute mutual information from a joint histogram.

counts = double(counts);
total = sum(counts(:));

if total <= 0
    MI = 0;
    return
end

pxy = counts / total;
px = sum(pxy, 2);
py = sum(pxy, 1);

mask = pxy > 0;
[iIdx, jIdx] = find(mask);

vals = zeros(numel(iIdx), 1);
for k = 1:numel(iIdx)
    vals(k) = pxy(iIdx(k), jIdx(k)) * log2( ...
        pxy(iIdx(k), jIdx(k)) / (px(iIdx(k)) * py(jIdx(k))));
end

MI = sum(vals);

end

function [translationXY, rotationDeg, scaleXY, determinantVal] = iExtractTformMetrics(tform)
%IEXTRACTTFORMMETRICS Extract QC metrics from a 2D transform object.

translationXY = [NaN NaN];
rotationDeg = NaN;
scaleXY = [NaN NaN];
determinantVal = NaN;

if isa(tform, 'affine2d') || isa(tform, 'projective2d') || ...
        isa(tform, 'images.geotrans.GeometricTransformation2D')
    T = tform.T;
    A = T(1:2,1:2);
    translationXY = T(3,1:2);
elseif isa(tform, 'rigidtform2d') || isa(tform, 'simtform2d') || isa(tform, 'affinetform2d')
    A = tform.A(1:2,1:2);
    translationXY = tform.A(3,1:2);
else
    return
end

sx = norm(A(:,1));
sy = norm(A(:,2));
scaleXY = [sx sy];
determinantVal = det(A);
rotationDeg = atan2d(A(2,1), A(1,1));

end

function [qcStatus, qcWarning] = iBuildQCStatus(MIDelta, translationXY, rotationDeg, scaleXY, imageSizeYX)
%IBUILDQCSTATUS Build simple advisory QC status and warning text.

qcStatus = 'estimated_not_reviewed';
warnings = {};

maxShift = max(abs(translationXY));
fovDiag = hypot(imageSizeYX(1), imageSizeYX(2));

if ~isfinite(MIDelta) || MIDelta <= 0
    warnings{end+1} = 'Mutual information did not improve after registration.'; %#ok<AGROW>
end

if isfinite(maxShift) && maxShift > 0.20 * fovDiag
    warnings{end+1} = 'Estimated translation is large relative to the field of view.'; %#ok<AGROW>
end

if isfinite(rotationDeg) && abs(rotationDeg) > 15
    warnings{end+1} = 'Estimated rotation is large.'; %#ok<AGROW>
end

if all(isfinite(scaleXY)) && any(abs(scaleXY - 1) > 0.10)
    warnings{end+1} = 'Estimated scale deviates substantially from 1.'; %#ok<AGROW>
end

if isempty(warnings)
    qcWarning = '';
else
    qcWarning = strjoin(warnings, ' ');
    qcStatus = 'review_recommended';
end

end

