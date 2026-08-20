function outFile = createRegistrationTform(SaveFolder, varargin)
%CREATEREGISTRATIONTFORM Estimate and store folder registration transform.
%
%   outFile = createRegistrationTform(SaveFolder)
%   outFile = createRegistrationTform(SaveFolder, 'Name', Value, ...)
%   info    = createRegistrationTform('pipelineInfo')
%
%   This function resolves SaveFolder's UMIT project binding and estimates a
%   2D geometric transform that aligns one image from the current folder to
%   the subject's active project-managed Image Reference. The exact fixed
%   image and managed-resource UUID/checksum provenance used for the estimate
%   are stored in DataParams.registration with the transform and QC metadata.
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
%                      'auto'  : .dat filename recorded in the active
%                                ImageReference.sourceFile provenance
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
%       outFile    - SaveFolder-relative file manifest of generated/updated
%                    files: {DataParams.mat, QC .fig, QC .png}
%
%   Notes:
%       - SaveFolder must have a valid UMITProjectBinding.umitlink whose
%         subject has an active, available managed Image Reference.
%       - Folder-local reference images are not used as a fallback.
%       - Managed Image Reference resources are read but never modified.
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
% Resolve and validate the active managed Image Reference
% -------------------------------------------------------------------------
dataParamsFile = fullfile(SaveFolder, 'DataParams.mat');
assert(isfile(dataParamsFile), ...
    'Umitoolbox:createRegistrationTform:MissingDataParams', ...
    'DataParams.mat was not found in "%s".', SaveFolder);

[~, managedReference, ImageReference] = ...
    iResolveManagedImageReference(SaveFolder);
refFr = single(ImageReference.image);

% Load only after managed-reference validation. loadDataParams normalizes
% older DataParams files in memory without writing them.
DataParams = loadDataParams(SaveFolder);
if ~isempty(DataParams.view.imageSizeYX) && ...
        ~isequal(size(refFr), double(DataParams.view.imageSizeYX(:).'))
    error('Umitoolbox:createRegistrationTform:ReferenceSizeMismatch', ...
        ['The active managed Image Reference size [%s] does not match ' ...
         'DataParams.view.imageSizeYX [%s].'], ...
        num2str(size(refFr)), ...
        num2str(double(DataParams.view.imageSizeYX(:).')));
end

% -------------------------------------------------------------------------
% Resolve moving-image source file
% -------------------------------------------------------------------------
targetFile = iResolveTargetFile(SaveFolder, UseFile, ImageReference);

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
% Band-pass via difference-of-Gaussians: a narrow 0.5-px blur removes pixel
% noise, and a wide blur (5% of the largest frame dimension, empirically
% chosen) approximates and removes the slow illumination gradient.
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
% Sweep of imregconfig('multimodal') optimizer settings from coarse/fast
% (index 1) to fine/slow (index 4); the best-scoring result across the sweep
% is kept below. GrowthFactor/Epsilon/InitialRadius values are empirically
% tuned step sizes, not derived from a formula.
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
tformQC = getTformQCMetrics(tform);
translationXY = tformQC.translationXY;
rotationDeg = tformQC.rotationDeg;
scaleXY = tformQC.scaleXY;
determinantVal = tformQC.determinant;

[qcStatus, qcWarning] = iBuildQCStatus( ...
    MIDelta, translationXY, rotationDeg, scaleXY, size(refFrMask));

% -------------------------------------------------------------------------
% Save QC artifacts
% -------------------------------------------------------------------------
ts = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
qcFigName = ['registrationQC_' ts '.fig'];
qcPngName = ['registrationQC_' ts '.png'];
qcFigFile = fullfile(SaveFolder, qcFigName);
qcPngFile = fullfile(SaveFolder, qcPngName);

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
DataParams.registration.referenceDescription = iReferenceDescription( ...
    managedReference, ImageReference);
DataParams.registration.referenceFile = managedReference.absolutePath;
DataParams.registration.referenceImage = refFr;
DataParams.registration.imageReferenceUUID = managedReference.uuid;
DataParams.registration.imageReferenceChecksum = managedReference.checksum;
DataParams.registration.createdOn = char( ...
    datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
DataParams.registration.source = targetFile;
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
DataParams.registration.qcFigureFile = qcFigName;
DataParams.registration.qcPreviewImageFile = qcPngName;
DataParams.registration.appliedOn = '';
DataParams.registration.appliedBy = '';
DataParams.registration.confirmationMode = '';
DataParams.registration.notes = '';
DataParams.registration.resourceUUID = '';

validateDataParams(DataParams);
saveDataParams(SaveFolder, DataParams);

outFile = {'DataParams.mat', qcFigName, qcPngName};

function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for this function.

info = PipelineManager.createPipelineInfo( ...
    mfilename, ...
    ['Estimate and store a folder registration transform using the bound ' ...
     'subject''s active project-managed Image Reference.']);

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    ['Project-bound folder containing DataParams.mat, AcqInfos.mat, and ' ...
     'source .dat files. Its subject must have an active managed Image ' ...
     'Reference.'], ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'UseFile', ...
    'parameter', ...
    ['Specific .dat file used to build the moving image. With ''auto'', ' ...
     'use the same .dat filename recorded in the active ' ...
     'ImageReference.sourceFile; do not fall back to another .dat file.'], ...
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
    'UnknownDataType', ...
    'file', ...
    ['Updated DataParams.mat saved in SaveFolder. The run also writes ' ...
     'timestamped registrationQC_*.fig and registrationQC_*.png artifacts ' ...
     'whose exact names are recorded in DataParams.registration.'], ...
    {'DataParams.mat'}, ...
    1, ...
    'isData', false, ...
    'saveFileName', '');

end

end

% =========================================================================
% Local helpers
% =========================================================================

function [bindingContext, resource, ImageReference] = ...
        iResolveManagedImageReference(saveFolder)
%IRESOLVEMANAGEDIMAGEREFERENCE Resolve and validate the fixed image source.

try
    [bindingContext, store] = ...
        UMITProjectStore.resolveProjectBinding(saveFolder);
catch ME
    error('Umitoolbox:createRegistrationTform:BindingResolutionFailed', ...
        ['Could not resolve the project binding for SaveFolder "%s". ' ...
         'No folder-local reference was used. %s'], saveFolder, ME.message);
end

try
    resource = store.getActiveImageReference(bindingContext.subjectID);
catch ME
    error( ...
        'Umitoolbox:createRegistrationTform:ActiveImageReferenceResolutionFailed', ...
        ['Could not resolve the active managed Image Reference for subject ' ...
         'UUID %s. %s'], bindingContext.subjectUUID, ME.message);
end

if isempty(resource)
    error('Umitoolbox:createRegistrationTform:NoActiveImageReference', ...
        ['Subject "%s" (%s) has no active managed Image Reference. ' ...
         'Select one in Image Reference Manager before estimating registration.'], ...
        bindingContext.subjectID, bindingContext.subjectUUID);
end

if ~strcmp(resource.type, 'imageReference') || ...
        ~strcmp(resource.status, 'active') || ...
        ~strcmpi(resource.ownerUUID, bindingContext.subjectUUID)
    error('Umitoolbox:createRegistrationTform:InvalidActiveImageReference', ...
        ['The resolved resource is not an active Image Reference owned by ' ...
         'the SaveFolder subject.']);
end

if ~resource.fileExists || ~isfile(resource.absolutePath)
    error('Umitoolbox:createRegistrationTform:MissingImageReferenceFile', ...
        'The active managed Image Reference file is missing: %s', ...
        resource.absolutePath);
end

try
    actualChecksum = computeFileChecksum(resource.absolutePath);
catch ME
    error('Umitoolbox:createRegistrationTform:InvalidImageReferenceFile', ...
        'Could not checksum the active managed Image Reference file: %s', ...
        ME.message);
end
if ~strcmp(actualChecksum, resource.checksum)
    error('Umitoolbox:createRegistrationTform:ImageReferenceChecksumMismatch', ...
        ['The active managed Image Reference file does not match its ' ...
         'registered SHA-256 checksum: %s'], resource.absolutePath);
end

try
    loaded = load(resource.absolutePath, 'ImageReference', '-mat');
catch ME
    error('Umitoolbox:createRegistrationTform:InvalidImageReferenceFile', ...
        'Could not load the active managed Image Reference file: %s', ...
        ME.message);
end
if ~isfield(loaded, 'ImageReference')
    error('Umitoolbox:createRegistrationTform:InvalidImageReferencePayload', ...
        ['The active managed Image Reference file does not contain the ' ...
         'required ImageReference variable: %s'], resource.absolutePath);
end

ImageReference = loaded.ImageReference;
try
    validateImageReferenceStruct(ImageReference);
catch ME
    error('Umitoolbox:createRegistrationTform:InvalidImageReferencePayload', ...
        'The active managed Image Reference payload is invalid: %s', ...
        ME.message);
end

if ~strcmpi(ImageReference.projectUUID, bindingContext.projectUUID) || ...
        ~strcmpi(ImageReference.subjectUUID, bindingContext.subjectUUID)
    error('Umitoolbox:createRegistrationTform:ImageReferenceIdentityMismatch', ...
        ['The active Image Reference payload project/subject UUIDs do not ' ...
         'match the resolved SaveFolder binding.']);
end

end


function description = iReferenceDescription(resource, ImageReference)
%IREFERENCEDESCRIPTION Choose the best managed-reference description.

candidates = { ...
    resource.description, ...
    ImageReference.description, ...
    resource.displayName, ...
    ImageReference.name};

description = '';
for k = 1:numel(candidates)
    candidate = strtrim(char(string(candidates{k})));
    if ~isempty(candidate)
        description = candidate;
        return
    end
end

end

function targetFile = iResolveTargetFile(saveFolder, useFile, ImageReference)
%IRESOLVETARGETFILE Resolve the .dat file used to estimate registration.

if strcmpi(useFile, 'auto')
    sourceFile = strtrim(char(string(ImageReference.sourceFile)));
    if isempty(sourceFile)
        error( ...
            'Umitoolbox:createRegistrationTform:AutoSourceFileUnavailable', ...
            ['UseFile=''auto'' requires the active Image Reference to record ' ...
             'the source channel filename in ImageReference.sourceFile. ' ...
             'Recreate the reference with source provenance or pass UseFile explicitly.']);
    end

    [~, sourceName, sourceExt] = fileparts(sourceFile);
    sourceFileName = [sourceName sourceExt];
    if isempty(sourceFileName) || ~strcmpi(sourceExt, '.dat')
        error( ...
            'Umitoolbox:createRegistrationTform:AutoSourceFileNotDat', ...
            ['UseFile=''auto'' requires the active Image Reference source ' ...
             'to be a .dat channel, but its sourceFile is "%s". Pass a ' ...
             'specific .dat filename with UseFile.'], sourceFile);
    end

    targetFile = fullfile(saveFolder, sourceFileName);
    if ~isfile(targetFile)
        error( ...
            'Umitoolbox:createRegistrationTform:AutoSourceFileNotFound', ...
            ['UseFile=''auto'' selected "%s" from the active Image ' ...
             'Reference provenance, but that file was not found in "%s".'], ...
            sourceFileName, saveFolder);
    end

    fprintf(['createRegistrationTform: UseFile=''auto'' selected "%s" ' ...
        'from the active Image Reference sourceFile provenance.\n'], ...
        sourceFileName);
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

function [qcStatus, qcWarning] = iBuildQCStatus(MIDelta, translationXY, rotationDeg, scaleXY, imageSizeYX)
%IBUILDQCSTATUS Build simple advisory QC status and warning text.

qcStatus = 'estimated_not_reviewed';
warnings = {};

maxShift = max(abs(translationXY));
fovDiag = hypot(imageSizeYX(1), imageSizeYX(2));

if ~isfinite(MIDelta) || MIDelta <= 0
    warnings{end+1} = 'Mutual information did not improve after registration.';
end

% QC thresholds below are empirically chosen review triggers, not derived
% limits: 20% of the frame diagonal for translation, 15 degrees for
% rotation, and 10% deviation from unity for scale.
if isfinite(maxShift) && maxShift > 0.20 * fovDiag
    warnings{end+1} = 'Estimated translation is large relative to the field of view.';
end

if isfinite(rotationDeg) && abs(rotationDeg) > 15
    warnings{end+1} = 'Estimated rotation is large.';
end

if all(isfinite(scaleXY)) && any(abs(scaleXY - 1) > 0.10)
    warnings{end+1} = 'Estimated scale deviates substantially from 1.';
end

if isempty(warnings)
    qcWarning = '';
else
    qcWarning = strjoin(warnings, ' ');
    qcStatus = 'review_recommended';
end

end

