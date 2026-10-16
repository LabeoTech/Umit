function outFile = createRegistrationTform(SaveFolder, varargin)
%CREATEREGISTRATIONTFORM Estimate and store folder registration transform.
%
%   outFile = createRegistrationTform(SaveFolder)
%   outFile = createRegistrationTform(SaveFolder, 'Name', Value, ...)
%   info    = createRegistrationTform('pipelineInfo')
%
%   This function resolves SaveFolder's UMIT project binding and estimates a
%   2D geometric transform that aligns one image from the current folder to
%   the subject's active project-managed Image Reference, resolved by the
%   shared resolveImageReferenceSource (source kind 'project'). A copy of
%   the reference image, its content hash and its source descriptor are
%   stored in DataParams.registration.reference; the transform, the QC
%   metadata and the hash it was estimated against in
%   DataParams.registration.estimate (version-2 record).
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
%       ReferenceSource - Where the fixed reference comes from (brief
%                    section 1.4): a reference source struct for
%                    resolveImageReferenceSource (kinds 'project', 'umitDat',
%                    'file') or a file path (a .dat is a 'umitDat' source,
%                    any other file a 'file' source: an exported reference
%                    or a TIF/MAT/UMT image). Resolution order:
%                      1) this parameter when given;
%                      2) else, for a project-bound folder, the bound
%                         subject's active Image Reference;
%                      3) else an error (Umitoolbox:ImageRegistration:
%                         noReferenceSource) explaining how to supply one.
%                    Projects are optional. Default: '' (steps 2-3)
%
%       UseFile    - File used to build the moving image. Use:
%                      'auto'  : the .dat named in the active
%                                ImageReference.sourceFile provenance when
%                                the folder has it; otherwise the folder's
%                                only Y-X-T .dat file (several candidates
%                                are an error: pass one explicitly)
%                      <name>  : specific .dat file name
%                    Default: 'auto'
%
%       ShowFigure - Logical scalar. If true, show the QC figure after
%                    creation. If false, create and save the figure
%                    invisibly.
%                    Default: true
%
%       RefStatistic - Statistic used to collapse the Y-X-T file into a 2D
%                    moving image (computed in X slabs, Low-RAM).
%                    Allowed: 'mean', 'median', 'first'
%                    Default: 'mean'
%
%   Output:
%       outFile    - SaveFolder-relative file manifest of generated/updated
%                    files: {DataParams.mat, QC .fig, QC .png}
%
%   Notes:
%       - SaveFolder must be bound (UMITProjectLedger.json) to a session whose
%         subject has an active, available managed Image Reference.
%       - Folder-local reference images are not used as a fallback.
%       - An aligned folder is refused
%         (Umitoolbox:ImageRegistration:alreadyAligned): only a fresh import
%         reopens alignment. The apply record is never modified here.
%       - The stored transform maps the folder's own image coordinates into
%         reference coordinates (the moving image is resized to the
%         reference for the estimate and the resize is composed back in),
%         which is the convention of applyImageAlignmentToFolder.
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
addParameter(p, 'ReferenceSource', '', ...
    @(x) isstruct(x) || ischar(x) || (isstring(x) && isscalar(x)));

parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
UseFile = char(string(p.Results.UseFile));
ShowFigure = p.Results.ShowFigure;
RefStatistic = lower(char(string(p.Results.RefStatistic)));
ReferenceSource = p.Results.ReferenceSource;

clear p

% -------------------------------------------------------------------------
% Resolve and validate the active managed Image Reference
% -------------------------------------------------------------------------
dataParamsFile = fullfile(SaveFolder, 'DataParams.mat');
assert(isfile(dataParamsFile), ...
    'Umitoolbox:createRegistrationTform:MissingDataParams', ...
    'DataParams.mat was not found in "%s".', SaveFolder);

% loadDataParams normalizes older DataParams files in memory without
% writing them. An aligned folder is committed to its alignment: estimating
% again is refused before anything is resolved or written.
DataParams = loadDataParams(SaveFolder);
operation = 'Estimating the image alignment';
if DataParams.registration.apply.isRegistered
    % Name the requested reference in the refusal when it can be resolved.
    requested = [];
    try
        requested = iResolveReference(SaveFolder, ReferenceSource);
    catch
    end
    assertFolderCanBeAligned(SaveFolder, DataParams, operation, requested);
end

% The reference, through the shared resolver (resolveImageReferenceSource):
% the explicit ReferenceSource, else the bound subject's active reference.
reference = iResolveReference(SaveFolder, ReferenceSource);
refFr = reference.image;

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
targetFile = iResolveTargetFile(SaveFolder, UseFile, reference.sourceFile);

% -------------------------------------------------------------------------
% Load moving image from file
% -------------------------------------------------------------------------
targetFr = iBuildMovingImageFromDat(targetFile, RefStatistic);

assert(~isempty(targetFr) && isnumeric(targetFr) && ismatrix(targetFr), ...
    'Umitoolbox:createRegistrationTform:InvalidMovingImage', ...
    'Failed to build a valid 2D moving image from "%s".', targetFile);

targetFr = single(targetFr);
movingSizeYX = size(targetFr);

% Resize moving image if needed so registration steps are well-defined. The
% stored transform is composed with this resize afterwards, so it always
% maps the folder's own image coordinates into reference coordinates (the
% convention of the shared apply backend applyImageAlignmentToFolder).
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

% QC metrics describe the estimated (resized-frame) transform; the stored
% transform maps native moving coordinates.
qcTform = tform;
tform = iComposeMovingResize(tform, movingSizeYX, size(refFr));

% -------------------------------------------------------------------------
% Derive QC metrics
% -------------------------------------------------------------------------
tformQC = getTformQCMetrics(qcTform);
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
% Update the estimate fields of DataParams.registration. The apply record
% (isRegistered, isReviewed, applied*, confirmationMode) is never touched
% here: the guard above guarantees the folder is not aligned.
% -------------------------------------------------------------------------
[~, referenceTemplate, estimateTemplate] = newRegistrationRecord();
createdOn = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

referenceRecord = referenceTemplate;
for fieldName = {'image', 'contentHash', 'name', 'description', 'axisConvention', ...
        'imageSizeYX', 'source'}
    referenceRecord.(fieldName{1}) = reference.(fieldName{1});
end
referenceRecord.selectedOn = createdOn;

estimate = estimateTemplate;
estimate.tform = tform;
estimate.transformType = 'similarity';
estimate.method = 'imregcorr initialization + multimodal imregtform optimization';
estimate.qcStatus = qcStatus;
estimate.qcWarning = qcWarning;
estimate.qcMetrics = struct( ...
    'MIBefore', MIBefore, ...
    'MIAfter', MIAfter, ...
    'MIDelta', MIDelta, ...
    'translationXY_px', translationXY, ...
    'rotationDeg', rotationDeg, ...
    'scaleXY', scaleXY, ...
    'determinant', determinantVal);
estimate.qcFigureFile = qcFigName;
estimate.qcPreviewImageFile = qcPngName;
estimate.createdOn = createdOn;
estimate.movingSource = targetFile;
estimate.referenceHash = reference.contentHash;

DataParams.registration.reference = referenceRecord;
DataParams.registration.estimate = estimate;

validateDataParams(DataParams);
saveDataParams(SaveFolder, DataParams);

outFile = {'DataParams.mat', qcFigName, qcPngName};

function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for this function.

info = PipelineManager.createPipelineInfo( ...
    mfilename, ...
    ['Estimate and store a folder registration transform against a reference: ' ...
     'the ReferenceSource given, else the bound subject''s active project-managed ' ...
     'Image Reference.']);

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    ['Folder containing DataParams.mat, AcqInfos.mat, and source .dat files. ' ...
     'Without ReferenceSource it must be project-bound and its subject must ' ...
     'have an active managed Image Reference.'], ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'UseFile', ...
    'parameter', ...
    ['Specific .dat file used to build the moving image. With ''auto'', ' ...
     'use the .dat named in the active ImageReference.sourceFile when the ' ...
     'folder has it, otherwise the folder''s only Y-X-T .dat file; several ' ...
     'candidates are an error.'], ...
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
    'ReferenceSource', ...
    'parameter', ...
    ['Reference file: a .dat of a UMIT save folder (mean of up to 100 evenly ' ...
     'spaced frames), an exported image reference (.mat) or a TIF/MAT/UMT image. ' ...
     'Empty: the bound subject''s active Image Reference. A .dat source from an ' ...
     'already aligned folder is allowed with a warning.'], ...
    'kind', 'parameter', ...
    'default', '', ...
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

function reference = iResolveReference(saveFolder, referenceSource)
%IRESOLVEREFERENCE Resolution order of the reference (brief section 3, Phase 4).

if ~(ischar(referenceSource) || isstring(referenceSource)) || strlength(string(referenceSource)) > 0
    reference = resolveImageReferenceSource( ...
        normalizeReferenceSourceInput(referenceSource), 'Interactive', false);
    return
end
switch UMITProjectStore.getProjectBindingState(saveFolder)
    case 'unbound'
        error('Umitoolbox:ImageRegistration:noReferenceSource', ...
            ['Estimating the image alignment was refused for folder "%s".\nWhy: no reference ' ...
             'source was given and the folder is not bound to a project.\nHow to fix:\n' ...
             '  (1) Pass ''ReferenceSource'': a .dat of a UMIT save folder, an exported ' ...
             'image reference (.mat), or a TIF/MAT/UMT image.\n' ...
             '  (2) Or bind the folder to a project session whose subject has an active ' ...
             'Image Reference.'], saveFolder);
    otherwise
        reference = iResolveProjectReference(saveFolder);
end
end

function reference = iResolveProjectReference(saveFolder)
%IRESOLVEPROJECTREFERENCE The bound subject's active library reference.

try
    [bindingContext, store] = ...
        UMITProjectStore.resolveProjectBinding(saveFolder);
catch ME
    error('Umitoolbox:createRegistrationTform:BindingResolutionFailed', ...
        ['Could not resolve the project binding for SaveFolder "%s". ' ...
         'No folder-local reference was used. %s'], saveFolder, ME.message);
end

reference = resolveImageReferenceSource(struct( ...
    'kind', 'project', ...
    'store', store, ...
    'subjectID', bindingContext.subjectID), 'Interactive', false);
end

function tform = iComposeMovingResize(tform, movingSizeYX, refSizeYX)
%ICOMPOSEMOVINGRESIZE Express TFORM in the native moving-image coordinates.
%
% The estimate runs on the moving image resized to the reference size.
% imresize maps intrinsic coordinate u to s*(u - 0.5) + 0.5 per axis, so the
% stored transform is TFORM applied after that scaling. It must stay a
% similarity transform (uniform scale): an aspect-ratio change cannot be
% applied by the shared backend.

movingSizeYX = double(movingSizeYX(1:2));
refSizeYX = double(refSizeYX(1:2));
if isequal(movingSizeYX, refSizeYX)
    return
end

sx = refSizeYX(2) / movingSizeYX(2);
sy = refSizeYX(1) / movingSizeYX(1);
if abs(sx - sy) > 1e-9 * max(sx, sy)
    error('Umitoolbox:createRegistrationTform:AnisotropicResize', ...
        ['The folder image [%s] and the reference [%s] have different aspect ' ...
         'ratios, so no similarity transform can map one onto the other.'], ...
        num2str(movingSizeYX), num2str(refSizeYX));
end

if isprop(tform, 'A')
    A = double(tform.A);
else
    A = double(tform.T).';
end
S = [sx 0 0.5 - 0.5 * sx; 0 sy 0.5 - 0.5 * sy; 0 0 1];
tform = affinetform2d(A * S);
end

function targetFile = iResolveTargetFile(saveFolder, useFile, sourceFile)
%IRESOLVETARGETFILE Resolve the .dat file used to estimate registration.

if strcmpi(useFile, 'auto')
    % 1) The .dat channel the reference was built from, when the folder has
    %    a file of that name.
    sourceFile = strtrim(char(string(sourceFile)));
    [~, sourceName, sourceExt] = fileparts(sourceFile);
    sourceFileName = [sourceName sourceExt];
    if strcmpi(sourceExt, '.dat') && isfile(fullfile(saveFolder, sourceFileName))
        targetFile = fullfile(saveFolder, sourceFileName);
        fprintf(['createRegistrationTform: UseFile=''auto'' selected "%s" ' ...
            'from the Image Reference source provenance.\n'], sourceFileName);
        return
    end

    % 2) Otherwise (reference built from a TIF/MAT/UMT file, or from a .dat
    %    of another name): the folder's only Y-X-T .dat. Several candidates
    %    are never guessed between.
    candidates = iYXTDatFiles(saveFolder);
    if isscalar(candidates)
        targetFile = fullfile(saveFolder, candidates{1});
        fprintf(['createRegistrationTform: UseFile=''auto'' selected "%s", the ' ...
            'only Y-X-T .dat file in the folder.\n'], candidates{1});
        return
    end
    if isempty(sourceFile)
        sourceText = 'records no source file';
    else
        sourceText = sprintf('was built from "%s"', sourceFile);
    end
    if isempty(candidates)
        error('Umitoolbox:createRegistrationTform:AutoSourceFileNotFound', ...
            ['UseFile=''auto'' found no Y-X-T .dat file in "%s" (the Image ' ...
             'Reference %s).'], saveFolder, sourceText);
    end
    error('Umitoolbox:createRegistrationTform:AutoSourceFileAmbiguous', ...
        ['UseFile=''auto'' cannot choose the moving image: the Image Reference %s, ' ...
         'and "%s" holds several Y-X-T .dat files (%s). Pass one with UseFile.'], ...
        sourceText, saveFolder, strjoin(candidates, ', '));
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

function names = iYXTDatFiles(saveFolder)
%IYXTDATFILES Names of the folder's .dat files with a Y-X-T layout.

listing = dir(fullfile(saveFolder, '*.dat'));
names = {};
for iFile = 1:numel(listing)
    try
        info = loadMetaData(fullfile(saveFolder, listing(iFile).name));
    catch
        continue
    end
    dimNames = cellstr(string(info.dimNames));
    if isequal(dimNames(:).', {'Y', 'X', 'T'})
        names{end+1} = listing(iFile).name; %#ok<AGROW>
    end
end
names = sort(names);
end

function img2D = iBuildMovingImageFromDat(datFile, refStatistic)
%IBUILDMOVINGIMAGEFROMDAT Build one 2D moving image from a Y-X-T .dat file.
%
% The file is read in X slabs sized from the available RAM, so the whole
% recording is never loaded. Only Y-X-T files are accepted.

info = loadMetaData(datFile);
assertDatLayout(info, {{'Y','X','T'}}, 'createRegistrationTform');
nY = datAxisSize(info, 'Y');
nX = datAxisSize(info, 'X');
nT = datAxisSize(info, 'T');

if ~ismember(refStatistic, {'mean', 'median', 'first'})
    error('Umitoolbox:createRegistrationTform:InvalidRefStatistic', ...
        'Unsupported RefStatistic "%s".', refStatistic);
end

slabIn = spatialSlabIO('open', datFile, 'Info', info);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

img2D = zeros(nY, nX, 'single');

if strcmp(refStatistic, 'first')
    img2D(:, :) = single(spatialSlabIO('read', slabIn, 1:nX, 1));
    return
end

% A slab and the reduction temporaries.
nChunks = calculateMaxChunkSize(double(nY) * nX * nT * 4, 3, 0.2);
chunkX = max(1, ceil(nX / nChunks));
for xStart = 1:chunkX:nX
    xIdx = xStart:min(xStart + chunkX - 1, nX);
    slab = single(spatialSlabIO('read', slabIn, xIdx));
    switch refStatistic
        case 'mean'
            img2D(:, xIdx) = mean(slab, 3, 'omitnan');
        case 'median'
            img2D(:, xIdx) = median(slab, 3, 'omitnan');
    end
end

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

