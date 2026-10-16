function report = applyImageAlignmentToFolder(targetFolder, finalTform, opts)
%APPLYIMAGEALIGNMENTTOFOLDER Apply an image alignment to a whole save folder.
%
%   report = applyImageAlignmentToFolder(targetFolder, finalTform, ...
%       'imageReference', ImageReference, ...
%       'referenceFile', referenceFileRel, ...
%       'referenceImage', referenceImage)
%
%   The single (canonical) apply path of UMIT image alignment. The GUI
%   (ImageAlignmentTool) calls it directly; the pipeline function
%   applyRegistrationTformOnFolder is a wrapper around it. Both therefore
%   produce the same folder state and share the same guards.
%
%   Processing steps:
%       1) Guard: an already-aligned folder is refused
%          (assertFolderCanBeAligned; Umitoolbox:ImageRegistration:
%          alreadyAligned). Aligning again would transform the data twice;
%          only a fresh import of the raw data reopens alignment.
%       2) Validate the folder, DataParams, every .dat, .umt and .roi file,
%          and the new DataParams, before any file is written.
%       3) Stage: write every transformed file (.dat, image .umt, current
%          .roi) and the alignment report as temporary siblings
%          (<name>.aligning.tmp). The originals are untouched.
%       4) Commit: move each original aside (<name>.aligning.bak), install
%          its transformed copy, then save DataParams.mat through
%          validateDataParams + atomic saveDataParams. Any failure restores
%          every original, so the folder is either fully aligned or
%          unchanged. The .bak files are deleted after a successful commit.
%       5) No backup zip is written: alignment is irreversible by design,
%          and a fresh import is the only way back.
%       AcqInfos.mat is not rewritten: every aligned .dat is headered and
%       describes its new frame size itself.
%
%   Transform convention:
%       finalTform maps the folder's current image coordinates (the .dat
%       frame) into reference coordinates; the output frame is the
%       reference image size.
%
%   Name-Value options:
%       imageReference       - ImageReference struct used as fixed reference.
%       referenceFile        - Reference file path, recorded as the path of
%                              the reference source descriptor when
%                              referenceSource is not given.
%       referenceSource      - Source descriptor of the reference
%                              (referenceSourceDescriptor, e.g. from
%                              resolveImageReferenceSource). Default: built
%                              from imageReference/referenceFile ('project'
%                              when imageReference carries a project
%                              identity, else 'file').
%       referenceImage       - Fixed reference image. If empty, this is read
%                              from imageReference.image.
%       baseTform            - Automatic/base transform before manual tuning.
%       previewTform         - Preview transform used in the GUI.
%       manualAdjustment     - Struct with dx/dy/rotation/scale.
%       movingSourceFile     - File or context used to estimate the transform.
%       preFlipHorizontally  - Whether the GUI pre-flipped the moving source.
%       createdBy            - Tool/app name; recorded as appliedBy (and as
%                              method when keepEstimate is false).
%       keepEstimate         - true when finalTform IS the transform of
%                              DataParams.registration.estimate (pipeline
%                              apply of createRegistrationTform's result):
%                              the reference and estimate records are kept
%                              and only the apply record is written; the
%                              estimate must have been made against the
%                              recorded reference. false (default, GUI): the
%                              reference and estimate records are rebuilt
%                              from these options.
%       confirmationMode     - Recorded confirmation mode.
%                              Default: 'user-confirmed-folder-transform'.
%       isReviewed           - Recorded review flag. Default: true.
%       resourceUUID         - RETIRED and IGNORED (brief section 4 Q5):
%                              accepted only until the ImageAlignmentTool
%                              edit that stops passing it is applied.
%       transformDat         - Transform .dat files. Default: true.
%       transformUmt         - Transform eligible image .umt files. Default: true.
%       transformRoi         - Transform current-format .roi files. Default: true.
%
%   Output:
%       report - Compact operation report (also saved as
%                manualAlignmentReport_<datetime>.mat).
%
%   Notes:
%       .dat files must have Y and X as their first two axes (Y-X, Y-X-T,
%       Y-X-E, Y-X-T-E, ...) and share one frame size; every Y-X plane is
%       warped with the same transform (linear interpolation, zero fill).
%       NaN pixels of floating-point data are zeroed for the interpolation
%       and restored through the warped NaN mask. Each aligned file is
%       written with a header: the new frame size, and the input's class,
%       trailing axes, frame rate, exposure, and name. Image .umt entries
%       and DataParams.mask must have the same frame size as the .dat files;
%       DataParams.view.imageSizeYX, when set, must too.

arguments
    targetFolder (1,1) string
    finalTform
    opts.imageReference = []
    opts.referenceFile (1,1) string = ""
    opts.referenceImage = []
    opts.baseTform = []
    opts.previewTform = []
    opts.manualAdjustment = struct()
    opts.movingSourceFile (1,1) string = ""
    opts.preFlipHorizontally (1,1) logical = false
    opts.createdBy (1,1) string = "ImageAlignmentTool"
    opts.referenceSource = []
    opts.keepEstimate (1,1) logical = false
    opts.confirmationMode (1,1) string = "user-confirmed-folder-transform"
    opts.isReviewed (1,1) logical = true
    opts.resourceUUID (1,1) string = ""
    opts.transformDat (1,1) logical = true
    opts.transformUmt (1,1) logical = true
    opts.transformRoi (1,1) logical = true
end

targetFolder = char(targetFolder);

if ~isfolder(targetFolder)
    error('ImageAlignmentTool:MissingTargetFolder', ...
        'Target folder does not exist: %s', targetFolder);
end

finalTform = iNormalizeTransform(finalTform);

referenceImage = opts.referenceImage;
if isempty(referenceImage) && isstruct(opts.imageReference) && ...
        isfield(opts.imageReference, 'image')
    referenceImage = opts.imageReference.image;
end
referenceImage = iValidateReferenceImage(referenceImage);

newSizeYX = size(referenceImage);
referenceView = imref2d(newSizeYX);

createdOn = datetime('now');
dateTag = char(datetime(createdOn, 'Format', 'yyyyMMdd_HHmmss'));

report = iInitializeReport( ...
    targetFolder, ...
    finalTform, ...
    opts, ...
    newSizeYX, ...
    createdOn);

DataParams = iLoadDataParams(targetFolder);
operation = 'Applying the image alignment';
if opts.keepEstimate
    requested = DataParams.registration.reference;
else
    requested = struct('name', 'the reference given to the tool', ...
        'contentHash', computeImageContentHash(referenceImage), ...
        'source', iReferenceSource(opts));
    if isstruct(opts.imageReference) && isfield(opts.imageReference, 'name') && ...
            ~isempty(opts.imageReference.name)
        requested.name = char(string(opts.imageReference.name));
    end
end
assertFolderCanBeAligned(targetFolder, DataParams, operation, requested);
if opts.keepEstimate
    verifyRegistrationReference(targetFolder, DataParams.registration, operation);
    iAssertEstimateMatchesReference(targetFolder, DataParams.registration, operation);
end
AcqInfoStream = iLoadAcqInfos(targetFolder);

datFiles = dir(fullfile(targetFolder, '*.dat'));
datInfos = iDescribeDatFiles(datFiles, opts.transformDat);

oldSizeYX = iResolveOldSizeYX(AcqInfoStream, DataParams, datFiles, datInfos);
report.oldSizeYX = oldSizeYX;

scaleInfo = iGetSimilarityScaleInfo(finalTform);
if ~scaleInfo.isSimilarity
    error('ImageAlignmentTool:NonSimilarityTransform', ...
        ['Registration transform is not a valid similarity transform. ' ...
         'Uniform scale/no-shear validation failed.']);
end
report.scaleInfo = scaleInfo;

umtFiles = dir(fullfile(targetFolder, '*.umt'));
roiFiles = dir(fullfile(targetFolder, '*.roi'));

if isempty(datFiles) && isempty(umtFiles)
    error('ImageAlignmentTool:NoImageFilesFound', ...
        'No .dat or .umt image files were found in target folder.');
end

% Build the complete plan before any write.
iPreflightDataParams(DataParams, oldSizeYX);

datPlan = iPreflightDatFiles(datFiles, datInfos, oldSizeYX, opts.transformDat);
[umtPlan, umtWarnings] = iPreflightUmtFiles(umtFiles, oldSizeYX, opts.transformUmt);
[roiPlan, roiWarnings] = iPreflightRoiFiles(roiFiles, oldSizeYX, opts.transformRoi);

report.warnings = [report.warnings; umtWarnings(:); roiWarnings(:)];

% The new DataParams is built and validated before any file is written.
precomputedMask = iPrecomputeTransformedMask( ...
    DataParams, oldSizeYX, finalTform, referenceView);
updatedView = iBuildUpdatedView(DataParams, opts.imageReference, newSizeYX, scaleInfo, finalTform);
newDataParams = iUpdateDataParams(DataParams, finalTform, referenceImage, opts, ...
    precomputedMask, scaleInfo, updatedView);
validateDataParams(newDataParams);

staged = iNewStagedList();

try
    % Stage every transformed file next to its original.
    for iFile = 1:numel(datPlan)
        tmpPath = iStagePath(datPlan(iFile).Path);
        staged = iAddStaged(staged, datPlan(iFile).Path, tmpPath);
        iTransformDatFile(datPlan(iFile), finalTform, referenceView, tmpPath);
        report.datFilesTransformed(end+1,1) = string(datPlan(iFile).Name);
    end

    for iFile = 1:numel(umtPlan)
        if ~umtPlan(iFile).Transform
            report.umtFilesSkipped(end+1,1) = ...
                string(umtPlan(iFile).Name) + ": " + string(umtPlan(iFile).Reason);
            continue
        end
        tmpPath = iStagePath(umtPlan(iFile).Path);
        staged = iAddStaged(staged, umtPlan(iFile).Path, tmpPath);
        iTransformUmtFile(umtPlan(iFile), finalTform, referenceView, tmpPath);
        report.umtFilesTransformed(end+1,1) = string(umtPlan(iFile).Name);
    end

    for iFile = 1:numel(roiPlan)
        if ~roiPlan(iFile).Transform
            report.roiFilesSkipped(end+1,1) = ...
                string(roiPlan(iFile).Name) + ": " + string(roiPlan(iFile).Reason);
            continue
        end
        tmpPath = iStagePath(roiPlan(iFile).Path);
        staged = iAddStaged(staged, roiPlan(iFile).Path, tmpPath);
        iTransformRoiFile(roiPlan(iFile), finalTform, referenceView, updatedView, tmpPath);
        report.roiFilesTransformed(end+1,1) = string(roiPlan(iFile).Name);
    end

    report.metadataFilesUpdated(end+1,1) = "DataParams.mat";
    reportPath = fullfile(targetFolder, ['manualAlignmentReport_' dateTag '.mat']);
    report.reportPath = reportPath;
    report.finishedOn = datetime('now');
    tmpPath = iStagePath(reportPath);
    staged = iAddStaged(staged, reportPath, tmpPath);
    save(tmpPath, 'report', '-mat');
catch ME
    iDiscardStaged(staged);
    iWriteFailureReport(targetFolder, iFailureInfo(ME, report, dateTag, ...
        'staging (no original file was changed)'));
    rethrow(ME);
end

% Commit: install every staged file, then DataParams. All or nothing.
installed = iNewStagedList();
try
    for iFile = 1:numel(staged)
        installed = iInstallStaged(installed, staged(iFile));
    end
    saveDataParams(targetFolder, newDataParams);
catch ME
    restoreMessage = iRollBack(installed);
    iDiscardStaged(staged);
    iWriteFailureReport(targetFolder, iFailureInfo(ME, report, dateTag, ...
        sprintf('commit; originals restored: %s', restoreMessage)));
    rethrow(ME);
end
iDeleteBackups(installed);

end

% =========================================================================
% Report / metadata loading
% =========================================================================

function report = iInitializeReport(targetFolder, finalTform, opts, newSizeYX, createdOn)
%IINITIALIZEREPORT Create compact operation report.

report = struct();
report.createdOn = createdOn;
report.finishedOn = NaT;
report.toolName = char(opts.createdBy);
report.targetFolder = targetFolder;
report.referenceFile = char(opts.referenceFile);
report.referenceImageSizeYX = newSizeYX;
report.oldSizeYX = [];
report.newSizeYX = newSizeYX;
report.transformType = 'similarity';
report.finalTform = finalTform;
report.baseTform = opts.baseTform;
report.previewTform = opts.previewTform;
report.manualAdjustment = opts.manualAdjustment;
report.movingSourceFile = char(opts.movingSourceFile);
report.preFlipHorizontally = logical(opts.preFlipHorizontally);
report.createdBy = char(opts.createdBy);
report.reportPath = '';
report.datFilesTransformed = strings(0,1);
report.umtFilesTransformed = strings(0,1);
report.umtFilesSkipped = strings(0,1);
report.roiFilesTransformed = strings(0,1);
report.roiFilesSkipped = strings(0,1);
report.metadataFilesUpdated = strings(0,1);
report.warnings = strings(0,1);
report.scaleInfo = struct();
end

function DataParams = iLoadDataParams(targetFolder)
%ILOADDATAPARAMS Load required DataParams.mat (normalized and validated).

if ~isfile(fullfile(targetFolder, 'DataParams.mat'))
    error('ImageAlignmentTool:MissingDataParams', ...
        'DataParams.mat is required before applying registration.');
end

DataParams = loadDataParams(targetFolder);
end

function AcqInfoStream = iLoadAcqInfos(targetFolder)
%ILOADACQINFOS Load required AcqInfos.mat.

acqInfosPath = fullfile(targetFolder, 'AcqInfos.mat');

if ~isfile(acqInfosPath)
    error('ImageAlignmentTool:MissingAcqInfos', ...
        'AcqInfos.mat is required before applying registration.');
end

S = load(acqInfosPath);

if ~isfield(S, 'AcqInfoStream')
    error('ImageAlignmentTool:InvalidAcqInfosFile', ...
        'AcqInfos.mat does not contain variable AcqInfoStream.');
end

AcqInfoStream = S.AcqInfoStream;
end

% =========================================================================
% Plan / preflight
% =========================================================================

function iPreflightDataParams(DataParams, oldSizeYX)
%IPREFLIGHTDATAPARAMS Validate DataParams fields affected by alignment.
%
%   DataParams.view and DataParams.mask describe the current image frame,
%   which alignment replaces; a stale size would give a wrong view or mask.

if isfield(DataParams, 'view') && isstruct(DataParams.view) && ...
        isfield(DataParams.view, 'imageSizeYX') && ~isempty(DataParams.view.imageSizeYX) && ...
        ~isequal(double(DataParams.view.imageSizeYX(:).'), double(oldSizeYX))
    error('ImageAlignmentTool:ViewSizeMismatch', ...
        ['DataParams.view.imageSizeYX [%s] does not match the folder''s image ' ...
         'data [%s]. The folder''s DataParams are inconsistent; no file was changed.'], ...
        num2str(double(DataParams.view.imageSizeYX(:).')), num2str(oldSizeYX));
end

if isfield(DataParams, 'mask') && isstruct(DataParams.mask) && ...
        isfield(DataParams.mask, 'logical') && ~isempty(DataParams.mask.logical)

    maskSize = size(DataParams.mask.logical);

    if numel(maskSize) < 2 || ~isequal(maskSize(1:2), oldSizeYX)
        error('ImageAlignmentTool:MaskSizeMismatch', ...
            ['DataParams.mask.logical size does not match the current folder ' ...
             'image size. Operation aborted to avoid warping a stale mask.']);
    end
end
end

function transformedMask = iPrecomputeTransformedMask(DataParams, oldSizeYX, finalTform, referenceView)
%IPRECOMPUTETRANSFORMEDMASK Warp DataParams mask before destructive writes.

transformedMask = [];

if ~isfield(DataParams, 'mask') || ~isstruct(DataParams.mask) || ...
        ~isfield(DataParams.mask, 'logical') || isempty(DataParams.mask.logical)
    return
end

oldMask = logical(DataParams.mask.logical);

if ~isequal(size(oldMask), oldSizeYX)
    error('ImageAlignmentTool:MaskSizeMismatch', ...
        'DataParams.mask.logical size does not match current folder image size.');
end

transformedMask = imwarp(single(oldMask), finalTform, ...
    'OutputView', referenceView, ...
    'InterpolationMethod', 'nearest', ...
    'FillValues', single(0)) > 0.5;
end

function datInfos = iDescribeDatFiles(datFiles, doTransform)
%IDESCRIBEDATFILES loadMetaData of every .dat that will be transformed.

datInfos = cell(1, numel(datFiles));
if ~doTransform
    return
end
for iFile = 1:numel(datFiles)
    datPath = fullfile(datFiles(iFile).folder, datFiles(iFile).name);
    if ~isfile(datPath)
        error('ImageAlignmentTool:MissingDatFile', ...
            'DAT file disappeared during preflight: %s', datPath);
    end
    info = loadMetaData(datPath);
    dimNames = cellstr(string(info.dimNames));
    dimNames = dimNames(:).';
    if numel(dimNames) < 2 || ~isequal(dimNames(1:2), {'Y', 'X'})
        error('ImageAlignmentTool:UnsupportedLayout', ...
            ['DAT file "%s" has dimensions {%s}. Alignment needs image data ' ...
             'whose first two axes are Y and X.'], ...
            datFiles(iFile).name, strjoin(dimNames, ','));
    end
    datInfos{iFile} = info;
end
end

function datPlan = iPreflightDatFiles(datFiles, datInfos, oldSizeYX, doTransform)
%IPREFLIGHTDATFILES Validate every .dat layout and size before any change.
%
% Any .dat whose first two axes are Y and X is aligned plane by plane; the
% trailing axes (T, E, F, ...) only enumerate planes and are kept.

datPlan = struct( ...
    'Name', {}, ...
    'Path', {}, ...
    'OldSizeYX', {}, ...
    'NPlanes', {}, ...
    'TrailingSizes', {}, ...
    'Info', {});

if ~doTransform
    return
end

for iFile = 1:numel(datFiles)
    datPath = fullfile(datFiles(iFile).folder, datFiles(iFile).name);
    info = datInfos{iFile};

    dimNames = cellstr(string(info.dimNames));
    dimNames = dimNames(:).';
    if numel(dimNames) < 2 || ~isequal(dimNames(1:2), {'Y', 'X'})
        error('ImageAlignmentTool:UnsupportedLayout', ...
            ['DAT file "%s" has dimensions {%s}. Alignment needs image data ' ...
             'whose first two axes are Y and X.'], ...
            datFiles(iFile).name, strjoin(dimNames, ','));
    end

    fileSizeYX = [datAxisSize(info, 'Y'), datAxisSize(info, 'X')];
    if ~isequal(fileSizeYX, oldSizeYX)
        error('ImageAlignmentTool:InvalidDatFileSize', ...
            ['DAT file "%s" has frames of %d x %d, but the folder''s image ' ...
             'size is %d x %d.'], ...
            datFiles(iFile).name, fileSizeYX(1), fileSizeYX(2), oldSizeYX(1), oldSizeYX(2));
    end

    dimSizes = double(info.dimSizes(:).');
    trailingSizes = dimSizes(3:end);
    nPlanes = prod(trailingSizes);
    if nPlanes < 1
        error('ImageAlignmentTool:InvalidDatFrameCount', ...
            'DAT file "%s" has no frames.', datFiles(iFile).name);
    end

    planEntry = struct( ...
        'Name', datFiles(iFile).name, ...
        'Path', datPath, ...
        'OldSizeYX', oldSizeYX, ...
        'NPlanes', nPlanes, ...
        'TrailingSizes', trailingSizes, ...
        'Info', info);
    datPlan(end+1) = planEntry; %#ok<AGROW>
end
end

function [umtPlan, warnings] = iPreflightUmtFiles(umtFiles, oldSizeYX, doTransform)
%IPREFLIGHTUMTFILES Validate and classify UMT files.

umtPlan = struct( ...
    'Name', {}, ...
    'Path', {}, ...
    'Transform', {}, ...
    'Reason', {}, ...
    'UmtVariableName', {}, ...
    'EntryNames', {});

warnings = strings(0,1);

if ~doTransform
    return
end

for iFile = 1:numel(umtFiles)
    umtPath = fullfile(umtFiles(iFile).folder, umtFiles(iFile).name);

    try
        S = load(umtPath,'-mat');
    catch ME
        error('ImageAlignmentTool:CorruptUmtFile', ...
            'Could not load UMT file "%s": %s', umtFiles(iFile).name, ME.message);
    end

    [umtVarName, umt] = iFindSingleUmtVariable(S, umtFiles(iFile).name);

    if isempty(umtVarName)
        planEntry = struct( ...
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', false, ...
            'Reason', "no UMT structure found", ...
            'UmtVariableName', '', ...
            'EntryNames', {{}});
        umtPlan(end+1) = planEntry; %#ok<AGROW>
        continue
    end

    validateUMTStruct(umt, 'requireEventInfo', false);

    if ~isfield(umt, 'kind') || ~strcmpi(char(string(umt.kind)), 'image')
        planEntry = struct( ...
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', false, ...
            'Reason', "kind is not image", ...
            'UmtVariableName', umtVarName, ...
            'EntryNames', {{}});
        umtPlan(end+1) = planEntry; %#ok<AGROW>
        continue
    end

    schema = getUMTSchema(umt.version);
    entryNames = fieldnames(umt.data);
    transformEntries = {};
    skippedEntries = strings(0,1);

    for iEntry = 1:numel(entryNames)
        entryName = entryNames{iEntry};
        entry = umt.data.(entryName);
        dimNames = iNormalizeDimNamesForPlan(entry.dimNames);

        if iIsAllowedImageSpatialPattern(dimNames, schema)
            entrySizeYX = double([size(entry.value, 1), size(entry.value, 2)]);
            if ~isequal(entrySizeYX, double(oldSizeYX))
                error('ImageAlignmentTool:UmtSizeMismatch', ...
                    ['UMT file "%s" entry "%s" is %d x %d, but the folder''s image ' ...
                     'size is %d x %d. No file was changed.'], umtFiles(iFile).name, ...
                    entryName, entrySizeYX(1), entrySizeYX(2), oldSizeYX(1), oldSizeYX(2));
            end
            transformEntries{end+1} = entryName; %#ok<AGROW>
        else
            skippedEntries(end+1,1) = string(entryName) + " (" + strjoin(string(dimNames), "/") + ")"; %#ok<AGROW>
        end
    end

    if isempty(transformEntries)
        planEntry = struct( ...
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', false, ...
            'Reason', "no eligible image entries", ...
            'UmtVariableName', umtVarName, ...
            'EntryNames', {{}});
        umtPlan(end+1) = planEntry; %#ok<AGROW>
    else
        if ~isempty(skippedEntries)
            warningText = sprintf( ...
                'UMT file "%s": skipped unsupported image entries: %s.', ...
                umtFiles(iFile).name, strjoin(skippedEntries, ', '));
            warnings(end+1,1) = warningText; %#ok<AGROW>
        end

        planEntry = struct( ...
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', true, ...
            'Reason', "", ...
            'UmtVariableName', umtVarName, ...
            'EntryNames', {transformEntries});
        umtPlan(end+1) = planEntry; %#ok<AGROW>
    end
end
end

function [roiPlan, warnings] = iPreflightRoiFiles(roiFiles, oldSizeYX, doTransform)
%IPREFLIGHTROIFILES Validate and classify ROI files.

roiPlan = struct( ...
    'Name', {}, ...
    'Path', {}, ...
    'Transform', {}, ...
    'Reason', {});

warnings = strings(0,1);

if ~doTransform
    return
end

for iFile = 1:numel(roiFiles)
    roiPath = fullfile(roiFiles(iFile).folder, roiFiles(iFile).name);

    try
        S = load(roiPath, 'ROIFile','-mat');
    catch ME
        error('ImageAlignmentTool:CorruptRoiFile', ...
            'Could not load ROI file "%s": %s', roiFiles(iFile).name, ME.message);
    end

    if ~isfield(S, 'ROIFile')
        error('ImageAlignmentTool:MissingRoiFileVariable', ...
            'ROI file "%s" does not contain variable ROIFile.', roiFiles(iFile).name);
    end

    ROIFile = S.ROIFile;

    if ~iIsCurrentRoiFile(ROIFile)
        warningText = sprintf( ...
            'ROI file "%s" is not a current UMIT ROI file and will be skipped.', ...
            roiFiles(iFile).name);
        warnings(end+1,1) = warningText; %#ok<AGROW>

        planEntry = struct( ...
            'Name', roiFiles(iFile).name, ...
            'Path', roiPath, ...
            'Transform', false, ...
            'Reason', "not a current UMIT ROI file");
        roiPlan(end+1) = planEntry; %#ok<AGROW>
        continue
    end

    roiSizeYX = double(ROIFile.imageInfo.imageSizeYX(:).');
    if ~isequal(roiSizeYX, double(oldSizeYX))
        error('ImageAlignmentTool:RoiSizeMismatch', ...
            ['ROI file "%s" was drawn on a %d x %d image, but the folder''s image ' ...
             'size is %d x %d. No file was changed.'], roiFiles(iFile).name, ...
            roiSizeYX(1), roiSizeYX(2), oldSizeYX(1), oldSizeYX(2));
    end

    ROIs = ROIFile.ROIs;
    for iROI = 1:numel(ROIs)
        ROI = ROIs(iROI);

        if ~isfield(ROI, 'geometry') || ...
                ~isfield(ROI.geometry, 'verticesXY_px') || ...
                isempty(ROI.geometry.verticesXY_px)
            error('ImageAlignmentTool:UnsupportedRoiGeometry', ...
                'ROI file "%s" contains ROI %d without geometry.verticesXY_px.', ...
                roiFiles(iFile).name, iROI);
        end

        vertices = ROI.geometry.verticesXY_px;
        if ~isnumeric(vertices) || size(vertices,2) ~= 2 || size(vertices,1) < 3
            error('ImageAlignmentTool:InvalidRoiVertices', ...
                'ROI file "%s" contains ROI %d with invalid verticesXY_px.', ...
                roiFiles(iFile).name, iROI);
        end
    end

    planEntry = struct( ...
        'Name', roiFiles(iFile).name, ...
        'Path', roiPath, ...
        'Transform', true, ...
        'Reason', "");
    roiPlan(end+1) = planEntry; %#ok<AGROW>
end
end

% =========================================================================
% DAT transformation
% =========================================================================

function iTransformDatFile(plan, tform, referenceView, tmpPath)
%ITRANSFORMDATFILE Write the aligned copy of one .dat file to tmpPath.
%
% Plane by plane through spatialSlabIO. The aligned file is headered: the
% reference frame size, the input's trailing axes, class, frame rate,
% exposure, and name. Planes stay in the input's class through imwarp.

datPath = plan.Path;
oldSizeYX = plan.OldSizeYX;
newSizeYX = referenceView.ImageSize;

[~, channelName] = fileparts(datPath);
outHeader = datHeaderFromInfo(plan.Info, channelName, ...
    'dimSizes', [newSizeYX(1), newSizeYX(2), plan.TrailingSizes]);

slabIn = spatialSlabIO('open', datPath, 'Info', plan.Info);
cleanupIn = onCleanup(@() spatialSlabIO('close', slabIn));
slabOut = spatialSlabIO('create', tmpPath, outHeader);
cleanupOut = onCleanup(@() spatialSlabIO('close', slabOut));

for iPlane = 1:plan.NPlanes
    plane = reshape(spatialSlabIO('read', slabIn, 1:oldSizeYX(2), iPlane), oldSizeYX);
    outPlane = iWarpPlane(plane, tform, referenceView);
    spatialSlabIO('write', slabOut, 1:newSizeYX(2), outPlane, iPlane);
end

spatialSlabIO('finalize', slabOut);
clear cleanupIn cleanupOut
end

function out = iWarpPlane(plane, tform, referenceView)
%IWARPPLANE Warp one Y-X plane (linear, zero fill, class kept).
%
% NaN pixels of floating-point data (outside the field of view, padding of
% short trials) are zeroed for the interpolation and restored through the
% nearest-neighbour warp of the NaN mask, so they do not spread.

if islogical(plane)
    out = imwarp(single(plane), tform, ...
        'OutputView', referenceView, ...
        'InterpolationMethod', 'nearest', ...
        'FillValues', single(0)) > 0.5;
    return
end

nanMask = [];
if isfloat(plane)
    nanMask = isnan(plane);
    if any(nanMask, 'all')
        plane(nanMask) = 0;
    else
        nanMask = [];
    end
end

out = imwarp(plane, tform, ...
    'OutputView', referenceView, ...
    'InterpolationMethod', 'linear', ...
    'FillValues', cast(0, class(plane)));

if ~isempty(nanMask)
    warpedNaN = imwarp(single(nanMask), tform, ...
        'OutputView', referenceView, ...
        'InterpolationMethod', 'nearest', ...
        'FillValues', single(0)) > 0.5;
    out(warpedNaN) = NaN;
end
end

% =========================================================================
% UMT transformation
% =========================================================================

function iTransformUmtFile(plan, tform, referenceView, tmpPath)
%ITRANSFORMUMTFILE Write the aligned copy of one UMT MAT file to tmpPath.

S = load(plan.Path,'-mat');
isStructLayout = strcmp(plan.UmtVariableName, iStructLayout());
if isStructLayout
    umt = S;
else
    umt = S.(plan.UmtVariableName);
end

validateUMTStruct(umt, 'requireEventInfo', false);

for iEntry = 1:numel(plan.EntryNames)
    entryName = plan.EntryNames{iEntry};
    entry = umt.data.(entryName);
    dimNames = iNormalizeDimNamesForPlan(entry.dimNames);

    entry.value = iTransformSpatialArrayGeneric(entry.value, dimNames, tform, referenceView);
    umt.data.(entryName) = entry;
end

umt = iRemoveResizedSpatialLabels(umt);

validateUMTStruct(umt, 'requireEventInfo', false);

if isStructLayout
    % Same layout and format as saveData (save2umt).
    save(tmpPath, '-struct', 'umt', '-mat', '-v7.3');
else
    S.(plan.UmtVariableName) = umt;
    save(tmpPath, '-struct', 'S', '-mat');
end
end

function name = iStructLayout()
%ISTRUCTLAYOUT Marker for a .umt whose variables are the UMT fields.

name = '<struct>';
end

function out = iTransformSpatialArrayGeneric(value, dimNames, tform, referenceView)
%ITRANSFORMSPATIALARRAYGENERIC Warp [Y X trailingDims...] arrays.
%
% Same plane warp as the .dat files (iWarpPlane). Numeric values are warped
% in single and cast back to their class; logical values use nearest.

valueClass = class(value);
isLogical = islogical(value);
if ~isLogical
    value = single(value);
end

newSizeYX = [referenceView.ImageSize(1), referenceView.ImageSize(2)];

oldSize = size(value);
nDims = numel(dimNames);

if numel(oldSize) < nDims
    oldSize(end+1:nDims) = 1;
end

trailingSize = oldSize(3:nDims);
if isempty(trailingSize)
    trailingSize = 1;
end

nPages = prod(trailingSize);
value3D = reshape(value, oldSize(1), oldSize(2), nPages);

out3D = zeros([newSizeYX, nPages], 'like', value);

for iPage = 1:nPages
    out3D(:, :, iPage) = iWarpPlane(value3D(:, :, iPage), tform, referenceView);
end

if nDims == 2
    out = out3D(:, :, 1);
else
    out = reshape(out3D, [newSizeYX, trailingSize]);
end

if ~isLogical
    out = cast(out, valueClass);
end
end

function umt = iRemoveResizedSpatialLabels(umt)
%IREMOVERESIZEDSPATIALLABELS Remove stale Y/X labels after spatial resizing.

if ~isfield(umt, 'labels') || ~isstruct(umt.labels)
    return
end

if isfield(umt.labels, 'Y')
    umt.labels = rmfield(umt.labels, 'Y');
end

if isfield(umt.labels, 'X')
    umt.labels = rmfield(umt.labels, 'X');
end

if isempty(fieldnames(umt.labels))
    umt = rmfield(umt, 'labels');
end
end

function [varName, umt] = iFindSingleUmtVariable(S, fileName)
%IFINDSINGLEUMTVARIABLE Find the UMT structure in loaded MAT content.
%
%   saveData writes a .umt with save -struct, so its variables ARE the UMT
%   fields; varName is then StructLayout. Older files hold one UMT-looking
%   variable; varName is its name. '' when no UMT structure is found.

varName = '';
umt = [];

if all(ismember({'version', 'kind', 'data'}, fieldnames(S)))
    varName = iStructLayout();
    umt = S;
    return
end

names = fieldnames(S);
isUMT = false(numel(names), 1);

for iName = 1:numel(names)
    x = S.(names{iName});

    isUMT(iName) = isstruct(x) && isscalar(x) && ...
        all(ismember({'version', 'kind', 'data'}, fieldnames(x)));
end

idx = find(isUMT);

if isempty(idx)
    return
end

if numel(idx) > 1
    error('ImageAlignmentTool:MultipleUmtVariables', ...
        ['UMT file "%s" contains multiple UMT variables. Multi-variable UMT ' ...
         'files are not supported by ImageAlignmentTool yet.'], ...
        fileName);
end

varName = names{idx};
umt = S.(varName);
end

function tf = iIsAllowedImageSpatialPattern(dimNames, schema)
%IISALLOWEDIMAGESPATIALPATTERN True for schema-allowed image pattern with Y/X.

if numel(dimNames) < 2 || ~strcmp(dimNames{1}, 'Y') || ~strcmp(dimNames{2}, 'X')
    tf = false;
    return
end

tf = any(cellfun(@(pattern) isequal(dimNames, pattern), schema.allowedPatterns.image));
end

function dimNames = iNormalizeDimNamesForPlan(dimNamesIn)
%INORMALIZEDIMNAMESFORPLAN Normalize dimNames for local plan decisions.

if isempty(dimNamesIn)
    dimNames = {};
elseif isstring(dimNamesIn)
    dimNames = cellstr(dimNamesIn(:).');
elseif iscell(dimNamesIn)
    dimNames = cellstr(string(dimNamesIn(:).'));
else
    error('ImageAlignmentTool:InvalidDimNames', ...
        'UMT entry dimNames must be a string vector or cell array of text.');
end
end

% =========================================================================
% ROI transformation
% =========================================================================

function iTransformRoiFile(plan, tform, referenceView, updatedView, tmpPath)
%ITRANSFORMROIFILE Write the aligned copy of a current-format ROI file to tmpPath.
%
% Every ROI becomes a polygon through its transformed vertices; its mask and
% geometric stats are recomputed at the new frame size. Intensity stats and
% the stats image refer to the old frame and are cleared. The file follows
% getROISchema and is checked with validateROIFile before it is written.

S = load(plan.Path, 'ROIFile','-mat');
ROIFile = S.ROIFile;

newSizeYX = [referenceView.ImageSize(1), referenceView.ImageSize(2)];
ROIs = ROIFile.ROIs;

fieldsToRemove = intersect(fieldnames(ROIs), {'ID', 'runtime'});
if ~isempty(fieldsToRemove)
    ROIs = rmfield(ROIs, fieldsToRemove);
end

pixelSize = updatedView.pixelSize_px_per_mm;
stampOn = datetime('now');
for iROI = 1:numel(ROIs)
    ROI = ROIs(iROI);
    vertices = double(ROI.geometry.verticesXY_px);

    [xNew, yNew] = transformPointsForward(tform, vertices(:,1), vertices(:,2));
    newVertices = [xNew(:), yNew(:)];

    ROI.type = 'polygon';
    ROI.modifiedOn = stampOn;
    ROI.geometry.verticesXY_px = newVertices;
    ROI.geometry.polyshape = polyshape(newVertices(:,1), newVertices(:,2), ...
        'Simplify', true);
    ROI.geometry.ROIType = 'polygon';
    ROI.geometry.ROIParameters = struct('ROIType', 'polygon', ...
        'Position', newVertices, 'Vertices', newVertices);
    ROI.mask = poly2mask(newVertices(:,1), newVertices(:,2), ...
        newSizeYX(1), newSizeYX(2));
    ROI.stats = iBuildRoiStats(ROI.stats, ROI.mask, pixelSize, stampOn);

    ROIs(iROI) = ROI;
end

ROIFile.ROIs = ROIs;
ROIFile.imageInfo.imageSizeYX = newSizeYX;
ROIFile.view.originXY_px = double(updatedView.origin_xy_px(:).');
ROIFile.view.pixelSize_px_per_mm = pixelSize;
ROIFile.statsImage.image = [];
ROIFile.modifiedOn = stampOn;

validateROIFile(ROIFile);
save(tmpPath, 'ROIFile', '-mat');
end

function tf = iIsCurrentRoiFile(ROIFile)
%IISCURRENTROIFILE True for a current UMIT ROI file (getROISchema).

schema = getROISchema();
tf = isstruct(ROIFile) && isscalar(ROIFile) && ...
    isfield(ROIFile, 'schemaName') && strcmp(char(string(ROIFile.schemaName)), schema.schemaName) && ...
    isfield(ROIFile, 'imageInfo') && isfield(ROIFile, 'ROIs') && isstruct(ROIFile.ROIs);
end

function stats = iBuildRoiStats(stats, mask, pixelSize, computedOn)
%IBUILDROISTATS Geometric ROI stats at the new frame size.
%
% Intensity stats (spatial*) and distances from the origin are cleared;
% the tools that own them recompute them from the aligned data.

stats.computedOn = computedOn;
stats.NPixels = nnz(mask);
stats.areaPx2 = nnz(mask);
stats.areaMM2 = [];
if isnumeric(pixelSize) && isscalar(pixelSize) && pixelSize > 0
    stats.areaMM2 = stats.areaPx2 / pixelSize^2;
end
stats.centroidXY_px = [];
if any(mask(:))
    [y, x] = find(mask);
    stats.centroidXY_px = [mean(x), mean(y)];
end
for fieldName = {'centroidXY_mm', 'distanceFromOrigin_px', 'distanceFromOrigin_mm', ...
        'spatialMean', 'spatialStd', 'spatialMedian', 'spatialMin', 'spatialMax'}
    stats.(fieldName{1}) = [];
end
end

% =========================================================================
% DataParams update
% =========================================================================

function updatedView = iBuildUpdatedView(DataParams, ImageReference, newSizeYX, scaleInfo, finalTform)
%IBUILDUPDATEDVIEW Build registered view metadata.
%
%   The image size is set to the reference image size. Spatial calibration is
%   propagated by the similarity scale. The view origin is transformed from
%   the old moving-image coordinate system into the registered/reference
%   coordinate system.

if isstruct(DataParams) && isfield(DataParams, 'view') && isstruct(DataParams.view)
    updatedView = DataParams.view;
else
    updatedView = struct();
end

updatedView.imageSizeYX = newSizeYX;

% Transform the existing origin instead of resetting it.
if isfield(updatedView, 'origin_xy_px') && ...
        ~isempty(updatedView.origin_xy_px) && ...
        isnumeric(updatedView.origin_xy_px) && ...
        numel(updatedView.origin_xy_px) >= 2

    oldOriginXY = double(updatedView.origin_xy_px(1:2));
else
    oldOriginXY = [1 1];
end

[xNew, yNew] = transformPointsForward(finalTform, oldOriginXY(1), oldOriginXY(2));
updatedView.origin_xy_px = [xNew yNew];

if isstruct(ImageReference) && isfield(ImageReference, 'axisConvention')
    updatedView.axisConvention = ImageReference.axisConvention;
elseif ~isfield(updatedView, 'axisConvention')
    updatedView.axisConvention = 'imageXY_topLeft';
end

if isfield(updatedView, 'pixelSize_px_per_mm') && ...
        ~isempty(updatedView.pixelSize_px_per_mm) && ...
        isnumeric(updatedView.pixelSize_px_per_mm)

    updatedView.pixelSize_px_per_mm = updatedView.pixelSize_px_per_mm .* scaleInfo.scale;
end
end

function DataParams = iUpdateDataParams(DataParams, finalTform, referenceImage, opts, ...
    precomputedMask, scaleInfo, updatedView)
%IUPDATEDATAPARAMS DataParams after alignment (validated, not yet saved).

DataParams.view = updatedView;

if ~isempty(precomputedMask)
    DataParams.mask.logical = logical(precomputedMask);
    DataParams.mask.space = 'registered';
    DataParams.mask.source = char(opts.createdBy);
end

DataParams.registration = iUpdateRegistrationStruct( ...
    DataParams.registration, finalTform, referenceImage, opts, scaleInfo);
end

function registration = iUpdateRegistrationStruct(registration, finalTform, ...
    referenceImage, opts, scaleInfo)
%IUPDATEREGISTRATIONSTRUCT Version-2 record after the apply.
%
%   The apply record is always written, committing the folder to the hash
%   of the reference it was aligned to (appliedReferenceHash). With
%   keepEstimate (pipeline apply of createRegistrationTform's result), the
%   reference and estimate records are kept: they describe exactly the
%   applied transform. Otherwise (GUI) both are rebuilt from the options;
%   nothing is carried over from an earlier estimate.

appliedOn = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
createdBy = char(opts.createdBy);
[~, referenceTemplate, estimateTemplate] = newRegistrationRecord();

if opts.keepEstimate
    % Checked before any write by iAssertEstimateMatchesReference.
else
    reference = referenceTemplate;
    [reference.contentHash, reference.image] = computeImageContentHash(referenceImage);
    reference.imageSizeYX = double(size(reference.image));
    reference.selectedOn = appliedOn;
    reference.source = iReferenceSource(opts);
    reference.name = 'Reference';
    if isstruct(opts.imageReference)
        if isfield(opts.imageReference, 'name') && ~isempty(opts.imageReference.name)
            reference.name = char(string(opts.imageReference.name));
        end
        if isfield(opts.imageReference, 'description')
            reference.description = char(string(opts.imageReference.description));
        end
        if isfield(opts.imageReference, 'axisConvention') && ...
                ~isempty(opts.imageReference.axisConvention)
            reference.axisConvention = char(string(opts.imageReference.axisConvention));
        end
    elseif ~isempty(reference.source.path)
        [~, fileName, fileExt] = fileparts(reference.source.path);
        reference.name = [fileName fileExt];
    end
    registration.reference = reference;

    estimate = estimateTemplate;
    estimate.tform = finalTform;
    estimate.transformType = 'similarity';
    estimate.method = createdBy;
    estimate.qcMetrics = iBuildTransformQcMetrics(finalTform, scaleInfo);
    estimate.createdOn = appliedOn;
    estimate.movingSource = char(opts.movingSourceFile);
    estimate.referenceHash = reference.contentHash;
    registration.estimate = estimate;
end

apply = registration.apply;
apply.isRegistered = true;
apply.isReviewed = logical(opts.isReviewed);
apply.appliedOn = appliedOn;
apply.appliedBy = createdBy;
apply.confirmationMode = char(opts.confirmationMode);
apply.appliedReferenceHash = registration.estimate.referenceHash;
if ~opts.keepEstimate
    apply.manualAdjustment = opts.manualAdjustment;
end
apply.notes = sprintf('Folder-level registration applied by %s.', createdBy);
registration.apply = apply;
end

function iAssertEstimateMatchesReference(folder, registration, operation)
%IASSERTESTIMATEMATCHESREFERENCE Apply guard (brief section 1.5).
%   The stored estimate must have been computed against the currently
%   selected reference; the hashes are compared, never the sources.

if isempty(registration.reference) || isempty(registration.estimate)
    message = UMITLedgerCore.buildErrorMessage(operation, folder, ...
        'the folder has no selected reference and estimate to apply.', ...
        {'Select a reference and estimate the alignment first (createRegistrationTform).'});
    throwAsCaller(MException('Umitoolbox:ImageRegistration:noEstimate', '%s', message));
end
if strcmp(registration.estimate.referenceHash, registration.reference.contentHash)
    return
end
why = sprintf(['the stored transform was estimated against %s, but the selected ' ...
    'reference is now %s.'], ...
    describeRegistrationReference([], registration.estimate.referenceHash), ...
    describeRegistrationReference(registration.reference));
message = UMITLedgerCore.buildErrorMessage(operation, folder, why, ...
    {['Estimate the alignment again against the selected reference (free before ' ...
    'apply), then apply.'], 'Or select again the reference the transform was estimated against.'});
throwAsCaller(MException('Umitoolbox:ImageRegistration:estimateReferenceMismatch', '%s', message));
end

function descriptor = iReferenceSource(opts)
%IREFERENCESOURCE Source descriptor of the GUI-supplied reference.

if isstruct(opts.referenceSource)
    descriptor = opts.referenceSource;
    validateReferenceSourceDescriptor(descriptor, 'referenceSource option');
    return
end
ImageReference = opts.imageReference;
path = char(opts.referenceFile);
if isstruct(ImageReference) && isfield(ImageReference, 'absolutePath') && ...
        ~isempty(ImageReference.absolutePath)
    path = char(string(ImageReference.absolutePath));
end
if isstruct(ImageReference) && isfield(ImageReference, 'projectUUID') && ...
        ~isempty(ImageReference.projectUUID)
    descriptor = referenceSourceDescriptor('project', path);
    for fieldName = {'projectUUID', 'projectName', 'subjectUUID', 'subjectID', ...
            'sessionUUID', 'sessionID', 'resourceUUID'}
        if isfield(ImageReference, fieldName{1})
            descriptor.(fieldName{1}) = char(string(ImageReference.(fieldName{1})));
        end
    end
else
    descriptor = referenceSourceDescriptor('file', path);
end
end

function qc = iBuildTransformQcMetrics(tform, scaleInfo)
%IBUILDTRANSFORMQCMETRICS Basic transform metrics.

T = iGetTransformMatrix(tform);

qc = struct();
qc.MIBefore = [];
qc.MIAfter = [];
qc.MIDelta = [];
qc.translationXY_px = [T(3,1), T(3,2)];
qc.rotationDeg = scaleInfo.rotationDeg;
qc.scale = scaleInfo.scale;
qc.scaleXY = [scaleInfo.scaleX, scaleInfo.scaleY];
qc.isUniformScale = scaleInfo.isUniformScale;
qc.hasNoShear = scaleInfo.hasNoShear;
qc.determinant = det(T(1:2, 1:2));
end

% =========================================================================
% Staged commit / failure report
% =========================================================================

function staged = iNewStagedList()
%INEWSTAGEDLIST Empty list of staged (or installed) files.

staged = struct('finalPath', {}, 'tmpPath', {}, 'backupPath', {}, 'hadOriginal', {});
end

function staged = iAddStaged(staged, finalPath, tmpPath)
%IADDSTAGED Record a staged file before it is written (so it is cleaned up).

staged(end+1) = struct('finalPath', finalPath, 'tmpPath', tmpPath, ...
    'backupPath', [finalPath '.aligning.bak'], 'hadOriginal', isfile(finalPath));
end

function tmpPath = iStagePath(finalPath)
%ISTAGEPATH Temporary sibling of finalPath. Its extension is .tmp, so a
%leftover is never mistaken for a .dat, .umt or .roi file.

tmpPath = [finalPath '.aligning.tmp'];
end

function installed = iInstallStaged(installed, entry)
%IINSTALLSTAGED Move the original aside, then install the staged copy.

if entry.hadOriginal
    [ok, message] = movefile(entry.finalPath, entry.backupPath, 'f');
    if ~ok
        error('ImageAlignmentTool:CommitFailed', ...
            'Could not move "%s" aside before installing its aligned copy: %s', ...
            entry.finalPath, message);
    end
end
installed(end+1) = entry;
[ok, message] = movefile(entry.tmpPath, entry.finalPath, 'f');
if ~ok
    error('ImageAlignmentTool:CommitFailed', ...
        'Could not install the aligned copy of "%s": %s', entry.finalPath, message);
end
end

function message = iRollBack(installed)
%IROLLBACK Restore every original moved aside by iInstallStaged.

problems = {};
for iFile = numel(installed):-1:1
    entry = installed(iFile);
    try
        if entry.hadOriginal
            if isfile(entry.backupPath)
                movefile(entry.backupPath, entry.finalPath, 'f');
            end
        elseif isfile(entry.finalPath)
            delete(entry.finalPath);
        end
    catch ME
        problems{end+1} = sprintf('%s (%s)', entry.finalPath, ME.message); %#ok<AGROW>
    end
end
if isempty(problems)
    message = 'yes';
else
    message = sprintf('NO, failed for: %s', strjoin(problems, '; '));
end
end

function iDiscardStaged(staged)
%IDISCARDSTAGED Delete staged temporary files that were not installed.

for iFile = 1:numel(staged)
    if isfile(staged(iFile).tmpPath)
        delete(staged(iFile).tmpPath);
    end
end
end

function iDeleteBackups(installed)
%IDELETEBACKUPS Drop the originals once the commit has succeeded.

for iFile = 1:numel(installed)
    if installed(iFile).hadOriginal && isfile(installed(iFile).backupPath)
        delete(installed(iFile).backupPath);
    end
end
end

function failureInfo = iFailureInfo(ME, report, dateTag, stage)
%IFAILUREINFO Fields of the failure report.

failureInfo = struct('error', ME, 'report', report, 'dateTag', dateTag, ...
    'stage', stage);
end

function reportPath = iWriteFailureReport(targetFolder, failureInfo)
%IWRITEFAILUREREPORT Write a human-readable failure report.

dateTag = failureInfo.dateTag;
reportPath = fullfile(targetFolder, ['ImageAlignmentFailedReport_' dateTag '.txt']);

fid = fopen(reportPath, 'w');
if fid < 0
    warning('ImageAlignmentTool:CouldNotWriteFailureReport', ...
        'Could not write failure report: %s', reportPath);
    return
end

cleanupObj = onCleanup(@() safeFclose(fid));

ME = failureInfo.error;

fprintf(fid, 'Image Alignment Failed Report\n');
fprintf(fid, '=============================\n\n');

fprintf(fid, 'Created on: %s\n', char(datetime('now')));
fprintf(fid, 'Target folder: %s\n', failureInfo.report.targetFolder);
fprintf(fid, 'Failed during: %s\n\n', failureInfo.stage);

fprintf(fid, 'Reference file: %s\n', failureInfo.report.referenceFile);
fprintf(fid, 'Old image size [Y X]: %s\n', mat2str(failureInfo.report.oldSizeYX));
fprintf(fid, 'New image size [Y X]: %s\n\n', mat2str(failureInfo.report.newSizeYX));

fprintf(fid, 'Files staged before the failure\n');
fprintf(fid, '-------------------------------\n');
iWriteStringList(fid, 'DAT files', failureInfo.report.datFilesTransformed);
iWriteStringList(fid, 'UMT files', failureInfo.report.umtFilesTransformed);
iWriteStringList(fid, 'ROI files', failureInfo.report.roiFilesTransformed);

fprintf(fid, '\nWarnings\n');
fprintf(fid, '--------\n');
iWriteStringList(fid, 'Warnings', failureInfo.report.warnings);

fprintf(fid, '\nError\n');
fprintf(fid, '-----\n');
fprintf(fid, 'Identifier: %s\n', ME.identifier);
fprintf(fid, 'Message: %s\n\n', ME.message);

fprintf(fid, 'Full MATLAB report\n');
fprintf(fid, '------------------\n');
fprintf(fid, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end

function iWriteStringList(fid, titleText, values)
%IWRITESTRINGLIST Write a simple string list to text report.

fprintf(fid, '%s:\n', titleText);

if isempty(values)
    fprintf(fid, '  <none>\n');
    return
end

values = string(values(:));
for iValue = 1:numel(values)
    fprintf(fid, '  - %s\n', values(iValue));
end
end

% =========================================================================
% Validation / utility helpers
% =========================================================================

function referenceImage = iValidateReferenceImage(referenceImage)
%IVALIDATEREFERENCEIMAGE Validate and convert reference image.

if ~(isnumeric(referenceImage) || islogical(referenceImage)) || ...
        ~ismatrix(referenceImage) || isempty(referenceImage)
    error('ImageAlignmentTool:InvalidReferenceImage', ...
        'referenceImage must be a non-empty 2D numeric/logical image.');
end

referenceImage = single(referenceImage);
referenceImage(~isfinite(referenceImage)) = 0;
end

function tform = iNormalizeTransform(tform)
%INORMALIZETRANSFORM Validate supported 2D geometric transform input.

if isnumeric(tform)
    if ~isequal(size(tform), [3 3]) || any(~isfinite(tform(:)))
        error('ImageAlignmentTool:InvalidTransformMatrix', ...
            'Numeric finalTform must be a finite 3 x 3 matrix.');
    end

    tform = affine2d(tform);
    return
end

supportedClasses = { ...
    'affine2d', ...
    'affinetform2d', ...
    'simtform2d', ...
    'rigidtform2d', ...
    'transltform2d'};

isSupported = any(cellfun(@(c) isa(tform, c), supportedClasses));

if ~isSupported
    error('ImageAlignmentTool:InvalidTransform', ...
        ['finalTform must be a supported 2D geometric transform. ' ...
         'Received class: %s'], class(tform));
end

T = iGetTransformMatrix(tform);

if ~isequal(size(T), [3 3]) || any(~isfinite(T(:)))
    error('ImageAlignmentTool:InvalidTransformMatrix', ...
        'Transform matrix must be a finite 3 x 3 matrix.');
end
end

function T = iGetTransformMatrix(tform)
%IGETTRANSFORMMATRIX Return 3 x 3 matrix from old/new MATLAB transform.

if isnumeric(tform)
    T = tform;
    return
end

if isprop(tform, 'T')
    T = tform.T;
    return
end

if isprop(tform, 'A')
    T = tform.A;
    return
end

error('ImageAlignmentTool:MissingTransformMatrix', ...
    'Could not extract transform matrix from object of class "%s".', class(tform));
end

function scaleInfo = iGetSimilarityScaleInfo(tform)
%IGETSIMILARITYSCALEINFO Extract and validate uniform similarity scale.

T = iGetTransformMatrix(tform);
A = T(1:2, 1:2);

scaleX = hypot(A(1,1), A(1,2));
scaleY = hypot(A(2,1), A(2,2));
scale = mean([scaleX, scaleY]);

tol = 1e-6 * max([1, scaleX, scaleY]);
isUniformScale = abs(scaleX - scaleY) <= tol;

colDot = dot(A(:,1), A(:,2));
hasNoShear = abs(colDot) <= 1e-6 * max(1, norm(A(:,1)) * norm(A(:,2)));

rotationDeg = atan2d(A(1,2), A(1,1));

scaleInfo = struct();
scaleInfo.scale = scale;
scaleInfo.scaleX = scaleX;
scaleInfo.scaleY = scaleY;
scaleInfo.isUniformScale = isUniformScale;
scaleInfo.hasNoShear = hasNoShear;
scaleInfo.isSimilarity = isUniformScale && hasNoShear && isfinite(scale) && scale > 0;
scaleInfo.rotationDeg = rotationDeg;
end

function oldSizeYX = iResolveOldSizeYX(~, DataParams, datFiles, datInfos)
%IRESOLVEOLDSIZEYX Resolve old [Y X] image size before alignment.
%
% The .dat files' own frame size comes first (they must all share it),
% then DataParams.view.imageSizeYX. AcqInfos.mat is not used: its
% Height/Width describe the raw acquisition, not the imported data (.dat
% header Phase 7a).

described = ~cellfun(@isempty, datInfos);
if any(described)
    sizes = cellfun(@(info) [datAxisSize(info, 'Y'), datAxisSize(info, 'X')], ...
        datInfos(described), 'UniformOutput', false);
    sizes = vertcat(sizes{:});
    if size(unique(sizes, 'rows'), 1) > 1
        names = {datFiles(described).name};
        error('ImageAlignmentTool:InconsistentDatSizes', ...
            ['The .dat files in the folder have different frame sizes (%s). ' ...
             'No file was changed.'], ...
            strjoin(cellfun(@(n, k) sprintf('%s: %d x %d', n, sizes(k, 1), sizes(k, 2)), ...
            names, num2cell(1:numel(names)), 'UniformOutput', false), '; '));
    end
    oldSizeYX = sizes(1, :);
    return
end

if isstruct(DataParams) && isfield(DataParams, 'view') && ...
        isfield(DataParams.view, 'imageSizeYX') && ~isempty(DataParams.view.imageSizeYX)
    oldSizeYX = double(DataParams.view.imageSizeYX(:).');
    return
end

error('ImageAlignmentTool:CouldNotResolveOldSize', ...
    ['Could not resolve old image size from the .dat files or ' ...
     'DataParams.view.imageSizeYX.']);
end
