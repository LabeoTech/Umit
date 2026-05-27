function report = applyImageAlignmentToFolder(targetFolder, finalTform, opts)
%APPLYIMAGEALIGNMENTTOFOLDER Apply ImageReference registration to a folder.
%
%   report = applyImageAlignmentToFolder(targetFolder, finalTform, ...
%       'imageReference', ImageReference, ...
%       'referenceFile', referenceFileRel, ...
%       'referenceImage', referenceImage)
%
%   Applies a folder-level geometric registration to all supported image-space
%   files in targetFolder.
%
%   Processing steps:
%       1) Validate target folder, metadata, .dat, .umt, and .roi files.
%       2) Create a managed-folder backup using genBackupFolder.
%       3) Transform all .dat files using RAM-safe frame-by-frame I/O.
%       4) Transform eligible image-kind .umt files.
%       5) Transform current-format .roi files and force transformed ROIs to
%          polygon geometry.
%       6) Update DataParams.mat.
%       7) Update AcqInfos.mat.
%       8) Save manualAlignmentReport_<datetime>.mat.
%
%   Transform convention:
%       finalTform maps the current raw moving-data coordinate space into
%       ImageReference coordinate space.
%
%   Name-Value options:
%       imageReference       - ImageReference struct used as fixed reference.
%       referenceFile        - Relative ImageReference file path.
%       referenceImage       - Fixed reference image. If empty, this is read
%                              from imageReference.image.
%       baseTform            - Automatic/base transform before manual tuning.
%       previewTform         - Preview transform used in the GUI.
%       manualAdjustment     - Struct with dx/dy/rotation/scale.
%       movingSourceFile     - File or context used to estimate the transform.
%       preFlipHorizontally  - Whether the GUI pre-flipped the moving source.
%       createdBy            - Tool/app name.
%       backupChoice         - Backup name or 'GENBACKUP'.
%       transformDat         - Transform .dat files. Default: true.
%       transformUmt         - Transform eligible image .umt files. Default: true.
%       transformRoi         - Transform current-format .roi files. Default: true.
%
%   Output:
%       report - Compact operation report.
%
%   Notes:
%       This function expects .dat files to be raw single-precision image
%       time series with dimensions [Y X T]. The frame count is inferred from
%       file size and validated before any file is modified.
%
%       genBackupFolder intentionally ignores raw acquisition files such as
%       .bin and .tif. This function therefore creates a managed-folder
%       backup, not a raw-data backup.

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
    opts.backupChoice (1,1) string = "GENBACKUP"
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

[DataParams, dataParamsPath] = iLoadDataParams(targetFolder);
[AcqInfoStream, acqInfosPath, acqInfosStruct] = iLoadAcqInfos(targetFolder);

oldSizeYX = iResolveOldSizeYX(AcqInfoStream, DataParams, targetFolder);
report.oldSizeYX = oldSizeYX;

scaleInfo = iGetSimilarityScaleInfo(finalTform);
if ~scaleInfo.isSimilarity
    error('ImageAlignmentTool:NonSimilarityTransform', ...
        ['Registration transform is not a valid similarity transform. ' ...
         'Uniform scale/no-shear validation failed.']);
end
report.scaleInfo = scaleInfo;

datFiles = dir(fullfile(targetFolder, '*.dat'));
umtFiles = dir(fullfile(targetFolder, '*.umt'));
roiFiles = dir(fullfile(targetFolder, '*.roi'));

if isempty(datFiles) && isempty(umtFiles)
    error('ImageAlignmentTool:NoImageFilesFound', ...
        'No .dat or .umt image files were found in target folder.');
end

% Build the complete plan before backup or mutation. This keeps temporary
% files generated during execution out of the operation list.
iPreflightDataParams(DataParams, oldSizeYX);

datPlan = iPreflightDatFiles(datFiles, oldSizeYX, opts.transformDat);
[umtPlan, umtWarnings] = iPreflightUmtFiles(umtFiles, opts.transformUmt);
[roiPlan, roiWarnings] = iPreflightRoiFiles(roiFiles, opts.transformRoi);

report.warnings = [report.warnings; umtWarnings(:); roiWarnings(:)];

% Precompute mask transformation before touching data files. This catches
% imwarp/mask errors before any destructive write.
precomputedMask = iPrecomputeTransformedMask( ...
    DataParams, oldSizeYX, finalTform, referenceView);

backupChoice = char(string(opts.backupChoice));
if strcmpi(backupChoice, 'GENBACKUP')
    backupChoice = ['alignment_bkp_' dateTag];
end

try
    backupPath = genBackupFolder(targetFolder, backupChoice, 'eraseFolder', false);
catch ME
    error('ImageAlignmentTool:BackupFailed', ...
        'Backup failed. No files were modified. Original error: %s', ...
        ME.message);
end

report.backupPath = backupPath;
report.backupCreatedOn = datetime('now');

completed = iInitializeCompletedTracker();

try
    % Transform .dat files with RAM-safe direct binary I/O.
    if opts.transformDat
        for iFile = 1:numel(datPlan)
            iTransformDatFile(datPlan(iFile), finalTform, referenceView);
            report.datFilesTransformed(end+1,1) = string(datPlan(iFile).Name); %#ok<AGROW>
            completed.datFiles(end+1,1) = string(datPlan(iFile).Name); %#ok<AGROW>
        end
    end

    % Transform eligible image-kind .umt files.
    if opts.transformUmt
        for iFile = 1:numel(umtPlan)
            if ~umtPlan(iFile).Transform
                report.umtFilesSkipped(end+1,1) = ...
                    string(umtPlan(iFile).Name) + ": " + string(umtPlan(iFile).Reason); %#ok<AGROW>
                continue
            end

            iTransformUmtFile(umtPlan(iFile), finalTform, referenceView);
            report.umtFilesTransformed(end+1,1) = string(umtPlan(iFile).Name); %#ok<AGROW>
            completed.umtFiles(end+1,1) = string(umtPlan(iFile).Name); %#ok<AGROW>
        end
    end
    
    updatedView = iBuildUpdatedView(DataParams, opts.imageReference, newSizeYX, scaleInfo, finalTform);

    % Transform current-format .roi files.
    if opts.transformRoi
        for iFile = 1:numel(roiPlan)
            if ~roiPlan(iFile).Transform
                report.roiFilesSkipped(end+1,1) = ...
                    string(roiPlan(iFile).Name) + ": " + string(roiPlan(iFile).Reason); %#ok<AGROW>
                continue
            end

            iTransformRoiFile(roiPlan(iFile), finalTform, referenceView, updatedView);
            report.roiFilesTransformed(end+1,1) = string(roiPlan(iFile).Name); %#ok<AGROW>
            completed.roiFiles(end+1,1) = string(roiPlan(iFile).Name); %#ok<AGROW>
        end
    end

    % Update DataParams after all file transforms succeed.
    DataParams = iUpdateDataParams( ...
        DataParams, ...
        finalTform, ...
        referenceImage, ...
        opts.imageReference, ...
        opts.referenceFile, ...
        opts.movingSourceFile, ...
        opts.manualAdjustment, ...
        opts.createdBy, ...
        precomputedMask, ...
        scaleInfo);

    save(dataParamsPath, 'DataParams', '-mat');
    report.metadataFilesUpdated(end+1,1) = "DataParams.mat";
    completed.metadataFiles(end+1,1) = "DataParams.mat";

    % Update AcqInfos after all file transforms succeed.
    AcqInfoStream.Height = newSizeYX(1);
    AcqInfoStream.Width = newSizeYX(2);
    acqInfosStruct.AcqInfoStream = AcqInfoStream;
    save(acqInfosPath, '-struct', 'acqInfosStruct');

    report.metadataFilesUpdated(end+1,1) = "AcqInfos.mat";
    completed.metadataFiles(end+1,1) = "AcqInfos.mat";

    report.finishedOn = datetime('now');
    reportPath = fullfile(targetFolder, ['manualAlignmentReport_' dateTag '.mat']);
    report.reportPath = reportPath;
    save(reportPath, 'report', '-mat');

catch ME
    failureInfo = struct();
    failureInfo.error = ME;
    failureInfo.completed = completed;
    failureInfo.backupPath = backupPath;
    failureInfo.dateTag = dateTag;
    failureInfo.report = report;

    [restoreSucceeded, restoreMessage] = iRestoreFromBackup(targetFolder, backupPath);

    failureInfo.restoreAttempted = true;
    failureInfo.restoreSucceeded = restoreSucceeded;
    failureInfo.restoreMessage = restoreMessage;

    iWriteFailureReport(targetFolder, failureInfo);

    rethrow(ME);
end

end

% =========================================================================
% Report / metadata loading
% =========================================================================

function report = iInitializeReport(targetFolder, finalTform, opts, newSizeYX, createdOn)
%IINITIALIZEREPORT Create compact operation report.

report = struct();
report.createdOn = createdOn;
report.finishedOn = NaT;
report.toolName = 'ImageAlignmentTool';
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
report.backupPath = '';
report.backupCreatedOn = NaT;
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

function [DataParams, dataParamsPath] = iLoadDataParams(targetFolder)
%ILOADDATAPARAMS Load required DataParams.mat.

dataParamsPath = fullfile(targetFolder, 'DataParams.mat');

if ~isfile(dataParamsPath)
    error('ImageAlignmentTool:MissingDataParams', ...
        'DataParams.mat is required before applying registration.');
end

S = load(dataParamsPath, 'DataParams');

if ~isfield(S, 'DataParams')
    error('ImageAlignmentTool:InvalidDataParamsFile', ...
        'DataParams.mat does not contain variable DataParams.');
end

DataParams = S.DataParams;
end

function [AcqInfoStream, acqInfosPath, S] = iLoadAcqInfos(targetFolder)
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

function datPlan = iPreflightDatFiles(datFiles, oldSizeYX, doTransform)
%IPREFLIGHTDATFILES Validate DAT file size and infer frame count.

datPlan = struct( ...
    'Name', {}, ...
    'Path', {}, ...
    'OldSizeYX', {}, ...
    'NFrames', {}, ...
    'FileBytes', {});

if ~doTransform
    return
end

for iFile = 1:numel(datFiles)
    datPath = fullfile(datFiles(iFile).folder, datFiles(iFile).name);
    info = dir(datPath);

    if isempty(info)
        error('ImageAlignmentTool:MissingDatFile', ...
            'DAT file disappeared during preflight: %s', datPath);
    end

    bytesPerFrame = prod(oldSizeYX) * 4; % single precision
    fileBytes = info.bytes;

    if bytesPerFrame <= 0 || mod(fileBytes, bytesPerFrame) ~= 0
        error('ImageAlignmentTool:InvalidDatFileSize', ...
            ['DAT file "%s" is incompatible with expected [Y X T] single layout. ' ...
             'File bytes: %d. Expected bytes per frame: %d.'], ...
            datFiles(iFile).name, fileBytes, bytesPerFrame);
    end

    nFrames = fileBytes / bytesPerFrame;

    if nFrames < 1 || mod(nFrames, 1) ~= 0
        error('ImageAlignmentTool:InvalidDatFrameCount', ...
            'DAT file "%s" has invalid frame count inferred from file size.', ...
            datFiles(iFile).name);
    end

    datPlan(end+1) = struct( ... %#ok<AGROW>
        'Name', datFiles(iFile).name, ...
        'Path', datPath, ...
        'OldSizeYX', oldSizeYX, ...
        'NFrames', nFrames, ...
        'FileBytes', fileBytes);
end
end

function [umtPlan, warnings] = iPreflightUmtFiles(umtFiles, doTransform)
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
        umtPlan(end+1) = struct( ... %#ok<AGROW>
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', false, ...
            'Reason', "no UMT structure found", ...
            'UmtVariableName', '', ...
            'EntryNames', {{}});
        continue
    end

    validateUMTStruct(umt, 'requireEventInfo', false);

    if ~isfield(umt, 'kind') || ~strcmpi(char(string(umt.kind)), 'image')
        umtPlan(end+1) = struct( ... %#ok<AGROW>
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', false, ...
            'Reason', "kind is not image", ...
            'UmtVariableName', umtVarName, ...
            'EntryNames', {{}});
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
            transformEntries{end+1} = entryName; %#ok<AGROW>
        else
            skippedEntries(end+1,1) = string(entryName) + " (" + strjoin(string(dimNames), "/") + ")"; %#ok<AGROW>
        end
    end

    if isempty(transformEntries)
        umtPlan(end+1) = struct( ... %#ok<AGROW>
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', false, ...
            'Reason', "no eligible image entries", ...
            'UmtVariableName', umtVarName, ...
            'EntryNames', {{}});
    else
        if ~isempty(skippedEntries)
            warnings(end+1,1) = sprintf( ... %#ok<AGROW>
                'UMT file "%s": skipped unsupported image entries: %s.', ...
                umtFiles(iFile).name, strjoin(skippedEntries, ', '));
        end

        umtPlan(end+1) = struct( ... %#ok<AGROW>
            'Name', umtFiles(iFile).name, ...
            'Path', umtPath, ...
            'Transform', true, ...
            'Reason', "", ...
            'UmtVariableName', umtVarName, ...
            'EntryNames', {transformEntries});
    end
end
end

function [roiPlan, warnings] = iPreflightRoiFiles(roiFiles, doTransform)
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
        warnings(end+1,1) = sprintf( ... %#ok<AGROW>
            'ROI file "%s" has unsupported legacy schema and will be skipped.', ...
            roiFiles(iFile).name);

        roiPlan(end+1) = struct( ... %#ok<AGROW>
            'Name', roiFiles(iFile).name, ...
            'Path', roiPath, ...
            'Transform', false, ...
            'Reason', "unsupported legacy schema");
        continue
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

    roiPlan(end+1) = struct( ... %#ok<AGROW>
        'Name', roiFiles(iFile).name, ...
        'Path', roiPath, ...
        'Transform', true, ...
        'Reason', "");
end
end

% =========================================================================
% DAT transformation
% =========================================================================

function iTransformDatFile(plan, tform, referenceView)
%ITRANSFORMDATFILE Transform one .dat file with RAM-safe direct binary I/O.

datPath = plan.Path;
oldSizeYX = plan.OldSizeYX;
nFrames = plan.NFrames;

tmpPath = iMakeTempSiblingPath(datPath);

fidIn = fopen(datPath, 'r');
if fidIn < 0
    error('ImageAlignmentTool:CouldNotOpenDatInput', ...
        'Could not open DAT input file: %s', datPath);
end
cleanupIn = onCleanup(@() safeFclose(fidIn));

fidOut = fopen(tmpPath, 'w');
if fidOut < 0
    error('ImageAlignmentTool:CouldNotOpenDatOutput', ...
        'Could not open temporary DAT output file: %s', tmpPath);
end
cleanupOut = onCleanup(@() safeFclose(fidOut));

try
    for iFrame = 1:nFrames
        frame = fread(fidIn, oldSizeYX, 'single=>single');

        if ~isequal(size(frame), oldSizeYX)
            error('ImageAlignmentTool:UnexpectedDatReadSize', ...
                'Unexpected read size in DAT file "%s" at frame %d.', ...
                datPath, iFrame);
        end

        outFrame = imwarp(frame, tform, ...
            'OutputView', referenceView, ...
            'InterpolationMethod', 'linear', ...
            'FillValues', single(0));

        nWritten = fwrite(fidOut, single(outFrame), 'single');

        if nWritten ~= numel(outFrame)
            error('ImageAlignmentTool:DatWriteFailed', ...
                'Could not write full transformed frame %d to "%s".', ...
                iFrame, tmpPath);
        end
    end

    safeFclose(fidIn);
    safeFclose(fidOut);
    clear cleanupIn cleanupOut

    movefile(tmpPath, datPath, 'f');

catch ME
    safeFclose(fidIn);
    safeFclose(fidOut);

    if isfile(tmpPath)
        delete(tmpPath);
    end

    rethrow(ME);
end
end

% =========================================================================
% UMT transformation
% =========================================================================

function iTransformUmtFile(plan, tform, referenceView)
%ITRANSFORMUMTFILE Transform planned entries in one UMT MAT file.

S = load(plan.Path,'-mat');
umt = S.(plan.UmtVariableName);

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

S.(plan.UmtVariableName) = umt;

tmpPath = iMakeTempSiblingPath(plan.Path);

try
    save(tmpPath, '-struct', 'S', '-mat');
    movefile(tmpPath, plan.Path, 'f');

catch ME
    if isfile(tmpPath)
        delete(tmpPath);
    end
    rethrow(ME);
end
end

function out = iTransformSpatialArrayGeneric(value, dimNames, tform, referenceView)
%ITRANSFORMSPATIALARRAYGENERIC Warp [Y X trailingDims...] arrays.

valueClass = class(value);
isLogical = islogical(value);

newSizeYX = [referenceView.ImageSize(1), referenceView.ImageSize(2)];

valueSingle = single(value);
oldSize = size(valueSingle);
nDims = numel(dimNames);

if numel(oldSize) < nDims
    oldSize(end+1:nDims) = 1;
end

trailingSize = oldSize(3:nDims);
if isempty(trailingSize)
    trailingSize = 1;
end

nPages = prod(trailingSize);
value3D = reshape(valueSingle, oldSize(1), oldSize(2), nPages);

out3D = zeros([newSizeYX, nPages], 'single');

if isLogical
    interpMethod = 'nearest';
else
    interpMethod = 'linear';
end

for iPage = 1:nPages
    out3D(:, :, iPage) = imwarp(value3D(:, :, iPage), tform, ...
        'OutputView', referenceView, ...
        'InterpolationMethod', interpMethod, ...
        'FillValues', single(0));
end

if nDims == 2
    out = out3D(:, :, 1);
else
    out = reshape(out3D, [newSizeYX, trailingSize]);
end

if isLogical
    out = out > 0.5;
else
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
%IFINDSINGLEUMTVARIABLE Find one UMT-looking variable in loaded MAT content.

varName = '';
umt = [];

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

function iTransformRoiFile(plan, tform, referenceView, updatedView)
%ITRANSFORMROIFILE Transform current-format ROIFile MAT file.

S = load(plan.Path, 'ROIFile','-mat');
ROIFile = S.ROIFile;

newSizeYX = [referenceView.ImageSize(1), referenceView.ImageSize(2)];
ROIs = ROIFile.ROIs;

fieldsToRemove = intersect(fieldnames(ROIs), {'runtime', 'persistent'});
if ~isempty(fieldsToRemove)
    ROIs = rmfield(ROIs, fieldsToRemove);
end

for iROI = 1:numel(ROIs)
    ROI = ROIs(iROI);
    vertices = double(ROI.geometry.verticesXY_px);

    [xNew, yNew] = transformPointsForward(tform, vertices(:,1), vertices(:,2));
    newVertices = [xNew(:), yNew(:)];

    ROI.geometry.verticesXY_px = newVertices;
    ROI.geometry.polyshape = polyshape(newVertices(:,1), newVertices(:,2), ...
        'Simplify', true);
    ROI.geometry.ROIType = 'polygon';
    ROI.geometry.ROIParameters = [];

    ROI.mask = poly2mask(newVertices(:,1), newVertices(:,2), ...
        newSizeYX(1), newSizeYX(2));

    ROI.stats = iBuildRoiStats(ROI.mask, newVertices);

    if isfield(ROI, 'type')
        ROI.type = 'polygon';
    end

    if isfield(ROI, 'modifiedOn')
        ROI.modifiedOn = datetime('now');
    end

    ROIs(iROI) = ROI;
end

ROIFile.ROIs = ROIs;
ROIFile.imageSizeYX = newSizeYX;

if isfield(ROIFile, 'DataParamsView')
    ROIFile.DataParamsView = updatedView;
end

if isfield(ROIFile, 'statsImage')
    ROIFile.statsImage = [];
end

if isfield(ROIFile, 'statsImageViewMode')
    ROIFile.statsImageViewMode = 'registered';
end

tmpPath = iMakeTempSiblingPath(plan.Path);

try
    save(tmpPath, 'ROIFile', '-mat');
    movefile(tmpPath, plan.Path, 'f');

catch ME
    if isfile(tmpPath)
        delete(tmpPath);
    end
    rethrow(ME);
end
end

function tf = iIsCurrentRoiFile(ROIFile)
%IISCURRENTROIFILE True for current supported ROIFile schema.

tf = isstruct(ROIFile) && isscalar(ROIFile) && ...
    isfield(ROIFile, 'fileType') && strcmpi(char(string(ROIFile.fileType)), 'UMIT_ROI_FILE') && ...
    isfield(ROIFile, 'ROIs') && isstruct(ROIFile.ROIs);
end

function stats = iBuildRoiStats(mask, verticesXY)
%IBUILDROISTATS Minimal ROI stats after geometric transform.

stats = struct();
stats.area_px = nnz(mask);
stats.verticesXY_px = verticesXY;

if any(mask(:))
    [y, x] = find(mask);
    stats.centroidXY_px = [mean(x), mean(y)];
    stats.boundingBoxXYWH = [min(x), min(y), max(x)-min(x)+1, max(y)-min(y)+1];
else
    stats.centroidXY_px = [NaN NaN];
    stats.boundingBoxXYWH = [NaN NaN NaN NaN];
end
end

% =========================================================================
% DataParams / AcqInfos update
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


function DataParams = iUpdateDataParams(DataParams, finalTform, referenceImage, ...
    ImageReference, referenceFile, movingSourceFile, manualAdjustment, createdBy, precomputedMask, scaleInfo)
%IUPDATEDATAPARAMS Update folder-level DataParams after registration.

newSizeYX = size(referenceImage);

DataParams.view = iBuildUpdatedView(DataParams, ImageReference, newSizeYX, scaleInfo, finalTform);

if ~isempty(precomputedMask)
    if ~isfield(DataParams, 'mask') || ~isstruct(DataParams.mask)
        DataParams.mask = struct();
    end

    DataParams.mask.logical = logical(precomputedMask);
    DataParams.mask.space = 'registered';
    DataParams.mask.source = 'ImageAlignmentTool';
end

DataParams.registration = iUpdateRegistrationStruct( ...
    DataParams, ...
    finalTform, ...
    referenceImage, ...
    ImageReference, ...
    referenceFile, ...
    movingSourceFile, ...
    manualAdjustment, ...
    createdBy, ...
    scaleInfo);

DataParams.lastModified = datetime('now');
end

function registration = iUpdateRegistrationStruct(DataParams, finalTform, referenceImage, ...
    ImageReference, referenceFile, movingSourceFile, manualAdjustment, createdBy, scaleInfo)
%IUPDATEREGISTRATIONSTRUCT Preserve and update DataParams.registration.

if isfield(DataParams, 'registration') && isstruct(DataParams.registration)
    registration = DataParams.registration;
else
    registration = struct();
end

registration.isRegistered = true;
registration.isReviewed = true;
registration.tform = finalTform;
registration.transformType = 'similarity';
registration.method = 'ImageAlignmentTool';
registration.referenceFile = char(string(referenceFile));
registration.referenceImage = referenceImage;

if isstruct(ImageReference) && isfield(ImageReference, 'description')
    registration.referenceDescription = ImageReference.description;
else
    registration.referenceDescription = '';
end

registration.source = char(string(movingSourceFile));
registration.createdOn = char(datetime('now'));
registration.appliedOn = char(datetime('now'));
registration.appliedBy = char(string(createdBy));
registration.confirmationMode = 'user-confirmed-folder-transform';
registration.notes = 'Folder-level registration applied by ImageAlignmentTool.';
registration.manualAdjustment = manualAdjustment;
registration.qcMetrics = iBuildTransformQcMetrics(finalTform, scaleInfo);
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
% Restore / failure report
% =========================================================================

function completed = iInitializeCompletedTracker()
%IINITIALIZECOMPLETEDTRACKER Track files modified before failure.

completed = struct();
completed.datFiles = strings(0,1);
completed.umtFiles = strings(0,1);
completed.roiFiles = strings(0,1);
completed.metadataFiles = strings(0,1);
end

function [restoreSucceeded, restoreMessage] = iRestoreFromBackup(targetFolder, backupPath)
%IRESTOREFROMBACKUP Restore managed files from backup zip.

restoreSucceeded = false;
restoreMessage = "";

if isempty(backupPath) || ~isfile(backupPath)
    restoreMessage = "Backup path is empty or backup file was not found.";
    return
end

try
    iDeleteManagedFilesForRestore(targetFolder, backupPath);
    unzip(backupPath, targetFolder);

    restoreSucceeded = true;
    restoreMessage = "Restore completed from backup zip.";

catch ME
    restoreSucceeded = false;
    restoreMessage = string(ME.message);
end
end

function iDeleteManagedFilesForRestore(targetFolder, backupPath)
%IDELETEMANAGEDFILESFORRESTORE Delete non-raw managed files before restore.
%
%   The backup zip itself and all other zip files are protected.

backupPath = char(string(backupPath));
backupPath = strrep(backupPath, '/', filesep);
backupPath = strrep(backupPath, '\', filesep);

allFiles = dir(targetFolder);
allFiles([allFiles.isdir]) = [];

protectedFiles = false(numel(allFiles), 1);

for iFile = 1:numel(allFiles)
    filePath = fullfile(allFiles(iFile).folder, allFiles(iFile).name);
    filePath = strrep(filePath, '/', filesep);
    filePath = strrep(filePath, '\', filesep);

    name = allFiles(iFile).name;

    protectedFiles(iFile) = ...
        strcmpi(filePath, backupPath) || ...
        endsWith(name, '.zip', 'IgnoreCase', true) || ...
        endsWith(name, '.bin', 'IgnoreCase', true) || ...
        endsWith(name, '.tif', 'IgnoreCase', true) || ...
        endsWith(name, '.tiff', 'IgnoreCase', true) || ...
        endsWith(name, 'info.txt', 'IgnoreCase', true) || ...
        strcmpi(name, 'Comments.txt') || ...
        (startsWith(name, 'Snapshot', 'IgnoreCase', true) && ...
         endsWith(name, '.png', 'IgnoreCase', true)) || ...
        startsWith(name, 'ImageAlignmentFailedReport_', 'IgnoreCase', true);
end

deleteFiles = allFiles(~protectedFiles);

for iFile = 1:numel(deleteFiles)
    filePath = fullfile(deleteFiles(iFile).folder, deleteFiles(iFile).name);

    if isfile(filePath)
        delete(filePath);
    end
end
end

function reportPath = iWriteFailureReport(targetFolder, failureInfo)
%IWRITEFAILUREREPORT Write human-readable failure report after restore.

dateTag = failureInfo.dateTag;
reportPath = fullfile(targetFolder, ['ImageAlignmentFailedReport_' dateTag '.txt']);

fid = fopen(reportPath, 'w');
if fid < 0
    warning('ImageAlignmentTool:CouldNotWriteFailureReport', ...
        'Could not write failure report: %s', reportPath);
    return
end

cleanupObj = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

ME = failureInfo.error;

fprintf(fid, 'Image Alignment Failed Report\n');
fprintf(fid, '=============================\n\n');

fprintf(fid, 'Created on: %s\n', char(datetime('now')));
fprintf(fid, 'Target folder: %s\n', failureInfo.report.targetFolder);
fprintf(fid, 'Backup path: %s\n', failureInfo.backupPath);
fprintf(fid, 'Restore attempted: %d\n', failureInfo.restoreAttempted);
fprintf(fid, 'Restore succeeded: %d\n', failureInfo.restoreSucceeded);
fprintf(fid, 'Restore message: %s\n\n', char(string(failureInfo.restoreMessage)));

fprintf(fid, 'Reference file: %s\n', failureInfo.report.referenceFile);
fprintf(fid, 'Old image size [Y X]: %s\n', mat2str(failureInfo.report.oldSizeYX));
fprintf(fid, 'New image size [Y X]: %s\n\n', mat2str(failureInfo.report.newSizeYX));

fprintf(fid, 'Completed before failure\n');
fprintf(fid, '------------------------\n');
iWriteStringList(fid, 'DAT files', failureInfo.completed.datFiles);
iWriteStringList(fid, 'UMT files', failureInfo.completed.umtFiles);
iWriteStringList(fid, 'ROI files', failureInfo.completed.roiFiles);
iWriteStringList(fid, 'Metadata files', failureInfo.completed.metadataFiles);

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
        ndims(referenceImage) ~= 2 || isempty(referenceImage)
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

function oldSizeYX = iResolveOldSizeYX(AcqInfoStream, DataParams, targetFolder)
%IRESOLVEOLDSIZEYX Resolve old [Y X] image size before alignment.

if isstruct(AcqInfoStream) && ...
        isfield(AcqInfoStream, 'Height') && isfield(AcqInfoStream, 'Width')
    oldSizeYX = [double(AcqInfoStream.Height), double(AcqInfoStream.Width)];
    return
end

if isstruct(DataParams) && isfield(DataParams, 'view') && ...
        isfield(DataParams.view, 'imageSizeYX')
    oldSizeYX = double(DataParams.view.imageSizeYX(:).');
    return
end

datFiles = dir(fullfile(targetFolder, '*.dat'));
if isempty(datFiles)
    error('ImageAlignmentTool:CouldNotResolveOldSize', ...
        'Could not resolve old image size from AcqInfos, DataParams, or DAT files.');
end

error('ImageAlignmentTool:CouldNotResolveOldSize', ...
    ['Could not resolve old image size. AcqInfoStream.Height/Width or ' ...
     'DataParams.view.imageSizeYX is required.']);
end

function tmpPath = iMakeTempSiblingPath(filePath)
%IMAKETEMPSIBLINGPATH Create unique sibling temp path preserving extension.

[folderName, baseName, ext] = fileparts(filePath);

dateTag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
randomTag = sprintf('%06d', randi(999999));

tmpPath = fullfile(folderName, sprintf('%s_%s_%s%s', ...
    baseName, dateTag, randomTag, ext));

while isfile(tmpPath)
    randomTag = sprintf('%06d', randi(999999));
    tmpPath = fullfile(folderName, sprintf('%s_%s_%s%s', ...
        baseName, dateTag, randomTag, ext));
end
end
