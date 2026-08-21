function [classifiedFiles, rigResolution] = run_ImagesClassification(RawFolder, SaveFolder, varargin)
%RUN_IMAGESCLASSIFICATION Classify raw interlaced IOS binaries into channels.
%
%   classifiedFiles = run_ImagesClassification(RawFolder, SaveFolder)
%   classifiedFiles = run_ImagesClassification(RawFolder, SaveFolder, 'Name', Value, ...)
%   [classifiedFiles, rigResolution] = run_ImagesClassification(...)
%   info = run_ImagesClassification('pipelineInfo')
%
%   This wrapper calls ImagesClassification to split raw interlaced binary
%   data into separate channel .dat files and a shared AcqInfos.mat file.
%   For dual-camera acquisitions, it automatically applies the resolved
%   Rig's active camera-coregistration resource when one is configured.
%
%   Inputs:
%       RawFolder  - Path to the folder containing the raw binary
%                    acquisition.
%       SaveFolder - Path to the folder where processed channel files will
%                    be saved.
%
%   Name-Value parameters:
%       BinningSpatial - Spatial binning factor.
%                        Default: 1
%
%       BinningTemp    - Temporal binning factor.
%                        Default: 1
%
%       backupOpts     - Backup handling option passed to
%                        ImagesClassification before writing outputs into
%                        SaveFolder.
%                        Values:
%                            ''           resolve interactively
%                            'ERASE'      delete without a backup
%                            'GENBACKUP'  zip, then delete
%                            <char>       zip under that base name, then
%                                         delete
%                        Default: '' for a direct call; 'GENBACKUP' when the
%                        function is driven by PipelineManager, which cannot
%                        answer an interactive prompt and must not default to
%                        a destructive action.
%
%                        Notes:
%                        - All values delete the managed existing files in
%                          SaveFolder; only 'ERASE' does so without first
%                          creating a .zip backup.
%                        - Raw acquisition files (.bin, .tif, info.txt,
%                          Comments.txt, Snapshot*.png) are never managed.
%
%   Outputs:
%       classifiedFiles - SaveFolder-relative file manifest of outputs saved
%                        in SaveFolder.
%       rigResolution  - Struct containing rigUUID, rigID, wasCreated, the
%                        default-Rig resolution path, and a nested
%                        cameraCoregistration struct. The nested struct
%                        reports status ('applied' or 'skipped'),
%                        wasApplied, skipReason, resourceUUID, displayName,
%                        fileName, managedFilePath, rigUUID, rigID, and
%                        whether DataParams provenance was recorded.
%
%   Notes:
%       - This function does not expose SubROI selection. It always calls
%         ImagesClassification with b_SubROI = 0.
%       - The returned file manifest uses SaveFolder-relative names.
%       - An active resource is never inferred from available or archived
%         Rig resources. A Rig with no active pointer imports without
%         coregistration.
%       - Invalid active-resource metadata or a missing managed file is an
%         error; import never silently substitutes another resource.
%
%   Examples:
%       classifiedFiles = run_ImagesClassification(rawFolder, saveFolder);
%
%       classifiedFiles = run_ImagesClassification( ...
%           rawFolder, saveFolder, ...
%           'BinningSpatial', 2, ...
%           'BinningTemp', 4);
%
%       classifiedFiles = run_ImagesClassification( ...
%           rawFolder, saveFolder, ...
%           'backupOpts', 'GENBACKUP');

default_Output = {'fluo_475.dat', 'fluo_567.dat', 'fluo.dat', ...
    'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat', 'AcqInfos.mat'};
rigResolution = [];
allowedBinning = 1:8;

if nargin == 1 && (ischar(RawFolder) || (isstring(RawFolder) && isscalar(RawFolder))) ...
        && strcmpi(strtrim(char(string(RawFolder))), 'pipelineInfo')
    classifiedFiles = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'BinningSpatial', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && any(x == allowedBinning));
addParameter(p, 'BinningTemp', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && any(x == allowedBinning));
addParameter(p, 'backupOpts', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, RawFolder, SaveFolder, varargin{:});

RawFolder = char(string(p.Results.RawFolder));
SaveFolder = char(string(p.Results.SaveFolder));
BinningSpatial = p.Results.BinningSpatial;
BinningTemp = p.Results.BinningTemp;
backupOpts = char(string(p.Results.backupOpts));

% Resolve an active Rig before the legacy importer mutates SaveFolder.
% ImagesClassification itself remains independent of UMITRigStore so the
% IOIAnalysis folder retains its legacy execution path.
[defaultRigStore, wasCreated, resolution] = UMITRigStore.getOrCreateDefaultRig();
defaultRigInfo = defaultRigStore.getRigInfo();
rigResolution = struct( ...
    'rigUUID', defaultRigInfo.uuid, ...
    'rigID', defaultRigInfo.rigID, ...
    'wasCreated', wasCreated, ...
    'resolution', resolution, ...
    'cameraCoregistration', localEmptyCoregistrationResolution(defaultRigInfo));

% Resolve and validate the active pointer before classification changes the
% destination folder. getActiveCameraCoregistration deliberately returns
% empty when no pointer is configured and throws for inconsistent pointers
% or missing managed files; those errors must not be converted into a skip.
try
    activeResource = defaultRigStore.getActiveCameraCoregistration();
    if ~isempty(activeResource)
        tformFile = defaultRigStore.resolveResourcePath(activeResource.uuid);
        tformPayload = localLoadCoregistrationResource(tformFile);
        rigResolution.cameraCoregistration = localResourceResolution( ...
            rigResolution.cameraCoregistration, activeResource, tformFile);
    end
catch ME
    coregError = MException( ...
        'Umitoolbox:run_ImagesClassification:InvalidActiveCameraCoregistration', ...
        ['The active camera-coregistration resource for Rig "%s" (%s) ' ...
         'could not be validated: %s'], ...
        defaultRigInfo.rigID, defaultRigInfo.uuid, ME.message);
    coregError = addCause(coregError, ME);
    throwAsCaller(coregError)
end

classifiedFiles = ImagesClassification( ...
    RawFolder, ...
    SaveFolder, ...
    BinningSpatial, ...
    BinningTemp, ...
    0, ...
    'backupOpts', backupOpts);
defaultRigID = defaultRigInfo.rigID;

if ~iscell(classifiedFiles)
    classifiedFiles = {};
end

% Bind the imported acquisition to the resolved Rig, then apply its active
% camera-coregistration resource to compatible dual-camera data.
acqInfoPath = fullfile(SaveFolder, 'AcqInfos.mat');
if ~isfile(acqInfoPath)
    error('Umitoolbox:run_ImagesClassification:InvalidAcqInfos', ...
        'ImagesClassification did not produce AcqInfos.mat.');
end

info = load(acqInfoPath, 'AcqInfoStream');
if ~isfield(info, 'AcqInfoStream') || ~isstruct(info.AcqInfoStream) || ...
        ~isscalar(info.AcqInfoStream)
    error('Umitoolbox:run_ImagesClassification:InvalidAcqInfos', ...
        'ImagesClassification did not produce a scalar AcqInfoStream.');
end
info.AcqInfoStream.rigUUID = defaultRigInfo.uuid;
info.AcqInfoStream.rigID = defaultRigInfo.rigID;
saveMatAtomic(acqInfoPath, 'AcqInfoStream', info.AcqInfoStream);

isMultiCamera = isfield(info.AcqInfoStream, 'MultiCam') && ...
    isscalar(info.AcqInfoStream.MultiCam) && logical(info.AcqInfoStream.MultiCam);
if isMultiCamera
    disp('Dual Camera data found!')
    if isempty(activeResource)
        rigResolution.cameraCoregistration.status = 'skipped';
        rigResolution.cameraCoregistration.skipReason = 'noActiveResource';
        fprintf(['No active camera-coregistration resource is configured for ' ...
            'Rig "%s" (%s); import will continue without coregistration.\n'], ...
            defaultRigInfo.rigID, defaultRigInfo.uuid);
    else
        fprintf(['Applying camera coregistration "%s" from managed file "%s" ' ...
            '(resource UUID %s), owned by Rig "%s" (UUID %s).\n'], ...
            activeResource.displayName, activeResource.fileName, activeResource.uuid, ...
            defaultRigInfo.rigID, defaultRigInfo.uuid);
        [status, warnmsg] = applyTform2Cams( ...
            SaveFolder, tformPayload.tform, tformPayload.tformInfo);
        if ~status
            rigResolution.cameraCoregistration.status = 'skipped';
            rigResolution.cameraCoregistration.skipReason = 'applicationFailed';
            warning('Umitoolbox:run_ImagesClassification:cameraCoregistrationFailed', ...
                ['Coregistration failed; data import will resume without ' ...
                 'coregistration. %s'], warnmsg);
        else
            rigResolution.cameraCoregistration.status = 'applied';
            rigResolution.cameraCoregistration.wasApplied = true;
            disp('Done!')

            % Record this application in the imported SaveFolder's
            % DataParams.mat. Best-effort: a failure here must not undo an
            % otherwise-successful coregistration.
            try
                if ~isfile(fullfile(SaveFolder, 'DataParams.mat'))
                    createDataParams(SaveFolder);
                end

                DataParams = loadDataParams(SaveFolder);
                cc = DataParams.cameraCoregistration;

                cc.isCoregistered = true;
                cc.isReviewed = false;
                cc.tform = tformPayload.tform;
                cc.resourceUUID = activeResource.uuid;
                cc.rigID = defaultRigID;
                cc.method = 'run_ImagesClassification';
                cc.sourceFile = tformFile;
                if isfield(tformPayload.tformInfo, 'SavedOn')
                    cc.sourceFileTimestamp = char(tformPayload.tformInfo.SavedOn);
                else
                    cc.sourceFileTimestamp = '';
                end
                cc.createdOn = char(datetime('now'));
                cc.source = SaveFolder;
                cc.appliedOn = char(datetime('now'));
                cc.appliedBy = 'run_ImagesClassification';
                cc.confirmationMode = 'automatic-import-application';
                cc.notes = ['Dual-camera coregistration automatically applied during ' ...
                    'raw data import using the rig''s active camera coregistration.'];

                % getTformQCMetrics normalises the premultiply/postmultiply
                % convention, so the recorded metrics do not depend on which
                % transform class the saved resource happens to hold.
                tformQC = getTformQCMetrics(tformPayload.tform);
                cc.transformType = tformQC.transformType;
                cc.qcMetrics.translationXY_px = tformQC.translationXY;
                cc.qcMetrics.rotationDeg = tformQC.rotationDeg;
                cc.qcMetrics.scaleXY = tformQC.scaleXY;
                cc.qcMetrics.determinant = tformQC.determinant;

                DataParams.cameraCoregistration = cc;
                DataParams.lastModified = datetime('now');

                saveDataParams(SaveFolder, DataParams);
                rigResolution.cameraCoregistration.dataParamsRecorded = true;
            catch ME
                warning('Umitoolbox:run_ImagesClassification:cameraCoregistrationDataParamsFailed', ...
                    ['Coregistration succeeded, but DataParams.cameraCoregistration could not ' ...
                    'be recorded: %s'], ME.message);
            end
        end
    end
else
    rigResolution.cameraCoregistration.status = 'skipped';
    rigResolution.cameraCoregistration.skipReason = 'singleCameraAcquisition';
end

% Classification has fully succeeded by this point (any failure above throws
% or asserts, aborting the function before reaching here). Ensure DataParams.mat
% exists and records the RawFolder actually used for this import, so dual-camera
% coregistration and other RawFolder-dependent features are available without a
% separate manual "Set Raw Folder" step. Best-effort: a failure here must not
% undo an otherwise-successful import.
try
    if ~isfile(fullfile(SaveFolder, 'DataParams.mat'))
        createDataParams(SaveFolder);
    end
    updateDataParam(SaveFolder, 'folders.RawFolder', RawFolder);
catch ME
    warning('Umitoolbox:run_ImagesClassification:dataParamsRawFolderFailed', ...
        'Import succeeded, but DataParams.mat RawFolder could not be recorded: %s', ME.message);
end

% Return SaveFolder-relative names so the manifest remains relocatable,
% consistent with funcTemplateFileManifest's convention.
classifiedFiles = unique(classifiedFiles, 'stable');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Classify raw interlaced IOS binary data into channel .dat files.');
        info.version = '1.0.0';
        info.freshSaveFolderRole = 'acquisition-initializer';

        info = PipelineManager.addInput(info, 'RawFolder', 'RawFolder', ...
            'Path to the folder containing the raw binary acquisition.', ...
            'kind', 'input', 'position', 1, 'callType', 'positional', 'isData', false);

        info = PipelineManager.addInput(info, 'SaveFolder', 'SaveFolder', ...
            'Path to the folder where processed channel files will be saved.', ...
            'kind', 'input', 'position', 2, 'callType', 'positional', 'isData', false);

        info = PipelineManager.addInput(info, 'BinningSpatial', 'parameter', ...
            'Spatial binning factor.', ...
            'kind', 'parameter', 'default', 1, 'allowed', allowedBinning, 'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'BinningTemp', 'parameter', ...
            'Temporal binning factor.', ...
            'kind', 'parameter', 'default', 1, 'allowed', allowedBinning, 'callType', 'namevalue');

        info = PipelineManager.addInput(info, 'backupOpts', 'parameter', ...
            ['Backup handling for managed files already present in SaveFolder. ' ...
             '''GENBACKUP'' zips them before the import; ''ERASE'' deletes them ' ...
             'without a backup.'], ...
            'kind', 'parameter', 'default', 'GENBACKUP', 'allowed', {'GENBACKUP','ERASE'}, 'callType', 'namevalue');

        info = PipelineManager.addOutput(info, 'classifiedFiles', 'ImageTimeSeries', 'file', ...
            'Generated file manifest saved in SaveFolder.', ...
            default_Output, 1, 'isData', true, 'saveFileName', '');

        info = PipelineManager.addOutput(info, 'rigResolution', 'UnknownDataType', 'data', ...
            ['Rig/optical resolution metadata for the classified data (rigUUID, rigID, ' ...
             'wasCreated, resolution, cameraCoregistration).'], ...
            '', 2, 'isData', false);
    end

    function T = localGetTransformMatrix(tform)
        %LOCALGETTRANSFORMMATRIX Return 3 x 3 transform matrix from a tform object.

        if isprop(tform, 'T')
            T = tform.T;
        elseif isprop(tform, 'A')
            T = tform.A;
        else
            T = tform;
        end
    end

    function result = localEmptyCoregistrationResolution(rigInfo)
        %LOCALEMPTYCOREGISTRATIONRESOLUTION Create a stable caller-facing result.

        result = struct( ...
            'status', 'pending', ...
            'wasApplied', false, ...
            'skipReason', '', ...
            'resourceUUID', '', ...
            'displayName', '', ...
            'fileName', '', ...
            'managedFilePath', '', ...
            'rigUUID', rigInfo.uuid, ...
            'rigID', rigInfo.rigID, ...
            'dataParamsRecorded', false);
    end

    function result = localResourceResolution(result, resource, resourcePath)
        %LOCALRESOURCERESOLUTION Copy active managed-resource identity.

        result.resourceUUID = resource.uuid;
        result.displayName = resource.displayName;
        result.fileName = resource.fileName;
        result.managedFilePath = resourcePath;
    end

    function payload = localLoadCoregistrationResource(resourcePath)
        %LOCALLOADCOREGISTRATIONRESOURCE Validate the managed MAT payload.

        payload = load(resourcePath, '-mat');
        if ~isfield(payload, 'tform') || isempty(payload.tform)
            error('Umitoolbox:run_ImagesClassification:InvalidCoregistrationPayload', ...
                'Managed file must contain a nonempty variable named "tform".');
        end
        payload.tform = localNormalizeCoregistrationTransform(payload.tform);
        if ~isfield(payload, 'tformInfo')
            payload.tformInfo = struct();
        elseif ~isstruct(payload.tformInfo) || ~isscalar(payload.tformInfo)
            error('Umitoolbox:run_ImagesClassification:InvalidCoregistrationPayload', ...
                'Variable "tformInfo" must be a scalar struct when present.');
        end

        T = localGetTransformMatrix(payload.tform);
        if ~isnumeric(T) || ~isequal(size(T), [3 3]) || any(~isfinite(T), 'all')
            error('Umitoolbox:run_ImagesClassification:InvalidCoregistrationPayload', ...
                'The active camera-coregistration transform must have a finite 3-by-3 matrix.');
        end
    end

    function tform = localNormalizeCoregistrationTransform(tform)
        %LOCALNORMALIZECOREGISTRATIONTRANSFORM Convert supported transforms for legacy use.
        %   applyTform2Cams uses affine2d's postmultiply T convention. Newer
        %   MATLAB transform classes expose a premultiply A matrix, so that
        %   matrix must be transposed during conversion.

        if isa(tform, 'affine2d')
            return
        end

        if isnumeric(tform)
            if ~isequal(size(tform), [3 3]) || any(~isfinite(tform), 'all')
                error('Umitoolbox:run_ImagesClassification:InvalidCoregistrationPayload', ...
                    'A numeric camera-coregistration transform must be a finite 3-by-3 matrix.');
            end
            tform = affine2d(double(tform));
            return
        end

        modernTransformClasses = { ...
            'affinetform2d', 'simtform2d', 'rigidtform2d', 'transltform2d'};
        isModernTransform = any(cellfun( ...
            @(className) isa(tform, className), modernTransformClasses));
        if ~isModernTransform || ~isprop(tform, 'A')
            error('Umitoolbox:run_ImagesClassification:InvalidCoregistrationPayload', ...
                ['Variable "tform" must be affine2d, a supported modern 2-D ' ...
                 'affine transform, or a finite numeric 3-by-3 matrix.']);
        end

        matrix = double(tform.A);
        if ~isequal(size(matrix), [3 3]) || any(~isfinite(matrix), 'all')
            error('Umitoolbox:run_ImagesClassification:InvalidCoregistrationPayload', ...
                'The modern camera-coregistration transform must have a finite 3-by-3 A matrix.');
        end
        tform = affine2d(matrix.');
    end
end
