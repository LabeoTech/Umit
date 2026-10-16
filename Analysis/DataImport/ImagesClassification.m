function [classifiedFiles, rigResolution] = ImagesClassification(RawFolder, SaveFolder, varargin)
%IMAGESCLASSIFICATION Classify raw interlaced IOS binaries into channels.
%
%   classifiedFiles = ImagesClassification(RawFolder, SaveFolder)
%   classifiedFiles = ImagesClassification(RawFolder, SaveFolder, 'Name', Value, ...)
%   [classifiedFiles, rigResolution] = ImagesClassification(...)
%   info = ImagesClassification('pipelineInfo')
%
%   This function reads raw interlaced .bin files produced by the LabeoTech
%   IOS systems, separates the illumination sequence into one .dat file per
%   unique output channel, and writes one shared AcqInfos.mat file in
%   SaveFolder. Rig information is recorded only in the folder rig ledger
%   (UMITRigLedger.json; see UMITRigLedger), never in AcqInfos.mat:
%       - Single-camera data resolve no rig. A first import records an
%         explicit "no rig"; a re-import keeps the existing ledger state.
%       - Dual-camera data use the Active Rig and its active
%         camera-coregistration resource. The transform is applied, the
%         folder is bound to the Rig's snapshots and the coregistration is
%         committed. Without an Active Rig or an active resource (or if the
%         transform cannot be applied), the data are imported WITHOUT
%         coregistration, with a warning: a new folder records "no rig", a
%         bound folder keeps its optics binding without a coregistration.
%       - Re-importing into a SaveFolder that already holds an import is a
%         FRESH IMPORT (prepareFreshImport): every artefact derived from the
%         previous data is deleted (planFreshImport lists them; raw files,
%         ledgers and unknown files are kept), and freshImport events are
%         recorded in the rig ledger (ALL rig commitments reset, binding
%         kept) and, when the folder is bound to a project, in the project
%         ledger. backupOpts is the only decision: GENBACKUP (default in
%         pipelines) zips the deleted files first, ERASE does not; a direct
%         call with backupOpts = '' asks on a desktop. A dual-camera
%         re-import then rebinds to the Active Rig, whatever rig the folder
%         used before.
%
%   Repeated illumination entries are supported within each camera sequence.
%   For example, Red-Green-Red creates red.dat and green.dat; red.dat stores
%   both Red sequence positions in chronological order and therefore has a
%   higher effective frame rate than green.dat. Repeated output channel tags
%   across cameras are rejected because the current file naming scheme cannot
%   represent the same output tag independently for multiple cameras.
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
%       backupOpts     - What to do when the SaveFolder already holds an
%                        import (a fresh import); ignored for an empty
%                        SaveFolder. Values:
%                            'GENBACKUP'  zip the files to delete, then
%                                         delete them and import
%                            'ERASE'      delete them without a zip, import
%                            <char>       zip under that base name, then
%                                         delete and import
%                            ''           direct call: confirmation dialog
%                                         listing the deletions (zip, delete
%                                         or cancel); GENBACKUP in batch
%                        Default: '' for a direct call; 'GENBACKUP' when the
%                        function is driven by PipelineManager.
%
%       b_SubROI       - Interactive sub-region selection before sorting
%                        (pre-defined ROI file or drawn rectangle). Not
%                        exposed to PipelineManager.
%                        Default: false
%
%       ApplyRig       - false classifies the raw data only: no Rig is
%                        resolved, the rig ledger is not read or written,
%                        no camera coregistration is applied, DataParams.mat
%                        is not touched, and rigResolution is empty. Used by
%                        calibration workflows that need unregistered source
%                        data. Not exposed to PipelineManager.
%                        Default: true
%
%   Outputs:
%       classifiedFiles - SaveFolder-relative file manifest of outputs saved
%                        in SaveFolder.
%       rigResolution  - Folder rig-ledger summary after the import:
%                        ledgerState, rigUUID, rigID, rigDisplayName,
%                        opticsHash and coregistrationHash (empty when the
%                        folder is not bound), and a nested
%                        cameraCoregistration struct reporting status
%                        ('applied' or 'skipped'), wasApplied, skipReason
%                        ('singleCameraAcquisition', 'noActiveRig',
%                        'noActiveResource' or 'applicationFailed'),
%                        resourceUUID, displayName, fileName,
%                        managedFilePath, and whether DataParams provenance
%                        was recorded.
%
%   Notes:
%       - The returned file manifest uses SaveFolder-relative names.
%       - Imported channel temporal metadata are stored in
%         AcqInfoStream.ImportedChannels. Spatial metadata remain at the top
%         level of AcqInfoStream.
%       - An active resource is never inferred from available or archived
%         Rig resources. A Rig with no active pointer imports without
%         coregistration.
%       - Invalid active-resource metadata or a missing managed file is an
%         error; import never silently substitutes another resource.
%       - Several Active Rigs, or an invalid folder rig ledger, are errors
%         raised before SaveFolder is touched. No Rig is ever created.
%       - AcqInfos.mat is written atomically.
%
%   Examples:
%       classifiedFiles = ImagesClassification(rawFolder, saveFolder);
%
%       classifiedFiles = ImagesClassification( ...
%           rawFolder, saveFolder, ...
%           'BinningSpatial', 2, ...
%           'BinningTemp', 4);
%
%       classifiedFiles = ImagesClassification( ...
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
addParameter(p, 'b_SubROI', false, ...
    @(x) isscalar(x) && (islogical(x) || ismember(x, [0 1])));
addParameter(p, 'ApplyRig', true, ...
    @(x) isscalar(x) && (islogical(x) || ismember(x, [0 1])));
parse(p, RawFolder, SaveFolder, varargin{:});

RawFolder = char(string(p.Results.RawFolder));
SaveFolder = char(string(p.Results.SaveFolder));
BinningSpatial = p.Results.BinningSpatial;
BinningTemp = p.Results.BinningTemp;
backupOpts = char(string(p.Results.backupOpts));
b_SubROI = logical(p.Results.b_SubROI);

if ~p.Results.ApplyRig
    % Calibration data into temporary folders: no ledger reads or writes.
    prepareFreshImport(SaveFolder, 'backupOpts', backupOpts, ...
        'Origin', 'ImagesClassification', 'RecordLedgerEvents', false);
    classifiedFiles = iClassifyRawBinaries(RawFolder, SaveFolder, ...
        BinningSpatial, BinningTemp, b_SubROI);
    return
end

% Rig information is recorded only in the folder rig ledger
% (UMITRigLedger.json), never in AcqInfos.mat. Everything rig-related is
% resolved and checked before classification touches SaveFolder.
ledgerOrigin = 'ImagesClassification';
rawInfo = ReadInfoFile(RawFolder);
isMultiCamera = isfield(rawInfo, 'MultiCam') && isscalar(rawInfo.MultiCam) && ...
    logical(rawInfo.MultiCam);
isReimport = UMITRigLedger.read(SaveFolder).exists; % throws for an invalid ledger
rigResolution = localNewRigResolution();

binding = [];
uncoregisteredReason = '';
if isMultiCamera
    [binding, activeResource, tformFile, tformPayload, uncoregisteredReason, ...
        rigResolution.cameraCoregistration.skipReason] = localResolveDualCameraRig();
end

% A re-import is a fresh import (confirmed or refused here, before any file
% is touched): derived artefacts are deleted and freshImport events reset
% every rig commitment and are logged in the project ledger when bound.
prepareFreshImport(SaveFolder, 'backupOpts', backupOpts, 'Origin', ledgerOrigin);

classifiedFiles = iClassifyRawBinaries(RawFolder, SaveFolder, ...
    BinningSpatial, BinningTemp, b_SubROI);
if ~iscell(classifiedFiles)
    classifiedFiles = {};
end
if ~isfile(fullfile(SaveFolder, 'AcqInfos.mat'))
    error('Umitoolbox:ImagesClassification:InvalidAcqInfos', ...
        'Classification did not produce AcqInfos.mat.');
end

if isMultiCamera
    disp('Dual Camera data found!')
    if ~isempty(binding)
        fprintf(['Applying camera coregistration "%s" (resource UUID %s, managed file ' ...
            '"%s") from Rig "%s" (%s, UUID %s).\n'], ...
            activeResource.displayName, activeResource.uuid, activeResource.fileName, ...
            binding.sourceRig.displayName, binding.sourceRig.rigID, binding.sourceRig.uuid);
        [status, warnmsg] = applyTform2Cams( ...
            SaveFolder, tformPayload.tform, tformPayload.tformInfo);
        if status
            disp('Done!')
            localBindCoregisteredFolder(binding);
            rigResolution.cameraCoregistration.status = 'applied';
            rigResolution.cameraCoregistration.wasApplied = true;
            rigResolution.cameraCoregistration.resourceUUID = activeResource.uuid;
            rigResolution.cameraCoregistration.displayName = activeResource.displayName;
            rigResolution.cameraCoregistration.fileName = activeResource.fileName;
            rigResolution.cameraCoregistration.managedFilePath = tformFile;
            rigResolution.cameraCoregistration.dataParamsRecorded = ...
                localRecordCoregistrationProvenance(binding, activeResource, tformFile, tformPayload);
        else
            rigResolution.cameraCoregistration.skipReason = 'applicationFailed';
            uncoregisteredReason = sprintf( ...
                'applying the camera-coregistration transform failed (%s)', warnmsg);
        end
    end
    if ~rigResolution.cameraCoregistration.wasApplied
        rigResolution.cameraCoregistration.status = 'skipped';
        localRecordUncoregisteredFolder(uncoregisteredReason);
    end
else
    rigResolution.cameraCoregistration.status = 'skipped';
    rigResolution.cameraCoregistration.skipReason = 'singleCameraAcquisition';
    % A first single-camera import records an explicit "no rig"; a
    % re-import keeps the existing binding (only the freshImport event above).
    if ~isReimport
        UMITRigLedger.writeNoRig(SaveFolder, 'Origin', ledgerOrigin);
    end
end
rigResolution = localSummarizeLedger(rigResolution);

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
    warning('Umitoolbox:ImagesClassification:dataParamsRawFolderFailed', ...
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
            ['When the SaveFolder already holds an import (a fresh import deletes every ' ...
             'artefact derived from the previous data and resets the rig commitments): ' ...
             '''GENBACKUP'' (default) zips the files it deletes, then re-imports; ' ...
             '''ERASE'' deletes them without a zip, then re-imports. Ignored for an ' ...
             'empty SaveFolder.'], ...
            'kind', 'parameter', 'default', 'GENBACKUP', 'allowed', {'GENBACKUP','ERASE'}, ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput(info, 'classifiedFiles', 'ImageTimeSeries', 'file', ...
            'Generated file manifest saved in SaveFolder.', ...
            default_Output, 1, 'isData', true, 'saveFileName', '');

        info = PipelineManager.addOutput(info, 'rigResolution', 'UnknownDataType', 'data', ...
            ['Folder rig-ledger summary after the import (ledgerState, rigUUID, rigID, ' ...
             'rigDisplayName, opticsHash, coregistrationHash, cameraCoregistration).'], ...
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

    function result = localNewRigResolution()
        %LOCALNEWRIGRESOLUTION Create a stable caller-facing ledger summary.

        result = struct( ...
            'ledgerState', '', ...
            'rigUUID', '', ...
            'rigID', '', ...
            'rigDisplayName', '', ...
            'opticsHash', '', ...
            'coregistrationHash', '', ...
            'cameraCoregistration', struct( ...
                'status', 'pending', ...
                'wasApplied', false, ...
                'skipReason', '', ...
                'resourceUUID', '', ...
                'displayName', '', ...
                'fileName', '', ...
                'managedFilePath', '', ...
                'dataParamsRecorded', false));
    end

    function [rigBinding, resource, resourcePath, payload, reasonText, skipCode] = ...
            localResolveDualCameraRig()
        %LOCALRESOLVEDUALCAMERARIG Snapshot the Active Rig and its transform.
        %   An empty rigBinding means the data are imported without
        %   coregistration (no Active Rig, or no active transform). Several
        %   Active Rigs or an invalid active resource are errors.

        rigBinding = [];
        resource = [];
        resourcePath = '';
        payload = [];
        reasonText = '';
        skipCode = '';
        if UMITRigStore.validateLifecycle().activeCount == 0
            skipCode = 'noActiveRig';
            reasonText = 'no Rig is Active in the rig library';
            return
        end
        activeStore = UMITRigStore.getActiveRig();
        activeInfo = activeStore.getRigInfo();
        % getActiveCameraCoregistration returns empty when no pointer is
        % configured and throws for inconsistent pointers or missing managed
        % files; those errors must not be converted into a skip.
        try
            resource = activeStore.getActiveCameraCoregistration();
            if ~isempty(resource)
                resourcePath = activeStore.resolveResourcePath(resource.uuid);
                payload = localLoadCoregistrationResource(resourcePath);
            end
        catch ME
            coregError = MException( ...
                'Umitoolbox:ImagesClassification:InvalidActiveCameraCoregistration', ...
                ['The active camera-coregistration resource for Rig "%s" (%s) ' ...
                 'could not be validated: %s'], activeInfo.rigID, activeInfo.uuid, ME.message);
            coregError = addCause(coregError, ME);
            throw(coregError)
        end
        if isempty(resource)
            skipCode = 'noActiveResource';
            reasonText = sprintf(['the Active Rig "%s" (%s) has no active ' ...
                'camera-coregistration resource'], activeInfo.displayName, activeInfo.rigID);
            return
        end
        rigBinding = UMITRigLedger.bindingFromRig(activeStore);
    end

    function localBindCoregisteredFolder(rigBinding)
        %LOCALBINDCOREGISTEREDFOLDER Bind the folder and commit the transform.

        ledgerState = UMITRigLedger.read(SaveFolder).state;
        if ismember(ledgerState, {'Unknown', 'NoRig'})
            UMITRigLedger.bind(SaveFolder, rigBinding, 'Origin', ledgerOrigin);
        else
            % The freshImport event reset every commitment: no conflict.
            rebindConflict = UMITRigLedger.rebind(SaveFolder, rigBinding, 'Origin', ledgerOrigin);
            if rebindConflict.hasConflict
                error('Umitoolbox:ImagesClassification:RigConflict', ...
                    'The re-import created a rig conflict in "%s":\n%s', ...
                    SaveFolder, rebindConflict.message);
            end
        end
        UMITRigLedger.commit(SaveFolder, 'coregistration', rigBinding.coregistration.hash, ...
            'Origin', ledgerOrigin, 'FirstUsedBy', ledgerOrigin);
        fprintf(['Folder rig ledger: bound to Rig "%s" (%s) and committed its camera ' ...
            'coregistration.\n'], rigBinding.sourceRig.displayName, rigBinding.sourceRig.rigID);
    end

    function localRecordUncoregisteredFolder(reasonText)
        %LOCALRECORDUNCOREGISTEREDFOLDER Ledger state for uncoregistered data.
        %   A new folder records "no rig"; a bound folder keeps its optics
        %   binding without a coregistration.

        ledgerState = UMITRigLedger.read(SaveFolder).state;
        switch ledgerState
            case 'Unknown'
                UMITRigLedger.writeNoRig(SaveFolder, 'Origin', ledgerOrigin);
                ledgerNote = 'an explicit "no rig"';
            case 'NoRig'
                ledgerNote = 'an explicit "no rig"';
            otherwise
                dropConflict = UMITRigLedger.dropCoregistrationBinding( ...
                    SaveFolder, 'Origin', ledgerOrigin);
                if dropConflict.hasConflict
                    error('Umitoolbox:ImagesClassification:RigConflict', ...
                        'The re-import created a rig conflict in "%s":\n%s', ...
                        SaveFolder, dropConflict.message);
                end
                ledgerNote = 'its existing optics binding, without a camera coregistration';
        end
        warning('Umitoolbox:ImagesClassification:ImportedWithoutCoregistration', ...
            ['Dual-camera data were imported into "%s" WITHOUT camera coregistration ' ...
            'because %s. The folder rig ledger records %s. To coregister, set an ' ...
            'Active Rig with an active camera-coregistration resource in the Rig ' ...
            'Manager, then re-import the raw data.'], SaveFolder, reasonText, ledgerNote);
    end

    function recorded = localRecordCoregistrationProvenance(rigBinding, resource, resourcePath, payload)
        %LOCALRECORDCOREGISTRATIONPROVENANCE Record the application in
        %DataParams.mat. Best-effort: a failure here must not undo an
        %otherwise-successful coregistration.

        recorded = false;
        try
            if ~isfile(fullfile(SaveFolder, 'DataParams.mat'))
                createDataParams(SaveFolder);
            end

            DataParams = loadDataParams(SaveFolder);
            cc = DataParams.cameraCoregistration;

            cc.isCoregistered = true;
            cc.isReviewed = false;
            cc.tform = payload.tform;
            cc.resourceUUID = resource.uuid;
            cc.rigID = rigBinding.sourceRig.rigID;
            cc.method = 'ImagesClassification';
            cc.sourceFile = resourcePath;
            if isfield(payload.tformInfo, 'SavedOn')
                cc.sourceFileTimestamp = char(payload.tformInfo.SavedOn);
            else
                cc.sourceFileTimestamp = '';
            end
            cc.createdOn = char(datetime('now'));
            cc.source = SaveFolder;
            cc.appliedOn = char(datetime('now'));
            cc.appliedBy = 'ImagesClassification';
            cc.confirmationMode = 'automatic-import-application';
            cc.notes = ['Dual-camera coregistration automatically applied during ' ...
                'raw data import using the rig''s active camera coregistration.'];

            % getTformQCMetrics normalises the premultiply/postmultiply
            % convention, so the recorded metrics do not depend on which
            % transform class the saved resource happens to hold.
            tformQC = getTformQCMetrics(payload.tform);
            cc.transformType = tformQC.transformType;
            cc.qcMetrics.translationXY_px = tformQC.translationXY;
            cc.qcMetrics.rotationDeg = tformQC.rotationDeg;
            cc.qcMetrics.scaleXY = tformQC.scaleXY;
            cc.qcMetrics.determinant = tformQC.determinant;

            DataParams.cameraCoregistration = cc;
            DataParams.lastModified = datetime('now');

            saveDataParams(SaveFolder, DataParams);
            recorded = true;
        catch provenanceError
            warning('Umitoolbox:ImagesClassification:cameraCoregistrationDataParamsFailed', ...
                ['Coregistration succeeded, but DataParams.cameraCoregistration could not ' ...
                'be recorded: %s'], provenanceError.message);
        end
    end

    function result = localSummarizeLedger(result)
        %LOCALSUMMARIZELEDGER Copy the folder ledger state into the result.

        ledgerInfo = UMITRigLedger.read(SaveFolder);
        result.ledgerState = ledgerInfo.state;
        if ~isempty(ledgerInfo.binding)
            result.rigUUID = ledgerInfo.binding.sourceRig.uuid;
            result.rigID = ledgerInfo.binding.sourceRig.rigID;
            result.rigDisplayName = ledgerInfo.binding.sourceRig.displayName;
            result.opticsHash = ledgerInfo.binding.opticsHash;
            result.coregistrationHash = ledgerInfo.binding.coregistrationHash;
        end
    end

    function payload = localLoadCoregistrationResource(resourcePath)
        %LOCALLOADCOREGISTRATIONRESOURCE Validate the managed MAT payload.

        payload = load(resourcePath, '-mat');
        if ~isfield(payload, 'tform') || isempty(payload.tform)
            error('Umitoolbox:ImagesClassification:InvalidCoregistrationPayload', ...
                'Managed file must contain a nonempty variable named "tform".');
        end
        payload.tform = localNormalizeCoregistrationTransform(payload.tform);
        if ~isfield(payload, 'tformInfo')
            payload.tformInfo = struct();
        elseif ~isstruct(payload.tformInfo) || ~isscalar(payload.tformInfo)
            error('Umitoolbox:ImagesClassification:InvalidCoregistrationPayload', ...
                'Variable "tformInfo" must be a scalar struct when present.');
        end

        T = localGetTransformMatrix(payload.tform);
        if ~isnumeric(T) || ~isequal(size(T), [3 3]) || any(~isfinite(T), 'all')
            error('Umitoolbox:ImagesClassification:InvalidCoregistrationPayload', ...
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
                error('Umitoolbox:ImagesClassification:InvalidCoregistrationPayload', ...
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
            error('Umitoolbox:ImagesClassification:InvalidCoregistrationPayload', ...
                ['Variable "tform" must be affine2d, a supported modern 2-D ' ...
                 'affine transform, or a finite numeric 3-by-3 matrix.']);
        end

        matrix = double(tform.A);
        if ~isequal(size(matrix), [3 3]) || any(~isfinite(matrix), 'all')
            error('Umitoolbox:ImagesClassification:InvalidCoregistrationPayload', ...
                'The modern camera-coregistration transform must have a finite 3-by-3 A matrix.');
        end
        tform = affine2d(matrix.');
    end

end

function outFile = iClassifyRawBinaries(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_SubROI)
%ICLASSIFYRAWBINARIES Split interlaced raw binaries into channel .dat files.
%   Rig-independent classification. Writes one headered .dat file per
%   output channel and AcqInfos.mat in SaveFolder, and returns the
%   SaveFolder-relative .dat manifest. The caller has already reset a
%   previously imported SaveFolder (prepareFreshImport).

if ~isfolder(SaveFolder)
    mkdir(SaveFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outFile = {};
AcqInfoStream = ReadInfoFile(DataFolder);
rawFrameRateHz = AcqInfoStream.FrameRateHz;

if ~strcmp(DataFolder(end), filesep)
    DataFolder = strcat(DataFolder, filesep);
end
if ~strcmp(SaveFolder(end), filesep)
    SaveFolder = strcat(SaveFolder, filesep);
end

if ~isfield(AcqInfoStream, 'Camera_Model')
    warning('Umitoolbox:ImagesClassification:CameraModelMissing', ...
        ['info.txt has no "Camera Model" entry; assuming "D1024". The camera ' ...
         'model selects how dropped frames are detected and filled.']);
    AcqInfoStream.Camera_Model = 'D1024';
end
if ~exist(SaveFolder, 'dir')
    mkdir(SaveFolder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Format and Header Information:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hWima = 5;
imgFilesList = dir([DataFolder 'img*.bin']);
if isempty(imgFilesList)
    error('No image binary files found! Classification aborted.')
end

% Validate camera-specific image file sequences. In dual-camera recordings,
% both img_*.bin and imgCam2_*.bin sequences are validated independently.
imgFilesListCam1 = dir([DataFolder 'img_*.bin']);
imgFilesListCam2 = dir([DataFolder 'imgCam2_*.bin']);

if isempty(imgFilesListCam1)
    error('Image binary files (img_xxxxx.bin) missing! Classification aborted.')
end
if AcqInfoStream.MultiCam && isempty(imgFilesListCam2)
    error('Camera #2 image binary files (imgCam2_xxxxx.bin) missing! Classification aborted.')
end

imgFileSet = {imgFilesListCam1};
if AcqInfoStream.MultiCam
    imgFileSet{end+1} = imgFilesListCam2;
end

for ii = 1:length(imgFileSet)
    thisImgFilesList = imgFileSet{ii};
    imgFileNames = sort({thisImgFilesList.name})';
    imgFileIndx = str2double(regexprep(imgFileNames, '^(img_|imgCam2_)|\.bin$', ''));
    if any(isnan(imgFileIndx)) || imgFileIndx(1) ~= 0 || any(diff(imgFileIndx) ~= 1)
        error('Image binary files missing! Classification aborted.')
    end
end

headerFiles = imgFilesListCam1;
header = memmapfile([DataFolder headerFiles(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
nx = double(header.Data.header(2));
ny = double(header.Data.header(3));
frameFormat = {'uint64', 3, 'framej'; 'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx * ny * 2 + 3 * 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_SubROI
    fprintf('Redefining Region Of Interest post-process:\n');
    ButtonName = questdlg('Would you like to use a pre-defined ROI?', ...
        'ROI', ...
        'Pre-defined', 'Draw', 'Cancel', 'Draw');
    switch ButtonName
        case 'Pre-defined'
            [filename, pathname] = uigetfile('*.mat', 'Select ROI file');
            if isequal(filename,0) || isequal(pathname,0)
                disp('User pressed cancel')
                Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
            else
                load([pathname filesep filename]); %#ok<LOAD>
            end
        case 'Draw'
            dat = memmapfile([DataFolder headerFiles(1).name], ...
                'Offset', hWima * 4 + 5 * SizeImage, ...
                'Format', frameFormat, 'repeat', 1);
            dat = rot90(fliplr(dat.Data.imgj));
            fig = figure('Name', 'Draw ROI', 'CloseRequestFcn', @closeFig);
            imagesc(dat);
            axis image;
            drawrectangle('Deletable', false, 'Tag', 'myRectangle');
            title('Close figure to confirm')
            waitfor(fig)
        case 'Cancel'
            disp('User pressed cancel')
            Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
    end

    LimX = [round(Pos(1)) round(Pos(1)+Pos(3)) - 1];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4)) - 1];
    save([SaveFolder 'ROI.mat'], 'Pos');
    clear Pos
else
    LimX = [1 ImRes_XY(1)];
    LimY = [1 ImRes_XY(2)];
end

Rx = round((LimX(2) - LimX(1) + 1) / BinningSpatial);
Ry = round((LimY(2) - LimY(1) + 1) / BinningSpatial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How many colors and in which order?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Sorting images per channels.\n');
if AcqInfoStream.MultiCam
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'CamIdx', {}, 'FrameIdx', {}, 'Exposure', {}, 'Tag', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        eval(['Colors(' int2str(indC) ').CamIdx = AcqInfoStream.Illumination' int2str(indC) '.CamIdx;']);
        eval(['Colors(' int2str(indC) ').FrameIdx = AcqInfoStream.Illumination' int2str(indC) '.FrameIdx;']);
        Colors(indC).Exposure = localResolveExposure(AcqInfoStream, Colors(indC).Color);
        Colors(indC).Tag = localGetChannelTag(Colors(indC).Color);
    end

    % Repeated output tags are allowed within one camera sequence but not
    % across cameras because each output tag maps to one .dat file name.
    TagList = {Colors.Tag};
    UniqueTags = unique(TagList, 'stable');
    for indT = 1:length(UniqueTags)
        idxTag = find(strcmp(TagList, UniqueTags{indT}));
        if length(unique([Colors(idxTag).CamIdx])) > 1
            error(['Repeated illumination channel "' UniqueTags{indT} ...
                '" was detected across cameras. Repeated illuminations ' ...
                'are supported only within a single camera sequence.']);
        end
    end

    fprintf('Camera #1.\n');
    imgFilesList = dir([DataFolder 'img_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 1, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    channelInfo1 = ChannelsSort(imgFilesList, Colors(idx));
    fprintf('Camera #1 Done.\n');

    fprintf('Camera #2.\n');
    imgFilesList = dir([DataFolder 'imgCam2_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 2, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    channelInfo2 = ChannelsSort(imgFilesList, Colors(idx));
    channelInfoByCam = {channelInfo1, channelInfo2};
    channelInfo = [channelInfo1, channelInfo2];
    fprintf('Camera #2 Done.\n');
else
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'Exposure', {}, 'Tag', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        Colors(indC).Exposure = localResolveExposure(AcqInfoStream, Colors(indC).Color);
        Colors(indC).Tag = localGetChannelTag(Colors(indC).Color);
    end

    imgFilesList = dir([DataFolder 'img_*.bin']);
    channelInfo = ChannelsSort(imgFilesList, Colors);
    channelInfoByCam = {channelInfo};
    fprintf('Done.\n');
end

% Final output validation. The output contract is validated independently
% for each camera sequence.
for indCam = 1:length(channelInfoByCam)
    datLenList = [channelInfoByCam{indCam}.datLength];
    repeatCountList = [channelInfoByCam{indCam}.RepeatCount];
    baseLenList = datLenList ./ repeatCountList;
    assert(~(any(mod(baseLenList, 1) ~= 0) || ~isscalar(unique(baseLenList))), ...
        ['Channel classification failed final length validation for camera #' ...
        int2str(indCam) '. Expected datLength/RepeatCount to be equal ' ...
        'for all output channels from the same camera.']);
end

% AcqInfos.mat describes the raw acquisition: Width, Height, and
% FrameRateHz keep their info.txt values (.dat header Phase 7b). The
% imported data (frame size, length, rate, class) is described by each
% .dat header. Only the import settings and the channel manifest are added.
AcqInfoStream.BinningSpatial = BinningSpatial;
AcqInfoStream.BinningTemp = BinningTemp;

if isfield(AcqInfoStream, 'ImportedChannels')
    AcqInfoStream = rmfield(AcqInfoStream, 'ImportedChannels');
end

for ind = 1:numel(channelInfo)
    importedInfo = struct();
    importedInfo.DatFile = channelInfo(ind).datFile;
    importedInfo.Tag = channelInfo(ind).Tag;
    importedInfo.Color = channelInfo(ind).Color;
    importedInfo.Length = channelInfo(ind).datLength;
    importedInfo.FrameRateHz = channelInfo(ind).FrameRateHz;
    importedInfo.ExposureMsec = channelInfo(ind).ExposureMsec;
    if isempty(channelInfo(ind).CamIdx)
        importedInfo.CamIdx = 1;
    else
        importedInfo.CamIdx = channelInfo(ind).CamIdx;
    end
    AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, importedInfo, 'Overwrite', true);
end

saveMatAtomic([SaveFolder 'AcqInfos.mat'], 'AcqInfoStream', AcqInfoStream);
outFile = [outFile, {channelInfo.datFile}];

    function channelInfoOut = ChannelsSort(fList, colors)
        %CHANNELSSORT Split interleaved binary frames into grouped .dat files.

        subNbColors = size(colors, 2);
        OutputGroups = struct('Tag', {}, 'Color', {}, 'SeqIdx', {}, ...
            'RepeatCount', {}, 'Exposure', {}, 'dTag', {}, 'CamIdx', {});

        for indSeq = 1:subNbColors
            idxGroup = [];
            if ~isempty(OutputGroups)
                idxGroup = find(strcmp({OutputGroups.Tag}, colors(indSeq).Tag), 1, 'first');
            end

            if isempty(idxGroup)
                idxGroup = length(OutputGroups) + 1;
                OutputGroups(idxGroup).Tag = colors(indSeq).Tag;
                OutputGroups(idxGroup).Color = colors(indSeq).Color;
                OutputGroups(idxGroup).SeqIdx = indSeq;
                OutputGroups(idxGroup).RepeatCount = 1;
                OutputGroups(idxGroup).Exposure = colors(indSeq).Exposure;
                OutputGroups(idxGroup).dTag = [colors(indSeq).Tag '.dat'];
                if isfield(colors, 'CamIdx')
                    OutputGroups(idxGroup).CamIdx = colors(indSeq).CamIdx;
                else
                    OutputGroups(idxGroup).CamIdx = [];
                end
            else
                if colors(indSeq).Exposure ~= OutputGroups(idxGroup).Exposure
                    error(['Repeated illumination channel "' colors(indSeq).Tag ...
                        '" has inconsistent exposure values.']);
                end
                OutputGroups(idxGroup).SeqIdx = [OutputGroups(idxGroup).SeqIdx indSeq];
                OutputGroups(idxGroup).RepeatCount = length(OutputGroups(idxGroup).SeqIdx);
            end
        end

        channelInfoOut = repmat(struct( ...
            'datFile', '', ...
            'datSize', [Ry, Rx], ...
            'datLength', 0, ...
            'Tag', '', ...
            'Color', '', ...
            'FrameRateHz', [], ...
            'ExposureMsec', [], ...
            'CamIdx', [], ...
            'RepeatCount', []), ...
            1, length(OutputGroups));
        % Headered channel files grow as chunks are sorted: each is created
        % (growable) at its first chunk and finalized after the length checks.
        slabs = cell(1, length(OutputGroups));
        % The map is a handle object, so the cleanup sees files created
        % later; on an error the started files stay marked incomplete.
        openSlabs = containers.Map('KeyType', 'double', 'ValueType', 'any');
        cleanupSlabs = onCleanup(@() iCloseSlabs(openSlabs));

        for indG = 1:length(OutputGroups)
            dTag = OutputGroups(indG).dTag;
            channelInfoOut(indG).datFile = dTag;
            channelInfoOut(indG).datSize = [Ry, Rx];
            channelInfoOut(indG).datLength = 0;
            channelInfoOut(indG).Tag = OutputGroups(indG).Tag;
            channelInfoOut(indG).Color = OutputGroups(indG).Color;
            channelInfoOut(indG).FrameRateHz = rawFrameRateHz * ...
                OutputGroups(indG).RepeatCount / (subNbColors * BinningTemp);
            channelInfoOut(indG).ExposureMsec = OutputGroups(indG).Exposure;
            channelInfoOut(indG).CamIdx = OutputGroups(indG).CamIdx;
            channelInfoOut(indG).RepeatCount = OutputGroups(indG).RepeatCount;
            if exist([SaveFolder dTag], 'file')
                delete([SaveFolder dTag]);
            end
        end

        ImBinBuffer = cell(1, length(OutputGroups));
        ImWriteBuffer = cell(1, length(OutputGroups));
        datLengthList = zeros(1, length(OutputGroups));
        for indG = 1:length(OutputGroups)
            ImBinBuffer{indG} = [];
            ImWriteBuffer{indG} = [];
        end

        oIm = [];
        Cnt = 0;
        for indF = 1:size(fList,1)
            fprintf('Sorting %s.', fList(indF).name);
            data = memmapfile([DataFolder fList(indF).name], ...
                'Offset', hWima * 4, 'Format', frameFormat, 'repeat', inf);
            data = data.Data;
            hData = reshape([data.framej], 3, []);
            iData = reshape([data.imgj], ImRes_XY(1), ImRes_XY(2), []);
            iData = permute(iData, [2 1 3]);
            clear data;

            if contains(AcqInfoStream.Camera_Model, {'D1024', 'D1312'})
                SkipNFirst = sum(hData(1,:) == 0);
                MissingOffset = cumsum(hData(2,:));
                hData(1,:) = hData(1,:) + MissingOffset - hData(1,1) + 1;
                goodFrames = find(accumarray(hData(1, (SkipNFirst+1):end)',1) == 1)';
                ConseqFromLeft = [1 diff(goodFrames,1,2) == 1];
                ConseqFromRight = fliplr([true diff(fliplr(goodFrames),1,2) == -1]);
                goodFrames = goodFrames(ConseqFromLeft | ConseqFromRight);
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1), 'uint16');
                [~, goodFramesInIData] = ismember(goodFrames(1,:), hData(1,:));
                Images(:,:,goodFrames) = iData(:,:,goodFramesInIData);
                iData = Images;
            elseif contains(AcqInfoStream.Camera_Model, 'BFLY')
                iNbF = hData(2,1) - Cnt;
                if (hData(2,end) - hData(2,1)) > 0
                    fprintf('\t WARNING: %d missing frames.', (hData(2,end) - hData(2,1)));
                    Cnt = Cnt + (hData(2,end) - hData(2,1));
                end
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1 + iNbF), 'uint16');
                Images(:,:,(hData(1,:) - hData(1,1) + 1 + iNbF)) = iData;
                iData = Images;
                clear Images;
            elseif any(hData(2,:))
                fprintf('\t WARNING: %d missing frames.', sum(hData(2,:)));
                hData(1,:) = 1:size(hData,2);
                hData(1,:) = hData(1,:) + cumsum(hData(2,:));
                if hData(2,1) >= 1
                    iNbF = hData(2,1);
                else
                    iNbF = 0;
                end
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1 + iNbF), 'uint16');
                Images(:,:,(hData(1,:) - hData(1,1) + 1 + iNbF)) = iData;
                iData = Images;
                clear Images;
            end

            iData = cat(3, oIm, iData);
            overflow = mod(size(iData,3), subNbColors);
            if overflow > 0
                oIm = iData(:,:, size(iData,3) - (overflow:-1:1) + 1);
            else
                oIm = [];
            end
            Images = iData(:,:,1:(size(iData,3) - overflow));
            clear iData hData overflow;

            if isempty(Images)
                fprintf('\n');
                continue
            end
            Images = reshape(Images, ImRes_XY(2), ImRes_XY(1), subNbColors, []);

            for indG = 1:length(OutputGroups)
                seqIdx = OutputGroups(indG).SeqIdx;
                Ims = Images(:, :, seqIdx, :);
                Ims = reshape(Ims, ImRes_XY(2), ImRes_XY(1), []);

                if b_SubROI
                    Ims = Ims(round(LimY(1)):round(LimY(2)), round(LimX(1)):round(LimX(2)), :);
                end

                if ~isempty(ImBinBuffer{indG})
                    Ims = cat(3, ImBinBuffer{indG}, Ims);
                end

                nFramesAll = size(Ims,3);
                nFullRaw = floor(nFramesAll / BinningTemp) * BinningTemp;

                if nFullRaw < nFramesAll
                    ImBinBuffer{indG} = Ims(:,:,nFullRaw+1:end);
                else
                    ImBinBuffer{indG} = [];
                end

                if nFullRaw == 0
                    continue
                end

                imgSum = squeeze(sum(sum(Ims,1),2));
                imgSum = imgSum(:)';
                if any(imgSum == 0)
                    idxGood = find(imgSum > 1);
                    if length(idxGood) > 1
                        nY = size(Ims,1);
                        nX = size(Ims,2);
                        Ims = interp1(idxGood, single(reshape(Ims(:,:,idxGood),[], length(idxGood)))', ...
                            1:size(Ims,3), 'linear', 'extrap');
                        Ims = reshape(Ims', nY, nX, []);
                    elseif length(idxGood) == 1
                        Ims = repmat(single(Ims(:,:,idxGood)), 1, 1, size(Ims,3));
                    end
                end

                Ims = Ims(:,:,1:nFullRaw);

                if BinningTemp > 1
                    Ims = imresize3(Ims, [size(Ims,1), size(Ims,2), ...
                        size(Ims,3) / BinningTemp], 'linear');
                end

                if BinningSpatial > 1
                    Ims = imresize(Ims, 1 / BinningSpatial);
                end

                if isempty(ImWriteBuffer{indG})
                    ImWriteBuffer{indG} = single(Ims);
                else
                    ImWriteBuffer{indG} = cat(3, ImWriteBuffer{indG}, single(Ims));
                end
            end
            clear Images;

            baseLenList = zeros(1, length(OutputGroups));
            for indG = 1:length(OutputGroups)
                if isempty(ImWriteBuffer{indG})
                    nBufferedFrames = 0;
                else
                    nBufferedFrames = size(ImWriteBuffer{indG},3);
                end
                baseLenList(indG) = floor((datLengthList(indG) + nBufferedFrames) / ...
                    OutputGroups(indG).RepeatCount);
            end
            targetBaseLen = min(baseLenList);
            currentBaseLenList = datLengthList ./ [OutputGroups.RepeatCount];
            assert(~(any(mod(currentBaseLenList, 1) ~= 0) || ...
                    ~isscalar(unique(currentBaseLenList))), ...
                'Internal length bookkeeping failed during channel sorting.');
            currentBaseLen = currentBaseLenList(1);

            if targetBaseLen > currentBaseLen
                for indG = 1:length(OutputGroups)
                    nFramesToWrite = (targetBaseLen - currentBaseLen) * ...
                        OutputGroups(indG).RepeatCount;
                    if nFramesToWrite > 0
                        block = ImWriteBuffer{indG}(:,:,1:nFramesToWrite);
                        if isempty(slabs{indG})
                            slabs{indG} = createChannelFile(indG, size(block, 1), size(block, 2));
                            openSlabs(indG) = slabs{indG};
                        end
                        slabs{indG} = spatialSlabIO('append', slabs{indG}, block);
                        clear block

                        if nFramesToWrite < size(ImWriteBuffer{indG},3)
                            ImWriteBuffer{indG} = ImWriteBuffer{indG}(:,:,nFramesToWrite+1:end);
                        else
                            ImWriteBuffer{indG} = [];
                        end
                        datLengthList(indG) = datLengthList(indG) + nFramesToWrite;
                        channelInfoOut(indG).datLength = datLengthList(indG);
                    end
                end
            end
            fprintf('\n');
        end

        for indG = 1:length(OutputGroups)
            if ~isempty(ImBinBuffer{indG})
                fprintf('Discarding %d unbinned overflow frame(s) from %s.\n', ...
                    size(ImBinBuffer{indG},3), OutputGroups(indG).Tag);
            end
            if ~isempty(ImWriteBuffer{indG})
                fprintf('Discarding %d ratio-overflow frame(s) from %s.\n', ...
                    size(ImWriteBuffer{indG},3), OutputGroups(indG).Tag);
            end
        end

        finalBaseLenList = datLengthList ./ [OutputGroups.RepeatCount];
        assert(~(any(mod(finalBaseLenList, 1) ~= 0) || ...
                ~isscalar(unique(finalBaseLenList))), ...
            ['Channel classification failed final length validation. ' ...
            'Expected datLength/RepeatCount to be equal for all output channels.']);

        for indG = 1:length(OutputGroups)
            if isempty(slabs{indG})
                error('Umitoolbox:ImagesClassification:NoFramesWritten', ...
                    'No frames were sorted into "%s".', OutputGroups(indG).dTag);
            end
        end
        for indG = 1:length(OutputGroups)
            spatialSlabIO('finalize', slabs{indG});
        end
        clear cleanupSlabs

        function h = createChannelFile(indGroup, nY, nX)
            % Growable Y-X-T single file carrying the channel's own rate,
            % exposure, and tag.
            hdr = datHeaderFromInfo(struct('dataClass', 'single', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [nY, nX, 1], ...
                'frameRateHz', channelInfoOut(indGroup).FrameRateHz, ...
                'exposureMsec', channelInfoOut(indGroup).ExposureMsec), ...
                OutputGroups(indGroup).Tag);
            h = spatialSlabIO('create', [SaveFolder OutputGroups(indGroup).dTag], hdr, ...
                'Growable', true);
        end
    end

    function closeFig(src,~)
        % CloseRequest function for drawing ROI.
        rectH = findobj(src, 'Tag', 'myRectangle');
        Pos = rectH.Position;
        delete(src)
    end
end

function iCloseSlabs(openSlabs)
%ICLOSESLABS Close channel files that are still open (closing twice is harmless).
handles = values(openSlabs);
for k = 1:numel(handles)
    spatialSlabIO('close', handles{k});
end
end

function exposure = localResolveExposure(AcqInfoStream, colorName)
%LOCALRESOLVEEXPOSURE Resolve per-channel exposure metadata.

if contains(colorName, {'red','amber','green'}, 'IgnoreCase', true)
    exposure = AcqInfoStream.ExposureMsec;
elseif contains(colorName, {'speckle'}, 'IgnoreCase', true)
    if ~isfield(AcqInfoStream, 'ExposureSpeckleMsec')
        exposure = AcqInfoStream.ExposureMsec;
    else
        exposure = AcqInfoStream.ExposureSpeckleMsec;
    end
else
    if ~isfield(AcqInfoStream, 'ExposureFluoMsec')
        exposure = AcqInfoStream.ExposureMsec;
    else
        exposure = AcqInfoStream.ExposureFluoMsec;
    end
end
end

function tag = localGetChannelTag(colorName)
%LOCALGETCHANNELTAG Convert illumination color metadata to output tag.

if contains(colorName, 'red', 'IgnoreCase', true)
    tag = 'red';
elseif contains(colorName, 'green', 'IgnoreCase', true)
    tag = 'green';
elseif contains(colorName, 'amber', 'IgnoreCase', true)
    tag = 'yellow';
elseif contains(colorName, 'speckle', 'IgnoreCase', true)
    tag = 'speckle';
else
    waveTag = regexp(colorName, '[0-9]{3}', 'match');
    if ~isempty(waveTag)
        tag = ['fluo_' waveTag{:}];
    else
        tag = 'fluo';
    end
end
end

