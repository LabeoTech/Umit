function classifiedFiles = run_ImagesClassification(RawFolder, SaveFolder, varargin)
%RUN_IMAGESCLASSIFICATION Classify raw interlaced IOS binaries into channels.
%
%   outFile = run_ImagesClassification(RawFolder, SaveFolder)
%   outFile = run_ImagesClassification(RawFolder, SaveFolder, 'Name', Value, ...)
%   info    = run_ImagesClassification('pipelineInfo')
%
%   This wrapper calls ImagesClassification to split raw interlaced binary
%   data into separate channel .dat files and a shared AcqInfos.mat file.
%   For dual-camera acquisitions, it also attempts to apply the saved
%   camera coregistration transform when available.
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
%                        Allowed values:
%                            'ERASE'
%                            'GENBACKUP'
%                        Default: ''
%
%                        Notes:
%                        - If left empty, the backup handling is resolved
%                          interactively by the called function.
%                        - 'ERASE' deletes managed existing files from
%                          SaveFolder before import.
%                        - 'GENBACKUP' creates a timestamped .zip backup of
%                          managed existing files before import.
%
%   Output:
%       outFile        - File manifest of outputs saved in SaveFolder.
%
%   Notes:
%       - This function does not expose SubROI selection. It always calls
%         ImagesClassification with b_SubROI = 0.
%       - The returned file manifest uses full paths.
%       - Dual-camera coregistration is attempted only when:
%           1) AcqInfos.mat indicates MultiCam = true
%           2) a saved coregistration tform file is available
%
%   Examples:
%       outFile = run_ImagesClassification(rawFolder, saveFolder);
%
%       outFile = run_ImagesClassification( ...
%           rawFolder, saveFolder, ...
%           'BinningSpatial', 2, ...
%           'BinningTemp', 4);
%
%       outFile = run_ImagesClassification( ...
%           rawFolder, saveFolder, ...
%           'backupOpts', 'GENBACKUP');

default_Output = {'fluo_475.dat', 'fluo_567.dat', 'fluo.dat', ...
    'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat', 'AcqInfos.mat'};
allowedBinning = [1:8];

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

classifiedFiles = ImagesClassification( ...
    RawFolder, ...
    SaveFolder, ...
    BinningSpatial, ...
    BinningTemp, ...
    0, ...
    'backupOpts', backupOpts);

if ~iscell(classifiedFiles)
    classifiedFiles = {};
end

% For Dual-Camera Imaging systems, apply the coregistration using the tform
% file created in DataViewer's OiS Dual Cam Coregistration utility.
acqInfoPath = fullfile(SaveFolder, 'AcqInfos.mat');
if isfile(acqInfoPath)
    info = load(acqInfoPath);
    if isfield(info, 'AcqInfoStream') && isstruct(info.AcqInfoStream) && ...
            isfield(info.AcqInfoStream, 'MultiCam') && info.AcqInfoStream.MultiCam

        disp('Dual Camera data found!')

        % Ensure a rig exists for this machine. Import is the entry point of
        % data into the toolbox, so a rigless project/DataViewer instance
        % should never arise from normal use -- a default rig is silently
        % auto-created here if none exists yet (see UMITRigStore.
        % getOrCreateDefaultRig). This never blocks import: which rig ends
        % up owning centralized resources is resolved later, independently
        % of this call.
        defaultRigID = '';
        defaultRigStore = [];
        try
            defaultRigStore = UMITRigStore.getOrCreateDefaultRig();
            defaultRigID = defaultRigStore.getRigInfo().rigID;
        catch ME
            warning('Umitoolbox:run_ImagesClassification:defaultRigFailed', ...
                ['Could not resolve/create a default rig for this machine. ' ...
                'Import will continue without one. %s'], ME.message);
        end

        % Resolve the active camera coregistration transform. Prefer the
        % rig's managed resource (kept current by DataViewer_Coreg2Cams's
        % "Save Calibration" action) so imported data carries a real
        % resourceUUID; fall back to the legacy flat tform file for rigs
        % that have never saved a camera coregistration through the
        % managed-resource path yet.
        tformFile = '';
        resourceUUID = '';
        if ~isempty(defaultRigStore)
            try
                activeResource = defaultRigStore.getActiveCameraCoregistration();
                if ~isempty(activeResource)
                    tformFile = defaultRigStore.resolveResourcePath(activeResource.uuid);
                    resourceUUID = activeResource.uuid;
                end
            catch ME
                warning('Umitoolbox:run_ImagesClassification:activeCoregistrationLookupFailed', ...
                    ['Could not resolve the rig''s active camera coregistration resource. ' ...
                    'Falling back to the legacy tform file. %s'], ME.message);
                tformFile = '';
                resourceUUID = '';
            end
        end

        if isempty(tformFile)
            LabeoFolder = getUmitFolder('tformFiles');
            tformFile = fullfile(LabeoFolder, 'coreg2cam_tform.mat');
        end

        if isfile(tformFile)
            tf = load(tformFile);
            disp('Applying coregistration to data from camera #2...');
            [status, warnmsg] = applyTform2Cams(SaveFolder, tf.tform, tf.tformInfo);
            if ~status
                warning(['Coregistration failed! Data import will resume without coregistration. ', warnmsg]);
            else
                disp('Done!')

                % Record this application in the imported SaveFolder's
                % DataParams.mat. Only DataViewer_Coreg2Cams's own "Save
                % Calibration" action recorded DataParams.cameraCoregistration
                % before now -- not the many downstream SaveFolders that
                % actually receive the transform via import, which is the
                % common case. Best-effort: a failure here must not undo an
                % otherwise-successful coregistration.
                try
                    if ~isfile(fullfile(SaveFolder, 'DataParams.mat'))
                        createDataParams(SaveFolder);
                    end

                    DataParams = loadDataParams(SaveFolder);
                    cc = DataParams.cameraCoregistration;

                    cc.isCoregistered = true;
                    cc.isReviewed = false;
                    cc.tform = tf.tform;
                    cc.resourceUUID = resourceUUID;
                    cc.rigID = defaultRigID;
                    cc.transformType = 'similarity';
                    cc.method = 'run_ImagesClassification';
                    cc.sourceFile = tformFile;
                    if isfield(tf.tformInfo, 'SavedOn')
                        cc.sourceFileTimestamp = char(tf.tformInfo.SavedOn);
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

                    T = localGetTransformMatrix(tf.tform);
                    A = T(1:2, 1:2);
                    cc.qcMetrics.translationXY_px = [T(3,1), T(3,2)];
                    cc.qcMetrics.rotationDeg = atan2d(A(1,2), A(1,1));
                    cc.qcMetrics.scaleXY = [hypot(A(1,1), A(1,2)), hypot(A(2,1), A(2,2))];
                    cc.qcMetrics.determinant = det(A);

                    DataParams.cameraCoregistration = cc;
                    DataParams.lastModified = datetime('now');

                    validateDataParams(DataParams);
                    save(fullfile(SaveFolder, 'DataParams.mat'), 'DataParams', '-mat');
                catch ME
                    warning('Umitoolbox:run_ImagesClassification:cameraCoregistrationDataParamsFailed', ...
                        ['Coregistration succeeded, but DataParams.cameraCoregistration could not ' ...
                        'be recorded: %s'], ME.message);
                end
            end
        else
            warning(['Coregistration TFORM file not found! Data from camera #2 will not be ' ...
                'coregistered with camera #1. Data import will resume without coregistration.'])
        end
    end
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

% Return full saved paths, consistent with wrapper-level file manifest style.
classifiedFiles = unique(cellfun(@(x) fullfile(SaveFolder, x), classifiedFiles, 'UniformOutput', false), 'stable');

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Classify raw interlaced IOS binary data into channel .dat files.');
        info.version = '1.0.0';

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
            ['Backup handling option passed to ImagesClassification. Use '''', ''ERASE'', ' ...
             '''GENBACKUP'', or a custom zip base name.'], ...
            'kind', 'parameter', 'default', 'ERASE','allowed',{'ERASE','GENBACKUP'}, 'callType', 'namevalue');

        info = PipelineManager.addOutput(info, 'classifiedFiles', 'ImageTimeSeries', 'file', ...
            'Generated file manifest saved in SaveFolder.', ...
            default_Output, 1, 'isData', true, 'saveFileName', '');
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
end