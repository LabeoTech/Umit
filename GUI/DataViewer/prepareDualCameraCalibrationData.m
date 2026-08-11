function temporarySaveFolder = prepareDualCameraCalibrationData(RawFolder, options)
%PREPAREDUALCAMERACALIBRATIONDATA Classify raw dual-camera calibration data.
%   temporarySaveFolder = prepareDualCameraCalibrationData(RawFolder)
%   validates the acquisition metadata through ReadInfoFile, creates a
%   workflow-owned temporary SaveFolder, and runs the Rig-independent legacy
%   ImagesClassification importer with default classification settings. This
%   deliberately produces unregistered source data even when the default Rig
%   has an active transform. The caller owns the returned folder and must
%   remove it when calibration succeeds or is cancelled.
%
%   This helper never writes to RawFolder. If metadata validation or
%   classification fails, any temporary folder created by this call is
%   removed before the error is rethrown.

arguments
    RawFolder (1, 1) string {mustBeFolder}
    options.InfoReader (1, 1) function_handle = @ReadInfoFile
    options.ClassificationFcn (1, 1) function_handle = @classifyUnregisteredData
    options.BeforeClassificationFcn (1, 1) function_handle = @(~) []
    options.TempRoot (1, 1) string = string(tempdir)
end

RawFolder = string(char(RawFolder));

try
    acquisitionInfo = options.InfoReader(char(RawFolder));
catch ME
    metadataError = MException( ...
        'DataViewerCoreg2Cams:UnreadableAcquisitionMetadata', ...
        'Could not read acquisition metadata through ReadInfoFile: %s', ME.message);
    metadataError = addCause(metadataError, ME);
    throw(metadataError)
end

if ~isstruct(acquisitionInfo) || ~isscalar(acquisitionInfo) || ...
        ~isfield(acquisitionInfo, 'MultiCam') || ...
        ~isPositiveMultiCameraValue(acquisitionInfo.MultiCam)
    error('DataViewerCoreg2Cams:NotMultiCameraAcquisition', ...
        ['The selected raw acquisition could not be positively validated as ' ...
        'multi-camera through ReadInfoFile. Select a dual-camera acquisition.']);
end

tempRoot = char(options.TempRoot);
if ~isfolder(tempRoot)
    [created, message] = mkdir(tempRoot);
    if ~created
        error('DataViewerCoreg2Cams:TemporaryFolderFailed', ...
            'Could not create the temporary-folder root: %s', message);
    end
end

temporarySaveFolder = string(tempname(tempRoot));
[created, message] = mkdir(temporarySaveFolder);
if ~created
    error('DataViewerCoreg2Cams:TemporaryFolderFailed', ...
        'Could not create a temporary calibration SaveFolder: %s', message);
end

% The default classification adapter receives only the required folders;
% binning and other classification settings remain at their defaults.
try
    options.BeforeClassificationFcn(char(temporarySaveFolder));
    options.ClassificationFcn(char(RawFolder), char(temporarySaveFolder));
catch ME
    removeFolderIfPresent(temporarySaveFolder);
    rethrow(ME)
end
end

function tf = isPositiveMultiCameraValue(value)
%ISPOSITIVEMULTICAMERAVALUE Accept only an unambiguous scalar true value.

tf = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    ~ismissing(value) && isfinite(double(value)) && double(value) == 1;
end

function removeFolderIfPresent(folder)
%REMOVEFOLDERIFPRESENT Remove a workflow-owned temporary folder.

if isfolder(folder)
    rmdir(folder, 's');
end
end

function classifiedFiles = classifyUnregisteredData(rawFolder, saveFolder)
%CLASSIFYUNREGISTEREDDATA Run the Rig-independent importer without a tform.

classifiedFiles = ImagesClassification(rawFolder, saveFolder, 1, 1, 0);
end
