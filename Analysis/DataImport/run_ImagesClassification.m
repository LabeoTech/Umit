function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
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
    outFile = localPipelineInfo();
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

outFile = ImagesClassification( ...
    RawFolder, ...
    SaveFolder, ...
    BinningSpatial, ...
    BinningTemp, ...
    0, ...
    'backupOpts', backupOpts);

if ~iscell(outFile)
    outFile = {};
end

% For Dual-Camera Imaging systems, apply the coregistration using the tform
% file created in DataViewer's OiS Dual Cam Coregistration utility.
acqInfoPath = fullfile(SaveFolder, 'AcqInfos.mat');
if isfile(acqInfoPath)
    info = load(acqInfoPath);
    if isfield(info, 'AcqInfoStream') && isstruct(info.AcqInfoStream) && ...
            isfield(info.AcqInfoStream, 'MultiCam') && info.AcqInfoStream.MultiCam

        disp('Dual Camera data found!')

        if ispc
            root = getenv('USERPROFILE');
        else
            root = getenv('HOME');
        end

        LabeoFolder = fullfile(root, 'Documents', 'LabeoTech', 'Config', 'umIT', 'tformFiles');
        tformFile = fullfile(LabeoFolder, 'coreg2cam_tform.mat');

        if isfile(tformFile)
            tf = load(tformFile);
            disp('Applying coregistration to data from camera #2...');
            [status, warnmsg] = applyTform2Cams(SaveFolder, tf.tform, tf.tformInfo);
            if ~status
                warning(['Coregistration failed! Data import will resume without coregistration. ', warnmsg]);
            else
                disp('Done!')
            end
        else
            warning(['Coregistration TFORM file not found! Data from camera #2 will not be ' ...
                'coregistered with camera #1. Data import will resume without coregistration.'])
        end
    end
end

% Return full saved paths, consistent with wrapper-level file manifest style.
outFile = unique(cellfun(@(x) fullfile(SaveFolder, x), outFile, 'UniformOutput', false), 'stable');

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

        info = PipelineManager.addOutput(info, 'outFile', 'ImageTimeSeries', 'file', ...
            'Generated file manifest saved in SaveFolder.', ...
            default_Output, 1, 'isData', true, 'saveFileName', '');
    end
end