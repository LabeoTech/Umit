function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
% RUN_IMAGESCLASSIFICATION calls the function
% IMAGESCLASSIFICATION from the IOI library (LabeoTech).
% In brief, this function classifies the imaging channels from the raw data
% (img_00000.bin) into separate .dat files for each illumination color.

% Defaults:
default_Output = {'fluo_475.dat', 'fluo_567.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; % This is here only as a reference for PIPELINEMANAGER.m. The real outputs will be stored in OUTFILE.
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1);
opts_values = struct('BinningSpatial', 2.^[0:4], 'BinningTemp',2.^[0:4]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
% Arguments validation:
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p, RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
% outFile = {};
clear p
%%%%
% Get existing ImagesClassification files in directory:
existing_ChanList  = dir(fullfile(SaveFolder,'*.dat'));
idxName = ismember({existing_ChanList.name}, default_Output);
existing_ChanList = existing_ChanList(idxName);
% Calls function from IOI library. Temporary for now.
ImagesClassification(RawFolder, SaveFolder, opts.BinningSpatial, opts.BinningTemp,0);

% For Dual-Camera Imaging systems, apply the coregistration using the tform
% file created in DataViewer's OiS Dual Cam Coregistration utility:
if isfile(fullfile(SaveFolder,'AcqInfos.mat'))
    info = load(fullfile(SaveFolder, 'AcqInfos.mat'));
    if info.AcqInfoStream.MultiCam
        disp('Dual Camera data found!')
        % If this is a dual-cam system, look for tform data:
        % Find tform file:
        if ispc
            %     For windows:
            root = getenv('USERPROFILE');
        else
            %     For unix:
            root = getenv('HOME');
        end
        % If the folder doesn't exist, create one:
        LabeoFolder = fullfile(root,'Documents', 'LabeoTech', 'Config','umIT','tformFiles');
        % Get the list of tform files:        
        if isfile(fullfile(LabeoFolder, 'coreg2cam_tform.mat'))
            % Load tform file and perform coregistration
            tf = load(fullfile(LabeoFolder, 'coreg2cam_tform.mat'));
            disp('Applying coregistration to data from camera #2...');
            [status,warnmsg] = applyTform2Cams(SaveFolder,tf.tform, tf.tformInfo);
            if ~status
                warning(['Coregistration failed! Data import will resume without coregistration. ',warnmsg]);
            else
                disp('Done!')
            end
        else
            warning('Coregistration TFORM file not found! Data from camera #2 will not be coregistered with camera #1. Data import will resume without coregistration.')
        end                
    end
end

% Get only new files created during ImagesClassification:
chanList = dir(fullfile(SaveFolder,'*.dat'));
idx = ismember({chanList.name}, default_Output);
chanList = chanList(idx);
idxName = ismember({chanList.name}, {existing_ChanList.name});
idxDate = ismember([chanList.datenum], [existing_ChanList.datenum]);
idxNew = ~all([idxName; idxDate],1);
chanList = {chanList(idxNew).name};
outFile = fullfile(SaveFolder, chanList);
%
end