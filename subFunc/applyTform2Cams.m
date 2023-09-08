function [status, warnmsg] = applyTform2Cams(DataFolder, tform)
% APPLYTFORM2CAMS -  Applies the geometric transformation (tform) to data
% files from Camera 2 (2-Camera systems from LabeoTech) in the given data 
% folder (DataFolder). The registered data is saved back to the original files.
%
%   Input:
%       - DataFolder(char): The path to the folder containing the data files.
%       - tform (affine2d): The geometric transformation to be applied to the data.
%   Output:
%       - status(bool): A status indicator. It is set to true if the operation is successful,
%                and false otherwise.
%       - warnmsg (char): A warning message indicating the status of the operation.

warnmsg = ''; status = false;
% Get Acquisition info:
if ~exist(fullfile(DataFolder, 'AcqInfos.mat'), 'file')
    % If the 'AcqInfos.mat' file doesn't exist, generate an error message:
    warnmsg = ['AcqInfos.mat file not found in "' DataFolder '". Run ImagesClassification and try again.'];
    return;
end
Info = matfile(fullfile(DataFolder, 'AcqInfos.mat'));
AcqInfo = Info.AcqInfoStream;
clear Info;
% Get list of channels for each camera:
if (AcqInfo.MultiCam)
    NbIllum = sum(cellfun(@(X) contains(X, 'Illumination'), fieldnames(AcqInfo)));
    Cam1List = {};
    Cam2List = {};
    for ind = 1:NbIllum
        idx = AcqInfo.("Illumination" + int2str(ind)).CamIdx;
        chan = lower(AcqInfo.("Illumination" + int2str(ind)).Color);
        if (contains(chan, 'fluo'))
            chan = ['fluo_' chan(9:11)];
        end
        if (contains(chan, 'amber'))
            chan = 'yellow';
        end
        if (idx == 1)
            Cam1List{end + 1} = [chan '.dat']; %#ok
        else
            Cam2List{end + 1} = [chan '.dat']; %#ok
        end
    end
    clear NbIllum ind idx chan
else
    % If only one camera was used, print a message and exit:
    disp('Only one camera was used. No need to coregister images');
    return;
end

% Apply tform to data from Camera 2:
for ii = 1:length(Cam2List)
    % Load Data and apply TFORM:
    disp(['Processing file: ' Cam2List{ii} '...'])
    disp('Loading data...')
    [dat, info] = loadDatFile(fullfile(DataFolder, Cam2List{ii}));
    disp('Applying geometric transformation...')
    dat = imwarp(dat, tform, 'OutputView', imref2d(info.datSize));
    % Replace data with registered one:
    disp('Overwriting data in .DAT file...');
    fid = fopen(fullfile(DataFolder, Cam2List{ii}), 'w');
    fwrite(fid, dat, 'single');
    fclose(fid);
    disp('Done')    
end
status = true; % Set status to true if the operation is successful.
% Save copy of tform in Data folder:
save(fullfile(DataFolder,'tform.mat'), 'tform');
end
