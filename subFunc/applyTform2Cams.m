function [status, warnmsg] = applyTform2Cams(DataFolder, tform,tformInfo)
% APPLYTFORM2CAMS -  Applies the geometric transformation (tform) to data
% files from Camera 2 (2-Camera systems from LabeoTech) in the given data 
% folder (DataFolder). The registered data is saved back to the original files.
%
%   Input:
%       - DataFolder(char): The path to the folder containing the data files.
%       - tform (affine2d): The geometric transformation to be applied to the data.
%       - tformInfo (struct): Structure containing extra parameters from
%       tform file.
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

% Check for existence of rotation and X/Y offset fields
fNames = {'Rotation', 'X_Offset', 'Y_Offset'};
fn = fieldnames(AcqInfo);
for ii = 1:length(fNames)
    % Add defaults:
    if ~ismember(fNames{ii}, fn)
        AcqInfo.(fNames{ii}) = 0;
    end
    if ~ismember(fNames{ii}, fn)
        tformInfo.(fNames{ii}) = 0;
    end
end
% Update spatial reference object:
fileInfo = load(fullfile(DataFolder,strrep(Cam2List{1},'.dat','.mat')));
% Account for rotation in acquisition software:
rot_diff = 90*(tformInfo.Rotation) - 90*(AcqInfo.Rotation);
RA = updateRA(tformInfo, fileInfo, AcqInfo);
clear fileInfo
% Apply tform to data from Camera 2:
for ii = 1:length(Cam2List)
    % Load Data and apply TFORM:
    fprintf('----------------------------------\n')    
    fprintf('Coregistration of file: "%s"\n', Cam2List{ii})
    fprintf('\t- Loading data...\n')       
    dat = loadDatFile(fullfile(DataFolder, Cam2List{ii})); 
    fprintf('\t- Applying geometric transformation...\n')
    if rot_diff
        dat = imrotate(dat,rot_diff);     
        dat = imwarp(dat,RA, tform, 'nearest','OutputView', RA);        
        dat = imrotate(dat,-rot_diff);
    else
        dat = imwarp(dat,RA, tform, 'nearest','OutputView', RA);
    end
    % Replace data with registered one:
    fprintf('\t- Overwriting data in .DAT file...\n')    
    fid = fopen(fullfile(DataFolder, Cam2List{ii}), 'w');
    fwrite(fid, dat, 'single');
    fclose(fid);
    fprintf('Done.\n')    
    fprintf('----------------------------------\n')    
end
status = true; % Set status to true if the operation is successful.
% Save copy of tform in Data folder:
save(fullfile(DataFolder,'tform.mat'), 'tform','rot_diff');
end

% Local function

function RA = updateRA(tf_info, img_info, acqInfo)
    % Update the spatial reference object (RA) based on transformation and acquisition information.
    %
    % This function takes transformation (tform), transformation information (tf_info),
    % image information (img_info), and acquisition information (acqInfo) as input and
    % updates the spatial reference object (RA) accordingly.
    %
    % Parameters:
    %   - tf_info (struct): Transformation information containing fields like 'X_Offset', and 'Y_Offset'.
    %   - img_info (struct): Image information containing data size and other attributes.
    %   - acqInfo (struct): Acquisition information containing binning, offsets, and other details.
    %
    % Returns:
    %   - RA (imref2d): The updated spatial reference object with adjusted X and Y offsets and pixel sizes.    
        
    % Process spatial binning:
    % Get binning from acquisition:
    AcqBinFactor = acqInfo.Binning;
    % Update binning from umIT:
    AcqBinFactor = AcqBinFactor * acqInfo.Width / img_info.datSize(2);
    % Get tform binning factor:
    binFactor = AcqBinFactor/tf_info.Binning; % Binning relative to tform file binning.    
    RA = imref2d(img_info.datSize, binFactor, binFactor);
    % Process ROIs:
    % Compensate for X and Y offsets:
    RA.XWorldLimits = RA.XWorldLimits + acqInfo.X_Offset - tf_info.X_Offset;
    RA.YWorldLimits = RA.YWorldLimits + acqInfo.Y_Offset - tf_info.Y_Offset;
end
