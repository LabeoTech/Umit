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

% For Retrocompatibility:
% Check for existence of Binning, rotation and X/Y offset fields
fNames = {'Binning','Rotation', 'X_Offset', 'Y_Offset'};
defaults = [1 0 0 0];
fn = fieldnames(AcqInfo);
for ii = 1:length(fNames)
    % Add defaults:
    if ~ismember(fNames{ii}, fn)
        AcqInfo.(fNames{ii}) = defaults(ii);
    end
    if ~ismember(fNames{ii}, fn)
        tformInfo.(fNames{ii}) = defaults(ii);
    end
end
% Update spatial reference object:
fileInfo = load(fullfile(DataFolder,strrep(Cam2List{1},'.dat','.mat')));
% Account for rotation in acquisition software:
rot_diff = 90*(AcqInfo.Rotation) - 90*(tformInfo.Rotation);
% Update TFORM to account for differences in rotation, binning and ROI
% offset:
tform = updateTForm(tform,tformInfo,fileInfo,AcqInfo,rot_diff);
RA = imref2d(fileInfo.datSize);
clear fileInfo

% Apply tform to data from Camera 2:
for ii = 1:length(Cam2List)
    % Load Data and apply TFORM:
    fprintf('----------------------------------\n')
    fprintf('Coregistration of file: "%s"\n', Cam2List{ii})
    fprintf('\t- Loading data...\n')
    dat = loadDatFile(fullfile(DataFolder, Cam2List{ii}));
    fprintf('\t- Applying geometric transformation...\n')
    dat = imwarp(dat,RA, tform, 'nearest','OutputView', RA);
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
save(fullfile(DataFolder,'tform.mat'), 'tform');
end

% Local function

function newtform = updateTForm(tform,tf_info, img_info, acqInfo, ang)
% UPDATETFORM Updates a geometric transformation matrix (tform) to account for
% spatial binning, rotation, and region of interest (ROI) offset.
%
% INPUTS:
%   tform       - Original geometric transformation matrix.
%   tf_info     - Transformation information.
%   img_info    - Image information.
%   acqInfo     - Acquisition information.
%   ang         - Rotation angle in degrees.
%
% OUTPUT:
%   newtform    - Updated geometric transformation matrix.
%
% This function performs a series of transformations on the input tform to
% accommodate differences in spatial binning, rotation, and ROI offset between
% acquired images and the desired transformation. It returns a new geometric
% transformation matrix that combines these adjustments:

% 1. Process Spatial Binning
AcqBinFactor = acqInfo.Binning;
AcqBinFactor = AcqBinFactor * acqInfo.Width / img_info.datSize(2);
binFactor = AcqBinFactor / tf_info.Binning;
binningMat = [binFactor 0 0; 0 binFactor 0; 0 0 1];

% 2. Process XY Offset
Xoffset = acqInfo.X_Offset - tf_info.X_Offset;
Yoffset = acqInfo.Y_Offset - tf_info.Y_Offset;
offsetMat = [1 0 0; 0 1 0; Xoffset Yoffset 1];

% 3. Create Rotation and Centering Matrices
frSize = img_info.datSize;
centerImg = [1 0 0; 0 1 0; -frSize(2)/2 -frSize(1)/2 1];
centerImgFlip = [1 0 0; 0 1 0; -frSize(1)/2 -frSize(2)/2 1];
rot = [cosd(ang) sind(ang) 0; -sind(ang) cosd(ang) 0; 0 0 1];

% 4. Create New tform 
switch abs(ang)
    case 0
        newMat = binningMat * offsetMat * tform.T * inv(offsetMat) * inv(binningMat);        
    case {90, 270}
        
        newMat = centerImg * rot * inv(centerImgFlip) * binningMat * offsetMat * tform.T * inv(offsetMat) * inv(binningMat) * centerImgFlip * inv(rot) * inv(centerImg);        
    case 180
        newMat = centerImg * rot * inv(centerImg) * binningMat * offsetMat * tform.T * inv(offsetMat) * inv(binningMat) * centerImg * inv(rot) * inv(centerImg);
end
newMat(:,3) = [0;0;1];
newtform = affine2d(newMat);
end
