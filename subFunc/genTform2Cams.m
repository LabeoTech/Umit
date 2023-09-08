function [tform, registeredImages, warnmsg] = genTform2Cams(DataFolder, b_enhanceContrast, b_isFlipped)
% GENTFORM2CAMS performs coregistration of data from a 2-camera imaging system (LabeoTech).
% This function generates the geometric transformations (TFORM) necessary to align the 
% data from the two cameras and produces registered images.
%
% Inputs:
%   - DataFolder (str): Full path to the folder containing .dat files from the 2-camera
%                     imaging system.
%   - b_enhanceContrast (logical): Set to TRUE to enhance image contrast before coregistration.
%   - b_isFlipped (logical): Set to TRUE if images from camera 2 are horizontally flipped.
%  
% Outputs:
%   - tform (affine2d): The computed geometric transformation to align the images.
%   - registeredImages (cell array): A cell array containing the preprocessed and aligned images
%                     from camera 1 and camera 2.
%   - warnmsg (char): A warning message indicating the status of the operation.

tform = []; % Initialize the transformation matrix
registeredImages = {}; % Initialize the cell array for registered images
warnmsg = ''; % Initialize the warning message

% Get Acquisition info:
if ~exist(fullfile(DataFolder, 'AcqInfos.mat'),'file')
    % If the 'AcqInfos.mat' file doesn't exist, generate an error message:
    warnmsg = ['AcqInfos.mat file not found in "' DataFolder '". Run ImagesClassification and try again'];
    return;
end

Info = matfile(fullfile(DataFolder,'AcqInfos.mat'));
AcqInfo = Info.AcqInfoStream;
clear Info;

% Get list of channels for each camera:
if( AcqInfo.MultiCam )
    NbIllum = sum(cellfun(@(X) contains(X, 'Illumination'), fieldnames(AcqInfo)));
    Cam1List = {};
    Cam2List = {};
    for ind = 1:NbIllum
        idx = AcqInfo.("Illumination" + int2str(ind)).CamIdx;
        chan = lower(AcqInfo.("Illumination" + int2str(ind)).Color);
        if( contains(chan, 'fluo') )
            chan = ['fluo_' chan(9:11)];
        end
        if( contains(chan, 'amber') )
            chan = 'yellow';
        end
        if( idx == 1 )
            Cam1List{end+1} = [chan '.dat'];%#ok
        else
            Cam2List{end+1} = [chan  '.dat'];%#ok
        end
    end
    clear NbIllum ind idx chan
else
    % If only one camera was used, print a message and exit:
    disp('Only one camera was used. No need to coregister images');
    return;
end

% Get data size from the first channel.
datInfo = mapDatFile(fullfile(DataFolder, Cam1List{1}));
cam1Img = zeros([length(Cam1List), size(datInfo.Data.data,1), size(datInfo.Data.data,2)], 'single');
cam2Img = zeros([length(Cam2List), size(datInfo.Data.data,1), size(datInfo.Data.data,2)], 'single');
clear datInfo

% Create a combined image from each camera to be used in the coregistration:
% Camera 1:
for ii = 1:length(Cam1List)
    dat = loadDatFile(fullfile(DataFolder,Cam1List{ii}));
    dat = mean(dat,3);
    if mean(dat(:) > 1e3)
        dat = dat./mean(dat(:));
    else
        dat = zeros(size(dat));
    end
    cam1Img(ii,:,:) = dat;
end
cam1Img = squeeze(sum(cam1Img,1));

% Check if signals recorded on camera 1 are too low to perform coregistration.
if ( mean(cam1Img(:) < 0.1) )
    warnmsg = 'Signals recorded on camera 1 are too low to perform coregistration.';
    return;
end

% Camera 2:
for ii = 1:length(Cam2List)
    dat = loadDatFile(fullfile(DataFolder,Cam2List{ii}));
    dat = mean(dat,3);
    if mean(dat(:) > 1e3)
        dat = dat./mean(dat(:));
    else
        dat = zeros(size(dat));
    end
    cam2Img(ii,:,:) = dat;
end
clear dat

% Camera 2:
cam2Img = squeeze(sum(cam2Img,1));

% Check if signals recorded on camera 2 are too low to perform coregistration.
if ( mean(cam2Img(:) < 0.1) )
    warnmsg = 'Signals recorded on camera 2 are too low to perform coregistration.';
    return;
end

% Normalize images:
cam1Img = (cam1Img - min(cam1Img(:)))./(max(cam1Img(:)) - min(cam1Img(:)));
cam2Img = (cam2Img - min(cam2Img(:)))./(max(cam2Img(:)) - min(cam2Img(:)));

if b_enhanceContrast
    % Enhance image contrast if enabled
    cam1Img = adapthisteq(cam1Img);
    cam2Img = adapthisteq(cam2Img);
end

initialTform = affine2d([1 0 0;0 1 0; 0 0 1]);
if b_isFlipped
    % If camera 2 image is flipped (mirrored) compared to camera 1 image, adjust the transformation.
    initialTform = affine2d( [-1.0 0 0; 0 1.0 0; size(cam2Img,1) 0 1]);
end

% Set coregistration parameters:
[opt, met] = imregconfig("multimodal");
opt.GrowthFactor = 1.05;
opt.Epsilon = 1.0e-8;
opt.InitialRadius = 1e-3;
opt.MaximumIterations = 500;

% Execute coregistration:
disp('Performing coregistration...')
tform = imregtform(cam2Img,cam1Img,'similarity',opt,met,'InitialTransformation',initialTform);
disp('Done!')

% Store registered images in a cell array.
registeredImages{1} = cam1Img;
registeredImages{2} = imwarp(cam2Img,tform,'OutputView',imref2d(size(cam1Img)));
end
