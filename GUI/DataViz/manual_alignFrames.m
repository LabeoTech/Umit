function b_redo = manual_alignFrames(object, applyToFile)
% MANUAL_ALIGNFRAMES uses cpselect, cpcorr, and fitgeotrans to select,
% refine and align user-selected points to a reference frame. 
% Inputs: 
%   - object: Imaging object.
%   - applyToFile: file containing data that will be aligned.
% Output:
%   - b_redo : boolean indicating if the user wants to rerun the function (used with the GUI lsaToolboxGUI.

%Initialize:
b_redo = false;
SaveFolder = object.SaveFolder;
figName = object.ID;
% Map reference frame:
try
    idx = false;
    ParentObj = object.MyParent;
    while ~idx
        ParentObj = ParentObj.MyParent;
        figName = strjoin({ParentObj.ID, figName}, '-');
        if isa(ParentObj, 'Subject')
            idx = true;
        end
    end
        ref_frame_info = matfile(fullfile(ParentObj.SaveFolder, 'ImagingReferenceFrame.mat'));
catch ME
    causeException = MException('MATLAB:UMIToolbox:FileNotFound', 'Imaging Reference Frame file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% Load Reference Frame;
refFr = ref_frame_info.reference_frame;
% Get Reference frame file name:
[~,ref_filename, ext] = fileparts(ref_frame_info.datFile);
ref_filename = [ref_filename ext];
% Load file from SaveFolder with the same name of Reference Frame file:
mData = mapDatFile(fullfile(SaveFolder, ref_filename));
% Load Data:
if size(mData.Data.data,3) < 100
    targetFr = mean(mData.Data.data);
else
    targetFr = mean(mData.Data.data(:,:,1:100),3);
end
% Apply unsharp mask to data:
targetFr = imsharpen(flipud(rot90(targetFr)), 'Radius', 1.5, 'Amount',1);
% Preprocessing:
% Normalize images:
refFr = (refFr - min(refFr(:)))./(max(refFr(:)) - min(refFr(:)));
targetFr = (targetFr - min(targetFr(:)))./(max(targetFr(:)) - min(targetFr(:)));
% Aling images' centers:
refFr_center(1) = mean(1:size(refFr,1));
refFr_center(2) = mean(1:size(refFr,2));
targetFr_center(1) = mean(1:size(targetFr,1));
targetFr_center(2) = mean(1:size(targetFr,2));
translation = (targetFr_center - refFr_center);
targetFr = imtranslate(targetFr, translation, 'FillValues', 0, 'OutputView','full');
% Lauch control points selector
[movPts, fxPts] = cpselect(targetFr, refFr,'Wait', true);
% Adjust moving points:
movPts_adj = cpcorr(movPts,fxPts,targetFr,refFr);
% Perform image registration:
tform = fitgeotrans(movPts_adj,fxPts,'nonreflectivesimilarity');
Rfixed = imref2d(size(refFr));
targetFr = imwarp(targetFr, tform, 'OutputView', Rfixed);
% Show results:
str = split(object.SaveFolder, filesep);
str = str(2:end); % Remove root folder.
fig = figure('Name', strjoin(str, '-'), 'WindowState', 'maximized');
imshowpair(refFr, targetFr);
waitfor(fig)
%%%%%%
answer = questdlg('Are you satisfied with the registration?', 'Image Registration result', 'Yes, proceed', 'No, redo', 'Cancel', 'Cancel');
switch answer
    case 'Yes, proceed'
        %%%%%
        mask_question = questdlg('Optional Operations:', 'Options', 'Apply Mask only', 'Apply Mask and crop', 'None', 'None');
        switch mask_question
            case 'Apply Mask only'
                b_applyMask = true;
                b_crop2Mask = false;
            case 'Apply Mask and crop'
                b_applyMask = true;
                b_crop2Mask = true;
            otherwise
                b_applyMask = false;
                b_crop2Mask = false;
        end
        warpData(tform,refFr,ref_frame_info,Rfixed,SaveFolder, applyToFile, b_applyMask, b_crop2Mask)
    case 'No, redo'
        b_redo = true;
end
close all;
end

function warpData(tform,refFr,ref_frame_info,Rfixed,SaveFolder, applyToFile, b_applyMask, b_crop2Mask)
try
    mData = mapDatFile(fullfile(SaveFolder, applyToFile));
    metaData_source = matfile(strrep(fullfile(SaveFolder, applyToFile), '.dat', '_info.mat'));
catch ME
    causeException = MException('MATLAB:UMIToolbox:FileNotFound', 'DAT file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% This section may be too greedy on RAM and not efficient...
data = flipud(rot90(mData.Data.data(:,:,:)));
warp_data = zeros(size(refFr,1),size(refFr,2), size(data,3), 'single');
disp(['Performing alignment in data from ' applyToFile '...']);
for i = 1:size(warp_data,3)
    frame = imwarp(data(:,:,i), tform, 'OutputView', Rfixed);
    warp_data(:,:,i) = frame;
end
disp('Alignment finished.')
%%%%%
BregmaXY = ref_frame_info.BregmaXY;
LambdaXY = ref_frame_info.LambdaXY;
%%%%%
if b_applyMask
    mask = ref_frame_info.logical_mask;
    if isempty(mask)
        msg = 'Logical Mask not found in ImagingReferenceFrame.mat';
        errID = 'MATLAB:UMIToolbox:VariableNotFound';
        error(errID, msg);
    end
    warp_data = bsxfun(@times, warp_data, mask);
    disp('Mask applied')
    if b_crop2Mask
        [r,c] = find(mask);
        warp_data = warp_data(min(r):max(r), min(c):max(c), :);
        BregmaXY(1) = BregmaXY(1) - min(c)+1;
        BregmaXY(2) = BregmaXY(2) - min(r)+1;
        LambdaXY(1) = LambdaXY(1) - min(c)+1;
        LambdaXY(2) = LambdaXY(2) - min(r)+1;
        disp('Frames cropped.')
    end
end
% Un-flip data before saving...
warp_data = flipud(rot90(warp_data));
%%%
% outFile = [outFile '.dat'];
% datFile = fullfile(SaveFolder, outFile);
datFile = fullfile(SaveFolder, 'mov_aligned.dat');
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, warp_data);
% Add Bregma and Lambda coordinates to file meta data:
% Add Bregma and Lambda coordinates to file meta data:
metaData_target = matfile(strrep(datFile, '.dat', '_info.mat'),'Writable', true);
metaData_target.BregmaXY = BregmaXY;
metaData_target.LambdaXY = LambdaXY;
% Inherit properties from applyToFile metadata (quick fix for the loophole on PIPELINEMANAGER when alignFrames is used):
props = setdiff(properties(metaData_source), properties(metaData_target));
for k = 1:length(props)
    eval(['metaData_target.' props{k} '= metaData_source.' props{k} ';'])
end
uiwait(msgbox(['Aligned Frames saved in ' datFile],'Data Saved', 'help'));
end
