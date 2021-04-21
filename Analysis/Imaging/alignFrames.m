function outFile = alignFrames(object, SaveFolder, varargin)
% ALIGNFRAMES uses phase correlation to align anatomical images (green
% channel) to a reference frame created using the app
% mouse_ref_frame_creator.
% Inputs:
% object: Imaging object pointing to data to be aligned.
% SaveFolder: Folder where data are stored.
% Output: name of .DAT file saved in SAVEFOLDER.
% Parameter:
%   -ApplyToFile: Target Data file. Name of .DAT file to which the image warping will be
% performed.
%   -ApplyMask: Apply logical mask to target Data file. (~ROIs == 0).
%   -Crop2Maks: Crop frames to fit ROI.
% Outputs:
% outFile: name of aligned file.
% 

% Defaults:
default_Output = 'mov_aligned.dat'; 
default_opts = struct('ApplyToFile', 'fChan_475.dat', 'ApplyMask', false, 'Crop2Mask', false);
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'object', @(x) isa(x,'Modality'));
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,object, SaveFolder, varargin{:});
% Initialize Variables:
object = p.Results.object;
SaveFolder = p.Results.SaveFolder;
Output = p.Results.Output;
opts = p.Results.opts;
% outFile = erase(Output, '.dat')
%%%%
% Map reference frame:
try
    idx = false;
    ParentObj = object.MyParent;
    while ~idx
        ParentObj = ParentObj.MyParent;
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
% MV method (use the unsharp mask to do the registration:)
refFr_mask = imgaussfilt(refFr, 1) - imgaussfilt(refFr, 3);
% Get Reference frame file name:
[~,ref_filename, ext] = fileparts(ref_frame_info.datFile);
ref_filename = [ref_filename ext];
% Load file from SaveFolder with the same name of Reference Frame file:
mData = mapDatFile(fullfile(SaveFolder, ref_filename));
% Load Data:
if size(mData.Data.data,3) < 100
    targetFr = flipud(rot90(mean(mData.Data.data,3)));
else
    targetFr = flipud(rot90(mean(mData.Data.data(:,:,1:100),3)));
end
% Apply unsharp mask to data:
targetFr = imsharpen(targetFr,'Radius', 1.5, 'Amount', 1);
targetFr_mask = imgaussfilt(targetFr, 1) - imgaussfilt(targetFr, 3);
% Preprocessing:
% Normalize images:
refFr_mask = (refFr_mask - min(refFr_mask(:)))./(max(refFr_mask(:)) - min(refFr_mask(:)));
targetFr_mask = (targetFr_mask - min(targetFr_mask(:)))./(max(targetFr_mask(:)) - min(targetFr_mask(:)));
% Aling images' centers:
refFr_center(1) = size(refFr,1)/2;
refFr_center(2) = size(refFr,2)/2;
targetFr_center(1) = size(targetFr,1)/2;
targetFr_center(2) = size(targetFr,2)/2;
translation = refFr_center - targetFr_center;
targetFr_mask = imtranslate(targetFr_mask, fliplr(translation), 'FillValues', 0, 'OutputView','same');
% Perform image registration:
try
[tform, peak] = imregcorr(targetFr_mask,refFr_mask, 'similarity');
catch ME
    causeException = MException('MATLAB:UMIToolbox:MissingOutput', 'your version of he built-in MATLAB function "imregcorr" does not provide "peak" as output. You need to add it to the function and try again.');
    addCause(ME, causeException);
    rethrow(ME)
end
if peak < 0.05
    disp('Phase correlation yielded a weak peak correlation value. Trying to apply intensity-based image registration...')
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 0.000000625;
    optimizer.MaximumIterations = 80000;
    tform = imregtform(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask, imref2d(size(refFr_mask)),'similarity',optimizer,metric);
end
Rfixed = imref2d(size(refFr));
targetFr = imwarp(targetFr,tform, 'OutputView',Rfixed);
%%%%%
% For debugging:
figure('Name', strjoin({object.MyParent.MyParent.ID object.MyParent.ID object.ID}, '-'));
subplot(211);imshowpair(refFr, targetFr);
subplot(212);imshowpair(refFr, targetFr, 'montage');
%%%%%%
% Apply mask to data file:
try
    mData = mapDatFile(fullfile(SaveFolder, opts.ApplyToFile));
    metaData_source = matfile(strrep(fullfile(SaveFolder, opts.ApplyToFile), '.dat', '_info.mat'));
catch ME
    causeException = MException('MATLAB:UMIToolbox:FileNotFound', 'DAT file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% This section may be too greedy on RAM and not efficient...
data = flipud(rot90(mData.Data.data(:,:,:)));
warp_data = zeros(size(refFr_mask,1),size(refFr_mask,2), size(data,3), 'single');
disp(['Performing alignment in data from ' opts.ApplyToFile '...']);
tic;
for i = 1:size(warp_data,3)
    frame = data(:,:,i);
    frame = imtranslate(frame, fliplr(translation), 'FillValues', 0, 'OutputView','same');
    frame = imwarp(frame, tform, 'OutputView', Rfixed);
    warp_data(:,:,i) = frame;
end
toc
disp('Alignment finished.')
%%%%%
BregmaXY = ref_frame_info.BregmaXY;
LambdaXY = ref_frame_info.LambdaXY;
%%%%%
if opts.ApplyMask
    mask = ref_frame_info.logical_mask;
    if isempty(mask)
        msg = 'Logical Mask not found in ImagingReferenceFrame.mat';
        errID = 'MATLAB:UMIToolbox:VariableNotFound';
        error(errID, msg);
    end
    warp_data = bsxfun(@times, warp_data, mask);
    disp('Mask applied')
    if opts.Crop2Mask
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
datFile = fullfile(SaveFolder, Output);
outFile = Output;
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, warp_data);
% Add Bregma and Lambda coordinates to file meta data:
metaData_target = matfile(strrep(datFile, '.dat', '_info.mat'),'Writable', true);
metaData_target.BregmaXY = BregmaXY;
metaData_target.LambdaXY = LambdaXY;
% Inherit properties from opts.Apply2File metadata (quick fix for the loophole on PIPELINEMANAGER when alignFrames is used):
props = setdiff(properties(metaData_source), properties(metaData_target));
for k = 1:length(props)
    eval(['metaData_target.' props{k} '= metaData_source.' props{k} ';'])
end

end