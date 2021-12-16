function outFile = alignFrames(object, SaveFolder, varargin)
% ALIGNFRAMES uses phase correlation to align images 
% to a reference frame created using the ROImanager app.

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
default_opts = struct('ApplyToFile', 'fluo.dat', 'ApplyMask', false, 'Crop2Mask', false);
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
clear p

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
% refFr_mask = imsharpen(refFr ,'Radius', 1.5, 'Amount', 1);
% MV method (use the unsharp mask to do the registration:)
refFr_mask = imgaussfilt(refFr, .5) - imgaussfilt(refFr, 8);
% Get Reference frame file name:
[~,ref_filename, ext] = fileparts(ref_frame_info.datFile);
ref_filename = [ref_filename ext];
% Load file from SaveFolder with the same name of Reference Frame file:
mData = mapDatFile(fullfile(SaveFolder, ref_filename));
% Load Data:
if size(mData.Data.data,3) < 100
    targetFr = mean(mData.Data.data,3);
else
    targetFr = mean(mData.Data.data(:,:,1:100),3);
end
% Apply unsharp mask to data:
% targetFr_mask = imsharpen(targetFr,'Radius', 1.5, 'Amount', 1);
targetFr_mask = imgaussfilt(targetFr, .5) - imgaussfilt(targetFr, 8);

% Preprocessing:
% Normalize images:
% refFr_mask = (refFr_mask - min(refFr_mask(:)))./(max(refFr_mask(:)) - min(refFr_mask(:)));
% targetFr_mask = (targetFr_mask - min(targetFr_mask(:)))./(max(targetFr_mask(:)) - min(targetFr_mask(:)));
% Aling images' centers:
% refFr_center(1) = size(refFr,1)/2;
% refFr_center(2) = size(refFr,2)/2;
% targetFr_center(1) = size(targetFr,1)/2;
% targetFr_center(2) = size(targetFr,2)/2;
% translation = refFr_center - targetFr_center;
% targetFr_mask = imtranslate(targetFr_mask, fliplr(translation), 'FillValues', 0, 'OutputView','same');
% Perform image registration:
try
[tform, peak] = imregcorr(targetFr_mask,refFr_mask, 'similarity', 'Window', true);
Rfixed = imref2d(size(refFr));
catch ME
    causeException = MException('MATLAB:UMIToolbox:MissingOutput', 'your version of he built-in MATLAB function "imregcorr" does not provide "peak" as output. You need to add it to the function and try again.');
    addCause(ME, causeException);
    rethrow(ME)
end
if peak < 0.05
    disp('Phase correlation yielded a weak peak correlation value. Trying to apply intensity-based image registration instead ...')
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 0.000000625;
    optimizer.MaximumIterations = 80000;
    tform = imregtform(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask, imref2d(size(refFr_mask)),'similarity',optimizer,metric);
end
targetFr_mask = imwarp(targetFr_mask ,tform,'cubic', 'OutputView',Rfixed); % Chose "cubic" because "nearest" was showing stripes when rotating the target.
%%%%%
% For debugging:
figure('Name', strjoin({object.MyParent.MyParent.ID object.MyParent.ID object.ID}, '-'));
subplot(211);imshowpair(refFr_mask, targetFr_mask);
subplot(212);imshowpair(refFr_mask, targetFr_mask, 'montage');drawnow;
%%%%%%
% Apply mask to data file:
try
    [mData, metaData_source] = mapDatFile(fullfile(SaveFolder, opts.ApplyToFile));
catch ME
    causeException = MException('MATLAB:UMIToolbox:FileNotFound', 'DAT file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% This section may be too greedy on RAM and not efficient...
data = mData.Data.data(:,:,:);
warp_data = zeros(size(refFr_mask,1),size(refFr_mask,2), size(data,3), 'single');
disp(['Performing alignment in data from ' opts.ApplyToFile '...']);
tic;
for i = 1:size(warp_data,3)
    frame = data(:,:,i);
%     frame = imtranslate(frame, fliplr(translation), 'FillValues', 0, 'OutputView','same');
    frame = imwarp(frame, tform, 'nearest', 'OutputView', Rfixed); % Interp method = "nearest" to be "safe"...
    warp_data(:,:,i) = frame;
end
toc
disp('Alignment finished.')
%%%%%
refPt = ref_frame_info.refPt;

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
        refPt(1) = refPt(1) - min(c)+1;
        refPt(2) = refPt(2) - min(r)+1;
        disp('Frames cropped.')
    end
end
% Un-flip data before saving...
% warp_data = flipud(rot90(warp_data));
%%%
outFile = [ 'aligned_' opts.ApplyToFile];
datFile = fullfile(SaveFolder, outFile);
% Save to .DAT file and create .MAT file with metaData:
save2Dat(datFile, warp_data, metaData_source.dim_names);
% Add Bregma and Lambda coordinates to file meta data:
[~, metaData_target] = mapDatFile(datFile);
metaData_target.Properties.Writable = true;
metaData_target.refPt = refPt;
% metaData_target.LambdaXY = LambdaXY;
% Inherit properties from opts.Apply2File metadata (quick fix for the loophole on PIPELINEMANAGER when alignFrames is used):
props = setdiff(properties(metaData_source), properties(metaData_target));
for k = 1:length(props)
    eval(['metaData_target.' props{k} '= metaData_source.' props{k} ';'])
end

end