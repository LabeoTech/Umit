function [outData, metaData] = alignFrames(data, metaData, object, varargin)
% ALIGNFRAMES uses phase correlation to align images 
% to a reference frame created using the ROImanager app.

% Inputs:
%   data: 3D numerical matrix containing image time series with dimensions {Y, X, T}.
%   metaData: .mat file with meta data associated with "data".
%   object: umIT's imaging object handle.
%   opts (optional): structure containing extra parameters.
% Outputs:
%   outData: 3D numerical matrix with dimensions {Y,X,T} containing aligned frames.   
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_Output = 'mov_aligned.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct('UseFile', 'self');

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(3) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p,'object', @(x) isa(x,'Modality') || isa(x,'Acquisition'));
addOptional(p,'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ischar(x.UseFile));
% Parse inputs:
parse(p,data, metaData, object, varargin{:});
% Initialize Variables:
data = p.Results.data;
metaData = p.Results.metaData;
object = p.Results.object;
opts = p.Results.opts;
clear p

% Further validation of optional parameter "opts.UseFile":
errID = 'MATLAB:UMIToolbox:InvalidInput';
errMsg = 'Invalid entry for "UseFile" field. Input must be "self", "auto" or a name of a .dat file';
validFcn = @(x) ismember(x, {'self','auto'}) || endsWith(x, '.dat');
assert(validFcn(opts.UseFile, errID, errMsg));
%%%%

% Look for reference image in Subject's folder:
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
% Load frame from file to compare with reference image:
switch opts.UseFile
    case 'self'
        % Use first frame from "data"
        targetFr = data(:,:,1);
    case 'auto'
        % Look for a file with the same name as the one used to create the
        % "ImagingReferenceFrame.mat" file:
        [~,filename,ext] = fileparts(ref_frame_info.datFile);
        try
            targetDat = mapDatFile(fullfile(object.SaveFolder, [filename,ext]));
            targetFr = targetDat.Data.data(:,:,1);
        catch ME
            causeException = MException('MATLAB:UMIToolbox:FileNotFound', ['Cannot find "' filename '" in object''s SaveFolder']);
            addCause(ME, causeException);
            rethrow(ME)
        end
    otherwise
        % Load the filename in "opts.UseFile"
        try
            targetDat = mapDatFile(fullfile(object.SaveFolder, opts.UseFile));
            targetFr = targetDat.Data.data(:,:,1);
        catch ME
            causeException = MException('MATLAB:UMIToolbox:FileNotFound', ['Cannot find "' opts.UseFile '" in object''s SaveFolder']);
            addCause(ME, causeException);
            rethrow(ME)
        end
end

% Load Reference Frame;
refFr = ref_frame_info.reference_frame;
% refFr_mask = imsharpen(refFr ,'Radius', 1.5, 'Amount', 1);

% MV method (use the unsharp mask to do the registration:)
refFr_mask = imgaussfilt(refFr, .5) - imgaussfilt(refFr, 8);

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
    causeException = MException('MATLAB:UMIToolbox:MissingOutput',...
        'your version of the built-in MATLAB function "imregcorr" does not provide "peak" as output. You need to add it to the function and try again.');
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
% For Visual quality control of alignment:
figure('Name', strjoin({object.MyParent.MyParent.ID object.MyParent.ID object.ID}, '-'));
subplot(211);imshowpair(refFr_mask, targetFr_mask);
subplot(212);imshowpair(refFr_mask, targetFr_mask, 'montage');drawnow;
%%%%%%

% Apply mask to data file:
h = waitbar(0,'Initiating alignment...');
outData = zeros(size(refFr_mask,1),size(refFr_mask,2), size(data,3), 'single');
for i = 1:size(outData,3)
    waitbar(h, i/size(outData,3), 'Performing alignment...')
    frame = data(:,:,i);
    frame = imwarp(frame, tform, 'nearest', 'OutputView', Rfixed); 
    outData(:,:,i) = frame;
end
waitbar(h, 1, 'Alignment finished!'); 
close(h);
%%%%%
% Create new metaData and add image parameters from "ImagingReferenceFrame.mat" file:
extraParams.refPt = ref_frame_info.refPt;
extraParams.pxPermm = ref_frame_info.pxPermm;
metaData = genMetaData(outData, metaData.dim_names, extraParams);
end