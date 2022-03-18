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
default_opts = struct('UseFile', 'auto');

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
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
assert(validFcn(opts.UseFile), errID, errMsg);
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
    causeException = MException('MATLAB:UMIToolbox:alignFrame:FileNotFound',...
        'Imaging Reference Frame file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% Load frame from file to compare with reference image:
switch opts.UseFile
    case 'self'
        % Use first frame from "data"
        targetFr = data(:,:,1);
        targetMetaData = metaData;
    case 'auto'
        % Look for a file with the same name as the one used to create the
        % "ImagingReferenceFrame.mat" file:
        if isempty(ref_frame_info.datFile)
            error('MATLAB:UMIToolbox:alignFrames:MissingInput',...
                ['Failed to locate reference file.' ...
                'The path to dat file in ImagingReferenceFrame.mat file is empty.' ...
                'Try again without the "auto" option.']);
        end
        [~,filename,ext] = fileparts(ref_frame_info.datFile);
        try
            [targetDat, targetMetaData]= mapDatFile(fullfile(object.SaveFolder, [filename,ext]));
            targetFr = targetDat.Data.data(:,:,1);
        catch ME
            causeException = MException('MATLAB:UMIToolbox:alignFrames:FileNotFound',...
                ['Cannot find "' filename '" in object''s SaveFolder']);
            addCause(ME, causeException);
            rethrow(ME)
        end
    otherwise
        % Load the filename in "opts.UseFile"
        try
            [targetDat, targetMetaData] = mapDatFile(fullfile(object.SaveFolder, opts.UseFile));
            targetFr = targetDat.Data.data(:,:,1);
        catch ME
            causeException = MException('MATLAB:UMIToolbox:alignFrame:FileNotFound',...
                ['Cannot find "' opts.UseFile '" in object''s SaveFolder']);
            addCause(ME, causeException);
            rethrow(ME)
        end
end
targetMetaData.dim_names;

% Load Reference Frame;
refFr = ref_frame_info.reference_frame;
% refFr_mask = imsharpen(refFr ,'Radius', 1.5, 'Amount', 1);

% MV method (use the unsharp mask to do the registration:)
radius = 0.05*max(size(refFr));
refFr_mask = imgaussfilt(refFr, .5) - imgaussfilt(refFr, radius);

% Apply unsharp mask to data:
% targetFr_mask = imsharpen(targetFr,'Radius', 1.5, 'Amount', 1);

targetFr_mask = imgaussfilt(targetFr, .5) - imgaussfilt(targetFr, radius);

% Perform image registration:
tform_init = imregcorr(targetFr_mask,refFr_mask, 'similarity', 'Window', true);
Rfixed = imref2d(size(refFr));
% Check if the first approximation is better than no registration at all:
counts = histcounts2(refFr_mask(:), targetFr_mask(:),50);
MIbefore = mutual_information(counts);
tmpFr = imwarp(targetFr_mask, tform_init, 'nearest', 'OutputView', Rfixed);
counts = histcounts2(refFr_mask(:), tmpFr(:),50);
MIafter = mutual_information(counts);
if MIafter < MIbefore
    disp('Phase correlation yielded a poor registration. Applying intensity-based image registration directly...')
    clear tform_init
end
% Set of HyperParameters for image registration:
GF = [1.10, 1.05, 1.02, 1.01];
Eps = [1e-10, 1e-15, 1e-20,1e-25];
IR = [6.25e-3, 6.25e-5, 6.25e-8, 6.25e-10];
MaxIter = 10000;
MI = -1000;
[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = MaxIter;
% Iterate over optimization parameters to obtain highest mutual
% information:
disp('Optimizing image registration parameters...')
for i = 1:4
    
    optimizer.GrowthFactor = GF(i);
    optimizer.Epsilon = Eps(i);
    optimizer.InitialRadius = IR(i);
    if exist('tform_init', 'var')
        tmpFr = imregister(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask,...
            imref2d(size(refFr_mask)),'similarity',optimizer,metric,...
            'DisplayOptimization', false, 'InitialTransformation', tform_init);
    else
        tmpFr = imregister(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask,...
            imref2d(size(refFr_mask)),'similarity',optimizer,metric, 'DisplayOptimization', false);
    end
    counts = histcounts2(refFr_mask(:), tmpFr(:),metric.NumberOfHistogramBins);
    tmpMI = mutual_information(counts);
    if tmpMI<=MI
        idx = i-1;
        break
    else
        MI = tmpMI;
        idx = i;
    end
end
fprintf('Maximum Mutual Information obtained: %.4f\n',MI)
% Re-calculate tform from best optimizer params:
optimizer.GrowthFactor = GF(idx);
optimizer.Epsilon = Eps(idx);
optimizer.InitialRadius = IR(idx);
disp('Calculating geometric transformation...')
if exist('tform_init','var')
    tform = imregtform(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask,...
        imref2d(size(refFr_mask)),'similarity',optimizer,metric, 'InitialTransformation', tform_init);
else
    tform = imregtform(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask,...
        imref2d(size(refFr_mask)),'similarity',optimizer,metric);
end
targetFr_mask = imwarp(targetFr_mask ,tform,'nearest', 'OutputView',Rfixed);
disp('Done.')
%%%%%

% For Visual quality control of alignment:
fig = figure('Name', strjoin({object.MyParent.MyParent.ID object.MyParent.ID object.ID}, '-'),...
    'WindowButtonMotionFcn', @moveDot, 'Visible', 'off');
s1=subplot(2,2,(1:2));imshowpair(refFr_mask, targetFr_mask);
s2=subplot(223); imagesc(s2,refFr_mask); colormap(s2,'gray');axis(s2,'off')
s3=subplot(224); imagesc(s3,targetFr_mask); colormap(s3,'gray');axis(s3,'off')
set(s2, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual');
set(s3, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual');
title(s1,['Mutual Information = ' num2str(MI)]);
title(s2, 'Reference');
title(s3, 'Registered');
% draw dots
hold(s2,'on');
plot(s2,1,1,'g+', 'Tag', 'gDot'); hold(s2,'off');
hold(s3,'on');
plot(s3,1,1,'rx', 'Tag', 'rDot'); hold(s3,'off');

%%%%%%

% Apply mask to data file:
h = waitbar(0,'Initiating alignment...');
outData = zeros(size(refFr_mask,1),size(refFr_mask,2), size(data,3), 'single');
for i = 1:size(outData,3)
    waitbar(i/size(outData,3), h, 'Performing alignment...')
    outData(:,:,i) = imwarp(data(:,:,i), tform, 'nearest', 'OutputView', Rfixed);
end
waitbar(1,h, 'Alignment finished!'); pause(.5);
disp('Check figure to validate alignment.')
waitbar(1,h, 'Check the Figure!'); pause(.5);
close(h);
%%%%%
% Create new metaData and add image parameters from "ImagingReferenceFrame.mat" file
% as well as the previous metaData fields:
extraParams = metaData;
extraParams.refPt = ref_frame_info.refPt;
extraParams.pxPermm = ref_frame_info.pxPermm;
metaData = genMetaData(outData, extraParams.dim_names, extraParams);
% Show Figure
fig.Visible = 'on';
end


% Figure callbacks:

function moveDot(src,~)
% 
% disp('moving...')
for i = 1:length(src.Children)
    coords = get(src.Children(i), 'CurrentPoint');
    coords = round(coords(1,1:2));
    b_in = getBounds(coords, src.Children(i));
    if b_in
        break        
    end    
end
if ~b_in
    return
end
% Update dot positions:
dot1 = findall(src, 'Tag', 'gDot');
dot2 = findall(src, 'Tag', 'rDot');

dot1.XData = coords(1);
dot2.XData = coords(1);

dot1.YData = coords(2);
dot2.YData = coords(2);
end

function [b_in_bounds,ax] = getBounds(pt,ax)
% GETBOUNDS verifies if the current position of the mouse
% Output:
% get axis limits:
x_lims = get(ax, 'XLim');
y_lims = get(ax, 'YLim');

% Check if cursor is inside the axis limits:
b_in_bounds = (pt(1) >= x_lims(1)) && ...
    (pt(1) <= x_lims(2)) && ...
    (pt(2) >= y_lims(1)) && ...
    (pt(2) <= y_lims(2));
end

