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
        [~,filename,ext] = fileparts(ref_frame_info.datFile);
        try
            [targetDat, targetMetaData]= mapDatFile(fullfile(object.SaveFolder, [filename,ext]));
            targetFr = targetDat.Data.data(:,:,1);           
        catch ME
            causeException = MException('MATLAB:UMIToolbox:alignFrame:FileNotFound',...
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
refFr_mask = imgaussfilt(refFr, .5) - imgaussfilt(refFr, 8);

% Apply unsharp mask to data:
% targetFr_mask = imsharpen(targetFr,'Radius', 1.5, 'Amount', 1);

targetFr_mask = imgaussfilt(targetFr, .5) - imgaussfilt(targetFr, 8);

% Perform image registration:
try
[tform, peak] = imregcorr(targetFr_mask,refFr_mask, 'similarity', 'Window', true);
Rfixed = imref2d(size(refFr));
catch ME
    causeException = MException('MATLAB:UMIToolbox:alignFrame:MissingOutput',...
        'your version of the built-in MATLAB function "imregcorr" does not provide "peak" as output. You need to add it to the function and try again.');
    addCause(ME, causeException);
    rethrow(ME)
end

% Make initial geometric transformation of target image if phase
% correlation was satisfactory. If not, try intensity-based registration
% directly:
if peak > 0.05
    targetFr_mask = imwarp(targetFr_mask ,tform,'cubic', 'OutputView',Rfixed); % Chose "cubic" because "nearest" was showing stripes when rotating the target.
else
    disp('Phase correlation yielded a weak peak correlation value. Applying intensity-based image registration directly...')
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
    tmpFr = imregister(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask,...
    imref2d(size(refFr_mask)),'similarity',optimizer,metric, 'DisplayOptimization', false);    
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
tform = imregtform(targetFr_mask, imref2d(size(targetFr_mask)),refFr_mask,...
    imref2d(size(refFr_mask)),'similarity',optimizer,metric);    
targetFr_mask = imwarp(targetFr_mask ,tform,'nearest', 'OutputView',Rfixed);
disp('Done.')
%%%%%
disp('Check figure to validate alignment.')
% For Visual quality control of alignment:
figure('Name', strjoin({object.MyParent.MyParent.ID object.MyParent.ID object.ID}, '-'));
subplot(211);imshowpair(refFr_mask, targetFr_mask);
subplot(212);imshowpair(refFr_mask, targetFr_mask, 'montage');drawnow;
%%%%%%

% Apply mask to data file:
h = waitbar(0,'Initiating alignment...');
outData = zeros(size(refFr_mask,1),size(refFr_mask,2), size(data,3), 'single');
for i = 1:size(outData,3)
    waitbar(i/size(outData,3), h, 'Performing alignment...')    
    outData(:,:,i) = imwarp(data(:,:,i), tform, 'nearest', 'OutputView', Rfixed);
end
waitbar(1,h, 'Alignment finished!'); pause(.5);
close(h);
%%%%%
% Create new metaData and add image parameters from "ImagingReferenceFrame.mat" file
% as well as the previous metaData fields:
extraParams = metaData;
extraParams.refPt = ref_frame_info.refPt;
extraParams.pxPermm = ref_frame_info.pxPermm;
metaData = genMetaData(outData, extraParams.dim_names, extraParams);
end
