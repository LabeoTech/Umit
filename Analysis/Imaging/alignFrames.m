function [outData, metaData] = alignFrames(data, metaData, SaveFolder, object, varargin)
% ALIGNFRAMES uses phase correlation to align images
% to a reference frame created using the ROImanager app.

% Inputs:
%   data: 2D or 3D numerical matrix containing image time series with dimensions {Y, X} or {Y, X, T}.
%   metaData: .mat file with meta data associated with "data".
%   SaveFolder: Folder containing the .dat files.
%   object: umIT's imaging object handle.
%   opts (optional): structure containing extra parameters.
% Outputs:
%   outData: "data" aligned.
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_Output = 'mov_aligned.dat'; %#ok This line is here just for Pipeline management.
default_opts = struct('UseFile', 'auto', 'RefFile','ImagingReferenceFrame.mat');
opts_values = struct('UseFile',{{'auto'}}, 'RefFile',{{'ImagingReferenceFrame.mat'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ismember(ndims(x),[2 3])); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p,'SaveFolder',@isfolder);
addRequired(p,'object', @(x) isa(x,'Modality') || isa(x,'Acquisition'));
addOptional(p,'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ischar(x.UseFile));
% Parse inputs:
parse(p,data, metaData, SaveFolder,object, varargin{:});
% Initialize Variables:
opts = p.Results.opts;
clear p
% Further validation of input data dimensions:
errID = 'umIToolbox:alignFrames:InvalidInput';
errMsg = 'Data must be have dimensions X,Y or X,Y,T';
assert(all(ismember(metaData.dim_names, {'X','Y','T'})),errID,errMsg);
% Further validation of optional parameter "opts.UseFile":
errID = 'MATLAB:UMIToolbox:InvalidInput';
errMsg = 'Invalid entry for "UseFile" field. Input must be "auto" or a name of a .dat file';
validFcn = @(x) ischar(x) || strcmpi(x, 'auto');
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
    refFile = fullfile(ParentObj.SaveFolder, opts.RefFile);
    if ~endsWith(refFile, '.mat')
        refFile = [refFile, '.mat'];
    end
    if ~isfile(refFile)
        error('umIToolbox:alignFrames:FileNotFound', 'Imaging reference file not found in Subject folder!');
    end
    ref_frame_info = matfile(refFile);
catch ME
    causeException = MException('umIToolbox:alignFrame:FileNotFound',...
        'Imaging Reference Frame file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% Load frame from file to compare with reference image:
switch lower(opts.UseFile)
    case 'auto'
        % Look for a file with the same name as the one used to create the
        % "ImagingReferenceFrame.mat" file:
        if isempty(ref_frame_info.datFile)
            error('umIToolbox:alignFrames:MissingInput',...
                ['Failed to locate reference file. Was the reference frame generated from an existing .dat file?' ...
                'Type the file name to use instead of "auto".']);
        end
        [~,filename,ext] = fileparts(ref_frame_info.datFile);
        try
            targetDat = mapDatFile(fullfile(object.SaveFolder, [filename,ext]));
            targetFr = targetDat.Data.data(:,:,1);
        catch ME
            causeException = MException('umIToolbox:alignFrames:FileNotFound',...
                ['Cannot find "' filename '" in object''s SaveFolder']);
            addCause(ME, causeException);
            rethrow(ME)
        end
    otherwise
        if ~endsWith(opts.UseFile, '.dat')
            opts.UseFile = [opts.UseFile , '.dat'];
        end
        % Load the filename in "opts.UseFile"
        try
            targetDat = mapDatFile(fullfile(object.SaveFolder, opts.UseFile));
            if ndims(targetDat.Data.data) == 3
                targetFr = targetDat.Data.data(:,:,1);
            else
                targetFr = targetDat.Data.data;
            end
        catch ME
            causeException = MException('umIToolbox:alignFrame:FileNotFound',...
                ['Cannot find "' opts.UseFile '" in object''s SaveFolder']);
            addCause(ME, causeException);
            rethrow(ME)
        end
end

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
if any(size(refFr_mask) ~= size(targetFr_mask))
    tmp = imresize(targetFr_mask, size(refFr_mask));
    counts = histcounts2(refFr_mask(:), tmp(:),50);
    clear tmp
else
    counts = histcounts2(refFr_mask(:), targetFr_mask(:),50);
end
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
MaxIter = 1000;
MI = -1000;
[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = MaxIter;
% Iterate over optimization parameters to obtain highest mutual
% information:
disp('Optimizing image registration parameters...')
for i = 1:4
    disp(['Iteration ' num2str(i) '/4']);
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
    if tmpMI <= MI
        idx = i-1;
        disp('Optimization stopped here!');
        break
    else
        MI = tmpMI;
        idx = i;
    end
end
% Re-calculate tform from best optimizer params:
optimizer.GrowthFactor = GF(idx);
optimizer.Epsilon = Eps(idx);
optimizer.InitialRadius = IR(idx);
disp('Calculating geometric transformation using optimal parameters...')
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
myID = strjoin({object.MyParent.MyParent.ID object.MyParent.ID object.ID}, '-');
fig = figure('Name', myID, 'WindowButtonMotionFcn', @moveDot, 'Visible', 'off',...
    'UserData',struct('ID',myID, 'Folder',object.SaveFolder, 'ListFile',''), 'Tag','alignFig');
uicontrol('Style','pushbutton','String','Mark as misaligned','Tooltip','Save the name of the recording to a text file', ...
    'Parent',fig,'Position',[10 10 120,30], 'BackgroundColor',[.8 0.2 0.2],...
    'FontSize',10,'Callback', @save2List)
s1 = subplot(2,2,(1:2));imshowpair(refFr_mask, targetFr_mask);
s2 = subplot(223); imagesc(s2,refFr_mask); colormap(s2,'gray');axis(s2,'off')
s3 = subplot(224); imagesc(s3,targetFr_mask); colormap(s3,'gray');axis(s3,'off')
set(s2, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual');
set(s3, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual');
title(s1,'Merged');
title(s2, 'Reference');
title(s3, 'Registered');
% draw dots
hold(s2,'on');
plot(s2,1,1,'g+', 'Tag', 'gDot'); hold(s2,'off');
hold(s3,'on');
plot(s3,1,1,'rx', 'Tag', 'rDot'); hold(s3,'off');
% Link axes
linkaxes([s1,s2,s3], 'xy');
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
% Save "tform" and the frames used in the aligment in the savefolder:
save(fullfile(SaveFolder,'alignmentParams.mat'),'tform','refFr','targetFr');
% Show Figure
fig.Visible = 'on';
end
% Figure callbacks:

function moveDot(src,~)
%
% disp('moving...')
myAx = findobj(src, 'Type','axes');
for i = 1:length(myAx)
    coords = get(myAx(i), 'CurrentPoint');
    coords = round(coords(1,1:2));
    b_in = getBounds(coords, myAx(i));
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

function save2List(src,~)
% SAVE2LIST appends the name of the folder to a .TXT file. This is useful
% for the user to make a list of failed alignments.

% Check for other opened figures:
h = findall(0,'Tag','alignFig');
% Look for a text file stored in the "UserData" property of those figures:
txtFileName = arrayfun(@(x) x.UserData.ListFile,h,'UniformOutput',false);
txtFileName(cellfun(@isempty,txtFileName)) = [];
% Open dialog to save file:
if isempty(txtFileName)
    defName = fullfile(src.Parent.UserData.Folder,['alignFrameList_' datestr(now(),'yyyy-mm-dd') '.csv']);
    [file,path] = uiputfile('*.csv','Create or choose file to save to list',defName);
else
    defName = fullfile(src.Parent.UserData.Folder,txtFileName{1});
    [file,path] = uigetfile('*.csv','Choose existing file to save to list',defName);
end

if file == 0
    return
end
info = cell2table({src.Parent.UserData.ID, src.Parent.UserData.Folder, datestr(now(),'yyyy-mm-dd HH:MM:SS')});
info.Properties.VariableNames = {'ID','Folder','DateAdded'};
txt = [];
% Read the selected file or create a new one:
saveFile = fullfile(path,file);
if exist(saveFile, 'file')
    txt = readtable(saveFile);
end
% Create/append the folder and name:
if isempty(txt)
    txt = info;
else
    txt = [txt;info]; % Append info.
end
writetable(txt,saveFile);
% Change button color to indicate that the data was saved:
figure(src.Parent);
src.BackgroundColor = [.2 .8 .2];
src.String = 'Added to list!';
pause(1);
src.BackgroundColor = [.7 .7 .7];
% Put file name into UserData:
src.Parent.UserData.ListFile = fullfile(path,file);
end