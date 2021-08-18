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
if ~isempty(ref_frame_info.datFile)
    [~,ref_filename, ext] = fileparts(ref_frame_info.datFile);
    ref_filename = [ref_filename ext];
    % Try to load file from SaveFolder with the same name of Reference Frame file:
    try
        mData = mapDatFile(fullfile(SaveFolder, ref_filename));
    catch
        answer = questdlg(['Cannot find file with name ' ref_filename ' in ' SaveFolder '.'],...
            'Failed to find target file', 'Load manually', 'Cancel', 'Cancel');
        if strcmp(answer, 'Load manually')
            cd(SaveFolder);
            [file, path] = uigetfile('*.dat', 'Select File to compare with reference frame', 'green.dat');
            if file == 0
                disp('Operation cancelled by user')
                return
            end
            mData = mapDatFile(fullfile(path,file));
        else
            disp('Operation cancelled by user')
            return
        end
    end
else
    cd(SaveFolder);
    [file, path] = uigetfile('*.dat', 'Select File to compare with reference frame', 'green.dat');
    if file == 0
        disp('Operation cancelled by user')
        return
    end
    mData = mapDatFile(fullfile(path,file));
end

% Load Data:
if size(mData.Data.data,3) < 100
    targetFr = mean(mData.Data.data);
else
    targetFr = mean(mData.Data.data(:,:,1:100),3);
end
% Preprocessing:
% Apply unsharp mask to data:
refFr = imsharpen(refFr, 'Radius', 1.5, 'Amount',1);
targetFr = imsharpen(targetFr, 'Radius', 1.5, 'Amount',1);
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
if isempty(movPts)
    disp('Operation cancelled by User')
    return
end
% Adjust moving points:
movPts_adj = cpcorr(movPts,fxPts,targetFr,refFr);
% Perform image registration:
tform = fitgeotrans(movPts_adj,fxPts,'nonreflectivesimilarity');
Rfixed = imref2d(size(refFr));
targetFr = imwarp(targetFr, tform,'nearest', 'OutputView', Rfixed);
% Show results:
str = split(object.SaveFolder, filesep);
str = str(2:end); % Remove root folder.
fig = figure('Name', strjoin(str, '-'), 'WindowState', 'maximized');
subplot(211);imshowpair(refFr, targetFr);
subplot(212);imshowpair(refFr, targetFr, 'montage');
waitfor(fig)
%%%%%%
answer = questdlg('Are you satisfied with the registration?', 'Image Registration result', 'Yes, proceed', 'No, redo', 'Cancel', 'Cancel');
switch answer
    case 'Yes, proceed'
        %%%%%
        mask_question = questdlg('Optional Operations:', 'Options', 'Apply Mask','None', 'None');
        switch mask_question
            case 'Apply Mask only'
                b_applyMask = true;
                b_crop2Mask = false;
            otherwise
                b_applyMask = false;
                b_crop2Mask = false;
        end
        warpData(object,tform,refFr,ref_frame_info,translation,Rfixed,SaveFolder, applyToFile, b_applyMask, b_crop2Mask)
    case 'No, redo'
        b_redo = true;
        disp('Operation cancellec by User')
end
disp('Done!');
close all;
end

function warpData(obj,tform,refFr,ref_frame_info, translation,Rfixed,SaveFolder, applyToFile, b_applyMask, b_crop2Mask)
try
    [mData, metaData_source] = mapDatFile(fullfile(SaveFolder, applyToFile));
catch ME
    causeException = MException('MATLAB:UMIToolbox:manual_alignFrames:FileNotFound', 'DAT file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% This section may be too greedy on RAM and not efficient...
data = mData.Data.data;
warp_data = zeros(size(refFr,1),size(refFr,2), size(data,3), 'single');
disp(['Performing alignment in data from ' applyToFile '...']);
tic;
for i = 1:size(warp_data,3)
    frame = data(:,:,i);
    frame = imtranslate(frame, translation, 'FillValues', 0, 'OutputView','full');
    frame = imwarp(frame, tform, 'nearest', 'OutputView', Rfixed);
    warp_data(:,:,i) = frame;
end
toc
disp('Alignment finished.')
%%%%%
% BregmaXY = ref_frame_info.BregmaXY;
% LambdaXY = ref_frame_info.LambdaXY;
%%%%%
if b_applyMask
    mask = ref_frame_info.logical_mask;
    if isempty(mask)
        msg = 'Logical Mask not found in ImagingReferenceFrame.mat';
        errID = 'MATLAB:UMIToolbox:manual_alignFrames:VariableNotFound';
        error(errID, msg);
    end
    sz = size(warp_data);
    mask = repmat(mask,1,1,sz(3));
    warp_data = warp_data(:);
    warp_data(~mask(:)) = nan;
    warp_data = reshape(warp_data,sz);
    disp('Mask applied')
%     if b_crop2Mask
%         [r,c] = find(mask);
%         warp_data = warp_data(min(r):max(r), min(c):max(c), :);
%         BregmaXY(1) = BregmaXY(1) - min(c)+1;
%         BregmaXY(2) = BregmaXY(2) - min(r)+1;
%         LambdaXY(1) = LambdaXY(1) - min(c)+1;
%         LambdaXY(2) = LambdaXY(2) - min(r)+1;
%         disp('Frames cropped.')
%     end
end
%%%
[~,filename,ext] = fileparts(mData.Filename);
outFile = [ 'aligned_' filename ext];
datFile = fullfile(SaveFolder, outFile);
save2Dat(datFile, warp_data, metaData_source.dim_names);
% Add Bregma and Lambda coordinates to file meta data:
% [~, metaData_target] = mapDatFile(datFile);
% metaData_target.Properties.Writable = true;
% metaData_target.refPt = ref_frame_info.refPt;
% % metaData_target.BregmaXY = BregmaXY;
% % metaData_target.LambdaXY = LambdaXY;
% % Inherit properties from applyToFile metadata (quick fix for the loophole on PIPELINEMANAGER when alignFrames is used):
% props = setdiff(properties(metaData_source), properties(metaData_target));
% for k = 1:length(props)
%     eval(['metaData_target.' props{k} '= metaData_source.' props{k} ';'])
% end
% Write to FilePtr:
write2FilePtr(obj, datFile, fullfile(SaveFolder, applyToFile))
uiwait(msgbox(['Aligned Frames saved in ' datFile],'Data Saved', 'help'));
end

function write2FilePtr(obj, datFile, inputFile)
% WRITE2FILEPTR writes the File information stored in structure FILEINFO
% in OBJ.TMP_TARGETOBJ.FILEPTR.

%Initialize
[~,metaData] = mapDatFile(datFile);
[~,input_metaData] = mapDatFile(inputFile);
[Path, filename,ext] = fileparts(datFile);
[~, input_filename, input_ext] = fileparts(inputFile);
root = getenv('Umitoolbox');
funcInfo = dir(fullfile(root, 'GUI','DataViz', 'manual_alignFrames.m'));
FileInfo = struct('Name', [filename ext], 'UUID', metaData.fileUUID, 'Folder', tokenizePath(Path, obj),...
    'InputFile_Path', tokenizePath(inputFile, obj), 'InputFile_UUID', input_metaData.fileUUID, ...
    'creationDateTime', datestr(now), 'FunctionInfo', ...
    struct('Name', 'manual_alignFrames', 'DateNum', funcInfo.datenum,...
    'Job', ['manual_alignFrames(object, ''' input_filename input_ext '''):'], 'opts', ''));
FileList = readFilePtr(obj);
% Check for Files already logged on FilePtr
idx = false(length(FileList),2);
for i = 1:length(FileList.Files)
    idx(i,1) = strcmp(FileInfo.Name, FileList.Files(i).Name);
    idx(i,2) = strcmp(FileInfo.FunctionInfo.Name, FileList.Files(i).FunctionInfo.Name);
end
idx = all(idx,2);
% If there are no Files
if isempty(FileList.Files)
    FileList.Files = FileInfo;
    % If there are files and one identical, replace it.
elseif ~isempty(FileList.Files) && any(idx)
    FileList.Files(idx) = FileInfo;
    % If there are files and none identical: Append
else
    FileList.Files = [FileList.Files; FileInfo];
end
for i = 1:numel(FileList.Files)
    FileList.Files(i).Folder = tokenizePath(FileList.Files(i).Folder, obj);
    FileList.Files(i).InputFile_Path = tokenizePath(FileList.Files(i).InputFile_Path, obj);
end
txt = jsonencode(FileList);
fid = fopen(obj.FilePtr, 'w');
fprintf(fid, '%s', txt);
fclose(fid);
end
function filePtr_struct = readFilePtr(obj)
% READFILEPTR loads the content of FILEPTR.JSON in a structure
% stored in OBJ.TMP_FILEPTR.
txt = fileread(obj.FilePtr);
a = jsondecode(txt);
for i = 1:numel(a.Files)
    a.Files(i).Folder = tokenizePath(a.Files(i).Folder, obj, 'detokenize');
    a.Files(i).InputFile_Path = tokenizePath(a.Files(i).InputFile_Path, obj, 'detokenize');
end
filePtr_struct = a;
end
