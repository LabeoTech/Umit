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
    causeException = MException('MATLAB:UMIToolbox:FileNotFound',...
        'Imaging Reference Frame file not found.');
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
targetFr = mData.Data.data(:,:,1);

% Preprocessing:
% Create unsharp mask for each image:
refFr = imgaussfilt(refFr, .5) - imgaussfilt(refFr, 8);
targetFr = imgaussfilt(targetFr, .5) - imgaussfilt(targetFr, 8);
% Normalize:
refFr = (refFr - min(refFr(:))) ./ (max(refFr(:)) - min(refFr(:)));
targetFr = (targetFr - min(targetFr(:))) ./ (max(targetFr(:)) - min(targetFr(:)));
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
answer = questdlg('Are you satisfied with the registration?',...
    'Image Registration result', 'Yes, proceed', 'No, redo', 'Cancel', 'Cancel');
switch answer
    case 'Yes, proceed'
        [~,outFileName,ext] = fileparts(applyToFile);
        outFileName = inputdlg('Save File as...', 'Type file name to save aligned data',...
            [1 35], {[outFileName '_aligned' ext]});        
        
        if isempty(outFileName)
            disp('Operation Cancelled by User')
            return
        end
        outFileName = outFileName{:};
        if ~endsWith(outFileName,ext)
            outFileName = [outFileName, ext];
        end
        warpData(tform,refFr,ref_frame_info,Rfixed,SaveFolder,applyToFile, outFileName)       
    case 'No, redo'
        b_redo = true;
        disp('Operation cancellec by User')
end
disp('Done!');
close all;
end

function warpData(tform,refFr,ref_frame_info, Rfixed,SaveFolder, applyToFile, SaveFileName)
try
    [mData, metaData_source] = mapDatFile(fullfile(SaveFolder, applyToFile));
catch ME
    causeException = MException('MATLAB:umIToolbox:manual_alignFrames:FileNotFound', 'DAT file not found.');
    addCause(ME, causeException);
    rethrow(ME)
end
% This section may be too greedy on RAM and not efficient...
warp_data = zeros(size(refFr,1),size(refFr,2), size(mData.Data.data,3), 'single');
h = waitbar(0,['Performing alignment in data from ' applyToFile ' ...']);
for i = 1:size(warp_data,3)        
    warp_data(:,:,i) = imwarp(mData.Data.data(:,:,i), tform, 'nearest', 'OutputView', Rfixed);
    waitbar(i/size(warp_data,3),h);
end
waitbar(1,h,'Alignment finished.')
pause(1)
close(h)
% Create metaData for aligned data:
% Add image parameters from "ImagingReferenceFrame.mat" file
% as well as the previous metaData fields:
extraParams = load(metaData_source.Properties.Source);
extraParams.refPt = ref_frame_info.refPt;
extraParams.pxPermm = ref_frame_info.pxPermm;
metaData = genMetaData(warp_data, extraParams.dim_names, extraParams);
% Add this function to "dataHistory" variable inside metaData:

fcnInfo = dir([mfilename('fullpath') '.m']);
fcnInfo.name = fcnInfo.name(1:end-2); % remove .m;
dtHist = genDataHistory(fcnInfo,['out = manual_alignFrames(object,''' applyToFile ''');'],...
    [],'none');
if isfield(metaData,'dataHistory')
    metaData.dataHistory = [metaData.dataHistory; dtHist];    
else
    metaData.dataHistory = dtHist; 
end
% Save Data to file:
disp('Saving aligned data to file...')
datFile = fullfile(SaveFolder, SaveFileName);
save2Dat(datFile, warp_data, metaData);
uiwait(msgbox(['Data saved in: ' datFile],'Data Saved', 'help'));
end
