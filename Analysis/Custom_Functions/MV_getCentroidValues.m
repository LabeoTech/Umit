function outFile = MV_getCentroidValues(File, SaveFolder, varargin)
% MV_GETCENTROIDVALUES is a custom function from Matthieu Vanni's lab.
% This function finds the centroid pixels of visual areas and extract
% values from RESTING STATE AND SFTF data files.
%
% Inputs:
% File: Imaging data file containing Resting State OR SFTF data.
% SaveFolder: Folder where data are stored.
% Output: name of .DAT file saved in SAVEFOLDER.
% Outputs:
% outFile: name of aligned file.
% 

% Defaults:
default_Output = 'centroid_Data.dat'; 
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'File', @isfile);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
% Initialize Variables:
File= p.Results.File;
SaveFolder = p.Results.SaveFolder;
Output = p.Results.Output;

% Control for types of data
mDat = mapDatFile(File);
fields = fieldnames(mDat.Data);
errID = 'MATLAB:UMIToolbox:WrongFile';
errMsg = 'Wrong .dat file. Accepts only Resting State or SFTF experients';
assert(any(ismember(fields, {'CM', 'AVG'})), errID, errMsg)
%%%
% Load meta data and check for anatomical landmarks:
metaData = matfile(strrep(File, '.dat', '_info.mat'));
errID = 'MATLAB:UMIToolbox:MissingVariable';
errMsg = 'File meta data does not containt Bregma and Lambda coordinates. Run "alignFrames" function and try again.';
assert(isprop(metaData, 'BregmaXY'), errID, errMsg);
% Load Allen Brain mask:
root = getenv('UMIToolbox');
a = load(fullfile(root, 'GUI', 'DataViz', 'mouse_ctx_borders.mat'));

% SECTION TO BE CREATED:
% Data transformation for Resting state experiments:
if isfield(mDat.Data, 'CM')
    %%% CODE HERE %%%
    
else
    data = flipud(rot90(mDat.Data.AVG));
    data_std = flipud(rot90(mDat.Data.STD));
end
BregmaXY = metaData.BregmaXY;
LambdaXY = metaData.LambdaXY;
data_size = size(data);
% Resize image to have max dimension = 64 pixels:
if max(data_size([1 2])) > 64
    dataAspectRatio = 64/max(data_size([1 2]));
    tmp = imresize(squeeze(data(:,:,1,1)), dataAspectRatio);
    tmp_std = imresize(squeeze(data_std(:,:,1,1)), dataAspectRatio);
    new_sz = [size(tmp) data_size([3 4])];
    data_shrinked = zeros(new_sz, 'single');
    data_std_shrinked = data_shrinked;
    for i = 1:data_size(4)
        data_shrinked(:,:,:,i) = imresize3(data(:,:,:,i), new_sz(1:3));
        data_std_shrinked(:,:,:,i) = imresize3(data_std(:,:,:,i), new_sz(1:3));
    end
    newBregma = round(BregmaXY*dataAspectRatio);
    newLambda = round(LambdaXY*dataAspectRatio);
    clear tmp* data data_std
    
% Fit Allen Brain mask to data:
tform = fitgeotrans([a.atlas.BX a.atlas.BY ; a.atlas.LX a.atlas.LY],[newBregma ; newLambda], 'nonreflectivesimilarity');
areasID = {'V1_R', 'V1_L', 'AL_R', 'AL_L', 'PM_R', 'PM_L', 'L_R', 'L_L', 'BC_L', 'BC_R', 'HL_L', 'HL_R', 'FL_L', 'FL_R'};
idx = cellfun(@(x) any(strcmp(x, areasID)), a.atlas.areatag);
areas = a.atlas.areatag(idx);
[Xnew, Ynew] = transformPointsForward(tform, a.atlas.sx(idx)', a.atlas.sy(idx)');
Xnew = round(Xnew);
Ynew = round(Ynew);
% Exclude points outside image
offPx_x = Xnew < 1 | Xnew > new_sz(1);
offPx_y = Ynew < 1 | Ynew > new_sz(2);
idx = any([offPx_x, offPx_y],2);
% Reassing X,Y coordinates of points outside figure to one (hoping that
% they will fall outside the ROI...)
Xnew(idx) = 1;
Ynew(idx) = 1;
% Calculate signal amplitude:
centralFr = metaData.eventFrame; centralFr = centralFr{1};
baseline = mean(data_shrinked(Xnew, Ynew, 1:centralFr-1,:),3);
Amp = max(data_shrinked(Xnew, Ynew, centralFr:end,:),[],3) - baseline;
Amp(idx) = 0;
AmpStd = max(data_std_shrinked(Xnew, Ynew, centralFr:end,:),[],3) - baseline;
Amp_avg_array = zeros(length(Xnew),new_sz(4), 'single');
Amp_std_array = Amp_avg_array;
for i = 1:length(Xnew)
    Amp_avg_array(i,:) = squeeze(Amp(i,i,1,:));
    Amp_std_array(i,:) = squeeze(AmpStd(i,i,1,:));
end
idx = all(Amp_avg_array==0,2);
Amp_avg_array(idx,:) = nan;
Amp_std_array(idx,:) = nan;
% Save array to .dat file:
datFile = fullfile(SaveFolder, Output);
outFile = Output;
% Create MetaData structure:
szArr = size(Amp_avg_array);
metaDat = struct('datName', {'AVG', 'STD'}, 'datSize', {szArr([1 2]), szArr([1 2])},...
    'datLength', {szArr(3:end) , szArr(3:end)}, 'Datatype', {'single', 'single'}, 'datFile', datFile);
% Save AVG and METADAT to DATFILE:
save2Dat(datFile, Amp_avg_array,'-w', metaDat)
% Append "STD" to DATFILE:
save2Dat(datFile, Amp_std_array, '-a')
% Add Bregma and Lambda coordinates to file meta data:
metaData = matfile(strrep(datFile, '.dat', '_info.mat'),'Writable', true);
metaData.areaTags = areas;
end