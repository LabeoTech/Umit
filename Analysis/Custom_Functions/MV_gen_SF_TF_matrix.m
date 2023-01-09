function outFile = MV_gen_SF_TF_matrix(SaveFolder, varargin)
% MV_GEN_SF_TF_MATRIX is a custom function from Matthieu Vanni's lab.
% This function finds the centroid pixels of visual areas and
% values from SF-TF experiments for statistical analysis.
%
% Inputs:
%   SaveFolder: Folder where data are stored.
%   Output: name of .DAT file saved in SAVEFOLDER.
%   datFileName: name of the .DAT file containing the SF-TF data
% Outputs:
%   outFile: name of aligned file.
% 

% Defaults:
default_opts = struct('datFileName', 'data_splitByEvent');
opts_values = struct('datFileName', {{'data_splitByEvent'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
default_Output = 'SF_TF_matrix.mat'; 
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
% Imaging Object:
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,SaveFolder, varargin{:});
% Initialize Variables:
File= p.Results.opts.datFileName;
SaveFolder = p.Results.SaveFolder;
outFile = p.Results.Output;

if ~endsWith(File,'.dat')
    File = [File '.dat'];
end

% Control for types of data
[mDat, metaData] = mapDatFile(fullfile(SaveFolder, File));
fields = fieldnames(mDat.Data);
errID = 'MATLAB:UMIToolbox:WrongFile';
errMsg = 'Wrong .dat file. Accepts only Resting State or SFTF experients';
assert(any(ismember(fields, {'CM', 'AVG'})), errID, errMsg)
%%%
% Check for anatomical landmarks:
errID = 'MATLAB:UMIToolbox:MissingVariable';
errMsg = 'File meta data does not containt Bregma and Lambda coordinates. Run "alignFrames" function and try again.';
assert(isprop(metaData, 'BregmaXY'), errID, errMsg);
% Load Allen Brain mask:
root = mfilename('fullpath');
str = strfind(root, 'Analysis');
root = root(1:str-1);
a = load(fullfile(root, 'GUI', 'DataViz', 'mouse_ctx_borders.mat'));

% SECTION TO BE CREATED:
% Data transformation for Resting state experiments:
if isfield(mDat.Data, 'CM')
    %%% CODE HERE %%%
    
else
    data = permute(mDat.Data.AVG,[2 3 1 4]);
    data = flipud(rot90(data));
    data = permute(data,[3 1 2 4]);
end
BregmaXY = metaData.BregmaXY;
LambdaXY = metaData.LambdaXY;
sz_data = size(data); 
% Resize image to have max dimension = 128 pixels:
mxSz = 128;
if max(sz_data([2 4])) > mxSz 
    dataAspectRatio = mxSz/max(sz_data([2 3]));
    tmp = imresize(squeeze(data(:,:,1,1)), dataAspectRatio);
    new_sz = [sz_data(1) size(tmp) sz_data(4)];
    data_shrinked = zeros(new_sz, 'single');
    for i = 1:sz_data(4)
        data_shrinked(:,:,:,i) = imresize3(data(:,:,:,i), new_sz(1:3), 'nearest');
    end
    BregmaXY = round(BregmaXY*dataAspectRatio);
    LambdaXY = round(LambdaXY*dataAspectRatio);
    clear tmp* data
else
    new_sz = sz_data;
    data_shrinked = data;
end

% Fit Allen Brain mask to data:
tform = fitgeotrans([a.atlas.BX a.atlas.BY ; a.atlas.LX a.atlas.LY],[BregmaXY ; LambdaXY], 'nonreflectivesimilarity');
areasID = {'V1_R', 'V1_L', 'AL_R', 'AL_L', 'PM_R', 'PM_L', 'L_R', 'L_L', 'BC_L', 'BC_R', 'HL_L', 'HL_R', 'FL_L', 'FL_R'};
idx = cellfun(@(x) any(strcmp(x, areasID)), a.atlas.areatag);
areas = a.atlas.areatag(idx);
[Xnew, Ynew] = transformPointsForward(tform, a.atlas.sx(idx)', a.atlas.sy(idx)');
Xnew = round(Xnew);
Ynew = round(Ynew);
% Exclude points outside image
offPx_x = Xnew < 1 | Xnew > new_sz(2);
offPx_y = Ynew < 1 | Ynew > new_sz(3);
idx = any([offPx_x, offPx_y],2);
% Reassign X,Y coordinates of points outside figure to one (hoping that
% they will fall outside the ROI...)
Xnew(idx) = 1;
Ynew(idx) = 1;
% Calculate signal amplitude:
centralFr = ceil(metaData.preEventTime_sec*metaData.Freq);
sub_data = zeros(length(Xnew),new_sz(1), new_sz(4), 'single');
for i = 1:length(Xnew)
    tmp = data_shrinked(:,Ynew(i),Xnew(i),:);
    tmp = squeeze(tmp);
    sub_data(i,:,:) = tmp;
end
baseline = mean(sub_data(:,:,1:centralFr-1),3);
Amp_avg_array = max(sub_data(:,:,centralFr:end,:),[],3, 'omitnan') - baseline;
Amp_avg_array(idx,:) = 0;
idx = all(Amp_avg_array==0,2);
Amp_avg_array(idx,:) = nan;
% Generate SF-TF matrix:
eventNames = metaData.eventList;
eventNames = eventNames{1}; eventNames = cell2mat(eventNames);
sf_list = unique(eventNames(:,1));
tf_list = unique(eventNames(:,2));
data = zeros(length(tf_list), length(sf_list), numel(areasID), 'single');
for i = 1:length(sf_list)
    idx = eventNames(:,1) == sf_list(i);
    data(i,:,:) = Amp_avg_array(:,idx)';
end
% Create Meta Data variables:
labels{1} = tf_list;
labels{2} = sf_list;
labels{3} = areas;
dim_names = {'X','Y','O'};
% Save to .MAT file:
% save2Mat(fullfile(SaveFolder, outFile), data, labels{3}, dim_names) % for testing...
save(fullfile(SaveFolder, outFile), 'data', 'labels', 'dim_names');
end