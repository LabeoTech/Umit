function outFile = MVcalculate_SF_TF_average(File, SaveFolder, varargin)
% MVCALCULATE_SF_TF_AVERAGE is a custom function used in the experimental
% protocol of Stroke project 2021 at Matthieu Vanni's. 
% This function averages specific conditions in the SF/TF experiments.
%
% Inputs:
%   File: fullpath of data_splitByEvent.DAT file (output of event_triggered_average function).
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
% Output:
%   outFile: name of Output file.
%   .DAT file containing average (AVG) and standard deviation (STD) movies
%   per event with dimensions (X,Y,T,Event).

% Defaults:
default_Output = 'SF_TF_AVG.dat'; 

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'File',@isfile)
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
Output = p.Results.Output;
%%%%
% Map movie and metadata to memory:
[mData, metaDat] = mapDatFile(File);
szdat = size(mData.Data.data);
% Check if file has 4-D:
if numel(szdat) ~= 4
    error('Wrong input file. Input must have 4 dimensions.');
end
% Combine Directions and calculate average and Standard deviation movies.
eventNameList = metaDat.eventNameList; 
eventNameList = cell2mat(eventNameList);
a = unique(eventNameList(:,1:3), 'rows');
AVG = zeros(size(a,1), szdat(2), szdat(3),szdat(4), 'single');
STD = zeros(size(a,1), szdat(2), szdat(3),szdat(4), 'single');
for i = 1:size(a,1)
    idx = find(all(eventNameList(:,1:3) == a(i,:),2));
    indx_trials = arrayfun(@(x) metaDat.eventID == x, idx, 'UniformOutput', false); indx_trials = any(cell2mat(indx_trials'),2);
    raw_data_subset = mData.Data.data(indx_trials,:,:,:);
    AVG(i,:,:,:)= mean(raw_data_subset,1, 'omitnan');
    STD(i,:,:,:)= std(raw_data_subset,0,1,'omitnan');
end
datFile = fullfile(SaveFolder, Output);
% Create MetaData structure:
szAVG = size(AVG);
metaDat = struct('datName', {'AVG', 'STD'}, 'datSize', {szAVG([1 2]), szAVG([1 2])},...
    'datLength', {szAVG(3:end) , szAVG(3:end)}, 'Datatype', {'single', 'single'},...
    'datFile', datFile);
metaDat(1).eventList = num2cell(a,2); % Added event description from events.mat file to the metadata.
% Save AVG and METADAT to DATFILE:
dim_names = {'E','X','Y', 'T'};
save2Dat(datFile, AVG, dim_names,'-w', metaDat)
% Append "STD" to DATFILE:
save2Dat(datFile, STD, dim_names, '-a')
% Output file names
outFile = Output;
end



