function outFile = MVcalculate_SF_TF_average(File, SaveFolder, varargin)
% MVCALCULATE_SF_TF_AVERAGE is a custom function used in the experimental
% protocol of Stroke project 2021 at Matthieu Vanni's. 
% This function averages specific conditions in the SF/TF experiments.
%
% Inputs:
%   File: fullpath of eventTrigAVG.DAT file (output of event_triggered_average function).
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
mData = mapDatFile(File);
metaDat = matfile(strrep(File, '.dat', '_info.mat'));
szdat = size(mData.Data.AVG);
% Combine Directions and calculate average and Standard deviation movies.
eventList = metaDat.eventList; 
eventList = cell2mat(eventList{1});
a = unique(eventList(:,1:3), 'rows');
AVG = zeros(szdat(1), szdat(2),szdat(3), size(a,1), 'single');
STD = zeros(szdat(1), szdat(2),szdat(3), size(a,1), 'single');
for i = 1:size(a,1)
    indx = find(all(eventList(:,1:3) == a(i,:),2));
    raw_data_subset = mData.Data.AVG(:,:,:,indx);
    AVG(:,:,:,i)= nanmean(raw_data_subset,4);
    STD(:,:,:,i)= std(raw_data_subset,0,4,'omitnan');
end
datFile = fullfile(SaveFolder, Output);
% Create MetaData structure:
szAVG = size(AVG);
metaDat = struct('datName', {'AVG', 'STD'}, 'datSize', {szAVG([1 2]), szAVG([1 2])},...
    'datLength', {szAVG(3:end) , szAVG(3:end)}, 'Datatype', {'single', 'single'}, 'datFile', datFile);
metaDat(1).eventList = num2cell(a,2); % Added event description from events.mat file to the metadata.
% Save AVG and METADAT to DATFILE:
save2Dat(datFile, AVG,'-w', metaDat)
% Append "STD" to DATFILE:
save2Dat(datFile, STD, '-a')
% Output file names
outFile = Output;
end



