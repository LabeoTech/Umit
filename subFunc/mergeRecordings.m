function mergeRecordings(SaveFilename,folderList,filename,varargin)
% MERGERECORDINGS concatenates Image time series (Y,X,T) obtained with a
% LabeoTech Optical Imaging system. The data is stored in a .dat file ("filename")
% located across the folders listed in "folderList". An "events.mat" file
% will be created for the merged data.
% Inputs:
%   SaveFilename(char): full path of the ".dat" file with the merged data.
%   folderList (cell): array of full paths of the folders containing the
%       files to be merged.
%   filename (char): name of the file to be merged.
% Optional inputs:
%
%   merge_order (int): array of integers with indices of "folderList". The data will be
%       merged in this order. If not provided, the data will be merged in
%       the ascending order of "folderList".
%   trialNames (cell array): provide a list of trial names (char)
%       to tag triggers. The length of this array must be equal to the number
%       of files to be merged.
%   b_IgnoreEvents (bool | default = false): If TRUE, the function will ignore all events
%       stored in the "events.mat" file from the source folders (folderList)
%       and will create new timestamps marking the first and last frames from
%       each source file.


%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFilename', @(x) ischar(x) & ~isempty(x));
addRequired(p, 'folderList', @(x) iscell(x) & ~isempty(x));
addRequired(p, 'filename', @(x) ischar(x) & ~isempty(x));
addOptional(p, 'merge_order',[], @isnumeric);
addOptional(p, 'trialNames',{}, @iscell);
addOptional(p, 'b_IgnoreEvents',false, @islogical);
% Parse inputs:
parse(p,SaveFilename,folderList,filename, varargin{:})

%%%%%% Further input validation %%%%%%
% Set optional variables:
b_IgnoreEvents = p.Results.b_IgnoreEvents(1,1);
merge_order = round(p.Results.merge_order);
trialNames = cellfun(@num2str, p.Results.trialNames,'UniformOutput',false);% Force data to strings.
if isempty(merge_order)
    % If not provided, the order of merging will be the order of
    % "folderList":
    merge_order = 1:length(folderList);
end
if ~isempty(trialNames)
    assert(isequaln(numel(trialNames),numel(folderList)),'umIToolbox:mergeRecordings:WrongInput',...
        'The number of trial IDs must be the same as the number of input folders!');
end

% Append ".dat" to filenames if not done yet:
if ~endsWith(SaveFilename,'.dat')
    SaveFilename = [SaveFilename '.dat'];
end
if ~endsWith(filename,'.dat')
    filename = [filename '.dat'];
end
SaveFolder = fileparts(SaveFilename);
if isempty(SaveFolder)
    SaveFolder = pwd;
    SaveFilename = fullfile(SaveFolder,SaveFilename);
end

metaData_filename = strrep(filename, '.dat', '.mat');
% Check if merge_order contains all the indices of folderList:
assert(isequal(sort(merge_order),1:length(folderList)),'umIToolbox:mergeRecordings:MissingInput',...
    'The merge order list is incompatible with the list of folders');
% Reorder folderList following merge_order:
folderList = folderList(merge_order);
% check if the file exists in all folders:
idx = cellfun(@isfile, fullfile(folderList,filename));
if all(~idx)
    error(['The file ' filename ' was not found in any of the folders provided!']);
elseif ~all(idx)
    disp(repmat('-',1,100))
    warning('The following folders do not contain the file %s and will be ignored:\n%s\n',...
        filename, folderList{~idx});
    disp(repmat('-',1,100))
    % Update folderList:
    folderList = folderList(idx);
end
% Get full path for input data and meta data files:
metaDatNames = fullfile(folderList, metaData_filename);
datNames = fullfile(folderList, filename);
% Check if the associated .mat file exists in the folder:
assert(all(cellfun(@(x) isfile(x), metaDatNames)), 'umIToolbox:mergeRecordings:MissingInput',...
    'One or more associated .mat file are missing!')
% Check if the input files are image time series with dimensions {Y,X,T}:
idxDim = false(size(datNames));
idxSz = idxDim;
refSz = load(metaDatNames{1}, 'datSize');
mD = cell(1,length(datNames));
for ii = 1:length(datNames)
    mD{ii} = matfile(metaDatNames{ii});
    idxDim(ii) = all(ismember(mD{ii}.dim_names, {'Y','X','T'}));
    idxSz(ii) = isequaln(mD{ii}.datSize,refSz.datSize);
end
assert(all(idxDim), 'umIToolbox:mergeRecordings:WrongInput',...
    'This function accepts only image time series with dimensions {"Y", "X","T"}!');
% Check if all data have the same Y,X dimensions sizes:
assert(all(idxSz), 'umIToolbox:mergeRecordings:WrongInput','All input files must have the same Y,X sizes!');
% Check if input data is valid (i.e. if it is from LabeoTech)
reqFields = {'Freq', 'datName', 'datLength', 'FirstDim', 'dim_names','Datatype', 'datSize'};
assert(all(cellfun(@(x) all(ismember(reqFields, fieldnames(matfile(x)))), metaDatNames)),...
    'umIToolbox:mergeRecordings:WrongInput', 'Input data must be generated from LabeoTech Imaging systems');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data merge:
% Open new .mat (meta data) file:
mOut = struct();
% Get required fields meta data from the first file:
mOut.Freq = mD{1}.Freq;
mOut.datName = mD{1}.datName;
mOut.FirstDim = mD{1}.FirstDim;
mOut.dim_names = mD{1}.dim_names;
mOut.datFile = SaveFilename;
mOut.Datatype = mD{1}.Datatype;
mOut.datLength = [];
mOut.datSize = mD{1}.datSize;

Stim_trialNames = [];
evFileList = cell(size(folderList));

w = waitbar(0,'Merging data...', 'Name', ['Merging ' filename], 'WindowStyle','modal');
w.Resize = 'on';
w.Children.Title.Interpreter = 'none';
% Open new .dat file:
fidOut = fopen(SaveFilename,'w');
% Concatenate data in time domain:
for ii = 1:length(datNames)
    waitbar(ii/length(datNames),w);    
    % Merge data:
    w.Children.Title.String = {datNames{ii} ; ['[' num2str(ii) '/' num2str(length(datNames)) ' Reading data...]']};drawnow;
    fidIn = fopen(datNames{ii},'r');
    data = fread(fidIn,inf,['*' mD{ii}.Datatype]);
    w.Children.Title.String = {datNames{ii}; ['[' num2str(ii) '/' num2str(length(datNames)) ' Merging data...]']};drawnow;
    fwrite(fidOut,data, mD{ii}.Datatype);
    fclose(fidIn);    
    % Check if an "events.mat" file exists:
    if isfile(fullfile(folderList{ii},'events.mat'))
        evFileList{ii} = matfile(fullfile(folderList{ii},'events.mat'));
    end
    % Update meta data:
    mOut.datLength = sum([mOut.datLength, mD{ii}.datLength]); % Update data length
    % Populate "Stim_trialNames":
    SubStim = ii*(ones(1,mD{ii}.datLength)); SubStim([1,end]) = 0; % Use the first/last frame to mark the onset and offset of the Trial.
    Stim_trialNames = [Stim_trialNames, SubStim];    
end
fclose(fidOut);
w.Children.Title.String = 'Creating "events.mat" file...';pause(1);
% Create "events.mat" file:
if all(~cellfun(@isempty,evFileList)) && ~b_IgnoreEvents
    % If there are already "events.mat" files in ALL source folders,
    % concatenate them into a single file.
    % Replace existing "events.mat" file in the SaveFolder.
    warning('off')
    delete(fullfile(SaveFolder,'events.mat')); % Clear existing events.mat file
    warning('on')
    % Merge events:
    eventID = {};
    timestamps = [];
    state = [];
    eventNameList = {};  
    datLen = zeros(size(datNames));
    for ii = 1:length(metaDatNames)
        eventID{ii} = evFileList{ii}.eventID;
        state = [state; evFileList{ii}.state];
        eventNameList{ii,1} = evFileList{ii}.eventNameList;
        datLen(ii) = mD{ii}.datLength/mD{ii}.Freq;
        % Shift timestamps:
        timestamps = [timestamps; evFileList{ii}.timestamps + sum(datLen) - datLen(1)];
    end
    timestamps = single(timestamps); state = logical(state);
    allEventID = [];
    if isempty(trialNames)
        % Update eventID to match merged eventNameLists:
        allEventNameList = unique([eventNameList{:}],'stable');
        allEventNames = {};
        for ii = 1:length(eventNameList)
            allEventNames = [allEventNames; arrayfun(@(x) eventNameList{ii}(x), eventID{ii})];            
        end
        [~,allEventID] = cellfun(@(x) ismember(x,allEventNameList),allEventNames);
    else
        % Overwrite event IDs with trialNamess:
        for ii = 1:length(eventID)
            allEventID = [allEventID; repmat(ii, numel(eventID{ii}),1)];
        end
        allEventNameList = trialNames;
    end   
    allEventID = uint16(allEventID);    
else
    % Create new "events.mat" file using timestamps "Stim_trialNames".
    [allEventID, state, timestamps] = getEventFromStim(Stim_trialNames,mOut.Freq);
    if isempty(trialNames)
        allEventNameList = arrayfun(@num2str,unique(allEventID),'UniformOutput',false);
    else
        allEventNameList = trialNames;
    end
end

% Save event info to file:
saveEventsFile(SaveFolder,allEventID,timestamps,state,allEventNameList)
% Copy AcqInfo file from one of the original files to get some experiment info. This is used by some IOI_ana functions.
copyfile(fullfile(folderList{end}, 'AcqInfos.mat'), fullfile(SaveFolder,'AcqInfos.mat'));
close(w)
% Add data history to meta data file with info of this function:
myInfo = dir([mfilename('fullpath') '.m']);
opts = struct;
% opts.sourceFile = filename;
opts.sourceFolderList = folderList;
opts.merge_order = merge_order;
dH = genDataHistory(myInfo, opts, {filename},SaveFilename);
mOut.dataHistory = dH;
% Save meta data to file:
save(strrep(SaveFilename, '.dat', '.mat'), '-struct', 'mOut');
disp('Done')
end

function [ID,state,timestamps] = getEventFromStim(data, FrameRateHz)
ID = [];state = [];timestamps = [];
id_list = unique(data(:)); id_list(id_list == 0) = [];
for i = 1:length(id_list)
    on_indx = find(data(1:end-1)<.5 & data(2:end)>.5 & data(2:end) == id_list(i));
    off_indx = find(data(1:end-1)>.5 & data(2:end)<.5 & data(1:end-1) == id_list(i));
    timestamps =[timestamps; (sort([on_indx;off_indx]))./FrameRateHz];
    state =[state; repmat([true;false], numel(on_indx),1)];
    ID = [ID; repmat(id_list(i),numel([on_indx,off_indx]),1)];
end
% Rearrange arrays by chronological order:
[timestamps,idxTime] = sort(timestamps);
state = state(idxTime);
ID = ID(idxTime);
end
