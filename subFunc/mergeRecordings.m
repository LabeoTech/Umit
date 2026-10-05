function mergeRecordings(SaveFilename,folderList,filename,varargin)
% MERGERECORDINGS concatenates image time series (Y,X,T) along T. The data
% is stored in a .dat file ("filename") located across the folders listed in
% "folderList". An "events.mat" file will be created for the merged data.
%
% Each input is described with loadMetaData (headered or legacy .dat), must
% have axes {Y,X,T}, and must share Y/X, class, frame rate, and exposure
% with the others (error Umitoolbox:mergeRecordings:incompatibleInputs
% otherwise). The merged file is written with the self-describing .dat
% header; no metadata .mat file is written next to it.
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
% Get full path for input data files:
datNames = fullfile(folderList, filename);
% Describe every input with its own metadata (headered or legacy .dat):
mD = cellfun(@loadMetaData, datNames, 'UniformOutput', false);
% Check if the input files are image time series with dimensions {Y,X,T}:
idxDim = cellfun(@(md) isequal(cellstr(string(md.dimNames)), {'Y','X','T'}), mD);
assert(all(idxDim), 'umIToolbox:mergeRecordings:WrongInput',...
    'This function accepts only image time series with dimensions {"Y", "X","T"}!');
% Check that all inputs can be concatenated in time:
iAssertCompatibleInputs(mD, datNames);
nFrames = cellfun(@(md) datAxisSize(md, 'T'), mD);
frameRateHz = double(mD{1}.frameRateHz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data merge:
Stim_trialNames = [];
evFileList = cell(size(folderList));

w = waitbar(0,'Merging data...', 'Name', ['Merging ' filename], 'WindowStyle','modal');
w.Resize = 'on';
w.Children.Title.Interpreter = 'none';
% Create the merged .dat file (headered: the inputs' common class, frame
% rate, and exposure; the summed length; the output's own name):
[~, outName] = fileparts(SaveFilename);
outHeader = datHeaderFromInfo(mD{1}, outName, ...
    'dimSizes', [datAxisSize(mD{1}, 'Y'), datAxisSize(mD{1}, 'X'), sum(nFrames)]);
try
    slabOut = spatialSlabIO('create', SaveFilename, outHeader);
    cleanupOut = onCleanup(@() spatialSlabIO('close', slabOut));
    % Concatenate data in time domain:
    for ii = 1:length(datNames)
        waitbar(ii/length(datNames),w);
        % Merge data:
        w.Children.Title.String = {datNames{ii} ; ['[' num2str(ii) '/' num2str(length(datNames)) ' Reading data...]']};drawnow;
        data = loadData(datNames{ii});
        w.Children.Title.String = {datNames{ii}; ['[' num2str(ii) '/' num2str(length(datNames)) ' Merging data...]']};drawnow;
        t0 = sum(nFrames(1:ii-1));
        spatialSlabIO('write', slabOut, 1:slabOut.Nx, data, t0 + (1:nFrames(ii)));
        clear data
        % Check if an "events.mat" file exists:
        if isfile(fullfile(folderList{ii},'events.mat'))
            evFileList{ii} = matfile(fullfile(folderList{ii},'events.mat'));
        end
        % Populate "Stim_trialNames":
        SubStim = ii*(ones(1,nFrames(ii))); SubStim([1,end]) = 0; % Use the first/last frame to mark the onset and offset of the Trial.
        Stim_trialNames = [Stim_trialNames, SubStim]; %#ok<AGROW>
    end
    spatialSlabIO('finalize', slabOut);
    clear cleanupOut
catch ME
    clear cleanupOut
    close(w)
    if isfile(SaveFilename)
        delete(SaveFilename);
    end
    rethrow(ME);
end
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
    for ii = 1:length(datNames)
        eventID{ii} = evFileList{ii}.eventID;
        state = [state; evFileList{ii}.state];
        eventNameList{ii,1} = evFileList{ii}.eventNameList;
        datLen(ii) = nFrames(ii)/frameRateHz;
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
    [allEventID, state, timestamps] = getEventFromStim(Stim_trialNames,frameRateHz);
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
% Keep the merged channel's manifest entry truthful (.dat header Phase 7b):
iUpdateMergedChannelLength(fullfile(SaveFolder,'AcqInfos.mat'), SaveFilename, sum(nFrames));
close(w)
% The merged .dat is headered, so no metadata .mat file is written.
disp('Done')
end

function iUpdateMergedChannelLength(acqInfoPath, mergedFile, mergedLength)
%IUPDATEMERGEDCHANNELLENGTH Set the merged file's ImportedChannels Length.
% The copied AcqInfos.mat describes the last recording. Its entry for the
% merged file (if any) gets the merged file's frame count; nothing is added.

S = load(acqInfoPath, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream') || ~isfield(S.AcqInfoStream, 'ImportedChannels') || ...
        isempty(S.AcqInfoStream.ImportedChannels)
    return
end
AcqInfoStream = S.AcqInfoStream;
[~, mergedName, mergedExt] = fileparts(mergedFile);
idx = find(strcmpi(cellstr(string({AcqInfoStream.ImportedChannels.DatFile})), ...
    [mergedName, mergedExt]));
if isempty(idx)
    return
end
[AcqInfoStream.ImportedChannels(idx).Length] = deal(double(mergedLength));
save(acqInfoPath, 'AcqInfoStream');
end

function iAssertCompatibleInputs(mD, datNames)
%IASSERTCOMPATIBLEINPUTS Inputs must share Y/X, class, frame rate, and exposure.
props = {'frame size (Y, X)', @(md) [datAxisSize(md, 'Y'), datAxisSize(md, 'X')]; ...
    'data class', @(md) char(md.dataClass); ...
    'frame rate', @(md) double(md.frameRateHz); ...
    'exposure', @(md) double(md.exposureMsec)};
for iProp = 1:size(props, 1)
    ref = props{iProp, 2}(mD{1});
    for ii = 2:numel(mD)
        if ~isequaln(props{iProp, 2}(mD{ii}), ref)
            error('Umitoolbox:mergeRecordings:incompatibleInputs', ...
                'Cannot merge "%s" with "%s": the %s differs.', ...
                datNames{ii}, datNames{1}, props{iProp, 1});
        end
    end
end
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
