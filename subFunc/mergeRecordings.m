function mergeRecordings(SaveFilename,folderList,filename, merge_type, varargin)
% MERGERECORDINGS concatenates Image time series (Y,X,T) obtained with a
% LabeoTech Optical Imaging system. The data is stored in a .dat file ("filename")
% located across the folders listed in "folderList".
% Inputs:
%   SaveFilename(char): full path of the ".dat" file with the merged data.
%   folderList (cell): array of full paths of the folders containing the
%       files to be merged.
%   filename (char): name of the file to be merged.
%   merge_type(char){'trial','movie'}: indicates how to merge the data.
%       If "trial" the output data is organized into a 4D array with
%       dimensions {E(vents), Y,X,T}. If "movie", the data will be
%       concatenated in the Time dimension.
%   merge_order (int, optional): array of integers with indices of "folderList". The data will be
%       merged in this order. If not provided, the data will be merged in
%       the ascending order of "folderList".
%   b_IgnoreStim (bool, optional): set to TRUE to skip concatenation of
%   "Stim" data from file's meta data.
%   trialMarker (cell array, optional): provide a list of trial IDs (char) to
%       create a "Stim" array containing the onset/offset timestamps for each trial/segment.
%       !! Beware that if a list is provided, the existing "Stim" data will be
%       overwritten!!


%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFilename', @(x) ischar(x) & ~isempty(x));
addRequired(p, 'folderList', @(x) iscell(x) & ~isempty(x));
addRequired(p, 'filename', @(x) ischar(x) & ~isempty(x));
addRequired(p, 'merge_type', @(x) ischar(x) & ismember(lower(x), {'trial', 'movie'}));
addOptional(p, 'merge_order',[], @isnumeric);
addOptional(p, 'b_IgnoreStim',false, @islogical);
addOptional(p, 'trialMarker',{}, @iscell);
% Parse inputs:
parse(p,SaveFilename,folderList,filename, merge_type,varargin{:})
% Set optional variables:
merge_order = p.Results.merge_order;
b_IgnoreStim = p.Results.b_IgnoreStim;
trialMarker = p.Results.trialMarker;
if isempty(merge_order)
    % If not provided, the order of merging will be the order of
    % "folderList":
    merge_order = 1:length(folderList);
else
    merge_order = round(merge_order);
end
if ~isempty(trialMarker)
    trialMarker = cellfun(@num2str, trialMarker, 'UniformOutput',false);    
    assert(isequaln(numel(trialMarker),numel(folderList)),'umIToolbox:mergeRecordings:WrongInput',...
        'The number of trial IDs must be the same as the number of input folders!');
end
    
% Append ".dat" to filenames if not done yet:
if ~endsWith(SaveFilename,'.dat')
    SaveFilename = [SaveFilename '.dat'];
end
if ~endsWith(filename,'.dat')
    filename = [filename '.dat'];
end
metaData_filename = strrep(filename, '.dat', '.mat');
%%%%%% Further input validation %%%%%%
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
metaDataPath = fullfile(folderList, metaData_filename);
dataPath = fullfile(folderList, filename);

% Check if the associated .mat file exists in the folder:
idx = cellfun(@(x) isfile(x), metaDataPath);
assert(all(idx), 'umIToolbox:mergeRecordings:MissingInput',...
    'One or more associated .mat file are missing!')
% Check if the input files are image time series with dimensions {Y,X,T}:
idxDim = false(size(dataPath));
idxSz = idxDim;
refSz = load(metaDataPath{1}, 'datSize');

for i = 1:length(dataPath)
    md = matfile(metaDataPath{i});
    idxDim(i) = all(ismember(md.dim_names, {'Y','X','T'}));
    idxSz(i) = isequaln(md.datSize,refSz.datSize);
end
assert(all(idxDim), 'umIToolbox:mergeRecordings:WrongInput',...
    'This function accepts only image time series with dimensions {"Y", "X","T"}!');
% Check if all data have the same Y,X dimensions sizes:
assert(all(idxSz), 'umIToolbox:mergeRecordings:WrongInput','All input files must have the same Y,X sizes!');
% Check if input data is valid (i.e. if it is from LabeoTech)

reqFields = {'Freq', 'datName', 'datLength', 'FirstDim', 'dim_names',...
    'datFile', 'Datatype', 'datSize'};
idx = cellfun(@(x) all(ismember(reqFields, fieldnames(matfile(x)))), metaDataPath);
assert(all(idx), 'umIToolbox:mergeRecordings:WrongInput', 'Input data must be generated from LabeoTech Imaging systems');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data merge:

% Open new .mat (meta data) file:
mOut = struct();
mIn = matfile(metaDataPath{1});
fn = setdiff(fieldnames(mIn),[{'Properties'}, reqFields]);
% Generate "Stim" fields:
stim_fn = fn(startsWith(fn, 'Stim'));
if ( b_IgnoreStim || isempty(stim_fn) || sum(mIn.Stim) == 0 )&& isempty(trialMarker)
    mOut.Stim = 0;
    b_hasStim = false;
elseif isempty(trialMarker)
    b_hasStim = true;
    for i  = 1:length(stim_fn)
        mOut.(stim_fn{i}) = [];
    end
elseif ~isempty(trialMarker)
    b_hasStim = false;
    mOut.Stim_TrialMarker = [];    
else
    error('This error should not be reached!')
end
% Get extra fields from the first file as well. Here, we assume that these
% fields are the same across all input files and fixed!
fn = setdiff(fn, fn(startsWith(fn, 'Stim')));
for i = 1:length(fn)
    mOut.(fn{i}) = mIn.(fn{i});
end
% Get required fields meta data from the first file:
mOut.Freq = mIn.Freq;
mOut.datName = mIn.datName;
mOut.FirstDim = mIn.FirstDim;
mOut.dim_names = mIn.dim_names;
mOut.datFile = SaveFilename;
mOut.Datatype = mIn.Datatype;
mOut.datLength = [];
mOut.datSize = mIn.datSize;

w = waitbar(0,'Merging data...', 'Name', ['Merging ' filename]);
if strcmpi(merge_type, 'movie')
    % Open new .dat file:
    fidOut = fopen(SaveFilename,'w');
    % Concatenate data in time domain:
    for i = 1:length(dataPath)
        % Merge data:
        fidIn = fopen(dataPath{i},'r');
        mIn = matfile(metaDataPath{i});
        data = fread(fidIn,inf,['*' mIn.Datatype]);
        fwrite(fidOut,data, mIn.Datatype);
        fclose(fidIn);
        % Update meta data:
        mOut.datLength = sum([mOut.datLength, mIn.datLength]); % Update data length         
        if b_hasStim 
            % Concatenate Stim-related variables:       
            for j = 1:length(stim_fn)
                mOut.(stim_fn{j}) = [mOut.(stim_fn{j}), mIn.(stim_fn{j})];
            end
        end
        if ~isempty(trialMarker)
            % Populate "Stim_trialMarker":
            SubStim = i*(ones(1,mIn.datLength)); SubStim([1,end]) = 0; % Use the first/last frame to mark the onset and offset of the Trial.
            mOut.Stim_TrialMarker = [mOut.Stim_TrialMarker, SubStim];
            
%             % Add Stim_TrialMarker info, if applicable:            
%             mOut.preEventTime_sec = single(1/mIn.Freq);
%             mOut.postEventTime_sec = single((mIn.datLength-1)/mIn.Freq);
        end
        waitbar(i/length(dataPath),w);
    end
    fclose(fidOut);   
    % Create info for "events.mat" file:
    
else
    % Concatenate data in "E"vent domain:
    % Create empty 4D matrix with dimensions {Trials,Y,X,T}:
    data = nan([length(dataPath), mOut.datSize, mIn.datLength], 'single');
    % Populate each trial with the input data
    for i = 1:length(dataPath)
        [datIn, mIn] = loadDatFile(dataPath{i});
        data(i,:,:,:) = datIn(:,:,1:size(data,4));        
        if b_hasStim
            % Concatenate Stim-related variables:       
            for j = 1:length(stim_fn)
                mOut.(stim_fn{j}) = [mOut.(stim_fn{j}), mIn.(stim_fn{j})];
            end
        end  
        if ~isempty(trialMarker)
            % Populate "Stim_trialMarker":
            SubStim = ones(1,mIn.datLength); SubStim([1,end]) = 0; % Use the first/last frame to mark the onset and offset of the Trial.
            mOut.Stim_TrialMarker = [mOut.Stim_TrialMarker, SubStim];
        end
        waitbar(i/length(dataPath),w);
    end
    % Update dimensions in metadata:
    mOut.dim_names = [{'E'}, mOut.dim_names];
    mOut.datSize = [length(dataPath) mIn.datSize(1)];
    mOut.datLength = [mIn.datSize(2) mIn.datLength(1)];
    % Add event info to meta data:
    if isempty(trialMarker)
        mOut.eventID = ones(mOut.datSize(1),1,'uint16');
        mOut.eventNameList = {'1'};
        mOut.preEventTime_sec = single(1/mOut.Freq);
        mOut.postEventTime_sec = single((mOut.datLength(2)-1)/mOut.Freq);
    end        
    waitbar(1,w,'Please wait. Writing data to file...')
    save2Dat(SaveFilename, data, mOut);
end
% Create "events.mat" file from "Stim" parameters:
genEventsFromStimParameters(fileparts(SaveFilename), mOut);
if ~isempty(trialMarker)
    % Open "events.mat" file and update eventNameList var:
    evntFile = matfile(fullfile(fileparts(SaveFilename), 'events.mat'), 'Writable', true);    
    evntFile.eventNameList = trialMarker;
    
end
% Copy AcqInfo file from one of the original files to get some experiment info. This is used by some IOI_ana functions.
copyfile(fullfile(fileparts(dataPath{end}), 'AcqInfos.mat'), fullfile(fileparts(SaveFilename),'AcqInfos.mat'));
close(w)
% Add data history to meta data file with info of this function:
myInfo = dir([mfilename('fullpath') '.m']);
opts = struct;
opts.filename = filename;
opts.folderList = folderList;
opts.merge_type = merge_type;
opts.merge_order = merge_order;
opts.b_IgnoreStim = b_IgnoreStim;
dH = genDataHistory(myInfo,[SaveFilename ' = mergeRecordings(opts);'], opts, {SaveFilename});
mOut.dataHistory = dH;
% Save meta data to file:
save(strrep(SaveFilename, '.dat', '.mat'), '-struct', 'mOut');   
disp('Done')
end



