function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
% RUN_IMAGESCLASSIFICATION calls the function
% IMAGESCLASSIFICATION from the IOI library (LabeoTech).

% Defaults:
default_Output = {'fluo_475.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; % This is here only as a reference for PIPELINEMANAGER.m . The real outputs will be stored in OUTFILE.
%%% Arguments parsing and validation %%%
p = inputParser;
% Raw folder:
addRequired(p, 'RawFolder', @isfolder);
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_SubROI', false, 'b_IgnoreStim', true);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));

parse(p, RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = {};
%%%%
% Get existing ImagesClassification files in directory:
existing_ChanList  = dir(fullfile(SaveFolder,'*.dat'));
idxName = ismember({existing_ChanList.name}, default_Output);
existing_ChanList = existing_ChanList(idxName);
% Calls function from IOI library. Temporary for now.
ImagesClassification(RawFolder, SaveFolder, opts.BinningSpatial, opts.BinningTemp,...
    opts.b_IgnoreStim, opts.b_SubROI);
% Get only new files created during ImagesClassification:
chanList = dir(fullfile(SaveFolder,'*.dat'));
idx = ismember({chanList.name}, default_Output);
chanList = chanList(idx);
idxName = ismember({chanList.name}, {existing_ChanList.name});
idxDate = ismember([chanList.datenum], [existing_ChanList.datenum]);
idxNew = ~all([idxName; idxDate],1);
chanList = {chanList(idxNew).name};
%%%%%% TEMPORARY FIX UNTIL CHANGES TO IMAGESCLASSIFICATION ARE APPLIED %%%%
for i = 1:numel(chanList)
    disp(['Flipping X and Y dimensions of ' chanList{i}])
    [dataMap, metaMap]= mapDatFile(fullfile(SaveFolder,chanList{i}));
    dataMap.Writable = true;
    dataMap.Data.data = permute(dataMap.Data.data,[2 1 3]);
    metaMap.Properties.Writable = true;
    metaMap.dim_names = {'Y','X','T'};
end

% If there is Stimulation, add "eventID" and "eventNameList" to the output
% files of ImagesClassification.
if ~opts.b_IgnoreStim
    disp('Creating events file.');    
    % Here the first channel fom "chanList" is chosen to retrieve the
    % "Stim" data:
    chan = matfile(fullfile(SaveFolder, strrep(chanList{1}, '.dat', '.mat'))); 
    [eventID, state, timestamps] = getEventsFromTTL(chan.Stim, chan.Freq, .5);    
    eventNameList = {num2str(unique(eventID))};    
    if isempty(eventID)
        warning('Stim signal not found! Skipped Event file creation.')
    else 
        saveEventsFile(SaveFolder, eventID, timestamps, state, eventNameList)
    end
end
outFile = fullfile(SaveFolder, chanList);
end