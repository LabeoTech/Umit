function outFile = event_triggered_average(File, SaveFolder, varargin)
% EVENT_TRIGGERED_AVERAGE calculates the average across repetitions of events 
% listed in  EVENTS.MAT inside the path of FILE.
%
% Inputs:
%   File: fullpath of functional imaging .DAT file.
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
%   opts (optional) : structure containing extra parameters.
% Output:
%   outFile: name of Output file.
%   .DAT file containing average (AVG) and standard deviation (STD) movies
%   per event with dimensions (X,Y,T,Event).

% Defaults:
default_opts = struct('preEventTime_sec',2, 'postEventTime_sec',4, 'PadWith', 'mean');
default_Output = 'eventTrigAVG.dat'; 

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'File',@isfile)
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ismember(x.PadWith, {'mean', 'NaN', 'zero'})); % Padding options for cases where movie snippets dont have the same length.
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%
% Map movie and metadata to memory:
[folder,~,~] = fileparts(File); % Assuming that File is in the same folder as events.mat.
mData = mapDatFile(File);
metaDat = matfile(strrep(File, '.dat', '_info.mat'));
szdat = size(mData.Data.data);
evDat = load(fullfile(folder, 'events.mat'));
% Create empty matrix:
sr = metaDat.Freq;
preFr = round(sr*opts.preEventTime_sec);
postFr = round(sr*opts.postEventTime_sec);
centralFr = preFr + 1;
len_trial = preFr + postFr;
n_reps = (sum(evDat.state == 1))/numel(evDat.eventNameList);
n_trial = numel(evDat.eventNameList);
data = nan(n_reps, n_trial, szdat(1), szdat(2), len_trial, 'single'); % (repetition, trialIndex, X, Y, T).
% Fill empty matrix with data segments
for i = 1:n_trial
    repFr = round(sr.*evDat.timestamps(evDat.eventID == i & evDat.state == 1));
    for j = 1:numel(repFr)
        start = repFr(j) - preFr;
        stop = repFr(j) + postFr - 1;
        if start < 1
            start = szdat(3);
        elseif stop > szdat(3)
            stop = szdat(3);
        end
        snippet = mData.Data.data(:,:,start:stop);
        startFr = centralFr - (repFr(j) - start);
       	stopFr = centralFr + (stop - repFr(j));
        data(j,i,:,:,startFr:stopFr) = snippet;
    end
end
% Check for NaNs and replace values using method specified by opts.PadWith:
idx = isnan(data);
if any(idx,'all')
    switch opts.PadWith
        case 'zero'
            data(idx) = 0;
        otherwise
            % The following lines are highly inefficient. To be improved.
            % Fill NaNs with the mean value of the trial
            [rep,trial,x,y,time] = ind2sub(size(data),find(idx));
            data(rep,trial,x,y,time) = nanmean(data(rep,trial,x,y,:),5);
    end
end
% Calculate average and Standard deviation movies.
AVG = nanmean(data,1); 
AVG = permute(AVG,[3 4 5 2 1]); % X,Y,T,Event
STD = std(data,0,1,'omitnan');
STD = permute(STD,[3 4 5 2 1]); % X,Y,T,Event
%
datFile = fullfile(SaveFolder, Output);
% Create MetaData structure:
szAVG = size(AVG);
metaDat = struct('datName', {'AVG', 'STD'}, 'datSize', {szAVG([1 2]), szAVG([1 2])},...
    'datLength', {szAVG(3:end) , szAVG(3:end)}, 'Datatype', {'single', 'single'}, 'datFile', datFile);
metaDat(1).eventList = evDat.eventNameList; % Added event description from events.mat file to the metadata... Not sure if I'll keep this.
metaDat(1).eventFrame = centralFr;
% Save AVG and METADAT to DATFILE:
save2Dat(datFile, AVG,'-w', metaDat)
% Append "STD" to DATFILE:
save2Dat(datFile, STD, '-a')
% Output file names
outFile = Output;
end



