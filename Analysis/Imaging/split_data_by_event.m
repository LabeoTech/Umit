function outFile = split_data_by_event(File, SaveFolder, varargin)
% SPLIT_DATA_BY_EVENT reshapes an image time series dataset in a 5D matrix of 
% dimensions  repetition (R), condition (C) , X,Y,T.%
% Inputs:
%   File: fullpath of functional imaging .DAT file.
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
%   opts (optional) : structure containing extra parameters.
% Output:
%   outFile: name of Output file.

% Defaults:
default_opts = struct('preEventTime_sec',2, 'postEventTime_sec',4, 'PadWith', 'mean');
default_Output = 'data_splitByEvent.dat'; 
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
data = nan(n_reps, n_trial, szdat(1), szdat(2), len_trial, 'single'); % (repetition, condition, X, Y, T).
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
datFile = fullfile(SaveFolder, Output);
% Save AVG and METADAT to DATFILE:
save2Dat(datFile, data)
metaData = matfile(strrep(datFile, '.dat', '_info.mat'));
metaData.Properties.Writable = true;
metaData.preEventTime_sec = opts.preEventTime_sec;
metaData.postEventTime_sec = opts.postEventTime_sec;
% Output file names
outFile = Output;
end



