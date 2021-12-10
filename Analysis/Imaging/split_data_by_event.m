function outFile = split_data_by_event(File, SaveFolder, varargin)
% SPLIT_DATA_BY_EVENT reshapes an image time series dataset in a 4D matrix of 
% dimensions  Event (E), X,Y,T.%
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
    ismember(x.PadWith, {'mean', 'NaN'}) || isnumeric(x)); % Padding options for cases where movie snippets dont have the same length.
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
[mData, metaDat] = mapDatFile(File);
szdat = size(mData.Data.data);
evDat = load(fullfile(folder, 'events.mat'));
% Create empty matrix:
sr = metaDat.Freq;
preFr = round(sr*opts.preEventTime_sec);
postFr = round(sr*opts.postEventTime_sec);
centralFr = preFr + 1;
len_trial = preFr + postFr;

n_trial = sum(evDat.state == 1);
timestamps = evDat.timestamps(evDat.state == 1);
new_dims = {'E', 'Y', 'X','T'};
[~, locB]= ismember(new_dims([2,3]), metaDat.dim_names);
data = nan([n_trial, szdat(locB), len_trial], 'single');
% Fill empty matrix with data segments
fix_snippet = false;
for i = 1:n_trial
    trialFr = round(sr*timestamps(i));
    start = trialFr - preFr;
    stop = trialFr + postFr - 1;
    if start < 1
        start = 1;
        fix_snippet = true;
    elseif stop > szdat(3)
        stop = szdat(3);
        fix_snippet = true;
    end
    snippet = mData.Data.data(:,:,start:stop);
    startFr = centralFr - (trialFr - start);
    stopFr = centralFr + (stop - trialFr);
    if fix_snippet
        warning(['Snippet size is out of bounds from Data.'...
            ' Missing data points will be replaced with ' opts.PadWith]);
        switch opts.PadWith
            case 'mean'
                avg = nanmean(snippet, 3);
                avg = repmat(avg,1,1,size(data,4));
                data(i,:,:,:) = avg;
            case 'NaN'
                % empty
            otherwise
                data(i,:,:,:) = opts.PadWith;
        end
    end
    data(i,:,:,startFr:stopFr) = snippet;
end

datFile = fullfile(SaveFolder, Output);
% Save DATA and METADATA to DATFILE:
save2Dat(datFile, data, new_dims)
[~, metaData] = mapDatFile(datFile);
metaData.Properties.Writable = true;
metaData.preEventTime_sec = opts.preEventTime_sec;
metaData.postEventTime_sec = opts.postEventTime_sec;
metaData.eventID = evDat.eventID(evDat.state == 1);
metaData.eventNameList = evDat.eventNameList;
% Output file names
outFile = Output;
end