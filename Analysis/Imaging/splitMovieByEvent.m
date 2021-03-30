function outFile = splitMovieByEvent(File, SaveFolder, varargin)
% SPLITMOVIEBYEVENT splits the imaging data in FILE by events stored in
% EVENTS.MAT inside the path of FILE.
%
% Inputs:
%   File: fullpath of functional imaging .DAT file.
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
%   opts (optional) : structure containing extra parameters.
% Output:
%   outFile: name of Output file.

% Defaults:
default_opts = struct('preEventTime_sec',2, 'postEventTime_sec',4, 'PadWith', 'mean');
default_Output = 'movSplitByEvent.dat'; 

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'File',@isfile)
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ismember(x.PadWith, {'mean', 'NaN', 'zero'})); % Padding options for cases where movie snippets dont have the same length.
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = opts.Output;
%%%%
% Map movie and metadata to memory:
[folder,~,~] = fileparts(File); % Assume that File is in the same folder as events.mat.
mData = mapDatFile(File);
metaDat = matfile(strrep(File, '.dat', '_info.mat'));
szdat = size(mData.Data.data);
evDat = load(fullfile(folder, 'events.mat'));
% Create empty matrix:
sr = metaDat.Freq;
preFr = round(sr*opts.preEventTime_sec);
postFr = round(sr*opts.postEventTime_sec);
len_trial = preFr + postFr;
n_reps = (sum(evDat.state == 1))/numel(evDat.eventNameList);
n_trial = numel(evDat.eventNameList);
Mdat = zeros(n_reps, n_trial, szdat(1), szdat(2), len_trial, 'single'); % (nReps, nTrials, X, Y, T).
% Populate Mdat with data snippets.

for i = 1:n_trial
    repFr = round(sr.*timestamps(evDat.eventID == i & evDat.state == 1));
    for j = 1:numel(repFr)
        start = repFr(j) - preFr;
        stop = repFr(j) + postFr - 1;
        if start < szdat(3)
            start = szdat(3);
            pad = -1;
        elseif stop > szdat(3)
            stop = szdat(3);
            pad = 1;
        else
            pad = 0;
        end
        snippet = mData.Data.data(:,:,start:stop);
            
               
    end
    
    
end







for indS = 1:nStim
    RepIdx = find(all(TrialList==StimList(indS,:),2));
    for indR = 1:nReps
        try
            pretm = sort(idx(RepIdx(indR)) - (1:(t_interStim*Infos.Freq)));
            postm = idx(RepIdx(indR)) + (0:(t_stim*Infos.Freq - 1));
            DeltaF_F = (dat(:,:,[pretm postm]) - mean(dat(:,:,pretm),3))./mean(dat(:,:,pretm),3);
            Mdat(indR, indS, :, :, :) = DeltaF_F;
        catch
            %Mdat(indR, indS, :, :, :) = nan(sdat(1),sdat(2),TrialLength*Infos.Freq);
            disp('Length Error in:');
            disp(['Repetition : ' num2str(indR)]);
            disp(['Condition: ' num2str(indS)]);            
        end
    end
end


















%Save data using save2dat.m function
save2Dat(datFile, data);
% Output file names
outFile = Output;
end



