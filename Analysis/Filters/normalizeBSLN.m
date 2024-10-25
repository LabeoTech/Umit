function outData = normalizeBSLN(data, SaveFolder, varargin)
% NORMALIZEBSLN normalizes event-triggered and image time series data by a baseline.
% Here, a baseline corresponds to the median of the pixel values before the event
% (pre-event time) or a value in seconds stored in the variable
% "baseline_sec" for each trial.

% Inputs:
%   data (3D numerical matrix): Image time series data.
%   SaveFolder(char): Folder containing the "events.mat" file.
%   opts (optional): structure containing extra parameters:
%       - b_centerAtOne (bool): Set to TRUE to center the normalized data at
%           one. Otherwise, the data will be centered at zero.
% Output:
%   outData(3D numerical matrix): "data" with values transformed to express
%       the normalized values by the baseline.

% Defaults:
default_Output = 'normBSLN.dat'; %#ok This is here only as a reference for PIPELINEMANAGER.m. 
default_opts = struct('b_centerAtOne', false);
opts_values = struct('b_centerAtOne', [true,false]);%#ok. This is here only as a reference for PIPELINEMANAGER.m. 
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a numerical 3D matrix.
addRequired(p,'SaveFolder',@(x) isfolder(x)); % Validate if the SaveFolder exists.
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x)); 
% Parse inputs:
parse(p,data, SaveFolder,varargin{:});
% Initialize Variables:
opts = p.Results.opts;
clear p
%%%%

% Check and load the "events.mat" file:
assert(isfile(fullfile(SaveFolder, 'events.mat')),'The "events.mat" file was not found in "%s"!',SaveFolder);
% Create EventsManager object and load the event info from the "events.mat" file:
ev = EventsManager(SaveFolder);ev.loadEvents(SaveFolder);
% Set output array:
outData = zeros(size(data),'single');
% 
bslnFrames = 1:round(ev.baselinePeriod*ev.AcqInfo.FrameRateHz);
for ii = 1:length(ev.eventNameList)
    frIndexMat = ev.getFrameMatrix(ev.eventNameList{ii});
    for jj = 1:size(frIndexMat,1)
        bslnTrialFrames = frIndexMat(jj,bslnFrames);bslnTrialFrames(isnan(bslnTrialFrames)) = [];
        allTrialFrames = frIndexMat(jj,~isnan(frIndexMat(jj,:)));
        % Calculate baseline (median from pre-event time):
        bsln = median(data(:,:,bslnTrialFrames),3,'omitnan');
        % Normalize by baseline:
        outData(:,:,allTrialFrames) = (data(:,:,allTrialFrames) - bsln)./bsln;
    end     
end

% Center data at ONE:
if opts.b_centerAtOne
    outData = outData + 1;
end
end