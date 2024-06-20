function [dataByEv, condIndx,repIndx] = reshape_data_by_event(data,SaveFolder)
% RESHAPE_DATA_BY_EVENT uses the information in the "events.mat" file to
% reshape a time series dataset (with dimensions Y,X,T) into a 4D array
% with dimensions E(vent), Y,X,T. 
% Inputs:
%   data (3D num array): Image time series to be reshaped.
%   SaveFolder (char): Path to folder containing the "events.mat" file.
% Outputs:
%   dataByEv (4D num array): Image time series split by events. Note that
%       the "selectedEvents" variable from the "events.mat" file is used to
%       output only selected conditions and/or repetitions. (See DataViewer
%       documentation for more information on how to ignore specific
%       conditions/repetitions).
%   condIndx (num vector): Vector of event index (from "eventNameList" 
%       variable in "events.mat" file).
%   repindx (num vector): Vector of repetition indices for each condition.

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p, 'SaveFolder', @isfolder);
parse(p,data,SaveFolder);
%%%

% Instantiate an events manager object:
evObj = EventsManager(SaveFolder);
% Load the events.mat info:
evObj.loadEvents(SaveFolder);
% Create the 4D array with dimensions E,Y,X,T:
nameIndx = find(any(evObj.selectedEvents,2));
condList = evObj.eventNameList(nameIndx);
frMat = [];
repIndx = [];
condIndx = [];
for ii = 1:length(condList)
    [frm,rpi] = evObj.createTrialFrameMatrix(condList{ii});
    frMat = [frMat;frm];
    repIndx = [repIndx;rpi];
    condIndx = [condIndx; repmat(nameIndx(ii),length(rpi),1)];
end
clear frm rpi
% Split data by events    
b_validFrames = ~isnan(frMat);    
% ! Missing frames will be set to NaN
dataByEv = nan(size(frMat,1),size(data,1),size(data,2),size(frMat,2),'single');
for ii = 1:size(frMat,1)
    dataByEv(ii,:,:,b_validFrames(ii,:)) = data(:,:,frMat(ii,b_validFrames(ii,:)));
end

end
