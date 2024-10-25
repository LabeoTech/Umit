function outData = split_data_by_event(data,SaveFolder)
% SPLIT_DATA_BY_EVENT uses the information in the "events.mat" file to
% reshape 3D arrays (Y,X,T) into 4D array with dimensions E,Y,X,T.
% The data must be an Image time series with dimensions
% {Y,X,T}.
%
% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T").
%   SaveFolder (chr): full path to folder containing the "events.mat" file.
%
% Outputs: 
%   outData: structure containing stats-ready data extracted from ROIs.
%
% Defaults:
default_Output = 'dataByEv.mat'; %#ok This line is here just for Pipeline management.

% Arguments validation:
assert(isnumeric(data) & ndims(data) == 3, 'Data has to be a 3D array with dimensions Y,X,T!'); % Validate if the input is a 3-D numerical matrix
assert(isfolder(SaveFolder),'The folder "%s" doesn''t exist!',SaveFolder);
% Instantiate EventsManager:
ev = EventsManager(SaveFolder);
% Load the info in the 'events.mat' file:
ev.loadEvents(SaveFolder);
% Split the data by events: Here, the trials will be cropped to the
% shortest length for uniformity.
dataByEv = ev.splitDataByEvents(data);
% Package the data in structure and add event info:
outData = genDataMetaStructure(dataByEv,'hasEvents',true,'extraInfo',ev.exportEventInfo);

end
