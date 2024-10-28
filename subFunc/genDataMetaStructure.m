function out = genDataMetaStructure(data,varargin)
% GENDATAMETASTRUCTURE validates all variables necessary for the statistical module
% and merges the DATA and METADATA into a single structure ("out").
%
% Inputs:
%   data = numeric array (e.g. Y,X,T image time series) OR struct OR cell array
%       with length equal to the number of elements of observations (obsID).
%   Optional:
%   obsID = 1D cell array of characters containing the description of each
%       observation.
%   hasEvents(bool | default= false): Set to TRUE, to indicate that the 
%       data is split by events.
%   extraInfo (struct): other meta data associated with "data".
% Output:
%   out = structure containing data and extraInfo.

% Arguments validation
p = inputParser;
addRequired(p, 'data',@(x) isstruct(x) || isnumeric(x) || iscell(x));
addOptional(p, 'obsID',{'genericObs'}, @(x) iscell(x) & ischar([x{:}]));
addParameter(p, 'hasEvents',false,@islogical);
addParameter(p, 'extraInfo',struct(), @isstruct);
parse(p, data, varargin{:});
% Instantiate input variables:
data = p.Results.data;
obsID = p.Results.obsID;
b_hasEvents = p.Results.hasEvents;
metaData = p.Results.extraInfo;
clear p
% Put data in structure if it is numeric:
if isnumeric(data)
    data = struct('data',data);
elseif iscell(data)
    data = cell2struct(data,'data');
end
% Check if the data length is the same as the number of observation IDs:
errID = 'Umitoolbox:genDataMetaStructure:IncompatibleSize';
errMsg = 'The length of data is different from the number of observations.';
assert(isequaln(length(data),length(obsID)), errID, errMsg);
% Categorize the data based on its dimensions:
fn = fieldnames(data);
dataCategory = struct();
for ii = 1:length(fn)
    this_data = data(1).(fn{ii});
    if b_hasEvents
        % Remove event dimension to classify the data:
        eval(['this_data = squeeze(this_data(1' repmat(',:',1,ndims(this_data)-1) '));'])
    end
    dataDims = ndims(this_data);
    if isscalar(this_data)
        dataCat = 'scalar';
    elseif isvector(this_data)
        dataCat = 'time-vector';
    elseif dataDims == 2
        dataCat = 'Map';
    elseif dataDims == 3
        dataCat = 'Image-time-series';
    else
        dataCat = 'unknown';
    end
    dataCategory.(fn{ii}) = dataCat;
end      
% Set output structure:
out = metaData;
out.obsID = obsID;
out.data = data;
out.b_hasEvents = b_hasEvents;
out.b_hasMultipleMeasures = length(fn)>1;
out.dataCategory = dataCategory;
end
