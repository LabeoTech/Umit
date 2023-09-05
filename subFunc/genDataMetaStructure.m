function out = genDataMetaStructure(data, obsID, dim_names, metaData, varargin)
% GENDATAMETASTRUCTURE validates all variables necessary for the statistical module
% and merges the DATA and METADATA into a single structure ("out").
%
% Inputs:
%   data = cell array OR struct with length equal to the number of elements of
%       observations "O".
%   obsID = 1D cell array of characters containing the description of each
%       observation.
%   dim_names = cell array of characters with the dimension description. See
%       documentation on how to label dimension names. Examples of dimension
%       names are listed below:
%           "O" = observation;
%           "X" = x axis;
%           "Y" = y axis;
%           "Z" = z axis;
%           "T" = time;
%           "E" = events; *In this case, the function will validate if event IDs and
%           Names are contained in "data" file.
%   metaData (struct): meta data associated with "data".
%   label (cell | ): array of characters
% Output:
%   out = structure containing data and metaData.

% Arguments validation
p = inputParser;
addRequired(p, 'data',@(x) iscell(x) || isstruct(x));
addRequired(p, 'obsID', @iscell);
addRequired(p, 'dim_names', @iscell);
addRequired(p, 'metaData', @isstruct);
addParameter(p, 'label', 'val', @(x) (iscell(x) && ischar([x{:}])) || (ischar(x)));
parse(p, data, obsID, dim_names, metaData, varargin{:});
% Instantiate input variables:
data = p.Results.data;
obsID = p.Results.obsID;
dim_names = upper(p.Results.dim_names);
metaData = p.Results.metaData;
label = p.Results.label;
clear p

% Further validate data:
errID = 'Umitoolbox:genDataMetaStructure:IncompatibleSize';
errMsg = 'The length of data is different from the number of observations.';
assert(isequaln(length(data),length(obsID)), errID, errMsg);
% Get the size of the largest dataset from "data":
if isstruct(data)
    % If the data has multiple measures:
    fn = fieldnames(data);
    max_data_size = arrayfun(@(x) size(x.(fn{1})), data,'UniformOutput',false);
else
    max_data_size = cellfun(@(x) size(x), data, 'UniformOutput',false);
end
max_data_size = vertcat(max_data_size{:});
max_data_size = max(max_data_size,[],1);
% Verify if the number of non-singleton dimensions in data match the
% number of dimensions names, exept "O":
errID = 'Umitoolbox:genDataMetaStructure:IncompatibleSize';
errMsg = 'The number of dimensions of data is different from the number of dimension names.';
assert(isequaln(sum(max_data_size~=1),length(dim_names(2:end))), errID, errMsg);
% Further validate dim_names:
root = fileparts(mfilename('fullpath'));
dim_names_info = load(fullfile(root, 'dimension_names.mat'));
%
errID = 'Umitoolbox:genDataMetaStructure:InvalidName';
errMsg = 'List of dimension names contain invalid values.';
assert(all(ismember(dim_names, dim_names_info.dims_dict)), errID, errMsg);
% Merge data and meta data:
out = metaData;
out.data = data;
% Add fields:
out.obsID = obsID;
out.dim_names = dim_names;
% Set dimension of data that will be labeled. Default = 2;
dim_label = 2;
if ismember('E', dim_names)
    % Create labels based in event names:        
    % Find "Events" dimension in data:
    idx_dim = zeros(length(out.data),ndims(out.data{1}));
    for ii = 1:length(out.data)
        idx_dim(ii,:) = (size(out.data{ii}) == length(out.eventID));
    end
    dim_label = find(all(idx_dim,1));
    fn = fieldnames(out);
    errID = 'Umitoolbox:genDataMetaStructure:MissingInfo';
    errMsg = 'An event dimension name ("E") was detected but no event info ("eventID" and "eventNameList") was found in metaData.';
    assert(all(ismember({'eventID', 'eventNameList'}, fn)), errID, errMsg);
    disp('Creating labels for Events...')
    % Overwrite "label" with the event names instead:
    label = cell(1,length(out.eventID));
    for i = 1:numel(out.eventNameList)
        indx = find(out.eventID == i);
        for j = 1:numel(indx)
            label{indx(j)} = strjoin({out.eventNameList{i}, 'rep', num2str(j)}, '_');
        end
    end
elseif strcmpi(label,'val')
    % Use generic labels:
    label = arrayfun(@(x) strjoin({label, num2str(x)}, '_'), 1:size(out.data{1},dim_label), 'UniformOutput', 0);
end
% Check if "label" has the same length of data:
errID = 'Umitoolbox:genDataMetaStructure:IncompatibleSize';
errMsg = 'The length of "labels" is different from the length of "data".';
assert(isequaln(size(out.data{1},dim_label),length(label)), errID, errMsg);
% Add "label":
out.label = label;
end
