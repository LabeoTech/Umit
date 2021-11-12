function metaData = genMetaData(data, dim_names, varargin)
% This function creates a strucure containing the metadata of "data" in
% order to work with PipelineManager class.
% Inputs:
%   data (numerical array) : Data which we want to retrieve the metaData.
%   dim_names (cell array of char) : names of dimensions from "data";
%    Examples of dimension names are listed below:
%       "O" = observation;
%       "X" = x axis;
%       "Y" = y axis;
%       "Z" = z axis;
%       "T" = time;
%       "E" = events; 
%   extraParams(struct) Optional: structure containing extra fields to add
%       to "metaData" or existing fields that you wish to overwrite. 
    
%   Here are the default fields and values for a metaData file:
%       datName(str) : 'data';
%       dim_names(cell array of char): dim_names;
%       datFile(str) : '';
%       datSize(double) 1x2 : size of 1st and 2nd dimension of "data";
%       datLength(double) 1xN : size of the 3rd + dimensions of "data";
%       Datatype(str) : class of "data";
%       
% Output:
%   metaData (struct): structure containing the minimal variables to work
%   with the umIToolbox.

p = inputParser;
addRequired(p, 'data', @isnumeric);
addRequired(p, 'dim_names', @iscell);
addOptional(p, 'extraParams', struct, @isstruct);
parse(p, data, dim_names, varargin{:});
% Further validate dim_names:
root = getenv('Umitoolbox');
dim_names_info = load(fullfile(root, 'subFunc','dimension_names.mat'));
errID = 'Umitoolbox:save2Dat:InvalidName';
errMsg = 'List of dimension names contain invalid values.';
assert(all(ismember(p.Results.dim_names, dim_names_info.dims_dict)), errID, errMsg);

% Extract info from data
metaData = struct;
metaData.datName = 'data';
metaData.dim_names = p.Results.dim_names;
metaData.datFile = '';
ds = size(p.Results.data);
metaData.datSize = ds([1 2]);
metaData.datLength = ds(3:end); % Accounts for 3+ dimensions.
metaData.Datatype = class(p.Results.data);

% Merge extrParams struct with metaData struct:
fNames = fieldnames(p.Results.extraParams);
for i = 1:numel(fNames)
    metaData.(fNames{i}) = p.Results.extraParams.(fNames{i});
end
end