function save2Mat(MatFileName, data, labels, dim_names)
% SAVE2MAT validates all variables necessary for the statistical module
% and saves DATA and metadata variables to MATFILENAME.
%
% Inputs:
%
% MatFileName = full path of .MAT filename
%
% data = 1 to 5-D matrix containing single-precision numerical data 
% to statistical analysis. 
%
% labels = 1+ cell array of characters containing the description of each
% observation. Must have the same length of the dimension "O" in "data".
%
% dim_names = cell array of characters with the dimension description. See
% documentation on how to label dimension names. Examples of dimension
% names are listed below:
% "O" = observation;
% "X" = x axis;
% "Y" = y axis;
% "Z" = z axis;
% "T" = time;

% Arguments validation
p = inputParser;
addRequired(p,'MatFileName',...
    @(x) isfolder(fileparts(x)) && endsWith(x,'.mat')); % Validates if MatFileName is a fullpath and .MAT.
validData = @(x) isa(x, 'single') & ~isempty(x)  & ndims(x) <= 5; % Data validation function.
addRequired(p, 'data', validData);
addRequired(p, 'labels');
addRequired(p, 'dim_names');
parse(p, MatFileName, data, labels, dim_names);
mFile = p.Results.MatFileName;
data = p.Results.data;
labels = p.Results.labels;
dim_names = upper(p.Results.dim_names);

% Validate dim_names:
root = getenv('Umitoolbox');
dim_names_info = load(fullfile(root, 'subFunc','dimension_names.mat'));

errID = 'Umitoolbox:InvalidName';
errMsg = 'Invalid dimension name. List of dimension names contain invalid values.';
assert(all(ismember(dim_names, dim_names_info.dims_dict)), errID, errMsg);

errID = 'Umitoolbox:InvalidSize';
errMsg = 'Invalid length of dimensions. The length of dimension description list do not match the size of Data.';
assert(isequal(length(dim_names), ndims(data)), errID, errMsg);

errID = 'Umitoolbox:InvalidInput';
errMsg = {'Invalid dimension names.',...
    'The number of unique dimension names is different from the length of dimension list.'};
assert(isequal(length(unique(dim_names)), length(dim_names)), errID, errMsg);

% Validate labels:
obs_idx = find(strcmp('O', dim_names));
errID = 'Umitoolbox:InvalidSize';
errMsg = 'The length of Labels do not match the length of observations in Data.';
assert(isequal(length(labels{obs_idx}),size(data,obs_idx)), errID, errMsg);
errID = 'Umitoolbox:InvalidType';
errMsg = 'Invalid Labels. Labels must be a cell array of strings/characters.';
for i = 1:numel(labels)
assert(iscell(labels) && all(cellfun(@(x) isa(x, 'char') || isa(x, 'string'),...
    labels{i})), errID, errMsg);
end
% Generate File UUID:
fileUUID = char(java.util.UUID.randomUUID);
% Save everything to .MAT file:
save(mFile, 'data', 'labels', 'dim_names', 'fileUUID');
[p,n,ext] = fileparts(mFile);
disp(['Stats data saved in ' p ' as ' [n ext]]);

end