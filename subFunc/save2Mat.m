function save2Mat(MatFileName, data, obsID, dim_names, varargin)
% SAVE2MAT validates all variables necessary for the statistical module
% and saves DATA and metadata variables to MATFILENAME.
%
% Inputs:
%
% MatFileName = full path of .MAT filename
%
% data = 1D cell array with length equal to the number of elements of 
% observations "O".
%
% obsID = 1D cell array of characters containing the description of each
% observation.
%
% dim_names = cell array of characters with the dimension description. See
% documentation on how to label dimension names. Examples of dimension
% names are listed below:
% "O" = observation;
% "X" = x axis;
% "Y" = y axis;
% "Z" = z axis;
% "T" = time;
% "E" = events; *In this case, the function will validate if event IDs and
% Names are contained in "data" file.

% appendMetaData (optional Name-value parameter):
%  handle to .MAT (matfile) containing metaData that will be added the
%  output .MAT file.
% Arguments validation
p = inputParser;
addRequired(p,'MatFileName',...
    @(x) isfolder(fileparts(x)) && endsWith(x,'.mat')); % Validates if MatFileName is a fullpath and .MAT.
% validData = @(x) iscell(x) & ~isempty(x) & isa(x{:}, 'single'); % Data validation function.
% addRequired(p, 'data', validData);
addRequired(p, 'data');
addRequired(p, 'obsID', @iscell);
addRequired(p, 'dim_names', @iscell);
addParameter(p, 'label', 'val', @(x) (iscell(x) && ischar([x{:}])) || (ischar(x)));
addParameter(p, 'appendMetaData', [], @(x) isempty(x) || isa(x, 'matlab.io.MatFile') || isstruct(x));
parse(p, MatFileName, data, obsID, dim_names, varargin{:});
% Instantiate input variables:
mFile = p.Results.MatFileName;
data = p.Results.data;
obsID = p.Results.obsID;
label = p.Results.label;
dim_names = upper(p.Results.dim_names);
metaData = p.Results.appendMetaData;

% Further validate data:
errID = 'Umitoolbox:save2Mat:IncompatibleSize';
errMsg = 'The lenght of data is different from the number of observations.';
assert(isequaln(length(data),length(obsID)), errID, errMsg);
errID = 'Umitoolbox:save2Mat:IncompatibleSize';
errMsg = 'The number of dimensions of data is different from the number of dimension names.';
assert(isequaln(ndims(data{1}),numel(dim_names)), errID, errMsg);

% Further validate dim_names:
root = getenv('Umitoolbox');
dim_names_info = load(fullfile(root, 'subFunc','dimension_names.mat'));

errID = 'Umitoolbox:save2Mat:InvalidName';
errMsg = 'List of dimension names contain invalid values.';
assert(all(ismember(dim_names, dim_names_info.dims_dict)), errID, errMsg);

% Generate File UUID:
fileUUID = char(java.util.UUID.randomUUID);

% Create structure with data and metadata:
if ~isempty(metaData)
    % If a metadata file is provided, load it and overwrite variables:
    s = load(metaData.Properties.Source);
end
s.data = data;
s.obsID = obsID;
s.dim_names = dim_names;
s.fileUUID = fileUUID;

% Check if "E" exists in dim_names and verify if event info exists in "s"
% struct:
if ismember('E', dim_names)
    fn = fieldnames(s);
    errID = 'Umitoolbox:save2Mat:MissingInfo';
    errMsg = 'An event dimension name ("E") was detected but no event info ("eventID" and "eventNameList")was found in metaData.';
    assert(all(ismember({'eventID', 'eventNameList'}, fn)), errID, errMsg);
    disp('Creating labels for Events...')
    % Build Event Labels:
    label = cell(1,length(s.eventID));
    for i = 1:numel(s.eventNameList)
        indx = find(s.eventID == i);
        for j = 1:numel(indx)
            label{indx(j)} = strjoin({s.eventNameList{i}, 'rep', num2str(j)}, '_');
        end
    end
end
% 
disp('Prepping label...');
% Prepare "label" to save:
if ischar(label)
    label = arrayfun(@(x) strjoin({label, num2str(x)}, '_'), 1:size(s.data{1},2), 'UniformOutput', 0);
end
% Check if "label" has the same length of data:
errID = 'Umitoolbox:save2Mat:IncompatibleSize';
errMsg = 'The lenght of Labels is different from the length of data.';
assert(isequaln(size(s.data{1},2),length(label)), errID, errMsg);
% Add "label" to s:
s.label = label;
% Save "s" struct to file: 
save(mFile, '-struct', 's', '-v7.3');
[p,n,ext] = fileparts(mFile);
disp(['Stats data saved in ' p ' as ' [n ext]]);

end