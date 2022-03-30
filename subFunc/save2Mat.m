function varargout = save2Mat(MatFileName, data, obsID, dim_names, varargin)
% SAVE2MAT validates all variables necessary for the statistical module
% and saves DATA and metadata variables to MATFILENAME.
%
% Inputs:
%
%   MatFileName = full path of .MAT filename or an empty array. If empty,
%       the function will generate a structure as output, otherwise it will
%       save the variables to a .MAT file "MatFileName".
%   data = 1D cell array with length equal to the number of elements of
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
%   appendMetaData (optional Name-value parameter):
%       handle to .MAT (matfile) containing metaData that will be added the
%       output .MAT file.
%   genFile (bool): If TRUE, this function will save the data to a .MAT,
%       otherwise it will output it as "out".

% Output (optional):
%   out = structure containing data and metaData.

% Arguments validation
p = inputParser;
addRequired(p,'MatFileName',...
    @(x) isempty(x) || (isfolder(fileparts(x)) && endsWith(x,'.mat'))); % Validates if MatFileName is a fullpath and .MAT OR is empty.
% validData = @(x) iscell(x) & ~isempty(x) & isa(x{:}, 'single'); % Data validation function.
% addRequired(p, 'data', validData);
addRequired(p, 'data');
addRequired(p, 'obsID', @iscell);
addRequired(p, 'dim_names', @iscell);
addParameter(p, 'label', 'val', @(x) (iscell(x) && ischar([x{:}])) || (ischar(x)));
addParameter(p, 'appendMetaData', [], @(x) isempty(x) || isa(x, 'matlab.io.MatFile') ||...
    isstruct(x));
addParameter(p, 'genFile', true, @islogical)
addParameter(p,'appendObjectInfo',[], @(x) isempty(x) || isa(x,'Subject') ||....
    isa(x,'Acquisition') || isa(x,'Modality'));
parse(p, MatFileName, data, obsID, dim_names, varargin{:});
% Instantiate input variables:
mFile = p.Results.MatFileName;
data = p.Results.data;
obsID = p.Results.obsID;
label = p.Results.label;
dim_names = upper(p.Results.dim_names);
metaData = p.Results.appendMetaData;
b_genFile = p.Results.genFile;
objHandle = p.Results.appendObjectInfo;
clear p
% Further validate data:
errID = 'Umitoolbox:save2Mat:IncompatibleSize';
errMsg = 'The length of data is different from the number of observations.';
assert(isequaln(length(data),length(obsID)), errID, errMsg);
errID = 'Umitoolbox:save2Mat:IncompatibleSize';
errMsg = 'The number of dimensions of data is different from the number of dimension names.';
% Get the size of the largest dataset from "data":
max_data_size = cellfun(@(x) size(x), data, 'UniformOutput',false);
max_data_size = vertcat(max_data_size{:});
max_data_size = max(max_data_size,[],1);
% Verify if the number of non-singleton dimensions in data match the
% number of dimensions names, exept "O":
assert(isequaln(sum(max_data_size~=1),sum(~strcmp(dim_names, 'O'))), errID, errMsg);

% Further validate dim_names:
root = getenv('Umitoolbox');
dim_names_info = load(fullfile(root, 'subFunc','dimension_names.mat'));

errID = 'Umitoolbox:save2Mat:InvalidName';
errMsg = 'List of dimension names contain invalid values.';
assert(all(ismember(dim_names, dim_names_info.dims_dict)), errID, errMsg);

% Create structure with data and metadata:
if isa(metaData, 'matlab.io.MatFile')
    % If a metadata file is provided as a matfile, load it and overwrite variables:
    s = load(metaData.Properties.Source);
elseif isstruct(metaData)
    s = metaData;
end
s.data = data;
s.obsID = obsID;
s.dim_names = dim_names;

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

if ~isempty(objHandle)
    % Add experiment metaData:
    objectInfo = struct();
    while isprop(objHandle, 'MyParent')
        Field = class(objHandle);
        % Remove object's handles from structure:        
        Props = setdiff(properties(objHandle), {'MyParent', 'Array'});
        for i=1:numel(Props)
            objectInfo.(Field).(Props{i}) = objHandle.(Props{i});
        end
        objHandle = objHandle.MyParent;
    end
    s.experimentMetaData = objectInfo;
end

% Save data to file
if b_genFile
    % Add file unique identifier:
    % Save "s" struct to file:
    save(mFile, '-struct', 's', '-v7.3');
    [p,n,ext] = fileparts(mFile);
    disp(['Stats data saved in ' p ' as ' [n ext]]);
else
    varargout{1} = s;
end

end
