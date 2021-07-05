function save2Dat(DatFileName, data, dim_names, varargin)
% SAVE2DAT saves data to binary .dat files. If it is a new file, it creates
% a .MAT file containing the data's METADATA. FLAG variable indicates if it
% appends ('-A', '-APPEND') to an existing file.
% Inputs:
% DatFileName = fullpath of .DAT file.
% data = non-empty "single" numeric multi-dimensional matrix.
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

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) validateattributes(x, {'single'}, {'nonempty'}));
addRequired(p, 'dim_names', @iscell);
addOptional(p, 'flag', '-w', @(x) ischar(x) || isstring(x) && mustBeMember(x, {'-w', '-write', '-a', '-append'}));
addOptional(p, 'metaData', struct(), @isstruct);
parse(p, DatFileName, data, dim_names, varargin{:});
%%%%%%
DatFileName = p.Results.DatFileName;
data = p.Results.data;
dim_names = upper(p.Results.dim_names);
flag = p.Results.flag;
metaDat = p.Results.metaData;
metaDatFilename = strrep(DatFileName, '.dat', '_info.mat');
% Further validate dim_names:
root = getenv('Umitoolbox');
dim_names_info = load(fullfile(root, 'subFunc','dimension_names.mat'));
errID = 'Umitoolbox:save2Dat:InvalidName';
errMsg = 'List of dimension names contain invalid values.';
assert(all(ismember(dim_names, dim_names_info.dims_dict)), errID, errMsg);
% Validate if data size is equal to the length of dim_names:
errID = 'Umitoolbox:save2Dat:IncompatibleSize';
errMsg = 'The number of dimensions of data is different from the number of dimension names.';
assert(isequaln(ndims(data),numel(dim_names)), errID, errMsg);
% Check if "E" exists in dim_names and verify if event info exists in "s"
% struct:
% if ismember('E', dim_names)
%     fn = fieldnames(metaDat);
%     errID = 'Umitoolbox:save2Mat:MissingInfo';
%     errMsg = 'An event dimension name ("E") was detected but no event info ("eventID" and "eventNameList") was found in metaData.';
%     assert(all(ismember({'eventID', 'eventNameList'}, fn)), errID, errMsg);
% end
% Create metaData structure (if not provided):
if isempty(fieldnames(metaDat))
    metaDat = struct;
    metaDat.datName = 'data';
    metaDat.datFile = DatFileName;
    ds = size(data);
    metaDat.datSize = ds([1 2]);
    metaDat.datLength = ds(3:end); % Accounts for 3+ dimensions.
    metaDat.Datatype = class(data);
    metaDat.dim_names = dim_names;
else
    metaDat(1).dim_names = dim_names;
end


% Create an unique file identifier. To be used by class PIPELINEMANAGER.
metaDat(1).fileUUID = char(java.util.UUID.randomUUID);


% Save data to .dat file and metaData to .mat file:
if any(strcmp(flag, {'-a', 'append'})) && exist(DatFileName, 'file')
    fid = fopen(DatFileName, 'a');
elseif any(strcmp(flag, {'-a', 'append'})) && ~exist(DatFileName, 'file')
    errID = 'UMIToolbox:save2Dat:FileNotFound';
    msg = ['Cant save file. The file ' DatFileName ' does not exist in Matlab''s path'];
    throwAsCaller(MException(errID,msg))
else
    if exist(DatFileName, 'file') % Delete existing files if "-w" option is chosen.
        delete(DatFileName);
    end
    fid = fopen(DatFileName, 'w');
    saveMetaData(metaDatFilename, metaDat);
end
fwrite(fid, data, class(data));
fclose(fid);
end

function saveMetaData(metaDatFilename, metaDat)
% SAVEMETADATA saves METADAT fields to METADATFILENAME.

if length(metaDat) > 1
    fn = fieldnames(metaDat);
    s = struct;
    for i=1:length(fn)
        s.(fn{i}) = {metaDat.(fn{i})};
    end
    save(metaDatFilename, '-struct', 's', '-v7.3');
else
    save(metaDatFilename, '-struct', 'metaDat', '-v7.3');
end
end


