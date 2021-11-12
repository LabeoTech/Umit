function save2Dat(DatFileName, data, metaData)
% SAVE2DAT saves data to binary .dat files. If it is a new file, it creates
% a .MAT file containing the data's METADATA. FLAG variable indicates if it
% appends ('-A', '-APPEND') to an existing file.
% Inputs:
% DatFileName (str): fullpath of .DAT file.
% data (numerical array): non-empty "single" numeric multi-dimensional matrix.
% metaData (struct): meta data associated with "data".

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) validateattributes(x, {'single'}, {'nonempty'}));
addOptional(p, 'metaData', struct(), @isstruct);
parse(p, DatFileName, data, metaData);
%%%%%%
DatFileName = p.Results.DatFileName;
data = p.Results.data;
metaData = upper(p.Results.dim_names);
clear p
% Further validate dim_names:
root = getenv('Umitoolbox');
dim_names_info = load(fullfile(root, 'subFunc','dimension_names.mat'));
errID = 'Umitoolbox:save2Dat:InvalidName';
errMsg = 'List of dimension names contain invalid values.';
assert(all(ismember(metaData.dim_names, dim_names_info.dims_dict)), errID, errMsg);
% Validate if data size is equal to the length of dim_names:
errID = 'Umitoolbox:save2Dat:IncompatibleSize';
errMsg = 'The number of dimensions of data is different from the number of dimension names.';
assert(isequaln(ndims(data),numel(metaData.dim_names)), errID, errMsg);
% Create an unique file identifier. To be used by class PIPELINEMANAGER.
metaData(1).fileUUID = char(java.util.UUID.randomUUID);
% Save meta data file:
save(strrep(DatFileName, '.mat', '.dat'), '-struct', 'metaData', '-v7.3');
% Save data to .dat file and metaData to .mat file:
if exist(DatFileName, 'file') % Delete existing files if "-w" option is chosen.
    delete(DatFileName);
end
disp('Writing data to .DAT file...');
fid = fopen(DatFileName, 'w');
fwrite(fid, data, class(data));
fclose(fid);
disp('Finished saving .DAT file.');
end