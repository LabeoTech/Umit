function save2Dat(DatFileName, data, varargin)
% SAVE2DAT saves data to binary .dat files. If it is a new file, it creates
% a .MAT file containing the data's METADATA. FLAG variable indicates if it
% appends ('-A', '-APPEND') to an existing file.

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) validateattributes(x, {'single'}, {'nonempty'}));
addOptional(p, 'flag', '-w', @(x) ischar(x) || isstring(x) && mustBeMember(x, {'-w', '-write', '-a', '-append'}));
addOptional(p, 'metaData', struct(), @isstruct);
parse(p, DatFileName, data, varargin{:});
%%%%%%
metaDat = p.Results.metaData;
metaDatFilename = strrep(p.Results.DatFileName, '.dat', '_info.mat');
if isempty(metaDat)
    metaDat = struct;
    metaDat.datName = 'data';
    metaDat.datFile = DatFileName;
    ds = size(data);
    metaDat.datSize = ds([1 2]);
    metaDat.datLength = ds(3:end); % Accounts for 3+ dimensions.
    metaDat.Datatype = class(data);
end

% Create an unique file identifier. To be used by class PIPELINEMANAGER.
metaDat(1).fileUUID = char(java.util.UUID.randomUUID);

if any(strcmp(p.Results.flag, {'-a', 'append'})) && exist(p.Results.DatFileName, 'file')
    fid = fopen(p.Results.DatFileName, 'a');
elseif any(strcmp(p.Results.flag, {'-a', 'append'})) && ~exist(p.Results.DatFileName, 'file')
    errID = 'IsaToolbox:FileNotFound';
    msg = ['Cant save file. The file ' p.Results.DatFileName ' does not exist in Matlab''s path'];
    throwAsCaller(MException(errID,msg))
else
    if exist(p.Results.DatFileName, 'file') % Delete existing files if "-w" option is chosen.
        delete(p.Results.DatFileName);
    end
    fid = fopen(p.Results.DatFileName, 'w');
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
    save(metaDatFilename, '-struct', 's');
else
    save(metaDatFilename, '-struct', 'metaDat');
end
end


