function save2Dat(DatFileName, data, opts)
% SAVE2DAT saves data to binary .dat files. If it is a new file, it creates
% a .MAT file containing the data's METADATA. FLAG variable indicates if it
% appends ('-A', '-APPEND') to an existing file.

arguments
    DatFileName
    data
    opts.metaData = []
    opts.flag string {mustBeMember(opts.flag, {'-a', '-append', '-w', '-write'})} = '-w'
end

metaDat = opts.metaData;
metaDatFilename = strrep(DatFileName, '.dat', '_info.mat');
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
metaDat(1).fileUUID = convertStringsToChars(matlab.lang.internal.uuid);

if any(strcmp(opts.flag, {'-a', 'append'})) && exist(DatFileName, 'file')
    fid = fopen(DatFileName, 'a');
elseif any(strcmp(opts.flag, {'-a', 'append'})) && ~exist(DatFileName, 'file')
    errID = 'IsaToolbox:FileNotFound';
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
    save(metaDatFilename, '-struct', 's');
else
    save(metaDatFilename, '-struct', 'metaDat');
end
end


