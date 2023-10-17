function saveDat(DatFileName, data)
% SAVE2DAT saves data to binary .dat files.
% Inputs:
% DatFileName (str): fullpath of .DAT file.
% data (numerical array): non-empty "single" numeric multi-dimensional matrix.

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) validateattributes(x, {'single'}, {'nonempty'}));
parse(p, DatFileName, data);
%%%%%%
DatFileName = p.Results.DatFileName;
data = p.Results.data;
clear p
% Apply filename extension:
if ~endsWith(DatFileName, '.dat')
    DatFileName = [DatFileName, '.dat'];
end

disp('Writing data to .DAT file ...');
fid = fopen(DatFileName, 'w');
fwrite(fid, data, 'single'); % Write data as single precision.
fclose(fid);
disp(['Data saved in : "' DatFileName '"']);
end