function deleteDatFile(file, FilePtr)
% DELETEDATFILE deletes a file from a folder and from the file pointer.
% Inputs:
% file : full path of .DAT file to be deleted.
% FilePtr: full path to file pointer JSON file listing the .DAT file to be
% deleted.
info = readPtr(FilePtr);
[~, filename, ext] = fileparts(file);
idx = strcmp([filename,ext], {info.Files.Name});
info.Files(idx) = [];
if sum(idx) == 1
    delete(file);
    delete(strrep(file,'.dat', '_info.mat'));
    disp('File deleted!')
    write2Ptr(info,FilePtr);
elseif sum(idx) == 0
    disp('File not found in File pointer')
else
    disp('More than one file with the same name was found in File Pointer. Fix this and try again!')
end
end

function info = readPtr(filePtr)
txt = fileread(filePtr);
info = jsondecode(txt);
end

function write2Ptr(info, filePtr)
txtOut = jsonencode(info);
fid = fopen(filePtr, 'w');
fprintf(fid,'%s',txtOut);
fclose(fid);
disp('File Pointer updated.');
end