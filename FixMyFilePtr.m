function FixMyFilePtr(protObj)
% Temporary fix for FilePtr.json files.

for i = 1:numel(protObj.Array.ObjList)
    tmpS = protObj.Array.ObjList(i);
    info = readFilePtr(tmpS.FilePtr);
    if ~isempty(info.Files)
%         info = tokenizeIt(info,tmpS);
%         info = addCreateDateTime(info, tmpS);
        info = fillFileInfo(info, tmpS);
        writeFilePtr(tmpS.FilePtr, info)
    end
    for j = 1:numel(tmpS.Array.ObjList)
        tmpA = tmpS.Array.ObjList(j);
        info = readFilePtr(tmpA.FilePtr);
        if ~isempty(info.Files)
%             info = tokenizeIt(info,tmpA);
%             info = addCreateDateTime(info, tmpA);
            info = fillFileInfo(info, tmpA);
            writeFilePtr(tmpA.FilePtr, info)
        end
        for k = 1:numel(tmpA.Array.ObjList)
            tmpM = tmpA.Array.ObjList(k);
            info = readFilePtr(tmpM.FilePtr);
            if ~isempty(info.Files)
%                 info = tokenizeIt(info,tmpM);
%                 info = addCreateDateTime(info, tmpM);
            else
                info = fillFileInfo(info, tmpM);
            end
            writeFilePtr(tmpM.FilePtr, info)
        end
    end
end
disp('Finished fix!');
end

% Fix functions :
function info = fillFileInfo(info, myObj)
folder = myObj.SaveFolder;
files = dir([folder '\*.dat']);
for i =1:length(files)
    file_metaData = matfile(strrep(fullfile(folder, files(i).name), '.dat', '_info.mat'));
    newStruct = struct('Name', files(i).name, 'UUID', file_metaData.fileUUID, 'Folder', folder, 'InputFile_Path', 'missing',...
                'InputFile_UUID', 'missing', 'creationDateTime', datestr(files(i).datenum), 'FunctionInfo', struct('Name', 'missing', 'DateNum', files(i).datenum, 'Job', 'missing'));
    info.Files = [info.Files newStruct];
end
end
function info = tokenizeIt(info, myObj)
for i = 1:numel(info.Files)
    Folder = info.Files(i).Folder;
    if strcmp(Folder, 'SaveFolder')
        Folder = myObj.SaveFolder;
    end
    info.Files(i).Folder = tokenizePath(Folder, myObj);
    info.Files(i).InputFile_Path = tokenizePath(info.Files(i).InputFile_Path, myObj);
end
disp('File Tokenized!')
end
function info = addCreateDateTime(info, myObj)
if isfield(info, 'creationDateTime')
    return
end
for i = 1:numel(info.Files)
    file_info = dir(tokenizePath(fullfile(info.Files(i).Folder, info.Files(i).Name), myObj, 'detokenize'));
    creationDateTime = datestr(file_info.datenum);
    info.Files(i).creationDateTime = creationDateTime;
end
disp('DateTime created!')
end

% File Pointer read/write functions:
function out = readFilePtr(fileptr)
txtIn = fileread(fileptr);
out= jsondecode(txtIn);
end

function writeFilePtr(fileptr, info)
txtOut = jsonencode(info);
fid = fopen(fileptr, 'w');
fprintf(fid, '%s', txtOut);
fclose(fid);
end