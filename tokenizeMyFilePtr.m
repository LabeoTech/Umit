function tokenizeMyFilePtr(protObj)
 % Temporary fix to tokenize all paths inside FilePtr.json files.
 for i = 1:numel(protObj.Array.ObjList)
     tmpS = protObj.Array.ObjList(i);
     tokenizeIt(tmpS);
     for j = 1:numel(tmpS.Array.ObjList)
         tmpA = tmpS.Array.ObjList(j);
         tokenizeIt(tmpA);
         for k = 1:numel(tmpA.Array.ObjList)
             tmpM = tmpA.Array.ObjList(k);
             tokenizeIt(tmpM);
         end
     end
 end
end 

function tokenizeIt(myObj)
    txtIn = fileread(myObj.FilePtr);
    a = jsondecode(txtIn);
    if isempty(a.Files)
        disp('Files field is empty. Skipped!')
        return
    end
    for i = 1:numel(a.Files)
%         Folder = a.Files(i).Folder;
%         if strcmp(Folder, 'SaveFolder')
%             Folder = myObj.SaveFolder;
%         end
%         a.Files(i).Folder = tokenizePath(Folder, myObj);
%         a.Files(i).InputFile_Path = tokenizePath(a.Files(i).InputFile_Path, myObj);
       a.Files(i).creationDateTime = datestr(now);
    end
    txtOut = jsonencode(a);
    fid = fopen(myObj.FilePtr, 'w');
    fprintf(fid, '%s', txtOut);
    fclose(fid);
    disp('File Tokenized!')
end
