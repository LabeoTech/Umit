function FixMyInheritance(protObj)
% Temporary fix 

for i = 1:numel(protObj.Array.ObjList)
    tmpS = protObj.Array.ObjList(i);
    for j = 1:numel(tmpS.Array.ObjList)
        tmpA = tmpS.Array.ObjList(j);
        for k = 1:numel(tmpA.Array.ObjList)
            tmpM = tmpA.Array.ObjList(k);
            info = readFilePtr(tmpM.FilePtr);
            if ~isempty(info.Files)
                if any(strcmp('SF_TF_AVG.dat', {info.Files.Name}))
                    mD1 = matfile(fullfile(tmpM.SaveFolder, 'mov_aligned_info.mat'));
                    mD2 = matfile(fullfile(tmpM.SaveFolder, 'SF_TF_AVG_info.mat'));
                    mD2.Properties.Writable = true;
                    props = setdiff(properties(mD1), properties(mD2));
                    for m = 1:length(props)
                        eval(['mD2.' props{m} '= mD1.' props{m} ';'])
                    end
                    disp(['Data inheritance done for ' tmpM.MyParent.ID])
                end
            end
        end
    end
end
disp('Finished fix!');
end

function out = readFilePtr(fileptr)
txtIn = fileread(fileptr);
out= jsondecode(txtIn);
end