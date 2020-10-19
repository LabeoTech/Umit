function SubjectInfo = ProtoFunc3(MainDir)
% This function creates a structure (struct) with the information of the
% recordings for each Subject.

% Regular expressions:
expSubj = 'Mouse_\d*'; % For Subject ID
expAcq = 'Acq_\d*'; % For Acquisition ID
expLabeo = '(ai|img)_\d*.bin'; % For Labeo files
expDLC = 'DLCdata_\d*.txt'; % For DeepLabCut files
%
SubjFolders = getNamesFromDir(MainDir,1);
PathSubj = regexp(SubjFolders, expSubj, 'match');
PathSubj = [PathSubj{:}];
SubjectInfo(numel(PathSubj)) = struct;
for i = 1:numel(PathSubj)
    AcqFolders = getNamesFromDir(PathSubj{i}, 1);
    AcqFolders = regexp(AcqFolders, expAcq,'match'); AcqFolders = [AcqFolders{:}];
    Acquisition(numel(AcqFolders)) = struct;
    for j = 1:numel(AcqFolders)
        ModFolders = getNamesFromDir([PathSubj{i} filesep AcqFolders{j}],1);
        ModFolders = ModFolders(strcmp(ModFolders, 'DLC') | strcmp(ModFolders, 'Labeo'));
        Modality(numel(ModFolders)) = struct;
        for k = 1:numel(ModFolders)
            if strcmp(ModFolders{k},'DLC')
                ModalityType = 'EyeMovement';
                exp = expDLC;
            else
                ModalityType = 'Labeo';
                exp = expLabeo;
            end
            Files = getNamesFromDir([PathSubj{i} filesep AcqFolders{j} filesep ModFolders{k}], 0);
            Files = regexp(Files, exp,'match'); Files = [Files{:}];
            Modality(k).Type = ModalityType;
            Modality(k).Folder = [MainDir PathSubj{i} filesep AcqFolders{j} filesep ModFolders{k}];
            Modality(k).Files = Files;
        end
        Acquisition(j).ID = AcqFolders{j};
        Acquisition(j).Modality = Modality;
        clear Modality
    end
    SubjectInfo(i).ID = PathSubj{i};
    SubjectInfo(i).GroupID = 'test';
    SubjectInfo(i).Acquisition = Acquisition;
    clear Acquisition
    
end
end


    
    
    
    
    
