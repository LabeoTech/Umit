function out = ProtoFunc4(MainDir)
% User-defined function. This Function is used by "Protocol" object to generate the
% list of Subjects and Acquisitions. In addition, this function will be used to update
% the content of the MainDir (Folder where the Data is).


% Regular expressions:
expSubj = 'Mouse_\d*'; % For Subject ID
expAcq = 'Acq_\d*'; % For Acquisition ID
expLabeo = '(ai|img)_\d*.bin'; % For Labeo Data
expDLC = 'DLCdata_\d*.txt'; % For DeepLabCut Data
cd(MainDir);
namesSubj = getNamesFromDir(MainDir, expSubj,1); % Extract Subject Names

namesAcq = arrayfun(@(x) getNamesFromDir(x{:}, expAcq, 1), namesSubj, 'UniformOutput', false);
tmpAcq = cellfun(@(a) cellfun(@(x) Acquisition(x,[]), a), namesAcq, 'UniformOutput', false);
FolderLabeo = cellfun(@(a,b) cellfun(@(x) [MainDir b filesep x filesep 'Labeo'], a,...
     'UniformOutput', false),namesAcq, namesSubj, 'UniformOutput', false);
FolderDLC = cellfun(@(a,b) cellfun(@(x) [MainDir b filesep x filesep 'DLC'], a,...
     'UniformOutput', false),namesAcq, namesSubj, 'UniformOutput', false);
filesLabeo = cellfun(@(a,b) cellfun(@(x) getNamesFromDir([MainDir b filesep x filesep 'Labeo'], expLabeo, 0), a,...
     'UniformOutput', false),namesAcq, namesSubj, 'UniformOutput', false);
filesDLC = cellfun(@(a,b) cellfun(@(x) getNamesFromDir([MainDir b filesep x filesep 'DLC'], expDLC, 0), a,...
     'UniformOutput', false),namesAcq, namesSubj, 'UniformOutput', false);
cellfun(@(a) arrayfun(@(x) x.Array.addObj(Labeo),a),tmpAcq)
cellfun(@(a) arrayfun(@(x) setfield(x.Array.ObjList(1), 'ID', x.ID),a),tmpAcq, 'UniformOutput', false);
cellfun(@(a,b) arrayfun(@(x,y) setfield(x.Array.ObjList(1), 'Folder', y{:}),a,b),...
    tmpAcq, FolderLabeo, 'UniformOutput', false);
cellfun(@(a,b) arrayfun(@(x,y) setfield(x.Array.ObjList(1), 'FileName', y{:}),a,b),...
    tmpAcq, filesLabeo, 'UniformOutput', false);
cellfun(@(a) arrayfun(@(x) x.Array.addObj(EyeMovement),a),tmpAcq)
cellfun(@(a) arrayfun(@(x) setfield(x.Array.ObjList(2), 'ID', x.ID),a),tmpAcq, 'UniformOutput', false);
cellfun(@(a,b) arrayfun(@(x,y) setfield(x.Array.ObjList(2), 'Folder', y{:}),a,b),...
    tmpAcq, FolderDLC, 'UniformOutput', false);
cellfun(@(a,b) arrayfun(@(x,y) setfield(x.Array.ObjList(2), 'FileName', y{:}),a,b),...
    tmpAcq, filesDLC, 'UniformOutput', false);

for i = 1:length(namesSubj)
    tmpSubj{i} = Subject();
    tmpSubj{i}.ID = namesSubj{i};
    tmpSubj{i}.GroupID = 'test';
    tmpSubj{i}.Array.addObj(tmpAcq{i})
end
    out = [tmpSubj{:}];
end


    
    
    
    
    
