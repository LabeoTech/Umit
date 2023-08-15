function fList = checkFolder(fList, flag)
%   This Function checks if a Folder exists. If not, it creates it and
%   adds a file separator ("/" or "\") at the end of the string of a path.
if ~exist('flag','var')
    flag = '';
end
if ischar(fList)
    fList = check(fList,flag);
    return
end
for ii = 1:length(fList)
    fList{ii} = check(fList{ii},flag);
end
end

function Folder = check(Folder,flag)
if ~isfolder(Folder) && strcmp(flag, 'new')
    mkdir(Folder)
end
if Folder(end) ~= filesep
    Folder(end+1) = filesep;
end
end
