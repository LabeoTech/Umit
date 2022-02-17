function Folder = checkFolder(Folder, flag)
%   This Function checks if a Folder exists. If not, it creates it and
%   adds a file separator ("/" or "\") at the end of the string of a path.

if ~isfolder(Folder) && strcmp(flag, 'new')
    mkdir(Folder)
end
if Folder(end) ~= filesep
    Folder(end+1) = filesep;
end
end