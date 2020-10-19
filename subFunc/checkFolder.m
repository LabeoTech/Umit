function Folder = checkFolder(Folder)
%   This Function adds a file separator ("/" or "\") at the end of the
%   string of a path.

if Folder(end) ~= filesep
    Folder(end+1) = filesep;
end
end