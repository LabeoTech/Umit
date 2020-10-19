function names = getNamesFromDir(Folder, isDir)
% This function extracts the names of elements from a structure obtained
% using the built-in "dir" function.
%   Inputs are the path "Folder" where the elements are and a boolean 
%   "isDir", with true if the element is a directory.

names = dir(Folder);
names = setdiff({names([names.isdir] == isDir).name}, {'.', '..'});

end