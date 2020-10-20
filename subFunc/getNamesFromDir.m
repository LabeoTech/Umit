function names = getNamesFromDir(Folder, exp, isDir)
% This function extracts the names of elements from a structure obtained
% using the built-in "dir" function.
%   Inputs are the path "Folder" where the elements are, a regular expression "exp"
%   and a boolean "isDir", with true if the element is a directory.

names = dir(Folder);
names = {names([names.isdir] == isDir).name};
names = regexp(names, exp, 'match');
names = [names{:}];

end