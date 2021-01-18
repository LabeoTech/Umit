function names = getNamesFromDir(Folder, exp, isDir, flag)
% This function extracts the names of elements from a structure obtained
% using the built-in "dir" function.
%   Inputs are the path "Folder" where the elements are, a regular expression "exp"
%   and a boolean "isDir", with true if the element is a directory.
arguments
   Folder {mustBeFolder}
   exp {mustBeText}
   isDir {mustBeNumericOrLogical} = true;
   flag {mustBeMember(flag, {'match', 'tokens'})} = 'match'
end
names = dir(Folder);
names = {names([names.isdir] == isDir).name};
names = regexp(names, exp, flag);
names = [names{:}];

end