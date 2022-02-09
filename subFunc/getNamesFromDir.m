function names = getNamesFromDir(Folder, exp, isDir, varargin)
% This function extracts the names of elements from a structure obtained
% using the built-in "dir" function.
% Inputs are the path "Folder" where the elements are, a regular expression "exp"
% and a boolean "isDir", with true if the element is a directory.

% Arguments validation
p = inputParser;
addRequired(p, 'Folder', @isfolder);
validate_exp = @(x) ischar(x) || isstring(x) && ~isempty(x);
addRequired(p, 'exp', validate_exp);
validate_isDir = @(x) isnumeric(x) || islogical(x);
addRequired(p, 'isDir', validate_isDir)
validate_flag = @(x) mustBeMember(x,{'match', 'tokens'});
addOptional(p, 'flag', 'match', validate_flag);
addOptional(p,'fullpath', true, @islogical)
parse(p,Folder,exp,isDir,varargin{:});
%%%%%
names = dir(p.Results.Folder);
names = {names([names.isdir] == p.Results.isDir).name};
names = regexp(names, p.Results.exp, p.Results.flag);
names = [names{:}];
if p.Results.fullpath
    names = fullfile(Folder,names);
end
end