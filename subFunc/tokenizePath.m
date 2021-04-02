function out = tokenizePath(Path, object, varargin)
% TOKENIZEPATH replaces parts of PATH with tokens pointing paths saved in
% objects contained in Protocol class.
% Inputs:
%   Path: full path of a folder or file.
%   obj: object (Subject, Acquisition or Modality) associated with PATH.
%   translation(optional) : flag indicating the direction of translation:
%               'tokenize'(default): translates real path to tokenized path
%               'detokenize' : translates tokenized path to real path.
% Output:
%   out = tokenized/detokenized version of PATH.
% Argument Validation
p = inputParser;
addRequired(p, 'Path');
addRequired(p, 'object');
addOptional(p, 'translation', 'tokenize', @(x) ismember(x, {'tokenize', 'detokenize'}));
% Parse inputs:
parse(p,Path, object, varargin{:});
Path = p.Results.Path;
object = p.Results.object;
translation = p.Results.translation;
% Control for empty paths:
if isempty(Path)
    out = [];
    return
end
% Initialize list of tokens
tokens = struct('Key', '', 'Path', '');
tokens(1).Key = '$RootRaw$';
tokens(2).Key = '$RootSave$';
tokens(3).Key = '$SubjSave$';
tokens(4).Key = '$AcqSave$';

typeObj = class(object);
switch typeObj
    case 'Subject'
        tokens(1).Path = strip(object.MyParent.MainDir, filesep);
        tokens(2).Path = strip(object.MyParent.SaveDir, filesep);
    case 'Acquisition'
        tokens(1).Path = strip(object.MyParent.MyParent.MainDir, filesep);
        tokens(2).Path = strip(object.MyParent.MyParent.SaveDir, filesep);
        tokens(3).Path= [filesep object.MyParent.ID filesep];
    otherwise
        tokens(1).Path = strip(object.MyParent.MyParent.MyParent.MainDir, filesep);
        tokens(2).Path = strip(object.MyParent.MyParent.MyParent.SaveDir, filesep);
        tokens(3).Path = [filesep object.MyParent.MyParent.ID filesep];
        tokens(4).Path = [filesep object.MyParent.ID filesep];
end

if strcmp(translation, 'tokenize')
    indx = find(cellfun(@(x) ~isempty(strfind(Path, x)), {tokens.Path}));%#ok
    out = Path;
    for i = 1:numel(indx)
        out = strrep(out, strip(tokens(indx(i)).Path, filesep), tokens(indx(i)).Key);
    end
else
    str = split(Path, filesep);
    for i = 1:length(str)
        idx = strcmp(str{i}, {tokens.Key});
        if sum(idx) == 0 
            continue
        end
        str{i} = tokens(idx).Path;
    end
    out = fullfile(str{:});
end
end
