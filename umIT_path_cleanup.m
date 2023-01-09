function umIT_path_cleanup(currPath)
% This function removes all paths that contains the umIT folder structure and
% any folder containing the name "ioi_ana" from the current matlab path.
% Finally, this function adds the toolbox's folder and subfolders from "currPath".
if contains(computer, 'PC')
    separator = ';';
else
    separator = ':';
end
P = strsplit(path,separator)';
search_str = ['docs' filesep 'css']; % Use UMIT's subfolders to identify where the toolbox is (and any other versions).
myRoot = unique(cellfun(@(x) erase(x, search_str),...
    P(contains(P, search_str) | contains(P, 'ioi_ana','IgnoreCase',true)), 'UniformOutput',false));
% Disable warnings
warning('off')
% Remove all "Umit" paths:
for ii = 1:length(myRoot)
    rmpath(genpath(myRoot{ii}))
end
% Add current Path to Matlab's path:
addpath(genpath(currPath));
warning('on')
end