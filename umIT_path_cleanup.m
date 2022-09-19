function umIT_path_cleanup(currPath)
% This function removes all paths that contains the umIT folder structure and
% any folder containing the name "ioi_ana" from the current matlab path.
% Finally, this function adds the toolbox's folder and subfolders from "currPath".

P = strsplit(path,';')';
search_str = ['GUI' filesep 'StatsModule'];
myRoot = unique(cellfun(@(x) erase(x, search_str),...
    P(contains(P, search_str) | contains(P, 'ioi_ana','IgnoreCase',true)), 'UniformOutput',false));
% Disable warnings
warning('off')
% Remove all "Umit" paths:
for i = 1:length(myRoot)
    rmpath(genpath(myRoot{i}))
end
% Add current Path to Matlab's path:
addpath(genpath(currPath));
warning('on')
end