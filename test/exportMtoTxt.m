function exportMtoTxt(srcRoot, destRoot, funcList)
% exportMtoTxt Export .m files from srcRoot (and subfolders) to .txt files.
%   exportMtoTxt(srcRoot)                      -> writes .txt files next to .m files
%   exportMtoTxt(srcRoot, destRoot)            -> writes .txt files under destRoot,
%                                                 preserving relative paths
%   exportMtoTxt(srcRoot, destRoot, funcList)  -> only exports .m files whose
%                                                 primary function name is in funcList
%
% funcList can be a string or cell array of character vectors.

if nargin < 1 || isempty(srcRoot)
    error('Source root folder must be provided.');
end
if nargin < 2 || isempty(destRoot)
    destRoot = srcRoot;
end
if nargin < 3
    funcList = [];
end

% Normalize paths
srcRoot = char(java.io.File(srcRoot).getCanonicalPath);
destRoot = char(java.io.File(destRoot).getCanonicalPath);
% Find .m files (recursive)
files = dir(fullfile(srcRoot, '**', '*.m')); % requires R2016b+; otherwise use genpath
% Normalize funcList to cellstr of lower-case names for comparison
if ~isempty(funcList)    
    if ischar(funcList)
        funcSet = {lower(funcList)};
    elseif isstring(funcList)
        funcSet = cellstr(lower(funcList));
    elseif iscell(funcList)
        funcSet = lower(cellfun(@char, funcList, 'UniformOutput', false));
    else
        error('funcList must be a string or cell array of character vectors.');
    end

    files_full = lower({files.name});
    
    [idx,locb] = ismember(funcSet,files_full);

    assert(all(idx),'One or more functions do not exist. Check inputs for typos.');
    files = files(locb);


end


for k = 1:numel(files)
    if files(k).isdir
        continue
    end
    srcFull = fullfile(files(k).folder, files(k).name);

    % Read file and detect primary function name from 'function' line
    fileContent = fileread(srcFull);  
    if ~exist(destRoot, 'dir')
        mkdir(destRoot);
    end
    [~, name, ~] = fileparts(files(k).name);
    timestamp = datestr(now,'dd-mm-yyyy_hh-MM');
    destFull = fullfile(destRoot, [name, '_', timestamp, '.txt']);
    % Copy file (text content preserved)
    fid = fopen(destFull,'w');
    fprintf(fid,'%s',fileContent);
    fclose(fid);
end
end
