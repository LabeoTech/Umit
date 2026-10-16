function removeAnalysisFixtureCategory(categoryFolder, projectRoot)
%REMOVEANALYSISFIXTURECATEGORY Remove a temporary Analysis/ fixture category.
%
%   removeAnalysisFixtureCategory(categoryFolder, projectRoot)
%
%   Clears its functions from memory, removes it from the path and deletes
%   it. Refuses any folder that is not directly under Analysis/. An empty
%   CATEGORYFOLDER or a folder already gone is not an error.
%
%   See also createAnalysisFixtureCategory.

categoryFolder = char(string(categoryFolder));
if isempty(categoryFolder)
    return
end
analysisFolder = fullfile(projectRoot, 'Analysis');
if ~strcmpi(fileparts(categoryFolder), analysisFolder)
    error('Umitoolbox:TestHelpers:UnsafeCleanupPath', ...
        'Refusing to remove a folder that is not an Analysis category: %s', categoryFolder);
end
functionFiles = dir(fullfile(categoryFolder, '*.m'));
for iFile = 1:numel(functionFiles)
    [~, functionName] = fileparts(functionFiles(iFile).name);
    clear(functionName);
end
if any(strcmpi(strsplit(path, pathsep), categoryFolder))
    rmpath(categoryFolder);
end
if isfolder(categoryFolder)
    rmdir(categoryFolder, 's');
end
rehash;
end
