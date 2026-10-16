function categoryFolder = createAnalysisFixtureCategory(projectRoot, files)
%CREATEANALYSISFIXTURECATEGORY Put test functions where PipelineManager finds them.
%
%   categoryFolder = createAnalysisFixtureCategory(projectRoot, files)
%
%   PipelineManager discovers functions only in Analysis/<category>/*.m, so
%   PipelineManager tests write their fixture functions into a temporary
%   category there and add it to the path. FILES is a struct array with
%   fields:
%       name    - function name (the file is <name>.m)
%       source  - path of a .m file to copy, or
%       text    - the file contents (used when non-empty)
%
%   The folder is Analysis/pmTestFixture_<unique>. If writing fails, it is
%   removed before the error is rethrown. Callers remove it in their
%   teardown with removeAnalysisFixtureCategory. A killed MATLAB runs no
%   teardown, so each call first removes categories left by earlier runs
%   (older than an hour; a concurrent run's folder is never touched).
%
%   See also removeAnalysisFixtureCategory.

arguments
    projectRoot (1,:) char
    files struct
end

analysisFolder = fullfile(projectRoot, 'Analysis');
iRemoveStaleCategories(analysisFolder, projectRoot);

[~, uniqueName] = fileparts(tempname);
categoryFolder = fullfile(analysisFolder, [analysisFixturePrefix() uniqueName]);
mkdir(categoryFolder);
try
    for iFile = 1:numel(files)
        target = fullfile(categoryFolder, [files(iFile).name '.m']);
        if isfield(files, 'text') && ~isempty(files(iFile).text)
            fid = fopen(target, 'w');
            fwrite(fid, files(iFile).text, 'char');
            fclose(fid);
        else
            copyfile(files(iFile).source, target);
        end
    end
    addpath(categoryFolder, '-begin');
    for iFile = 1:numel(files)
        clear(files(iFile).name);
    end
    rehash;
catch ME
    removeAnalysisFixtureCategory(categoryFolder, projectRoot);
    rethrow(ME);
end
end

function iRemoveStaleCategories(analysisFolder, projectRoot)
% Categories from killed runs: the fixture prefix, or the older tempname
% form "tp<hex>_..." holding only .m files. Real categories never match.
listing = dir(analysisFolder);
listing = listing([listing.isdir] & ~startsWith({listing.name}, '.'));
for iFolder = 1:numel(listing)
    name = listing(iFolder).name;
    isFixture = startsWith(name, analysisFixturePrefix());
    isLegacy = ~isempty(regexp(name, '^tp[0-9a-f]{8}_', 'once'));
    if ~isFixture && ~isLegacy
        continue
    end
    age = datetime('now') - datetime(listing(iFolder).datenum, 'ConvertFrom', 'datenum');
    if age < hours(1)
        continue
    end
    folder = fullfile(analysisFolder, name);
    contents = dir(folder);
    contents = contents(~[contents.isdir]);
    if isLegacy && (isempty(contents) || ~all(endsWith({contents.name}, '.m')))
        continue
    end
    removeAnalysisFixtureCategory(folder, projectRoot);
end
end
