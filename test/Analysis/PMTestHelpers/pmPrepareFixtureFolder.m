function workFolder = pmPrepareFixtureFolder(parentFolder, caseName, sourceFolder, fileNames, varargin)
%PMPREPAREFIXTUREFOLDER Build an isolated PipelineManager-ready fixture folder.
%
%   workFolder = pmPrepareFixtureFolder(parentFolder, caseName, sourceFolder, fileNames)
%   workFolder = pmPrepareFixtureFolder(..., 'OptionalFiles', {path1, path2})
%   workFolder = pmPrepareFixtureFolder(..., 'RegisterChannels', {'green.dat'})
%
%   Creates PARENTFOLDER/CASENAME, copies FILENAMES into it from
%   SOURCEFOLDER, and registers the imported-channel entries PipelineManager
%   needs so executePipeline does not refuse the folder as legacy-schema.
%
%   Each RAM scenario must run in its own folder: comparing what two
%   scenarios produced is only meaningful if neither could see the other's
%   output files, and PipelineManager's overwrite-avoidance renames a new
%   output when a file of that name already exists. Passing a distinct
%   CASENAME per run is what keeps those runs independent.
%
%   Name-Value options:
%       'OptionalFiles'    - Cell array of full paths copied into the folder
%                            only when they exist. Use this for state a
%                            fixture generates elsewhere, such as an
%                            events.mat written into the parent SaveFolder by
%                            a TestMethodSetup helper.
%       'RegisterChannels' - Cell array of .dat names passed to
%                            ensurePMReadyAcqInfos. Default: every .dat
%                            copied in through FILENAMES. Registration
%                            happens before the pipeline runs, while only the
%                            input channels are present, so an analysis
%                            output can never be recorded as an imported
%                            channel.
%
%   See also ensurePMReadyAcqInfos, buildPMForScenario, pmCollectScenarioOutputs.

p = inputParser;
addRequired(p, 'parentFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'caseName', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'sourceFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'fileNames', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'OptionalFiles', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'RegisterChannels', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
parse(p, parentFolder, caseName, sourceFolder, fileNames, varargin{:});

parentFolder = char(string(p.Results.parentFolder));
caseName = char(string(p.Results.caseName));
sourceFolder = char(string(p.Results.sourceFolder));
fileNames = cellstr(string(p.Results.fileNames));

workFolder = fullfile(parentFolder, caseName);
assert(~isfolder(workFolder), ...
    'Umitoolbox:pmPrepareFixtureFolder:caseFolderExists', ...
    ['Case folder "%s" already exists. Each scenario run needs its own ' ...
    'folder, so caseName must be unique within the parent folder.'], workFolder);
mkdir(workFolder);

% -------------------------------------------------------------------------
% Required fixture files
% -------------------------------------------------------------------------
for iFile = 1:numel(fileNames)
    thisName = char(string(fileNames{iFile}));
    sourcePath = fullfile(sourceFolder, thisName);

    assert(isfile(sourcePath), ...
        'Umitoolbox:pmPrepareFixtureFolder:missingFixtureFile', ...
        'Required fixture file "%s" was not found.', sourcePath);

    copyfile(sourcePath, fullfile(workFolder, thisName));
end

% -------------------------------------------------------------------------
% Optional fixture files
% -------------------------------------------------------------------------
optionalFiles = p.Results.OptionalFiles;
if ~isempty(optionalFiles)
    optionalFiles = cellstr(string(optionalFiles));
    for iFile = 1:numel(optionalFiles)
        sourcePath = char(string(optionalFiles{iFile}));
        if ~isfile(sourcePath)
            continue
        end
        [~, baseName, ext] = fileparts(sourcePath);
        copyfile(sourcePath, fullfile(workFolder, [baseName ext]));
    end
end

% -------------------------------------------------------------------------
% Register imported channels
% -------------------------------------------------------------------------
registerChannels = p.Results.RegisterChannels;
if isempty(registerChannels)
    isDat = endsWith(lower(string(fileNames)), '.dat');
    registerChannels = fileNames(isDat);
else
    registerChannels = cellstr(string(registerChannels));
end

if ~isempty(registerChannels)
    ensurePMReadyAcqInfos(workFolder, registerChannels);
end
end
