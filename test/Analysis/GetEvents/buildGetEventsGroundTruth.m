function buildGetEventsGroundTruth()
% BUILDGETEVENTSGROUNDTRUTH Regenerate ground-truth events.mat files for all
% getEvents fixture cases.
%
% Each case must contain:
%   <caseFolder>/Raw/...
% or:
%   <caseFolder>/info.txt
%
% Optional:
%   <caseFolder>/caseConfig.mat with variable "cfg", where:
%       cfg.NameValue = {...};   % Name-Value pairs passed to getEvents
%
% Ground truth is written to:
%   <caseFolder>/groundTruth/events.mat

rootDir = fileparts(mfilename('fullpath'));
cases = localDiscoverCases(rootDir);

assert(~isempty(cases), 'No getEvents fixture cases were found.');

for iCase = 1:numel(cases)
    c = cases(iCase);

    gtFolder = fullfile(c.caseFolder, 'groundTruth');
    if ~isfolder(gtFolder)
        mkdir(gtFolder);
    end

    gtFile = fullfile(gtFolder, 'events.mat');
    if isfile(gtFile)
        delete(gtFile);
    end

    fprintf('Building ground truth for case: %s\n', c.name);
    outFile = getEvents(c.rawFolder, gtFolder, c.cfg.NameValue{:});
    outPath = fullfile(gtFolder, outFile);

    assert(isfile(outPath), ...
        'Ground-truth generation failed for case "%s".', c.name);

    fprintf('  -> %s\n', outPath);
end

fprintf('Ground-truth generation completed.\n');

end

% =========================================================================
% Local helpers
% =========================================================================
function cases = localDiscoverCases(rootDir)

d = dir(rootDir);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));

cases = struct( ...
    'name', {}, ...
    'caseFolder', {}, ...
    'rawFolder', {}, ...
    'cfg', {});

for ii = 1:numel(d)
    caseFolder = fullfile(rootDir, d(ii).name);

    if isfolder(fullfile(caseFolder, 'Raw'))
        rawFolder = fullfile(caseFolder, 'Raw');
    elseif isfile(fullfile(caseFolder, 'info.txt'))
        rawFolder = caseFolder;
    else
        continue
    end

    if ~isfile(fullfile(rawFolder, 'info.txt'))
        continue
    end

    cfg = localLoadCaseConfig(caseFolder);

    cases(end+1).name = d(ii).name; %#ok<AGROW>
    cases(end).caseFolder = caseFolder;
    cases(end).rawFolder = rawFolder;
    cases(end).cfg = cfg;
end
end

function cfg = localLoadCaseConfig(caseFolder)

cfg = struct();
cfg.NameValue = {};

cfgFile = fullfile(caseFolder, 'caseConfig.mat');
if ~isfile(cfgFile)
    return
end

S = load(cfgFile);
assert(isfield(S, 'cfg') && isstruct(S.cfg), ...
    'caseConfig.mat in "%s" must contain a struct variable named cfg.', caseFolder);

if isfield(S.cfg, 'NameValue')
    cfg.NameValue = S.cfg.NameValue;
end
end