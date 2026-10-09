function outputs = pmCollectScenarioOutputs(prepareFcn, funcName, varargin)
%PMCOLLECTSCENARIOOUTPUTS Run one node under each RAM scenario, collect outputs.
%
%   outputs = pmCollectScenarioOutputs(prepareFcn, funcName)
%   outputs = pmCollectScenarioOutputs(prepareFcn, funcName, 'Input', 'green.dat')
%
%   PREPAREFCN is a function handle taking no arguments and returning the path
%   to a freshly prepared, PipelineManager-ready SaveFolder. It is called once
%   per scenario so that no scenario observes another scenario's leftovers.
%
%   For each scenario, PMCOLLECTSCENARIOOUTPUTS snapshots the folder's data
%   files, executes the single-step pipeline, and records which data files are
%   new. The returned struct array has fields:
%       .scenario - scenario name.
%       .files    - sorted cellstr of newly created data file names.
%       .folder   - the SaveFolder that scenario ran in (kept for inspection).
%
%   The point of collecting outputs this way is that the set of files a node
%   produces is part of that node's contract and must NOT depend on which RAM
%   scenario PipelineManager happened to choose: RAM availability is a resource
%   decision, not a data-identity decision. Tests compare .files across
%   scenarios to assert that invariance.
%
%   Name-Value options:
%       'Input'      - Input reference forwarded to addStep. Always supply this
%                      for functions declaring a DATA input.
%       'SaveAs'     - saveas target(s) forwarded to addStep.
%       'Extensions' - Data-file extensions to track.
%                      Default: {'.dat', '.umt'}.
%       'Scenarios'  - Scenario names to run. Default: pmRamScenarioList().
%       'PMOptions'  - Cell array of extra Name-Value pairs forwarded to
%                      buildPMForScenario (e.g. {'RawFolder', rawPath}).
%
%   NOTE: this helper deliberately does not trap errors or warnings. Do not use
%   it for a scenario that is contracted to fail — a function whose DATA input
%   has supportsFile==false is refused under 'ramsafe-strict'. Restrict
%   'Scenarios' and assert the failure separately with verifyError.
%
%   See also buildPMForScenario, pmRamScenarioList, pmRamScenarioExpectation.

p = inputParser;
addRequired(p, 'prepareFcn', @(x) isa(x, 'function_handle'));
addRequired(p, 'funcName', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Input', []);
addParameter(p, 'SaveAs', '');
addParameter(p, 'Extensions', {'.dat', '.umt'}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'Scenarios', pmRamScenarioList(), @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'PMOptions', {}, @iscell);
parse(p, prepareFcn, funcName, varargin{:});

funcName = char(string(p.Results.funcName));
extensions = cellstr(string(p.Results.Extensions));
scenarios = cellstr(string(p.Results.Scenarios));

outputs = struct( ...
    'scenario', cell(1, numel(scenarios)), ...
    'files', cell(1, numel(scenarios)), ...
    'folder', cell(1, numel(scenarios)));

for iScenario = 1:numel(scenarios)
    scenario = scenarios{iScenario};

    workFolder = p.Results.prepareFcn();
    workFolder = char(string(workFolder));

    before = iListDataFiles(workFolder, extensions);

    pmArgs = {'Input', p.Results.Input, 'SaveAs', p.Results.SaveAs};
    pm = buildPMForScenario(workFolder, funcName, scenario, ...
        pmArgs{:}, p.Results.PMOptions{:});
    result = pm.executePipeline('PrintSummary', false);
    assert(strcmpi(char(string(result.status)), 'completed'), ...
        'Umitoolbox:PMTestHelpers:PipelineExecutionFailed', ...
        '%s failed under %s: %s', funcName, scenario, ...
        char(strjoin(string(result.globalPipeLog.Messages_short), ' | ')));

    after = iListDataFiles(workFolder, extensions);

    outputs(iScenario).scenario = scenario;
    outputs(iScenario).files = sort(setdiff(after, before));
    outputs(iScenario).folder = workFolder;
end
end

% =========================================================================
% Local helpers
% =========================================================================
function names = iListDataFiles(folderPath, extensions)
%ILISTDATAFILES List data file names in a folder, filtered by extension.

names = {};
for iExt = 1:numel(extensions)
    ext = char(string(extensions{iExt}));
    if ~startsWith(ext, '.')
        ext = ['.' ext]; %#ok<AGROW>
    end
    found = dir(fullfile(folderPath, ['*' ext]));
    found = found(~[found.isdir]);
    if ~isempty(found)
        names = [names, {found.name}]; %#ok<AGROW>
    end
end
names = unique(names);
end
