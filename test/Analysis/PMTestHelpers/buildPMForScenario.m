function pm = buildPMForScenario(saveFolder, funcName, scenario, varargin)
%BUILDPMFORSCENARIO Configure a PipelineManager for one RAM execution scenario.
%
%   pm = buildPMForScenario(saveFolder, funcName, scenario)
%   pm = buildPMForScenario(saveFolder, funcName, scenario, 'Input', 'green.dat')
%
%   BUILDPMFORSCENARIO constructs a PipelineManager over SAVEFOLDER, applies
%   the ramMode / ramSafePolicy pair named by SCENARIO, adds a single step for
%   FUNCNAME, and returns the manager ready for executePipeline. It deliberately
%   does not execute anything and does not trap errors or warnings, so callers
%   keep matlab.unittest's native verifyError / verifyWarning qualifications.
%
%   SCENARIO is one of the values returned by PMRAMSCENARIOLIST:
%       'auto', 'ramsafe-strict', 'ramsafe-bestEffort'.
%
%   Name-Value options:
%       'Input'            - Input reference forwarded to addStep. Char, string
%                            or cell. Always pass this for functions that
%                            declare a DATA input: addStep prompts the user
%                            interactively when it cannot resolve a source,
%                            which would hang an automated test.
%       'SaveAs'           - saveas target(s) forwarded to addStep.
%       'RawFolder'        - Raw folder path, for functions consuming RawFolder.
%       'ProjectFolder'    - Project folder. Default: PipelineManager's own
%                            default (the current working folder).
%       'SkipSteps'        - Value for pm.b_skipSteps. Default: false, so a
%                            step always runs instead of being skipped from a
%                            previous result.
%       'LeafOutputPolicy' - Value for pm.leafOutputPolicy. Default: leave the
%                            PipelineManager default ('saveLeaves') untouched,
%                            which persists DATA leaf outputs under their
%                            declared defOutfilename.
%       'Parameters'       - Scalar struct of parameter-name/value overrides.
%                            PipelineManager's public parameter editor is a UI
%                            dialog, so automated tests apply overrides through
%                            the supported save/load .pipe representation.
%
%   Example:
%       pm = buildPMForScenario(folder, 'GSR', 'ramsafe-strict', ...
%           'Input', 'green.dat');
%       pm.executePipeline();
%
%   See also pmRamScenarioList, pmRamScenarioExpectation, PipelineManager.

p = inputParser;
addRequired(p, 'saveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'funcName', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'scenario', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Input', []);
addParameter(p, 'SaveAs', '');
addParameter(p, 'RawFolder', '');
addParameter(p, 'ProjectFolder', '');
addParameter(p, 'SkipSteps', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'LeafOutputPolicy', '');
addParameter(p, 'Parameters', struct(), @(x) isstruct(x) && isscalar(x));
parse(p, saveFolder, funcName, scenario, varargin{:});

saveFolder = char(string(p.Results.saveFolder));
funcName = char(string(p.Results.funcName));
scenario = char(string(p.Results.scenario));
rawFolder = char(string(p.Results.RawFolder));
projectFolder = char(string(p.Results.ProjectFolder));
leafPolicy = char(string(p.Results.LeafOutputPolicy));

% -------------------------------------------------------------------------
% Construct the manager
% -------------------------------------------------------------------------
ctorArgs = {{saveFolder}};
if ~isempty(rawFolder)
    ctorArgs{end+1} = {rawFolder};
elseif ~isempty(projectFolder)
    % ProjectFolder is the third positional input, so a placeholder RawFolder
    % is required to reach it.
    ctorArgs{end+1} = {''};
end
if ~isempty(projectFolder)
    ctorArgs{end+1} = projectFolder;
end

pm = PipelineManager(ctorArgs{:});

% -------------------------------------------------------------------------
% Apply the RAM execution scenario
% -------------------------------------------------------------------------
switch lower(scenario)
    case 'auto'
        pm.ramMode = 'auto';

    case 'ramsafe-strict'
        pm.ramMode = 'ramsafe';
        pm.ramSafePolicy = 'strict';

    case 'ramsafe-besteffort'
        pm.ramMode = 'ramsafe';
        pm.ramSafePolicy = 'bestEffort';

    otherwise
        error('Umitoolbox:buildPMForScenario:UnknownScenario', ...
            ['Unknown RAM scenario "%s". Expected one of the values returned ' ...
            'by pmRamScenarioList: auto, ramsafe-strict, ramsafe-bestEffort.'], ...
            scenario);
end

pm.b_skipSteps = p.Results.SkipSteps;

if ~isempty(leafPolicy)
    pm.leafOutputPolicy = leafPolicy;
end

% -------------------------------------------------------------------------
% Add the single step under test
% -------------------------------------------------------------------------
stepArgs = {funcName};
if ~isempty(p.Results.Input)
    stepArgs = [stepArgs, {'input', p.Results.Input}];
end
if ~isempty(p.Results.SaveAs)
    stepArgs = [stepArgs, {'saveas', p.Results.SaveAs}];
end

pm.addStep(stepArgs{:});

if ~isempty(fieldnames(p.Results.Parameters))
    pm = iApplyParameterOverrides(pm, p.Results.Parameters, saveFolder, funcName);
end
end

function pm = iApplyParameterOverrides(pm, overrides, saveFolder, funcName)
%IAPPLYPARAMETEROVERRIDES Round-trip overrides through the public pipe API.

pipeBase = tempname(saveFolder);
pipeFile = [pipeBase '.pipe'];
pm.savePipe(pipeFile);
cleanupObj = onCleanup(@() iDeleteIfExists(pipeFile));

loaded = load(pipeFile, 'pipeStruct', '-mat');
pipeStruct = loaded.pipeStruct;

streamIdx = [];
for iNode = 1:numel(pipeStruct.nodes)
    node = pipeStruct.nodes(iNode);
    if strcmpi(node.kind, 'stream') && isfield(node, 'info') && ...
            isfield(node.info, 'name') && strcmpi(node.info.name, funcName)
        streamIdx = iNode;
    end
end
assert(~isempty(streamIdx), ...
    'Umitoolbox:buildPMForScenario:MissingStreamNode', ...
    'Could not find stream node "%s" in the saved pipeline.', funcName);

params = pipeStruct.nodes(streamIdx).info.parameters;
overrideNames = fieldnames(overrides);
for iOverride = 1:numel(overrideNames)
    overrideName = overrideNames{iOverride};
    paramIdx = find(strcmpi({params.name}, overrideName), 1, 'first');
    assert(~isempty(paramIdx), ...
        'Umitoolbox:buildPMForScenario:UnknownParameter', ...
        'Function "%s" has no pipeline parameter named "%s".', ...
        funcName, overrideName);
    params(paramIdx).value = overrides.(overrideName);
end
pipeStruct.nodes(streamIdx).info.parameters = params;
save(pipeFile, 'pipeStruct', '-mat');

pm.loadPipe(pipeFile);
end

function iDeleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
