function names = pmDeclaredOutputNames(funcName, varargin)
%PMDECLAREDOUTPUTNAMES Default output file names an analysis function declares.
%
%   names = pmDeclaredOutputNames(funcName)
%   names = pmDeclaredOutputNames(funcName, 'DataOnly', true)
%
%   Returns the flattened cellstr of defOutfilename entries declared by
%   FUNCNAME's pipelineInfo outputs. When no explicit saveas target is given,
%   PipelineManager persists a DATA leaf output under this declared name, so
%   these are the file names a single-step pipeline is expected to leave in the
%   SaveFolder.
%
%   Reading the names from pipelineInfo rather than hard-coding them in each
%   test keeps the expectation tied to the function's own declaration, so a
%   renamed output surfaces as a real behavioural failure instead of a stale
%   literal in a test.
%
%   Name-Value options:
%       'DataOnly' - Only include outputs with isData==true. Default: false.
%
%   See also buildPMForScenario, pmCollectScenarioOutputs.

p = inputParser;
addRequired(p, 'funcName', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'DataOnly', false, @(x) islogical(x) && isscalar(x));
parse(p, funcName, varargin{:});

funcName = char(string(p.Results.funcName));
info = feval(funcName, 'pipelineInfo');

assert(isstruct(info) && isfield(info, 'outputs'), ...
    'Umitoolbox:pmDeclaredOutputNames:invalidPipelineInfo', ...
    '"%s" did not return a pipelineInfo struct with an outputs field.', funcName);

names = {};
for iOutput = 1:numel(info.outputs)
    thisOutput = info.outputs(iOutput);

    if p.Results.DataOnly
        if ~isfield(thisOutput, 'isData') || ~thisOutput.isData
            continue
        end
    end

    if ~isfield(thisOutput, 'defOutfilename') || isempty(thisOutput.defOutfilename)
        continue
    end

    names = [names, cellstr(string(thisOutput.defOutfilename))]; %#ok<AGROW>
end

names = unique(names, 'stable');
end
