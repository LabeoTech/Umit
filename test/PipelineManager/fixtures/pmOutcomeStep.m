function outData = pmOutcomeStep(data, SaveFolder, varargin)
%PMOUTCOMESTEP Test fixture for the run-outcome tests: returns data, optionally warns or fails.
%
%   Used by TestPipelineManagerRunOutcomeSummary. The Mode parameter selects
%       'ok'       return the data unchanged,
%       'warn'     raise warning 'Umitoolbox:pmOutcomeStep:Test', then return the data,
%       'warnfail' raise that warning, then error 'Umitoolbox:pmOutcomeStep:Failed',
%       'fail'     error 'Umitoolbox:pmOutcomeStep:Failed'.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

mode = 'ok';
for k = 1:2:numel(varargin)
    if strcmpi(char(string(varargin{k})), 'Mode')
        mode = char(string(varargin{k + 1}));
    end
end

if ischar(data) || (isstring(data) && isscalar(data))
    fileName = char(string(data));
    if ~isfile(fileName)
        fileName = fullfile(SaveFolder, fileName);
    end
    data = loadData(fileName);
end

if any(strcmp(mode, {'warn', 'warnfail'}))
    warning('Umitoolbox:pmOutcomeStep:Test', 'Test warning raised by pmOutcomeStep.');
end
if any(strcmp(mode, {'fail', 'warnfail'}))
    error('Umitoolbox:pmOutcomeStep:Failed', 'Test failure raised by pmOutcomeStep.');
end

outData = single(data);
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo(mfilename, ...
    'Test fixture: returns its data, optionally warning or failing.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Input data.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, 'Mode', 'parameter', 'Outcome of the step.', ...
    'kind', 'parameter', 'position', 3, 'callType', 'namevalue', ...
    'default', 'ok', 'allowed', {'ok', 'warn', 'warnfail', 'fail'}, 'dataType', 'char');
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'The input data.', [mfilename '.dat'], 1, 'isData', true);
end
