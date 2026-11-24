function outData = pmReuseGainStep(data, SaveFolder, varargin)
%PMREUSEGAINSTEP Test fixture for the reuse tests: counts its calls, returns data * Gain.
%
%   Used by TestPipelineManagerReuseExistingFiles. The Gain parameter lets a test
%   change a step's parameters between two runs; each real execution increments
%   the counter getappdata(0, 'pmSkipLog').pmReuseGainStep.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

gain = 1;
for k = 1:2:numel(varargin)
    if strcmpi(char(string(varargin{k})), 'Gain')
        gain = varargin{k + 1};
    end
end

if ischar(data) || (isstring(data) && isscalar(data))
    fileName = char(string(data));
    if ~isfile(fileName)
        fileName = fullfile(SaveFolder, fileName);
    end
    data = loadData(fileName);
end

key = mfilename;
callLog = getappdata(0, 'pmSkipLog');
if isempty(callLog)
    callLog = struct();
end
if isfield(callLog, key)
    callLog.(key) = callLog.(key) + 1;
else
    callLog.(key) = 1;
end
setappdata(0, 'pmSkipLog', callLog);

outData = single(data) * gain;
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo(mfilename, ...
    'Test fixture: counts its calls and multiplies the data by Gain.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Input data.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, 'Gain', 'parameter', 'Gain factor.', ...
    'kind', 'parameter', 'position', 3, 'callType', 'namevalue', ...
    'default', 1, 'allowed', [0 10], 'dataType', 'double');
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'Input times Gain.', [mfilename '.dat'], 1, 'isData', true);
end
