function outData = pmSkipStepB(data, SaveFolder, varargin)
%PMSKIPSTEPB Test fixture for the skipped-step tests: counts its calls, returns data + 1.
%
%   Used by TestPipelineManagerSkipSteps. Each real execution increments the
%   counter getappdata(0, 'pmSkipLog').pmSkipStepB, so a skipped step is
%   visible as a counter that does not move. pmSkipStepA is an identical copy
%   (own name and counter) used to build two-step pipelines.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
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

outData = single(data) + 1;
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo(mfilename, ...
    'Test fixture: counts its calls and adds 1 to the data.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Input data.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'Input plus one.', [mfilename '.dat'], 1, 'isData', true);
end
