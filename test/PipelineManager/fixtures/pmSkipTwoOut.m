function [outData, outAlt] = pmSkipTwoOut(data, SaveFolder, varargin)
%PMSKIPTWOOUT Test fixture with two DATA outputs: counts its calls.
%
%   Used by TestPipelineManagerSkipSteps. outData is the input plus 1 and
%   outAlt the input plus 10, so a consumer reading one port or the other
%   gives a different result. Each real execution increments the counter
%   getappdata(0, 'pmSkipLog').pmSkipTwoOut.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    outAlt = [];
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
outAlt = single(data) + 10;
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo(mfilename, ...
    'Test fixture: two DATA outputs (input + 1 and input + 10).');
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
    'Input plus one.', 'pmSkipTwoOutData.dat', 1, 'isData', true);
info = PipelineManager.addOutput(info, ...
    'outAlt', 'ImageTimeSeries', 'data', ...
    'Input plus ten.', 'pmSkipTwoOutAlt.dat', 2, 'isData', true);
end
