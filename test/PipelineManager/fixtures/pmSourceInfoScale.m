function outData = pmSourceInfoScale(data, SaveFolder)
%PMSOURCEINFOSCALE Test fixture: multiply the input by 2 and return it in RAM.
%
%   Used by TestPipelineManagerSourceInfo. Accepts an array or a file name
%   (RAM-safe mode passes files).

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

outData = 2 .* localData(data, SaveFolder);
end

function data = localData(data, SaveFolder)
if ischar(data) || (isstring(data) && isscalar(data))
    fileName = char(string(data));
    if ~isfile(fileName)
        fileName = fullfile(SaveFolder, fileName);
    end
    data = loadData(fileName);
end
data = single(data);
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoScale', 'Test fixture multiplying the input by 2.');
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
    'Input multiplied by 2.', 'pmScaled.dat', 1, 'isData', true);
end
