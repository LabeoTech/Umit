function [outData, metaData] = pmSourceInfoRate(data, SaveFolder)
%PMSOURCEINFORATE Test fixture: return the input and a metaData update.
%
%   Used by TestPipelineManagerSourceInfo. metaData sets frameRateHz = 5
%   and exposureMsec = 9, and carries a field PipelineManager must ignore.

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
outData = single(data) + 1;
metaData = struct('frameRateHz', 5, 'exposureMsec', 9, 'channelName', 'ignored');
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoRate', 'Test fixture updating the frame rate through metaData.');
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
    'Input plus 1.', 'pmRated.dat', 1, 'isData', true);
info = PipelineManager.addOutput(info, ...
    'metaData', 'metaData', 'data', ...
    'Updated frame rate and exposure.', '', 2, 'isData', false);
end
