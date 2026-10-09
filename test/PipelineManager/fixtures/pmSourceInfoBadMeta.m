function [outData, metaData] = pmSourceInfoBadMeta(data, SaveFolder)
%PMSOURCEINFOBADMETA Test fixture: return the input and a non-struct metaData.
%
%   Used by TestPipelineManagerSourceInfo. PipelineManager must ignore the
%   metaData output because it is not a struct.

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
outData = single(data);
metaData = 'not a struct';
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoBadMeta', 'Test fixture returning a non-struct metaData.');
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
    'Unchanged input.', 'pmBadMeta.dat', 1, 'isData', true);
info = PipelineManager.addOutput(info, ...
    'metaData', 'metaData', 'data', ...
    'Invalid metadata.', '', 2, 'isData', false);
end
