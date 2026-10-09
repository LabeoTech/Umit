function [outData, metaData] = pmSourceInfoAxes(data, SaveFolder)
%PMSOURCEINFOAXES Test fixture: relabel a 3-D input as Y-X-E through metaData.
%
%   Used by TestPipelineManagerSourceInfo (.dat header Phase 6a). The values
%   are returned unchanged; the metaData output sets dimNames {'Y','X','E'},
%   which replaces the inherited Y-X-T axes of the source file.

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
metaData = struct('dimNames', {{'Y', 'X', 'E'}});
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoAxes', 'Test fixture setting the output axes through metaData.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Input data.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addOutput(info, ...
    'outData', 'Image', 'data', ...
    'Input relabelled as event-split images.', 'pmAxes.dat', 1, 'isData', true);
info = PipelineManager.addOutput(info, ...
    'metaData', 'metaData', 'data', ...
    'Output axes.', '', 2, 'isData', false);
end
