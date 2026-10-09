function outData = pmSourceInfoMean(data, SaveFolder)
%PMSOURCEINFOMEAN Test fixture: average a Y-X-T input over T (a Y-X image).
%
%   Used by TestPipelineManagerSourceInfo (.dat header Phase 6a). Its output
%   is declared 'Image', so PipelineManager saves it with axes {'Y','X'}
%   although the source file is Y-X-T.

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
outData = mean(single(data), 3);
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoMean', 'Test fixture averaging the input over time.');
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
    'Temporal mean image.', 'pmMean.dat', 1, 'isData', true);
end
