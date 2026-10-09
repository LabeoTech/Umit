function outData = pmSourceInfoEvents(data, SaveFolder)
%PMSOURCEINFOEVENTS Test fixture: split a Y-X-T input into 2 events (Y-X-T-E).
%
%   Used by TestPipelineManagerSourceInfo (.dat header Phase 6a). The T axis
%   (even length) is cut into two halves stacked along a fourth axis; the
%   output is declared 'ImageTimeSeries', so PipelineManager saves it with
%   axes {'Y','X','T','E'}.

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
data = single(data);
nT = size(data, 3) / 2;
outData = cat(4, data(:, :, 1:nT), data(:, :, nT+1:end));
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoEvents', 'Test fixture splitting the input into two events.');
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
    'Event-split data.', 'pmEvents.dat', 1, 'isData', true);
end
