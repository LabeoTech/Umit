function outData = pmInjEcho(data, SaveFolder, varargin)
%PMINJECHO Test fixture: record the injected per-data metadata, return the input.
%
%   Used by TestPipelineManagerSourceInfoInjection (.dat header Phase 6b-2).
%   Declares 'sourceInfo' inputs FrameRateHz (required), DimNames and
%   ExposureMsec (optional). Each call appends a record to
%   getappdata(0, 'pmInjLog'): the step function name, the values received
%   ([] when not passed), and the size of the data.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
addParameter(p, 'FrameRateHz', []);
addParameter(p, 'DimNames', []);
addParameter(p, 'ExposureMsec', []);
parse(p, varargin{:});

if ischar(data) || (isstring(data) && isscalar(data))
    fileName = char(string(data));
    if ~isfile(fileName)
        fileName = fullfile(SaveFolder, fileName);
    end
    data = loadData(fileName);
end

rec = struct('fcn', mfilename, 'FrameRateHz', {p.Results.FrameRateHz}, ...
    'DimNames', {p.Results.DimNames}, 'ExposureMsec', {p.Results.ExposureMsec}, ...
    'dataSize', size(data));
log = getappdata(0, 'pmInjLog');
if isempty(log)
    log = rec;
else
    log(end+1) = rec;
end
setappdata(0, 'pmInjLog', log);

outData = single(data);
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo(mfilename, ...
    'Test fixture recording injected per-data metadata.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Input data.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, 'FrameRateHz', 'sourceInfo', ...
    'Frame rate of the input data.', 'kind', 'sourceInfo', 'sourceField', 'frameRateHz');
info = PipelineManager.addInput(info, 'DimNames', 'sourceInfo', ...
    'Axes of the input data.', 'kind', 'sourceInfo', 'sourceField', 'dimNames', ...
    'required', false);
info = PipelineManager.addInput(info, 'ExposureMsec', 'sourceInfo', ...
    'Exposure of the input data.', 'kind', 'sourceInfo', 'sourceField', 'exposureMsec', ...
    'required', false);
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'Unchanged input.', 'pmInjEcho.dat', 1, 'isData', true);
end
