function outData = pmInjTwoInputs(data, ref, SaveFolder, varargin)
%PMINJTWOINPUTS Test fixture: per-data metadata from a non-primary input.
%
%   Used by TestPipelineManagerSourceInfoInjection (.dat header Phase 6b-2).
%   FrameRateHz comes from the second data input 'ref' (sourceInput), and
%   ExposureMsec from the primary input 'data'. Each call appends a record
%   to getappdata(0, 'pmInjLog') and returns DATA unchanged.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
addParameter(p, 'FrameRateHz', []);
addParameter(p, 'ExposureMsec', []);
parse(p, varargin{:});

data = localLoad(data, SaveFolder);
ref = localLoad(ref, SaveFolder);

rec = struct('fcn', mfilename, 'FrameRateHz', {p.Results.FrameRateHz}, ...
    'DimNames', {[]}, 'ExposureMsec', {p.Results.ExposureMsec}, ...
    'dataSize', size(data));
log = getappdata(0, 'pmInjLog');
if isempty(log)
    log = rec;
else
    log(end+1) = rec;
end
setappdata(0, 'pmInjLog', log);

outData = single(data) + 0 * mean(ref(:));
end

function v = localLoad(v, SaveFolder)
if ischar(v) || (isstring(v) && isscalar(v))
    fileName = char(string(v));
    if ~isfile(fileName)
        fileName = fullfile(SaveFolder, fileName);
    end
    v = loadData(fileName);
end
v = single(v);
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo(mfilename, ...
    'Test fixture taking per-data metadata from two inputs.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Primary input.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'ref', 'ImageTimeSeries', 'Reference input.', ...
    'position', 2, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 3, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, 'FrameRateHz', 'sourceInfo', ...
    'Frame rate of the reference input.', 'kind', 'sourceInfo', ...
    'sourceField', 'frameRateHz', 'sourceInput', 'ref');
info = PipelineManager.addInput(info, 'ExposureMsec', 'sourceInfo', ...
    'Exposure of the primary input.', 'kind', 'sourceInfo', ...
    'sourceField', 'exposureMsec', 'required', false);
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'Primary input, unchanged.', 'pmInjTwo.dat', 1, 'isData', true);
end
