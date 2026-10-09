function outData = pmLegacyMetaDataEcho(data, metaData, SaveFolder)
%PMLEGACYMETADATAECHO Test fixture recording the metaData input it receives.
%
%   Declares a non-data input named "metaData". Since .dat header Phase 8a
%   PipelineManager rejects this (MetaDataInputUnsupported); the fixture
%   is kept to test that. If it ever runs, it appends the received value to
%   getappdata(0, 'pmLegacyMetaDataEcho') and returns the data unchanged.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

received = getappdata(0, 'pmLegacyMetaDataEcho');
if isempty(received)
    received = {};
end
received{end+1} = metaData;
setappdata(0, 'pmLegacyMetaDataEcho', received);

if ischar(data) || (isstring(data) && isscalar(data))
    data = loadData(fullfile(SaveFolder, char(string(data))));
end
outData = data;
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmLegacyMetaDataEcho', ...
    'Test fixture recording the metaData input it receives.');
info = PipelineManager.addInput(info, ...
    'data', 'ImageTimeSeries', 'Input data.', ...
    'position', 1, 'callType', 'positional', ...
    'isData', true, 'supportsFile', true, 'dataMode', 'either');
info = PipelineManager.addInput(info, ...
    'metaData', 'metaData', 'Metadata of the input data.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 3, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'Unchanged input data.', 'pmEcho.dat', 1, 'isData', true);
end
