function outData = pmLeafMismatchedFile(data, SaveFolder)
%PMLEAFMISMATCHEDFILE Test fixture returning a non-declared backing file.

declaredFile = 'declared.dat';
producedFile = 'actual.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo(declaredFile);
    return
end

info = [];
if ischar(data) || (isstring(data) && isscalar(data))
    [data, info] = loadData(fullfile(SaveFolder, char(string(data))));
end
% Explicit axes and rate (.dat header Phase 6a): the input file's, else
% the test folder's Y-X-T at 20 Hz.
saveArgs = {'DimNames', {'Y', 'X', 'T'}, 'FrameRateHz', 20};
if ~isempty(info)
    saveArgs = {'DimNames', info.dimNames, 'Info', info};
end
saveData(fullfile(SaveFolder, producedFile), data, saveArgs{:});
outData = producedFile;
end

function info = localPipelineInfo(declaredFile)
info = PipelineManager.createPipelineInfo( ...
    'pmLeafMismatchedFile', ...
    'Test fixture returning a permanent file different from its declaration.');
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
    'File-backed output.', declaredFile, 1, 'isData', true);
end
