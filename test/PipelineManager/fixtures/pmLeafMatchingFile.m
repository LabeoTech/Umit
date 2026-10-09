function outData = pmLeafMatchingFile(data, SaveFolder)
%PMLEAFMATCHINGFILE Test fixture returning its declared backing file.

outputFile = 'declared.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo(outputFile);
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
saveData(fullfile(SaveFolder, outputFile), data, saveArgs{:});
outData = outputFile;
end

function info = localPipelineInfo(outputFile)
info = PipelineManager.createPipelineInfo( ...
    'pmLeafMatchingFile', ...
    'Test fixture returning the permanent file declared for its DATA output.');
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
    'File-backed output.', outputFile, 1, 'isData', true);
end
