function outData = pmSourceInfoMake(SaveFolder)
%PMSOURCEINFOMAKE Test fixture: create a 2 x 3 x 4 array with no data input.
%
%   Used by TestPipelineManagerSourceInfo: its output has no source file,
%   so saveData falls back to AcqInfos.mat for the frame rate.

if nargin == 1 && (ischar(SaveFolder) || (isstring(SaveFolder) && isscalar(SaveFolder))) && ...
        strcmpi(strtrim(char(string(SaveFolder))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

outData = reshape(single(1:24), [2 3 4]);
end

function info = localPipelineInfo()
info = PipelineManager.createPipelineInfo( ...
    'pmSourceInfoMake', 'Test fixture creating data without a data input.');
info = PipelineManager.addInput(info, ...
    'SaveFolder', 'SaveFolder', 'Output folder.', ...
    'kind', 'input', 'position', 1, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addOutput(info, ...
    'outData', 'ImageTimeSeries', 'data', ...
    'Created data.', 'pmMade.dat', 1, 'isData', true);
end
