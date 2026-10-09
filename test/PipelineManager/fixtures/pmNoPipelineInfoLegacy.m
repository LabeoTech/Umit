function outData = pmNoPipelineInfoLegacy(data, metaData, SaveFolder, opts) %#ok<INUSD>
%PMNOPIPELINEINFOLEGACY Test fixture: a legacy, signature-style function.
%
%   It does not answer pmNoPipelineInfoLegacy('pipelineInfo'). Since .dat
%   header Phase 8a PipelineManager skips such functions with the warning
%   PipelineManager:createFcnList:NoPipelineInfo.
%
% default_Output = 'pmNoPipelineInfoLegacy.dat'

outData = data;
end
