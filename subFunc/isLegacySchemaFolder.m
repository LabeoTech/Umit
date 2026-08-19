function [isLegacy, message] = isLegacySchemaFolder(saveFolder)
%ISLEGACYSCHEMAFOLDER Detect whether a SaveFolder holds legacy-schema metadata.
%
%   isLegacy = isLegacySchemaFolder(saveFolder)
%   [isLegacy, message] = isLegacySchemaFolder(saveFolder)
%
%   A legacy-schema dataset is one whose AcqInfos.mat predates dev's
%   centralized AcqInfoStream.ImportedChannels registry (see
%   appendImportedChannelInfo), which is the field the current pipeline
%   path relies on for per-.dat-file / per-illumination-channel metadata
%   (e.g. FrameRateHz). Older (Astrocyte-era) datasets carried this kind
%   of metadata in per-.dat sidecar .mat files instead (see
%   loadMetaData's iLoadLegacySidecar), or never wrote it at all.
%
%   Reinterpreted fields such as FrameRate (overall camera rate under the
%   legacy convention vs. per-channel rate under the current one) cannot
%   be safely read under current semantics for a legacy-schema dataset.
%   PipelineManager and DataViewer use this check to block/disable
%   processing while leaving viewing/inspection unaffected.
%
%   This check only inspects AcqInfos.mat field names (no .dat data is
%   read), so it is cheap to run at pipeline start-up and in UI refresh
%   callbacks.
%
%   isLegacy is true and message explains why whenever the centralized
%   ImportedChannels registry is not available, including when
%   AcqInfos.mat itself is missing.

acqInfoPath = fullfile(saveFolder, 'AcqInfos.mat');

if ~isfile(acqInfoPath)
    isLegacy = true;
    message = iBuildMessage('This dataset has no AcqInfos.mat metadata file (pre-centralization import).');
    return
end

S = load(acqInfoPath, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream') || ~isstruct(S.AcqInfoStream) || ~isscalar(S.AcqInfoStream)
    isLegacy = true;
    message = iBuildMessage('This dataset''s AcqInfos.mat does not contain a valid AcqInfoStream.');
    return
end

if ~isfield(S.AcqInfoStream, 'ImportedChannels')
    isLegacy = true;
    message = iBuildMessage(['This dataset uses a legacy metadata schema (its AcqInfos.mat predates ' ...
        'the centralized per-channel ImportedChannels registry).']);
    return
end

isLegacy = false;
message = '';

end

function message = iBuildMessage(prefix)
%IBUILDMESSAGE Compose the shared legacy-schema explanation text.

message = [prefix, ' Processing is not supported for legacy-schema datasets, because ', ...
    'reinterpreted fields (e.g. FrameRate) cannot be safely reused under current per-channel ', ...
    'semantics. Reprocess this dataset from raw/continuous data to enable pipeline processing. ', ...
    'Viewing and inspection remain fully available.'];

end
