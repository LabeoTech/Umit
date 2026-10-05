function AcqInfoStream = genAcqInfo(Width,Height,Length,FrameRateHz,ExposureMsec,extraInfo)
% GENACQINFO creates the 'AcqInfoStream' structure saved in AcqInfos.mat.
%
% AcqInfos.mat describes the raw acquisition: pass the raw (pre-binning)
% Width, Height, FrameRateHz, and ExposureMsec. The imported data is
% described by each .dat header. Pass an empty Length to leave the field
% out, which importers do (.dat header Phase 7b).

% Set AcqInfoStream structure with basic information:
AcqInfoStream.Width = Width;
AcqInfoStream.Height = Height;
if ~isempty(Length)
    AcqInfoStream.Length = Length;
end
AcqInfoStream.FrameRateHz = FrameRateHz;
AcqInfoStream.ExposureMsec = ExposureMsec;

if ~exist('extraInfo','var');return;end
% Add extra fields:
fn = fieldnames(extraInfo);
for ii = 1:length(fn)
    AcqInfoStream.(fn{ii}) = extraInfo.(fn{ii});
end

end
