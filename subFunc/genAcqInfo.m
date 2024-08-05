function AcqInfoStream = genAcqInfo(Width,Height,Length,FrameRateHz,ExposureMsec,extraInfo)
% GENACQINFO creates the 'AcqInfoStream' structure necessary to read all
% .dat files inside the saveFolder.

% Set AcqInfoStream structure with basic information:
AcqInfoStream.Width = Width;
AcqInfoStream.Height = Height;
AcqInfoStream.Length = Length;
AcqInfoStream.FrameRateHz = FrameRateHz;
AcqInfoStream.ExposureMsec = ExposureMsec;

if ~exist('extraInfo','var');return;end
% Add extra fields:
fn = fieldnames(extraInfo);
for ii = 1:length(fn)
    AcqInfoStream.(fn{ii}) = extraInfo.(fn{ii});
end

end