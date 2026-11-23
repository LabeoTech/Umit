function [frame, padFrame] = eventAverageFrame(getFrameFcn, frMat, col, nLead, evOnsetFrame, padFrame)
%EVENTAVERAGEFRAME Mean image of one column of the event frame matrix.
%
%   [frame, padFrame] = eventAverageFrame(getFrameFcn, frMat, col, nLead, ...
%       evOnsetFrame, padFrame)
%
%   Averages the source frames of column COL over the trials (rows of FRMAT)
%   that have one. When the first trial starts before the recording (NLEAD
%   leading NaN frames) and COL is one of them, that trial contributes the
%   mean image of its valid baseline frames NLEAD+1 to EVONSETFRAME, as
%   split_data_by_event fills it (padLeadingBaseline).
%
%   GETFRAMEFCN reads one source frame (e.g. @(k) app.DataSource.getFrame(k)).
%   PADFRAME is that baseline mean image, computed on first use and returned
%   so the caller can cache it across columns ([] to compute it).
%
%   Returns an empty FRAME when no trial has a frame at COL.

col = min(max(col, 1), size(frMat, 2));
frames = frMat(:, col);
useFill = nLead >= 1 && col <= nLead && isnan(frames(1));

sourceFrames = frames(isfinite(frames));
frame = [];
n = 0;
for k = 1:numel(sourceFrames)
    thisFrame = single(getFrameFcn(sourceFrames(k)));
    if isempty(frame)
        frame = zeros(size(thisFrame), 'single');
    end
    n = n + 1;
    frame = frame + (thisFrame - frame) ./ n;
end

if useFill
    if isempty(padFrame)
        padFrame = iBaselineMeanFrame(getFrameFcn, frMat(1, :), nLead, evOnsetFrame);
    end
    if ~isempty(padFrame)
        if isempty(frame)
            frame = padFrame;
        else
            n = n + 1;
            frame = frame + (padFrame - frame) ./ n;
        end
    end
end
end

function padFrame = iBaselineMeanFrame(getFrameFcn, row, nLead, evOnsetFrame)
%IBASELINEMEANFRAME Mean of the valid baseline frames of the first trial.

padFrame = [];
baseCols = nLead + 1:min(evOnsetFrame, numel(row));
baseFrames = row(baseCols);
baseFrames = baseFrames(isfinite(baseFrames));
n = 0;
for k = 1:numel(baseFrames)
    thisFrame = single(getFrameFcn(baseFrames(k)));
    if isempty(padFrame)
        padFrame = zeros(size(thisFrame), 'single');
    end
    n = n + 1;
    padFrame = padFrame + (thisFrame - padFrame) ./ n;
end
end
