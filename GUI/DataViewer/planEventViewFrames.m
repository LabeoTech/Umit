function [frMat, nLead, evOnsetFrame] = planEventViewFrames(frMat, baselinePeriod, frameRateHz, bAverage)
%PLANEVENTVIEWFRAMES Frame matrix of DataViewer's event view for a continuous .dat.
%
%   [frMat, nLead, evOnsetFrame] = planEventViewFrames(frMat, baselinePeriod, ...
%       frameRateHz, bAverage)
%
%   FRMAT is EventsManager.getFrameMatrix's output for the current condition
%   (NaN = frame outside the recording). The event view applies the same rule
%   as split_data_by_event: the FIRST trial (row 1) may start before the
%   recording; its NLEAD leading baseline frames are kept as NaN so the view
%   fills them with the mean of the baseline frames the trial has (see
%   padEventTrialTrace, eventAverageFrame).
%
%   BAVERAGE true (AVERAGE of the condition): every trial is cropped to the
%   first frame any trial lacks (splitTrialFramePlan). When the plan fails
%   (a trial other than the first starts outside the recording, or the first
%   has no baseline frame inside the recording) the old crop at the first NaN
%   column applies and NLEAD is 0.
%
%   BAVERAGE false (one repetition): FRMAT is not cropped; NLEAD is the
%   number of leading NaN frames when the trial has at least one valid
%   baseline frame, 0 otherwise.
%
%   EVONSETFRAME is round(baselinePeriod * frameRateHz): the baseline frames
%   before the onset, as in splitTrialFramePlan.

evOnsetFrame = round(double(baselinePeriod) * double(frameRateHz));
nLead = 0;

if isempty(frMat)
    return
end

if bAverage
    [frPlanned, nLeadPlanned, ok] = splitTrialFramePlan(frMat, evOnsetFrame, frameRateHz);
    if ok
        frMat = frPlanned;
        nLead = nLeadPlanned;
        return
    end
    firstInvalidCol = find(any(isnan(frMat), 1), 1, 'first');
    if ~isempty(firstInvalidCol)
        frMat = frMat(:, 1:firstInvalidCol - 1);
    end
    return
end

firstValid = find(~isnan(frMat(1, :)), 1, 'first');
if isempty(firstValid) || firstValid == 1
    return
end
% The fill needs a valid baseline frame: one before the onset (column
% evOnsetFrame + 1), inside the recording.
if firstValid - 1 < evOnsetFrame
    nLead = firstValid - 1;
end
end
