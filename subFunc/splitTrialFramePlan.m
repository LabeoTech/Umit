function [frMat, nLead, ok, msg] = splitTrialFramePlan(frMat, evOnsetFrame, frameRateHz)
%SPLITTRIALFRAMEPLAN Crop trial frames and detect a first trial that starts before the recording.
%
%   [frMat, nLead, ok, msg] = splitTrialFramePlan(frMat, evOnsetFrame, frameRateHz)
%
%   FRMAT is the repetition-by-frame matrix of EventsManager.getFrameMatrix
%   (NaN = frame outside the recording). Every trial is cropped to the first
%   column that any trial lacks, as event splitting always did, with one
%   exception: the FIRST trial (row 1) may start before the recording. Its
%   leading NaN frames (NLEAD of them) are baseline frames the recording does
%   not have; they are kept in FRMAT (still NaN) so the caller can fill them
%   (see padLeadingBaseline) and do not count as a crop point.
%
%   Inputs:
%       frMat        - Trial-by-frame matrix of frame indices (NaN outside).
%       evOnsetFrame - Number of baseline frames before the onset of a trial,
%                      round(baselinePeriod * frameRateHz); the onset is at
%                      column evOnsetFrame + 1.
%       frameRateHz  - Data frame rate, only used for the message.
%
%   Outputs:
%       frMat  - Cropped matrix (leading NaN of row 1 kept).
%       nLead  - Leading NaN frames of the first trial (0 when none).
%       ok     - False when the trials cannot be split: no trial, a trial
%                other than the first one starting outside the recording, or
%                a first trial without any valid baseline frame.
%       msg    - Warning text when NLEAD > 0, '' otherwise.
%
%   A trial with leading NaN needs at least one valid baseline frame (a
%   frame of the recording before its onset), to give the fill a value.

nLead = 0;
ok = false;
msg = '';

if isempty(frMat)
    return
end

firstValid = find(~isnan(frMat(1, :)), 1, 'first');
if isempty(firstValid)
    return                                  % first trial wholly outside
end
if any(isnan(frMat(2:end, 1)))
    return                                  % another trial starts outside
end
nLead = firstValid - 1;

nanMask = isnan(frMat);
nanMask(1, 1:nLead) = false;                % the padded frames are not a crop point
firstNaNCol = find(any(nanMask, 1), 1, 'first');
if ~isempty(firstNaNCol)
    frMat = frMat(:, 1:firstNaNCol - 1);
end

if nLead > 0
    nValid = min(evOnsetFrame, size(frMat, 2)) - nLead;
    if nValid < 1
        nLead = 0;
        return                              % no baseline frame to fill from
    end
    msg = sprintf(['The first event instance starts %d frames (%.3g s) before the ' ...
        'recording, so its baseline period is incomplete. The %d missing baseline ' ...
        'frames of that trial were filled with the mean of its %d valid baseline ' ...
        'frame(s). To avoid this, reduce the baselinePeriod parameter in the ' ...
        'EventsManager (setBaselinePeriod) and run this function again.'], ...
        nLead, nLead / frameRateHz, nLead, nValid);
end
ok = true;
end
