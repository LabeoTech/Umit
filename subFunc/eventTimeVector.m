function t = eventTimeVector(nSamples, frameRateHz, baselinePeriod)
%EVENTTIMEVECTOR Event-locked time axis of one trial, in seconds.
%
%   t = eventTimeVector(nSamples, frameRateHz, baselinePeriod)
%
%   Returns the times (s, relative to the event onset) of the NSAMPLES frames
%   of one event-split trial. Time 0 is the frame EventsManager uses as the
%   onset: column round(baselinePeriod * frameRateHz) of each getFrameMatrix
%   row (and of splitDataByEvents' T axis). Use this helper wherever an
%   event-split trace is plotted or exported so the time axis cannot drift
%   from the split (.dat header Phase 8b).
%
%   Inputs:
%       nSamples       - Number of frames of the trial (positive integer).
%       frameRateHz    - Frame rate of the data (positive scalar).
%       baselinePeriod - Pre-event baseline (s, non-negative). Empty or NaN
%                        gives a time axis starting at 0.
%
%   Output:
%       t - Row vector [1, nSamples].
%
%   See also EventsManager.getFrameMatrix.

validateattributes(nSamples, {'numeric'}, {'scalar', 'integer', 'positive'}, ...
    'eventTimeVector', 'nSamples');
validateattributes(frameRateHz, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
    'eventTimeVector', 'frameRateHz');

frameRateHz = double(frameRateHz);
if nargin < 3 || isempty(baselinePeriod) || ~isfinite(baselinePeriod)
    t = (0:double(nSamples) - 1) ./ frameRateHz;
    return
end
validateattributes(baselinePeriod, {'numeric'}, {'scalar', 'real', 'nonnegative'}, ...
    'eventTimeVector', 'baselinePeriod');

onsetColumn = round(double(baselinePeriod) * frameRateHz);
t = ((1:double(nSamples)) - onsetColumn) ./ frameRateHz;
end
