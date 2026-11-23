function warnIfBaselineNearZero(warningId, stepName, nChecked, nLow, advice)
%WARNIFBASELINENEARZERO Warn once when most traces have a near-zero baseline.
%
%   warnIfBaselineNearZero(warningId, stepName, nChecked, nLow, advice)
%
%   Companion of COUNTLOWBASELINE for steps that divide by a per-trace
%   baseline. The step counts its traces while it processes them (array
%   rows, file slabs, trials), sums the counts, and calls this function once
%   at the end. When at least 20 % of the examined traces are low the input
%   looks zero-centered, so the division is unstable, and a non-blocking
%   warning names the step, the share of affected traces and what to do.
%   The result of the step is not changed.
%
%   Inputs:
%       warningId - Warning identifier, e.g.
%                   'Umitoolbox:normalizeLPF:ZeroCenteredInput'. Users can
%                   silence it with warning('off', warningId).
%       stepName  - Name used in the message, e.g. 'normalizeLPF'.
%       nChecked  - Traces examined (sum of the COUNTLOWBASELINE outputs).
%       nLow      - Traces counted as low.
%       advice    - Step-specific remedy appended to the message.
%
%   Nothing is emitted when no trace was examined or the share is below 20 %.
%
%   Example:
%       warnIfBaselineNearZero('Umitoolbox:normalizeLPF:ZeroCenteredInput', ...
%           'normalizeLPF', nChecked, nLow, 'Set Normalize to false.');
%
%   See also: countLowBaseline

arguments
    warningId {mustBeTextScalar}
    stepName {mustBeTextScalar}
    nChecked (1,1) {mustBeNonnegative}
    nLow (1,1) {mustBeNonnegative}
    advice {mustBeTextScalar} = ''
end

minShare = 0.2;

if nChecked == 0 || nLow / nChecked < minShare
    return
end

message = sprintf(['%s: %.0f%% of the %d traces have a baseline smaller than 3 times ' ...
    'their own fluctuation, so the input looks zero-centered (already ' ...
    'baseline-normalized, z-scored or baseline-subtracted). Dividing by a baseline ' ...
    'near zero makes the result unstable, with values of arbitrary size.'], ...
    char(stepName), 100 * nLow / nChecked, nChecked);
if strlength(string(advice)) > 0
    message = sprintf('%s %s', message, char(advice));
end

warning(char(warningId), '%s', message);
end
