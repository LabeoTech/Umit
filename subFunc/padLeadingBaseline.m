function trial = padLeadingBaseline(trial, nLead, evOnsetFrame)
%PADLEADINGBASELINE Fill the missing leading baseline frames of a trial with its baseline mean.
%
%   trial = padLeadingBaseline(trial, nLead, evOnsetFrame)
%
%   TRIAL is a Y-by-X-by-T array whose first NLEAD frames are missing (NaN).
%   They are replaced, pixel by pixel, by the mean of the valid baseline
%   frames of the trial: frames NLEAD+1 to EVONSETFRAME (the baseline
%   columns that exist; the onset is column EVONSETFRAME + 1). See
%   splitTrialFramePlan, which guarantees at least one such frame.
%
%   NLEAD = 0 returns TRIAL unchanged.

if nLead < 1
    return
end
lastBaseline = min(evOnsetFrame, size(trial, 3));
assert(lastBaseline > nLead, 'Umitoolbox:padLeadingBaseline:noBaselineFrame', ...
    'The trial has no valid baseline frame to fill from.');

baselineMean = mean(trial(:, :, nLead + 1:lastBaseline), 3, 'omitnan');
trial(:, :, 1:nLead) = repmat(baselineMean, 1, 1, nLead);
end
