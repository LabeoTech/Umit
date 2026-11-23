function trace = padEventTrialTrace(trace, nLead, evOnsetFrame)
%PADEVENTTRIALTRACE Fill the leading baseline of a first-trial trace with its baseline mean.
%
%   trace = padEventTrialTrace(trace, nLead, evOnsetFrame)
%
%   TRACE is a vector (any orientation) whose first NLEAD samples are missing
%   (NaN) because the trial starts before the recording. They are replaced by
%   the mean of the valid baseline samples NLEAD+1 to EVONSETFRAME, the same
%   rule split_data_by_event applies to the images (padLeadingBaseline).
%
%   NLEAD = 0, or a trace without a valid baseline sample, returns TRACE
%   unchanged (NaN stays NaN).

if nLead < 1 || isempty(trace)
    return
end
lastBaseline = min(evOnsetFrame, numel(trace));
if lastBaseline <= nLead
    return
end
baselineMean = mean(double(trace(nLead + 1:lastBaseline)), 'omitnan');
trace(1:nLead) = baselineMean;
end
