function [nChecked, nLow] = countLowBaseline(baseline, fluctuation)
%COUNTLOWBASELINE Count traces whose baseline is small against their fluctuation.
%
%   [nChecked, nLow] = countLowBaseline(baseline, fluctuation)
%
%   A step that divides by a per-trace baseline (DeltaR/R) is meaningful
%   only for intensity-like data, whose baseline is much larger than its
%   own fluctuation. Zero-centered data (already baseline-normalized,
%   z-scored or baseline-subtracted) have a baseline within the noise, and
%   dividing by it gives results of arbitrary size.
%
%   A trace counts as low when abs(baseline) < 3 * fluctuation. On the
%   MultiChannel_SmallData fixture the ratio abs(baseline)/fluctuation is
%   32 to 74 (5th to 95th percentile) for raw data and below 1 for the
%   normalized, z-scored and baseline-subtracted variants, so the factor
%   3 leaves a wide margin on both sides.
%
%   Inputs:
%       baseline    - Numeric array, one baseline per trace (any shape).
%       fluctuation - Numeric array of the same size, the temporal spread
%                     of each trace (for example its standard deviation).
%
%   Outputs:
%       nChecked - Number of traces examined. Traces with a non-finite
%                  baseline or fluctuation, and masked traces (baseline
%                  and fluctuation both 0), are not examined.
%       nLow     - Number of examined traces counted as low.
%
%   Example:
%       [nChecked, nLow] = countLowBaseline(mean(x, 2), std(x, 0, 2));
%
%   See also: warnIfBaselineNearZero

arguments
    baseline {mustBeNumeric}
    fluctuation {mustBeNumeric}
end

if ~isequal(size(baseline), size(fluctuation))
    error('Umitoolbox:countLowBaseline:sizeMismatch', ...
        'baseline and fluctuation must have the same size.');
end

ratio = 3;

baseline = double(baseline(:));
fluctuation = double(fluctuation(:));

isExamined = isfinite(baseline) & isfinite(fluctuation) & ...
    ~(baseline == 0 & fluctuation == 0);

nChecked = nnz(isExamined);
nLow = nnz(isExamined & abs(baseline) < ratio .* fluctuation);
end
