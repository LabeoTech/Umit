function out = computeAggregates(vals,aggfcn,idxDim)
% computeAggregates - Aggregates values along a specified dimension using a selected aggregation function.
%
% Syntax: out = computeAggregates(vals, aggfcn, idxDim)
%
% Inputs:
%   vals - A numerical array of values to be aggregated.
%   aggfcn - A string specifying the aggregation function ('mean', 'median', 'mode', 'std', 'max', 'min', 'sum').
%   idxDim - An integer specifying the dimension along which to aggregate.
%
% Outputs:
%   out - The aggregated result based on the selected function and dimension.
%
% Example:
%   vals = [1, 2, 3; 4, 5, 6];
%   aggfcn = 'mean';
%   idxDim = 2;
%   out = computeAggregates(vals, aggfcn, idxDim)
%   % This returns the mean of each row.
%
% Supported aggregation functions:
%   - 'mean': Arithmetic mean, ignoring NaNs.
%   - 'median': Median value, ignoring NaNs.
%   - 'mode': Mode (most frequent value).
%   - 'std': Standard deviation, ignoring NaNs.
%   - 'max': Maximum value, ignoring NaNs.
%   - 'min': Minimum value, ignoring NaNs.
%   - 'sum': Sum of values, ignoring NaNs.
%
% If an unsupported function is provided, the input values are returned unchanged, and a warning is issued.
%
% See also: mean, median, mode, std, max, min, sum
%
% Note: NaNs are ignored in calculations where applicable.

switch aggfcn
    case 'mean'
        out = mean(vals, idxDim,'omitnan');
    case 'median'
        out = median(vals, idxDim, 'omitnan');
    case 'mode'
        out = mode(vals, idxDim);
    case 'std'
        out = std(vals, 0, idxDim, 'omitnan');
    case 'max'
        out = max(vals, [], idxDim, 'omitnan');
    case 'min'
        out = min(vals, [], idxDim, 'omitnan');
    case 'sum'
        out = sum(vals, idxDim, 'omitnan');
    otherwise
        warning('Unsupported aggregation function "%s". Returning original values.', aggfcn);
        out = vals;
end
end


switch aggfcn
    case 'mean'
        out = mean(vals, idxDim,'omitnan');
    case 'median'
        out = median(vals, idxDim, 'omitnan');
    case 'mode'
        out = mode(vals, idxDim);
    case 'std'
        out = std(vals, 0, idxDim, 'omitnan');
    case 'max'
        out = max(vals, [], idxDim, 'omitnan');
    case 'min'
        out = min(vals, [], idxDim, 'omitnan');
    case 'sum'
        out = sum(vals, idxDim, 'omitnan');
    otherwise
        warning('Unsupported aggreation function "%s". Returning original values.' aggfcn);
        out = vals;
end
end
