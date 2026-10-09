function signal = makePulseSignal(nSamples, pulseStarts, pulseWidth, varargin)
%MAKEPULSESIGNAL Create a synthetic 1-D pulse train.
%
% Syntax:
%   signal = makePulseSignal(nSamples, pulseStarts, pulseWidth)
%   signal = makePulseSignal(..., 'Amplitude', 5, 'Baseline', 0, ...)
%
% Inputs:
%   nSamples    - Total number of samples.
%   pulseStarts - Vector of 1-based pulse start indices.
%   pulseWidth  - Scalar pulse width in samples.
%
% Name-Value Pairs:
%   'Amplitude' - Pulse amplitude. Default = 5.
%   'Baseline'  - Baseline level. Default = 0.
%   'NoiseStd'  - Additive Gaussian noise std. Default = 0.
%   'Polarity'  - 1 for positive pulses, -1 for negative pulses.
%
% Output:
%   signal - Column vector.

    p = inputParser;
    addRequired(p, 'nSamples', @(x) isnumeric(x) && isscalar(x) && x == round(x) && x > 0);
    addRequired(p, 'pulseStarts', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'pulseWidth', @(x) isnumeric(x) && isscalar(x) && x == round(x) && x >= 1);
    addParameter(p, 'Amplitude', 5, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'Baseline', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NoiseStd', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Polarity', 1, @(x) isnumeric(x) && isscalar(x) && any(x == [-1 1]));
    parse(p, nSamples, pulseStarts, pulseWidth, varargin{:});

    signal = ones(nSamples, 1, 'single') .* single(p.Results.Baseline);
    pulseStarts = unique(p.Results.pulseStarts(:)');
    for ii = 1:numel(pulseStarts)
        idx1 = max(1, pulseStarts(ii));
        idx2 = min(nSamples, pulseStarts(ii) + p.Results.pulseWidth - 1);
        signal(idx1:idx2) = single(p.Results.Baseline + p.Results.Polarity * p.Results.Amplitude);
    end

    if p.Results.NoiseStd > 0
        signal = signal + single(p.Results.NoiseStd .* randn(size(signal), 'single'));
    end
end
