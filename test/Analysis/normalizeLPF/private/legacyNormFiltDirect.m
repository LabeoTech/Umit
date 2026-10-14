function OutData = legacyNormFiltDirect(OutData, lowFreq, highFreq, bDivide, bExpFit, Freq)
%LEGACYNORMFILTDIRECT Test oracle: the legacy in-RAM normalization filter.
%
%   OutData = legacyNormFiltDirect(OutData, lowFreq, highFreq, bDivide, bExpFit, Freq)
%
%   Frozen copy of NormFiltDirect from IOIAnalysis/NormalisationFiltering.m
%   as it was before the core was absorbed into Analysis/Filters/normalizeLPF.m
%   (2026-10-13). The normalizeLPF tests compare against it, so the
%   reference stays the legacy algorithm. Only the frame-rate inference
%   from a folder (never used by the tests) was removed: FREQ is required.
%   Do not change this file to follow normalizeLPF.

ExpFun = @(P,x) abs(P(1)).*exp(-abs(P(2)).*x) + abs(P(3)).*exp(-abs(P(4)).*x) ...
    - abs(P(5))*x + P(6);
Opt = optimset(@fminsearch);
Opt.Display = 'off';

% Temporal filtering
if lowFreq > 0
    UseLPFilt = 1;
    if (1/lowFreq) > (size(OutData,3)/Freq)
        lowFreq = 1 / (2 * (size(OutData,3)/Freq));
    end
    f = fdesign.lowpass('N,F3dB', 4, lowFreq, Freq);
    lpass = design(f, 'butter');
else
    UseLPFilt = 0;
end

if highFreq > 0
    UseHPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, highFreq, Freq);
    hpass = design(f, 'butter');
else
    UseHPFilt = 0;
end

dims = size(OutData);
if bExpFit
    rng('shuffle');
    S = mean(reshape(OutData, [], dims(3)), 1);
    initB = double(rand(1,6) .* [30 1 20 1 1 mean(double(S))]);
    B = fminsearch(@(P) sum((double(S) - ExpFun(P, (1:size(S,2)))).^2), ...
        initB, Opt);
    Approx = ExpFun([B(1:4) 0 0], 1:size(S,2));
    Pred = [ones(1, size(S,2)); linspace(0,1,size(S,2)); Approx]';
end

for ind = 1:dims(1)
    Signal = double(squeeze(OutData(ind,:,:)));

    if bExpFit
        B = Pred \ Signal';
        Approx = (Pred * B)';
        Signal = Signal ./ Approx;
    end

    if UseLPFilt
        LP_lowCutOff = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
    else
        LP_lowCutOff = ones(size(Signal));
    end

    if UseHPFilt
        LP_highCutOff = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
    else
        LP_highCutOff = Signal;
    end

    if bDivide
        OutData(ind,:,:) = single(LP_highCutOff ./ LP_lowCutOff);
    else
        OutData(ind,:,:) = single(LP_highCutOff - LP_lowCutOff);
    end
end
end
