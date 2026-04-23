function outData = calculateResponseFeatures(data, varargin)
%CALCULATERESPONSEFEATURES Extract response features from ROI event data.
%
%   outData = calculateResponseFeatures(data)
%   outData = calculateResponseFeatures(data, 'STD_threshold', 1.5, ...)
%
%   This function calculates response features from ROI-organized temporal
%   data split by event. It is intended to operate on the roi UMT output of
%   getDataFromROI when the selected entry uses dimensions {'ROI','T','E'}.
%
%   Inputs:
%       data - ROI UMT structure.
%
%   Name-Value parameters:
%       STD_threshold    - Positive scalar multiplier used to define the
%                          response threshold from the baseline standard
%                          deviation. Default: 1
%       ResponsePolarity - 'positive' or 'negative'. Default: 'positive'
%       TimeWindow_sec   - 'all' or [start end] in seconds relative to
%                          stimulus onset. Default: 'all'
%
%   Output:
%       outData - ROI UMT structure with one entry using dimensions
%                 {'ROI','Measure','E'}.
%
%   Notes:
%       - The input must have kind = 'roi'.
%       - The selected entry must use dimensions {'ROI','T','E'}.
%       - Shared top-level eventInfo must exist and use
%         eventAxisMode = 'instances'.
%       - baselinePeriod must exist in top-level eventInfo.
%       - FrameRateHz must exist in the selected entry meta.
%       - ROI entries that retain the Pixel dimension are not supported by
%         this function in the current version.


dependency = 'getDataFromROI'; %#ok<NASGU>
default_Output = 'RespFeatures.umt'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addParameter(p, 'STD_threshold', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ResponsePolarity', 'positive', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'positive','negative'}));
addParameter(p, 'TimeWindow_sec', 'all', @(x) iValidateTimeWindowInput(x));

parse(p, data, varargin{:});

STD_threshold = double(p.Results.STD_threshold);
responsePolarity = lower(char(string(p.Results.ResponsePolarity)));
timeWindowSec = p.Results.TimeWindow_sec;

assert(isstruct(data) && isscalar(data), ...
    'Umitoolbox:calculateResponseFeatures:InvalidInputType', ...
    'Input data must be a scalar roi UMT struct.');

validateUMTStruct(data, 'requireEventInfo', true);

assert(strcmpi(char(string(data.kind)), 'roi'), ...
    'Umitoolbox:calculateResponseFeatures:InvalidUMTKind', ...
    'Input UMT must have kind = "roi".');

assert(isfield(data, 'eventInfo') && isstruct(data.eventInfo), ...
    'Umitoolbox:calculateResponseFeatures:MissingEventInfo', ...
    'Input roi UMT must contain shared top-level eventInfo.');

assert(strcmpi(char(string(data.eventInfo.eventAxisMode)), 'instances'), ...
    'Umitoolbox:calculateResponseFeatures:InvalidEventAxisMode', ...
    'Input eventInfo.eventAxisMode must be "instances".');

assert(isfield(data.eventInfo, 'baselinePeriod') && ~isempty(data.eventInfo.baselinePeriod), ...
    'Umitoolbox:calculateResponseFeatures:MissingBaselinePeriod', ...
    'Input eventInfo must contain baselinePeriod.');

entryNames = fieldnames(data.data);
assert(~isempty(entryNames), ...
    'Umitoolbox:calculateResponseFeatures:EmptyInput', ...
    'Input roi UMT contains no data entries.');

entryName = entryNames{1};
entry = data.data.(entryName);
dimNames = cellstr(string(entry.dimNames));

assert(~any(strcmp(dimNames, 'Pixel')), ...
    'Umitoolbox:calculateResponseFeatures:PixelDimensionNotSupported', ...
    ['Inputs that retain the Pixel dimension are not supported in the ' ...
     'current version. Aggregate ROI pixels first in getDataFromROI.']);

assert(isequal(dimNames, {'ROI','T','E'}), ...
    'Umitoolbox:calculateResponseFeatures:InvalidInputDims', ...
    'Input ROI entry must use dimensions {''ROI'',''T'',''E''}.');

assert(isfield(entry, 'meta') && isstruct(entry.meta) && isfield(entry.meta, 'FrameRateHz') && ...
    ~isempty(entry.meta.FrameRateHz), ...
    'Umitoolbox:calculateResponseFeatures:MissingFrameRate', ...
    'Input ROI entry meta must contain FrameRateHz.');

frameRateHz = double(entry.meta.FrameRateHz);
baselinePeriodSec = double(data.eventInfo.baselinePeriod);
roiData = single(entry.value);

nROI = size(roiData, 1);
nT = size(roiData, 2);
nE = size(roiData, 3);

baselineFrames = round(baselinePeriodSec * frameRateHz);
baselineFrames = max(1, min(baselineFrames, nT));

if ischar(timeWindowSec) || (isstring(timeWindowSec) && isscalar(timeWindowSec))
    frOn = baselineFrames + 1;
    frOff = nT;
else
    timeWindowSec = double(timeWindowSec(:).');
    frOn = round(timeWindowSec(1) * frameRateHz) + baselineFrames;
    frOff = round(timeWindowSec(2) * frameRateHz) + baselineFrames;

    assert(frOn >= baselineFrames + 1 && frOff <= nT && frOn <= frOff, ...
        'Umitoolbox:calculateResponseFeatures:InvalidTimeWindow', ...
        ['TimeWindow_sec is out of range for the available response ' ...
         'duration.']);
end

featureNames = { ...
    'PeakAmplitude', ...
    'PeakLatency', ...
    'TimeWindowAverageAmplitude', ...
    'AUCamplitude', ...
    'OnsetAmplitude', ...
    'OnsetLatency'};

nMeasure = numel(featureNames);
outVals = nan(nROI, nMeasure, nE, 'single');

for iEvent = 1:nE
    avgResp = roiData(:, :, iEvent);

    if strcmp(responsePolarity, 'negative')
        avgResp = -1 .* avgResp;
    end

    avgBsln = mean(avgResp(:, 1:baselineFrames), 2, 'omitnan');
    thr = avgBsln + STD_threshold .* std(avgResp(:, 1:baselineFrames), 0, 2, 'omitnan');

    windowResp = avgResp(:, frOn:frOff);
    [peakValue, peakIdxLocal] = max(windowResp, [], 2);
    peakIdx = frOn + peakIdxLocal - 1;

    avgAmp = mean(windowResp, 2, 'omitnan') - avgBsln;
    aucAmp = trapz(windowResp, 2) - trapz(avgResp(:, 1:baselineFrames), 2);
    peakAmp = peakValue - avgBsln;

    onsetAmp = nan(nROI, 1, 'single');
    onsetLat = nan(nROI, 1, 'single');
    peakLat = nan(nROI, 1, 'single');

    for iROI = 1:nROI
        if peakValue(iROI) > thr(iROI)
            onsetLocal = find(avgResp(iROI, frOn:frOff) > thr(iROI), 1, 'first');
            if ~isempty(onsetLocal)
                onsetIdx = frOn + onsetLocal - 1;
                onsetAmp(iROI) = avgResp(iROI, onsetIdx) - avgBsln(iROI);
                onsetLat(iROI) = (onsetIdx - baselineFrames) / frameRateHz;
                peakLat(iROI) = (peakIdx(iROI) - baselineFrames) / frameRateHz;
            end
        end
    end

    outVals(:, 1, iEvent) = peakAmp;
    outVals(:, 2, iEvent) = peakLat;
    outVals(:, 3, iEvent) = avgAmp;
    outVals(:, 4, iEvent) = aucAmp;
    outVals(:, 5, iEvent) = onsetAmp;
    outVals(:, 6, iEvent) = onsetLat;
end

labels = struct();
if isfield(data, 'labels') && isfield(data.labels, 'ROI')
    labels.ROI = data.labels.ROI;
else
    labels.ROI = arrayfun(@num2str, 1:nROI, 'UniformOutput', false);
end
labels.Measure = featureNames;

outData = genUMTStruct( ...
    outVals, ...
    'kind', 'roi', ...
    'entryName', entryName, ...
    'dimNames', {'ROI','Measure','E'}, ...
    'labels', labels, ...
    'eventInfo', data.eventInfo);

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Calculate ROI response features from roi UMT event data.');

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ProcessedData'}, ...
            'ROI UMT input with dimensions {ROI,T,E}.', ...
            'kind', 'input', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', false, ...
            'dataMode', 'either');

        info = PipelineManager.addInput( ...
            info, ...
            'STD_threshold', ...
            'parameter', ...
            'Baseline standard-deviation multiplier used for thresholding.', ...
            'kind', 'parameter', ...
            'default', 1, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'ResponsePolarity', ...
            'parameter', ...
            'Response polarity: positive or negative.', ...
            'kind', 'parameter', ...
            'default', 'positive', ...
            'allowed', {'positive','negative'}, ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'TimeWindow_sec', ...
            'parameter', ...
            'Response time window in seconds or ''all''.', ...
            'kind', 'parameter', ...
            'default', 'all', ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'ROI feature UMT output.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

function tf = iValidateTimeWindowInput(x)
%IVALIDATETIMEWINDOWINPUT Validate TimeWindow_sec input.

tf = false;

if ischar(x) || (isstring(x) && isscalar(x))
    tf = strcmpi(char(string(x)), 'all');
    return
end

if isnumeric(x) && isvector(x) && numel(x) == 2
    x = double(x(:).');
    tf = all(isfinite(x)) && all(x >= 0) && x(1) <= x(2);
end
end
