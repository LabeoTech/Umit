function fields = sourceInfoFields()
%SOURCEINFOFIELDS Closed whitelist of per-data metadata fields that can be injected.
%
%   fields = sourceInfoFields()
%
%   Frame rate, axes, and exposure describe the data, not the pipeline.
%   This is the single definition of which source-Info fields a function
%   may receive as an explicit Name-Value parameter, under which name, and
%   what counts as a valid value. PipelineManager ('sourceInfo' inputs,
%   their injection and generated scripts), resolveDataInfoValue, and
%   EventsManager read it.
%
%   Output:
%       fields - 1x3 struct array with fields:
%           field   - source-Info / .dat Info schema field name
%           nvName  - Name-Value parameter name in function calls
%           isValid - function handle, true for a resolved value. NaN or
%                     empty values are unresolved.
%
%       field          nvName          valid value
%       frameRateHz    FrameRateHz     finite numeric scalar > 0
%       dimNames       DimNames        .dat layout: 'Y','X', then distinct
%                                      'T','E','F' in that order
%       exposureMsec   ExposureMsec    finite real numeric scalar
%
%   See also: resolveDataInfoValue, PipelineManager, loadMetaData

fields = struct( ...
    'field', {'frameRateHz', 'dimNames', 'exposureMsec'}, ...
    'nvName', {'FrameRateHz', 'DimNames', 'ExposureMsec'}, ...
    'isValid', {@iValidRate, @iValidLayout, @iValidExposure});
end

function tf = iValidRate(value)
tf = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value) && value > 0;
end

function tf = iValidExposure(value)
tf = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end

function tf = iValidLayout(value)
tf = false;
if ~(iscell(value) || isstring(value)) || isempty(value)
    return
end
try
    names = cellstr(string(value(:).'));
catch
    return
end
if numel(names) < 2 || numel(names) > 5 || ~strcmp(names{1}, 'Y') || ~strcmp(names{2}, 'X')
    return
end
[isKnown, slot] = ismember(names(3:end), {'T', 'E', 'F'});
tf = all(isKnown) && all(diff(slot) > 0);
end
