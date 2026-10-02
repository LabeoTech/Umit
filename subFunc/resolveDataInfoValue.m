function value = resolveDataInfoValue(field, explicitValue, data, fnName, varargin)
%RESOLVEDATAINFOVALUE Per-data metadata value with the explicit > header > error precedence.
%
%   value = resolveDataInfoValue(field, explicitValue, data, fnName)
%   value = resolveDataInfoValue(..., 'OwnValue', v, 'OwnSource', text)
%
%   Shared rule for functions that take frame rate, axes, or exposure as an
%   explicit Name-Value parameter (FrameRateHz, DimNames, ExposureMsec; see
%   sourceInfoFields). PipelineManager injects these parameters from the
%   data flowing into the step; direct callers pass them by hand.
%   AcqInfos.mat is never used.
%
%   Inputs:
%       field         - sourceInfoFields field: 'frameRateHz', 'dimNames',
%                       or 'exposureMsec'.
%       explicitValue - Value of the function's Name-Value parameter ([]
%                       when not given).
%       data          - The function's data input. When it is a .dat file
%                       name (char or string), its header is the data's own
%                       metadata.
%       fnName        - Function name, used in identifiers and messages.
%
%   Name-Value options:
%       OwnValue  - The data's own metadata value when the data is not a
%                   .dat file, e.g. a .umt entry's meta.FrameRateHz. Used
%                   in place of the header. Default: [] (none).
%       OwnSource - Text naming where OwnValue comes from, for messages.
%
%   Precedence:
%       1) A valid explicitValue. If the data's own value (header or
%          OwnValue) is valid and different, warning
%          Umitoolbox:<fnName>:sourceInfoConflict; the explicit value wins.
%       2) The data's own value (.dat header, or OwnValue), when valid.
%       3) Error Umitoolbox:<fnName>:missing<NVName>, telling direct
%          callers which Name-Value pair to pass.
%   An invalid explicit value raises Umitoolbox:<fnName>:invalid<NVName>.
%
%   See also: sourceInfoFields, loadMetaData

p = inputParser;
p.FunctionName = 'resolveDataInfoValue';
addRequired(p, 'field', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'explicitValue');
addRequired(p, 'data');
addRequired(p, 'fnName', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && strlength(string(x)) > 0);
addParameter(p, 'OwnValue', []);
addParameter(p, 'OwnSource', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, field, explicitValue, data, fnName, varargin{:});

field = char(p.Results.field);
fnName = char(p.Results.fnName);
spec = sourceInfoFields();
idx = find(strcmp({spec.field}, field), 1);
if isempty(idx)
    error('Umitoolbox:resolveDataInfoValue:invalidInput', ...
        'Unknown source-Info field "%s". Allowed: %s.', field, strjoin({spec.field}, ', '));
end
spec = spec(idx);

[ownValue, ownSource] = iOwnValue(data, field, spec, p.Results.OwnValue, char(p.Results.OwnSource));

if ~iIsUnset(explicitValue)
    if ~spec.isValid(explicitValue)
        error(sprintf('Umitoolbox:%s:invalid%s', fnName, spec.nvName), ...
            '%s: ''%s'' has an invalid value.', fnName, spec.nvName);
    end
    value = iNormalize(explicitValue, field);
    if ~isempty(ownValue) && ~isequal(value, ownValue)
        warning(sprintf('Umitoolbox:%s:sourceInfoConflict', fnName), ...
            ['%s: the explicit ''%s'' (%s) differs from %s (%s). The explicit ' ...
             'value is used.'], fnName, spec.nvName, iText(value), ownSource, iText(ownValue));
    end
    return
end

if ~isempty(ownValue)
    value = ownValue;
    return
end

error(sprintf('Umitoolbox:%s:missing%s', fnName, spec.nvName), ...
    ['%s: the %s of the input data is unknown. Pass it as ''%s'', value ' ...
     '(PipelineManager injects it automatically); AcqInfos.mat is not used.'], ...
    fnName, iDescribe(field), spec.nvName);
end

function [value, source] = iOwnValue(data, field, spec, ownValueIn, ownSourceIn)
%IOWNVALUE The data's own value: the .dat header, else OwnValue; [] if none.
value = [];
source = '';
if (ischar(data) || (isstring(data) && isscalar(data))) && ...
        endsWith(char(data), '.dat', 'IgnoreCase', true) && isfile(char(data))
    Info = loadMetaData(char(data));
    if isfield(Info, field) && spec.isValid(Info.(field))
        value = iNormalize(Info.(field), field);
        source = sprintf('the header of "%s"', char(data));
    end
    return
end
if ~iIsUnset(ownValueIn) && spec.isValid(ownValueIn)
    value = iNormalize(ownValueIn, field);
    source = ownSourceIn;
    if isempty(source)
        source = 'the data''s own metadata';
    end
end
end

function tf = iIsUnset(value)
tf = isempty(value) || (isnumeric(value) && isscalar(value) && isnan(value));
end

function value = iNormalize(value, field)
if strcmp(field, 'dimNames')
    value = cellstr(string(value(:).'));
else
    value = double(value);
end
end

function text = iText(value)
if iscell(value)
    text = ['{' strjoin(value, ',') '}'];
else
    text = sprintf('%g', value);
end
end

function text = iDescribe(field)
switch field
    case 'frameRateHz'
        text = 'frame rate';
    case 'dimNames'
        text = 'axis layout';
    otherwise
        text = 'exposure';
end
end
