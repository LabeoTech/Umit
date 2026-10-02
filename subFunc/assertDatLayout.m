function assertDatLayout(Info, acceptedLayouts, fnName)
%ASSERTDATLAYOUT Refuse a .dat input whose axes a function does not support.
%
%   assertDatLayout(Info, acceptedLayouts, fnName)
%
%   Consumers check the layout of every .dat input before reading data or
%   writing outputs (.dat header Phase 6b-1): since writers emit any
%   supported layout (Y-X, Y-X-T, Y-X-T-E, Y-X-E, Y-X-F), a function that
%   assumes Y-X-T would otherwise mis-shape other layouts silently.
%
%   Inputs:
%       Info            - .dat Info from loadMetaData (or a spatialSlabIO
%                         handle's Info): uses dimNames and filePath.
%       acceptedLayouts - Cell array of axis-name lists, e.g.
%                         {{'Y','X','T'}} or {{'Y','X','T'}, {'Y','X'}}.
%       fnName          - Name of the calling function, used in the error
%                         identifier and message.
%
%   Raises Umitoolbox:<fnName>:unsupportedLayout when Info.dimNames equals
%   none of acceptedLayouts. The message names the file, its axes, and the
%   accepted layouts; for event-split files it says that event-split .dat
%   input is not supported by that function yet.
%
%   See also: loadMetaData, datAxisSize

p = inputParser;
p.FunctionName = 'assertDatLayout';
addRequired(p, 'Info', @(x) isstruct(x) && isscalar(x) && isfield(x, 'dimNames'));
addRequired(p, 'acceptedLayouts', @(x) iscell(x) && ~isempty(x) && ...
    all(cellfun(@(c) iscell(c) || isstring(c), x)));
addRequired(p, 'fnName', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && strlength(string(x)) > 0);
parse(p, Info, acceptedLayouts, fnName);

fnName = char(fnName);
dimNames = cellstr(string(Info.dimNames(:).'));
for k = 1:numel(acceptedLayouts)
    if isequal(dimNames, cellstr(string(acceptedLayouts{k}(:).')))
        return
    end
end

fileText = '';
if isfield(Info, 'filePath') && ~isempty(Info.filePath)
    fileText = sprintf(' "%s"', char(string(Info.filePath)));
end
acceptedText = strjoin(cellfun(@(c) ['{' strjoin(cellstr(string(c(:).')), ',') '}'], ...
    acceptedLayouts, 'UniformOutput', false), ' or ');
eventText = '';
if any(strcmp(dimNames, 'E'))
    eventText = ' Event-split .dat input is not supported by this function yet.';
end
error(sprintf('Umitoolbox:%s:unsupportedLayout', fnName), ...
    '%s: the .dat input%s has axes {%s}; supported layouts: %s.%s', ...
    fnName, fileText, strjoin(dimNames, ','), acceptedText, eventText);
end
