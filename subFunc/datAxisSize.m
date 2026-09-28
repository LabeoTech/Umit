function n = datAxisSize(Info, axisName)
%DATAXISSIZE Size of a named axis of a .dat file, or 0 when absent.
%
%   n = datAxisSize(Info, axisName)
%
%   Looks up axisName in Info.dimNames and returns the matching entry of
%   Info.dimSizes. Returns 0 when the file has no such axis.
%
%   Inputs:
%       Info     - .dat metadata from loadMetaData (needs dimNames and
%                  dimSizes).
%       axisName - Axis name, e.g. 'Y', 'X', 'T', 'E'.
%
%   Output:
%       n        - Axis size (double), 0 when absent.
%
%   Axis names are not checked against the header vocabulary, so legacy
%   axis orders such as E-Y-X-T are supported.
%
%   Example:
%       Info = loadMetaData('green.dat');
%       nFrames = datAxisSize(Info, 'T');
%
%   See also: loadMetaData

if ~isstruct(Info) || ~isscalar(Info) || ~isfield(Info, 'dimNames') || ...
        ~isfield(Info, 'dimSizes')
    error('Umitoolbox:datAxisSize:invalidInput', ...
        'Info must be a scalar struct with fields dimNames and dimSizes.');
end
if ~(ischar(axisName) && (isrow(axisName) || isempty(axisName))) && ...
        ~(isstring(axisName) && isscalar(axisName))
    error('Umitoolbox:datAxisSize:invalidInput', ...
        'axisName must be a character vector or string scalar.');
end

idx = find(strcmp(cellstr(Info.dimNames), char(axisName)), 1);
if isempty(idx)
    n = 0;
else
    n = double(Info.dimSizes(idx));
end
end
