function hdr = datHeaderFromInfo(Info, channelName, varargin)
%DATHEADERFROMINFO Header description of an output written from a source Info.
%
%   hdr = datHeaderFromInfo(Info, channelName)
%   hdr = datHeaderFromInfo(Info, channelName, 'Name', Value, ...)
%
%   Returns the header description that spatialSlabIO('create', file, hdr)
%   and encodeDatHeader expect, for an output derived from the .dat whose
%   Info (as returned by loadMetaData) is given: the output inherits the
%   source's class, axes, sizes, frame rate, and exposure, and carries
%   CHANNELNAME (normally the output file's base name).
%
%   Inputs:
%       Info        - .dat Info struct with fields dataClass, dimNames,
%                     dimSizes, and optionally frameRateHz and exposureMsec
%                     (NaN when absent).
%       channelName - Text. Characters outside printable ASCII are replaced
%                     by '_', and the name is truncated to the header field
%                     (datHeaderSchema channelNameMaxChars), silently.
%
%   Name-Value options (override the value taken from Info):
%       dataClass, dimNames, dimSizes, frameRateHz, exposureMsec
%
%   Output:
%       hdr - Struct with fields dataClass, frameRateHz, exposureMsec,
%             channelName, dimNames, dimSizes. It is not validated here;
%             spatialSlabIO('create') and encodeDatHeader validate it.
%
%   See also: spatialSlabIO, encodeDatHeader, saveData, loadMetaData

p = inputParser;
p.FunctionName = 'datHeaderFromInfo';
addRequired(p, 'Info', @(x) isstruct(x) && isscalar(x));
addRequired(p, 'channelName', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'dataClass', []);
addParameter(p, 'dimNames', []);
addParameter(p, 'dimSizes', []);
addParameter(p, 'frameRateHz', []);
addParameter(p, 'exposureMsec', []);
parse(p, Info, channelName, varargin{:});
opts = p.Results;

required = {'dataClass', 'dimNames', 'dimSizes'};
for k = 1:numel(required)
    if ~isfield(Info, required{k}) && isempty(opts.(required{k}))
        error('Umitoolbox:datHeaderFromInfo:invalidInfo', ...
            'Info must contain dataClass, dimNames, and dimSizes (missing %s).', required{k});
    end
end

hdr = struct();
hdr.dataClass = char(string(iPick(opts.dataClass, Info, 'dataClass', '')));
hdr.frameRateHz = double(iPick(opts.frameRateHz, Info, 'frameRateHz', NaN));
hdr.exposureMsec = double(iPick(opts.exposureMsec, Info, 'exposureMsec', NaN));
hdr.channelName = iSanitizeName(channelName);
hdr.dimNames = cellstr(string(iPick(opts.dimNames, Info, 'dimNames', {})));
hdr.dimSizes = double(iPick(opts.dimSizes, Info, 'dimSizes', []));
end

function value = iPick(override, Info, fieldName, default)
if ~isempty(override)
    value = override;
elseif isfield(Info, fieldName) && ~isempty(Info.(fieldName))
    value = Info.(fieldName);
else
    value = default;
end
end

function name = iSanitizeName(name)
[~, codes] = datHeaderSchema(1);
name = char(name);
name(name < 32 | name > 126) = '_';
name = name(1:min(end, codes.constants.channelNameMaxChars));
end
