function validateDatHeader(hdr, varargin)
%VALIDATEDATHEADER Check a .dat header description against the version-1 rules.
%
%   validateDatHeader(hdr)
%   validateDatHeader(hdr, 'FileBytes', n)
%   validateDatHeader(hdr, 'ReservedBytes', bytes)
%
%   Errors with identifier Umitoolbox:validateDatHeader:invalidHeader and a
%   message naming the problem when the descriptive fields of HDR break a
%   rule of docs/dev/dat-header-spec.md. Returns nothing when HDR is valid.
%
%   Input:
%       hdr - Scalar struct using the field names of decodeDatHeader.
%             Required fields: dataClass, frameRateHz, exposureMsec,
%             channelName, dimNames, dimSizes. Optional fields used by the
%             'FileBytes' check: writeComplete (default false) and
%             dataOffset (default: the version-1 header length).
%
%   Name-Value arguments:
%       'FileBytes'     - Total file size in bytes. Must equal
%                         dataOffset + prod(dimSizes) * bytesPerValue. A
%                         mismatch errors (Umitoolbox:validateDatHeader:
%                         fileSizeMismatch) when writeComplete is true and
%                         warns (Umitoolbox:validateDatHeader:incompleteFile)
%                         when it is false.
%       'ReservedBytes' - Bytes of the channel-metadata reserved area. Any
%                         nonzero byte warns (Umitoolbox:validateDatHeader:
%                         reservedBytesNotZero).
%
%   Rules checked (version 1):
%       - dataClass has a data class code.
%       - dimNames and dimSizes have the same length.
%       - Every name in dimNames is an assigned axis (Y, X, T, E, F).
%       - Y and X are present.
%       - No axis appears more than once.
%       - dimNames is in slot order (Y, X, T, E, F), the data's memory order.
%       - Every size is an integer between 1 and intmax('uint32').
%       - With a T axis, frameRateHz is finite and greater than 0.
%       - Without a T axis, frameRateHz is NaN.
%
%   See also: encodeDatHeader, decodeDatHeader, datHeaderSchema

p = inputParser;
addParameter(p, 'FileBytes', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'ReservedBytes', [], @(x) isempty(x) || isnumeric(x));
parse(p, varargin{:});
fileBytes = p.Results.FileBytes;
reservedBytes = p.Results.ReservedBytes;

[~, codes] = datHeaderSchema(1);

if ~isstruct(hdr) || ~isscalar(hdr)
    iFail('The header must be a scalar struct.');
end

requiredFields = {'dataClass', 'frameRateHz', 'exposureMsec', ...
    'channelName', 'dimNames', 'dimSizes'};
missingFields = setdiff(requiredFields, fieldnames(hdr), 'stable');
if ~isempty(missingFields)
    iFail(sprintf('Missing header field(s): %s.', strjoin(missingFields, ', ')));
end

% --- data class --------------------------------------------------------
dataClass = hdr.dataClass;
if isstring(dataClass) && isscalar(dataClass)
    dataClass = char(dataClass);
end
classIdx = [];
if ischar(dataClass)
    classIdx = find(strcmp({codes.dataClass.className}, dataClass), 1);
end
if isempty(classIdx)
    iFail(sprintf('dataClass "%s" has no valid data class code.', iToText(dataClass)));
end

% --- dimensions --------------------------------------------------------
dimNames = hdr.dimNames;
if isstring(dimNames)
    dimNames = cellstr(dimNames);
end
if ~iscellstr(dimNames)
    iFail('dimNames must be a cell array of axis names.');
end
dimNames = dimNames(:).';

dimSizes = hdr.dimSizes;
if ~isnumeric(dimSizes) || ~(isvector(dimSizes) || isempty(dimSizes))
    iFail('dimSizes must be a numeric vector.');
end
dimSizes = double(dimSizes(:).');

if numel(dimNames) ~= numel(dimSizes)
    iFail(sprintf('dimNames has %d entries but dimSizes has %d.', ...
        numel(dimNames), numel(dimSizes)));
end

assignedMask = ~[codes.axis.reserved];
assignedNames = {codes.axis(assignedMask).name};
assignedSlots = [codes.axis(assignedMask).slot];

[isAssigned, nameIdx] = ismember(dimNames, assignedNames);
if ~all(isAssigned)
    iFail(sprintf('"%s" is not an assigned axis. Assigned axes: %s.', ...
        dimNames{find(~isAssigned, 1)}, strjoin(assignedNames, ', ')));
end

for requiredAxis = {'Y', 'X'}
    if ~any(strcmp(dimNames, requiredAxis{1}))
        iFail(sprintf('Axis %s is missing; every .dat file needs both Y and X.', ...
            requiredAxis{1}));
    end
end

[uniqueNames, ~, whichUnique] = unique(dimNames, 'stable');
if numel(uniqueNames) ~= numel(dimNames)
    counts = accumarray(whichUnique(:), 1);
    iFail(sprintf('Axis %s appears more than once.', uniqueNames{find(counts > 1, 1)}));
end

slots = assignedSlots(nameIdx);
if any(diff(slots) <= 0)
    iFail(sprintf(['dimNames {%s} is not in slot order (%s), which is the ' ...
        'memory order of the data.'], strjoin(dimNames, ','), strjoin(assignedNames, ', ')));
end

maxSize = double(intmax('uint32'));
badSize = ~isreal(dimSizes) | ~isfinite(dimSizes) | dimSizes ~= round(dimSizes) | ...
    dimSizes < 1 | dimSizes > maxSize;
if any(badSize)
    k = find(badSize, 1);
    iFail(sprintf('Size of axis %s is %s; sizes must be integers from 1 to %d.', ...
        dimNames{k}, num2str(dimSizes(k)), maxSize));
end

% --- frame rate --------------------------------------------------------
frameRateHz = hdr.frameRateHz;
if ~isnumeric(frameRateHz) || ~isscalar(frameRateHz) || ~isreal(frameRateHz)
    iFail('frameRateHz must be a real numeric scalar.');
end
hasT = any(strcmp(dimNames, 'T'));
if hasT && ~(isfinite(frameRateHz) && frameRateHz > 0)
    iFail(sprintf(['The file has a T axis, so frameRateHz must be finite and ' ...
        'greater than 0 (got %s).'], num2str(frameRateHz)));
end
if ~hasT && ~isnan(frameRateHz)
    iFail(sprintf('The file has no T axis, so frameRateHz must be NaN (got %s).', ...
        num2str(frameRateHz)));
end

if ~isnumeric(hdr.exposureMsec) || ~isscalar(hdr.exposureMsec) || ~isreal(hdr.exposureMsec)
    iFail('exposureMsec must be a real numeric scalar (NaN when unknown).');
end

% --- optional checks ---------------------------------------------------
if ~isempty(reservedBytes) && any(reservedBytes(:) ~= 0)
    warning('Umitoolbox:validateDatHeader:reservedBytesNotZero', ...
        'The channel-metadata reserved area contains %d nonzero byte(s).', ...
        nnz(reservedBytes));
end

if ~isempty(fileBytes)
    writeComplete = isfield(hdr, 'writeComplete') && logical(hdr.writeComplete);
    if isfield(hdr, 'dataOffset') && ~isempty(hdr.dataOffset)
        dataOffset = double(hdr.dataOffset);
    else
        dataOffset = codes.constants.headerLength;
    end
    expectedBytes = dataOffset + prod(dimSizes) * codes.dataClass(classIdx).bytesPerValue;
    if double(fileBytes) ~= expectedBytes
        msg = sprintf(['File size is %d bytes, but the header describes %d bytes ' ...
            '(%d header + %s data).'], double(fileBytes), expectedBytes, dataOffset, ...
            iToText(expectedBytes - dataOffset));
        if writeComplete
            error('Umitoolbox:validateDatHeader:fileSizeMismatch', '%s', msg);
        else
            warning('Umitoolbox:validateDatHeader:incompleteFile', ...
                '%s The write-complete flag is not set, so the file may be incomplete.', msg);
        end
    end
end

end

% =========================================================================
function iFail(msg)
error('Umitoolbox:validateDatHeader:invalidHeader', 'Invalid .dat header: %s', msg);
end

function txt = iToText(value)
if ischar(value)
    txt = value;
elseif isstring(value) && isscalar(value)
    txt = char(value);
elseif isnumeric(value) && isscalar(value)
    txt = num2str(value);
else
    txt = sprintf('<%s>', class(value));
end
end
