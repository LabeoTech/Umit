function eventInfo = normalizeOptionalEventInfoFields(eventInfo, eventInfoIn, errID)
%NORMALIZEOPTIONALEVENTINFOFIELDS Validate and copy optional UMT eventInfo fields.
%
%   eventInfo = normalizeOptionalEventInfoFields(eventInfo, eventInfoIn, errID)
%
%   Copies the optional per-row fields of a UMT eventInfo from EVENTINFOIN
%   to EVENTINFO (whose eventID already defines the number of rows),
%   validating them; raises ERRID on invalid values. Shared by genUMTStruct,
%   appendUMTEventInfo, and validateUMTStruct so every path keeps them:
%       selected    - logical, false for ignored instances (.dat header 8b);
%                     for aggregated rows, true when the slice has data
%       durationSec - ON-to-OFF duration (s), non-negative or NaN (8b)
%       nInstances  - aggregated rows only: number of selected instances
%                     reduced into the slice, non-negative integer (8c)
%   Pass EVENTINFO = struct() with eventID set, or the same struct as
%   EVENTINFOIN to validate in place.

nE = numel(eventInfo.eventID);

if isfield(eventInfoIn, 'selected') && ~isempty(eventInfoIn.selected)
    sel = eventInfoIn.selected;
    if ~(islogical(sel) || isnumeric(sel)) || numel(sel) ~= nE
        error(errID, ...
            'Operation aborted. eventInfo.selected must be a logical vector with one value per row.');
    end
    eventInfo.selected = logical(sel(:));
end

if isfield(eventInfoIn, 'durationSec') && ~isempty(eventInfoIn.durationSec)
    dur = eventInfoIn.durationSec;
    if ~isnumeric(dur) || numel(dur) ~= nE || any(dur(:) < 0)
        error(errID, ...
            ['Operation aborted. eventInfo.durationSec must be a numeric vector ' ...
             'with one non-negative value (or NaN) per row.']);
    end
    eventInfo.durationSec = double(dur(:));
end

if isfield(eventInfoIn, 'nInstances') && ~isempty(eventInfoIn.nInstances)
    nInst = eventInfoIn.nInstances;
    if ~isnumeric(nInst) || numel(nInst) ~= nE || any(nInst(:) < 0) || ...
            any(mod(nInst(:), 1) ~= 0)
        error(errID, ...
            ['Operation aborted. eventInfo.nInstances must be a vector of ' ...
             'non-negative integers with one value per row.']);
    end
    eventInfo.nInstances = double(nInst(:));
end
end
