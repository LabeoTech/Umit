function mapping = resolveDatEventMapping(Info, saveFolder)
%RESOLVEDATEVENTMAPPING Map the E axis of an event-split .dat onto events.mat.
%
%   mapping = resolveDatEventMapping(Info, saveFolder)
%
%   An event-split .dat (layout Y-X-T-E or Y-X-E) stores no event labels.
%   Its E axis is matched to the SaveFolder's events.mat by size:
%     - E = number of event instances: one slice per instance, ignored ones
%       included, in EventsManager split order (.dat header Phase 8b);
%     - E = number of conditions (fewer than instances): one slice per
%       condition, aggregated over events, in order of first appearance
%       (EventsManager.conditionAggregationPlan, Phase 8c);
%     - with one repetition per condition both counts are equal; the file is
%       read per instance, which is the same data.
%   Any other size is not matched.
%
%   Inputs:
%       Info       - .dat Info from loadMetaData (dimNames must contain 'E').
%       saveFolder - Folder holding events.mat (normally the file's folder).
%
%   Output struct:
%       status    - 'matched' (per instance), 'aggregated' (per condition),
%                   'mismatch', or 'noEvents'
%       message   - '' when matched, else a user-facing explanation
%       nE        - size of the E axis
%       nEvents   - number of event instances in events.mat (NaN without)
%       eventInfo - UMT-style eventInfo with one row per E slice:
%                   eventID, repetitionIndex, eventName (string),
%                   eventAxisMode ('instances'), selected, durationSec, and
%                   baselinePeriod when known.
%                   When matched it comes from events.mat. When aggregated
%                   it is the aggregated eventInfo of the conditions, with
%                   nInstances and selected recomputed from the CURRENT
%                   events.mat flags (they may differ from the flags used
%                   when the file was computed). Otherwise all
%                   slices are repetitions 1..nE of one condition ("All
%                   events"), all selected, durations unknown (NaN); the
%                   baselinePeriod of events.mat is kept when it exists.
%
%   See also EventsManager.getEventInstances, DatImageSource.

dimNames = cellstr(string(Info.dimNames));
idxE = find(strcmpi(dimNames, 'E'), 1, 'first');
if isempty(idxE)
    error('Umitoolbox:resolveDatEventMapping:noEventAxis', ...
        'The file "%s" has no E axis.', char(string(Info.filePath)));
end
nE = double(Info.dimSizes(idxE));

mapping = struct('status', 'noEvents', 'message', '', 'nE', nE, ...
    'nEvents', NaN, 'eventInfo', struct());
baselinePeriod = [];

eventsFile = fullfile(saveFolder, 'events.mat');
if ~isfile(eventsFile)
    mapping.message = sprintf(['No events.mat was found next to this file, so its %d ' ...
        'event slices cannot be labeled. All slices are shown as repetitions of ' ...
        'one condition.'], nE);
    mapping.eventInfo = iSingleCondition(nE, baselinePeriod);
    return
end

try
    ev = EventsManager(saveFolder);
    inst = ev.getEventInstances();
    if ~isempty(ev.baselinePeriod)
        baselinePeriod = double(ev.baselinePeriod);
    end
catch ME
    mapping.message = sprintf(['events.mat could not be read (%s). All %d event ' ...
        'slices are shown as repetitions of one condition.'], ME.message, nE);
    mapping.eventInfo = iSingleCondition(nE, baselinePeriod);
    return
end

mapping.nEvents = numel(inst.eventID);
nConditions = numel(EventsManager.conditionOrder(inst.eventID));
if mapping.nEvents ~= nE && nConditions == nE
    instInfo = struct('eventID', inst.eventID, 'repetitionIndex', inst.repetitionIndex, ...
        'eventName', inst.eventName, 'eventAxisMode', 'instances', ...
        'selected', inst.selected, 'durationSec', inst.durationSec);
    if ~isempty(baselinePeriod)
        instInfo.baselinePeriod = baselinePeriod;
    end
    plan = EventsManager.conditionAggregationPlan(instInfo);
    mapping.status = 'aggregated';
    mapping.eventInfo = plan.eventInfoOut;
    return
end
if mapping.nEvents ~= nE
    mapping.status = 'mismatch';
    mapping.message = sprintf(['This file has %d event slices but events.mat lists %d ' ...
        'event instances in %d conditions, so the slices cannot be matched to the ' ...
        'events. All slices are shown as repetitions of one condition. (Event-split ' ...
        'data saved since .dat header Phase 8b keeps every event; older files or ' ...
        'files split from another events.mat may not.)'], nE, mapping.nEvents, nConditions);
    mapping.eventInfo = iSingleCondition(nE, baselinePeriod);
    return
end

eventInfo = struct();
eventInfo.eventID = inst.eventID(:);
eventInfo.repetitionIndex = inst.repetitionIndex(:);
eventInfo.eventName = inst.eventName(:);
eventInfo.eventAxisMode = 'instances';
eventInfo.selected = logical(inst.selected(:));
eventInfo.durationSec = inst.durationSec(:);
if ~isempty(baselinePeriod)
    eventInfo.baselinePeriod = baselinePeriod;
end
mapping.status = 'matched';
mapping.eventInfo = eventInfo;
end

function eventInfo = iSingleCondition(nE, baselinePeriod)
%ISINGLECONDITION All slices as repetitions 1..nE of one condition.
eventInfo = struct();
eventInfo.eventID = ones(nE, 1);
eventInfo.repetitionIndex = (1:nE).';
eventInfo.eventName = repmat("All events", nE, 1);
eventInfo.eventAxisMode = 'instances';
eventInfo.selected = true(nE, 1);
eventInfo.durationSec = nan(nE, 1);
if ~isempty(baselinePeriod)
    eventInfo.baselinePeriod = baselinePeriod;
end
end
