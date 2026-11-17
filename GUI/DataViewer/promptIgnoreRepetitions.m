function result = promptIgnoreRepetitions(parentFigure, evObj, conditionName, opts)
%PROMPTIGNOREREPETITIONS Ignore or restore single repetitions of one condition.
%
%   result = promptIgnoreRepetitions(parentFigure, evObj, conditionName)
%   result = promptIgnoreRepetitions(..., 'Interactive', false, ...
%       'Ignored', repetitionIndices)
%
%   Shows the repetitions of CONDITIONNAME (EventsManager.listRepetitions:
%   index and onset time) with a checkbox "Ignore" per repetition, the
%   currently ignored ones ticked. Apply ignores the newly ticked repetitions
%   (EventsManager.removeRepetition) and includes the unticked ones again
%   (EventsManager.restoreRepetition); nothing else is touched. The change is
%   in the EventsManager object only: the caller saves events.mat (the Events
%   Manager's Save events.mat button, DataViewer's File > Save Event
%   Selection).
%
%   RESULT has fields changed (true when any repetition changed state),
%   message and severity ('info' or 'error'). Errors are shown in a uialert
%   and reported in RESULT, never rethrown.
%
%   Name-Value options:
%       Interactive - false: no dialog (tests); Ignored then names the
%                     repetitions that must end up ignored. Default true.
%       Ignored     - with Interactive false, the repetition indices that must
%                     be ignored afterwards (every other one is included).

arguments
    parentFigure
    evObj (1,1) EventsManager
    conditionName (1,:) char
    opts.Interactive (1,1) logical = true
    opts.Ignored (1,:) double = []
    opts.RequireIncludedEvents (1,1) logical = false
end

result = struct('changed', false, 'message', '', 'severity', 'info');
try
    repetitions = evObj.listRepetitions(conditionName);
    if isempty(repetitions)
        result.message = sprintf('Condition "%s" has no repetitions.', conditionName);
        return
    end
    % With RequireIncludedEvents, ignoring every repetition is refused when no
    % other condition still has included events (it would leave nothing).
    mustKeepOne = opts.RequireIncludedEvents && ~evObj.hasOtherSelectedEvents(conditionName);
    keepOneMessage = ['At least one repetition must stay included: ignoring every repetition ' ...
        'of the only condition with included events would leave no events.'];

    if opts.Interactive
        [ignoredAfter, accepted] = iAsk(parentFigure, conditionName, repetitions, ...
            mustKeepOne, keepOneMessage);
        if ~accepted
            result.message = 'Ignore / restore repetitions cancelled.';
            return
        end
    else
        ignoredAfter = unique(opts.Ignored);
        unknown = setdiff(ignoredAfter, repetitions.Repetition');
        if ~isempty(unknown)
            error('Umitoolbox:promptIgnoreRepetitions:unknownRepetition', ...
                'Condition "%s" has no repetition %s.', conditionName, num2str(unknown));
        end
    end

    if mustKeepOne && all(ismember(repetitions.Repetition', ignoredAfter))
        result.message = keepOneMessage;
        result.severity = 'warning';
        return
    end

    ignoredBefore = repetitions.Repetition(repetitions.Ignored)';
    toIgnore = setdiff(ignoredAfter, ignoredBefore);
    toRestore = setdiff(ignoredBefore, ignoredAfter);
    if isempty(toIgnore) && isempty(toRestore)
        result.message = sprintf('Condition "%s": no repetition changed.', conditionName);
        return
    end

    % Restore first: ignoring a repetition of a condition that is ignored
    % as a whole needs at least one included repetition.
    if ~isempty(toRestore)
        evObj.restoreRepetition(conditionName, toRestore);
    end
    if ~isempty(toIgnore)
        evObj.removeRepetition(conditionName, toIgnore);
    end

    result.changed = true;
    result.message = sprintf(['Condition "%s": %d repetition(s) ignored, %d included again. ' ...
        'Save events.mat to keep this change.'], conditionName, numel(toIgnore), numel(toRestore));
    if numel(ignoredAfter) == height(repetitions)
        result.message = [result.message ' All repetitions are ignored: the whole condition is ignored.'];
    end
catch ME
    result.message = sprintf('Could not ignore or restore the repetitions: %s', ME.message);
    result.severity = 'error';
    if opts.Interactive
        uialert(parentFigure, ME.message, 'Ignore repetitions failed', 'Icon', 'error');
    end
end
end

function [ignoredAfter, accepted] = iAsk(parentFigure, conditionName, repetitions, ...
    mustKeepOne, keepOneMessage)
%IASK Modal table with an Ignore checkbox per repetition.
ignoredAfter = [];
accepted = false;

dialog = uifigure('Name', ['Ignore repetitions: ' conditionName], ...
    'WindowStyle', 'modal', 'Position', [100 100 440 420], 'Resize', 'off');
closeDialog = onCleanup(@() iDeleteIfValid(dialog));
if ~isempty(parentFigure) && isvalid(parentFigure) && ...
        exist('placeAppInsideCaller', 'file') == 2
    try
        placeAppInsideCaller(parentFigure, dialog, 'center');
    catch
        movegui(dialog, 'center');
    end
else
    movegui(dialog, 'center');
end
try
    initializeUmitTheme(dialog);
catch
    % Theme is cosmetic.
end

grid = uigridlayout(dialog, [3 2], 'RowHeight', {'fit', '1x', 30}, ...
    'ColumnWidth', {'1x', '1x'});
message = uilabel(grid, 'Text', sprintf(['Condition "%s". Tick the repetitions to ignore ' ...
    '(for example a faulty trigger, or a trial you decided not to use); untick to ' ...
    'include them again. Both the start and the end of a repetition are affected.'], ...
    conditionName), 'WordWrap', 'on');
message.Layout.Column = [1 2];
data = table(repetitions.Repetition, repetitions.Onset_s, repetitions.Ignored, ...
    'VariableNames', {'Repetition', 'Onset_s', 'Ignore'});
tbl = uitable(grid, 'Data', data, 'ColumnEditable', [false false true], ...
    'ColumnName', {'Repetition', 'Onset (s)', 'Ignore'}, 'RowName', {});
tbl.Layout.Column = [1 2];
uibutton(grid, 'Text', 'Apply', 'ButtonPushedFcn', @(~, ~) iClose(true));
uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) iClose(false));
dialog.CloseRequestFcn = @(~, ~) iClose(false);

uiwait(dialog);
% Delete now: the nested-function callbacks keep this workspace (and so
% closeDialog) alive after return, and a lingering modal dialog would hide
% the caller's next dialog.
iDeleteIfValid(dialog);

    function iClose(isAccepted)
        accepted = isAccepted;
        if isAccepted
            current = tbl.Data;
            ignoredAfter = current.Repetition(logical(current.Ignore))';
            if mustKeepOne && all(logical(current.Ignore))
                % Stay in the dialog: the user unticks at least one repetition.
                ignoredAfter = [];
                accepted = false;
                uialert(dialog, keepOneMessage, 'Cannot ignore every repetition', ...
                    'Icon', 'warning');
                return
            end
        end
        uiresume(dialog);
    end
end

function iDeleteIfValid(handle)
if isvalid(handle)
    delete(handle);
end
end
