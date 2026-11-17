function result = saveEventSelection(evObj, saveFolder, opts)
%SAVEEVENTSELECTION Write the loaded event selection to events.mat.
%
%   result = saveEventSelection(evObj, saveFolder)
%   result = saveEventSelection(evObj, saveFolder, 'ParentFigure', fig)
%
%   With 'ParentFigure' (DataViewer's figure) a uiconfirm asks before events.mat
%   is overwritten; Cancel leaves the file unchanged. Without it nothing is
%   asked (tests).
%
%   DataViewer's Delete Condition, Delete Repetition and Restore change the
%   selection (which events are ignored) of the loaded EventsManager object
%   only. This is the explicit save: it writes the object to SAVEFOLDER's
%   events.mat with EventsManager.saveEvents, which keeps every other
%   variable (detected events, condition names, analog signal, baseline).
%   Event-split data follow the saved selection: an ignored repetition stays
%   in the file, flagged "(ignored)" and left out of AVERAGE.
%
%   RESULT has fields changed (true when events.mat was written), message and
%   severity ('info', 'warning' or 'error'). When the file already holds the
%   same selection nothing is written. Errors are reported in RESULT, never
%   rethrown.

arguments
    evObj
    saveFolder
    opts.ParentFigure = []
end

result = struct('changed', false, 'message', '', 'severity', 'info');
try
    if ~isa(evObj, 'EventsManager') || ~isscalar(evObj) || ~isvalid(evObj)
        error('Umitoolbox:saveEventSelection:noEvents', ...
            'There is no event source to save.');
    end
    saveFolder = char(string(saveFolder));
    if isempty(saveFolder) || ~isfolder(saveFolder)
        error('Umitoolbox:saveEventSelection:noFolder', ...
            'The SaveFolder of the data is not accessible: "%s".', saveFolder);
    end

    eventsFile = fullfile(saveFolder, 'events.mat');
    if ~isfile(eventsFile)
        result.message = ['This folder has no events.mat: detect the events with the ' ...
            'Events Manager first.'];
        result.severity = 'warning';
        return
    end

    saved = load(eventsFile, 'selectedEvents');
    if isfield(saved, 'selectedEvents') && ...
            isequal(logical(saved.selectedEvents(:)), logical(evObj.selectedEvents(:)))
        result.message = 'The event selection is already saved in events.mat.';
        return
    end

    nIgnored = nnz(evObj.state & ~evObj.selectedEvents);
    if ~isempty(opts.ParentFigure) && isvalid(opts.ParentFigure)
        choice = uiconfirm(opts.ParentFigure, sprintf(['Write the event selection (%d ignored ' ...
            'repetition(s)) to events.mat of this folder?\n\nIgnored trials stay excluded from ' ...
            'DataViewer, from the average of a condition and from analysis functions that use ' ...
            'events.mat, until they are restored and saved again.'], nIgnored), ...
            'Save event selection', 'Options', {'Save', 'Cancel'}, ...
            'DefaultOption', 2, 'CancelOption', 2, 'Icon', 'question');
        if ~strcmp(choice, 'Save')
            result.message = 'Saving the event selection was cancelled; events.mat is unchanged.';
            return
        end
    end

    evObj.saveEvents(saveFolder);
    result.changed = true;
    result.message = sprintf(['Event selection saved to events.mat (%d ignored ' ...
        'repetition(s)).'], nIgnored);
catch ME
    result.message = sprintf('Could not save the event selection: %s', ME.message);
    result.severity = 'error';
end
end
