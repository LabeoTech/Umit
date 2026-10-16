function [context, store] = bindFolderToProjectSession(dataFolder, selection)
%BINDFOLDERTOPROJECTSESSION Bind a SaveFolder to a project session, all or nothing.
%
%   [context, store] = bindFolderToProjectSession(dataFolder, selection)
%
%   The back end of ProjectBindingManager's Bind button (audit item 11).
%   It creates the project and/or subject when asked, then adds a new
%   session for DATAFOLDER or binds DATAFOLDER to an existing session. If
%   any step fails, what this call created is removed again (a new
%   subject, then a new project), so a retry never finds a half-created
%   project or subject and never creates a duplicate project.
%
%   SELECTION fields:
%       project  projectUUID (char/string) of an existing project, or a
%                struct for UMITProjectStore.create (a new project)
%       subject  subjectID (char/string) of an existing subject, or a
%                struct for addSubject (a new subject)
%       session  sessionID (char/string) of an existing session (the folder
%                is bound to it, or relocated onto it), or a struct for
%                addSession (a new session; processedDataFolder is set to
%                DATAFOLDER)
%       replaceOrphanBinding  optional logical (default false): DATAFOLDER
%                may carry a binding to a project that is not on this
%                machine; it is replaced (a rebind event in its ledger)
%
%   Outputs: CONTEXT from UMITProjectStore.resolveProjectBinding, and the
%   project STORE.
%
%   Errors: Umitoolbox:UMITProjectStore:bindFolderFailed, with the failing
%   step's error as cause; the compensation outcome is in the message.
%
%   See also describeProjectBinding, UMITProjectStore.

arguments
    dataFolder (1,1) string
    selection (1,1) struct
end
dataFolder = char(dataFolder);
replaceOrphan = isfield(selection, 'replaceOrphanBinding') && ...
    logical(selection.replaceOrphanBinding);
for name = {'project', 'subject', 'session'}
    if ~isfield(selection, name{1})
        error('Umitoolbox:UMITProjectStore:bindFolderFailed', ...
            'selection.%s is required.', name{1});
    end
end

store = [];
createdProject = false;
createdSubjectUUID = '';
step = 'project';
try
    if isstruct(selection.project)
        store = UMITProjectStore.create(selection.project);
        createdProject = true;
    else
        store = UMITProjectStore.open(char(string(selection.project)));
    end

    step = 'subject';
    if isstruct(selection.subject)
        createdSubjectUUID = store.addSubject(selection.subject);
        subjectID = store.getSubjectInfoByUUID(createdSubjectUUID).subjectID;
    else
        subjectID = store.getSubjectInfo(char(string(selection.subject))).subjectID;
    end

    step = 'session';
    if isstruct(selection.session)
        sessionInfo = selection.session;
        sessionInfo.processedDataFolder = dataFolder;
        store.addSession(subjectID, sessionInfo, 'ReplaceOrphanBinding', replaceOrphan);
    else
        store.bindProcessedDataFolder(subjectID, char(string(selection.session)), ...
            dataFolder, 'ReplaceOrphanBinding', replaceOrphan);
    end
catch ME
    outcome = iCompensate(store, createdProject, createdSubjectUUID);
    failure = MException('Umitoolbox:UMITProjectStore:bindFolderFailed', ...
        'Binding "%s" failed at the %s step: %s %s', dataFolder, step, ME.message, outcome);
    failure = addCause(failure, ME);
    throw(failure)
end
context = UMITProjectStore.resolveProjectBinding(dataFolder);
end

function outcome = iCompensate(store, createdProject, createdSubjectUUID)
% Undo what this call created, newest first. Never throws.
outcome = 'Nothing was created.';
if isempty(store)
    return
end
try
    if createdProject
        ProjectInfo = store.getProjectInfo();
        store.deleteProject(ProjectInfo.projectUUID);
        outcome = 'The new project was removed again.';
    elseif ~isempty(createdSubjectUUID)
        store.removeSubjectFromProject(createdSubjectUUID);
        outcome = 'The new subject was removed again.';
    end
catch cleanupError
    outcome = sprintf(['Removing what was created also failed (%s); remove it in ' ...
        'the Project Manager.'], cleanupError.message);
end
end
