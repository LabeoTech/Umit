classdef UMITProjectStore < handle
    %UMITPROJECTSTORE Manage a centralized UMIT project metadata folder.
    %
    %   store = UMITProjectStore.create(projectInfo)
    %   store = UMITProjectStore.open(projectUUID)
    %   projects = UMITProjectStore.listProjects()
    %
    %   UMITProjectStore is the only supported write interface for the
    %   centralized UMIT project folder. All projects are stored below
    %   getUmitFolder('projects'). It creates canonical folders, writes
    %   metadata atomically, imports managed resources, archives/restores
    %   resources, performs transactional ID renames, and validates the project.
    %
    %   The centralized folder stores metadata and image-reference files
    %   only. Imaging data remains external and may be referenced by
    %   session metadata. Rigs are independent of any project -- see
    %   RigManagement/UMITRigStore for rig metadata and rig-owned resources
    %   (camera coregistration, calibration files). Sessions keep no rig
    %   pointer: a folder's rig is recorded in its rig ledger (UMITRigLedger).
    %   A folder's image reference and alignment are recorded only in its
    %   DataParams.registration; the project keeps no transform resources.
    %
    %   Main operations:
    %       getProjectsRoot
    %       create
    %       open
    %       resolveProjectRoot
    %       listProjects
    %       getProjectInfoByUUID
    %       readProjectBinding
    %       isOrphanProjectBinding
    %       removeOrphanProjectBinding
    %       resolveProjectBinding
    %       addSubject
    %       addSession
    %       addImageReference
    %       archiveResource
    %       restoreResource
    %       setActiveImageReference
    %       updateSubjectMetadata
    %       updateSessionMetadata
    %       updateProjectMetadata
    %       getMetadataFieldDescriptions
    %       updateResourceMetadata
    %       renameSubjectID
    %       renameSessionID
    %       getLockInfo
    %       clearStaleLock
    %       listSubjectResources
    %       listSessionResources
    %       getResource
    %       resolveResourcePath
    %       getActiveImageReference
    %       getSubjectInfoByUUID
    %       getSessionInfoByUUID
    %       getSessionContextByUUID
    %       bindRawDataFolder
    %       bindProcessedDataFolder
    %       getRawDataFolderBinding
    %       getProcessedDataFolderBinding
    %       getSaveFolderBindingStatus
    %       retrySaveFolderAvailability
    %       relocateSaveFolder
    %       relocateRawDataFolder
    %       repairSaveFolderBinding
    %       removeSessionFromProject
    %       removeSubjectFromProject
    %       deleteProject
    %       unbindRawDataFolder
    %       findSessionByDataFolder
    %       validate
    %
    %   Design rules:
    %       - UUIDs are immutable identities.
    %       - Subject and session IDs are validated filesystem names.
    %       - Display names are editable without renaming physical folders.
    %       - Resource filenames are immutable after import.
    %       - Archived resources remain registered and restorable.
    %       - Every mutation uses a project lock and atomic metadata writes.
    %       - Validity is per session: a project opens read-only only for
    %         store-level errors. Errors of one subject or session block the
    %         mutations of that subject or session only, except narrowly
    %         scoped SaveFolder repair/relocation/session-removal operations
    %         when the selected binding is the only inconsistency. Mutations
    %         validate what they touch (quick mode), not the whole project.
    %       - UMITProjectStore is the only supported creator, reader,
    %         validator and writer of folder project ledgers (UMITProjectLedger.json).

    properties (SetAccess = private)
        ProjectRoot
        Schema
        IsReadOnly = false
        LastValidationReport
    end

    methods (Static)
        function projectsRoot = getProjectsRoot()
            %GETPROJECTSROOT Return the single UMIT projects root folder.
            %
            %   projectsRoot = UMITProjectStore.getProjectsRoot()
            %
            %   The root is always resolved through getUmitFolder. Individual
            %   project locations are not user-configurable.

            errID = 'Umitoolbox:UMITProjectStore:projectsRootFailed';

            if exist('getUmitFolder', 'file') ~= 2
                error(errID, ...
                    'Required function getUmitFolder.m was not found.');
            end

            try
                projectsRoot = getUmitFolder('projects');
            catch ME
                error(errID, ...
                    'Could not resolve getUmitFolder(''projects''): %s', ...
                    ME.message);
            end

            projectsRoot = UMITProjectStore.iAbsolutePath(projectsRoot);

            if isfile(projectsRoot)
                error(errID, ...
                    'The UMIT projects root resolves to a file: %s', ...
                    projectsRoot);
            end

            if ~isfolder(projectsRoot)
                [ok, message] = mkdir(projectsRoot);
                if ~ok
                    error(errID, ...
                        'Could not create the UMIT projects root: %s', ...
                        message);
                end
            end
        end

        function obj = create(projectInfo)
            %CREATE Create a new project under the static UMIT projects root.
            %
            %   store = UMITProjectStore.create(projectInfo)
            %
            %   Required projectInfo field:
            %       projectName
            %
            %   Optional projectInfo field:
            %       description
            %
            %   The project directory is generated internally under:
            %
            %       getUmitFolder('projects')
            %
            %   Users do not choose an individual project path.

            errID = 'Umitoolbox:UMITProjectStore:createFailed';

            if nargin < 1
                projectInfo = struct();
            end

            if ~isstruct(projectInfo) || ~isscalar(projectInfo)
                error(errID, '"projectInfo" must be a scalar struct.');
            end

            schema = getUMITProjectSchema();
            projectsRoot = UMITProjectStore.getProjectsRoot();

            projectName = UMITProjectStore.iGetTextField( ...
                projectInfo, 'projectName', '', false, errID);
            projectName = strtrim(projectName);
            if isempty(projectName)
                error(errID, 'Field "projectName" cannot be empty.');
            end

            description = UMITProjectStore.iGetTextField( ...
                projectInfo, 'description', '', true, errID);

            projectUUID = UMITProjectStore.iGenerateUUID();
            projectFolderName = ...
                UMITProjectStore.iMakeProjectFolderName( ...
                projectName, projectUUID);
            projectRoot = fullfile(projectsRoot, projectFolderName);

            if isfolder(projectRoot) || isfile(projectRoot)
                error(errID, ...
                    'Generated project destination already exists: %s', ...
                    projectRoot);
            end

            stagingRoot = fullfile(projectsRoot, ...
                ['.__creating_', projectUUID]);

            if isfolder(stagingRoot) || isfile(stagingRoot)
                error(errID, ...
                    'Could not allocate staging folder: %s', stagingRoot);
            end

            cleanupObj = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(stagingRoot));

            mkdir(stagingRoot);
            mkdir(fullfile(stagingRoot, schema.folders.internal));
            mkdir(fullfile(stagingRoot, schema.folders.internal, ...
                schema.folders.transactions));
            mkdir(fullfile(stagingRoot, schema.folders.internal, ...
                schema.folders.recovery));
            mkdir(fullfile(stagingRoot, schema.folders.internal, ...
                schema.folders.logs));
            mkdir(fullfile(stagingRoot, schema.folders.subjects));

            nowTime = datetime('now');
            ProjectInfo = struct();
            ProjectInfo.schemaVersion = schema.version;
            ProjectInfo.projectUUID = projectUUID;
            ProjectInfo.projectName = projectName;
            ProjectInfo.description = description;
            ProjectInfo = UMITProjectStore.iAddMetadataDefaults( ...
                ProjectInfo, 'project');
            ProjectInfo = UMITProjectStore.iApplyMetadataInput( ...
                ProjectInfo, projectInfo, schema.editableFields.project, errID);
            ProjectInfo.createdOn = nowTime;
            ProjectInfo.modifiedOn = nowTime;
            ProjectInfo.subjectRegistry = ...
                UMITProjectStore.iEmptyRegistry();

            projectFile = fullfile( ...
                stagingRoot, schema.files.projectMetadata);
            saveMatAtomic(projectFile, ...
                schema.metadataVariables.project, ProjectInfo);

            stagedStore = UMITProjectStore( ...
                stagingRoot, schema, false);
            report = stagedStore.validate('Mode', 'full');
            if ~report.isValid
                error(errID, ...
                    'New project failed validation: %s', ...
                    UMITProjectStore.iJoinIssueMessages( ...
                    report.errors));
            end

            [ok, message] = movefile( ...
                stagingRoot, projectRoot, 'f');
            if ~ok
                error(errID, ...
                    'Could not finalize project folder: %s', message);
            end

            clear cleanupObj

            obj = UMITProjectStore(projectRoot, schema, false);
            obj.LastValidationReport = ...
                obj.validate('Mode', 'quick');
            obj.iAppendLog( ...
                'createProject', projectUUID, 'completed');
        end

        function obj = open(projectUUID)
            %OPEN Open an existing project by immutable project UUID.
            %
            %   store = UMITProjectStore.open(projectUUID)
            %
            %   The project path is resolved only under:
            %
            %       getUmitFolder('projects')
            %
            %   A project with a store-level error (project metadata, schema,
            %   required folders, unregistered top-level folders, duplicate
            %   UUIDs) is returned in read-only mode. Errors of one subject or
            %   session do not make the project read-only: they block only
            %   the mutations of that subject or session (per-session
            %   validity). The validation result, with each issue's
            %   storeLevel flag, is available from store.LastValidationReport.

            errID = 'Umitoolbox:UMITProjectStore:openFailed';
            projectUUID = ...
                UMITProjectStore.iNormalizeProjectUUID(projectUUID);
            projectRoot = ...
                UMITProjectStore.resolveProjectRoot(projectUUID);

            defaultSchema = getUMITProjectSchema();
            projectFile = fullfile(projectRoot, ...
                defaultSchema.files.projectMetadata);

            if ~isfile(projectFile)
                error(errID, ...
                    'Project metadata file is missing: %s', projectFile);
            end

            loaded = load(projectFile, ...
                defaultSchema.metadataVariables.project, '-mat');
            if ~isfield(loaded, ...
                    defaultSchema.metadataVariables.project)
                error(errID, ...
                    'Project metadata variable "%s" is missing.', ...
                    defaultSchema.metadataVariables.project);
            end

            ProjectInfo = loaded.( ...
                defaultSchema.metadataVariables.project);
            if ~isstruct(ProjectInfo) || ~isscalar(ProjectInfo) || ...
                    ~isfield(ProjectInfo, 'schemaVersion')
                error(errID, ...
                    ['Project metadata does not define a valid ' ...
                     'schemaVersion.']);
            end

            if ~isfield(ProjectInfo, 'projectUUID') || ...
                    ~strcmpi(ProjectInfo.projectUUID, projectUUID)
                error(errID, ...
                    ['Resolved project metadata does not match requested ' ...
                     'project UUID: %s'], projectUUID);
            end

            schema = getUMITProjectSchema( ...
                ProjectInfo.schemaVersion);
            obj = UMITProjectStore(projectRoot, schema, false);
            report = obj.validate('Mode', 'quick');
            obj.LastValidationReport = report;
            obj.IsReadOnly = any([report.errors.storeLevel]);
        end

        function projectRoot = resolveProjectRoot(projectUUID)
            %RESOLVEPROJECTROOT Resolve one project UUID under the static root.
            %
            %   projectRoot = UMITProjectStore.resolveProjectRoot(projectUUID)

            errID = ...
                'Umitoolbox:UMITProjectStore:resolveProjectFailed';
            projectUUID = ...
                UMITProjectStore.iNormalizeProjectUUID(projectUUID);

            matchingRoots = UMITProjectStore.iFindProjectRoots(projectUUID);

            if isempty(matchingRoots)
                error(errID, ...
                    ['Project UUID was not found under ' ...
                     'getUmitFolder(''projects''): %s'], projectUUID);
            end

            if numel(matchingRoots) > 1
                error(errID, ...
                    ['Project UUID is duplicated under the static ' ...
                     'projects root: %s'], projectUUID);
            end

            projectRoot = char(matchingRoots(1));
        end

        function projects = listProjects()
            %LISTPROJECTS List projects found under the static UMIT root.
            %
            %   projects = UMITProjectStore.listProjects()
            %
            %   Returns one table row per folder containing project.mat.
            %   Malformed metadata remains visible with IsReadable=false.

            projectsRoot = UMITProjectStore.getProjectsRoot();
            schema = getUMITProjectSchema();

            ProjectUUID = strings(0, 1);
            ProjectName = strings(0, 1);
            Description = strings(0, 1);
            SchemaVersion = zeros(0, 1);
            CreatedOn = datetime.empty(0, 1);
            ModifiedOn = datetime.empty(0, 1);
            ProjectFolder = strings(0, 1);
            ProjectRoot = strings(0, 1);
            IsReadable = false(0, 1);
            Diagnostic = strings(0, 1);

            folderList = dir(projectsRoot);
            folderList = folderList([folderList.isdir]);

            for iFolder = 1:numel(folderList)
                folderName = folderList(iFolder).name;

                if ismember(folderName, {'.', '..'}) || ...
                        startsWith(folderName, '.') % internal staging folders
                    continue
                end

                projectRoot = fullfile(projectsRoot, folderName);
                projectFile = fullfile( ...
                    projectRoot, schema.files.projectMetadata);

                if ~isfile(projectFile)
                    continue
                end

                rowUUID = "";
                rowName = string(folderName);
                rowDescription = "";
                rowSchemaVersion = NaN;
                rowCreatedOn = NaT;
                rowModifiedOn = NaT;
                rowReadable = false;
                rowDiagnostic = "";

                try
                    loaded = load(projectFile, ...
                        schema.metadataVariables.project, '-mat');

                    if ~isfield(loaded, ...
                            schema.metadataVariables.project)
                        error('Missing metadata variable "%s".', ...
                            schema.metadataVariables.project);
                    end

                    ProjectInfo = loaded.( ...
                        schema.metadataVariables.project);

                    if ~isstruct(ProjectInfo) || ...
                            ~isscalar(ProjectInfo)
                        error('ProjectInfo is not a scalar structure.');
                    end

                    if isfield(ProjectInfo, 'projectUUID')
                        rowUUID = string(ProjectInfo.projectUUID);
                    end
                    if isfield(ProjectInfo, 'projectName')
                        rowName = string(ProjectInfo.projectName);
                    end
                    if isfield(ProjectInfo, 'description')
                        rowDescription = ...
                            string(ProjectInfo.description);
                    end
                    if isfield(ProjectInfo, 'schemaVersion')
                        rowSchemaVersion = ...
                            double(ProjectInfo.schemaVersion);
                    end
                    if isfield(ProjectInfo, 'createdOn') && ...
                            isdatetime(ProjectInfo.createdOn) && ...
                            isscalar(ProjectInfo.createdOn)
                        rowCreatedOn = ProjectInfo.createdOn;
                    end
                    if isfield(ProjectInfo, 'modifiedOn') && ...
                            isdatetime(ProjectInfo.modifiedOn) && ...
                            isscalar(ProjectInfo.modifiedOn)
                        rowModifiedOn = ProjectInfo.modifiedOn;
                    end

                    UMITProjectStore.iNormalizeProjectUUID(rowUUID);
                    rowReadable = true;

                catch ME
                    rowDiagnostic = string(ME.message);
                end

                ProjectUUID(end+1, 1) = rowUUID; %#ok<AGROW>
                ProjectName(end+1, 1) = rowName; %#ok<AGROW>
                Description(end+1, 1) = rowDescription; %#ok<AGROW>
                SchemaVersion(end+1, 1) = ...
                    rowSchemaVersion; %#ok<AGROW>
                CreatedOn(end+1, 1) = rowCreatedOn; %#ok<AGROW>
                ModifiedOn(end+1, 1) = rowModifiedOn; %#ok<AGROW>
                ProjectFolder(end+1, 1) = ...
                    string(folderName); %#ok<AGROW>
                ProjectRoot(end+1, 1) = ...
                    string(projectRoot); %#ok<AGROW>
                IsReadable(end+1, 1) = rowReadable; %#ok<AGROW>
                Diagnostic(end+1, 1) = rowDiagnostic; %#ok<AGROW>
            end

            projects = table( ...
                ProjectUUID, ...
                ProjectName, ...
                Description, ...
                SchemaVersion, ...
                CreatedOn, ...
                ModifiedOn, ...
                ProjectFolder, ...
                ProjectRoot, ...
                IsReadable, ...
                Diagnostic);

            if ~isempty(projects)
                projects = sortrows(projects, ...
                    {'IsReadable', 'ModifiedOn', 'ProjectName'}, ...
                    {'descend', 'descend', 'ascend'});
            end
        end

        function fileName = getProjectBindingFileName()
            %GETPROJECTBINDINGFILENAME Name of the folder's project ledger file.
            %
            %   The file exists in every folder that was ever bound (its
            %   history is kept after an unbind): use getProjectBindingState
            %   to know whether a folder is bound.

            schema = getUMITProjectSchema();
            fileName = schema.files.projectBinding;
        end

        function state = getProjectBindingState(dataFolder)
            %GETPROJECTBINDINGSTATE 'unbound', 'bound' or 'invalid'. Never writes.
            %
            %   state = UMITProjectStore.getProjectBindingState(dataFolder)
            %
            %   'unbound': no project ledger, or its binding was removed.
            %   'bound': the ledger binds the folder (resolveProjectBinding
            %   then decides orphan, mismatch or conflict against the store).
            %   'invalid': the ledger fails integrity validation.

            state = UMITProjectStore.iBindingState(dataFolder);
            if strcmp(state, 'none')
                state = 'unbound';
            end
        end

        function [ProjectInfo, projectRoot] = getProjectInfoByUUID(projectUUID)
            %GETPROJECTINFOBYUUID Resolve and return project metadata by UUID.
            %
            %   ProjectInfo = UMITProjectStore.getProjectInfoByUUID(projectUUID)
            %   [ProjectInfo, projectRoot] = ...

            store = UMITProjectStore.open(projectUUID);
            ProjectInfo = store.getProjectInfo();
            projectRoot = store.ProjectRoot;
        end

        function [ProjectBinding, bindingPath] = readProjectBinding(dataFolder)
            %READPROJECTBINDING Current binding of a folder from its project ledger.
            %
            %   ProjectBinding = UMITProjectStore.readProjectBinding(dataFolder)
            %   [ProjectBinding, ledgerPath] = ...
            %
            %   Reads UMITProjectLedger.json (UMITProjectLedger) and returns
            %   the binding record {version, bindingUUID, projectUUID,
            %   subjectUUID, sessionUUID, folderRole, createdOn,
            %   modifiedOn}. Only the identity fields are read from the
            %   ledger; its snapshots are never read. A folder without a
            %   ledger, or whose binding was removed, has no binding.
            %   Never writes.

            errID = 'Umitoolbox:UMITProjectStore:bindingReadFailed';

            if ~(ischar(dataFolder) || ...
                    (isstring(dataFolder) && isscalar(dataFolder)))
                error(errID, ...
                    '"dataFolder" must be a character vector or string scalar.');
            end

            dataFolder = UMITProjectStore.iAbsolutePath(dataFolder);
            if ~isfolder(dataFolder)
                error(errID, ...
                    'Data folder does not exist: %s', dataFolder);
            end

            schema = getUMITProjectSchema();
            bindingPath = UMITProjectLedger.filePath(dataFolder);
            try
                info = UMITProjectLedger.read(dataFolder);
            catch ME
                error(errID, 'Could not read the project ledger: %s', ME.message);
            end
            if ~strcmp(info.state, 'Bound')
                error(errID, ...
                    'The folder is not bound to a project (project ledger state %s): %s', ...
                    info.state, bindingPath);
            end

            ProjectBinding = UMITProjectStore.iBindingRecordFromLedger(info.binding, schema);
            UMITProjectStore.iValidateProjectBindingStruct( ...
                ProjectBinding, schema, errID);
        end

        function [tf, ProjectBinding] = isOrphanProjectBinding(dataFolder)
            %ISORPHANPROJECTBINDING Test whether a valid link references no project.
            %
            %   tf = UMITProjectStore.isOrphanProjectBinding(dataFolder)
            %   [tf, ProjectBinding] = ...
            %
            %   A binding is orphaned only when its payload is structurally valid and
            %   its project UUID is not present below getUmitFolder('projects').
            %   Duplicate project UUIDs and malformed binding files are errors rather
            %   than orphan states.

            ProjectBinding = ...
                UMITProjectStore.readProjectBinding(dataFolder);

            matchingRoots = UMITProjectStore.iFindProjectRoots( ...
                ProjectBinding.projectUUID);
            if numel(matchingRoots) > 1
                % A duplicated project UUID is an error, not an orphan.
                UMITProjectStore.resolveProjectRoot(ProjectBinding.projectUUID);
            end
            tf = isempty(matchingRoots);
        end

        function ProjectBinding = removeStaleProjectBinding(dataFolder)
            %REMOVESTALEPROJECTBINDING Clear a link no session matches any more.
            %
            %   ProjectBinding = UMITProjectStore.removeStaleProjectBinding(dataFolder)
            %
            %   A stale link (audit item 9) binds the folder to a project that
            %   IS on this machine, but none of its sessions binds this folder
            %   any more (the session was removed, relocated or rebound while
            %   the folder was offline). describeProjectBinding reports it as
            %   'stale'. Appends an unbind event with reason 'staleBinding';
            %   the ledger and its history stay. Any other state is refused
            %   (Umitoolbox:UMITProjectStore:bindingNotStale).

            info = describeProjectBinding(dataFolder);
            if ~strcmp(info.state, 'stale')
                error('Umitoolbox:UMITProjectStore:bindingNotStale', ...
                    ['The folder''s project binding is "%s", not stale; only a ' ...
                     'stale binding can be cleared this way: %s'], info.state, ...
                    char(string(dataFolder)));
            end
            ProjectBinding = info.binding;
            UMITProjectLedger.unbind(UMITProjectStore.iAbsolutePath(dataFolder), ...
                'Reason', 'staleBinding', ...
                'Origin', 'UMITProjectStore.removeStaleProjectBinding');
        end

        function ProjectBinding = removeOrphanProjectBinding(dataFolder)
            %REMOVEORPHANPROJECTBINDING Remove a binding to a project absent here.
            %
            %   Appends an orphanRemoval event to the folder's project
            %   ledger (the history and the file stay); the folder becomes
            %   unbound.
            %
            %   ProjectBinding = ...
            %       UMITProjectStore.removeOrphanProjectBinding(dataFolder)
            %
            %   This operation refuses to remove a binding when its project UUID is
            %   currently resolvable. Available-project bindings must be detached with
            %   the owning store's normal unbind method so SessionInfo remains
            %   reciprocal.

            errID = ...
                'Umitoolbox:UMITProjectStore:removeOrphanBindingFailed';

            [isOrphan, ProjectBinding] = ...
                UMITProjectStore.isOrphanProjectBinding(dataFolder);

            if ~isOrphan
                error('Umitoolbox:UMITProjectStore:bindingProjectAvailable', ...
                    ['The referenced project is available. Use the owning ' ...
                     'UMITProjectStore unbind method instead.']);
            end

            dataFolder = UMITProjectStore.iAbsolutePath(dataFolder);
            try
                UMITProjectLedger.recordOrphanRemoval(dataFolder, ...
                    'Origin', 'UMITProjectStore.removeOrphanProjectBinding');
            catch ME
                error(errID, ...
                    'Could not record the orphan binding removal: %s', ME.message);
            end
        end

        function [context, store] = resolveProjectBinding(dataFolder)
            %RESOLVEPROJECTBINDING Resolve a SaveFolder link to project context.
            %
            %   context = UMITProjectStore.resolveProjectBinding(dataFolder)
            %   [context, store] = ...
            %
            %   The binding UUIDs are authoritative. The current folder path must
            %   also match the path registered in SessionInfo. A mismatch indicates
            %   a moved or copied bound folder and is not repaired silently.

            errID = 'Umitoolbox:UMITProjectStore:bindingResolveFailed';
            [ProjectBinding, bindingPath] = ...
                UMITProjectStore.readProjectBinding(dataFolder);

            store = UMITProjectStore.open(ProjectBinding.projectUUID);
            ProjectInfo = store.getProjectInfo();
            [SubjectInfo, ~] = store.getSubjectInfoByUUID( ...
                ProjectBinding.subjectUUID);
            [SessionInfo, ParentSubjectInfo, ~] = ...
                store.getSessionInfoByUUID(ProjectBinding.sessionUUID);

            if ~strcmpi(ParentSubjectInfo.uuid, SubjectInfo.uuid)
                error(errID, ...
                    'Binding subject UUID does not own the bound session UUID.');
            end

            [pathField, bindingField] = store.iGetBindingRoleFields( ...
                ProjectBinding.folderRole);

            if ~strcmp(SessionInfo.(bindingField), ...
                    ProjectBinding.bindingUUID)
                error(errID, ...
                    ['SessionInfo.%s does not match the binding UUID in ' ...
                     'the SaveFolder.'], bindingField);
            end

            registeredFolder = SessionInfo.(pathField);
            if isempty(registeredFolder)
                error(errID, ...
                    'The session has no registered %s.', pathField);
            end

            currentFolder = UMITProjectStore.iAbsolutePath(dataFolder);
            if ~strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath(currentFolder), ...
                    UMITProjectStore.iNormalizeComparisonPath(registeredFolder))
                error('Umitoolbox:UMITProjectStore:bindingPathMismatch', ...
                    ['The project ledger UUIDs are valid, but the current folder path ' ...
                     'does not match SessionInfo.%s. The folder may have been ' ...
                     'moved or copied.'], pathField);
            end

            context = struct();
            context.projectUUID = ProjectInfo.projectUUID;
            context.projectName = ProjectInfo.projectName;
            context.projectRoot = store.ProjectRoot;
            context.subjectUUID = SubjectInfo.uuid;
            context.subjectID = SubjectInfo.subjectID;
            context.sessionUUID = SessionInfo.uuid;
            context.sessionID = SessionInfo.sessionID;
            context.folderRole = ProjectBinding.folderRole;
            context.matchedField = ProjectBinding.folderRole;
            context.bindingUUID = ProjectBinding.bindingUUID;
            context.bindingPath = bindingPath;
            context.dataFolder = currentFolder;
        end

    end
    methods
        function ProjectInfo = getProjectInfo(obj)
            %GETPROJECTINFO Return current project metadata.

            ProjectInfo = obj.iLoadMetadata( ...
                fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), ...
                obj.Schema.metadataVariables.project);
        end

        function descriptions = getMetadataFieldDescriptions(obj, entityType)
            %GETMETADATAFIELDDESCRIPTIONS Return UI-agnostic field descriptions.
            %
            %   descriptions = store.getMetadataFieldDescriptions(entityType)
            %   returns a scalar structure whose field names are metadata
            %   fields and whose values are concise descriptions. entityType is
            %   'project', 'subject', or 'session'. Unknown entity types return
            %   an empty structure so callers can safely omit a tooltip.

            entityType = lower(strtrim(char(string(entityType))));
            descriptions = struct();
            if isfield(obj.Schema.metadataFieldDescriptions, entityType)
                descriptions = obj.Schema.metadataFieldDescriptions.(entityType);
            end
        end

        function SubjectInfo = getSubjectInfo(obj, subjectID)
            %GETSUBJECTINFO Return metadata for one subject ID.

            [~, SubjectInfo] = obj.iResolveSubject(subjectID);
        end

        function SessionInfo = getSessionInfo(obj, subjectID, sessionID)
            %GETSESSIONINFO Return metadata for one session ID.

            [~, ~, ~, SessionInfo] = ...
                obj.iResolveSession(subjectID, sessionID);
        end


        function [SubjectInfo, subjectRecord] = getSubjectInfoByUUID(obj, subjectUUID)
            %GETSUBJECTINFOBYUUID Return subject metadata by immutable UUID.

            subjectUUID = UMITProjectStore.iNormalizeUUIDInput(subjectUUID);
            ProjectInfo = obj.getProjectInfo();

            if isempty(ProjectInfo.subjectRegistry)
                error('Umitoolbox:UMITProjectStore:subjectNotFound', ...
                    'Subject UUID was not found: %s', subjectUUID);
            end

            idx = find(strcmpi( ...
                {ProjectInfo.subjectRegistry.uuid}, subjectUUID));

            if isempty(idx)
                error('Umitoolbox:UMITProjectStore:subjectNotFound', ...
                    'Subject UUID was not found: %s', subjectUUID);
            end
            if numel(idx) > 1
                error('Umitoolbox:UMITProjectStore:duplicateSubjectUUID', ...
                    'Subject UUID is registered more than once: %s', ...
                    subjectUUID);
            end

            subjectRecord = ProjectInfo.subjectRegistry(idx);
            subjectPath = obj.iResolveRelativePath( ...
                subjectRecord.relativePath);
            SubjectInfo = obj.iLoadMetadata( ...
                fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                obj.Schema.metadataVariables.subject);
        end

        function [SessionInfo, SubjectInfo, sessionRecord] = ...
                getSessionInfoByUUID(obj, sessionUUID)
            %GETSESSIONINFOBYUUID Return session metadata by immutable UUID.

            [~, SubjectInfo, ~, SessionInfo, ~, ~, sessionRecord] = ...
                obj.iResolveSessionByUUID(sessionUUID);
        end

        function context = getSessionContextByUUID(obj, sessionUUID)
            %GETSESSIONCONTEXTBYUUID Return project/subject/session context.

            ProjectInfo = obj.getProjectInfo();
            [SessionInfo, SubjectInfo] = ...
                obj.getSessionInfoByUUID(sessionUUID);

            context = struct();
            context.projectUUID = ProjectInfo.projectUUID;
            context.projectName = ProjectInfo.projectName;
            context.projectRoot = obj.ProjectRoot;
            context.subjectUUID = SubjectInfo.uuid;
            context.subjectID = SubjectInfo.subjectID;
            context.sessionUUID = SessionInfo.uuid;
            context.sessionID = SessionInfo.sessionID;
            context.rawDataFolder = SessionInfo.rawDataFolder;
            context.processedDataFolder = ...
                SessionInfo.processedDataFolder;
        end

        function resources = listSubjectResources(obj, subjectID, varargin)
            %LISTSUBJECTRESOURCES List resources owned by one subject.
            %
            %   resources = store.listSubjectResources(subjectID)
            %   resources = store.listSubjectResources(subjectID, ...
            %       'Type', resourceType, ...
            %       'Status', statusFilter, ...
            %       'VerifyFiles', true)
            %
            %   By default, archived resources are excluded. Each returned
            %   record includes ownerType, ownerUUID, absolutePath, and
            %   fileExists in addition to the stored registry fields.

            [~, SubjectInfo] = obj.iResolveSubject(subjectID);
            resources = obj.iListResources( ...
                SubjectInfo, 'subject', SubjectInfo.uuid, varargin{:});
        end

        function resources = listSessionResources(obj, subjectID, sessionID, varargin)
            %LISTSESSIONRESOURCES List resources owned by one session.
            %
            %   resources = store.listSessionResources(subjectID, sessionID)
            %   resources = store.listSessionResources(subjectID, sessionID, ...
            %       'Type', resourceType, ...
            %       'Status', statusFilter, ...
            %       'VerifyFiles', true)

            [~, ~, ~, SessionInfo] = ...
                obj.iResolveSession(subjectID, sessionID);
            resources = obj.iListResources( ...
                SessionInfo, 'session', SessionInfo.uuid, varargin{:});
        end

        function resource = getResource(obj, resourceUUID)
            %GETRESOURCE Return one managed resource by immutable UUID.
            %
            %   resource = store.getResource(resourceUUID)
            %
            %   The returned record includes its owner type, owner UUID,
            %   absolute path, and current file-existence state.

            resourceUUID = UMITProjectStore.iNormalizeUUIDInput(resourceUUID);
            locations = obj.iFindResourceLocations(resourceUUID);

            if isempty(locations)
                error('Umitoolbox:UMITProjectStore:resourceNotFound', ...
                    'Resource UUID was not found: %s', resourceUUID);
            end

            if numel(locations) > 1
                error('Umitoolbox:UMITProjectStore:duplicateResourceUUID', ...
                    ['Resource UUID is registered more than once and cannot ' ...
                    'be resolved safely: %s'], resourceUUID);
            end

            location = locations(1);
            record = location.metadata.resourceRegistry(location.recordIndex);
            resource = obj.iEnrichResourceRecord( ...
                record, location.ownerType, location.metadata.uuid);
        end

        function filePath = resolveResourcePath(obj, resourceUUID, varargin)
            %RESOLVERESOURCEPATH Resolve one managed resource to an absolute path.
            %
            %   filePath = store.resolveResourcePath(resourceUUID)
            %   filePath = store.resolveResourcePath(resourceUUID, ...
            %       'RequireFile', false)
            %
            %   RequireFile is true by default.

            p = inputParser;
            p.FunctionName = 'UMITProjectStore.resolveResourcePath';
            addParameter(p, 'RequireFile', true, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            resource = obj.getResource(resourceUUID);
            if p.Results.RequireFile && ~resource.fileExists
                error('Umitoolbox:UMITProjectStore:missingResourceFile', ...
                    'Registered resource file is missing: %s', ...
                    resource.absolutePath);
            end

            filePath = resource.absolutePath;
        end

        function resource = getActiveImageReference(obj, subjectID)
            %GETACTIVEIMAGEREFERENCE Return the active subject image reference.
            %
            %   Returns [] when the subject has no active image reference.

            owner = obj.iResolveOwner('subject', {subjectID});
            resource = obj.iGetActiveResource(owner, 'imageReference');
        end

        function context = findSessionByDataFolder(obj, dataFolder)
            %FINDSESSIONBYDATAFOLDER Resolve context from a folder's project ledger.
            %
            %   context = store.findSessionByDataFolder(dataFolder)
            %
            %   Returns [] when the folder is not bound (no project ledger, or its binding was removed).
            %   The binding must resolve to this store's project UUID.

            if ~(ischar(dataFolder) || ...
                    (isstring(dataFolder) && isscalar(dataFolder)))
                error('Umitoolbox:UMITProjectStore:invalidPath', ...
                    '"dataFolder" must be a character vector or string scalar.');
            end

            dataFolder = UMITProjectStore.iAbsolutePath(dataFolder);
            obj.iAssertExternalDataPath(dataFolder, 'dataFolder');

            if strcmp(UMITProjectStore.iBindingState(dataFolder), 'none')
                context = [];
                return
            end

            ProjectBinding = ...
                UMITProjectStore.readProjectBinding(dataFolder);
            ProjectInfo = obj.getProjectInfo();

            if ~strcmpi(ProjectBinding.projectUUID, ...
                    ProjectInfo.projectUUID)
                error('Umitoolbox:UMITProjectStore:bindingProjectMismatch', ...
                    ['The folder is bound to project UUID %s, not the ' ...
                     'currently open project UUID %s.'], ...
                    ProjectBinding.projectUUID, ProjectInfo.projectUUID);
            end

            [SessionInfo, SubjectInfo] = ...
                obj.getSessionInfoByUUID(ProjectBinding.sessionUUID);

            if ~strcmpi(SessionInfo.subjectUUID, ...
                    ProjectBinding.subjectUUID) || ...
                    ~strcmpi(SubjectInfo.uuid, ...
                    ProjectBinding.subjectUUID)
                error('Umitoolbox:UMITProjectStore:bindingOwnerMismatch', ...
                    'Binding subject UUID does not own the bound session.');
            end

            [pathField, bindingField] = obj.iGetBindingRoleFields( ...
                ProjectBinding.folderRole);
            obj.iAssertBindingMatchesSession( ...
                ProjectBinding, SessionInfo, SubjectInfo, pathField, ...
                bindingField, dataFolder);

            context = struct();
            context.projectUUID = ProjectInfo.projectUUID;
            context.projectName = ProjectInfo.projectName;
            context.projectRoot = obj.ProjectRoot;
            context.subjectID = SubjectInfo.subjectID;
            context.subjectUUID = SubjectInfo.uuid;
            context.sessionID = SessionInfo.sessionID;
            context.sessionUUID = SessionInfo.uuid;
            context.matchedField = pathField;
            context.folderRole = pathField;
            context.bindingUUID = ProjectBinding.bindingUUID;
            context.dataFolder = dataFolder;
        end

        function subjectUUID = addSubject(obj, subjectInfo)
            %ADDSUBJECT Add a subject and its canonical metadata folders.
            %
            %   subjectUUID = store.addSubject(subjectInfo)
            %
            %   Required subjectInfo field:
            %       subjectID
            %
            %   Optional fields:
            %       displayName
            %       description

            errID = 'Umitoolbox:UMITProjectStore:addSubjectFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('addSubject');
            obj.iAssertHealthyForMutation(UMITProjectStore.iScope());

            if ~isstruct(subjectInfo) || ~isscalar(subjectInfo)
                error(errID, '"subjectInfo" must be a scalar struct.');
            end
            UMITProjectStore.iAssertKnownInputFields(subjectInfo, ...
                [{'subjectID'}, obj.Schema.editableFields.subject], errID);

            subjectID = obj.iNormalizeManagedID( ...
                UMITProjectStore.iGetTextField(subjectInfo, 'subjectID', '', false, errID), ...
                'subjectID', false);
            displayName = UMITProjectStore.iGetTextField( ...
                subjectInfo, 'displayName', subjectID, true, errID);
            if isempty(strtrim(displayName))
                displayName = subjectID;
            end
            description = UMITProjectStore.iGetTextField( ...
                subjectInfo, 'description', '', true, errID);

            ProjectInfo = obj.getProjectInfo();
            originalProjectInfo = ProjectInfo;
            obj.iAssertUniqueRegistryID( ...
                ProjectInfo.subjectRegistry, subjectID, 'subject');

            subjectUUID = UMITProjectStore.iGenerateUUID();
            subjectRel = UMITProjectStore.iJoinRelative( ...
                obj.Schema.folders.subjects, subjectID);
            subjectPath = obj.iResolveRelativePath(subjectRel);

            if isfolder(subjectPath) || isfile(subjectPath)
                error(errID, ...
                    'Subject destination already exists: %s', subjectPath);
            end

            transactionPath = obj.iCreateTransactionFolder('addSubject');
            stagedSubjectPath = fullfile(transactionPath, subjectID);
            cleanupStage = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(transactionPath));

            obj.iCreateSubjectFolders(stagedSubjectPath);

            nowTime = datetime('now');
            SubjectInfo = struct();
            SubjectInfo.schemaVersion = obj.Schema.version;
            SubjectInfo.uuid = subjectUUID;
            SubjectInfo.subjectID = subjectID;
            SubjectInfo.displayName = displayName;
            SubjectInfo.description = description;
            SubjectInfo = UMITProjectStore.iAddMetadataDefaults( ...
                SubjectInfo, 'subject');
            SubjectInfo = UMITProjectStore.iApplyMetadataInput( ...
                SubjectInfo, subjectInfo, obj.Schema.editableFields.subject, errID);
            if isempty(strtrim(SubjectInfo.displayName))
                SubjectInfo.displayName = subjectID;
            end
            displayName = SubjectInfo.displayName;
            SubjectInfo.createdOn = nowTime;
            SubjectInfo.modifiedOn = nowTime;
            SubjectInfo.activeImageReferenceUUID = '';
            SubjectInfo.sessionRegistry = UMITProjectStore.iEmptyRegistry();
            SubjectInfo.resourceRegistry = UMITProjectStore.iEmptyResourceRegistry();

            saveMatAtomic( ...
                fullfile(stagedSubjectPath, obj.Schema.files.subjectMetadata), ...
                obj.Schema.metadataVariables.subject, SubjectInfo);

            [ok, message] = movefile(stagedSubjectPath, subjectPath, 'f');
            if ~ok
                error(errID, ...
                    'Could not install subject folder: %s', message);
            end

            record = UMITProjectStore.iNewRegistryRecord( ...
                subjectUUID, subjectID, displayName, subjectRel);
            ProjectInfo.subjectRegistry(end+1) = record;
            ProjectInfo.modifiedOn = datetime('now');

            try
                obj.iSaveProjectInfo(ProjectInfo);
                obj.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {subjectUUID}));
            catch ME
                UMITProjectStore.iRemoveFolderIfPresent(subjectPath);
                obj.iSaveProjectInfo(originalProjectInfo);
                rethrow(ME)
            end

            obj.iAppendLog('addSubject', subjectUUID, 'completed');
            clear cleanupStage lockCleanup
        end

        function sessionUUID = addSession(obj, subjectID, sessionInfo, varargin)
            %ADDSESSION Add a dataset-backed session under an existing subject.
            %
            %   sessionUUID = store.addSession(subjectID, sessionInfo)
            %   sessionUUID = store.addSession(..., 'ReplaceOrphanBinding', true)
            %
            %   ReplaceOrphanBinding (default false): the SaveFolder may carry
            %   a binding to a project that is not on this machine (an
            %   orphan); it is replaced by a rebind event in the folder's
            %   ledger (audit item 12). Any other existing binding is refused.
            %
            %   Required sessionInfo fields:
            %       sessionID
            %       processedDataFolder  Existing SaveFolder for the dataset.
            %
            %   Optional fields:
            %       displayName
            %       description
            %       acquisitionDate
            %
            %   Any other field is rejected (unsupportedField). Sessions keep
            %   no rig: the folder's rig is recorded in its rig ledger
            %   (UMITRigLedger) and read from there.
            %
            %   The session registry entry, SessionInfo, and SaveFolder
            %   the project ledger bind event are installed as one rollback-safe
            %   operation. A persisted session is therefore never created
            %   without its authoritative SaveFolder reference.

            errID = 'Umitoolbox:UMITProjectStore:addSessionFailed';
            p = inputParser;
            p.FunctionName = 'UMITProjectStore.addSession';
            addParameter(p, 'ReplaceOrphanBinding', false, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('addSession');

            if ~isstruct(sessionInfo) || ~isscalar(sessionInfo)
                error(errID, '"sessionInfo" must be a scalar struct.');
            end

            forbiddenFolderFields = intersect(fieldnames(sessionInfo), ...
                {'rawDataFolder', 'rawDataBindingUUID', ...
                 'processedDataBindingUUID'});
            if ~isempty(forbiddenFolderFields)
                error('Umitoolbox:UMITProjectStore:folderBindingRequired', ...
                    ['Raw-folder and binding UUID fields cannot be assigned ' ...
                     'through addSession. The store owns binding identities.']);
            end

            UMITProjectStore.iAssertKnownInputFields(sessionInfo, ...
                [{'sessionID', 'processedDataFolder'}, ...
                obj.Schema.editableFields.session], errID);

            if ~isfield(sessionInfo, 'processedDataFolder') || ...
                    isempty(strtrim(char(string( ...
                    sessionInfo.processedDataFolder))))
                error('Umitoolbox:UMITProjectStore:saveFolderRequired', ...
                    ['A session must reference an existing SaveFolder. Set ' ...
                     'sessionInfo.processedDataFolder when adding the session.']);
            end

            saveFolder = UMITProjectStore.iAbsolutePath( ...
                sessionInfo.processedDataFolder);
            if ~isfolder(saveFolder)
                error('Umitoolbox:UMITProjectStore:saveFolderUnavailable', ...
                    'SaveFolder does not exist: %s', saveFolder);
            end
            obj.iAssertExternalDataPath(saveFolder, ...
                'processedDataFolder');
            obj.iAssertValidSaveFolderDataset(saveFolder);

            existingMatches = obj.iFindSessionPathMatches(saveFolder);
            if ~isempty(existingMatches)
                error('Umitoolbox:UMITProjectStore:dataFolderAlreadyBound', ...
                    ['The SaveFolder is already registered to session %s ' ...
                     'through %s.'], existingMatches(1).sessionID, ...
                    existingMatches(1).matchedField);
            end

            if ~strcmp(UMITProjectStore.iBindingState(saveFolder), 'none') && ...
                    ~(p.Results.ReplaceOrphanBinding && ...
                    strcmp(UMITProjectStore.iBindingState(saveFolder), 'bound') && ...
                    UMITProjectStore.isOrphanProjectBinding(saveFolder))
                error('Umitoolbox:UMITProjectStore:bindingFileExists', ...
                    ['The SaveFolder is already bound to a project (%s). Resolve ' ...
                     'that binding before adding a session (an orphan binding can be ' ...
                     'replaced with ''ReplaceOrphanBinding'', true).'], ...
                    obj.Schema.files.projectBinding);
            end

            [subjectRecord, SubjectInfo, subjectPath] = ...
                obj.iResolveSubject(subjectID);
            obj.iAssertHealthyForMutation(UMITProjectStore.iScope( ...
                'subjects', {SubjectInfo.uuid}));
            originalSubjectInfo = SubjectInfo;

            sessionID = obj.iNormalizeManagedID( ...
                UMITProjectStore.iGetTextField( ...
                sessionInfo, 'sessionID', '', false, errID), ...
                'sessionID', false);
            displayName = UMITProjectStore.iGetTextField( ...
                sessionInfo, 'displayName', sessionID, true, errID);
            if isempty(strtrim(displayName))
                displayName = sessionID;
            end
            description = UMITProjectStore.iGetTextField( ...
                sessionInfo, 'description', '', true, errID);
            acquisitionDate = UMITProjectStore.iGetDateField( ...
                sessionInfo, 'acquisitionDate', NaT, errID);

            obj.iAssertUniqueRegistryID( ...
                SubjectInfo.sessionRegistry, sessionID, 'session');

            sessionUUID = UMITProjectStore.iGenerateUUID();
            sessionRel = UMITProjectStore.iJoinRelative( ...
                subjectRecord.relativePath, ...
                obj.Schema.folders.sessions, ...
                sessionID);
            sessionPath = obj.iResolveRelativePath(sessionRel);

            if isfolder(sessionPath) || isfile(sessionPath)
                error(errID, ...
                    'Session destination already exists: %s', sessionPath);
            end

            transactionPath = obj.iCreateTransactionFolder('addSession');
            stagedSessionPath = fullfile(transactionPath, sessionID);
            cleanupStage = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(transactionPath));

            obj.iCreateSessionFolders(stagedSessionPath);

            nowTime = datetime('now');
            SessionInfo = struct();
            SessionInfo.schemaVersion = obj.Schema.version;
            SessionInfo.uuid = sessionUUID;
            SessionInfo.sessionID = sessionID;
            SessionInfo.subjectUUID = SubjectInfo.uuid;
            SessionInfo.subjectID = SubjectInfo.subjectID;
            SessionInfo.displayName = displayName;
            SessionInfo.description = description;
            SessionInfo.acquisitionDate = acquisitionDate;
            SessionInfo = UMITProjectStore.iAddMetadataDefaults( ...
                SessionInfo, 'session');
            SessionInfo = UMITProjectStore.iApplyMetadataInput( ...
                SessionInfo, sessionInfo, obj.Schema.editableFields.session, errID);
            if isempty(strtrim(SessionInfo.displayName))
                SessionInfo.displayName = sessionID;
            end
            displayName = SessionInfo.displayName;
            SessionInfo.rawDataFolder = '';
            SessionInfo.rawDataBindingUUID = '';
            SessionInfo.processedDataFolder = saveFolder;
            SessionInfo.processedDataBindingUUID = ...
                UMITProjectStore.iGenerateUUID();
            SessionInfo.createdOn = nowTime;
            SessionInfo.modifiedOn = nowTime;
            SessionInfo.resourceRegistry = ...
                UMITProjectStore.iEmptyResourceRegistry();

            ProjectInfo = obj.getProjectInfo();
            ProjectBinding = struct();
            ProjectBinding.version = obj.Schema.projectBinding.version;
            ProjectBinding.bindingUUID = ...
                SessionInfo.processedDataBindingUUID;
            ProjectBinding.projectUUID = ProjectInfo.projectUUID;
            ProjectBinding.subjectUUID = SubjectInfo.uuid;
            ProjectBinding.sessionUUID = sessionUUID;
            ProjectBinding.folderRole = 'processedDataFolder';
            ProjectBinding.createdOn = nowTime;
            ProjectBinding.modifiedOn = nowTime;
            UMITProjectStore.iValidateProjectBindingStruct( ...
                ProjectBinding, obj.Schema, errID);

            saveMatAtomic( ...
                fullfile(stagedSessionPath, ...
                    obj.Schema.files.sessionMetadata), ...
                obj.Schema.metadataVariables.session, SessionInfo);

            ledgerBackup = UMITProjectStore.iBackupLedger(saveFolder);

            [ok, message] = movefile(stagedSessionPath, sessionPath, 'f');
            if ~ok
                error(errID, ...
                    'Could not install session folder: %s', message);
            end

            sessionRecord = UMITProjectStore.iNewRegistryRecord( ...
                sessionUUID, sessionID, displayName, sessionRel);
            SubjectInfo.sessionRegistry(end+1) = sessionRecord;
            SubjectInfo.modifiedOn = datetime('now');

            try
                obj.iSaveSubjectInfo(subjectPath, SubjectInfo);

                obj.iInstallBinding(saveFolder, ProjectBinding, ...
                    SubjectInfo.subjectID, sessionID);

                obj.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {SubjectInfo.uuid}, 'sessions', {sessionUUID}));
            catch ME
                UMITProjectStore.iRestoreLedger(ledgerBackup);
                UMITProjectStore.iRemoveFolderIfPresent(sessionPath);
                obj.iSaveSubjectInfo(subjectPath, originalSubjectInfo);
                rethrow(ME)
            end

            obj.iAppendLog('addSession', sessionUUID, 'completed');
            clear cleanupStage lockCleanup
        end

        function ProjectBinding = bindRawDataFolder( ...
                obj, subjectID, sessionID, rawDataFolder, varargin)
            %BINDRAWDATAFOLDER Bind one raw-data folder to a session.
            %
            %   store.bindRawDataFolder(subjectID, sessionID, rawDataFolder)
            %   store.bindRawDataFolder(subjectID, sessionID, '')
            %
            %   The RawFolder may be the same folder as the session
            %   SaveFolder. It must contain a .bin, .tif, or .tiff file.
            %
            %   Single RawFolder authority (audit item 7): the SaveFolder's
            %   DataParams.folders.RawFolder, written by the importers. An
            %   empty RAWDATAFOLDER takes it from there; a different folder
            %   is refused (Umitoolbox:UMITProjectStore:rawFolderMismatch).
            %   SessionInfo.rawDataFolder only records where the store placed
            %   the RawFolder's ledger. A SaveFolder whose DataParams has no
            %   RawFolder (older imports) accepts any RawFolder.
            %
            %   Name-Value option:
            %       ReplaceOrphanBinding - Replace a structurally valid binding
            %                              whose project UUID is unavailable.

            if isempty(rawDataFolder) || strlength(string(rawDataFolder)) == 0
                rawDataFolder = obj.iDataParamsRawFolder( ...
                    obj.getSessionInfo(subjectID, sessionID));
                if isempty(rawDataFolder)
                    error('Umitoolbox:UMITProjectStore:rawFolderUnknown', ...
                        ['The SaveFolder''s DataParams records no RawFolder; give ' ...
                         'the RawFolder explicitly.']);
                end
            end
            ProjectBinding = obj.iBindSessionDataFolder( ...
                subjectID, sessionID, rawDataFolder, 'rawDataFolder', ...
                varargin{:});
        end

        function ProjectBinding = bindProcessedDataFolder( ...
                obj, subjectID, sessionID, processedDataFolder, varargin)
            %BINDPROCESSEDDATAFOLDER Bind or relocate a session SaveFolder.
            %
            %   Name-Value option:
            %       ReplaceOrphanBinding - Replace a structurally valid binding
            %                              whose project UUID is unavailable.

            SessionInfo = obj.getSessionInfo(subjectID, sessionID);
            if isempty(SessionInfo.processedDataFolder)
                % Compatibility for schema-v1 projects created by older
                % releases. New sessions can never enter this state.
                ProjectBinding = obj.iBindSessionDataFolder( ...
                    subjectID, sessionID, processedDataFolder, ...
                    'processedDataFolder', varargin{:});
            else
                ProjectBinding = obj.relocateSaveFolder( ...
                    subjectID, sessionID, processedDataFolder, varargin{:});
            end
        end

        function ProjectBinding = getRawDataFolderBinding( ...
                obj, subjectID, sessionID)
            %GETRAWDATAFOLDERBINDING Return the validated raw-folder binding.

            ProjectBinding = obj.iGetSessionDataFolderBinding( ...
                subjectID, sessionID, 'rawDataFolder');
        end

        function ProjectBinding = getProcessedDataFolderBinding( ...
                obj, subjectID, sessionID)
            %GETPROCESSEDDATAFOLDERBINDING Return the validated SaveFolder binding.

            ProjectBinding = obj.iGetSessionDataFolderBinding( ...
                subjectID, sessionID, 'processedDataFolder');
        end

        function ProjectBinding = unbindRawDataFolder( ...
                obj, subjectID, sessionID)
            %UNBINDRAWDATAFOLDER Remove a reciprocal raw-folder binding.

            ProjectBinding = obj.iUnbindSessionDataFolder( ...
                subjectID, sessionID, 'rawDataFolder');
        end

        function ProjectBinding = unbindProcessedDataFolder( ...
                ~, ~, ~)
            %UNBINDPROCESSEDDATAFOLDER Deprecated: sessions require a SaveFolder.
            %
            %   Use relocateSaveFolder to preserve the session, or
            %   removeSessionFromProject for explicit ledger removal.

            ProjectBinding = [];
            error('Umitoolbox:UMITProjectStore:unbindSaveFolderDeprecated', ...
                ['A persisted session must retain a SaveFolder. Use ' ...
                 'relocateSaveFolder or removeSessionFromProject instead.']);
        end

        function status = getSaveFolderBindingStatus( ...
                obj, subjectID, sessionID)
            %GETSAVEFOLDERBINDINGSTATUS Inspect runtime SaveFolder state.
            %
            %   The stored path is never changed by this read-only check.
            %   status.state is available, missing, invalid, or conflicting.

            [~, SubjectInfo, ~, SessionInfo] = ...
                obj.iResolveSession(subjectID, sessionID);
            status = obj.iGetBindingRuntimeStatus( ...
                SessionInfo, SubjectInfo, 'processedDataFolder');
        end

        function status = retrySaveFolderAvailability( ...
                obj, subjectID, sessionID)
            %RETRYSAVEFOLDERAVAILABILITY Re-run the read-only runtime check.

            status = obj.getSaveFolderBindingStatus( ...
                subjectID, sessionID);
        end

        function ProjectBinding = relocateSaveFolder( ...
                obj, subjectID, sessionID, newSaveFolder, varargin)
            %RELOCATESAVEFOLDER Rebind a session to a different SaveFolder.
            %
            %   Session UUID, ID, and metadata are preserved. The target link
            %   and ledger update are rollback-safe. An unavailable old folder
            %   does not prevent relocation.

            ProjectBinding = obj.iRelocateSessionDataFolder( ...
                subjectID, sessionID, newSaveFolder, ...
                'processedDataFolder', varargin{:});
        end

        function ProjectBinding = relocateRawDataFolder( ...
                obj, subjectID, sessionID, newRawDataFolder, varargin)
            %RELOCATERAWDATAFOLDER Rebind a session to one RawFolder.
            %
            %   Session identity and metadata are preserved. The target link
            %   and ledger update are rollback-safe, and an unavailable old
            %   RawFolder does not prevent relocation. The new RawFolder may
            %   equal the SaveFolder and must contain recognizable raw data.
            %   It must equal the SaveFolder's DataParams.folders.RawFolder
            %   when that is recorded (audit item 7): after moving raw data,
            %   update DataParams first.

            ProjectBinding = obj.iRelocateSessionDataFolder( ...
                subjectID, sessionID, newRawDataFolder, ...
                'rawDataFolder', varargin{:});
        end

        function ProjectBinding = repairSaveFolderBinding( ...
                obj, subjectID, sessionID)
            %REPAIRSAVEFOLDERBINDING Synchronize the stored SaveFolder link.

            [~, ~, ~, SessionInfo] = ...
                obj.iResolveSession(subjectID, sessionID);
            if isempty(SessionInfo.processedDataFolder)
                error('Umitoolbox:UMITProjectStore:saveFolderRequired', ...
                    'The session has no stored SaveFolder to repair.');
            end
            ProjectBinding = obj.iRelocateSessionDataFolder( ...
                subjectID, sessionID, ...
                SessionInfo.processedDataFolder, ...
                'processedDataFolder', 'RepairExisting', true);
        end

        function removedSession = removeSessionFromProject( ...
                obj, subjectID, sessionID)
            %REMOVESESSIONFROMPROJECT Remove a session from the project ledger.
            %
            %   The project-owned session metadata and registry entry are
            %   removed. External scientific data is never deleted. When the
            %   SaveFolder is accessible, its matching project-ledger binding is removed (unbind event);
            %   malformed links that still occupy the managed filename are
            %   invalidated without touching other dataset files.

            removedSession = obj.iRemoveSessionFromProject( ...
                subjectID, sessionID);
        end

        function result = removeSubjectFromProject(obj, subjectUUID)
            %REMOVESUBJECTFROMPROJECT Remove one subject and its sessions safely.
            %
            %   result = store.removeSubjectFromProject(subjectUUID)
            %
            %   Subject identity is resolved from its immutable UUID.  All
            %   registered sessions are removed from the project ledger and
            %   their accessible raw/processed-folder bindings are removed in
            %   one rollback-safe operation.  External scientific data is
            %   never deleted.  The operation rejects unavailable or invalid
            %   bindings before changing project metadata, rather than leaving
            %   an orphaned binding that could later reappear with its folder.

            errID = 'Umitoolbox:UMITProjectStore:removeSubjectFailed';
            subjectUUID = UMITProjectStore.iNormalizeUUIDInput(subjectUUID);
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('removeSubjectFromProject'); %#ok<NASGU>
            obj.iAssertHealthyForMutation(UMITProjectStore.iScope( ...
                'subjectTrees', {subjectUUID}));

            [SubjectInfo, subjectRecord] = obj.getSubjectInfoByUUID(subjectUUID);
            ProjectInfo = obj.getProjectInfo();
            subjectIndex = find(strcmpi( ...
                {ProjectInfo.subjectRegistry.uuid}, subjectUUID), 1, 'first');
            if isempty(subjectIndex)
                error(errID, 'Subject UUID was not found: %s', subjectUUID);
            end
            subjectPath = obj.iResolveRelativePath(subjectRecord.relativePath);

            recoveryPath = obj.iCreateRecoveryFolder('removeSubjectFromProject');
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));
            linkChanges = obj.iPreflightSubjectRemoval( ...
                SubjectInfo, subjectPath, recoveryPath);
            backupProjectPath = fullfile(recoveryPath, ...
                obj.Schema.files.projectMetadata);
            copyfile(fullfile(obj.ProjectRoot, ...
                obj.Schema.files.projectMetadata), backupProjectPath, 'f');
            stagedSubjectPath = fullfile(recoveryPath, 'subject');

            try
                [ok, message] = movefile(subjectPath, stagedSubjectPath, 'f');
                if ~ok
                    error(errID, ...
                        'Could not stage subject metadata removal: %s', message);
                end

                ProjectInfo.subjectRegistry(subjectIndex) = [];
                ProjectInfo.modifiedOn = datetime('now');
                obj.iSaveProjectInfo(ProjectInfo);
                obj.iDeletePreflightedBindings(linkChanges);
                obj.iAssertValidAfterMutation(UMITProjectStore.iScope());
            catch ME
                obj.iRestorePreflightedBindings(linkChanges);
                if ~isfolder(subjectPath) && isfolder(stagedSubjectPath)
                    movefile(stagedSubjectPath, subjectPath, 'f');
                end
                copyfile(backupProjectPath, fullfile(obj.ProjectRoot, ...
                    obj.Schema.files.projectMetadata), 'f');
                rethrow(ME)
            end

            result = struct('subjectUUID', SubjectInfo.uuid, ...
                'subjectID', SubjectInfo.subjectID, ...
                'sessionsRemoved', numel(SubjectInfo.sessionRegistry));
            obj.iAppendLog('removeSubjectFromProject', SubjectInfo.uuid, ...
                sprintf('%d session(s) removed', result.sessionsRemoved));
            clear cleanupRecovery lockCleanup
        end

        function result = deleteProject(obj, projectUUID)
            %DELETEPROJECT Remove a centralized project after unbinding sessions.
            %
            %   result = store.deleteProject(projectUUID)
            %
            %   The project UUID prevents a stale GUI selection from deleting
            %   a different open store.  All affected data-folder bindings are
            %   preflighted, then removed with rollback support before the
            %   project-owned metadata and managed resources are deleted.
            %   Deletion is deliberately available for a project opened
            %   read-only after validation failure, so corrupted managed
            %   metadata does not become undeletable (audit item 23): an
            %   unreadable subject or session, an offline data folder, or a
            %   folder whose ledger no longer binds it to this project is
            %   skipped instead of blocking the deletion. Skipped folders
            %   keep their ledger, which then reads as an orphan binding
            %   (removeOrphanProjectBinding clears it). RESULT.skipped lists
            %   them (fields folder/node and reason).
            %   External raw and processed imaging folders are never deleted.

            errID = 'Umitoolbox:UMITProjectStore:deleteProjectFailed';
            projectUUID = UMITProjectStore.iNormalizeUUIDInput(projectUUID);
            lockCleanup = obj.iAcquireWriteLock('deleteProject'); %#ok<NASGU>

            ProjectInfo = obj.getProjectInfo();
            if ~strcmpi(ProjectInfo.projectUUID, projectUUID)
                error(errID, ['The requested project UUID does not match ' ...
                    'this open project store.']);
            end
            sessionCount = obj.iCountProjectSessions(ProjectInfo);

            projectsRoot = UMITProjectStore.getProjectsRoot();
            recoveryPath = fullfile(projectsRoot, ...
                ['.deleteProjectRecovery_' UMITProjectStore.iGenerateUUID()]);
            stagedProjectPath = fullfile(projectsRoot, ...
                ['.deleteProjectStage_' UMITProjectStore.iGenerateUUID()]);
            mkdir(recoveryPath);
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));
            [linkChanges, skipped] = obj.iPreflightProjectRemoval( ...
                ProjectInfo, recoveryPath);

            try
                [ok, message] = movefile( ...
                    obj.ProjectRoot, stagedProjectPath, 'f');
                if ~ok
                    error(errID, 'Could not stage project deletion: %s', message);
                end

                obj.iDeletePreflightedBindings(linkChanges);
                UMITProjectStore.iRemoveFolderIfPresent(stagedProjectPath);
            catch ME
                obj.iRestorePreflightedBindings(linkChanges);
                if ~isfolder(obj.ProjectRoot) && isfolder(stagedProjectPath)
                    movefile(stagedProjectPath, obj.ProjectRoot, 'f');
                end
                rethrow(ME)
            end

            result = struct('projectUUID', ProjectInfo.projectUUID, ...
                'subjectsRemoved', numel(ProjectInfo.subjectRegistry), ...
                'sessionsRemoved', sessionCount, ...
                'skipped', skipped);
            clear cleanupRecovery lockCleanup
        end

        function resourceUUID = addImageReference(obj, subjectID, sourceFile, resourceInfo)
            %ADDIMAGEREFERENCE Import a managed subject image-reference file.

            if nargin < 4
                resourceInfo = struct();
            end
            resourceUUID = obj.iAddManagedResource( ...
                'subject', {subjectID}, 'imageReference', ...
                sourceFile, resourceInfo);
        end

        function archiveResource(obj, resourceUUID, varargin)
            %ARCHIVERESOURCE Archive a managed resource without deleting it.
            %
            %   store.archiveResource(resourceUUID)
            %   store.archiveResource(resourceUUID, ...
            %       'ReplacementUUID', replacementUUID)
            %
            %   An active resource requires a valid non-archived replacement.

            errID = 'Umitoolbox:UMITProjectStore:archiveFailed';

            p = inputParser;
            p.FunctionName = 'archiveResource';
            addRequired(p, 'resourceUUID', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));
            addParameter(p, 'ReplacementUUID', '', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));
            parse(p, resourceUUID, varargin{:});

            resourceUUID = char(string(p.Results.resourceUUID));
            replacementUUID = char(string(p.Results.ReplacementUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('archiveResource');
            location = obj.iLocateResource(resourceUUID);
            scope = UMITProjectStore.iLocationScope(location);
            obj.iAssertHealthyForMutation(scope);
            metadata = location.metadata;
            originalMetadata = metadata;
            record = metadata.resourceRegistry(location.recordIndex);

            if strcmp(record.status, 'archived')
                error(errID, 'Resource is already archived: %s', resourceUUID);
            end

            resourceDef = obj.Schema.resourceTypes.(record.type);
            pointerField = resourceDef.activePointerField;

            if strcmp(record.status, 'active')
                if isempty(replacementUUID)
                    error(errID, ...
                        ['Active resource "%s" cannot be archived without ' ...
                        'a replacement resource.'], resourceUUID);
                end

                replacementIndex = obj.iFindResourceIndex( ...
                    metadata.resourceRegistry, replacementUUID);
                if isempty(replacementIndex)
                    error(errID, ...
                        'Replacement resource is not owned by the same node.');
                end

                replacement = metadata.resourceRegistry(replacementIndex);
                if ~strcmp(replacement.type, record.type) || ...
                        strcmp(replacement.status, 'archived') || ...
                        strcmpi(replacement.uuid, record.uuid)
                    error(errID, ...
                        'Replacement must be a different non-archived resource of the same type.');
                end

                metadata.resourceRegistry(replacementIndex).status = 'active';
                metadata.resourceRegistry(replacementIndex).modifiedOn = datetime('now');
                metadata.(pointerField) = replacement.uuid;
            end

            oldPath = obj.iResolveRelativePath(record.relativePath);
            archiveRel = obj.iBuildResourceRelativePath( ...
                location.ownerBaseRelativePath, record.type, ...
                obj.Schema.folders.archive, record.fileName);
            archivePath = obj.iResolveRelativePath(archiveRel);

            if isfile(archivePath)
                error(errID, ...
                    'Archive destination already exists: %s', archivePath);
            end

            [ok, message] = movefile(oldPath, archivePath, 'f');
            if ~ok
                error(errID, 'Could not archive resource: %s', message);
            end

            metadata.resourceRegistry(location.recordIndex).status = 'archived';
            metadata.resourceRegistry(location.recordIndex).archivedOn = datetime('now');
            metadata.resourceRegistry(location.recordIndex).modifiedOn = datetime('now');
            metadata.resourceRegistry(location.recordIndex).relativePath = archiveRel;
            metadata.modifiedOn = datetime('now');

            try
                obj.iSaveLocatedMetadata(location, metadata);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                movefile(archivePath, oldPath, 'f');
                obj.iSaveLocatedMetadata(location, originalMetadata);
                rethrow(ME)
            end

            obj.iAppendLog('archiveResource', resourceUUID, 'completed');
            clear lockCleanup
        end

        function restoreResource(obj, resourceUUID)
            %RESTORERESOURCE Restore an archived resource as available.
            %
            %   Restoring does not automatically make the resource active.

            errID = 'Umitoolbox:UMITProjectStore:restoreFailed';
            resourceUUID = char(string(resourceUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('restoreResource');
            location = obj.iLocateResource(resourceUUID);
            scope = UMITProjectStore.iLocationScope(location);
            obj.iAssertHealthyForMutation(scope);
            metadata = location.metadata;
            originalMetadata = metadata;
            record = metadata.resourceRegistry(location.recordIndex);

            if ~strcmp(record.status, 'archived')
                error(errID, 'Resource is not archived: %s', resourceUUID);
            end

            oldPath = obj.iResolveRelativePath(record.relativePath);
            activeRel = obj.iBuildResourceRelativePath( ...
                location.ownerBaseRelativePath, record.type, ...
                obj.Schema.folders.active, record.fileName);
            activePath = obj.iResolveRelativePath(activeRel);

            if isfile(activePath)
                error(errID, ...
                    'Restore destination already exists: %s', activePath);
            end

            [ok, message] = movefile(oldPath, activePath, 'f');
            if ~ok
                error(errID, 'Could not restore resource: %s', message);
            end

            metadata.resourceRegistry(location.recordIndex).status = 'available';
            metadata.resourceRegistry(location.recordIndex).archivedOn = NaT;
            metadata.resourceRegistry(location.recordIndex).modifiedOn = datetime('now');
            metadata.resourceRegistry(location.recordIndex).relativePath = activeRel;
            metadata.modifiedOn = datetime('now');

            try
                obj.iSaveLocatedMetadata(location, metadata);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                movefile(activePath, oldPath, 'f');
                obj.iSaveLocatedMetadata(location, originalMetadata);
                rethrow(ME)
            end

            obj.iAppendLog('restoreResource', resourceUUID, 'completed');
            clear lockCleanup
        end

        function setActiveImageReference(obj, subjectID, resourceUUID)
            %SETACTIVEIMAGEREFERENCE Select the active subject image reference.

            obj.iSetActiveResource('subject', {subjectID}, ...
                'imageReference', resourceUUID);
        end

        function updateProjectMetadata(obj, updates)
            %UPDATEPROJECTMETADATA Update editable project metadata fields.

            errID = 'Umitoolbox:UMITProjectStore:updateProjectFailed';
            obj.iAssertUpdateStruct(updates, errID);
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('updateProjectMetadata');
            obj.iAssertHealthyForMutation(UMITProjectStore.iScope());

            ProjectInfo = obj.getProjectInfo();
            originalProjectInfo = ProjectInfo;
            ProjectInfo = obj.iApplyEditableUpdates( ...
                ProjectInfo, updates, obj.Schema.editableFields.project, errID);
            ProjectInfo.modifiedOn = datetime('now');

            backupPath = obj.iCreateRecoveryFolder('updateProjectMetadata');
            cleanupBackup = onCleanup(@() UMITProjectStore.iRemoveFolderIfPresent(backupPath));
            copyfile(fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), ...
                fullfile(backupPath, obj.Schema.files.projectMetadata), 'f');

            try
                obj.iSaveProjectInfo(ProjectInfo);
                obj.iAssertValidAfterMutation(UMITProjectStore.iScope());
            catch ME
                obj.iSaveProjectInfo(originalProjectInfo);
                rethrow(ME)
            end

            obj.iAppendLog('updateProjectMetadata', ProjectInfo.projectUUID, 'completed');
            clear cleanupBackup lockCleanup
        end

        function updateSubjectMetadata(obj, subjectID, updates)
            %UPDATESUBJECTMETADATA Update editable subject metadata fields.

            errID = 'Umitoolbox:UMITProjectStore:updateSubjectFailed';
            obj.iAssertUpdateStruct(updates, errID);
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('updateSubjectMetadata');

            [subjectRecord, SubjectInfo, subjectPath, projectIndex, ProjectInfo] = ...
                obj.iResolveSubject(subjectID);
            scope = UMITProjectStore.iScope('subjects', {SubjectInfo.uuid});
            obj.iAssertHealthyForMutation(scope);
            SubjectInfo = obj.iApplyEditableUpdates( ...
                SubjectInfo, updates, obj.Schema.editableFields.subject, errID);
            SubjectInfo.modifiedOn = datetime('now');

            ProjectInfo.subjectRegistry(projectIndex).displayName = ...
                SubjectInfo.displayName;
            ProjectInfo.modifiedOn = datetime('now');

            backupPath = obj.iCreateRecoveryFolder('updateSubjectMetadata');
            cleanupBackup = onCleanup(@() UMITProjectStore.iRemoveFolderIfPresent(backupPath));
            copyfile(fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                fullfile(backupPath, obj.Schema.files.subjectMetadata), 'f');
            copyfile(fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), ...
                fullfile(backupPath, obj.Schema.files.projectMetadata), 'f');

            try
                obj.iSaveSubjectInfo(subjectPath, SubjectInfo);
                obj.iSaveProjectInfo(ProjectInfo);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                copyfile(fullfile(backupPath, obj.Schema.files.subjectMetadata), ...
                    fullfile(subjectPath, obj.Schema.files.subjectMetadata), 'f');
                copyfile(fullfile(backupPath, obj.Schema.files.projectMetadata), ...
                    fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), 'f');
                rethrow(ME)
            end

            obj.iAppendLog('updateSubjectMetadata', subjectRecord.uuid, 'completed');
            clear cleanupBackup lockCleanup
        end

        function updateSessionMetadata(obj, subjectID, sessionID, updates)
            %UPDATESESSIONMETADATA Update editable non-binding session metadata.
            %
            %   Data-folder paths and binding UUIDs are managed only through
            %   bindRawDataFolder, bindProcessedDataFolder, and the matching
            %   unbind methods. Any field that is not editable is rejected.

            errID = 'Umitoolbox:UMITProjectStore:updateSessionFailed';
            obj.iAssertUpdateStruct(updates, errID);

            forbiddenFields = intersect(fieldnames(updates), ...
                {'rawDataFolder', 'processedDataFolder', ...
                 'rawDataBindingUUID', 'processedDataBindingUUID'});
            if ~isempty(forbiddenFields)
                error('Umitoolbox:UMITProjectStore:folderBindingRequired', ...
                    ['Session data folders and binding UUIDs cannot be changed ' ...
                     'through updateSessionMetadata. Use the binding methods.']);
            end

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('updateSessionMetadata');

            [~, SubjectInfo, subjectPath, SessionInfo, ...
                sessionPath, sessionIndex] = ...
                obj.iResolveSession(subjectID, sessionID);
            scope = UMITProjectStore.iScope('subjects', {SubjectInfo.uuid}, ...
                'sessions', {SessionInfo.uuid});
            obj.iAssertHealthyForMutation(scope);

            SessionInfo = obj.iApplyEditableUpdates( ...
                SessionInfo, updates, obj.Schema.editableFields.session, errID);

            SessionInfo.modifiedOn = datetime('now');
            SubjectInfo.sessionRegistry(sessionIndex).displayName = ...
                SessionInfo.displayName;
            SubjectInfo.modifiedOn = datetime('now');

            backupPath = ...
                obj.iCreateRecoveryFolder('updateSessionMetadata');
            cleanupBackup = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(backupPath));
            copyfile(fullfile(sessionPath, ...
                obj.Schema.files.sessionMetadata), ...
                fullfile(backupPath, ...
                obj.Schema.files.sessionMetadata), 'f');
            copyfile(fullfile(subjectPath, ...
                obj.Schema.files.subjectMetadata), ...
                fullfile(backupPath, ...
                obj.Schema.files.subjectMetadata), 'f');

            try
                obj.iSaveSessionInfo(sessionPath, SessionInfo);
                obj.iSaveSubjectInfo(subjectPath, SubjectInfo);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                copyfile(fullfile(backupPath, ...
                    obj.Schema.files.sessionMetadata), ...
                    fullfile(sessionPath, ...
                    obj.Schema.files.sessionMetadata), 'f');
                copyfile(fullfile(backupPath, ...
                    obj.Schema.files.subjectMetadata), ...
                    fullfile(subjectPath, ...
                    obj.Schema.files.subjectMetadata), 'f');
                rethrow(ME)
            end

            obj.iAppendLog('updateSessionMetadata', ...
                SessionInfo.uuid, 'completed');
            clear cleanupBackup lockCleanup
        end

        function updateResourceMetadata(obj, resourceUUID, updates)
            %UPDATERESOURCEMETADATA Update a resource display name or description.

            errID = 'Umitoolbox:UMITProjectStore:updateResourceFailed';
            obj.iAssertUpdateStruct(updates, errID);
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('updateResourceMetadata');
            location = obj.iLocateResource(resourceUUID);
            scope = UMITProjectStore.iLocationScope(location);
            obj.iAssertHealthyForMutation(scope);
            metadata = location.metadata;
            originalMetadata = metadata;
            record = metadata.resourceRegistry(location.recordIndex);
            record = obj.iApplyEditableUpdates( ...
                record, updates, obj.Schema.editableFields.resource, errID);
            record.modifiedOn = datetime('now');
            metadata.resourceRegistry(location.recordIndex) = record;
            metadata.modifiedOn = datetime('now');

            try
                obj.iSaveLocatedMetadata(location, metadata);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                obj.iSaveLocatedMetadata(location, originalMetadata);
                rethrow(ME)
            end
            obj.iAppendLog('updateResourceMetadata', resourceUUID, 'completed');
            clear lockCleanup
        end

        function renameSubjectID(obj, oldSubjectID, newSubjectID)
            %RENAMESUBJECTID Transactionally rename a subject ID and folder.

            errID = 'Umitoolbox:UMITProjectStore:renameSubjectFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('renameSubjectID');

            [~, SubjectInfo, oldPath, projectIndex, ProjectInfo] = ...
                obj.iResolveSubject(oldSubjectID);
            scope = UMITProjectStore.iScope('subjectTrees', {SubjectInfo.uuid});
            obj.iAssertHealthyForMutation(scope);
            oldSubjectID = SubjectInfo.subjectID;
            newSubjectID = obj.iNormalizeManagedID(newSubjectID, 'subjectID', false);
            if strcmp(oldSubjectID, newSubjectID)
                clear lockCleanup
                return
            end
            isCaseOnlyRename = strcmpi(oldSubjectID, newSubjectID);

            obj.iAssertUniqueRegistryID( ...
                ProjectInfo.subjectRegistry, newSubjectID, 'subject', ...
                SubjectInfo.uuid);

            oldRel = ProjectInfo.subjectRegistry(projectIndex).relativePath;
            newRel = UMITProjectStore.iJoinRelative( ...
                obj.Schema.folders.subjects, newSubjectID);
            newPath = obj.iResolveRelativePath(newRel);

            if (isfolder(newPath) || isfile(newPath)) && ...
                    ~isCaseOnlyRename
                error(errID, ...
                    'Destination subject folder already exists: %s', newPath);
            end

            recoveryPath = obj.iCreateRecoveryFolder('renameSubjectID');
            backupSubject = fullfile(recoveryPath, 'subject');
            backupProject = fullfile(recoveryPath, obj.Schema.files.projectMetadata);
            copyfile(oldPath, backupSubject, 'f');
            copyfile(fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), ...
                backupProject, 'f');

            temporaryPath = fullfile(fileparts(oldPath), ...
                ['.__rename_', SubjectInfo.uuid(1:8)]);

            try
                [ok, message] = movefile(oldPath, temporaryPath, 'f');
                if ~ok
                    error(errID, 'Could not stage subject rename: %s', message);
                end

                SubjectInfo.subjectID = newSubjectID;
                SubjectInfo.modifiedOn = datetime('now');
                SubjectInfo.resourceRegistry = obj.iReplaceResourcePrefix( ...
                    SubjectInfo.resourceRegistry, oldRel, newRel);
                SubjectInfo.sessionRegistry = obj.iReplaceRegistryPrefix( ...
                    SubjectInfo.sessionRegistry, oldRel, newRel);

                for iSession = 1:numel(SubjectInfo.sessionRegistry)
                    sessionRel = SubjectInfo.sessionRegistry(iSession).relativePath;
                    sessionPath = obj.iResolveRelativePathFromBase( ...
                        temporaryPath, obj.iRelativeTail(sessionRel, newRel));
                    SessionInfo = obj.iLoadMetadata( ...
                        fullfile(sessionPath, obj.Schema.files.sessionMetadata), ...
                        obj.Schema.metadataVariables.session);
                    SessionInfo.subjectID = newSubjectID;
                    SessionInfo.resourceRegistry = obj.iReplaceResourcePrefix( ...
                        SessionInfo.resourceRegistry, oldRel, newRel);
                    SessionInfo.modifiedOn = datetime('now');
                    obj.iSaveSessionInfo(sessionPath, SessionInfo);
                end

                obj.iSaveSubjectInfo(temporaryPath, SubjectInfo);

                [ok, message] = movefile(temporaryPath, newPath, 'f');
                if ~ok
                    error(errID, 'Could not finalize subject rename: %s', message);
                end

                ProjectInfo.subjectRegistry(projectIndex).id = newSubjectID;
                ProjectInfo.subjectRegistry(projectIndex).relativePath = newRel;
                ProjectInfo.modifiedOn = datetime('now');
                obj.iSaveProjectInfo(ProjectInfo);
                obj.iAssertValidAfterMutation(scope);

            catch ME
                UMITProjectStore.iRemoveFolderIfPresent(temporaryPath);
                UMITProjectStore.iRemoveFolderIfPresent(newPath);
                copyfile(backupSubject, oldPath, 'f');
                copyfile(backupProject, ...
                    fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), 'f');
                warning(errID, ...
                    'Subject rename rolled back. Recovery snapshot: %s', recoveryPath);
                rethrow(ME)
            end

            UMITProjectStore.iRemoveFolderIfPresent(recoveryPath);
            obj.iAppendLog('renameSubjectID', SubjectInfo.uuid, ...
                sprintf('%s -> %s', oldSubjectID, newSubjectID));
            clear lockCleanup
        end

        function renameSessionID(obj, subjectID, oldSessionID, newSessionID)
            %RENAMESESSIONID Transactionally rename a session ID and folder.

            errID = 'Umitoolbox:UMITProjectStore:renameSessionFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('renameSessionID');

            [subjectRecord, SubjectInfo, subjectPath, SessionInfo, oldPath, sessionIndex] = ...
                obj.iResolveSession(subjectID, oldSessionID);
            scope = UMITProjectStore.iScope('subjects', {SubjectInfo.uuid}, ...
                'sessions', {SessionInfo.uuid});
            obj.iAssertHealthyForMutation(scope);
            oldSessionID = SessionInfo.sessionID;
            newSessionID = obj.iNormalizeManagedID(newSessionID, 'sessionID', false);
            if strcmp(oldSessionID, newSessionID)
                clear lockCleanup
                return
            end
            isCaseOnlyRename = strcmpi(oldSessionID, newSessionID);

            obj.iAssertUniqueRegistryID( ...
                SubjectInfo.sessionRegistry, newSessionID, 'session', ...
                SessionInfo.uuid);

            oldRel = SubjectInfo.sessionRegistry(sessionIndex).relativePath;
            newRel = UMITProjectStore.iJoinRelative( ...
                subjectRecord.relativePath, ...
                obj.Schema.folders.sessions, newSessionID);
            newPath = obj.iResolveRelativePath(newRel);

            if (isfolder(newPath) || isfile(newPath)) && ...
                    ~isCaseOnlyRename
                error(errID, ...
                    'Destination session folder already exists: %s', newPath);
            end

            recoveryPath = obj.iCreateRecoveryFolder('renameSessionID');
            backupSession = fullfile(recoveryPath, 'session');
            backupSubject = fullfile(recoveryPath, obj.Schema.files.subjectMetadata);
            copyfile(oldPath, backupSession, 'f');
            copyfile(fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                backupSubject, 'f');

            temporaryPath = fullfile(fileparts(oldPath), ...
                ['.__rename_', SessionInfo.uuid(1:8)]);

            try
                [ok, message] = movefile(oldPath, temporaryPath, 'f');
                if ~ok
                    error(errID, 'Could not stage session rename: %s', message);
                end

                SessionInfo.sessionID = newSessionID;
                SessionInfo.resourceRegistry = obj.iReplaceResourcePrefix( ...
                    SessionInfo.resourceRegistry, oldRel, newRel);
                SessionInfo.modifiedOn = datetime('now');
                obj.iSaveSessionInfo(temporaryPath, SessionInfo);

                [ok, message] = movefile(temporaryPath, newPath, 'f');
                if ~ok
                    error(errID, 'Could not finalize session rename: %s', message);
                end

                SubjectInfo.sessionRegistry(sessionIndex).id = newSessionID;
                SubjectInfo.sessionRegistry(sessionIndex).relativePath = newRel;
                SubjectInfo.modifiedOn = datetime('now');
                obj.iSaveSubjectInfo(subjectPath, SubjectInfo);
                obj.iAssertValidAfterMutation(scope);

            catch ME
                UMITProjectStore.iRemoveFolderIfPresent(temporaryPath);
                UMITProjectStore.iRemoveFolderIfPresent(newPath);
                copyfile(backupSession, oldPath, 'f');
                copyfile(backupSubject, ...
                    fullfile(subjectPath, obj.Schema.files.subjectMetadata), 'f');
                warning(errID, ...
                    'Session rename rolled back. Recovery snapshot: %s', recoveryPath);
                rethrow(ME)
            end

            UMITProjectStore.iRemoveFolderIfPresent(recoveryPath);
            obj.iAppendLog('renameSessionID', SessionInfo.uuid, ...
                sprintf('%s -> %s', oldSessionID, newSessionID));
            clear lockCleanup
        end

        function renameEntity(obj, entityType, entityUUID, newName)
            %RENAMEENTITY Rename a subject or session identified by UUID.
            %
            %   store.renameEntity('subject', subjectUUID, newSubjectID)
            %   store.renameEntity('session', sessionUUID, newSessionID)
            %
            %   A thin UUID-based dispatcher over renameSubjectID/
            %   renameSessionID. Those methods are keyed by the current
            %   (mutable) ID string and, for sessions, require the caller to
            %   already know the owning subjectID. This method exists for
            %   callers that only hold an entity's immutable UUID (e.g. a
            %   GUI list selection) and resolves the current ID(s) first.
            %
            %   Preconditions:
            %     - entityType is 'subject' or 'session'.
            %     - entityUUID resolves to exactly one such entity in this
            %       project.
            %     - All preconditions of the delegated rename method also
            %       apply (project writable/healthy, newName is a valid
            %       managed ID, no ID collision at the destination).
            %
            %   Postconditions / failure modes:
            %     - Identical to renameSubjectID/renameSessionID. This
            %       method performs no mutation of its own; it only
            %       resolves entityUUID to the arguments those methods
            %       expect, so any failure there is reported unchanged.

            errID = 'Umitoolbox:UMITProjectStore:renameEntityFailed';
            entityType = lower(char(string(entityType)));

            switch entityType
                case 'subject'
                    SubjectInfo = obj.getSubjectInfoByUUID(entityUUID);
                    obj.renameSubjectID(SubjectInfo.subjectID, newName);
                case 'session'
                    [SessionInfo, SubjectInfo] = ...
                        obj.getSessionInfoByUUID(entityUUID);
                    obj.renameSessionID(SubjectInfo.subjectID, ...
                        SessionInfo.sessionID, newName);
                otherwise
                    error(errID, ...
                        ['Unsupported entityType: %s. Use ''subject'' ' ...
                         'or ''session''.'], entityType);
            end
        end

        function moveSessionToSubject(obj, sessionUUID, targetSubjectUUID)
            %MOVESESSIONTOSUBJECT Re-parent a session to another subject.
            %
            %   store.moveSessionToSubject(sessionUUID, targetSubjectUUID)
            %
            %   Transactionally relocates a session folder from its current
            %   subject to a different subject WITHIN THE SAME project. This
            %   is distinct from renameSessionID, which only changes the
            %   sessionID string in place; here the session's parent subject
            %   changes while sessionUUID, sessionID, and every resource
            %   UUID under the session are preserved.
            %
            %   Preconditions:
            %     - Project is writable and passes full validation.
            %     - sessionUUID resolves to exactly one session in this
            %       project.
            %     - targetSubjectUUID resolves to exactly one subject in
            %       this project.
            %     - The target subject does not already have a session
            %       whose sessionID collides (case-insensitively) with the
            %       session being moved (sessionIDs are only unique within
            %       one subject, so a same-ID session can already exist
            %       under the destination).
            %
            %   Postconditions:
            %     - The session folder is moved under the target subject's
            %       "sessions" folder; SessionInfo.subjectUUID/subjectID are
            %       updated to the new owner; the session's resourceRegistry
            %       relative paths are rewritten to the new location.
            %     - The source and target subjects' sessionRegistry entries
            %       are updated together under one project write lock.
            %     - If the session has a bound rawDataFolder and/or
            %       processedDataFolder, each bound folder's project-ledger binding
            %       is rewritten so its subjectUUID matches the new owner
            %       (projectUUID and sessionUUID are unchanged, since the
            %       project itself does not change).
            %     - Moving a session to the subject it already belongs to is
            %       a no-op (returns without error or mutation).
            %
            %   Failure modes (all leave the project byte-for-byte
            %   unchanged; a recovery snapshot is restored automatically):
            %     - sessionUUID or targetSubjectUUID not found in this
            %       project.
            %     - The target subject already has a session with the same
            %       sessionID.
            %     - A destination folder collision, a filesystem move
            %       failure, or a post-mutation validation failure at any
            %       step triggers a full rollback of the folder move, both
            %       subject metadata files, and any rewritten bindings.

            errID = 'Umitoolbox:UMITProjectStore:moveSessionFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('moveSessionToSubject');

            targetSubjectUUID = ...
                UMITProjectStore.iNormalizeUUIDInput(targetSubjectUUID);

            [~, OldSubjectInfo, oldSubjectPath, SessionInfo, ...
                oldSessionPath, oldSessionIndex, oldSessionRecord] = ...
                obj.iResolveSessionByUUID(sessionUUID);
            % The session and both subject nodes (registries change).
            scope = UMITProjectStore.iScope( ...
                'subjects', {OldSubjectInfo.uuid, targetSubjectUUID}, ...
                'sessions', {SessionInfo.uuid});
            obj.iAssertHealthyForMutation(scope);

            if strcmpi(OldSubjectInfo.uuid, targetSubjectUUID)
                clear lockCleanup
                return
            end

            ProjectInfo = obj.getProjectInfo();
            targetIndex = find(strcmpi( ...
                {ProjectInfo.subjectRegistry.uuid}, targetSubjectUUID));
            if isempty(targetIndex)
                error(errID, 'Target subject UUID was not found: %s', ...
                    targetSubjectUUID);
            end

            targetSubjectRecord = ProjectInfo.subjectRegistry(targetIndex);
            targetSubjectPath = obj.iResolveRelativePath( ...
                targetSubjectRecord.relativePath);
            TargetSubjectInfo = obj.iLoadMetadata( ...
                fullfile(targetSubjectPath, obj.Schema.files.subjectMetadata), ...
                obj.Schema.metadataVariables.subject);

            obj.iAssertUniqueRegistryID( ...
                TargetSubjectInfo.sessionRegistry, SessionInfo.sessionID, ...
                'session');

            oldRel = oldSessionRecord.relativePath;
            newRel = UMITProjectStore.iJoinRelative( ...
                targetSubjectRecord.relativePath, obj.Schema.folders.sessions, ...
                SessionInfo.sessionID);
            newSessionPath = obj.iResolveRelativePath(newRel);

            if isfolder(newSessionPath) || isfile(newSessionPath)
                error(errID, 'Destination session folder already exists: %s', ...
                    newSessionPath);
            end

            recoveryPath = obj.iCreateRecoveryFolder('moveSessionToSubject');
            backupSession = fullfile(recoveryPath, 'session');
            backupOldSubject = fullfile(recoveryPath, 'old_subject.mat');
            backupTargetSubject = fullfile(recoveryPath, 'target_subject.mat');
            copyfile(oldSessionPath, backupSession, 'f');
            copyfile(fullfile(oldSubjectPath, obj.Schema.files.subjectMetadata), ...
                backupOldSubject, 'f');
            copyfile(fullfile(targetSubjectPath, obj.Schema.files.subjectMetadata), ...
                backupTargetSubject, 'f');

            boundRoles = {'rawDataFolder', 'processedDataFolder'};
            boundRoles = boundRoles(~cellfun(@(role) ...
                isempty(SessionInfo.(role)), boundRoles));
            bindingBackups = struct('role', {}, 'path', {}, 'backup', {});
            for iRole = 1:numel(boundRoles)
                role = boundRoles{iRole};
                bindingPath = fullfile(SessionInfo.(role), ...
                    obj.Schema.files.projectBinding);
                backupPath = fullfile(recoveryPath, ...
                    sprintf('%s_%s', role, obj.Schema.files.projectBinding));
                copyfile(bindingPath, backupPath, 'f');
                bindingBackups(end+1) = struct( ...
                    'role', {role}, 'path', {bindingPath}, ...
                    'backup', {backupPath}); %#ok<AGROW>
            end

            try
                [ok, message] = movefile(oldSessionPath, newSessionPath, 'f');
                if ~ok
                    error(errID, 'Could not relocate session folder: %s', ...
                        message);
                end

                SessionInfo.subjectUUID = TargetSubjectInfo.uuid;
                SessionInfo.subjectID = TargetSubjectInfo.subjectID;
                SessionInfo.resourceRegistry = obj.iReplaceResourcePrefix( ...
                    SessionInfo.resourceRegistry, oldRel, newRel);
                SessionInfo.modifiedOn = datetime('now');
                obj.iSaveSessionInfo(newSessionPath, SessionInfo);

                for iRole = 1:numel(bindingBackups)
                    role = bindingBackups(iRole).role;
                    ProjectBinding = UMITProjectStore.readProjectBinding( ...
                        SessionInfo.(role));
                    ProjectBinding.subjectUUID = TargetSubjectInfo.uuid;
                    obj.iInstallBinding(SessionInfo.(role), ProjectBinding, ...
                        TargetSubjectInfo.subjectID, SessionInfo.sessionID);
                end

                OldSubjectInfo.sessionRegistry(oldSessionIndex) = [];
                OldSubjectInfo.modifiedOn = datetime('now');
                obj.iSaveSubjectInfo(oldSubjectPath, OldSubjectInfo);

                newSessionRecord = oldSessionRecord;
                newSessionRecord.relativePath = newRel;
                TargetSubjectInfo.sessionRegistry(end+1) = newSessionRecord;
                TargetSubjectInfo.modifiedOn = datetime('now');
                obj.iSaveSubjectInfo(targetSubjectPath, TargetSubjectInfo);

                obj.iAssertValidAfterMutation(scope);

            catch ME
                UMITProjectStore.iRemoveFolderIfPresent(newSessionPath);
                if ~isfolder(oldSessionPath)
                    copyfile(backupSession, oldSessionPath, 'f');
                end
                copyfile(backupOldSubject, ...
                    fullfile(oldSubjectPath, obj.Schema.files.subjectMetadata), 'f');
                copyfile(backupTargetSubject, ...
                    fullfile(targetSubjectPath, obj.Schema.files.subjectMetadata), 'f');
                for iRole = 1:numel(bindingBackups)
                    copyfile(bindingBackups(iRole).backup, ...
                        bindingBackups(iRole).path, 'f');
                end
                warning(errID, ...
                    'Session move rolled back. Recovery snapshot: %s', ...
                    recoveryPath);
                rethrow(ME)
            end

            UMITProjectStore.iRemoveFolderIfPresent(recoveryPath);
            obj.iAppendLog('moveSessionToSubject', SessionInfo.uuid, ...
                sprintf('%s -> %s', OldSubjectInfo.subjectID, ...
                TargetSubjectInfo.subjectID));
            clear lockCleanup
        end

        function rebindSessionToProject(obj, sessionUUID, targetProjectUUID)
            %REBINDSESSIONTOPROJECT Move a session to a different project.
            %
            %   store.rebindSessionToProject(sessionUUID, targetProjectUUID)
            %
            %   Moves a session, and its bound data-folder links, from this
            %   project to a different UMITProjectStore project. A subject
            %   with the same subjectID is created in the target project if
            %   one does not already exist there (a fresh subjectUUID is
            %   assigned in that case: subject identity is NOT preserved
            %   across projects, only sessionUUID is).
            %
            %   *** WEAKER GUARANTEE THAN moveSessionToSubject ***
            %   Each UMITProjectStore instance owns its own write lock and
            %   its own atomic, rollback-safe transaction (see
            %   iAcquireWriteLock); there is no cross-store lock and no
            %   two-phase commit spanning both projects. This method is
            %   composed of two independently-atomic phases:
            %     1) DETACH the session from the source project (fully
            %        rolled back on failure; the target project is never
            %        touched during this phase).
            %     2) Only once (1) has committed, ATTACH the session under
            %        the target project.
            %   This ordering is required, not just a preference: a bound
            %   rawDataFolder/processedDataFolder's project-ledger binding can only
            %   ever satisfy ONE project's full validation at a time (its
            %   projectUUID/subjectUUID must match whichever project's
            %   SessionInfo currently references that folder). Attaching to
            %   the target before detaching from the source would make both
            %   projects' SessionInfo reference the same external folder at
            %   once, and neither project's own iAssertHealthyForMutation
            %   full-validate would then pass. Detach-first avoids that.
            %   The cost is the opposite failure mode: if phase 1 (detach)
            %   fails, nothing has changed anywhere. If phase 2 (attach)
            %   fails AFTER phase 1 already committed, this method attempts
            %   to compensate by restoring the session back into the source
            %   project from an on-disk recovery snapshot. If that
            %   compensating restore ALSO fails, the session ends up
            %   referenced by NEITHER project's registry -- its folder and
            %   metadata still exist under the recovery snapshot reported
            %   in the error, but reattaching it requires manual action.
            %   There is no protection against another writer mutating
            %   either project during the gap between phase 1 and phase 2.
            %
            %   Preconditions:
            %     - targetProjectUUID resolves to an existing, writable
            %       project that passes full validation and is different
            %       from this project (rebinding to the same project is a
            %       no-op).
            %     - sessionUUID resolves to exactly one session in this
            %       project.
            %     - The (possibly newly-created) same-subjectID subject in
            %       the target project does not already have a session
            %       whose sessionID collides with the session being moved.
            %
            %   Postconditions (full success):
            %     - The session folder, sessionUUID, sessionID, and every
            %       resource UUID under it are preserved, now under a
            %       subject with the source's subjectID in the target
            %       project. SessionInfo.subjectUUID is updated to that
            %       target-project subject's UUID.
            %     - Any bound rawDataFolder/processedDataFolder project-ledger binding
            %       file is rewritten with the target project's projectUUID
            %       and the target subject's UUID; sessionUUID and
            %       bindingUUID are unchanged.
            %     - The source project no longer references the session.
            %
            %   Failure modes:
            %     - targetProjectUUID not found, read-only, or invalid:
            %       rejected before the source project is touched at all.
            %     - Phase-1 (detach) failure: fully rolled back; neither
            %       project is changed.
            %     - Phase-2 (attach) failure: target project is rolled back
            %       to not referencing the session; this method then
            %       attempts to restore the session into the source
            %       project. See the "referenced by neither project" risk
            %       above if that compensating restore also fails.

            errID = 'Umitoolbox:UMITProjectStore:rebindSessionFailed';
            sessionUUID = UMITProjectStore.iNormalizeUUIDInput(sessionUUID);
            targetProjectUUID = ...
                UMITProjectStore.iNormalizeUUIDInput(targetProjectUUID);

            ProjectInfo = obj.getProjectInfo();
            if strcmpi(ProjectInfo.projectUUID, targetProjectUUID)
                return
            end

            targetStore = UMITProjectStore.open(targetProjectUUID);
            targetStore.iAssertWritable();
            obj.iAssertWritable();

            % --- Ensure a same-subjectID subject exists in the target,
            %     before touching the source project at all ---
            [SessionInfo, SubjectInfo] = obj.getSessionInfoByUUID(sessionUUID); %#ok<ASGLU>


            TargetProjectInfo = targetStore.getProjectInfo();
            targetSubjectIdx = targetStore.iFindRegistryIndex( ...
                TargetProjectInfo.subjectRegistry, SubjectInfo.subjectID);
            if isempty(targetSubjectIdx)
                targetStore.addSubject(struct( ...
                    'subjectID', SubjectInfo.subjectID, ...
                    'displayName', SubjectInfo.displayName, ...
                    'description', SubjectInfo.description));
            end

            % --- Phase 1: detach from the source project (its own atomic
            %     transaction; target is not touched here) ---
            sourceLock = obj.iAcquireWriteLock( ...
                'rebindSessionToProject_detach');
            try
                [~, SubjectInfo, subjectPath, SessionInfo, sessionPath, ...
                    sessionIndex, sessionRecord] = ...
                    obj.iResolveSessionByUUID(sessionUUID);
                obj.iAssertHealthyForMutation(UMITProjectStore.iScope( ...
                    'subjects', {SubjectInfo.uuid}, 'sessions', {sessionUUID}));

                recoveryPath = obj.iCreateRecoveryFolder( ...
                    'rebindSessionToProject');
                stagingSessionPath = fullfile(recoveryPath, 'session');
                copyfile(sessionPath, stagingSessionPath, 'f');
                backupSubject = fullfile(recoveryPath, 'source_subject.mat');
                copyfile(fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                    backupSubject, 'f');

                UMITProjectStore.iRemoveFolderIfPresent(sessionPath);
                SubjectInfo.sessionRegistry(sessionIndex) = [];
                SubjectInfo.modifiedOn = datetime('now');
                obj.iSaveSubjectInfo(subjectPath, SubjectInfo);
                obj.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {SubjectInfo.uuid}));

            catch ME
                if exist('sessionPath', 'var') && exist('stagingSessionPath', 'var') ...
                        && isfile(fullfile(stagingSessionPath, ...
                        obj.Schema.files.sessionMetadata)) && ~isfolder(sessionPath)
                    copyfile(stagingSessionPath, sessionPath, 'f');
                end
                if exist('subjectPath', 'var') && exist('backupSubject', 'var') ...
                        && isfile(backupSubject)
                    copyfile(backupSubject, ...
                        fullfile(subjectPath, obj.Schema.files.subjectMetadata), 'f');
                end
                if exist('recoveryPath', 'var')
                    UMITProjectStore.iRemoveFolderIfPresent(recoveryPath);
                end
                clear sourceLock
                error(errID, ...
                    ['Could not detach session from the source project; ' ...
                     'no changes were made anywhere. Cause: %s'], ME.message);
            end
            clear sourceLock

            % --- Phase 2: attach under the target project, now that the
            %     source project no longer references the session ---
            newSessionPath = '';
            targetLock = targetStore.iAcquireWriteLock( ...
                'rebindSessionToProject_attach');
            try
                TargetProjectInfo = targetStore.getProjectInfo();
                targetSubjectIdx = targetStore.iFindRegistryIndex( ...
                    TargetProjectInfo.subjectRegistry, SubjectInfo.subjectID);
                targetSubjectRecord = ...
                    TargetProjectInfo.subjectRegistry(targetSubjectIdx);
                targetSubjectPath = targetStore.iResolveRelativePath( ...
                    targetSubjectRecord.relativePath);
                TargetSubjectInfo = targetStore.iLoadMetadata( ...
                    fullfile(targetSubjectPath, ...
                    targetStore.Schema.files.subjectMetadata), ...
                    targetStore.Schema.metadataVariables.subject);
                targetStore.iAssertHealthyForMutation(UMITProjectStore.iScope( ...
                    'subjects', {TargetSubjectInfo.uuid}));

                targetStore.iAssertUniqueRegistryID( ...
                    TargetSubjectInfo.sessionRegistry, ...
                    SessionInfo.sessionID, 'session');

                newRel = UMITProjectStore.iJoinRelative( ...
                    targetSubjectRecord.relativePath, ...
                    targetStore.Schema.folders.sessions, ...
                    SessionInfo.sessionID);
                newSessionPath = targetStore.iResolveRelativePath(newRel);

                if isfolder(newSessionPath) || isfile(newSessionPath)
                    error(errID, ...
                        'Destination session folder already exists: %s', ...
                        newSessionPath);
                end

                [ok, message] = movefile( ...
                    stagingSessionPath, newSessionPath, 'f');
                if ~ok
                    error(errID, ...
                        'Could not install session folder in target project: %s', ...
                        message);
                end

                AttachedSessionInfo = SessionInfo;
                AttachedSessionInfo.subjectUUID = TargetSubjectInfo.uuid;
                AttachedSessionInfo.subjectID = TargetSubjectInfo.subjectID;
                AttachedSessionInfo.resourceRegistry = ...
                    targetStore.iReplaceResourcePrefix( ...
                    AttachedSessionInfo.resourceRegistry, ...
                    sessionRecord.relativePath, newRel);
                AttachedSessionInfo.modifiedOn = datetime('now');
                targetStore.iSaveSessionInfo(newSessionPath, AttachedSessionInfo);

                boundRoles = {'rawDataFolder', 'processedDataFolder'};
                for iRole = 1:numel(boundRoles)
                    role = boundRoles{iRole};
                    if isempty(AttachedSessionInfo.(role))
                        continue
                    end
                    ProjectBinding = UMITProjectStore.readProjectBinding( ...
                        AttachedSessionInfo.(role));
                    ProjectBinding.projectUUID = TargetProjectInfo.projectUUID;
                    ProjectBinding.subjectUUID = TargetSubjectInfo.uuid;
                    targetStore.iInstallBinding(AttachedSessionInfo.(role), ...
                        ProjectBinding, TargetSubjectInfo.subjectID, ...
                        AttachedSessionInfo.sessionID);
                end

                newSessionRecord = sessionRecord;
                newSessionRecord.relativePath = newRel;
                TargetSubjectInfo.sessionRegistry(end+1) = newSessionRecord;
                TargetSubjectInfo.modifiedOn = datetime('now');
                targetStore.iSaveSubjectInfo( ...
                    targetSubjectPath, TargetSubjectInfo);

                targetStore.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {TargetSubjectInfo.uuid}, 'sessions', {sessionUUID}));

            catch ME
                if ~isempty(newSessionPath) && isfolder(newSessionPath)
                    % The staged copy was already moved here; move it back
                    % so the compensation step below has something to
                    % restore into the source project.
                    UMITProjectStore.iRemoveFolderIfPresent(stagingSessionPath);
                    movefile(newSessionPath, stagingSessionPath, 'f');
                end
                clear targetLock
                obj.iCompensateFailedRebindDetach( ...
                    sessionUUID, SubjectInfo.subjectID, SessionInfo, ...
                    stagingSessionPath, backupSubject, recoveryPath, ...
                    targetProjectUUID, ME);
                return
            end
            clear targetLock

            UMITProjectStore.iRemoveFolderIfPresent(recoveryPath);
            obj.iAppendLog('rebindSessionToProject', SessionInfo.uuid, ...
                sprintf('moved to project %s', targetProjectUUID));
            targetStore.iAppendLog('rebindSessionToProject', SessionInfo.uuid, ...
                sprintf('received from project %s', ProjectInfo.projectUUID));
        end

        function LockInfo = getLockInfo(obj)
            %GETLOCKINFO Return the current project-lock metadata.
            %
            %   Returns an empty struct when the project is not locked.

            lockFolder = fullfile(obj.ProjectRoot, ...
                obj.Schema.folders.internal, obj.Schema.folders.lock);
            lockMetadataPath = fullfile(lockFolder, ...
                obj.Schema.files.lockMetadata);
            LockInfo = struct();
            if ~isfolder(lockFolder)
                return
            end

            if ~isfile(lockMetadataPath)
                error('Umitoolbox:UMITProjectStore:invalidLockFile', ...
                    'Project lock metadata is missing: %s', lockMetadataPath);
            end

            loaded = load(lockMetadataPath, 'LockInfo', '-mat');
            if ~isfield(loaded, 'LockInfo') || ...
                    ~isstruct(loaded.LockInfo) || ~isscalar(loaded.LockInfo)
                error('Umitoolbox:UMITProjectStore:invalidLockFile', ...
                    'Project lock metadata is malformed: %s', lockMetadataPath);
            end
            LockInfo = loaded.LockInfo;
        end

        function clearStaleLock(obj, varargin)
            %CLEARSTALELOCK Archive and remove an abandoned project lock.
            %
            %   store.clearStaleLock()
            %   store.clearStaleLock('MinimumAgeMinutes', 60)
            %   store.clearStaleLock('Force', true)
            %
            %   Without Force=true, the lock must contain a valid creation time
            %   and be at least MinimumAgeMinutes old. The removed lock is moved
            %   into the project log folder for auditability.

            p = inputParser;
            p.FunctionName = 'clearStaleLock';
            addParameter(p, 'MinimumAgeMinutes', 60, ...
                @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
            addParameter(p, 'Force', false, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            lockFolder = fullfile(obj.ProjectRoot, ...
                obj.Schema.folders.internal, obj.Schema.folders.lock);
            if ~isfolder(lockFolder)
                return
            end

            force = p.Results.Force;
            if ~force
                LockInfo = obj.getLockInfo();
                if ~isfield(LockInfo, 'createdOn') || ...
                        ~(isdatetime(LockInfo.createdOn) && ...
                        isscalar(LockInfo.createdOn))
                    error('Umitoolbox:UMITProjectStore:invalidLockFile', ...
                        ['Lock age cannot be established. Use Force=true ' ...
                        'only after confirming that no writer is active.']);
                end

                ageMinutes = minutes(datetime('now') - LockInfo.createdOn);
                if ageMinutes < p.Results.MinimumAgeMinutes
                    error('Umitoolbox:UMITProjectStore:lockNotStale', ...
                        ['Project lock is only %.1f minutes old. The required ' ...
                        'minimum age is %.1f minutes.'], ...
                        ageMinutes, p.Results.MinimumAgeMinutes);
                end
            end

            timestamp = char(datetime('now', ...
                'Format', 'yyyyMMdd''T''HHmmss'));
            archiveName = sprintf('stale_lock_%s_%s', ...
                timestamp, UMITProjectStore.iGenerateUUID());
            archivePath = fullfile(obj.ProjectRoot, ...
                obj.Schema.folders.internal, obj.Schema.folders.logs, ...
                archiveName);

            [ok, message] = movefile(lockFolder, archivePath, 'f');
            if ~ok
                error('Umitoolbox:UMITProjectStore:lockRecoveryFailed', ...
                    'Could not archive stale lock: %s', message);
            end

            obj.iAppendLog('clearStaleLock', '', archiveName);
        end

        function report = validate(obj, varargin)
            %VALIDATE Validate project metadata, folders, registries, and files.
            %
            %   report = store.validate()
            %   report = store.validate('Mode', 'quick')
            %   report = store.validate('Mode', 'full')
            %   report = store.validate('Mode', 'quick', 'Scope', scope)
            %
            %   Quick mode performs all structural and registry checks but skips
            %   file checksum computation. Full mode also validates checksums.
            %
            %   Every issue carries a storeLevel flag: true for errors of the
            %   store itself (project metadata, schema, required folders,
            %   registries, folders below 'subjects' that belong to no
            %   registered subject, global UUIDs), false for errors that
            %   belong to one registered subject or session. Only store-level
            %   errors make a project read-only (per-session validity).
            %
            %   'Scope' (used by mutations, brief section 1.11) limits the
            %   traversal to the entities a mutation touches: a struct with
            %   cellstr fields subjects (subject nodes, without their
            %   sessions), sessions (session nodes; their subject node must
            %   be listed in subjects) and subjectTrees (subjects with all
            %   their sessions). Project-level checks always run; the global
            %   UUID scan is skipped. A scoped validation does not replace
            %   LastValidationReport.

            p = inputParser;
            p.FunctionName = 'UMITProjectStore.validate';
            addParameter(p, 'Mode', 'full', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));
            addParameter(p, 'Scope', [], @(x) isempty(x) || isstruct(x));
            parse(p, varargin{:});

            mode = lower(char(string(p.Results.Mode)));
            if ~ismember(mode, {'quick', 'full'})
                error('Umitoolbox:UMITProjectStore:invalidValidationMode', ...
                    'Validation mode must be "quick" or "full".');
            end
            scope = p.Results.Scope;
            isScoped = isstruct(scope);
            if isScoped
                scope = UMITProjectStore.iScope( ...
                    'subjects', scope.subjects, 'sessions', scope.sessions, ...
                    'subjectTrees', scope.subjectTrees);
            end

            report = obj.iNewValidationReport(mode);

            projectFile = fullfile(obj.ProjectRoot, ...
                obj.Schema.files.projectMetadata);
            if ~isfile(projectFile)
                report = obj.iAddIssue(report, 'error', ...
                    'missing_project_metadata', 'project', '', ...
                    'project.mat is missing.', false);
                report = obj.iFinishReport(report, [], isScoped);
                return
            end

            try
                ProjectInfo = obj.iLoadMetadata(projectFile, ...
                    obj.Schema.metadataVariables.project);
            catch ME
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_project_metadata', 'project', ...
                    obj.Schema.files.projectMetadata, ME.message, false);
                report = obj.iFinishReport(report, [], isScoped);
                return
            end

            report = obj.iValidateMetadataFields(report, ProjectInfo, ...
                obj.Schema.requiredProjectFields, 'project', ...
                obj.Schema.files.projectMetadata);

            if ~isfield(ProjectInfo, 'schemaVersion') || ...
                    ~isequal(ProjectInfo.schemaVersion, obj.Schema.version)
                report = obj.iAddIssue(report, 'error', ...
                    'schema_version_mismatch', 'project', ...
                    obj.Schema.files.projectMetadata, ...
                    'Project schemaVersion does not match the loaded schema.', ...
                    false);
            end

            % Required top-level and internal folders.
            for iFolder = 1:numel(obj.Schema.requiredTopFolders)
                rel = obj.Schema.requiredTopFolders{iFolder};
                if ~isfolder(obj.iResolveRelativePath(rel))
                    report = obj.iAddIssue(report, 'error', ...
                        'missing_required_folder', 'project', rel, ...
                        sprintf('Required project folder is missing: %s', rel), ...
                        true);
                end
            end

            for iFolder = 1:numel(obj.Schema.requiredInternalFolders)
                rel = UMITProjectStore.iJoinRelative( ...
                    obj.Schema.folders.internal, ...
                    obj.Schema.requiredInternalFolders{iFolder});
                if ~isfolder(obj.iResolveRelativePath(rel))
                    report = obj.iAddIssue(report, 'error', ...
                        'missing_internal_folder', 'project', rel, ...
                        sprintf('Required internal folder is missing: %s', rel), ...
                        true);
                end
            end

            report = obj.iValidateAllowedChildren(report, ...
                obj.ProjectRoot, '', 'project', { ...
                obj.Schema.files.projectMetadata, ...
                obj.Schema.folders.internal, ...
                obj.Schema.folders.subjects});

            if isfield(ProjectInfo, 'subjectRegistry')
                report = obj.iValidateRegistry(report, ...
                    ProjectInfo.subjectRegistry, 'subject', ...
                    obj.Schema.folders.subjects);
            end

            if isfield(ProjectInfo, 'subjectRegistry') && ...
                    obj.iIsRegistryStruct(ProjectInfo.subjectRegistry)
                for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                    subjectRecord = ProjectInfo.subjectRegistry(iSubject);
                    if isScoped && ~UMITProjectStore.iScopeHasSubject( ...
                            scope, subjectRecord.uuid)
                        continue
                    end
                    try
                        report = obj.iValidateSubjectNode(report, ...
                            subjectRecord, ProjectInfo, mode, scope);
                    catch ME
                        report = obj.iAddIssue(report, 'error', ...
                            'subject_validation_failed', 'subject', ...
                            subjectRecord.relativePath, ME.message, false);
                    end
                end
            end

            if isfield(ProjectInfo, 'subjectRegistry')
                report = obj.iValidateUnregisteredTopNodes(report, ProjectInfo);
            end

            if ~isScoped
                try
                    report = obj.iValidateGlobalUUIDs(report, ProjectInfo);
                catch ME
                    report = obj.iAddIssue(report, 'error', ...
                        'global_uuid_validation_failed', 'project', '', ...
                        ME.message, false);
                end
            end
            report = obj.iFinishReport(report, ProjectInfo, isScoped);
        end
    end

    methods (Access = private)
        function iCompensateFailedRebindDetach(obj, sessionUUID, ...
                subjectID, SessionInfo, stagingSessionPath, backupSubject, ...
                recoveryPath, targetProjectUUID, attachError)
            %ICOMPENSATEFAILEDREBINDDETACH Best-effort restore into source.
            %
            % Called only when rebindSessionToProject's phase 2 (attach to
            % the target project) fails after phase 1 (detach from this,
            % the source, project) already committed. Restores the session
            % folder and this project's subject registry entry from the
            % phase-1 recovery snapshot. SessionInfo is the pristine,
            % pre-attach struct (NOT the copy phase 2 mutated with the
            % target subject's UUID/ID) and is re-saved on top of the
            % restored folder so the session's own metadata cannot end up
            % claiming the target subject as its owner. Always errors: with
            % the original attach failure if the restore succeeds, or with
            % both errors combined (and the recovery path highlighted for
            % manual recovery) if the restore itself fails.

            restoreErrID = 'Umitoolbox:UMITProjectStore:rebindSessionFailed';
            sourceLock = obj.iAcquireWriteLock('rebindSessionToProject_undo');
            try
                % Phase 1 already removed the session from this subject, so
                % only the subject (not the session) can still be resolved
                % normally; the session's folder location is rebuilt from
                % the subject's own relativePath plus the known sessionID.
                [subjectRecord, ~, subjectPath] = ...
                    obj.iResolveSubject(subjectID);
                sessionRel = UMITProjectStore.iJoinRelative( ...
                    subjectRecord.relativePath, obj.Schema.folders.sessions, ...
                    SessionInfo.sessionID);
                sessionPath = obj.iResolveRelativePath(sessionRel);

                UMITProjectStore.iRemoveFolderIfPresent(sessionPath);
                copyfile(stagingSessionPath, sessionPath, 'f');
                obj.iSaveSessionInfo(sessionPath, SessionInfo);
                copyfile(backupSubject, ...
                    fullfile(subjectPath, obj.Schema.files.subjectMetadata), 'f');

                SourceProjectInfo = obj.getProjectInfo();
                boundRoles = {'rawDataFolder', 'processedDataFolder'};
                for iRole = 1:numel(boundRoles)
                    role = boundRoles{iRole};
                    if isempty(SessionInfo.(role)) || ...
                            ~isfolder(SessionInfo.(role))
                        continue
                    end
                    ProjectBinding = UMITProjectStore.readProjectBinding( ...
                        SessionInfo.(role));
                    ProjectBinding.projectUUID = ...
                        SourceProjectInfo.projectUUID;
                    ProjectBinding.subjectUUID = SessionInfo.subjectUUID;
                    obj.iInstallBinding(SessionInfo.(role), ProjectBinding, ...
                        SessionInfo.subjectID, SessionInfo.sessionID);
                end
                obj.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {SessionInfo.subjectUUID}, 'sessions', {sessionUUID}));

                clear sourceLock

            catch ME
                clear sourceLock
                warning(restoreErrID, ...
                    ['Could not attach session (UUID: %s) to target ' ...
                     'project %s, AND the automatic restore into this ' ...
                     'source project also failed. The session is not ' ...
                     'referenced by either project''s registry. Its files ' ...
                     'are preserved for manual recovery at: %s'], ...
                    sessionUUID, targetProjectUUID, recoveryPath);
                error(restoreErrID, ...
                    ['Rebind failed and automatic restore also failed. ' ...
                     'Attach cause: %s | Restore cause: %s | Recovery ' ...
                     'snapshot: %s'], attachError.message, ME.message, ...
                    recoveryPath);
            end

            UMITProjectStore.iRemoveFolderIfPresent(recoveryPath);
            error(restoreErrID, ...
                ['Could not attach session to target project %s. The ' ...
                 'session was restored to this (source) project. ' ...
                 'Cause: %s'], targetProjectUUID, attachError.message);
        end

        function obj = UMITProjectStore(projectRoot, schema, isReadOnly)
            %UMITPROJECTSTORE Construct an object bound to one project root.

            obj.ProjectRoot = UMITProjectStore.iAbsolutePath(projectRoot);
            obj.Schema = schema;
            obj.IsReadOnly = isReadOnly;
            obj.LastValidationReport = obj.iNewValidationReport('none');
        end

        function resourceUUID = iAddManagedResource(obj, ownerType, ownerIDs, resourceType, sourceFile, resourceInfo)
            %IADDMANAGEDRESOURCE Import one file into a managed resource folder.

            errID = 'Umitoolbox:UMITProjectStore:addResourceFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock(['add_', resourceType]);

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                error(errID, 'Unsupported resource type: %s', resourceType);
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            if ~strcmp(resourceDef.ownerType, ownerType)
                error(errID, ...
                    'Resource type "%s" cannot be owned by "%s".', ...
                    resourceType, ownerType);
            end

            if ~isstruct(resourceInfo) || ~isscalar(resourceInfo)
                error(errID, '"resourceInfo" must be a scalar struct.');
            end

            sourceFile = UMITProjectStore.iAbsolutePath(sourceFile);
            if ~isfile(sourceFile)
                error(errID, 'Source file does not exist: %s', sourceFile);
            end

            [~, ~, extension] = fileparts(sourceFile);
            if ~any(strcmpi(extension, resourceDef.allowedExtensions))
                error(errID, ...
                    'Resource type "%s" does not allow extension "%s".', ...
                    resourceType, extension);
            end

            % Verify that the managed MAT resource can be read and that its
            % domain payload matches the schema-defined resource contract.
            try
                sourceProbe = load(sourceFile, '-mat');
            catch ME
                error(errID, ...
                    'Source MAT file cannot be loaded: %s', ME.message);
            end

            obj.iValidateManagedResourcePayload( ...
                resourceType, sourceProbe, sourceFile);

            displayName = UMITProjectStore.iGetTextField( ...
                resourceInfo, 'displayName', resourceDef.filePrefix, true, errID);
            description = UMITProjectStore.iGetTextField( ...
                resourceInfo, 'description', '', true, errID);

            owner = obj.iResolveOwner(ownerType, ownerIDs);
            scope = UMITProjectStore.iLocationScope(owner);
            obj.iAssertHealthyForMutation(scope);
            metadata = owner.metadata;
            originalMetadata = metadata;

            if strcmp(resourceType, 'imageReference')
                obj.iAssertImageReferenceImport( ...
                    sourceProbe.ImageReference, owner, sourceFile);
            end

            resourceUUID = UMITProjectStore.iGenerateUUID();
            timestamp = char(datetime('now', ...
                'Format', 'yyyyMMdd''T''HHmmss'));
            compactUUID = strrep(resourceUUID, '-', '');
            fileName = sprintf('%s_%s_%s%s', ...
                resourceDef.filePrefix, timestamp, compactUUID(1:8), ...
                lower(extension));

            destinationRel = obj.iBuildResourceRelativePath( ...
                owner.baseRelativePath, resourceType, ...
                obj.Schema.folders.active, fileName);
            destinationPath = obj.iResolveRelativePath(destinationRel);

            transactionPath = obj.iCreateTransactionFolder(['add_', resourceType]);
            stagedFile = fullfile(transactionPath, fileName);
            cleanupStage = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(transactionPath));

            [ok, message] = copyfile(sourceFile, stagedFile, 'f');
            if ~ok
                error(errID, 'Could not stage resource: %s', message);
            end

            checksum = computeFileChecksum(stagedFile);
            nowTime = datetime('now');
            pointerField = resourceDef.activePointerField;

            if isempty(metadata.(pointerField))
                status = 'active';
                metadata.(pointerField) = resourceUUID;
            else
                status = 'available';
            end

            record = UMITProjectStore.iNewResourceRecord( ...
                resourceUUID, resourceType, fileName, destinationRel, ...
                displayName, description, nowTime, status, checksum, sourceFile);
            metadata.resourceRegistry(end+1) = record;
            metadata.modifiedOn = nowTime;

            [ok, message] = movefile(stagedFile, destinationPath, 'f');
            if ~ok
                error(errID, 'Could not install resource: %s', message);
            end

            try
                obj.iSaveOwnerMetadata(owner, metadata);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                if isfile(destinationPath)
                    delete(destinationPath);
                end
                obj.iSaveOwnerMetadata(owner, originalMetadata);
                rethrow(ME)
            end

            obj.iAppendLog(['add_', resourceType], resourceUUID, 'completed');
            clear cleanupStage lockCleanup
        end

        function iAssertImageReferenceImport(obj, ImageReference, owner, sourceFile)
            %IASSERTIMAGEREFERENCEIMPORT Version and identity of a new library reference.
            %
            %   Subject libraries store ImageReference version 2. Identity
            %   fields are optional in v2, but any that are present must
            %   name this project and the owner subject (audit item 21);
            %   a reference exported from another project or subject is
            %   imported only after its identity fields are cleared.

            if ImageReference.version ~= 2
                error('Umitoolbox:UMITProjectStore:imageReferenceVersion', ...
                    ['Subject libraries store ImageReference version 2, but "%s" ' ...
                     'holds version %d. Recreate it with genImageReferenceStruct.'], ...
                    sourceFile, ImageReference.version);
            end
            ProjectInfo = obj.getProjectInfo();
            payloadProject = strtrim(char(string(ImageReference.projectUUID)));
            if ~isempty(payloadProject) && ~strcmpi(payloadProject, ProjectInfo.projectUUID)
                error('Umitoolbox:UMITProjectStore:imageReferenceIdentityMismatch', ...
                    ['The image reference "%s" belongs to project "%s" (%s), not to ' ...
                     'project "%s" (%s).'], sourceFile, ImageReference.projectName, ...
                    payloadProject, ProjectInfo.projectName, ProjectInfo.projectUUID);
            end
            payloadSubject = strtrim(char(string(ImageReference.subjectUUID)));
            if ~isempty(payloadSubject) && ~strcmpi(payloadSubject, owner.metadata.uuid)
                error('Umitoolbox:UMITProjectStore:imageReferenceIdentityMismatch', ...
                    ['The image reference "%s" belongs to subject "%s" (%s), not to ' ...
                     'subject "%s" (%s).'], sourceFile, ImageReference.subjectID, ...
                    payloadSubject, owner.metadata.subjectID, owner.metadata.uuid);
            end
        end

        function iValidateManagedResourcePayload(obj, resourceType, loadedData, sourceFile)
            %IVALIDATEMANAGEDRESOURCEPAYLOAD Validate schema-defined MAT payload.
            %
            %   The project store performs this check before importing a file.
            %   Full project validation repeats the same check for managed
            %   resources whose schema defines a payload contract.

            errID = 'Umitoolbox:UMITProjectStore:invalidResourcePayload';

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                error(errID, 'Unsupported resource type: %s', resourceType);
            end

            if ~isstruct(loadedData) || ~isscalar(loadedData)
                error(errID, ...
                    'Loaded MAT-file content must be represented by a scalar struct.');
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            requiredVariable = '';
            validatorFunction = '';

            if isfield(resourceDef, 'requiredVariable')
                requiredVariable = char(string(resourceDef.requiredVariable));
            end
            if isfield(resourceDef, 'validatorFunction')
                validatorFunction = char(string(resourceDef.validatorFunction));
            end

            if ~isempty(requiredVariable) && ...
                    ~isfield(loadedData, requiredVariable)
                error(errID, ...
                    ['Resource type "%s" requires variable "%s" in MAT-file ' ...
                     '"%s".'], ...
                    resourceType, requiredVariable, sourceFile);
            end

            if isempty(validatorFunction)
                return
            end

            if exist(validatorFunction, 'file') ~= 2
                error(errID, ...
                    ['Validator "%s" required by resource type "%s" is not ' ...
                     'available on the MATLAB path.'], ...
                    validatorFunction, resourceType);
            end

            if isempty(requiredVariable)
                payload = loadedData;
            else
                payload = loadedData.(requiredVariable);
            end

            try
                feval(validatorFunction, payload);
            catch ME
                wrapped = MException(errID, ...
                    ['Resource type "%s" failed payload validation for "%s".'], ...
                    resourceType, sourceFile);
                wrapped = addCause(wrapped, ME);
                throw(wrapped);
            end
        end

        function tf = iResourceTypeHasPayloadContract(obj, resourceType)
            %IRESOURCETYPEHASPAYLOADCONTRACT True for validated MAT payloads.

            tf = false;

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                return
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);

            if isfield(resourceDef, 'requiredVariable')
                tf = tf || strlength(string(resourceDef.requiredVariable)) > 0;
            end
            if isfield(resourceDef, 'validatorFunction')
                tf = tf || strlength(string(resourceDef.validatorFunction)) > 0;
            end
        end

        function iSetActiveResource(obj, ownerType, ownerIDs, resourceType, resourceUUID)
            %ISETACTIVERESOURCE Select one non-archived resource as active.

            errID = 'Umitoolbox:UMITProjectStore:setActiveResourceFailed';
            resourceUUID = char(string(resourceUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock(['setActive_', resourceType]);

            owner = obj.iResolveOwner(ownerType, ownerIDs);
            scope = UMITProjectStore.iLocationScope(owner);
            obj.iAssertHealthyForMutation(scope);
            metadata = owner.metadata;
            originalMetadata = metadata;
            resourceIndex = obj.iFindResourceIndex( ...
                metadata.resourceRegistry, resourceUUID);

            if isempty(resourceIndex)
                error(errID, ...
                    'Resource is not owned by the requested project node.');
            end

            selectedRecord = metadata.resourceRegistry(resourceIndex);
            if ~strcmp(selectedRecord.type, resourceType)
                error(errID, ...
                    'Resource type does not match "%s".', resourceType);
            end
            if strcmp(selectedRecord.status, 'archived')
                error(errID, 'Archived resources cannot be made active.');
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            pointerField = resourceDef.activePointerField;
            oldUUID = metadata.(pointerField);

            if ~isempty(oldUUID)
                oldIndex = obj.iFindResourceIndex( ...
                    metadata.resourceRegistry, oldUUID);
                if ~isempty(oldIndex) && oldIndex ~= resourceIndex && ...
                        ~strcmp(metadata.resourceRegistry(oldIndex).status, 'archived')
                    metadata.resourceRegistry(oldIndex).status = 'available';
                    metadata.resourceRegistry(oldIndex).modifiedOn = datetime('now');
                end
            end

            metadata.resourceRegistry(resourceIndex).status = 'active';
            metadata.resourceRegistry(resourceIndex).modifiedOn = datetime('now');
            metadata.(pointerField) = resourceUUID;
            metadata.modifiedOn = datetime('now');

            try
                obj.iSaveOwnerMetadata(owner, metadata);
                obj.iAssertValidAfterMutation(scope);
            catch ME
                obj.iSaveOwnerMetadata(owner, originalMetadata);
                rethrow(ME)
            end
            obj.iAppendLog(['setActive_', resourceType], resourceUUID, 'completed');
            clear lockCleanup
        end


        function resources = iListResources(obj, metadata, ownerType, ownerUUID, varargin)
            %ILISTRESOURCES Filter and resolve resources owned by one node.

            p = inputParser;
            p.FunctionName = 'UMITProjectStore.listResources';
            addParameter(p, 'Type', '', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));
            addParameter(p, 'Status', {'active', 'available'});
            addParameter(p, 'VerifyFiles', false, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            if ~isfield(metadata, 'resourceRegistry') || ...
                    ~obj.iIsResourceRegistryStruct(metadata.resourceRegistry)
                error('Umitoolbox:UMITProjectStore:invalidResourceRegistry', ...
                    'Owner metadata contains an invalid resourceRegistry.');
            end

            typeFilter = char(string(p.Results.Type));
            resourceTypeNames = fieldnames(obj.Schema.resourceTypes);
            if ~isempty(typeFilter)
                typeIndex = find(strcmpi( ...
                    resourceTypeNames, typeFilter), 1, 'first');
                if isempty(typeIndex)
                    error('Umitoolbox:UMITProjectStore:invalidResourceFilter', ...
                        'Unsupported resource type filter: %s', typeFilter);
                end

                typeFilter = resourceTypeNames{typeIndex};
                if ~strcmp( ...
                        obj.Schema.resourceTypes.(typeFilter).ownerType, ...
                        ownerType)
                    error('Umitoolbox:UMITProjectStore:invalidResourceFilter', ...
                        ['Resource type "%s" cannot be owned by a %s ' ...
                        'project node.'], typeFilter, ownerType);
                end
            end

            statusFilter = UMITProjectStore.iNormalizeTextList( ...
                p.Results.Status, 'Status');
            if isempty(statusFilter)
                statusFilter = obj.Schema.resourceStates;
            else
                canonicalStatus = cell(size(statusFilter));
                for iStatus = 1:numel(statusFilter)
                    statusIndex = find(strcmpi( ...
                        obj.Schema.resourceStates, statusFilter{iStatus}), ...
                        1, 'first');
                    if isempty(statusIndex)
                        error( ...
                            'Umitoolbox:UMITProjectStore:invalidResourceFilter', ...
                            'Unsupported resource status filter: %s', ...
                            statusFilter{iStatus});
                    end
                    canonicalStatus{iStatus} = ...
                        obj.Schema.resourceStates{statusIndex};
                end
                statusFilter = unique(canonicalStatus, 'stable');
            end

            registry = metadata.resourceRegistry;
            selected = true(1, numel(registry));
            if ~isempty(registry)
                % Retired types (e.g. registration transforms) are never listed.
                selected = ~ismember({registry.type}, obj.Schema.retiredResourceTypes);
            end

            if ~isempty(typeFilter)
                selected = selected & strcmp({registry.type}, typeFilter);
            end
            selected = selected & ismember({registry.status}, statusFilter);
            selectedIndices = find(selected);

            resources = UMITProjectStore.iEmptyResolvedResourceRegistry();
            for iResource = 1:numel(selectedIndices)
                record = registry(selectedIndices(iResource));
                resource = obj.iEnrichResourceRecord( ...
                    record, ownerType, ownerUUID);

                if p.Results.VerifyFiles && ~resource.fileExists
                    error('Umitoolbox:UMITProjectStore:missingResourceFile', ...
                        'Registered resource file is missing: %s', ...
                        resource.absolutePath);
                end

                resources(end+1, 1) = resource; %#ok<AGROW>
            end
        end

        function resource = iEnrichResourceRecord(obj, record, ownerType, ownerUUID)
            %IENRICHRESOURCERECORD Add resolved owner and filesystem fields.

            resource = record;
            resource.ownerType = ownerType;
            resource.ownerUUID = ownerUUID;
            resource.absolutePath = ...
                obj.iResolveRelativePath(record.relativePath);
            resource.fileExists = isfile(resource.absolutePath);
        end

        function resource = iGetActiveResource(obj, owner, resourceType)
            %IGETACTIVERESOURCE Resolve and validate one active resource.

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                error('Umitoolbox:UMITProjectStore:invalidResourceType', ...
                    'Unsupported resource type: %s', resourceType);
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            if ~strcmp(resourceDef.ownerType, owner.type)
                error('Umitoolbox:UMITProjectStore:invalidResourceOwner', ...
                    'Resource type "%s" cannot belong to owner type "%s".', ...
                    resourceType, owner.type);
            end

            metadata = owner.metadata;
            pointerField = resourceDef.activePointerField;
            if ~isfield(metadata, pointerField)
                error('Umitoolbox:UMITProjectStore:invalidActiveResource', ...
                    'Owner metadata is missing active pointer "%s".', ...
                    pointerField);
            end

            pointerUUID = metadata.(pointerField);
            if ~(ischar(pointerUUID) || ...
                    (isstring(pointerUUID) && isscalar(pointerUUID)))
                error('Umitoolbox:UMITProjectStore:invalidActiveResource', ...
                    'Active pointer "%s" must be a text scalar.', ...
                    pointerField);
            end
            pointerUUID = char(string(pointerUUID));

            if isempty(pointerUUID)
                resource = [];
                return
            end

            registry = metadata.resourceRegistry;
            matchingIndices = find(strcmpi({registry.uuid}, pointerUUID));
            if isempty(matchingIndices)
                error('Umitoolbox:UMITProjectStore:invalidActiveResource', ...
                    ['Active pointer "%s" does not resolve to a registered ' ...
                    'resource.'], pointerField);
            end
            if numel(matchingIndices) > 1
                error('Umitoolbox:UMITProjectStore:duplicateResourceUUID', ...
                    'Active resource UUID is registered more than once: %s', ...
                    pointerUUID);
            end

            record = registry(matchingIndices);
            if ~strcmp(record.type, resourceType) || ...
                    ~strcmp(record.status, 'active')
                error('Umitoolbox:UMITProjectStore:invalidActiveResource', ...
                    ['Active pointer "%s" does not reference an active ' ...
                    'resource of type "%s".'], ...
                    pointerField, resourceType);
            end

            resource = obj.iEnrichResourceRecord( ...
                record, owner.type, metadata.uuid);
            if ~resource.fileExists
                error('Umitoolbox:UMITProjectStore:missingResourceFile', ...
                    'Active resource file is missing: %s', ...
                    resource.absolutePath);
            end
        end

        function locations = iFindResourceLocations(obj, resourceUUID)
            %IFINDRESOURCELOCATIONS Find all registrations of one resource UUID.
            %   Unreadable subject or session metadata is skipped: one broken
            %   node must not block resource lookups elsewhere in the project
            %   (per-session validity). Resources of a broken node are then
            %   not found.

            resourceUUID = UMITProjectStore.iNormalizeUUIDInput(resourceUUID);
            template = UMITProjectStore.iNewResourceLocation( ...
                '', struct(), 0, '', '', '', '');
            locations = repmat(template, 0, 1);
            ProjectInfo = obj.getProjectInfo();

            for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                subjectRecord = ProjectInfo.subjectRegistry(iSubject);
                subjectPath = obj.iResolveRelativePath( ...
                    subjectRecord.relativePath);
                try
                    SubjectInfo = obj.iLoadMetadata( ...
                        fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                        obj.Schema.metadataVariables.subject);
                    indices = find(strcmpi( ...
                        {SubjectInfo.resourceRegistry.uuid}, resourceUUID));
                catch
                    continue
                end
                for iIndex = 1:numel(indices)
                    locations(end+1, 1) = ...
                        UMITProjectStore.iNewResourceLocation( ...
                        'subject', SubjectInfo, indices(iIndex), ...
                        subjectPath, subjectRecord.relativePath, ...
                        fullfile(subjectPath, ...
                        obj.Schema.files.subjectMetadata), ...
                        obj.Schema.metadataVariables.subject); %#ok<AGROW>
                end

                for iSession = 1:numel(SubjectInfo.sessionRegistry)
                    sessionRecord = SubjectInfo.sessionRegistry(iSession);
                    sessionPath = obj.iResolveRelativePath( ...
                        sessionRecord.relativePath);
                    try
                        SessionInfo = obj.iLoadMetadata( ...
                            fullfile(sessionPath, ...
                            obj.Schema.files.sessionMetadata), ...
                            obj.Schema.metadataVariables.session);
                        indices = find(strcmpi( ...
                            {SessionInfo.resourceRegistry.uuid}, resourceUUID));
                    catch
                        continue
                    end
                    for iIndex = 1:numel(indices)
                        locations(end+1, 1) = ...
                            UMITProjectStore.iNewResourceLocation( ...
                            'session', SessionInfo, indices(iIndex), ...
                            sessionPath, sessionRecord.relativePath, ...
                            fullfile(sessionPath, ...
                            obj.Schema.files.sessionMetadata), ...
                            obj.Schema.metadataVariables.session); %#ok<AGROW>
                    end
                end
            end

        end

        function owner = iResolveOwner(obj, ownerType, ownerIDs)
            %IRESOLVEOWNER Resolve owner metadata and canonical base path.

            owner = struct();
            owner.type = ownerType;

            switch ownerType
                case 'subject'
                    [record, metadata, path] = obj.iResolveSubject(ownerIDs{1});
                    owner.record = record;
                    owner.metadata = metadata;
                    owner.path = path;
                    owner.metadataPath = fullfile(path, obj.Schema.files.subjectMetadata);
                    owner.variableName = obj.Schema.metadataVariables.subject;
                    owner.baseRelativePath = record.relativePath;

                case 'session'
                    [~, ~, ~, metadata, path, ~, record] = ...
                        obj.iResolveSession(ownerIDs{1}, ownerIDs{2});
                    owner.record = record;
                    owner.metadata = metadata;
                    owner.path = path;
                    owner.metadataPath = fullfile(path, obj.Schema.files.sessionMetadata);
                    owner.variableName = obj.Schema.metadataVariables.session;
                    owner.baseRelativePath = record.relativePath;

                otherwise
                    error('Umitoolbox:UMITProjectStore:invalidOwnerType', ...
                        'Unsupported owner type: %s', ownerType);
            end
        end

        function [record, SubjectInfo, subjectPath, index, ProjectInfo] = iResolveSubject(obj, subjectID)
            %IRESOLVESUBJECT Resolve a subject registry entry and metadata.

            subjectID = char(string(subjectID));
            ProjectInfo = obj.getProjectInfo();
            index = obj.iFindRegistryIndex(ProjectInfo.subjectRegistry, subjectID);
            if isempty(index)
                subjectID = obj.iNormalizeManagedID(subjectID, 'subjectID');
                index = obj.iFindRegistryIndex( ...
                    ProjectInfo.subjectRegistry, subjectID);
            end
            if isempty(index)
                error('Umitoolbox:UMITProjectStore:subjectNotFound', ...
                    'Subject ID was not found: %s', subjectID);
            end

            record = ProjectInfo.subjectRegistry(index);
            subjectPath = obj.iResolveRelativePath(record.relativePath);
            SubjectInfo = obj.iLoadMetadata( ...
                fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                obj.Schema.metadataVariables.subject);
        end

        function [subjectRecord, SubjectInfo, subjectPath, SessionInfo, sessionPath, sessionIndex, sessionRecord] = iResolveSession(obj, subjectID, sessionID)
            %IRESOLVESESSION Resolve a session registry entry and metadata.

            [subjectRecord, SubjectInfo, subjectPath] = ...
                obj.iResolveSubject(subjectID);
            sessionID = char(string(sessionID));
            sessionIndex = obj.iFindRegistryIndex( ...
                SubjectInfo.sessionRegistry, sessionID);
            if isempty(sessionIndex)
                sessionID = obj.iNormalizeManagedID(sessionID, 'sessionID');
                sessionIndex = obj.iFindRegistryIndex( ...
                    SubjectInfo.sessionRegistry, sessionID);
            end
            if isempty(sessionIndex)
                error('Umitoolbox:UMITProjectStore:sessionNotFound', ...
                    'Session ID was not found: %s', sessionID);
            end

            sessionRecord = SubjectInfo.sessionRegistry(sessionIndex);
            sessionPath = obj.iResolveRelativePath(sessionRecord.relativePath);
            SessionInfo = obj.iLoadMetadata( ...
                fullfile(sessionPath, obj.Schema.files.sessionMetadata), ...
                obj.Schema.metadataVariables.session);
        end


        function [subjectRecord, SubjectInfo, subjectPath, SessionInfo, ...
                sessionPath, sessionIndex, sessionRecord] = ...
                iResolveSessionByUUID(obj, sessionUUID)
            %IRESOLVESESSIONBYUUID Resolve a unique session UUID project-wide.

            sessionUUID = UMITProjectStore.iNormalizeUUIDInput(sessionUUID);
            ProjectInfo = obj.getProjectInfo();
            matches = struct( ...
                'subjectRecord', {}, ...
                'SubjectInfo', {}, ...
                'subjectPath', {}, ...
                'SessionInfo', {}, ...
                'sessionPath', {}, ...
                'sessionIndex', {}, ...
                'sessionRecord', {});

            for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                thisSubjectRecord = ProjectInfo.subjectRegistry(iSubject);
                thisSubjectPath = obj.iResolveRelativePath( ...
                    thisSubjectRecord.relativePath);
                thisSubjectInfo = obj.iLoadMetadata( ...
                    fullfile(thisSubjectPath, ...
                    obj.Schema.files.subjectMetadata), ...
                    obj.Schema.metadataVariables.subject);

                if isempty(thisSubjectInfo.sessionRegistry)
                    continue
                end

                idx = find(strcmpi( ...
                    {thisSubjectInfo.sessionRegistry.uuid}, ...
                    sessionUUID));
                for iMatch = idx
                    thisSessionRecord = ...
                        thisSubjectInfo.sessionRegistry(iMatch);
                    thisSessionPath = obj.iResolveRelativePath( ...
                        thisSessionRecord.relativePath);
                    thisSessionInfo = obj.iLoadMetadata( ...
                        fullfile(thisSessionPath, ...
                        obj.Schema.files.sessionMetadata), ...
                        obj.Schema.metadataVariables.session);

                    match = struct();
                    match.subjectRecord = thisSubjectRecord;
                    match.SubjectInfo = thisSubjectInfo;
                    match.subjectPath = thisSubjectPath;
                    match.SessionInfo = thisSessionInfo;
                    match.sessionPath = thisSessionPath;
                    match.sessionIndex = iMatch;
                    match.sessionRecord = thisSessionRecord;
                    matches(end+1) = match; %#ok<AGROW>
                end
            end

            if isempty(matches)
                error('Umitoolbox:UMITProjectStore:sessionNotFound', ...
                    'Session UUID was not found: %s', sessionUUID);
            end
            if numel(matches) > 1
                error('Umitoolbox:UMITProjectStore:duplicateSessionUUID', ...
                    'Session UUID is registered more than once: %s', ...
                    sessionUUID);
            end

            subjectRecord = matches.subjectRecord;
            SubjectInfo = matches.SubjectInfo;
            subjectPath = matches.subjectPath;
            SessionInfo = matches.SessionInfo;
            sessionPath = matches.sessionPath;
            sessionIndex = matches.sessionIndex;
            sessionRecord = matches.sessionRecord;
        end

        function ProjectBinding = iBindSessionDataFolder( ...
                obj, subjectID, sessionID, dataFolder, folderRole, varargin)
            %IBINDSESSIONDATAFOLDER Create reciprocal session/project-ledger state.

            errID = 'Umitoolbox:UMITProjectStore:bindDataFolderFailed';

            p = inputParser;
            p.FunctionName = 'UMITProjectStore.bindDataFolder';
            addParameter(p, 'ReplaceOrphanBinding', false, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});
            replaceOrphanBinding = p.Results.ReplaceOrphanBinding;

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock( ...
                ['bind_', folderRole]);

            if ~(ischar(dataFolder) || ...
                    (isstring(dataFolder) && isscalar(dataFolder)))
                error(errID, ...
                    '"dataFolder" must be a character vector or string scalar.');
            end

            dataFolder = UMITProjectStore.iAbsolutePath(dataFolder);
            if ~isfolder(dataFolder)
                error(errID, ...
                    'Data folder does not exist: %s', dataFolder);
            end
            obj.iAssertExternalDataPath(dataFolder, folderRole);
            if strcmp(folderRole, 'processedDataFolder')
                obj.iAssertValidSaveFolderDataset(dataFolder);
            else
                obj.iAssertValidRawDataFolder(dataFolder);
            end

            [~, SubjectInfo, ~, SessionInfo, sessionPath] = ...
                obj.iResolveSession(subjectID, sessionID);
            scope = UMITProjectStore.iScope('subjects', {SubjectInfo.uuid}, ...
                'sessions', {SessionInfo.uuid});
            obj.iAssertHealthyForMutation(scope);
            if strcmp(folderRole, 'rawDataFolder')
                obj.iAssertRawFolderAuthority(SessionInfo, dataFolder);
            end
            [pathField, bindingField] = ...
                obj.iGetBindingRoleFields(folderRole);

            hasPath = ~isempty(SessionInfo.(pathField));
            hasBindingUUID = ~isempty(SessionInfo.(bindingField));

            if hasPath || hasBindingUUID
                if hasPath && hasBindingUUID && strcmp( ...
                        UMITProjectStore.iNormalizeComparisonPath( ...
                        SessionInfo.(pathField)), ...
                        UMITProjectStore.iNormalizeComparisonPath( ...
                        dataFolder))
                    ProjectBinding = ...
                        obj.iGetSessionDataFolderBinding( ...
                        subjectID, sessionID, folderRole);
                    clear lockCleanup
                    return
                end

                error('Umitoolbox:UMITProjectStore:dataFolderAlreadyBound', ...
                    ['Session %s is already bound through %s. Unbind it ' ...
                     'before assigning another folder.'], ...
                    SessionInfo.sessionID, pathField);
            end

            [usesSharedBinding, sharedBinding] = ...
                obj.iGetCoLocatedSessionBinding( ...
                dataFolder, SessionInfo, SubjectInfo, pathField);
            if usesSharedBinding
                ProjectBinding = obj.iAdoptCoLocatedSessionBinding( ...
                    SessionInfo, sessionPath, pathField, bindingField, ...
                    dataFolder, sharedBinding, 'bindDataFolder');
                clear lockCleanup
                return
            end

            existingMatches = obj.iFindSessionPathMatches(dataFolder);
            if ~isempty(existingMatches)
                error('Umitoolbox:UMITProjectStore:dataFolderAlreadyBound', ...
                    ['The data folder is already registered to session %s ' ...
                     'through %s.'], ...
                    existingMatches(1).sessionID, ...
                    existingMatches(1).matchedField);
            end

            hasOrphanBinding = false;

            if UMITProjectStore.iHasBinding(dataFolder)
                if ~replaceOrphanBinding
                    error('Umitoolbox:UMITProjectStore:bindingFileExists', ...
                        ['The target folder is already bound to a project (%s). ' ...
                         'Resolve or detach that binding before assigning this folder.'], ...
                        obj.Schema.files.projectBinding);
                end

                [hasOrphanBinding, ~] = ...
                    UMITProjectStore.isOrphanProjectBinding(dataFolder);

                if ~hasOrphanBinding
                    error('Umitoolbox:UMITProjectStore:bindingProjectAvailable', ...
                        ['The existing binding references an available project. ' ...
                         'Use the owning store''s unbind method instead.']);
                end
            end

            ProjectInfo = obj.getProjectInfo();
            nowTime = datetime('now');
            ProjectBinding = struct();
            ProjectBinding.version = ...
                obj.Schema.projectBinding.version;
            ProjectBinding.bindingUUID = ...
                UMITProjectStore.iGenerateUUID();
            ProjectBinding.projectUUID = ...
                ProjectInfo.projectUUID;
            ProjectBinding.subjectUUID = SubjectInfo.uuid;
            ProjectBinding.sessionUUID = SessionInfo.uuid;
            ProjectBinding.folderRole = pathField;
            ProjectBinding.createdOn = nowTime;
            ProjectBinding.modifiedOn = nowTime;

            UMITProjectStore.iValidateProjectBindingStruct( ...
                ProjectBinding, obj.Schema, errID);

            recoveryPath = ...
                obj.iCreateRecoveryFolder(['bind_', folderRole]);
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));
            sessionMetadataPath = fullfile(sessionPath, ...
                obj.Schema.files.sessionMetadata);
            backupSessionPath = fullfile(recoveryPath, ...
                obj.Schema.files.sessionMetadata);
            copyfile(sessionMetadataPath, backupSessionPath, 'f');

            % An orphan binding is replaced by a rebind event (history kept).
            ledgerBackup = UMITProjectStore.iBackupLedger(dataFolder);

            try
                SessionInfo.(pathField) = dataFolder;
                SessionInfo.(bindingField) = ...
                    ProjectBinding.bindingUUID;
                SessionInfo.modifiedOn = datetime('now');
                obj.iSaveSessionInfo(sessionPath, SessionInfo);

                obj.iInstallBinding(dataFolder, ProjectBinding, ...
                    SubjectInfo.subjectID, SessionInfo.sessionID);

                obj.iAssertValidAfterMutation(scope);

            catch ME
                copyfile(backupSessionPath, ...
                    sessionMetadataPath, 'f');
                UMITProjectStore.iRestoreLedger(ledgerBackup);
                rethrow(ME)
            end

            if hasOrphanBinding
                logAction = 'replaceOrphanDataFolderBinding';
            else
                logAction = 'bindDataFolder';
            end

            obj.iAppendLog(logAction, ...
                ProjectBinding.bindingUUID, pathField);
            clear cleanupRecovery lockCleanup
        end

        function ProjectBinding = iGetSessionDataFolderBinding( ...
                obj, subjectID, sessionID, folderRole)
            %IGETSESSIONDATAFOLDERBINDING Load and verify one reciprocal binding.

            [~, SubjectInfo, ~, SessionInfo] = ...
                obj.iResolveSession(subjectID, sessionID);
            [pathField, bindingField] = ...
                obj.iGetBindingRoleFields(folderRole);

            hasPath = ~isempty(SessionInfo.(pathField));
            hasBindingUUID = ~isempty(SessionInfo.(bindingField));

            if ~hasPath && ~hasBindingUUID
                ProjectBinding = [];
                return
            end
            if xor(hasPath, hasBindingUUID)
                error('Umitoolbox:UMITProjectStore:incompleteBinding', ...
                    ['Session binding metadata is incomplete for %s.'], ...
                    pathField);
            end
            if ~isfolder(SessionInfo.(pathField))
                error('Umitoolbox:UMITProjectStore:bindingFolderUnavailable', ...
                    'Bound data folder is unavailable: %s', ...
                    SessionInfo.(pathField));
            end

            ProjectBinding = UMITProjectStore.readProjectBinding( ...
                SessionInfo.(pathField));
            obj.iAssertBindingMatchesSession( ...
                ProjectBinding, SessionInfo, SubjectInfo, ...
                pathField, bindingField, SessionInfo.(pathField));
        end

        function status = iGetBindingRuntimeStatus( ...
                obj, SessionInfo, SubjectInfo, folderRole)
            %IGETBINDINGRUNTIMESTATUS Classify availability without mutation.

            [pathField, bindingField] = ...
                obj.iGetBindingRoleFields(folderRole);
            status = struct( ...
                'state', 'invalid', ...
                'path', SessionInfo.(pathField), ...
                'isAvailable', false, ...
                'message', '', ...
                'binding', []);

            if isempty(SessionInfo.(pathField)) || ...
                    isempty(SessionInfo.(bindingField))
                status.message = ...
                    'The persisted binding metadata is incomplete.';
                return
            end
            if ~isfolder(SessionInfo.(pathField))
                status.state = 'missing';
                status.message = sprintf( ...
                    'SaveFolder is currently unavailable: %s', ...
                    SessionInfo.(pathField));
                return
            end

            if strcmp(UMITProjectStore.iBindingState(SessionInfo.(pathField)), 'none')
                status.message = sprintf( ...
                    'SaveFolder is not bound in its %s.', ...
                    obj.Schema.files.projectBinding);
                return
            end

            try
                ProjectBinding = UMITProjectStore.readProjectBinding( ...
                    SessionInfo.(pathField));
                status.binding = ProjectBinding;
            catch ME
                status.message = ME.message;
                return
            end

            try
                obj.iAssertBindingMatchesSession( ...
                    ProjectBinding, SessionInfo, SubjectInfo, ...
                    pathField, bindingField, SessionInfo.(pathField));
            catch ME
                status.state = 'conflicting';
                status.message = ME.message;
                return
            end

            status.state = 'available';
            status.isAvailable = true;
            status.message = 'SaveFolder and project binding are available.';
        end

        function ProjectBinding = iRelocateSessionDataFolder( ...
                obj, subjectID, sessionID, newDataFolder, folderRole, varargin)
            %IRELOCATESESSIONDATAFOLDER Atomically relocate or repair a binding.

            if strcmp(folderRole, 'rawDataFolder')
                errID = ...
                    'Umitoolbox:UMITProjectStore:relocateRawDataFolderFailed';
                functionName = 'UMITProjectStore.relocateRawDataFolder';
                logOperation = 'relocateRawDataFolder';
            else
                errID = ...
                    'Umitoolbox:UMITProjectStore:relocateSaveFolderFailed';
                functionName = 'UMITProjectStore.relocateSaveFolder';
                logOperation = 'relocateSaveFolder';
            end
            p = inputParser;
            p.FunctionName = functionName;
            addParameter(p, 'ReplaceOrphanBinding', false, ...
                @(x) islogical(x) && isscalar(x));
            addParameter(p, 'RepairExisting', false, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            lockCleanup = obj.iAcquireWriteLock( ...
                ['relocate_', folderRole]);
            [~, SubjectInfo, ~, SessionInfo, sessionPath, ...
                sessionIndex] = obj.iResolveSession(subjectID, sessionID);
            scope = UMITProjectStore.iScope('subjects', {SubjectInfo.uuid}, ...
                'sessions', {SessionInfo.uuid});
            obj.iAssertBindingMutationAllowed(scope, ...
                SubjectInfo.sessionRegistry(sessionIndex).relativePath);

            if ~(ischar(newDataFolder) || ...
                    (isstring(newDataFolder) && isscalar(newDataFolder)))
                error(errID, ...
                    '"newSaveFolder" must be a character vector or string scalar.');
            end
            newDataFolder = UMITProjectStore.iAbsolutePath(newDataFolder);
            if ~isfolder(newDataFolder)
                error('Umitoolbox:UMITProjectStore:saveFolderUnavailable', ...
                    'Target SaveFolder does not exist: %s', newDataFolder);
            end
            obj.iAssertExternalDataPath(newDataFolder, folderRole);
            if strcmp(folderRole, 'processedDataFolder')
                obj.iAssertValidSaveFolderDataset(newDataFolder);
            else
                obj.iAssertValidRawDataFolder(newDataFolder);
            end

            if strcmp(folderRole, 'rawDataFolder')
                obj.iAssertRawFolderAuthority(SessionInfo, newDataFolder);
            end

            [pathField, bindingField] = ...
                obj.iGetBindingRoleFields(folderRole);
            oldDataFolder = SessionInfo.(pathField);
            sameFolder = ~isempty(oldDataFolder) && strcmp( ...
                UMITProjectStore.iNormalizeComparisonPath(oldDataFolder), ...
                UMITProjectStore.iNormalizeComparisonPath(newDataFolder));

            matches = obj.iFindSessionPathMatches(newDataFolder);
            if ~isempty(matches)
                isOwnMatch = strcmpi({matches.sessionUUID}, SessionInfo.uuid);
                if any(~isOwnMatch)
                    conflict = matches(find(~isOwnMatch, 1));
                    error('Umitoolbox:UMITProjectStore:dataFolderAlreadyBound', ...
                        ['The target SaveFolder is registered to session %s ' ...
                         'through %s.'], conflict.sessionID, ...
                        conflict.matchedField);
                end
            end

            targetBindingExists = UMITProjectStore.iHasBinding(newDataFolder);
            targetBinding = [];
            targetBindingReadable = false;
            if targetBindingExists
                try
                    targetBinding = UMITProjectStore.readProjectBinding( ...
                        newDataFolder);
                    targetBindingReadable = true;
                catch ME
                    if ~(sameFolder && p.Results.RepairExisting)
                        error('Umitoolbox:UMITProjectStore:bindingFileExists', ...
                            ['Target SaveFolder contains an invalid %s: %s'], ...
                            obj.Schema.files.projectBinding, ME.message);
                    end
                end
            end

            ProjectInfo = obj.getProjectInfo();
            ownsTargetBinding = targetBindingReadable && ...
                strcmpi(targetBinding.projectUUID, ProjectInfo.projectUUID) && ...
                strcmpi(targetBinding.subjectUUID, SubjectInfo.uuid) && ...
                strcmpi(targetBinding.sessionUUID, SessionInfo.uuid) && ...
                strcmp(targetBinding.folderRole, pathField);

            [usesSharedBinding, sharedBinding] = ...
                obj.iGetCoLocatedSessionBinding( ...
                newDataFolder, SessionInfo, SubjectInfo, pathField);
            if usesSharedBinding
                ProjectBinding = obj.iAdoptCoLocatedSessionBinding( ...
                    SessionInfo, sessionPath, pathField, bindingField, ...
                    newDataFolder, sharedBinding, logOperation);
                clear lockCleanup
                return
            end

            if targetBindingReadable && ~ownsTargetBinding
                replaceOrphan = false;
                if p.Results.ReplaceOrphanBinding
                    [replaceOrphan, ~] = ...
                        UMITProjectStore.isOrphanProjectBinding(newDataFolder);
                end
                if ~replaceOrphan
                    error('Umitoolbox:UMITProjectStore:bindingConflict', ...
                        ['Target SaveFolder is bound to another project or ' ...
                         'session. Its link was not modified.']);
                end
                % Audit item 13: taking over an orphan folder must not pull
                % a session away from a folder that still binds it.
                if ~isempty(oldDataFolder) && ~sameFolder && isfolder(oldDataFolder) && ...
                        strcmp(describeProjectBinding(oldDataFolder).state, 'bound')
                    error('Umitoolbox:UMITProjectStore:sessionAlreadyBound', ...
                        ['Session %s is healthily bound to "%s". Unbind it there, or ' ...
                         'choose another session, before giving it this folder.'], ...
                        SessionInfo.sessionID, oldDataFolder);
                end
            end

            if sameFolder && ownsTargetBinding && ...
                    strcmpi(targetBinding.bindingUUID, ...
                    SessionInfo.(bindingField))
                ProjectBinding = targetBinding;
                clear lockCleanup
                return
            end

            nowTime = datetime('now');
            ProjectBinding = struct();
            ProjectBinding.version = obj.Schema.projectBinding.version;
            if sameFolder && ~isempty(SessionInfo.(bindingField))
                ProjectBinding.bindingUUID = SessionInfo.(bindingField);
            else
                ProjectBinding.bindingUUID = ...
                    UMITProjectStore.iGenerateUUID();
            end
            ProjectBinding.projectUUID = ProjectInfo.projectUUID;
            ProjectBinding.subjectUUID = SubjectInfo.uuid;
            ProjectBinding.sessionUUID = SessionInfo.uuid;
            ProjectBinding.folderRole = pathField;
            ProjectBinding.createdOn = nowTime;
            if ownsTargetBinding
                ProjectBinding.createdOn = targetBinding.createdOn;
            end
            ProjectBinding.modifiedOn = nowTime;
            UMITProjectStore.iValidateProjectBindingStruct( ...
                ProjectBinding, obj.Schema, errID);

            recoveryPath = obj.iCreateRecoveryFolder( ...
                ['relocate_', folderRole]);
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));
            sessionMetadataPath = fullfile(sessionPath, ...
                obj.Schema.files.sessionMetadata);
            backupSessionPath = fullfile(recoveryPath, ...
                obj.Schema.files.sessionMetadata);
            copyfile(sessionMetadataPath, backupSessionPath, 'f');

            targetLedgerBackup = UMITProjectStore.iBackupLedger(newDataFolder);

            oldLedgerBackup = [];
            removeOldBinding = false;
            [otherPathField, otherBindingField] = ...
                obj.iGetOtherBindingRoleFields(pathField);
            oldFolderStillUsed = ~isempty(oldDataFolder) && ...
                ~isempty(SessionInfo.(otherPathField)) && strcmp( ...
                UMITProjectStore.iNormalizeComparisonPath(oldDataFolder), ...
                UMITProjectStore.iNormalizeComparisonPath( ...
                SessionInfo.(otherPathField)));
            if ~sameFolder && ~isempty(oldDataFolder) && ...
                    ~oldFolderStillUsed && isfolder(oldDataFolder) && ...
                    UMITProjectStore.iHasBinding(oldDataFolder)
                oldLedgerBackup = UMITProjectStore.iBackupLedger(oldDataFolder);
                try
                    oldBinding = UMITProjectStore.readProjectBinding( ...
                        oldDataFolder);
                    removeOldBinding = strcmpi( ...
                        oldBinding.projectUUID, ProjectInfo.projectUUID) && ...
                        strcmpi(oldBinding.sessionUUID, SessionInfo.uuid);
                catch
                    % An invalid ledger cannot take an unbind event;
                    % iRemoveBinding leaves it as it is (with a warning).
                    removeOldBinding = true;
                end
            end

            remainingBinding = [];
            if ~sameFolder && oldFolderStillUsed && ...
                    isfolder(oldDataFolder)
                oldLedgerBackup = UMITProjectStore.iBackupLedger(oldDataFolder);
                remainingBinding = struct();
                remainingBinding.version = ...
                    obj.Schema.projectBinding.version;
                remainingBinding.bindingUUID = ...
                    SessionInfo.(otherBindingField);
                remainingBinding.projectUUID = ...
                    ProjectInfo.projectUUID;
                remainingBinding.subjectUUID = SubjectInfo.uuid;
                remainingBinding.sessionUUID = SessionInfo.uuid;
                remainingBinding.folderRole = otherPathField;
                remainingBinding.createdOn = nowTime;
                remainingBinding.modifiedOn = nowTime;
                UMITProjectStore.iValidateProjectBindingStruct( ...
                    remainingBinding, obj.Schema, errID);
            end

            try
                SessionInfo.(pathField) = newDataFolder;
                SessionInfo.(bindingField) = ProjectBinding.bindingUUID;
                SessionInfo.modifiedOn = datetime('now');
                obj.iSaveSessionInfo(sessionPath, SessionInfo);

                obj.iInstallBinding(newDataFolder, ProjectBinding, ...
                    SubjectInfo.subjectID, SessionInfo.sessionID);

                if removeOldBinding
                    UMITProjectStore.iRemoveBinding(oldDataFolder, 'relocated');
                elseif ~isempty(remainingBinding)
                    % The old folder stays bound for the other role.
                    obj.iInstallBinding(oldDataFolder, remainingBinding, ...
                        SubjectInfo.subjectID, SessionInfo.sessionID);
                end

                obj.iAssertValidAfterMutation(scope);
                obj.IsReadOnly = false;
            catch ME
                copyfile(backupSessionPath, sessionMetadataPath, 'f');
                UMITProjectStore.iRestoreLedger(targetLedgerBackup);
                UMITProjectStore.iRestoreLedger(oldLedgerBackup);
                rethrow(ME)
            end

            obj.iAppendLog(logOperation, ...
                SessionInfo.uuid, newDataFolder);
            clear cleanupRecovery lockCleanup
        end

        function removedSession = iRemoveSessionFromProject( ...
                obj, subjectID, sessionID)
            %IREMOVESESSIONFROMPROJECT Remove ledger metadata, never science data.

            lockCleanup = obj.iAcquireWriteLock('removeSessionFromProject');
            [~, SubjectInfo, subjectPath, SessionInfo, sessionPath, ...
                sessionIndex, sessionRecord] = ...
                obj.iResolveSession(subjectID, sessionID);
            obj.iAssertBindingMutationAllowed(UMITProjectStore.iScope( ...
                'subjects', {SubjectInfo.uuid}, 'sessions', {SessionInfo.uuid}), ...
                sessionRecord.relativePath);

            removedSession = SessionInfo;
            recoveryPath = obj.iCreateRecoveryFolder( ...
                'removeSessionFromProject');
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));
            backupSubjectPath = fullfile(recoveryPath, ...
                obj.Schema.files.subjectMetadata);
            backupSessionPath = fullfile(recoveryPath, 'session');
            copyfile(fullfile(subjectPath, ...
                obj.Schema.files.subjectMetadata), ...
                backupSubjectPath, 'f');

            linkChanges = struct( ...
                'path', {}, 'backupPath', {}, 'invalidatedPath', {});
            roles = {'rawDataFolder', 'processedDataFolder'};
            ProjectInfo = obj.getProjectInfo();

            try
                [ok, message] = movefile( ...
                    sessionPath, backupSessionPath, 'f');
                if ~ok
                    error('Umitoolbox:UMITProjectStore:removeSessionFailed', ...
                        'Could not stage session metadata removal: %s', message);
                end

                SubjectInfo.sessionRegistry(sessionIndex) = [];
                SubjectInfo.modifiedOn = datetime('now');
                obj.iSaveSubjectInfo(subjectPath, SubjectInfo);

                for iRole = 1:numel(roles)
                    dataFolder = SessionInfo.(roles{iRole});
                    if isempty(dataFolder) || ~isfolder(dataFolder)
                        continue
                    end
                    linkPath = fullfile(dataFolder, ...
                        obj.Schema.files.projectBinding);
                    if ~UMITProjectStore.iHasBinding(dataFolder)
                        continue
                    end

                    change = struct();
                    change.path = linkPath;
                    change.backupPath = fullfile(recoveryPath, ...
                        sprintf('link_%d_%s', iRole, obj.Schema.files.projectBinding));
                    change.invalidatedPath = '';
                    copyfile(linkPath, change.backupPath, 'f');

                    removeManagedLink = false;
                    try
                        binding = UMITProjectStore.readProjectBinding( ...
                            dataFolder);
                        removeManagedLink = strcmpi( ...
                            binding.projectUUID, ProjectInfo.projectUUID) && ...
                            strcmpi(binding.sessionUUID, SessionInfo.uuid);
                    catch
                        removeManagedLink = true;
                    end

                    if removeManagedLink
                        UMITProjectStore.iRemoveBinding(dataFolder, 'sessionRemoved');
                        linkChanges(end+1) = change; %#ok<AGROW>
                    end
                end

                obj.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {SubjectInfo.uuid}));
                obj.IsReadOnly = false;
            catch ME
                if ~isfolder(sessionPath) && isfolder(backupSessionPath)
                    movefile(backupSessionPath, sessionPath, 'f');
                end
                copyfile(backupSubjectPath, fullfile(subjectPath, ...
                    obj.Schema.files.subjectMetadata), 'f');
                for iChange = 1:numel(linkChanges)
                    copyfile(linkChanges(iChange).backupPath, ...
                        linkChanges(iChange).path, 'f');
                end
                rethrow(ME)
            end

            obj.iAppendLog('removeSessionFromProject', ...
                SessionInfo.uuid, SessionInfo.processedDataFolder);
            clear cleanupRecovery lockCleanup
        end

        function linkChanges = iPreflightSubjectRemoval( ...
                obj, SubjectInfo, subjectPath, recoveryPath)
            %IPREFLIGHTSUBJECTREMOVAL Verify every external link before mutation.

            linkChanges = struct('path', {}, 'backupPath', {});
            for iSession = 1:numel(SubjectInfo.sessionRegistry)
                sessionRecord = SubjectInfo.sessionRegistry(iSession);
                sessionPath = obj.iResolveRelativePath(sessionRecord.relativePath);
                SessionInfo = obj.iLoadMetadata( ...
                    fullfile(sessionPath, obj.Schema.files.sessionMetadata), ...
                    obj.Schema.metadataVariables.session);
                linkChanges = obj.iPreflightSessionBindings( ...
                    SessionInfo, SubjectInfo, recoveryPath, linkChanges);
            end

            if ~isfolder(subjectPath)
                error('Umitoolbox:UMITProjectStore:removeSubjectFailed', ...
                    'Subject folder is unavailable: %s', subjectPath);
            end
        end

        function [linkChanges, skipped] = iPreflightProjectRemoval( ...
                obj, ProjectInfo, recoveryPath)
            %IPREFLIGHTPROJECTREMOVAL Snapshot the links a deletion removes.
            %
            %   Tolerant (deleteProject, audit item 23): unreadable subject
            %   or session metadata, offline folders and folders whose
            %   ledger does not bind them to this session are skipped and
            %   reported in SKIPPED (struct array: folder, reason).

            linkChanges = struct('path', {}, 'backupPath', {});
            skipped = struct('folder', {}, 'reason', {});
            for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                subjectRecord = ProjectInfo.subjectRegistry(iSubject);
                subjectPath = obj.iResolveRelativePath(subjectRecord.relativePath);
                try
                    SubjectInfo = obj.iLoadMetadata( ...
                        fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                        obj.Schema.metadataVariables.subject);
                catch ME
                    skipped(end+1) = struct('folder', subjectPath, ...
                        'reason', ['unreadable subject metadata: ' ME.message]); %#ok<AGROW>
                    continue
                end
                for iSession = 1:numel(SubjectInfo.sessionRegistry)
                    sessionRecord = SubjectInfo.sessionRegistry(iSession);
                    sessionPath = obj.iResolveRelativePath( ...
                        sessionRecord.relativePath);
                    try
                        SessionInfo = obj.iLoadMetadata( ...
                            fullfile(sessionPath, obj.Schema.files.sessionMetadata), ...
                            obj.Schema.metadataVariables.session);
                    catch ME
                        skipped(end+1) = struct('folder', sessionPath, ...
                            'reason', ['unreadable session metadata: ' ME.message]); %#ok<AGROW>
                        continue
                    end
                    [linkChanges, skipped] = obj.iPreflightSessionBindings( ...
                        SessionInfo, SubjectInfo, recoveryPath, linkChanges, skipped);
                end
            end
        end

        function [linkChanges, skipped] = iPreflightSessionBindings( ...
                obj, SessionInfo, SubjectInfo, recoveryPath, linkChanges, skipped)
            %IPREFLIGHTSESSIONBINDINGS Validate and snapshot managed links.
            %
            %   Strict (called without SKIPPED): an offline or mismatched
            %   folder is an error. Tolerant (SKIPPED passed in): such
            %   folders are recorded in SKIPPED instead.

            tolerant = nargin >= 6;
            roles = {'rawDataFolder', 'processedDataFolder'};
            for iRole = 1:numel(roles)
                [pathField, bindingField] = ...
                    obj.iGetBindingRoleFields(roles{iRole});
                dataFolder = SessionInfo.(pathField);
                if isempty(dataFolder)
                    continue
                end
                if ~isfolder(dataFolder)
                    if tolerant
                        skipped(end+1) = struct('folder', dataFolder, ...
                            'reason', 'folder unavailable (offline or moved)'); %#ok<AGROW>
                        continue
                    end
                    error('Umitoolbox:UMITProjectStore:bindingFolderUnavailable', ...
                        ['Bound folder must be available before project ' ...
                         'removal: %s'], dataFolder);
                end

                try
                    ProjectBinding = UMITProjectStore.readProjectBinding(dataFolder);
                    obj.iAssertBindingMatchesSession( ...
                        ProjectBinding, SessionInfo, SubjectInfo, ...
                        pathField, bindingField, dataFolder);
                catch ME
                    if tolerant
                        skipped(end+1) = struct('folder', dataFolder, ...
                            'reason', ['not bound to this session: ' ME.message]); %#ok<AGROW>
                        continue
                    end
                    rethrow(ME)
                end
                bindingPath = fullfile(dataFolder, ...
                    obj.Schema.files.projectBinding);
                if any(strcmpi({linkChanges.path}, bindingPath))
                    continue
                end
                backupPath = fullfile(recoveryPath, sprintf( ...
                    'link_%d_%s', numel(linkChanges) + 1, obj.Schema.files.projectBinding));
                [ok, message] = copyfile(bindingPath, backupPath, 'f');
                if ~ok
                    error('Umitoolbox:UMITProjectStore:bindingBackupFailed', ...
                        'Could not snapshot project binding: %s', message);
                end
                linkChanges = obj.iAddPreflightedBinding( ...
                    bindingPath, backupPath, linkChanges);
            end
        end

        function linkChanges = iAddPreflightedBinding(~, ...
                bindingPath, backupPath, linkChanges)
            %IADDPREFLIGHTEDBINDING Add one unique external binding snapshot.

            if any(strcmpi({linkChanges.path}, bindingPath))
                return
            end
            linkChanges(end+1) = struct( ...
                'path', bindingPath, 'backupPath', backupPath);
        end

        function iDeletePreflightedBindings(~, linkChanges)
            %IDELETEPREFLIGHTEDBINDINGS Unbind only folders validated in preflight.
            %
            %   Each folder's project ledger gets an unbind event (reason
            %   'removedFromProject'); the ledger and its history stay.

            for iChange = 1:numel(linkChanges)
                dataFolder = fileparts(linkChanges(iChange).path);
                if ~UMITProjectLedger.isBound(dataFolder)
                    error('Umitoolbox:UMITProjectStore:bindingChanged', ...
                        ['A project binding changed after deletion preflight: ' ...
                         '%s'], linkChanges(iChange).path);
                end
                UMITProjectStore.iRemoveBinding(dataFolder, 'removedFromProject');
            end
        end

        function iRestorePreflightedBindings(~, linkChanges)
            %IRESTOREPREFLIGHTEDBINDINGS Restore binding snapshots on failure.

            for iChange = 1:numel(linkChanges)
                if isfile(linkChanges(iChange).backupPath)
                    try
                        copyfile(linkChanges(iChange).backupPath, ...
                            linkChanges(iChange).path, 'f');
                    catch
                        % Preserve the original mutation exception. Recovery
                        % artifacts remain available for a manual repair.
                    end
                end
            end
        end

        function count = iCountProjectSessions(obj, ProjectInfo)
            %ICOUNTPROJECTSESSIONS Number of registered sessions (readable subjects).

            count = 0;
            for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                subjectPath = obj.iResolveRelativePath( ...
                    ProjectInfo.subjectRegistry(iSubject).relativePath);
                try
                    SubjectInfo = obj.iLoadMetadata( ...
                        fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                        obj.Schema.metadataVariables.subject);
                    count = count + numel(SubjectInfo.sessionRegistry);
                catch
                    % Unreadable subject: counted as no sessions.
                end
            end
        end

        function ProjectBinding = iUnbindSessionDataFolder( ...
                obj, subjectID, sessionID, folderRole)
            %IUNBINDSESSIONDATAFOLDER Remove reciprocal binding transactionally.

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock( ...
                ['unbind_', folderRole]);

            [~, SubjectInfo, ~, SessionInfo, sessionPath] = ...
                obj.iResolveSession(subjectID, sessionID);
            scope = UMITProjectStore.iScope('subjects', {SubjectInfo.uuid}, ...
                'sessions', {SessionInfo.uuid});
            obj.iAssertHealthyForMutation(scope);
            [pathField, bindingField] = ...
                obj.iGetBindingRoleFields(folderRole);

            if isempty(SessionInfo.(pathField)) && ...
                    isempty(SessionInfo.(bindingField))
                ProjectBinding = [];
                clear lockCleanup
                return
            end

            if isempty(SessionInfo.(pathField)) || ...
                    isempty(SessionInfo.(bindingField))
                error('Umitoolbox:UMITProjectStore:incompleteBinding', ...
                    ['Session binding metadata is incomplete for %s.'], ...
                    pathField);
            end

            dataFolder = SessionInfo.(pathField);
            if ~isfolder(dataFolder)
                error('Umitoolbox:UMITProjectStore:bindingFolderUnavailable', ...
                    ['Bound folder must be available before it can be ' ...
                     'detached: %s'], dataFolder);
            end

            ProjectBinding = ...
                UMITProjectStore.readProjectBinding(dataFolder);
            obj.iAssertBindingMatchesSession( ...
                ProjectBinding, SessionInfo, SubjectInfo, ...
                pathField, bindingField, dataFolder);

            bindingPath = fullfile(dataFolder, ...
                obj.Schema.files.projectBinding);
            [otherPathField, otherBindingField] = ...
                obj.iGetOtherBindingRoleFields(pathField);
            folderStillUsed = ~isempty(SessionInfo.(otherPathField)) && ...
                strcmp( ...
                UMITProjectStore.iNormalizeComparisonPath(dataFolder), ...
                UMITProjectStore.iNormalizeComparisonPath( ...
                SessionInfo.(otherPathField)));
            recoveryPath = ...
                obj.iCreateRecoveryFolder(['unbind_', folderRole]);
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));

            sessionMetadataPath = fullfile(sessionPath, ...
                obj.Schema.files.sessionMetadata);
            backupSessionPath = fullfile(recoveryPath, ...
                obj.Schema.files.sessionMetadata);
            backupBindingPath = fullfile(recoveryPath, ...
                obj.Schema.files.projectBinding);
            copyfile(sessionMetadataPath, backupSessionPath, 'f');
            copyfile(bindingPath, backupBindingPath, 'f');

            replacementBinding = [];
            if folderStillUsed && ...
                    strcmp(ProjectBinding.folderRole, pathField)
                replacementBinding = ProjectBinding;
                replacementBinding.bindingUUID = ...
                    SessionInfo.(otherBindingField);
                replacementBinding.folderRole = otherPathField;
                replacementBinding.modifiedOn = datetime('now');
            end

            try
                SessionInfo.(pathField) = '';
                SessionInfo.(bindingField) = '';
                SessionInfo.modifiedOn = datetime('now');
                obj.iSaveSessionInfo(sessionPath, SessionInfo);

                if ~folderStillUsed
                    UMITProjectStore.iRemoveBinding(dataFolder, 'folderUnbound');
                elseif ~isempty(replacementBinding)
                    % The folder stays bound for the other role.
                    obj.iInstallBinding(dataFolder, replacementBinding, ...
                        SubjectInfo.subjectID, SessionInfo.sessionID);
                end
                obj.iAssertValidAfterMutation(scope);

            catch ME
                copyfile(backupSessionPath, ...
                    sessionMetadataPath, 'f');
                if (~folderStillUsed || ...
                        ~isempty(replacementBinding)) && ...
                        isfile(backupBindingPath)
                    copyfile(backupBindingPath, bindingPath, 'f');
                end
                rethrow(ME)
            end

            obj.iAppendLog('unbindDataFolder', ...
                ProjectBinding.bindingUUID, pathField);
            clear cleanupRecovery lockCleanup
        end

        function [pathField, bindingField] = ...
                iGetBindingRoleFields(~, folderRole)
            %IGETBINDINGROLEFIELDS Resolve folder and reciprocal UUID fields.

            folderRole = char(string(folderRole));
            switch folderRole
                case 'rawDataFolder'
                    pathField = 'rawDataFolder';
                    bindingField = 'rawDataBindingUUID';
                case 'processedDataFolder'
                    pathField = 'processedDataFolder';
                    bindingField = 'processedDataBindingUUID';
                otherwise
                    error('Umitoolbox:UMITProjectStore:invalidBindingRole', ...
                        'Unsupported binding folder role: %s', folderRole);
            end
        end

        function iInstallBinding(obj, dataFolder, ProjectBinding, subjectID, sessionID)
            %IINSTALLBINDING Record ProjectBinding in the folder's project ledger.
            %
            %   An unbound folder gets a bind event; a bound folder (an
            %   orphan being replaced, a session moved to another subject or
            %   project) gets a rebind event recording the previous binding.
            %   Both carry a frozen human snapshot of the project metadata
            %   and the subject and session IDs. Callers back the ledger up
            %   first (iBackupLedger) and restore it on failure.

            identity = struct( ...
                'bindingUUID', lower(ProjectBinding.bindingUUID), ...
                'projectUUID', lower(ProjectBinding.projectUUID), ...
                'subjectUUID', lower(ProjectBinding.subjectUUID), ...
                'sessionUUID', lower(ProjectBinding.sessionUUID), ...
                'folderRole', ProjectBinding.folderRole);
            snapshot = UMITProjectLedger.makeSnapshot(obj.getProjectInfo(), ...
                subjectID, sessionID);
            origin = 'UMITProjectStore';
            if UMITProjectLedger.isBound(dataFolder)
                current = UMITProjectLedger.read(dataFolder).binding;
                if isequal(rmfield(current, {'boundOn', 'modifiedOn'}), ...
                        orderfields(identity, UMITProjectLedger.BindingFields))
                    return
                end
                UMITProjectLedger.rebind(dataFolder, identity, snapshot, 'Origin', origin);
            else
                UMITProjectLedger.bind(dataFolder, identity, snapshot, 'Origin', origin);
            end
        end

        function iAssertBindingMatchesSession(obj, ProjectBinding, ...
                SessionInfo, SubjectInfo, pathField, bindingField, ...
                currentFolder)
            %IASSERTBINDINGMATCHESSESSION Verify reciprocal binding identity.

            ProjectInfo = obj.getProjectInfo();

            if ~strcmpi(ProjectBinding.projectUUID, ...
                    ProjectInfo.projectUUID)
                error('Umitoolbox:UMITProjectStore:bindingProjectMismatch', ...
                    'Binding project UUID does not match this project.');
            end
            if ~strcmpi(ProjectBinding.subjectUUID, SubjectInfo.uuid) || ...
                    ~strcmpi(SessionInfo.subjectUUID, SubjectInfo.uuid)
                error('Umitoolbox:UMITProjectStore:bindingOwnerMismatch', ...
                    'Binding subject UUID does not match the session owner.');
            end
            if ~strcmpi(ProjectBinding.sessionUUID, SessionInfo.uuid)
                error('Umitoolbox:UMITProjectStore:bindingSessionMismatch', ...
                    'Binding session UUID does not match SessionInfo.uuid.');
            end
            roleMatches = strcmp(ProjectBinding.folderRole, pathField);
            bindingMatches = strcmpi(ProjectBinding.bindingUUID, ...
                SessionInfo.(bindingField));
            if ~roleMatches
                [otherPathField, otherBindingField] = ...
                    obj.iGetOtherBindingRoleFields(pathField);
                roleMatches = strcmp(ProjectBinding.folderRole, ...
                    otherPathField) && ...
                    ~isempty(SessionInfo.(otherPathField)) && ...
                    strcmpi(ProjectBinding.bindingUUID, ...
                    SessionInfo.(otherBindingField)) && ...
                    strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath( ...
                    SessionInfo.(otherPathField)), ...
                    UMITProjectStore.iNormalizeComparisonPath( ...
                    currentFolder));
            end
            if ~roleMatches
                error('Umitoolbox:UMITProjectStore:bindingRoleMismatch', ...
                    ['Binding folderRole does not match %s or a co-located ' ...
                     'session data role.'], pathField);
            end
            if ~bindingMatches
                error('Umitoolbox:UMITProjectStore:bindingUUIDMismatch', ...
                    'Binding UUID does not match SessionInfo.%s.', ...
                    bindingField);
            end

            registeredFolder = SessionInfo.(pathField);
            if isempty(registeredFolder) || ~strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath( ...
                    registeredFolder), ...
                    UMITProjectStore.iNormalizeComparisonPath( ...
                    currentFolder))
                error('Umitoolbox:UMITProjectStore:bindingPathMismatch', ...
                    ['Current folder does not match SessionInfo.%s. The ' ...
                     'bound folder may have been moved or copied.'], ...
                    pathField);
            end
        end

        function matches = iFindSessionPathMatches(obj, dataFolder)
            %IFINDSESSIONPATHMATCHES Find sessions storing an exact data path.
            %   Unreadable subject or session metadata is skipped (per-session
            %   validity); the folder's own project ledger still guards
            %   against binding a folder twice.

            comparisonPath = ...
                UMITProjectStore.iNormalizeComparisonPath(dataFolder);
            ProjectInfo = obj.getProjectInfo();
            template = struct( ...
                'subjectID', '', ...
                'subjectUUID', '', ...
                'sessionID', '', ...
                'sessionUUID', '', ...
                'matchedField', '');
            matches = repmat(template, 0, 1);

            for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                subjectRecord = ProjectInfo.subjectRegistry(iSubject);
                subjectPath = obj.iResolveRelativePath( ...
                    subjectRecord.relativePath);
                try
                    SubjectInfo = obj.iLoadMetadata( ...
                        fullfile(subjectPath, ...
                        obj.Schema.files.subjectMetadata), ...
                        obj.Schema.metadataVariables.subject);
                catch
                    continue
                end

                for iSession = 1:numel(SubjectInfo.sessionRegistry)
                    sessionRecord = ...
                        SubjectInfo.sessionRegistry(iSession);
                    sessionPath = obj.iResolveRelativePath( ...
                        sessionRecord.relativePath);
                    try
                        SessionInfo = obj.iLoadMetadata( ...
                            fullfile(sessionPath, ...
                            obj.Schema.files.sessionMetadata), ...
                            obj.Schema.metadataVariables.session);
                    catch
                        continue
                    end

                    roles = {'rawDataFolder', ...
                        'processedDataFolder'};
                    for iRole = 1:numel(roles)
                        pathField = roles{iRole};
                        if isempty(SessionInfo.(pathField))
                            continue
                        end
                        if strcmp( ...
                                UMITProjectStore.iNormalizeComparisonPath( ...
                                SessionInfo.(pathField)), comparisonPath)
                            match = template;
                            match.subjectID = SubjectInfo.subjectID;
                            match.subjectUUID = SubjectInfo.uuid;
                            match.sessionID = SessionInfo.sessionID;
                            match.sessionUUID = SessionInfo.uuid;
                            match.matchedField = pathField;
                            matches(end+1) = match; %#ok<AGROW>
                        end
                    end
                end
            end
        end

        function location = iLocateResource(obj, resourceUUID)
            %ILOCATERESOURCE Find one uniquely registered resource UUID.

            resourceUUID = UMITProjectStore.iNormalizeUUIDInput(resourceUUID);
            locations = obj.iFindResourceLocations(resourceUUID);

            if isempty(locations)
                error('Umitoolbox:UMITProjectStore:resourceNotFound', ...
                    'Resource UUID was not found: %s', resourceUUID);
            end

            if numel(locations) > 1
                error('Umitoolbox:UMITProjectStore:duplicateResourceUUID', ...
                    ['Resource UUID is registered more than once and cannot ' ...
                    'be resolved safely: %s'], resourceUUID);
            end

            location = locations(1);
        end

        function iSaveLocatedMetadata(~, location, metadata)
            %ISAVELOCATEDMETADATA Save metadata described by a resource location.

            saveMatAtomic(location.metadataPath, ...
                location.variableName, metadata);
        end

        function iSaveOwnerMetadata(~, owner, metadata)
            %ISAVEOWNERMETADATA Save resolved owner metadata atomically.

            saveMatAtomic(owner.metadataPath, owner.variableName, metadata);
        end

        function iSaveProjectInfo(obj, ProjectInfo)
            %ISAVEPROJECTINFO Save project metadata atomically.

            saveMatAtomic( ...
                fullfile(obj.ProjectRoot, obj.Schema.files.projectMetadata), ...
                obj.Schema.metadataVariables.project, ProjectInfo);
        end

        function iSaveSubjectInfo(obj, subjectPath, SubjectInfo)
            %ISAVESUBJECTINFO Save subject metadata atomically.

            saveMatAtomic( ...
                fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                obj.Schema.metadataVariables.subject, SubjectInfo);
        end

        function iSaveSessionInfo(obj, sessionPath, SessionInfo)
            %ISAVESESSIONINFO Save session metadata atomically.

            saveMatAtomic( ...
                fullfile(sessionPath, obj.Schema.files.sessionMetadata), ...
                obj.Schema.metadataVariables.session, SessionInfo);
        end

        function value = iLoadMetadata(~, filePath, variableName)
            %ILOADMETADATA Load one required metadata variable from a MAT file.

            if ~isfile(filePath)
                error('Umitoolbox:UMITProjectStore:metadataMissing', ...
                    'Metadata file is missing: %s', filePath);
            end

            loaded = load(filePath, variableName, '-mat');
            if ~isfield(loaded, variableName)
                error('Umitoolbox:UMITProjectStore:metadataVariableMissing', ...
                    'Variable "%s" is missing from "%s".', ...
                    variableName, filePath);
            end

            value = loaded.(variableName);
            if ~isstruct(value) || ~isscalar(value)
                error('Umitoolbox:UMITProjectStore:invalidMetadata', ...
                    'Metadata variable "%s" must be a scalar struct.', ...
                    variableName);
            end

            switch variableName
                case 'ProjectInfo'
                    entityType = 'project';
                case 'SubjectInfo'
                    entityType = 'subject';
                case 'SessionInfo'
                    entityType = 'session';
                otherwise
                    entityType = '';
            end
            if ~isempty(entityType)
                % Older project folders remain readable: newly optional
                % fields are supplied in memory and persist on the next
                % supported metadata update.
                value = UMITProjectStore.iAddMetadataDefaults(value, entityType);
            end
        end

        function cleanupObj = iAcquireWriteLock(obj, operation)
            %IACQUIREWRITELOCK Acquire the exclusive project mutation lock.
            %
            % The lock is represented by a directory. Directory creation is an
            % atomic filesystem operation, preventing two writers from both
            % acquiring the project lock during the same check/create window.

            lockFolder = fullfile(obj.ProjectRoot, ...
                obj.Schema.folders.internal, obj.Schema.folders.lock);
            lockMetadataPath = fullfile(lockFolder, ...
                obj.Schema.files.lockMetadata);

            lockDirectory = java.io.File(lockFolder);
            if ~lockDirectory.mkdir()
                if ~isfolder(lockFolder)
                    error('Umitoolbox:UMITProjectStore:lockCreationFailed', ...
                        'Could not create project lock directory: %s', ...
                        lockFolder);
                end

                details = '';
                try
                    existing = obj.getLockInfo();
                    if isfield(existing, 'operation')
                        details = sprintf(' Existing operation: %s.', ...
                            existing.operation);
                    end
                catch
                end

                error('Umitoolbox:UMITProjectStore:projectLocked', ...
                    'Project is locked by another operation.%s', details);
            end

            try
                LockInfo = struct();
                LockInfo.operation = operation;
                LockInfo.createdOn = datetime('now');
                LockInfo.processID = feature('getpid');
                LockInfo.userName = UMITProjectStore.iGetEnvironmentValue( ...
                    {'USERNAME', 'USER'}, 'unknown');
                LockInfo.hostName = UMITProjectStore.iGetEnvironmentValue( ...
                    {'COMPUTERNAME', 'HOSTNAME'}, 'unknown');

                saveMatAtomic(lockMetadataPath, 'LockInfo', LockInfo);
            catch ME
                UMITProjectStore.iRemoveFolderIfPresent(lockFolder);
                rethrow(ME)
            end

            cleanupObj = onCleanup(@() obj.iReleaseWriteLock(lockFolder));
        end

        function iReleaseWriteLock(~, lockFolder)
            %IRELEASEWRITELOCK Remove the project mutation lock directory.

            try
                UMITProjectStore.iRemoveFolderIfPresent(lockFolder);
            catch ME
                warning('Umitoolbox:UMITProjectStore:lockReleaseFailed', ...
                    ['Project operation completed, but the lock directory ' ...
                    'could not be removed: %s'], ME.message);
            end
        end

        function iAssertWritable(obj)
            %IASSERTWRITABLE Reject writes when the project is read-only.

            if obj.IsReadOnly
                error('Umitoolbox:UMITProjectStore:readOnlyProject', ...
                    ['Project is open in read-only mode because validation ' ...
                    'failed. Repair the project before modifying it.']);
            end
        end

        function iAssertHealthyForMutation(obj, scope)
            %IASSERTHEALTHYFORMUTATION Validate what a mutation touches first.
            %
            %   obj.iAssertHealthyForMutation(scope)
            %
            %   Per-session validity (brief section 1.11): a quick, scoped
            %   validation (UMITProjectStore.iScope) of the project level and
            %   the subjects/sessions the mutation touches; no checksum of
            %   every resource. Any error found blocks the mutation. Errors of
            %   other subjects or sessions are not looked at and never block
            %   it. Only store-level errors make the store read-only.

            report = obj.validate('Mode', 'quick', 'Scope', scope);
            obj.iThrowBlockingErrors(report.errors);
        end

        function iThrowBlockingErrors(obj, blockingErrors)
            %ITHROWBLOCKINGERRORS Refuse a mutation blocked by validation errors.
            %
            %   Store-level errors: Umitoolbox:UMITProjectStore:invalidProject
            %   (and the store becomes read-only). Errors of the touched
            %   subject or session only: ...:invalidTarget (the rest of the
            %   project stays writable).

            if isempty(blockingErrors)
                return
            end
            if any([blockingErrors.storeLevel])
                obj.IsReadOnly = true;
                error('Umitoolbox:UMITProjectStore:invalidProject', ...
                    'Project validation failed: %s', ...
                    UMITProjectStore.iJoinIssueMessages(blockingErrors));
            end
            error('Umitoolbox:UMITProjectStore:invalidTarget', ...
                ['The subject or session to modify has validation errors ' ...
                 '(the rest of the project is unaffected): %s'], ...
                UMITProjectStore.iJoinIssueMessages(blockingErrors));
        end

        function [pathField, bindingField] = ...
                iGetOtherBindingRoleFields(obj, folderRole)
            %IGETOTHERBINDINGROLEFIELDS Return the complementary data role.

            folderRole = char(string(folderRole));
            if strcmp(folderRole, 'rawDataFolder')
                otherRole = 'processedDataFolder';
            elseif strcmp(folderRole, 'processedDataFolder')
                otherRole = 'rawDataFolder';
            else
                error('Umitoolbox:UMITProjectStore:invalidBindingRole', ...
                    'Unsupported binding folder role: %s', folderRole);
            end
            [pathField, bindingField] = ...
                obj.iGetBindingRoleFields(otherRole);
        end

        function [tf, ProjectBinding] = ...
                iGetCoLocatedSessionBinding(obj, dataFolder, ...
                SessionInfo, SubjectInfo, requestedPathField)
            %IGETCOLOCATEDSESSIONBINDING Find the other role's shared link.

            tf = false;
            ProjectBinding = [];
            [otherPathField, otherBindingField] = ...
                obj.iGetOtherBindingRoleFields(requestedPathField);
            if isempty(SessionInfo.(otherPathField)) || ...
                    isempty(SessionInfo.(otherBindingField)) || ...
                    ~strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath( ...
                    SessionInfo.(otherPathField)), ...
                    UMITProjectStore.iNormalizeComparisonPath(dataFolder))
                return
            end

            if ~UMITProjectLedger.isBound(dataFolder)
                return
            end

            try
                ProjectBinding = ...
                    UMITProjectStore.readProjectBinding(dataFolder);
                obj.iAssertBindingMatchesSession( ...
                    ProjectBinding, SessionInfo, SubjectInfo, ...
                    otherPathField, otherBindingField, dataFolder);
                tf = true;
            catch
                ProjectBinding = [];
            end
        end

        function ProjectBinding = iAdoptCoLocatedSessionBinding( ...
                obj, SessionInfo, sessionPath, pathField, bindingField, ...
                dataFolder, ProjectBinding, logOperation)
            %IADOPTCOLOCATEDSESSIONBINDING Share the other role's link atomically.

            oldDataFolder = SessionInfo.(pathField);
            if ~isempty(oldDataFolder) && strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath(oldDataFolder), ...
                    UMITProjectStore.iNormalizeComparisonPath(dataFolder)) && ...
                    strcmp(SessionInfo.(bindingField), ...
                    ProjectBinding.bindingUUID)
                return
            end

            recoveryPath = obj.iCreateRecoveryFolder( ...
                ['share_', pathField]);
            cleanupRecovery = onCleanup(@() ...
                UMITProjectStore.iRemoveFolderIfPresent(recoveryPath));
            sessionMetadataPath = fullfile(sessionPath, ...
                obj.Schema.files.sessionMetadata);
            backupSessionPath = fullfile(recoveryPath, ...
                obj.Schema.files.sessionMetadata);
            copyfile(sessionMetadataPath, backupSessionPath, 'f');

            oldBindingPath = '';
            backupOldPath = '';
            removeOldBinding = false;
            if ~isempty(oldDataFolder) && ...
                    ~strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath(oldDataFolder), ...
                    UMITProjectStore.iNormalizeComparisonPath(dataFolder)) && ...
                    isfolder(oldDataFolder)
                oldBindingPath = fullfile(oldDataFolder, ...
                    obj.Schema.files.projectBinding);
                if UMITProjectStore.iHasBinding(oldDataFolder)
                    backupOldPath = fullfile(recoveryPath, ...
                        ['old_', obj.Schema.files.projectBinding]);
                    copyfile(oldBindingPath, backupOldPath, 'f');
                    try
                        oldBinding = UMITProjectStore.readProjectBinding( ...
                            oldDataFolder);
                        ProjectInfo = obj.getProjectInfo();
                        removeOldBinding = strcmpi( ...
                            oldBinding.projectUUID, ...
                            ProjectInfo.projectUUID) && strcmpi( ...
                            oldBinding.sessionUUID, SessionInfo.uuid);
                    catch
                        removeOldBinding = true;
                    end
                end
            end

            try
                SessionInfo.(pathField) = dataFolder;
                SessionInfo.(bindingField) = ProjectBinding.bindingUUID;
                SessionInfo.modifiedOn = datetime('now');
                obj.iSaveSessionInfo(sessionPath, SessionInfo);
                if removeOldBinding
                    UMITProjectStore.iRemoveBinding(oldDataFolder, 'relocated');
                end
                obj.iAssertValidAfterMutation(UMITProjectStore.iScope( ...
                    'subjects', {SessionInfo.subjectUUID}, 'sessions', {SessionInfo.uuid}));
                obj.IsReadOnly = false;
            catch ME
                copyfile(backupSessionPath, sessionMetadataPath, 'f');
                if removeOldBinding && isfile(backupOldPath)
                    copyfile(backupOldPath, oldBindingPath, 'f');
                end
                rethrow(ME)
            end

            obj.iAppendLog(logOperation, ...
                SessionInfo.uuid, dataFolder);
            clear cleanupRecovery
        end

        function iAssertBindingMutationAllowed(obj, scope, sessionRelativePath)
            %IASSERTBINDINGMUTATIONALLOWED Permit narrowly scoped repair writes.
            %
            % Relocation, repair, and explicit session removal may proceed
            % when the only validation errors in their scope are
            % binding-consistency or dataset errors on the selected session.

            report = obj.validate('Mode', 'quick', 'Scope', scope);
            if report.isValid
                return
            end

            allowedCodes = { ...
                'missing_project_binding', ...
                'invalid_project_binding', ...
                'incomplete_data_folder_binding', ...
                'not_a_dataset', ...
                'invalid_acq_infos', ...
                'missing_raw_data_files'};
            errors = report.errors;
            isAllowed = ismember({errors.code}, allowedCodes) & ...
                strcmp({errors.relativePath}, sessionRelativePath);
            obj.iThrowBlockingErrors(errors(~isAllowed));
        end

        function iAssertValidAfterMutation(obj, scope)
            %IASSERTVALIDAFTERMUTATION Validate what a mutation touched.
            %   Quick, scoped validation (see iAssertHealthyForMutation).

            report = obj.validate('Mode', 'quick', 'Scope', scope);
            if ~report.isValid
                error('Umitoolbox:UMITProjectStore:postWriteValidationFailed', ...
                    'Project failed post-write validation: %s', ...
                    UMITProjectStore.iJoinIssueMessages(report.errors));
            end
        end

        function id = iNormalizeManagedID(obj, idIn, idLabel, allowRewrite)
            %INORMALIZEMANAGEDID Normalize one filesystem-backed MATLAB ID.
            %
            % For lookups (ALLOWREWRITE true, the default), IDs that fail the
            % managed-ID grammar are converted with matlab.lang.makeValidName
            % before the managed-ID checks, so "Mouse 104" finds "Mouse104".
            % For new and renamed IDs (ALLOWREWRITE false) a rewrite is never
            % applied silently (audit item 24): the ID is refused with
            % Umitoolbox:UMITProjectStore:invalidID and the message names
            % the valid form to use instead. The store's reserved-name and
            % case-insensitive uniqueness protections remain in force.

            if nargin < 4
                allowRewrite = true;
            end
            if ~(ischar(idIn) || (isstring(idIn) && isscalar(idIn)))
                error('Umitoolbox:UMITProjectStore:invalidID', ...
                    '"%s" must be a character vector or string scalar.', idLabel);
            end

            id = strtrim(char(string(idIn)));
            rules = obj.Schema.namingRules;

            if isempty(id)
                error('Umitoolbox:UMITProjectStore:invalidID', ...
                    '"%s" must contain 1 to %d characters.', ...
                    idLabel, rules.maxIDLength);
            end

            if isempty(regexp(id, rules.idPattern, 'once'))
                requested = id;
                id = matlab.lang.makeValidName(id, ...
                    'ReplacementStyle', 'underscore');
                if ~allowRewrite
                    error('Umitoolbox:UMITProjectStore:invalidID', ...
                        ['%s "%s" is not a valid name (letters, digits, "_" and "-", ' ...
                         'starting with a letter or digit). Use "%s" instead.'], ...
                        idLabel, requested, id);
                end
            end
            if strlength(string(id)) > rules.maxIDLength
                error('Umitoolbox:UMITProjectStore:invalidID', ...
                    '"%s" must normalize to at most %d characters.', ...
                    idLabel, rules.maxIDLength);
            end

            if isempty(regexp(id, rules.idPattern, 'once'))
                error('Umitoolbox:UMITProjectStore:invalidID', ...
                    '"%s" could not be converted to a valid managed name.', ...
                    idLabel);
            end

            if any(strcmpi(id, rules.reservedNames))
                error('Umitoolbox:UMITProjectStore:invalidID', ...
                    '"%s" uses a reserved filesystem name: %s', idLabel, id);
            end

            for iPrefix = 1:numel(rules.reservedPrefixes)
                if startsWith(id, rules.reservedPrefixes{iPrefix}, ...
                        'IgnoreCase', true)
                    error('Umitoolbox:UMITProjectStore:invalidID', ...
                        '"%s" uses reserved prefix "%s".', ...
                        idLabel, rules.reservedPrefixes{iPrefix});
                end
            end
        end

        function iAssertUniqueRegistryID(~, registry, newID, typeName, ignoredUUID)
            %IASSERTUNIQUEREGISTRYID Enforce case-insensitive ID uniqueness.
            %
            %   Excludes the registry entry identified by ignoredUUID when validating
            %   an ID change. This permits case-only renames while continuing to reject
            %   collisions with other registered items.

            if nargin < 5
                ignoredUUID = '';
            end

            for iRecord = 1:numel(registry)

                % Ignore the item currently being renamed.
                if ~isempty(ignoredUUID) && ...
                        strcmpi(registry(iRecord).uuid, ignoredUUID)
                    continue
                end

                if strcmpi(registry(iRecord).id, newID)
                    error('Umitoolbox:UMITProjectStore:duplicateID', ...
                        '%s ID already exists: %s', typeName, newID);
                end
            end
        end

        function metadata = iApplyEditableUpdates(~, metadata, updates, allowedFields, errID)
            %IAPPLYEDITABLEUPDATES Apply only explicitly allowed metadata fields.

            updateFields = fieldnames(updates);
            invalidFields = setdiff(updateFields, allowedFields);
            if ~isempty(invalidFields)
                error(errID, 'Unsupported editable field(s): %s', ...
                    strjoin(invalidFields, ', '));
            end

            for iField = 1:numel(updateFields)
                fieldName = updateFields{iField};
                metadata.(fieldName) = UMITProjectStore.iNormalizeMetadataValue( ...
                    fieldName, updates.(fieldName), errID);
            end
        end

        function iAssertExternalDataPath(obj, folderPath, fieldName)
            %IASSERTEXTERNALDATAPATH Reject data paths inside the project root.

            if isempty(folderPath)
                return
            end

            folderPath = UMITProjectStore.iAbsolutePath(folderPath);
            projectRoot = obj.ProjectRoot;
            comparePath = folderPath;
            compareRoot = projectRoot;
            if ispc || ismac
                comparePath = lower(comparePath);
                compareRoot = lower(compareRoot);
            end

            if strcmp(comparePath, compareRoot) || ...
                    startsWith(comparePath, [compareRoot, filesep])
                error('Umitoolbox:UMITProjectStore:internalDataPath', ...
                    ['Field "%s" must reference a folder outside the ' ...
                    'centralized project metadata folder.'], fieldName);
            end
        end

        function iAssertValidSaveFolderDataset(~, saveFolder)
            %IASSERTVALIDSAVEFOLDERDATASET Require imported imaging data.
            %
            %   Format-neutral (audit section 12.1 item 5): the folder must
            %   hold a .dat file with a UMIT header, or a valid AcqInfos.mat
            %   (isUmitDataset). AcqInfos.mat alone is no longer required.
            %   Errors: Umitoolbox:UMITProjectStore:notADataset (no marker),
            %   ...:invalidAcqInfos (only an unreadable AcqInfos.mat).

            [isDataset, reason] = isUmitDataset(saveFolder);
            if isDataset
                return
            end
            if isfile(fullfile(saveFolder, 'AcqInfos.mat'))
                error('Umitoolbox:UMITProjectStore:invalidAcqInfos', ...
                    'SaveFolder is not a UMIT dataset: %s.', reason);
            end
            error('Umitoolbox:UMITProjectStore:notADataset', ...
                ['SaveFolder is not a UMIT dataset: %s. Import the raw data ' ...
                 'into it first.'], reason);
        end

        function rawFolder = iDataParamsRawFolder(~, SessionInfo)
            %IDATAPARAMSRAWFOLDER The SaveFolder's DataParams.folders.RawFolder.
            %   '' when the SaveFolder, its DataParams or the field is
            %   unavailable or 'Missing'. Read-only.

            rawFolder = '';
            saveFolder = SessionInfo.processedDataFolder;
            if isempty(saveFolder) || ~isfile(fullfile(saveFolder, 'DataParams.mat'))
                return
            end
            try
                DataParams = loadDataParams(saveFolder);
                value = char(string(DataParams.folders.RawFolder));
            catch
                return
            end
            if ~isempty(value) && ~strcmpi(value, 'Missing')
                rawFolder = value;
            end
        end

        function iAssertRawFolderAuthority(obj, SessionInfo, rawFolder)
            %IASSERTRAWFOLDERAUTHORITY Raw binding must follow DataParams (item 7).

            authority = obj.iDataParamsRawFolder(SessionInfo);
            if isempty(authority)
                return
            end
            if ~isfolder(authority) || ~strcmp( ...
                    UMITProjectStore.iNormalizeComparisonPath(authority), ...
                    UMITProjectStore.iNormalizeComparisonPath(rawFolder))
                error('Umitoolbox:UMITProjectStore:rawFolderMismatch', ...
                    ['The SaveFolder''s DataParams.folders.RawFolder is "%s", the ' ...
                     'authority for its RawFolder; "%s" was given. Update DataParams ' ...
                     'first if the raw data moved.'], authority, rawFolder);
            end
        end

        function iAssertValidRawDataFolder(~, rawDataFolder)
            %IASSERTVALIDRAWDATAFOLDER Require recognizable raw image data.

            entries = dir(rawDataFolder);
            entries = entries(~[entries.isdir]);
            fileNames = {entries.name};
            hasRawData = any(~cellfun('isempty', regexpi( ...
                fileNames, '\.(bin|tiff?)$', 'once')));
            if ~hasRawData
                error( ...
                    'Umitoolbox:UMITProjectStore:missingRawDataFiles', ...
                    ['RawFolder must contain at least one .bin, .tif, or ' ...
                     '.tiff file: %s'], rawDataFolder);
            end
        end

        function iAssertUpdateStruct(~, updates, errID)
            %IASSERTUPDATESTRUCT Validate an update payload.

            if ~isstruct(updates) || ~isscalar(updates)
                error(errID, '"updates" must be a scalar struct.');
            end
        end

        function path = iCreateTransactionFolder(obj, operation)
            %ICREATETRANSACTIONFOLDER Create a unique internal staging folder.

            name = sprintf('%s_%s', operation, UMITProjectStore.iGenerateUUID());
            path = fullfile(obj.ProjectRoot, obj.Schema.folders.internal, ...
                obj.Schema.folders.transactions, name);
            mkdir(path);
        end

        function path = iCreateRecoveryFolder(obj, operation)
            %ICREATERECOVERYFOLDER Create a unique recovery snapshot folder.

            name = sprintf('%s_%s', operation, UMITProjectStore.iGenerateUUID());
            path = fullfile(obj.ProjectRoot, obj.Schema.folders.internal, ...
                obj.Schema.folders.recovery, name);
            mkdir(path);
        end

        function iCreateSubjectFolders(obj, subjectPath)
            %ICREATESUBJECTFOLDERS Create canonical subject folders.

            mkdir(subjectPath);
            base = fullfile(subjectPath, ...
                obj.Schema.folders.references, obj.Schema.folders.image);
            mkdir(fullfile(base, obj.Schema.folders.active));
            mkdir(fullfile(base, obj.Schema.folders.archive));
            mkdir(fullfile(subjectPath, obj.Schema.folders.sessions));
        end

        function iCreateSessionFolders(~, sessionPath)
            %ICREATESESSIONFOLDERS Create the canonical session folder.
            %   Sessions own no resource folders since the registration
            %   transforms were retired (section 4 Q5/Q9).

            mkdir(sessionPath);
        end

        function rel = iBuildResourceRelativePath(obj, ownerBaseRel, resourceType, stateFolder, fileName)
            %IBUILDRESOURCERELATIVEPATH Build a canonical resource path.

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            parts = [{ownerBaseRel}, resourceDef.relativeFolderParts, ...
                {stateFolder, fileName}];
            rel = UMITProjectStore.iJoinRelative(parts{:});
        end

        function path = iResolveRelativePath(obj, relativePath)
            %IRESOLVERELATIVEPATH Resolve a safe project-relative path.

            relativePath = char(string(relativePath));
            normalized = strrep(relativePath, '\', '/');

            if isempty(normalized) || startsWith(normalized, '/') || ...
                    ~isempty(regexp(normalized, '^[A-Za-z]:', 'once'))
                error('Umitoolbox:UMITProjectStore:unsafePath', ...
                    'Expected a non-empty project-relative path: %s', relativePath);
            end

            parts = strsplit(normalized, '/');
            if any(cellfun(@isempty, parts)) || ...
                    any(ismember(parts, {'.', '..'}))
                error('Umitoolbox:UMITProjectStore:unsafePath', ...
                    'Relative path contains an unsafe segment: %s', relativePath);
            end

            path = fullfile(obj.ProjectRoot, parts{:});
        end

        function path = iResolveRelativePathFromBase(~, basePath, relativeTail)
            %IRESOLVERELATIVEPATHFROMBASE Resolve a validated relative tail.

            parts = strsplit(strrep(relativeTail, '\', '/'), '/');
            if any(cellfun(@isempty, parts)) || any(ismember(parts, {'.', '..'}))
                error('Umitoolbox:UMITProjectStore:unsafePath', ...
                    'Relative tail contains an unsafe segment: %s', relativeTail);
            end
            path = fullfile(basePath, parts{:});
        end

        function tail = iRelativeTail(~, fullRelativePath, prefixRelativePath)
            %IRELATIVETAIL Return a relative path after a canonical prefix.

            fullPath = strrep(char(string(fullRelativePath)), '\', '/');
            prefix = strrep(char(string(prefixRelativePath)), '\', '/');
            expected = [prefix, '/'];
            if ~startsWith(fullPath, expected, 'IgnoreCase', true)
                error('Umitoolbox:UMITProjectStore:pathPrefixMismatch', ...
                    'Path "%s" is not below "%s".', fullPath, prefix);
            end
            tail = fullPath(numel(expected)+1:end);
        end

        function registry = iReplaceRegistryPrefix(obj, registry, oldPrefix, newPrefix)
            %IREPLACEREGISTRYPREFIX Update relative paths after parent rename.

            for iRecord = 1:numel(registry)
                registry(iRecord).relativePath = obj.iReplaceRelativePrefix( ...
                    registry(iRecord).relativePath, oldPrefix, newPrefix);
            end
        end

        function registry = iReplaceResourcePrefix(obj, registry, oldPrefix, newPrefix)
            %IREPLACERESOURCEPREFIX Update resource paths after parent rename.

            for iRecord = 1:numel(registry)
                registry(iRecord).relativePath = obj.iReplaceRelativePrefix( ...
                    registry(iRecord).relativePath, oldPrefix, newPrefix);
            end
        end

        function output = iReplaceRelativePrefix(~, input, oldPrefix, newPrefix)
            %IREPLACERELATIVEPREFIX Replace one case-insensitive path prefix.

            input = strrep(char(string(input)), '\', '/');
            oldPrefix = strrep(char(string(oldPrefix)), '\', '/');
            newPrefix = strrep(char(string(newPrefix)), '\', '/');

            if strcmpi(input, oldPrefix)
                output = newPrefix;
                return
            end

            expected = [oldPrefix, '/'];
            if ~startsWith(input, expected, 'IgnoreCase', true)
                error('Umitoolbox:UMITProjectStore:pathPrefixMismatch', ...
                    'Path "%s" does not start with "%s".', input, oldPrefix);
            end
            output = [newPrefix, input(numel(oldPrefix)+1:end)];
        end

        function report = iValidateSubjectNode(obj, report, subjectRecord, ProjectInfo, mode, scope)
            %IVALIDATESUBJECTNODE Validate one subject and its child sessions.
            %   SCOPE ([] = all) limits the sessions validated (see validate).

            if nargin < 6
                scope = [];
            end

            subjectRel = subjectRecord.relativePath;
            subjectPath = obj.iResolveRelativePath(subjectRel);

            if ~isfolder(subjectPath)
                report = obj.iAddIssue(report, 'error', ...
                    'missing_subject_folder', 'subject', subjectRel, ...
                    sprintf('Registered subject folder is missing: %s', subjectRel), ...
                    false);
                return
            end

            metadataRel = UMITProjectStore.iJoinRelative( ...
                subjectRel, obj.Schema.files.subjectMetadata);
            try
                SubjectInfo = obj.iLoadMetadata( ...
                    fullfile(subjectPath, obj.Schema.files.subjectMetadata), ...
                    obj.Schema.metadataVariables.subject);
            catch ME
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_subject_metadata', 'subject', metadataRel, ...
                    ME.message, false);
                return
            end

            report = obj.iValidateMetadataFields(report, SubjectInfo, ...
                obj.Schema.requiredSubjectFields, 'subject', metadataRel);

            if isfield(SubjectInfo, 'uuid') && ...
                    ~strcmpi(SubjectInfo.uuid, subjectRecord.uuid)
                report = obj.iAddIssue(report, 'error', ...
                    'subject_uuid_mismatch', 'subject', subjectRel, ...
                    'Subject UUID does not match the project registry.', false);
            end

            if isfield(SubjectInfo, 'subjectID') && ...
                    ~strcmp(SubjectInfo.subjectID, subjectRecord.id)
                report = obj.iAddIssue(report, 'error', ...
                    'subject_id_mismatch', 'subject', subjectRel, ...
                    'Subject ID does not match its registry record.', false);
            end

            try
                obj.iNormalizeManagedID(subjectRecord.id, 'subjectID');
            catch ME
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_subject_id', 'subject', subjectRel, ME.message, false);
            end

            requiredFolders = { ...
                UMITProjectStore.iJoinRelative(obj.Schema.folders.references, ...
                obj.Schema.folders.image, obj.Schema.folders.active), ...
                UMITProjectStore.iJoinRelative(obj.Schema.folders.references, ...
                obj.Schema.folders.image, obj.Schema.folders.archive), ...
                obj.Schema.folders.sessions};
            report = obj.iValidateRequiredNodeFolders( ...
                report, subjectPath, subjectRel, 'subject', requiredFolders);
            report = obj.iValidateAllowedChildren(report, ...
                subjectPath, subjectRel, 'subject', { ...
                obj.Schema.files.subjectMetadata, ...
                obj.Schema.folders.references, ...
                obj.Schema.folders.sessions});
            report = obj.iValidateAllowedChildren(report, ...
                fullfile(subjectPath, obj.Schema.folders.references), ...
                UMITProjectStore.iJoinRelative(subjectRel, ...
                obj.Schema.folders.references), 'subject', ...
                {obj.Schema.folders.image});
            report = obj.iValidateAllowedChildren(report, ...
                fullfile(subjectPath, obj.Schema.folders.references, ...
                obj.Schema.folders.image), ...
                UMITProjectStore.iJoinRelative(subjectRel, ...
                obj.Schema.folders.references, obj.Schema.folders.image), ...
                'subject', {obj.Schema.folders.active, ...
                obj.Schema.folders.archive});

            if isfield(SubjectInfo, 'resourceRegistry')
                report = obj.iValidateResourceRegistry(report, ...
                    SubjectInfo, subjectRel, 'subject', mode);
            end

            if isfield(SubjectInfo, 'sessionRegistry')
                report = obj.iValidateRegistry(report, ...
                    SubjectInfo.sessionRegistry, 'session', ...
                    UMITProjectStore.iJoinRelative(subjectRel, obj.Schema.folders.sessions));

                if obj.iIsRegistryStruct(SubjectInfo.sessionRegistry)
                    for iSession = 1:numel(SubjectInfo.sessionRegistry)
                        if ~UMITProjectStore.iScopeHasSession(scope, ...
                                subjectRecord.uuid, SubjectInfo.sessionRegistry(iSession).uuid)
                            continue
                        end
                        try
                            report = obj.iValidateSessionNode(report, ...
                                SubjectInfo.sessionRegistry(iSession), ...
                                SubjectInfo, subjectRecord, ProjectInfo, mode);
                        catch ME
                            report = obj.iAddIssue(report, 'error', ...
                                'session_validation_failed', 'session', ...
                                SubjectInfo.sessionRegistry(iSession).relativePath, ...
                                ME.message, false);
                        end
                    end
                end
            end

            report = obj.iValidateUnregisteredSessionFolders( ...
                report, SubjectInfo, subjectRel);
        end

        function report = iValidateSessionNode(obj, report, sessionRecord, SubjectInfo, subjectRecord, ~, mode)
            %IVALIDATESESSIONNODE Validate one session metadata node.

            sessionRel = sessionRecord.relativePath;
            sessionPath = obj.iResolveRelativePath(sessionRel);

            if ~isfolder(sessionPath)
                report = obj.iAddIssue(report, 'error', ...
                    'missing_session_folder', 'session', sessionRel, ...
                    sprintf('Registered session folder is missing: %s', sessionRel), ...
                    false);
                return
            end

            metadataRel = UMITProjectStore.iJoinRelative( ...
                sessionRel, obj.Schema.files.sessionMetadata);
            try
                SessionInfo = obj.iLoadMetadata( ...
                    fullfile(sessionPath, obj.Schema.files.sessionMetadata), ...
                    obj.Schema.metadataVariables.session);
            catch ME
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_session_metadata', 'session', metadataRel, ...
                    ME.message, false);
                return
            end

            report = obj.iValidateMetadataFields(report, SessionInfo, ...
                obj.Schema.requiredSessionFields, 'session', metadataRel);

            if isfield(SessionInfo, 'uuid') && ...
                    ~strcmpi(SessionInfo.uuid, sessionRecord.uuid)
                report = obj.iAddIssue(report, 'error', ...
                    'session_uuid_mismatch', 'session', sessionRel, ...
                    'Session UUID does not match the subject registry.', false);
            end

            if isfield(SessionInfo, 'sessionID') && ...
                    ~strcmp(SessionInfo.sessionID, sessionRecord.id)
                report = obj.iAddIssue(report, 'error', ...
                    'session_id_mismatch', 'session', sessionRel, ...
                    'Session ID does not match its registry record.', false);
            end

            if isfield(SessionInfo, 'subjectUUID') && ...
                    isfield(SubjectInfo, 'uuid') && ...
                    ~strcmpi(SessionInfo.subjectUUID, SubjectInfo.uuid)
                report = obj.iAddIssue(report, 'error', ...
                    'session_subject_uuid_mismatch', 'session', sessionRel, ...
                    'Session subjectUUID does not match its parent subject.', false);
            end

            if isfield(SessionInfo, 'subjectID') && ...
                    ~strcmp(SessionInfo.subjectID, subjectRecord.id)
                report = obj.iAddIssue(report, 'error', ...
                    'session_subject_id_mismatch', 'session', sessionRel, ...
                    'Session subjectID does not match its parent subject.', false);
            end

            try
                obj.iNormalizeManagedID(sessionRecord.id, 'sessionID');
            catch ME
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_session_id', 'session', sessionRel, ME.message, false);
            end

            % 'transforms' is a legacy, optional child (retired
            % registration-transform resources); its content is not checked.
            report = obj.iValidateAllowedChildren(report, ...
                sessionPath, sessionRel, 'session', { ...
                obj.Schema.files.sessionMetadata, ...
                obj.Schema.folders.transforms});

            if isfield(SessionInfo, 'resourceRegistry')
                report = obj.iValidateResourceRegistry(report, ...
                    SessionInfo, sessionRel, 'session', mode);
            end

            bindingRoles = { ...
                'rawDataFolder', 'processedDataFolder'};
            for iRole = 1:numel(bindingRoles)
                pathField = bindingRoles{iRole};
                [~, bindingField] = ...
                    obj.iGetBindingRoleFields(pathField);

                hasPath = isfield(SessionInfo, pathField) && ...
                    ~isempty(SessionInfo.(pathField));
                hasBindingUUID = isfield(SessionInfo, bindingField) && ...
                    ~isempty(SessionInfo.(bindingField));

                if ~hasPath && ~hasBindingUUID
                    if strcmp(pathField, 'processedDataFolder')
                        report = obj.iAddIssue(report, 'error', ...
                            'session_save_folder_required', 'session', ...
                            sessionRel, ...
                            ['A persisted session must retain a nonempty ' ...
                             'processedDataFolder SaveFolder binding.'], ...
                            false);
                    end
                    continue
                end

                if xor(hasPath, hasBindingUUID)
                    report = obj.iAddIssue(report, 'error', ...
                        'incomplete_data_folder_binding', 'session', ...
                        sessionRel, sprintf( ...
                        'Fields %s and %s must either both be empty or both be set.', ...
                        pathField, bindingField), true);
                    continue
                end

                if isempty(regexp(SessionInfo.(bindingField), ...
                        obj.Schema.namingRules.uuidPattern, 'once'))
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_binding_uuid', 'session', sessionRel, ...
                        sprintf('Field %s is not a valid UUID.', ...
                        bindingField), false);
                    continue
                end

                try
                    obj.iAssertExternalDataPath( ...
                        SessionInfo.(pathField), pathField);
                catch ME
                    report = obj.iAddIssue(report, 'error', ...
                        'internal_data_folder_reference', 'session', ...
                        sessionRel, ME.message, true);
                    continue
                end

                if ~isfolder(SessionInfo.(pathField))
                    report = obj.iAddIssue(report, 'warning', ...
                        'external_data_folder_missing', 'session', ...
                        sessionRel, sprintf( ...
                        'External folder does not exist: %s', ...
                        SessionInfo.(pathField)), true);
                    continue
                end

                if strcmp(pathField, 'processedDataFolder')
                    try
                        obj.iAssertValidSaveFolderDataset( ...
                            SessionInfo.(pathField));
                    catch ME
                        if strcmp(ME.identifier, ...
                                'Umitoolbox:UMITProjectStore:notADataset')
                            issueCode = 'not_a_dataset';
                        else
                            issueCode = 'invalid_acq_infos';
                        end
                        report = obj.iAddIssue(report, 'error', ...
                            issueCode, 'session', sessionRel, ...
                            ME.message, false);
                    end
                else
                    try
                        obj.iAssertValidRawDataFolder( ...
                            SessionInfo.(pathField));
                    catch ME
                        report = obj.iAddIssue(report, 'error', ...
                            'missing_raw_data_files', 'session', ...
                            sessionRel, ME.message, false);
                    end
                    authority = obj.iDataParamsRawFolder(SessionInfo);
                    if ~isempty(authority) && ~strcmp( ...
                            UMITProjectStore.iNormalizeComparisonPath(authority), ...
                            UMITProjectStore.iNormalizeComparisonPath(SessionInfo.(pathField)))
                        report = obj.iAddIssue(report, 'warning', ...
                            'raw_folder_mismatch', 'session', sessionRel, sprintf( ...
                            ['The session''s RawFolder "%s" differs from the ' ...
                             'SaveFolder''s DataParams.folders.RawFolder "%s" (the ' ...
                             'authority).'], SessionInfo.(pathField), authority), true);
                    end
                end

                if strcmp(UMITProjectStore.iBindingState(SessionInfo.(pathField)), 'none')
                    report = obj.iAddIssue(report, 'error', ...
                        'missing_project_binding', 'session', ...
                        sessionRel, sprintf( ...
                        'Bound folder is not bound in its %s: %s', ...
                        obj.Schema.files.projectBinding, ...
                        SessionInfo.(pathField)), true);
                    continue
                end

                try
                    ProjectBinding = ...
                        UMITProjectStore.readProjectBinding( ...
                        SessionInfo.(pathField));
                    obj.iAssertBindingMatchesSession( ...
                        ProjectBinding, SessionInfo, SubjectInfo, ...
                        pathField, bindingField, ...
                        SessionInfo.(pathField));
                catch ME
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_project_binding', 'session', ...
                        sessionRel, ME.message, true);
                end
            end
        end

        function report = iValidateRegistry(obj, report, registry, nodeType, expectedParentRel)
            %IVALIDATEREGISTRY Validate one registry struct array.

            if ~obj.iIsRegistryStruct(registry)
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_registry', nodeType, expectedParentRel, ...
                    sprintf('%s registry has an invalid structure.', nodeType), ...
                    false);
                return
            end

            if isempty(registry)
                return
            end

            validRecord = true(1, numel(registry));
            textFields = obj.Schema.registryFields;

            for iRecord = 1:numel(registry)
                record = registry(iRecord);
                for iField = 1:numel(textFields)
                    fieldName = textFields{iField};
                    value = record.(fieldName);
                    if ~(ischar(value) || ...
                            (isstring(value) && isscalar(value)))
                        report = obj.iAddIssue(report, 'error', ...
                            'invalid_registry_field', nodeType, ...
                            expectedParentRel, ...
                            sprintf('Registry field "%s" must be a text scalar.', ...
                            fieldName), false);
                        validRecord(iRecord) = false;
                    end
                end

                if ~validRecord(iRecord)
                    continue
                end

                try
                    obj.iNormalizeManagedID(record.id, [nodeType, 'ID']);
                catch ME
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_registry_id', nodeType, record.relativePath, ...
                        ME.message, false);
                end

                expectedRel = UMITProjectStore.iJoinRelative( ...
                    expectedParentRel, record.id);
                if ~strcmp(record.relativePath, expectedRel)
                    report = obj.iAddIssue(report, 'error', ...
                        'noncanonical_registry_path', nodeType, ...
                        record.relativePath, ...
                        sprintf('Expected canonical path: %s', expectedRel), ...
                        true);
                end
            end

            validRegistry = registry(validRecord);
            if isempty(validRegistry)
                return
            end

            ids = {validRegistry.id};
            uuids = {validRegistry.uuid};
            if numel(unique(lower(string(ids)))) ~= numel(ids)
                report = obj.iAddIssue(report, 'error', ...
                    'duplicate_registry_id', nodeType, expectedParentRel, ...
                    sprintf('%s registry contains duplicate IDs.', nodeType), false);
            end
            if numel(unique(string(uuids))) ~= numel(uuids)
                report = obj.iAddIssue(report, 'error', ...
                    'duplicate_registry_uuid', nodeType, expectedParentRel, ...
                    sprintf('%s registry contains duplicate UUIDs.', nodeType), false);
            end
        end

        function report = iValidateResourceRegistry(obj, report, metadata, ownerRel, ownerType, mode)
            %IVALIDATERESOURCEREGISTRY Validate resources owned by one node.

            registry = metadata.resourceRegistry;
            if ~obj.iIsResourceRegistryStruct(registry)
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_resource_registry', ownerType, ownerRel, ...
                    'resourceRegistry has an invalid structure.', false);
                return
            end

            resourceTypes = fieldnames(obj.Schema.resourceTypes);
            validRecord = true(1, numel(registry));
            textFields = { ...
                'uuid', 'type', 'fileName', 'relativePath', ...
                'displayName', 'description', 'status', ...
                'checksum', 'sourceFile'};

            for iRecord = 1:numel(registry)
                record = registry(iRecord);

                for iField = 1:numel(textFields)
                    fieldName = textFields{iField};
                    value = record.(fieldName);
                    if ~(ischar(value) || ...
                            (isstring(value) && isscalar(value)))
                        report = obj.iAddIssue(report, 'error', ...
                            'invalid_resource_field', ownerType, ownerRel, ...
                            sprintf('Resource field "%s" must be a text scalar.', ...
                            fieldName), false);
                        validRecord(iRecord) = false;
                    end
                end

                dateFields = {'createdOn', 'modifiedOn', 'archivedOn'};
                for iField = 1:numel(dateFields)
                    fieldName = dateFields{iField};
                    value = record.(fieldName);
                    if ~(isdatetime(value) && isscalar(value))
                        report = obj.iAddIssue(report, 'error', ...
                            'invalid_resource_date', ownerType, ownerRel, ...
                            sprintf('Resource field "%s" must be a datetime scalar.', ...
                            fieldName), false);
                        validRecord(iRecord) = false;
                    end
                end

                if ~validRecord(iRecord)
                    continue
                end

                if ismember(record.type, obj.Schema.retiredResourceTypes)
                    % Legacy record (e.g. a registration transform written
                    % before section 4 Q5): tolerated, never used.
                    report = obj.iAddIssue(report, 'warning', ...
                        'retired_resource_type', ownerType, ...
                        record.relativePath, ...
                        sprintf(['Resource %s has the retired type "%s"; it is ' ...
                        'kept on disk but no longer used.'], record.uuid, record.type), ...
                        false);
                    validRecord(iRecord) = false;
                    continue
                end

                if ~ismember(record.type, resourceTypes)
                    report = obj.iAddIssue(report, 'error', ...
                        'unsupported_resource_type', ownerType, ...
                        record.relativePath, ...
                        sprintf('Unsupported resource type: %s', record.type), ...
                        false);
                    validRecord(iRecord) = false;
                    continue
                end

                resourceDef = obj.Schema.resourceTypes.(record.type);
                if ~strcmp(resourceDef.ownerType, ownerType)
                    report = obj.iAddIssue(report, 'error', ...
                        'resource_owner_mismatch', ownerType, ...
                        record.relativePath, ...
                        sprintf('Resource type "%s" cannot belong to "%s".', ...
                        record.type, ownerType), false);
                end

                if ~ismember(record.status, obj.Schema.resourceStates)
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_resource_status', ownerType, ...
                        record.relativePath, ...
                        sprintf('Invalid resource status: %s', record.status), ...
                        false);
                    validRecord(iRecord) = false;
                    continue
                end

                [~, ~, extension] = fileparts(record.fileName);
                if ~any(strcmpi(extension, resourceDef.allowedExtensions))
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_resource_extension', ownerType, ...
                        record.relativePath, ...
                        sprintf('Invalid extension for resource type "%s".', ...
                        record.type), false);
                end

                if isempty(regexp(record.checksum, '^[0-9a-f]{64}$', 'once'))
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_resource_checksum', ownerType, ...
                        record.relativePath, ...
                        'Resource checksum must be a lowercase SHA-256 digest.', ...
                        false);
                end

                if strcmp(record.status, 'archived')
                    stateFolder = obj.Schema.folders.archive;
                    if isnat(record.archivedOn)
                        report = obj.iAddIssue(report, 'error', ...
                            'missing_archive_timestamp', ownerType, ...
                            record.relativePath, ...
                            'Archived resource is missing archivedOn.', false);
                    end
                else
                    stateFolder = obj.Schema.folders.active;
                    if ~isnat(record.archivedOn)
                        report = obj.iAddIssue(report, 'error', ...
                            'unexpected_archive_timestamp', ownerType, ...
                            record.relativePath, ...
                            'Non-archived resource must use archivedOn=NaT.', ...
                            false);
                    end
                end

                expectedRel = obj.iBuildResourceRelativePath( ...
                    ownerRel, record.type, stateFolder, record.fileName);
                if ~strcmp(record.relativePath, expectedRel)
                    report = obj.iAddIssue(report, 'error', ...
                        'noncanonical_resource_path', ownerType, ...
                        record.relativePath, ...
                        sprintf('Expected resource path: %s', expectedRel), true);
                end

                try
                    resourcePath = obj.iResolveRelativePath(record.relativePath);
                catch ME
                    report = obj.iAddIssue(report, 'error', ...
                        'unsafe_resource_path', ownerType, ...
                        record.relativePath, ME.message, false);
                    continue
                end

                if ~isfile(resourcePath)
                    report = obj.iAddIssue(report, 'error', ...
                        'missing_resource_file', ownerType, ...
                        record.relativePath, ...
                        'Registered resource file is missing.', false);
                elseif strcmp(mode, 'full')
                    try
                        checksum = computeFileChecksum(resourcePath);
                        if ~strcmp(checksum, record.checksum)
                            report = obj.iAddIssue(report, 'error', ...
                                'resource_checksum_mismatch', ownerType, ...
                                record.relativePath, ...
                                'Resource checksum does not match its registry.', ...
                                false);
                        end
                    catch ME
                        report = obj.iAddIssue(report, 'error', ...
                            'resource_checksum_failed', ownerType, ...
                            record.relativePath, ME.message, false);
                    end

                    try
                        if obj.iResourceTypeHasPayloadContract(record.type)
                            resourceProbe = load(resourcePath, '-mat');
                            obj.iValidateManagedResourcePayload( ...
                                record.type, resourceProbe, resourcePath);
                        end
                    catch ME
                        report = obj.iAddIssue(report, 'error', ...
                            'invalid_resource_payload', ownerType, ...
                            record.relativePath, ME.message, false);
                    end
                end
            end

            validRegistry = registry(validRecord);
            if ~isempty(validRegistry)
                uuids = {validRegistry.uuid};
                if numel(unique(string(uuids))) ~= numel(uuids)
                    report = obj.iAddIssue(report, 'error', ...
                        'duplicate_resource_uuid', ownerType, ownerRel, ...
                        'resourceRegistry contains duplicate UUIDs.', false);
                end
            end

            % Validate active pointers and status consistency for every
            % resource type owned by this metadata node.
            for iType = 1:numel(resourceTypes)
                typeName = resourceTypes{iType};
                def = obj.Schema.resourceTypes.(typeName);
                if ~strcmp(def.ownerType, ownerType)
                    continue
                end

                pointerField = def.activePointerField;
                if ~isfield(metadata, pointerField)
                    report = obj.iAddIssue(report, 'error', ...
                        'missing_active_pointer', ownerType, ownerRel, ...
                        sprintf('Metadata is missing "%s".', pointerField), false);
                    continue
                end

                pointerUUID = metadata.(pointerField);
                if ~(ischar(pointerUUID) || ...
                        (isstring(pointerUUID) && isscalar(pointerUUID)))
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_active_pointer_type', ownerType, ownerRel, ...
                        sprintf('Active pointer "%s" must be a text scalar.', ...
                        pointerField), false);
                    continue
                end
                pointerUUID = char(string(pointerUUID));

                typeIndices = [];
                for iRecord = find(validRecord)
                    if strcmp(registry(iRecord).type, typeName)
                        typeIndices(end+1) = iRecord; %#ok<AGROW>
                    end
                end

                activeIndices = [];
                for iRecord = typeIndices
                    if strcmp(registry(iRecord).status, 'active')
                        activeIndices(end+1) = iRecord; %#ok<AGROW>
                    end
                end

                if isempty(pointerUUID)
                    if ~isempty(activeIndices)
                        report = obj.iAddIssue(report, 'error', ...
                            'orphan_active_status', ownerType, ownerRel, ...
                            sprintf('Active %s exists without an active pointer.', ...
                            typeName), false);
                    end
                else
                    pointerIndex = obj.iFindResourceIndex(registry, pointerUUID);
                    if isempty(pointerIndex) || ...
                            ~validRecord(pointerIndex) || ...
                            ~strcmp(registry(pointerIndex).type, typeName) || ...
                            ~strcmp(registry(pointerIndex).status, 'active')
                        report = obj.iAddIssue(report, 'error', ...
                            'invalid_active_pointer', ownerType, ownerRel, ...
                            sprintf('Active pointer for %s is invalid.', typeName), ...
                            false);
                    end
                    if numel(activeIndices) ~= 1
                        report = obj.iAddIssue(report, 'error', ...
                            'invalid_active_count', ownerType, ownerRel, ...
                            sprintf('Exactly one active %s is required when its pointer is set.', ...
                            typeName), false);
                    end
                end
            end

            report = obj.iValidateUnregisteredResourceFiles( ...
                report, validRegistry, ownerRel, ownerType);
        end

        function report = iValidateRequiredNodeFolders(obj, report, nodePath, nodeRel, nodeType, requiredFolders)
            %IVALIDATEREQUIREDNODEFOLDERS Validate canonical child folders.

            for iFolder = 1:numel(requiredFolders)
                folderRel = requiredFolders{iFolder};
                parts = strsplit(folderRel, '/');
                fullPath = fullfile(nodePath, parts{:});
                if ~isfolder(fullPath)
                    rel = UMITProjectStore.iJoinRelative(nodeRel, folderRel);
                    report = obj.iAddIssue(report, 'error', ...
                        'missing_required_folder', nodeType, rel, ...
                        sprintf('Required folder is missing: %s', rel), true);
                end
            end
        end

        function report = iValidateUnregisteredTopNodes(obj, report, ProjectInfo)
            %IVALIDATEUNREGISTEREDTOPNODES Detect unregistered subject folders.

            report = obj.iCompareFolderNamesToRegistry(report, ...
                obj.Schema.folders.subjects, ProjectInfo.subjectRegistry, ...
                'subject');
        end

        function report = iValidateUnregisteredSessionFolders(obj, report, SubjectInfo, subjectRel)
            %IVALIDATEUNREGISTEREDSESSIONFOLDERS Detect unknown session folders.

            sessionsRel = UMITProjectStore.iJoinRelative( ...
                subjectRel, obj.Schema.folders.sessions);
            report = obj.iCompareFolderNamesToRegistry(report, ...
                sessionsRel, SubjectInfo.sessionRegistry, 'session');
        end

        function report = iCompareFolderNamesToRegistry(obj, report, parentRel, registry, nodeType)
            %ICOMPAREFOLDERNAMETOREGISTRY Compare physical folders to registry IDs.

            parentPath = obj.iResolveRelativePath(parentRel);
            if ~isfolder(parentPath) || ~obj.iIsRegistryStruct(registry)
                return
            end

            entries = dir(parentPath);
            entries = entries(~ismember({entries.name}, {'.', '..'}));
            registered = {registry.id};

            for iEntry = 1:numel(entries)
                name = entries(iEntry).name;
                if any(strcmpi(name, obj.Schema.allowedIncidentalNames))
                    continue
                end

                rel = UMITProjectStore.iJoinRelative(parentRel, name);
                if ~entries(iEntry).isdir
                    report.unregisteredItems{end+1} = rel;
                    report = obj.iAddIssue(report, 'error', ...
                        'unexpected_registry_file', nodeType, rel, ...
                        sprintf('Unexpected file in managed registry folder: %s', ...
                        rel), false);
                    continue
                end

                matchIndex = find(strcmpi(name, registered), 1, 'first');
                if isempty(matchIndex)
                    report.unregisteredItems{end+1} = rel;
                    report = obj.iAddIssue(report, 'error', ...
                        'unregistered_folder', nodeType, rel, ...
                        sprintf('Unregistered managed folder: %s', rel), false);
                elseif ~strcmp(name, registered{matchIndex})
                    report = obj.iAddIssue(report, 'error', ...
                        'noncanonical_folder_case', nodeType, rel, ...
                        sprintf('Folder name case does not match ID "%s".', ...
                        registered{matchIndex}), true);
                end
            end
        end

        function report = iValidateUnregisteredResourceFiles(obj, report, registry, ownerRel, ownerType)
            %IVALIDATEUNREGISTEREDRESOURCEFILES Detect files absent from registry.

            resourceTypes = fieldnames(obj.Schema.resourceTypes);
            for iType = 1:numel(resourceTypes)
                typeName = resourceTypes{iType};
                def = obj.Schema.resourceTypes.(typeName);
                if ~strcmp(def.ownerType, ownerType)
                    continue
                end

                for stateFolder = {obj.Schema.folders.active, obj.Schema.folders.archive}
                    folderRel = UMITProjectStore.iJoinRelative( ...
                        ownerRel, def.relativeFolderParts{:}, stateFolder{1});
                    folderPath = obj.iResolveRelativePath(folderRel);
                    if ~isfolder(folderPath)
                        continue
                    end

                    entries = dir(folderPath);
                    entries = entries(~ismember({entries.name}, {'.', '..'}));
                    for iEntry = 1:numel(entries)
                        if any(strcmpi(entries(iEntry).name, ...
                                obj.Schema.allowedIncidentalNames))
                            continue
                        end

                        rel = UMITProjectStore.iJoinRelative(folderRel, entries(iEntry).name);
                        if entries(iEntry).isdir
                            report.unregisteredItems{end+1} = rel;
                            report = obj.iAddIssue(report, 'error', ...
                                'unexpected_resource_folder', ownerType, rel, ...
                                sprintf('Unexpected folder in resource location: %s', rel), ...
                                false);
                        elseif ~any(strcmp({registry.relativePath}, rel))
                            report.unregisteredItems{end+1} = rel;
                            report = obj.iAddIssue(report, 'error', ...
                                'unregistered_resource_file', ownerType, rel, ...
                                sprintf('Resource file is not registered: %s', rel), ...
                                false);
                        end
                    end
                end
            end
        end

        function report = iValidateAllowedChildren(obj, report, folderPath, folderRel, nodeType, allowedNames)
            %IVALIDATEALLOWEDCHILDREN Reject unexpected files and folders.

            if ~isfolder(folderPath)
                return
            end

            entries = dir(folderPath);
            names = setdiff({entries.name}, {'.', '..'});

            for iEntry = 1:numel(names)
                name = names{iEntry};
                if any(strcmpi(name, obj.Schema.allowedIncidentalNames))
                    continue
                end

                if any(strcmp(name, allowedNames))
                    continue
                end

                relativePath = name;
                if ~isempty(folderRel)
                    relativePath = UMITProjectStore.iJoinRelative(folderRel, name);
                end

                if any(strcmpi(name, allowedNames))
                    code = 'noncanonical_item_case';
                    message = sprintf( ...
                        'Managed item uses noncanonical letter case: %s', ...
                        relativePath);
                    repairable = true;
                else
                    code = 'unexpected_managed_item';
                    message = sprintf( ...
                        'Unexpected file or folder in managed location: %s', ...
                        relativePath);
                    repairable = false;
                end

                report.unregisteredItems{end+1} = relativePath;
                report = obj.iAddIssue(report, 'error', code, nodeType, ...
                    relativePath, message, repairable);
            end
        end

        function report = iValidateGlobalUUIDs(obj, report, ProjectInfo)
            %IVALIDATEGLOBALUUIDS Validate UUID format and project-wide uniqueness.

            uuids = {};
            labels = {};

            if isfield(ProjectInfo, 'projectUUID')
                uuids{end+1} = ProjectInfo.projectUUID; %#ok<AGROW>
                labels{end+1} = 'project'; %#ok<AGROW>
            end

            if isfield(ProjectInfo, 'subjectRegistry') && ...
                    obj.iIsRegistryStruct(ProjectInfo.subjectRegistry)
                for iSubject = 1:numel(ProjectInfo.subjectRegistry)
                    subjectRecord = ProjectInfo.subjectRegistry(iSubject);
                    uuids{end+1} = subjectRecord.uuid; %#ok<AGROW>
                    labels{end+1} = subjectRecord.relativePath; %#ok<AGROW>

                    try
                        subjectPath = obj.iResolveRelativePath( ...
                            subjectRecord.relativePath);
                        SubjectInfo = obj.iLoadMetadata( ...
                            fullfile(subjectPath, ...
                            obj.Schema.files.subjectMetadata), ...
                            obj.Schema.metadataVariables.subject);
                    catch
                        continue
                    end

                    if isfield(SubjectInfo, 'resourceRegistry') && ...
                            obj.iIsResourceRegistryStruct( ...
                            SubjectInfo.resourceRegistry)
                        for iResource = 1:numel(SubjectInfo.resourceRegistry)
                            uuids{end+1} = ...
                                SubjectInfo.resourceRegistry(iResource).uuid; %#ok<AGROW>
                            labels{end+1} = ...
                                SubjectInfo.resourceRegistry(iResource).relativePath; %#ok<AGROW>
                        end
                    end

                    if isfield(SubjectInfo, 'sessionRegistry') && ...
                            obj.iIsRegistryStruct(SubjectInfo.sessionRegistry)
                        for iSession = 1:numel(SubjectInfo.sessionRegistry)
                            sessionRecord = SubjectInfo.sessionRegistry(iSession);
                            uuids{end+1} = sessionRecord.uuid; %#ok<AGROW>
                            labels{end+1} = sessionRecord.relativePath; %#ok<AGROW>

                            try
                                sessionPath = obj.iResolveRelativePath( ...
                                    sessionRecord.relativePath);
                                SessionInfo = obj.iLoadMetadata( ...
                                    fullfile(sessionPath, ...
                                    obj.Schema.files.sessionMetadata), ...
                                    obj.Schema.metadataVariables.session);
                            catch
                                continue
                            end

                            bindingUUIDFields = { ...
                                'rawDataBindingUUID', ...
                                'processedDataBindingUUID'};
                            hasSharedDataFolderBinding = ...
                                isfield(SessionInfo, ...
                                'rawDataBindingUUID') && ...
                                isfield(SessionInfo, ...
                                'processedDataBindingUUID') && ...
                                isfield(SessionInfo, 'rawDataFolder') && ...
                                isfield(SessionInfo, ...
                                'processedDataFolder') && ...
                                ~isempty(SessionInfo.rawDataBindingUUID) && ...
                                strcmpi(SessionInfo.rawDataBindingUUID, ...
                                SessionInfo.processedDataBindingUUID) && ...
                                ~isempty(SessionInfo.rawDataFolder) && ...
                                strcmp( ...
                                UMITProjectStore.iNormalizeComparisonPath( ...
                                SessionInfo.rawDataFolder), ...
                                UMITProjectStore.iNormalizeComparisonPath( ...
                                SessionInfo.processedDataFolder));
                            if hasSharedDataFolderBinding
                                bindingUUIDFields = { ...
                                    'processedDataBindingUUID'};
                            end
                            for iBindingUUID = 1:numel(bindingUUIDFields)
                                bindingField = ...
                                    bindingUUIDFields{iBindingUUID};
                                if isfield(SessionInfo, bindingField) && ...
                                        ~isempty(SessionInfo.(bindingField))
                                    uuids{end+1} = ...
                                        SessionInfo.(bindingField); %#ok<AGROW>
                                    labels{end+1} = sprintf( ...
                                        '%s/%s', ...
                                        sessionRecord.relativePath, ...
                                        bindingField); %#ok<AGROW>
                                end
                            end

                            if isfield(SessionInfo, 'resourceRegistry') && ...
                                    obj.iIsResourceRegistryStruct( ...
                                    SessionInfo.resourceRegistry)
                                for iResource = 1:numel(SessionInfo.resourceRegistry)
                                    uuids{end+1} = ...
                                        SessionInfo.resourceRegistry(iResource).uuid; %#ok<AGROW>
                                    labels{end+1} = ...
                                        SessionInfo.resourceRegistry(iResource).relativePath; %#ok<AGROW>
                                end
                            end
                        end
                    end
                end
            end

            normalizedUUIDs = strings(size(uuids));
            validUUID = false(size(uuids));
            for iUUID = 1:numel(uuids)
                value = uuids{iUUID};
                if ischar(value) || (isstring(value) && isscalar(value))
                    textValue = char(string(value));
                    normalizedUUIDs(iUUID) = lower(string(textValue));
                    validUUID(iUUID) = ~isempty(regexp( ...
                        textValue, obj.Schema.namingRules.uuidPattern, 'once'));
                end

                if ~validUUID(iUUID)
                    report = obj.iAddIssue(report, 'error', ...
                        'invalid_uuid', 'project', labels{iUUID}, ...
                        sprintf('Invalid UUID at %s.', labels{iUUID}), false);
                end
            end

            validValues = normalizedUUIDs(validUUID);
            if numel(unique(validValues)) ~= numel(validValues)
                [uniqueValues, ~, groupIndex] = unique(validValues);
                counts = accumarray(groupIndex(:), 1);
                duplicates = uniqueValues(counts > 1);
                for iDuplicate = 1:numel(duplicates)
                    report = obj.iAddIssue(report, 'error', ...
                        'duplicate_global_uuid', 'project', '', ...
                        sprintf('UUID is used by multiple project items: %s', ...
                        char(duplicates(iDuplicate))), false);
                end
            end
        end

        function report = iValidateMetadataFields(obj, report, metadata, requiredFields, nodeType, relativePath)
            %IVALIDATEMETADATAFIELDS Validate required metadata field presence.

            if ~isstruct(metadata) || ~isscalar(metadata)
                report = obj.iAddIssue(report, 'error', ...
                    'invalid_metadata_struct', nodeType, relativePath, ...
                    'Metadata must be a scalar struct.', false);
                return
            end

            missingFields = setdiff(requiredFields, fieldnames(metadata));
            if ~isempty(missingFields)
                report = obj.iAddIssue(report, 'error', ...
                    'missing_metadata_fields', nodeType, relativePath, ...
                    sprintf('Missing metadata field(s): %s', ...
                    strjoin(missingFields, ', ')), false);
            end
        end

        function report = iAddIssue(~, report, severity, code, nodeType, relativePath, message, repairable)
            %IADDISSUE Append one structured validation issue.

            issue = struct();
            issue.code = code;
            issue.severity = severity;
            issue.nodeType = nodeType;
            issue.relativePath = relativePath;
            issue.message = message;
            issue.repairable = repairable;
            issue.storeLevel = false; % set by iFinishReport

            switch severity
                case 'error'
                    report.errors(end+1) = issue;
                case 'warning'
                    report.warnings(end+1) = issue;
            end
        end

        function report = iNewValidationReport(~, mode)
            %INEWVALIDATIONREPORT Create an empty validation report.

            issueTemplate = struct( ...
                'code', '', ...
                'severity', '', ...
                'nodeType', '', ...
                'relativePath', '', ...
                'message', '', ...
                'repairable', false, ...
                'storeLevel', false);

            report = struct();
            report.isValid = true;
            report.mode = mode;
            report.errors = repmat(issueTemplate, 0, 1);
            report.warnings = repmat(issueTemplate, 0, 1);
            report.unregisteredItems = {};
            report.checkedOn = NaT;
        end

        function report = iFinishReport(obj, report, ProjectInfo, isScoped)
            %IFINISHREPORT Classify issues, set isValid, keep unscoped reports.
            %
            %   An issue is entity-level (storeLevel false) when it is not a
            %   project node and its relativePath lies inside a registered
            %   subject's folder; every other issue is store-level.

            subjectRels = {};
            if isstruct(ProjectInfo) && isfield(ProjectInfo, 'subjectRegistry') && ...
                    obj.iIsRegistryStruct(ProjectInfo.subjectRegistry)
                subjectRels = {ProjectInfo.subjectRegistry.relativePath};
            end
            for iIssue = 1:numel(report.errors)
                report.errors(iIssue).storeLevel = ...
                    UMITProjectStore.iIsStoreLevelIssue(report.errors(iIssue), subjectRels);
            end
            for iIssue = 1:numel(report.warnings)
                report.warnings(iIssue).storeLevel = ...
                    UMITProjectStore.iIsStoreLevelIssue(report.warnings(iIssue), subjectRels);
            end
            report.isValid = isempty(report.errors);
            report.checkedOn = datetime('now');
            if ~isScoped
                obj.LastValidationReport = report;
            end
        end

        function tf = iIsRegistryStruct(obj, registry)
            %IISREGISTRYSTRUCT Return true for a valid empty/non-empty registry.

            tf = isstruct(registry) && ...
                all(ismember(obj.Schema.registryFields, fieldnames(registry)));
        end

        function tf = iIsResourceRegistryStruct(obj, registry)
            %IISRESOURCEREGISTRYSTRUCT Return true for a resource registry.

            tf = isstruct(registry) && ...
                all(ismember(obj.Schema.resourceRecordFields, ...
                fieldnames(registry)));
        end

        function index = iFindRegistryIndex(~, registry, id)
            %IFINDREGISTRYINDEX Find a registry ID case-insensitively.

            index = [];
            if isempty(registry)
                return
            end
            index = find(strcmpi({registry.id}, char(string(id))), 1, 'first');
        end

        function index = iFindResourceIndex(~, registry, uuid)
            %IFINDRESOURCEINDEX Find a resource UUID.

            index = [];
            if isempty(registry)
                return
            end
            index = find(strcmpi({registry.uuid}, char(string(uuid))), 1, 'first');
        end

        function iAppendLog(obj, operation, targetUUID, result)
            %IAPPENDLOG Append a compact mutation record to the project log.

            logPath = fullfile(obj.ProjectRoot, obj.Schema.folders.internal, ...
                obj.Schema.folders.logs, 'project.log');
            fid = fopen(logPath, 'a');
            if fid < 0
                warning('Umitoolbox:UMITProjectStore:logWriteFailed', ...
                    'Could not open project log: %s', logPath);
                return
            end

            cleanupObj = onCleanup(@() fclose(fid));
            timestamp = char(datetime('now', ...
                'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            fprintf(fid, '%s\t%s\t%s\t%s\n', ...
                timestamp, operation, targetUUID, result);
            clear cleanupObj
        end
    end

    methods (Static, Access = private)
        function state = iBindingState(dataFolder)
            %IBINDINGSTATE 'none' (no ledger, or its binding was removed),
            %'bound', or 'invalid' (a ledger that fails validation).

            report = UMITProjectLedger.validate(dataFolder);
            if ~report.isValid
                state = 'invalid';
            elseif strcmp(report.info.state, 'Bound')
                state = 'bound';
            else
                state = 'none';
            end
        end

        function tf = iHasBinding(dataFolder)
            %IHASBINDING True for a bound folder or an invalid ledger (both
            %block a new binding, as an existing .umitlink did).

            tf = ~strcmp(UMITProjectStore.iBindingState(dataFolder), 'none');
        end

        function backup = iBackupLedger(dataFolder)
            %IBACKUPLEDGER In-memory copy of the folder's ledger file bytes.

            backup = struct('path', UMITProjectLedger.filePath(dataFolder), ...
                'existed', false, 'bytes', uint8([]));
            if isfile(backup.path)
                fid = fopen(backup.path, 'r');
                backup.bytes = fread(fid, Inf, '*uint8');
                fclose(fid);
                backup.existed = true;
            end
        end

        function iRestoreLedger(backup)
            %IRESTORELEDGER Undo ledger writes of a failed store operation.
            %
            %   The failed operation never happened, so its events are
            %   removed by restoring the previous bytes (or the absence of
            %   the file). This is a rollback, not a history rewrite.

            if isempty(backup)
                return
            end
            if ~backup.existed
                UMITProjectStore.iDeleteFileIfPresent(backup.path);
                return
            end
            fid = fopen(backup.path, 'w');
            if fid == -1
                warning('Umitoolbox:UMITProjectStore:ledgerRestoreFailed', ...
                    'Could not restore the project ledger %s.', backup.path);
                return
            end
            fwrite(fid, backup.bytes, 'uint8');
            fclose(fid);
        end

        function iRemoveBinding(dataFolder, reason)
            %IREMOVEBINDING Record that the store removed a folder's binding.
            %
            %   Appends an unbind event (the ledger and its history stay).
            %   A folder that is not bound is left unchanged. An invalid
            %   ledger cannot take an event and is left untouched (warning).

            switch UMITProjectStore.iBindingState(dataFolder)
                case 'bound'
                    UMITProjectLedger.unbind(dataFolder, 'Reason', reason, ...
                        'Origin', 'UMITProjectStore');
                case 'invalid'
                    warning('Umitoolbox:UMITProjectStore:invalidLedgerKept', ...
                        ['The project ledger in "%s" is invalid and was left as it is; ' ...
                         'the folder is no longer registered to the session.'], dataFolder);
            end
        end

        function ProjectBinding = iBindingRecordFromLedger(binding, schema)
            %IBINDINGRECORDFROMLEDGER The store's binding record (identity
            %from the ledger; times from its events).

            ProjectBinding = struct( ...
                'version', schema.projectBinding.version, ...
                'bindingUUID', binding.bindingUUID, ...
                'projectUUID', binding.projectUUID, ...
                'subjectUUID', binding.subjectUUID, ...
                'sessionUUID', binding.sessionUUID, ...
                'folderRole', binding.folderRole, ...
                'createdOn', UMITLedgerCore.localTime(binding.boundOn), ...
                'modifiedOn', UMITLedgerCore.localTime(binding.modifiedOn));
        end

        function metadata = iAddMetadataDefaults(metadata, entityType)
            %IADDMETADATADEFAULTS Add optional metadata fields missing from old files.

            switch entityType
                case 'project'
                    defaults = struct( ...
                        'institution', '', 'lab', '', ...
                        'principalInvestigator', '', 'experimenters', {{}}, ...
                        'projectStartDate', NaT, 'projectEndDate', NaT, ...
                        'notes', '');
                case 'subject'
                    defaults = struct( ...
                        'species', '', 'strain', '', 'genotype', '', 'sex', '', ...
                        'dateOfBirth', NaT, 'ageAtProjectEntry', '', ...
                        'ageReference', '', 'weight_g', NaN, 'housing', '', ...
                        'diet', '', 'healthStatus', '', 'earTag', '', 'rfid', '', ...
                        'notes', '');
                case 'session'
                    defaults = struct( ...
                        'behavioralTask', '', 'stimulusNotes', '', ...
                        'pharmacology', '', 'surgery', '', 'virus', '', ...
                        'anesthesia', '', 'notes', '');
                otherwise
                    return
            end

            fields = fieldnames(defaults);
            for iField = 1:numel(fields)
                fieldName = fields{iField};
                if ~isfield(metadata, fieldName)
                    metadata.(fieldName) = defaults.(fieldName);
                end
            end
        end

        function iAssertKnownInputFields(inputInfo, allowedFields, errID)
            %IASSERTKNOWNINPUTFIELDS Reject creation fields that would be ignored.
            %   (audit item 24: addSubject/addSession used to drop them
            %   silently, while the update APIs rejected them.)

            unknownFields = setdiff(fieldnames(inputInfo), allowedFields);
            if ~isempty(unknownFields)
                error('Umitoolbox:UMITProjectStore:unsupportedField', ...
                    '%s: unsupported field(s) %s. Supported fields: %s.', ...
                    errID, strjoin(unknownFields, ', '), strjoin(unique(allowedFields), ', '));
            end
        end

        function metadata = iApplyMetadataInput(metadata, inputInfo, allowedFields, errID)
            %IAPPLYMETADATAINPUT Apply supplied editable creation fields only.

            suppliedFields = intersect(fieldnames(inputInfo), allowedFields);
            for iField = 1:numel(suppliedFields)
                fieldName = suppliedFields{iField};
                metadata.(fieldName) = UMITProjectStore.iNormalizeMetadataValue( ...
                    fieldName, inputInfo.(fieldName), errID);
            end
        end

        function value = iNormalizeMetadataValue(fieldName, value, errID)
            %INORMALIZEMETADATAVALUE Validate values accepted by update APIs.

            dateFields = {'acquisitionDate', 'projectStartDate', ...
                'projectEndDate', 'dateOfBirth'};
            textListFields = {'experimenters'};
            numericFields = {'weight_g'};

            if ismember(fieldName, dateFields)
                if ~(isdatetime(value) && isscalar(value))
                    error(errID, 'Field "%s" must be a datetime scalar.', fieldName);
                end
            elseif ismember(fieldName, textListFields)
                try
                    value = UMITProjectStore.iNormalizeTextList(value, fieldName);
                catch ME
                    error(errID, 'Field "%s" must be text or a text list: %s', ...
                        fieldName, ME.message);
                end
            elseif ismember(fieldName, numericFields)
                if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
                        (isfinite(value) || isnan(value)))
                    error(errID, ...
                        'Field "%s" must be a real scalar numeric value.', fieldName);
                end
                value = double(value);
            else
                if ~(ischar(value) || (isstring(value) && isscalar(value)))
                    error(errID, 'Field "%s" must be a text scalar.', fieldName);
                end
                value = char(string(value));
            end
        end

        function projectUUID = iNormalizeProjectUUID(projectUUID)
            %INORMALIZEPROJECTUUID Validate a project UUID input.

            if ~(ischar(projectUUID) || ...
                    (isstring(projectUUID) && isscalar(projectUUID)))
                error('Umitoolbox:UMITProjectStore:invalidProjectUUID', ...
                    ['Project UUID must be a character vector or ' ...
                     'string scalar.']);
            end

            projectUUID = lower(strtrim(char(string(projectUUID))));
            schema = getUMITProjectSchema();

            if isempty(regexp(projectUUID, ...
                    schema.namingRules.uuidPattern, 'once'))
                error('Umitoolbox:UMITProjectStore:invalidProjectUUID', ...
                    'Invalid project UUID: %s', projectUUID);
            end
        end

        function folderName = iMakeProjectFolderName( ...
                projectName, projectUUID)
            %IMAKEPROJECTFOLDERNAME Build an immutable project directory name.

            projectName = strtrim(char(string(projectName)));
            safeName = matlab.lang.makeValidName( ...
                projectName, 'ReplacementStyle', 'underscore');
            safeName = regexprep(safeName, '_+', '_');
            safeName = regexprep(safeName, '^_+|_+$', '');

            if isempty(safeName)
                safeName = 'Project';
            end

            maxNameLength = 48;
            if numel(safeName) > maxNameLength
                safeName = safeName(1:maxNameLength);
                safeName = regexprep(safeName, '_+$', '');
            end

            folderName = [safeName, '__', projectUUID];
        end

        function registry = iEmptyRegistry()
            %IEMPTYREGISTRY Create an empty canonical node registry.

            template = struct( ...
                'uuid', '', ...
                'id', '', ...
                'displayName', '', ...
                'relativePath', '');
            registry = repmat(template, 0, 1);
        end

        function registry = iEmptyResourceRegistry()
            %IEMPTYRESOURCEREGISTRY Create an empty resource registry.

            template = struct( ...
                'uuid', '', ...
                'type', '', ...
                'fileName', '', ...
                'relativePath', '', ...
                'displayName', '', ...
                'description', '', ...
                'createdOn', NaT, ...
                'modifiedOn', NaT, ...
                'status', '', ...
                'archivedOn', NaT, ...
                'checksum', '', ...
                'sourceFile', '');
            registry = repmat(template, 0, 1);
        end


        function registry = iEmptyResolvedResourceRegistry()
            %IEMPTYRESOLVEDRESOURCEREGISTRY Create an empty query result.

            template = UMITProjectStore.iNewResourceRecord( ...
                '', '', '', '', '', '', NaT, '', '', '');
            template.ownerType = '';
            template.ownerUUID = '';
            template.absolutePath = '';
            template.fileExists = false;
            registry = repmat(template, 0, 1);
        end

        function record = iNewRegistryRecord(uuid, id, displayName, relativePath)
            %INEWREGISTRYRECORD Create one canonical node registry record.

            record = struct( ...
                'uuid', uuid, ...
                'id', id, ...
                'displayName', displayName, ...
                'relativePath', relativePath);
        end

        function record = iNewResourceRecord(uuid, type, fileName, relativePath, displayName, description, createdOn, status, checksum, sourceFile)
            %INEWRESOURCERECORD Create one canonical resource record.

            record = struct( ...
                'uuid', uuid, ...
                'type', type, ...
                'fileName', fileName, ...
                'relativePath', relativePath, ...
                'displayName', displayName, ...
                'description', description, ...
                'createdOn', createdOn, ...
                'modifiedOn', createdOn, ...
                'status', status, ...
                'archivedOn', NaT, ...
                'checksum', checksum, ...
                'sourceFile', sourceFile);
        end

        function location = iNewResourceLocation(ownerType, metadata, recordIndex, ownerPath, ownerBaseRelativePath, metadataPath, variableName)
            %INEWRESOURCELOCATION Create a resolved resource-location record.

            location = struct();
            location.ownerType = ownerType;
            location.metadata = metadata;
            location.recordIndex = recordIndex;
            location.ownerPath = ownerPath;
            location.ownerBaseRelativePath = ownerBaseRelativePath;
            location.metadataPath = metadataPath;
            location.variableName = variableName;
        end

        function value = iGetTextField(info, fieldName, defaultValue, allowEmpty, errID)
            %IGETTEXTFIELD Read and normalize one optional text field.

            if ~isfield(info, fieldName)
                value = defaultValue;
                return
            end

            raw = info.(fieldName);
            if ~(ischar(raw) || (isstring(raw) && isscalar(raw)))
                error(errID, 'Field "%s" must be a text scalar.', fieldName);
            end

            value = char(string(raw));
            if ~allowEmpty && isempty(value)
                error(errID, 'Field "%s" cannot be empty.', fieldName);
            end
        end

        function value = iGetDateField(info, fieldName, defaultValue, errID)
            %IGETDATEFIELD Read one optional datetime scalar field.

            if ~isfield(info, fieldName)
                value = defaultValue;
                return
            end

            value = info.(fieldName);
            if ~(isdatetime(value) && isscalar(value))
                error(errID, 'Field "%s" must be a datetime scalar.', fieldName);
            end
        end


        function uuid = iNormalizeUUIDInput(uuidIn)
            %INORMALIZEUUIDINPUT Normalize one resource UUID input.

            if ~(ischar(uuidIn) || ...
                    (isstring(uuidIn) && isscalar(uuidIn)))
                error('Umitoolbox:UMITProjectStore:invalidUUID', ...
                    'Resource UUID must be a character vector or string scalar.');
            end

            uuid = char(string(uuidIn));
            if isempty(uuid)
                error('Umitoolbox:UMITProjectStore:invalidUUID', ...
                    'Resource UUID cannot be empty.');
            end
        end

        function values = iNormalizeTextList(value, fieldName)
            %INORMALIZETEXTLIST Normalize a text scalar or vector to cellstr.

            if isempty(value)
                values = {};
                return
            end

            if ischar(value)
                values = {value};
            elseif isstring(value) && isvector(value)
                values = cellstr(value(:));
            elseif iscell(value) && isvector(value) && ...
                    all(cellfun(@(x) ischar(x) || ...
                    (isstring(x) && isscalar(x)), value))
                values = cellstr(string(value(:)));
            else
                error('Umitoolbox:UMITProjectStore:invalidResourceFilter', ...
                    '"%s" must be a text scalar or text vector.', fieldName);
            end

            if any(cellfun(@isempty, values))
                error('Umitoolbox:UMITProjectStore:invalidResourceFilter', ...
                    '"%s" cannot contain empty values.', fieldName);
            end
        end

        function pathKey = iNormalizeComparisonPath(pathIn)
            %INORMALIZECOMPARISONPATH Normalize a path for equality testing.

            pathKey = UMITProjectStore.iAbsolutePath(pathIn);
            if ispc || ismac
                pathKey = lower(pathKey);
            end
        end

        function uuid = iGenerateUUID()
            %IGENERATEUUID Generate a lowercase UUID string.

            uuid = lower(char(java.util.UUID.randomUUID()));
        end

        function path = iAbsolutePath(pathIn)
            %IABSOLUTEPATH Return a canonical absolute filesystem path.

            if ~(ischar(pathIn) || (isstring(pathIn) && isscalar(pathIn)))
                error('Umitoolbox:UMITProjectStore:invalidPath', ...
                    'Path must be a character vector or string scalar.');
            end

            pathIn = char(string(pathIn));
            if isempty(pathIn)
                error('Umitoolbox:UMITProjectStore:invalidPath', ...
                    'Path cannot be empty.');
            end

            path = char(java.io.File(pathIn).getCanonicalPath());
        end

        function matchingRoots = iFindProjectRoots(projectUUID)
            %IFINDPROJECTROOTS Project roots whose project.mat has PROJECTUUID.
            %   Never throws for a missing project (empty result). Internal
            %   staging folders (names starting with '.') are skipped.

            projectUUID = UMITProjectStore.iNormalizeProjectUUID(projectUUID);
            projectsRoot = UMITProjectStore.getProjectsRoot();
            schema = getUMITProjectSchema();

            folderList = dir(projectsRoot);
            folderList = folderList([folderList.isdir]);
            folderNames = {folderList.name};
            folderList = folderList(~startsWith(folderNames, '.'));

            matchingRoots = strings(0, 1);
            % New projects use a full UUID suffix. Search those candidates
            % first, then fall back to all project folders to tolerate a
            % manually renamed directory inside the static root.
            suffix = ['__', projectUUID];
            candidateMask = endsWith( ...
                {folderList.name}, suffix, 'IgnoreCase', true);
            candidateOrder = [find(candidateMask), find(~candidateMask)];

            for iCandidate = candidateOrder
                candidateRoot = fullfile( ...
                    projectsRoot, folderList(iCandidate).name);
                projectFile = fullfile( ...
                    candidateRoot, schema.files.projectMetadata);
                if ~isfile(projectFile)
                    continue
                end
                try
                    loaded = load(projectFile, ...
                        schema.metadataVariables.project, '-mat');
                catch
                    continue
                end
                if ~isfield(loaded, schema.metadataVariables.project)
                    continue
                end
                ProjectInfo = loaded.(schema.metadataVariables.project);
                if isstruct(ProjectInfo) && isscalar(ProjectInfo) && ...
                        isfield(ProjectInfo, 'projectUUID') && ...
                        strcmpi(ProjectInfo.projectUUID, projectUUID)
                    matchingRoots(end+1, 1) = string(candidateRoot); %#ok<AGROW>
                end
            end
            matchingRoots = unique(matchingRoots, 'stable');
        end

        function scope = iScope(varargin)
            %ISCOPE Validation scope of a mutation (see validate 'Scope').
            %
            %   scope = UMITProjectStore.iScope('subjects', uuids, ...
            %       'sessions', uuids, 'subjectTrees', uuids)
            %   Omitted lists are empty; iScope() is project-level only.

            p = inputParser;
            addParameter(p, 'subjects', {});
            addParameter(p, 'sessions', {});
            addParameter(p, 'subjectTrees', {});
            parse(p, varargin{:});
            normalize = @(values) lower(cellstr(string(values(:)')));
            scope = struct( ...
                'subjects', {normalize(p.Results.subjects)}, ...
                'sessions', {normalize(p.Results.sessions)}, ...
                'subjectTrees', {normalize(p.Results.subjectTrees)});
            scope.subjects = unique([scope.subjects, scope.subjectTrees], 'stable');
        end

        function scope = iLocationScope(location)
            %ILOCATIONSCOPE Validation scope of a resource's owner.
            %   LOCATION is an iNewResourceLocation record or an
            %   iResolveOwner result (fields type/ownerType and metadata).

            if isfield(location, 'ownerType')
                ownerType = location.ownerType;
            else
                ownerType = location.type;
            end
            metadata = location.metadata;
            if strcmp(ownerType, 'session')
                scope = UMITProjectStore.iScope('subjects', {metadata.subjectUUID}, ...
                    'sessions', {metadata.uuid});
            else
                scope = UMITProjectStore.iScope('subjects', {metadata.uuid});
            end
        end

        function tf = iScopeHasSubject(scope, subjectUUID)
            tf = ~isstruct(scope) || ismember(lower(char(string(subjectUUID))), scope.subjects);
        end

        function tf = iScopeHasSession(scope, subjectUUID, sessionUUID)
            tf = ~isstruct(scope) || ...
                ismember(lower(char(string(sessionUUID))), scope.sessions) || ...
                ismember(lower(char(string(subjectUUID))), scope.subjectTrees);
        end

        function tf = iIsStoreLevelIssue(issue, subjectRels)
            %IISSTORELEVELISSUE True unless the issue lies inside a registered
            %subject's folder (and is not a project-node issue).

            tf = true;
            if strcmp(issue.nodeType, 'project')
                return
            end
            rel = strrep(char(string(issue.relativePath)), '\', '/');
            for iRel = 1:numel(subjectRels)
                subjectRel = strrep(char(string(subjectRels{iRel})), '\', '/');
                if ~isempty(subjectRel) && (strcmp(rel, subjectRel) || ...
                        startsWith(rel, [subjectRel '/']))
                    tf = false;
                    return
                end
            end
        end

        function rel = iJoinRelative(varargin)
            %IJOINRELATIVE Join path parts using canonical forward slashes.

            parts = cell(size(varargin));
            for iPart = 1:numel(varargin)
                value = strrep(char(string(varargin{iPart})), '\', '/');
                value = regexprep(value, '^/+|/+$', '');
                parts{iPart} = value;
            end
            rel = strjoin(parts, '/');
        end

        function iRemoveFolderIfPresent(path)
            %IREMOVEFOLDERIFPRESENT Recursively remove a folder when it exists.

            if isfolder(path)
                rmdir(path, 's');
            end
        end

        function value = iGetEnvironmentValue(names, defaultValue)
            %IGETENVIRONMENTVALUE Return the first populated environment value.

            value = defaultValue;
            for iName = 1:numel(names)
                candidate = getenv(names{iName});
                if ~isempty(candidate)
                    value = candidate;
                    return
                end
            end
        end

        function messages = iJoinIssueMessages(issues)
            %IJOINISSUEMESSAGES Combine validation issue messages compactly.

            if isempty(issues)
                messages = 'Unknown validation failure.';
                return
            end

            count = min(numel(issues), 5);
            messages = strjoin({issues(1:count).message}, ' | ');
            if numel(issues) > count
                messages = sprintf('%s | ... (%d additional errors)', ...
                    messages, numel(issues) - count);
            end
        end

        function iValidateProjectBindingStruct( ...
                ProjectBinding, schema, errID)
            %IVALIDATEPROJECTBINDINGSTRUCT Strictly validate the store's binding record.

            if ~isstruct(ProjectBinding) || ~isscalar(ProjectBinding)
                error(errID, ...
                    'ProjectBinding must be a scalar structure.');
            end

            requiredFields = schema.requiredProjectBindingFields;
            actualFields = fieldnames(ProjectBinding);
            missingFields = setdiff(requiredFields, actualFields);
            unknownFields = setdiff(actualFields, requiredFields);

            if ~isempty(missingFields)
                error(errID, ...
                    'ProjectBinding is missing field(s): %s', ...
                    strjoin(missingFields, ', '));
            end
            if ~isempty(unknownFields)
                error(errID, ...
                    'ProjectBinding contains unsupported field(s): %s', ...
                    strjoin(unknownFields, ', '));
            end

            if ~isnumeric(ProjectBinding.version) || ...
                    ~isscalar(ProjectBinding.version) || ...
                    ProjectBinding.version ~= ...
                    schema.projectBinding.version
                error(errID, ...
                    'Unsupported ProjectBinding version.');
            end

            uuidFields = {'bindingUUID', 'projectUUID', ...
                'subjectUUID', 'sessionUUID'};
            for iField = 1:numel(uuidFields)
                fieldName = uuidFields{iField};
                value = ProjectBinding.(fieldName);
                if ~(ischar(value) || ...
                        (isstring(value) && isscalar(value))) || ...
                        isempty(regexp(char(string(value)), ...
                        schema.namingRules.uuidPattern, 'once'))
                    error(errID, ...
                        'ProjectBinding.%s must be a valid lowercase UUID.', ...
                        fieldName);
                end
            end

            folderRole = char(string(ProjectBinding.folderRole));
            if ~ismember(folderRole, ...
                    schema.projectBinding.folderRoles)
                error(errID, ...
                    'Unsupported ProjectBinding.folderRole: %s', ...
                    folderRole);
            end

            dateFields = {'createdOn', 'modifiedOn'};
            for iField = 1:numel(dateFields)
                value = ProjectBinding.(dateFields{iField});
                if ~isdatetime(value) || ~isscalar(value) || isnat(value)
                    error(errID, ...
                        'ProjectBinding.%s must be a valid datetime scalar.', ...
                        dateFields{iField});
                end
            end
        end

        function iDeleteFileIfPresent(filePath)
            %IDELETEFILEIFPRESENT Delete one file during rollback/cleanup.

            if isfile(filePath)
                delete(filePath);
            end
        end

    end
end
