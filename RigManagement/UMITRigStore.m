classdef UMITRigStore < handle
    %UMITRIGSTORE Manage a centralized UMIT rig metadata folder.
    %
    %   store = UMITRigStore.create(rigInfo)
    %   store = UMITRigStore.open(rigUUID)
    %   store = UMITRigStore.openByRigID(rigID)
    %   store = UMITRigStore.getOrCreateDefaultRig()
    %   rigs  = UMITRigStore.listRigs()
    %
    %   UMITRigStore is the only supported write interface for centralized
    %   rig metadata. A rig describes one physical imaging setup (cameras,
    %   filters, dual-camera coregistration) and is completely
    %   independent of any project -- projects/sessions reference a rig by
    %   rigUUID/rigID, they never own its folder or metadata. All rigs are
    %   stored below getUmitFolder('rigs'). It creates canonical folders,
    %   writes metadata atomically, imports/archives/restores managed
    %   resources, performs transactional ID renames, and validates the rig.
    %
    %   Main operations:
    %       getRigsRoot
    %       create
    %       open
    %       openByRigID
    %       listRigs
    %       rigExists
    %       getOrCreateDefaultRig
    %       getRigInfo
    %       updateRigMetadata
    %       renameRigID
    %       archiveRig
    %       addCameraCoregistration
    %       addCalibrationFile
    %       archiveResource
    %       restoreResource
    %       purgeResource
    %       setActiveCameraCoregistration
    %       clearActiveCameraCoregistration
    %       setActiveCalibrationFile
    %       getActiveCameraCoregistration
    %       getActiveCalibrationFile
    %       updateResourceMetadata
    %       listResources
    %       getResource
    %       resolveResourcePath
    %       getLockInfo
    %       clearStaleLock
    %       validate
    %
    %   Design rules:
    %       - UUIDs are immutable identities.
    %       - Rig IDs are validated filesystem names and double as the rig's
    %         folder name (uniqueness enforced by folder existence, since
    %         there is no parent registry -- rigs are themselves top-level).
    %       - Display names are editable without renaming physical folders.
    %       - Resource filenames are immutable after import.
    %       - Archived resources remain registered and restorable.
    %       - Resource files are never deleted directly by callers; purgeResource is the only
    %         sanctioned deletion path, and only for already-archived resources.
    %       - A rig itself is soft-deleted via archiveRig (status ->
    %         'archived'), never removed -- sessions that reference a
    %         rigID/rigUUID by value keep resolving via open/openByRigID/
    %         rigExists. Archived rigs are excluded from
    %         getOrCreateDefaultRig's candidate pool. The default rig
    %         cannot be archived directly.
    %       - RigInfo.metadata is an open struct reserved for future
    %         rig-scoped fields not yet designed (camera names, per-camera
    %         info, filter spectral profiles, ...). updateRigMetadata
    %         shallow-merges new keys into it, so adding a field later
    %         never requires migrating existing rig records.
    %       - Every mutation uses a rig lock and atomic metadata writes.
    %       - Invalid rigs open read-only.

    properties (SetAccess = private)
        RigRoot
        Schema
        IsReadOnly = false
        LastValidationReport
    end

    methods (Static)
        function rigsRoot = getRigsRoot()
            %GETRIGSROOT Return the single UMIT rigs root folder.

            errID = 'Umitoolbox:UMITRigStore:rigsRootFailed';

            if exist('getUmitFolder', 'file') ~= 2
                error(errID, ...
                    'Required function getUmitFolder.m was not found.');
            end

            try
                rigsRoot = getUmitFolder('rigs');
            catch ME
                error(errID, ...
                    'Could not resolve getUmitFolder(''rigs''): %s', ...
                    ME.message);
            end

            rigsRoot = UMITRigStore.iAbsolutePath(rigsRoot);

            if isfile(rigsRoot)
                error(errID, ...
                    'The UMIT rigs root resolves to a file: %s', rigsRoot);
            end

            if ~isfolder(rigsRoot)
                [ok, message] = mkdir(rigsRoot);
                if ~ok
                    error(errID, ...
                        'Could not create the UMIT rigs root: %s', message);
                end
            end
        end

        function obj = create(rigInfo)
            %CREATE Create a new rig under the static UMIT rigs root.
            %
            %   store = UMITRigStore.create(rigInfo)
            %
            %   Required rigInfo field:
            %       rigID
            %
            %   Optional fields:
            %       displayName
            %       description
            %       isDefault (default false)

            errID = 'Umitoolbox:UMITRigStore:createFailed';

            if ~isstruct(rigInfo) || ~isscalar(rigInfo)
                error(errID, '"rigInfo" must be a scalar struct.');
            end

            schema = getUMITRigSchema();
            rigID = UMITRigStore.iNormalizeManagedID( ...
                schema, UMITRigStore.iGetTextField(rigInfo, 'rigID', '', false, errID), ...
                'rigID');
            displayName = UMITRigStore.iGetTextField( ...
                rigInfo, 'displayName', rigID, true, errID);
            description = UMITRigStore.iGetTextField( ...
                rigInfo, 'description', '', true, errID);
            isDefault = false;
            if isfield(rigInfo, 'isDefault')
                if ~islogical(rigInfo.isDefault) || ~isscalar(rigInfo.isDefault)
                    error(errID, '"isDefault" must be a scalar logical.');
                end
                isDefault = rigInfo.isDefault;
            end

            rigsRoot = UMITRigStore.getRigsRoot();
            rigPath = fullfile(rigsRoot, rigID);

            if isfolder(rigPath) || isfile(rigPath)
                error(errID, 'Rig already exists: %s', rigID);
            end

            % Case-insensitive collision check against existing rig folders,
            % since rigID uniqueness is not case sensitive on Windows.
            existing = dir(rigsRoot);
            existing = existing([existing.isdir]);
            existingNames = {existing.name};
            existingNames = existingNames(~ismember(existingNames, {'.', '..'}));
            if any(strcmpi(existingNames, rigID))
                error(errID, 'Rig ID already exists: %s', rigID);
            end

            transactionPath = fullfile(tempname(rigsRoot));
            mkdir(transactionPath);
            cleanupStage = onCleanup(@() UMITRigStore.iRemoveFolderIfPresent(transactionPath));

            rigUUID = UMITRigStore.iGenerateUUID();
            stagedRigPath = fullfile(transactionPath, rigID);
            UMITRigStore.iCreateRigFolders(schema, stagedRigPath);

            nowTime = datetime('now');
            RigInfo = struct();
            RigInfo.schemaVersion = schema.version;
            RigInfo.uuid = rigUUID;
            RigInfo.rigID = rigID;
            RigInfo.displayName = displayName;
            RigInfo.description = description;
            RigInfo.createdOn = nowTime;
            RigInfo.modifiedOn = nowTime;
            RigInfo.isDefault = isDefault;
            RigInfo.status = 'active';
            RigInfo.archivedOn = NaT;
            RigInfo.metadata = struct();
            RigInfo.activeCoregistrationUUID = '';
            RigInfo.activeCalibrationFileUUID = '';
            RigInfo.resourceRegistry = UMITRigStore.iEmptyResourceRegistry();

            saveMatAtomic( ...
                fullfile(stagedRigPath, schema.files.rigMetadata), ...
                schema.metadataVariables.rig, RigInfo);

            [ok, message] = movefile(stagedRigPath, rigPath, 'f');
            if ~ok
                error(errID, 'Could not install rig folder: %s', message);
            end

            obj = UMITRigStore(rigPath, schema, false);
            clear cleanupStage
        end

        function obj = open(rigUUID)
            %OPEN Open an existing rig by immutable rig UUID.

            errID = 'Umitoolbox:UMITRigStore:openFailed';
            rigUUID = UMITRigStore.iNormalizeUUIDInput(rigUUID);

            rigsRoot = UMITRigStore.getRigsRoot();
            schema = getUMITRigSchema();

            entries = dir(rigsRoot);
            entries = entries([entries.isdir]);
            entryNames = {entries.name};
            entryNames = entryNames(~ismember(entryNames, {'.', '..'}));

            for iEntry = 1:numel(entryNames)
                candidateRoot = fullfile(rigsRoot, entryNames{iEntry});
                rigFile = fullfile(candidateRoot, schema.files.rigMetadata);

                if ~isfile(rigFile)
                    continue
                end

                try
                    loaded = load(rigFile, schema.metadataVariables.rig, '-mat');
                catch
                    continue
                end

                if ~isfield(loaded, schema.metadataVariables.rig)
                    continue
                end

                RigInfo = loaded.(schema.metadataVariables.rig);
                if isstruct(RigInfo) && isscalar(RigInfo) && ...
                        isfield(RigInfo, 'uuid') && strcmpi(RigInfo.uuid, rigUUID)
                    obj = UMITRigStore(candidateRoot, getUMITRigSchema(), false);
                    report = obj.validate('Mode', 'quick');
                    obj.LastValidationReport = report;
                    obj.IsReadOnly = ~report.isValid;
                    return
                end
            end

            error(errID, 'Rig UUID was not found: %s', rigUUID);
        end

        function obj = openByRigID(rigID)
            %OPENBYRIGID Open an existing rig by its human-readable rig ID.

            errID = 'Umitoolbox:UMITRigStore:openFailed';
            schema = getUMITRigSchema();
            rigID = UMITRigStore.iNormalizeManagedID(schema, rigID, 'rigID');

            rigsRoot = UMITRigStore.getRigsRoot();
            rigPath = fullfile(rigsRoot, rigID);

            if ~isfolder(rigPath)
                % Fall back to a case-insensitive folder match.
                entries = dir(rigsRoot);
                entries = entries([entries.isdir]);
                entryNames = {entries.name};
                idx = find(strcmpi(entryNames, rigID), 1, 'first');
                if isempty(idx)
                    error(errID, 'Rig ID was not found: %s', rigID);
                end
                rigPath = fullfile(rigsRoot, entryNames{idx});
            end

            rigFile = fullfile(rigPath, schema.files.rigMetadata);
            if ~isfile(rigFile)
                error(errID, 'Rig metadata file is missing: %s', rigFile);
            end

            loaded = load(rigFile, schema.metadataVariables.rig, '-mat');
            if ~isfield(loaded, schema.metadataVariables.rig)
                error(errID, 'Rig metadata variable is missing.');
            end

            obj = UMITRigStore(rigPath, getUMITRigSchema(), false);
            report = obj.validate('Mode', 'quick');
            obj.LastValidationReport = report;
            obj.IsReadOnly = ~report.isValid;
        end

        function rigs = listRigs()
            %LISTRIGS Enumerate all rigs under the static rigs root.
            %
            %   rigs = UMITRigStore.listRigs()
            %
            %   Returns a table with columns: RigUUID, RigID, DisplayName,
            %   IsDefault, Status, IsReadable, RigRoot. Status is
            %   'active'/'archived' for readable rigs, "" otherwise.

            schema = getUMITRigSchema();
            rigsRoot = UMITRigStore.getRigsRoot();

            entries = dir(rigsRoot);
            entries = entries([entries.isdir]);
            entryNames = {entries.name};
            entryNames = entryNames(~ismember(entryNames, {'.', '..'}));

            RigUUID = strings(0, 1);
            RigID = strings(0, 1);
            DisplayName = strings(0, 1);
            IsDefault = false(0, 1);
            Status = strings(0, 1);
            IsReadable = false(0, 1);
            RigRoot = strings(0, 1);

            for iEntry = 1:numel(entryNames)
                candidateRoot = fullfile(rigsRoot, entryNames{iEntry});
                rigFile = fullfile(candidateRoot, schema.files.rigMetadata);

                if ~isfile(rigFile)
                    continue
                end

                try
                    loaded = load(rigFile, schema.metadataVariables.rig, '-mat');
                    RigInfo = loaded.(schema.metadataVariables.rig);
                    RigInfo = UMITRigStore.iUpgradeRigInfo(RigInfo);
                    RigUUID(end+1, 1) = string(RigInfo.uuid); %#ok<AGROW>
                    RigID(end+1, 1) = string(RigInfo.rigID); %#ok<AGROW>
                    DisplayName(end+1, 1) = string(RigInfo.displayName); %#ok<AGROW>
                    IsDefault(end+1, 1) = logical(RigInfo.isDefault); %#ok<AGROW>
                    Status(end+1, 1) = string(RigInfo.status); %#ok<AGROW>
                    IsReadable(end+1, 1) = true; %#ok<AGROW>
                catch
                    RigUUID(end+1, 1) = ""; %#ok<AGROW>
                    RigID(end+1, 1) = string(entryNames{iEntry}); %#ok<AGROW>
                    DisplayName(end+1, 1) = ""; %#ok<AGROW>
                    IsDefault(end+1, 1) = false; %#ok<AGROW>
                    Status(end+1, 1) = ""; %#ok<AGROW>
                    IsReadable(end+1, 1) = false; %#ok<AGROW>
                end
                RigRoot(end+1, 1) = string(candidateRoot); %#ok<AGROW>
            end

            rigs = table(RigUUID, RigID, DisplayName, IsDefault, Status, IsReadable, RigRoot);
        end

        function tf = rigExists(rigID)
            %RIGEXISTS Return true if rigID resolves to a readable rig record.
            %
            %   tf = UMITRigStore.rigExists(rigID)
            %
            %   Accepts either a human-readable rigID (what DataParams
            %   stores at the session level, e.g.
            %   DataParams.cameraCoregistration.rigID) or a rig UUID.
            %   Resolves both active and archived rigs -- intended for
            %   validation callers (e.g. a GUI checking whether a
            %   session's stored rigID is still valid) that want a plain
            %   boolean instead of a try/catch around open/openByRigID.
            %   Never throws: any resolution failure (not found,
            %   unreadable metadata, bad input) returns false.

            tf = false;

            if nargin < 1 || ~(ischar(rigID) || (isstring(rigID) && isscalar(rigID)))
                return
            end

            rigID = strtrim(char(string(rigID)));
            if isempty(rigID)
                return
            end

            schema = getUMITRigSchema();
            isUUID = ~isempty(regexp(rigID, schema.namingRules.uuidPattern, 'once'));

            try
                if isUUID
                    UMITRigStore.open(rigID);
                else
                    UMITRigStore.openByRigID(rigID);
                end
                tf = true;
            catch
                tf = false;
            end
        end

        function obj = getOrCreateDefaultRig()
            %GETORCREATEDEFAULTRIG Resolve "the" rig for this machine.
            %
            %   obj = UMITRigStore.getOrCreateDefaultRig()
            %
            %   Intended for entry points (e.g. raw-data import) that need a
            %   rig but have no way to ask the user which one. Resolution
            %   order:
            %       1) A rig explicitly flagged isDefault.
            %       2) If exactly one rig exists, it is promoted to default
            %          (isDefault is persisted) and returned -- a single rig
            %          per machine is the common case.
            %       3) If no rig exists at all, one is created automatically,
            %          named after this machine, and flagged isDefault.
            %       4) If multiple rigs exist and none is flagged isDefault,
            %          this is genuinely ambiguous and throws -- mark one as
            %          default via updateRigMetadata, or open a specific rig
            %          by rigID/rigUUID instead.
            %
            %   Archived rigs are excluded from every resolution branch --
            %   an archived rig can still be opened explicitly (openByRigID/
            %   open/rigExists), but is never auto-selected for new work.

            errID = 'Umitoolbox:UMITRigStore:noDefaultRig';
            rigs = UMITRigStore.listRigs();
            rigs = rigs(rigs.Status ~= "archived", :);

            defaultIdx = find(rigs.IsDefault & rigs.IsReadable, 1, 'first');
            if ~isempty(defaultIdx)
                obj = UMITRigStore.open(char(rigs.RigUUID(defaultIdx)));
                return
            end

            readableRigs = rigs(rigs.IsReadable, :);

            if height(readableRigs) == 1
                obj = UMITRigStore.open(char(readableRigs.RigUUID(1)));
                obj.updateRigMetadata(struct('isDefault', true));
                return
            end

            if height(readableRigs) == 0
                rigID = UMITRigStore.iMakeDefaultRigID();
                obj = UMITRigStore.create(struct( ...
                    'rigID', rigID, ...
                    'displayName', rigID, ...
                    'isDefault', true));
                return
            end

            error(errID, ...
                ['Multiple rigs exist and none is marked as default. Mark one ' ...
                'as default (updateRigMetadata), or open a specific rig by ' ...
                'rigID/rigUUID instead of requesting the default.']);
        end
    end

    methods
        function obj = UMITRigStore(rigRoot, schema, isReadOnly)
            %UMITRIGSTORE Construct an object bound to one rig root.

            obj.RigRoot = UMITRigStore.iAbsolutePath(rigRoot);
            obj.Schema = schema;
            obj.IsReadOnly = isReadOnly;
            obj.LastValidationReport = obj.iNewValidationReport('none');
        end

        function RigInfo = getRigInfo(obj)
            %GETRIGINFO Return current rig metadata.

            RigInfo = obj.iLoadRigInfo();
        end

        function updateRigMetadata(obj, updates)
            %UPDATERIGMETADATA Update editable rig metadata fields.

            errID = 'Umitoolbox:UMITRigStore:updateRigFailed';
            obj.iAssertUpdateStruct(updates, errID);
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('updateRigMetadata'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            RigInfo = obj.iApplyEditableUpdates( ...
                RigInfo, updates, obj.Schema.editableFields.rig, errID);
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog('updateRigMetadata', RigInfo.uuid, 'completed');
        end

        function renameRigID(obj, newRigID)
            %RENAMERIGID Transactionally rename this rig's ID and folder.

            errID = 'Umitoolbox:UMITRigStore:renameRigFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('renameRigID'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            newRigID = obj.iNormalizeManagedID(obj.Schema, newRigID, 'rigID');

            if strcmp(RigInfo.rigID, newRigID)
                return
            end

            rigsRoot = fileparts(obj.RigRoot);
            newRigPath = fullfile(rigsRoot, newRigID);

            if ~strcmpi(RigInfo.rigID, newRigID) && ...
                    (isfolder(newRigPath) || isfile(newRigPath))
                error(errID, 'Rig ID already exists: %s', newRigID);
            end

            oldRigPath = obj.RigRoot;
            isCaseOnly = strcmpi(RigInfo.rigID, newRigID);

            if isCaseOnly
                % Case-only renames require staging through a temporary name
                % because some filesystems treat case-only renames as no-ops.
                stagingPath = [oldRigPath, '__renaming_', UMITRigStore.iGenerateUUID()];
                [ok, message] = movefile(oldRigPath, stagingPath, 'f');
                if ~ok
                    error(errID, 'Could not stage rig rename: %s', message);
                end
                [ok, message] = movefile(stagingPath, newRigPath, 'f');
                if ~ok
                    movefile(stagingPath, oldRigPath, 'f');
                    error(errID, 'Could not complete rig rename: %s', message);
                end
            else
                [ok, message] = movefile(oldRigPath, newRigPath, 'f');
                if ~ok
                    error(errID, 'Could not rename rig folder: %s', message);
                end
            end

            RigInfo.rigID = newRigID;
            RigInfo.modifiedOn = datetime('now');

            % The lock directory lives inside RigRoot and has just moved along
            % with the rest of the rig folder to newRigPath. Reassigning
            % lockCleanup fires the old onCleanup callback immediately (a
            % harmless no-op against oldRigPath, which no longer exists) and
            % installs a new one tracking the lock at its real, current
            % location, so protection stays continuous and release still
            % fires correctly when this function returns.
            newLockFolder = fullfile(newRigPath, ...
                obj.Schema.folders.internal, obj.Schema.folders.lock);
            lockCleanup = onCleanup(@() obj.iReleaseWriteLock(newLockFolder));

            try
                saveMatAtomic( ...
                    fullfile(newRigPath, obj.Schema.files.rigMetadata), ...
                    obj.Schema.metadataVariables.rig, RigInfo);
            catch ME
                movefile(newRigPath, oldRigPath, 'f');
                rethrow(ME)
            end

            obj.RigRoot = UMITRigStore.iAbsolutePath(newRigPath);
            obj.iAppendLog('renameRigID', RigInfo.uuid, 'completed');
        end

        function archiveRig(obj)
            %ARCHIVERIG Soft-delete this rig: mark it archived, keep its data.
            %
            %   store.archiveRig()
            %
            %   Nothing is deleted or moved. Sessions that reference this
            %   rig's rigID/rigUUID (e.g. DataParams.cameraCoregistration.
            %   rigID) keep resolving via openByRigID/open/rigExists.
            %   Archived rigs are excluded from getOrCreateDefaultRig's
            %   candidate pool, so they stop being selected for new work.
            %   The default rig cannot be archived directly -- mark a
            %   different rig as default first (updateRigMetadata).

            errID = 'Umitoolbox:UMITRigStore:archiveRigFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('archiveRig'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;

            if strcmp(RigInfo.status, 'archived')
                error(errID, 'Rig is already archived: %s', RigInfo.rigID);
            end

            if RigInfo.isDefault
                error(errID, ...
                    ['Cannot archive the default rig "%s". Mark a different ' ...
                    'rig as default first (updateRigMetadata).'], RigInfo.rigID);
            end

            RigInfo.status = 'archived';
            RigInfo.archivedOn = datetime('now');
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog('archiveRig', RigInfo.uuid, 'completed');
        end

        function resourceUUID = addCameraCoregistration(obj, sourceFile, resourceInfo)
            %ADDCAMERACOREGISTRATION Import a camera-coregistration transform.

            if nargin < 3
                resourceInfo = struct();
            end
            resourceUUID = obj.iAddManagedResource('cameraCoregistration', sourceFile, resourceInfo);
        end

        function resourceUUID = addCalibrationFile(obj, sourceFile, resourceInfo)
            %ADDCALIBRATIONFILE Import a general managed rig calibration file.

            if nargin < 3
                resourceInfo = struct();
            end
            resourceUUID = obj.iAddManagedResource('calibrationFile', sourceFile, resourceInfo);
        end

        function archiveResource(obj, resourceUUID, varargin)
            %ARCHIVERESOURCE Archive a managed resource without deleting it.
            %
            %   store.archiveResource(resourceUUID)
            %   store.archiveResource(resourceUUID, 'ReplacementUUID', replacementUUID)
            %
            %   An active resource requires a valid non-archived replacement.

            errID = 'Umitoolbox:UMITRigStore:archiveFailed';

            p = inputParser;
            p.FunctionName = 'archiveResource';
            addRequired(p, 'resourceUUID', @(x) ischar(x) || (isstring(x) && isscalar(x)));
            addParameter(p, 'ReplacementUUID', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
            parse(p, resourceUUID, varargin{:});

            resourceUUID = char(string(p.Results.resourceUUID));
            replacementUUID = char(string(p.Results.ReplacementUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('archiveResource'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            recordIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, resourceUUID);
            if isempty(recordIndex)
                error(errID, 'Resource is not owned by this rig: %s', resourceUUID);
            end
            record = RigInfo.resourceRegistry(recordIndex);

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
                    RigInfo.resourceRegistry, replacementUUID);
                if isempty(replacementIndex)
                    error(errID, 'Replacement resource is not owned by this rig.');
                end

                replacement = RigInfo.resourceRegistry(replacementIndex);
                if ~strcmp(replacement.type, record.type) || ...
                        strcmp(replacement.status, 'archived') || ...
                        strcmp(replacement.uuid, record.uuid)
                    error(errID, ...
                        'Replacement must be a different non-archived resource of the same type.');
                end

                RigInfo.resourceRegistry(replacementIndex).status = 'active';
                RigInfo.resourceRegistry(replacementIndex).modifiedOn = datetime('now');
                RigInfo.(pointerField) = replacement.uuid;
            end

            oldPath = obj.iResolveRelativePath(record.relativePath);
            archiveRel = obj.iBuildResourceRelativePath(record.type, obj.Schema.folders.archive, record.fileName);
            archivePath = obj.iResolveRelativePath(archiveRel);

            if isfile(archivePath)
                error(errID, 'Archive destination already exists: %s', archivePath);
            end

            [ok, message] = movefile(oldPath, archivePath, 'f');
            if ~ok
                error(errID, 'Could not archive resource: %s', message);
            end

            RigInfo.resourceRegistry(recordIndex).status = 'archived';
            RigInfo.resourceRegistry(recordIndex).archivedOn = datetime('now');
            RigInfo.resourceRegistry(recordIndex).modifiedOn = datetime('now');
            RigInfo.resourceRegistry(recordIndex).relativePath = archiveRel;
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                movefile(archivePath, oldPath, 'f');
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog('archiveResource', resourceUUID, 'completed');
        end

        function restoreResource(obj, resourceUUID)
            %RESTORERESOURCE Restore an archived resource as available.
            %
            %   Restoring does not automatically make the resource active.

            errID = 'Umitoolbox:UMITRigStore:restoreFailed';
            resourceUUID = char(string(resourceUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('restoreResource'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            recordIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, resourceUUID);
            if isempty(recordIndex)
                error(errID, 'Resource is not owned by this rig: %s', resourceUUID);
            end
            record = RigInfo.resourceRegistry(recordIndex);

            if ~strcmp(record.status, 'archived')
                error(errID, 'Resource is not archived: %s', resourceUUID);
            end

            oldPath = obj.iResolveRelativePath(record.relativePath);
            activeRel = obj.iBuildResourceRelativePath(record.type, obj.Schema.folders.active, record.fileName);
            activePath = obj.iResolveRelativePath(activeRel);

            if isfile(activePath)
                error(errID, 'Restore destination already exists: %s', activePath);
            end

            [ok, message] = movefile(oldPath, activePath, 'f');
            if ~ok
                error(errID, 'Could not restore resource: %s', message);
            end

            RigInfo.resourceRegistry(recordIndex).status = 'available';
            RigInfo.resourceRegistry(recordIndex).archivedOn = NaT;
            RigInfo.resourceRegistry(recordIndex).modifiedOn = datetime('now');
            RigInfo.resourceRegistry(recordIndex).relativePath = activeRel;
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                movefile(activePath, oldPath, 'f');
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog('restoreResource', resourceUUID, 'completed');
        end

        function purgeResource(obj, resourceUUID)
            %PURGERESOURCE Permanently delete an archived resource's file and registry entry.
            %
            %   store.purgeResource(resourceUUID)
            %
            %   This is the only sanctioned way to delete a managed resource file -- callers
            %   must never delete resource files directly. Only resources with status
            %   'archived' can be purged (archive first); this keeps deletion a deliberate
            %   two-step action rather than something that can happen to an active or
            %   available resource by accident.

            errID = 'Umitoolbox:UMITRigStore:purgeFailed';
            resourceUUID = char(string(resourceUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('purgeResource'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            recordIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, resourceUUID);
            if isempty(recordIndex)
                error(errID, 'Resource is not owned by this rig: %s', resourceUUID);
            end
            record = RigInfo.resourceRegistry(recordIndex);

            if ~strcmp(record.status, 'archived')
                error(errID, ...
                    'Only archived resources can be purged (status is "%s"): %s', ...
                    record.status, resourceUUID);
            end

            resourcePath = obj.iResolveRelativePath(record.relativePath);

            RigInfo.resourceRegistry(recordIndex) = [];
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            if isfile(resourcePath)
                delete(resourcePath);
            end

            obj.iAppendLog('purgeResource', resourceUUID, 'completed');
        end

        function setActiveCameraCoregistration(obj, resourceUUID)
            %SETACTIVECAMERACOREGISTRATION Select an active camera transform.

            obj.iSetActiveResource('cameraCoregistration', resourceUUID);
        end

        function clearActiveCameraCoregistration(obj)
            %CLEARACTIVECAMERACOREGISTRATION Temporarily leave no active transform.
            %
            %   The previously active camera coregistration remains available and
            %   can be reactivated with setActiveCameraCoregistration. This is used
            %   by recalibration workflows that must classify unregistered source
            %   data without deleting or archiving resource history.

            obj.iClearActiveResource('cameraCoregistration');
        end

        function setActiveCalibrationFile(obj, resourceUUID)
            %SETACTIVECALIBRATIONFILE Select an active general calibration file.

            obj.iSetActiveResource('calibrationFile', resourceUUID);
        end

        function resource = getActiveCameraCoregistration(obj)
            %GETACTIVECAMERACOREGISTRATION Return the active camera transform.

            resource = obj.iGetActiveResource('cameraCoregistration');
        end

        function resource = getActiveCalibrationFile(obj)
            %GETACTIVECALIBRATIONFILE Return the active general calibration file.

            resource = obj.iGetActiveResource('calibrationFile');
        end

        function updateResourceMetadata(obj, resourceUUID, updates)
            %UPDATERESOURCEMETADATA Update a resource display name or description.

            errID = 'Umitoolbox:UMITRigStore:updateResourceFailed';
            obj.iAssertUpdateStruct(updates, errID);
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('updateResourceMetadata'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            recordIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, resourceUUID);
            if isempty(recordIndex)
                error(errID, 'Resource is not owned by this rig: %s', resourceUUID);
            end

            record = RigInfo.resourceRegistry(recordIndex);
            record = obj.iApplyEditableUpdates( ...
                record, updates, obj.Schema.editableFields.resource, errID);
            record.modifiedOn = datetime('now');
            RigInfo.resourceRegistry(recordIndex) = record;
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog('updateResourceMetadata', resourceUUID, 'completed');
        end

        function resources = listResources(obj, varargin)
            %LISTRESOURCES List this rig's managed resources.
            %
            %   resources = store.listResources()
            %   resources = store.listResources('Type', 'cameraCoregistration')
            %   resources = store.listResources('Status', 'active')

            p = inputParser;
            p.FunctionName = 'listResources';
            addParameter(p, 'Type', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
            addParameter(p, 'Status', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
            parse(p, varargin{:});

            RigInfo = obj.iLoadRigInfo();
            resources = RigInfo.resourceRegistry;

            typeFilter = char(string(p.Results.Type));
            if ~isempty(typeFilter)
                resources = resources(strcmp({resources.type}, typeFilter));
            end

            statusFilter = char(string(p.Results.Status));
            if ~isempty(statusFilter)
                resources = resources(strcmp({resources.status}, statusFilter));
            end
        end

        function resource = getResource(obj, resourceUUID)
            %GETRESOURCE Return one resource record by UUID.

            errID = 'Umitoolbox:UMITRigStore:getResourceFailed';
            RigInfo = obj.iLoadRigInfo();
            recordIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, resourceUUID);
            if isempty(recordIndex)
                error(errID, 'Resource is not owned by this rig: %s', resourceUUID);
            end
            resource = RigInfo.resourceRegistry(recordIndex);
        end

        function filePath = resolveResourcePath(obj, resourceUUID, varargin)
            %RESOLVERESOURCEPATH Return the absolute path of one resource.

            p = inputParser;
            p.FunctionName = 'resolveResourcePath';
            addParameter(p, 'MustExist', true, @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            record = obj.getResource(resourceUUID);
            filePath = obj.iResolveRelativePath(record.relativePath);

            if p.Results.MustExist && ~isfile(filePath)
                error('Umitoolbox:UMITRigStore:resourceFileMissing', ...
                    'Resource file does not exist on disk: %s', filePath);
            end
        end

        function LockInfo = getLockInfo(obj)
            %GETLOCKINFO Return metadata about the current rig lock, if any.

            LockInfo = struct();
            lockFolder = fullfile(obj.RigRoot, obj.Schema.folders.internal, obj.Schema.folders.lock);
            lockMetadataPath = fullfile(lockFolder, obj.Schema.files.lockMetadata);

            if ~isfile(lockMetadataPath)
                return
            end

            try
                loaded = load(lockMetadataPath, 'LockInfo', '-mat');
                if isfield(loaded, 'LockInfo')
                    LockInfo = loaded.LockInfo;
                end
            catch
            end
        end

        function clearStaleLock(obj, varargin)
            %CLEARSTALELOCK Forcibly remove the rig lock directory.
            %
            %   Use only when certain no other process holds the lock (e.g.
            %   after a crash). This does not validate that assumption.

            p = inputParser;
            p.FunctionName = 'clearStaleLock';
            addParameter(p, 'Force', false, @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            lockFolder = fullfile(obj.RigRoot, obj.Schema.folders.internal, obj.Schema.folders.lock);

            if ~isfolder(lockFolder)
                return
            end

            if ~p.Results.Force
                error('Umitoolbox:UMITRigStore:clearLockRequiresForce', ...
                    'Pass ''Force'', true to confirm clearing a rig lock.');
            end

            UMITRigStore.iRemoveFolderIfPresent(lockFolder);
        end

        function report = validate(obj, varargin)
            %VALIDATE Validate this rig's metadata and resource registry.
            %
            %   report = store.validate()
            %   report = store.validate('Mode', 'quick')
            %   report = store.validate('Mode', 'full')
            %
            %   'quick' checks structural presence only. 'full' additionally
            %   cross-checks resource files against the registry.

            p = inputParser;
            p.FunctionName = 'validate';
            addParameter(p, 'Mode', 'full', @(x) ismember(lower(string(x)), ["quick", "full"]));
            parse(p, varargin{:});
            mode = lower(char(string(p.Results.Mode)));

            report = obj.iNewValidationReport(mode);

            rigFile = fullfile(obj.RigRoot, obj.Schema.files.rigMetadata);
            if ~isfile(rigFile)
                report = obj.iAddIssue(report, 'error', 'missingRigMetadata', ...
                    sprintf('Rig metadata file is missing: %s', rigFile));
                report.isValid = isempty(report.errors);
                return
            end

            try
                RigInfo = obj.iLoadRigInfo();
            catch ME
                report = obj.iAddIssue(report, 'error', 'unreadableRigMetadata', ME.message);
                report.isValid = isempty(report.errors);
                return
            end

            missingFields = setdiff(obj.Schema.requiredRigFields, fieldnames(RigInfo));
            if ~isempty(missingFields)
                report = obj.iAddIssue(report, 'error', 'missingRigField', ...
                    sprintf('Missing rig field(s): %s', strjoin(missingFields, ', ')));
            end

            if isfield(RigInfo, 'status') && isfield(obj.Schema, 'rigStatuses') && ...
                    ~ismember(RigInfo.status, obj.Schema.rigStatuses)
                report = obj.iAddIssue(report, 'error', 'invalidRigStatus', ...
                    sprintf('Rig status "%s" is not one of the allowed values.', RigInfo.status));
            end

            if isfield(RigInfo, 'status') && isfield(RigInfo, 'isDefault') && ...
                    strcmp(RigInfo.status, 'archived') && RigInfo.isDefault
                report = obj.iAddIssue(report, 'error', 'archivedDefaultRig', ...
                    'Rig is archived but still flagged as the default rig.');
            end

            if strcmp(mode, 'full') && isfield(RigInfo, 'resourceRegistry')
                for iRecord = 1:numel(RigInfo.resourceRegistry)
                    record = RigInfo.resourceRegistry(iRecord);

                    if ~isfield(obj.Schema.resourceTypes, record.type)
                        report = obj.iAddIssue(report, 'error', 'unknownResourceType', ...
                            sprintf('Resource %s has unknown type "%s".', record.uuid, record.type));
                        continue
                    end

                    resourcePath = obj.iResolveRelativePath(record.relativePath);
                    if ~isfile(resourcePath)
                        report = obj.iAddIssue(report, 'error', 'missingResourceFile', ...
                            sprintf('Resource %s file is missing: %s', record.uuid, resourcePath));
                    end
                end

                for iType = fieldnames(obj.Schema.resourceTypes)'
                    resourceType = iType{1};
                    resourceDef = obj.Schema.resourceTypes.(resourceType);
                    pointerField = resourceDef.activePointerField;
                    activeUUID = RigInfo.(pointerField);

                    if isempty(activeUUID)
                        continue
                    end

                    idx = obj.iFindResourceIndex(RigInfo.resourceRegistry, activeUUID);
                    if isempty(idx)
                        report = obj.iAddIssue(report, 'error', 'danglingActivePointer', ...
                            sprintf('%s points to a resource not in the registry: %s', ...
                                pointerField, activeUUID));
                    elseif strcmp(RigInfo.resourceRegistry(idx).status, 'archived')
                        report = obj.iAddIssue(report, 'error', 'activePointsToArchived', ...
                            sprintf('%s points to an archived resource: %s', ...
                                pointerField, activeUUID));
                    end
                end
            end

            report.isValid = isempty(report.errors);
        end
    end

    methods (Access = private)
        function RigInfo = iLoadRigInfo(obj)
            %ILOADRIGINFO Load current rig metadata from disk.

            rigFile = fullfile(obj.RigRoot, obj.Schema.files.rigMetadata);
            loaded = load(rigFile, obj.Schema.metadataVariables.rig, '-mat');
            if ~isfield(loaded, obj.Schema.metadataVariables.rig)
                error('Umitoolbox:UMITRigStore:invalidRigMetadata', ...
                    'Rig metadata variable "%s" is missing.', obj.Schema.metadataVariables.rig);
            end
            RigInfo = loaded.(obj.Schema.metadataVariables.rig);
            RigInfo = UMITRigStore.iUpgradeRigInfo(RigInfo);
        end

        function iSaveRigInfo(obj, RigInfo)
            %ISAVERIGINFO Save rig metadata atomically.

            saveMatAtomic( ...
                fullfile(obj.RigRoot, obj.Schema.files.rigMetadata), ...
                obj.Schema.metadataVariables.rig, RigInfo);
        end

        function resourceUUID = iAddManagedResource(obj, resourceType, sourceFile, resourceInfo)
            %IADDMANAGEDRESOURCE Import one file into a managed resource folder.

            errID = 'Umitoolbox:UMITRigStore:addResourceFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock(['add_', resourceType]); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                error(errID, 'Unsupported resource type: %s', resourceType);
            end

            if ~isstruct(resourceInfo) || ~isscalar(resourceInfo)
                error(errID, '"resourceInfo" must be a scalar struct.');
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);

            sourceFile = UMITRigStore.iAbsolutePath(sourceFile);
            if ~isfile(sourceFile)
                error(errID, 'Source file does not exist: %s', sourceFile);
            end

            [~, ~, extension] = fileparts(sourceFile);
            if ~any(strcmpi(extension, resourceDef.allowedExtensions))
                error(errID, 'Resource type "%s" does not allow extension "%s".', ...
                    resourceType, extension);
            end

            try
                sourceProbe = load(sourceFile, '-mat'); %#ok<NASGU>
            catch ME
                error(errID, 'Source MAT file cannot be loaded: %s', ME.message);
            end

            displayName = UMITRigStore.iGetTextField( ...
                resourceInfo, 'displayName', resourceDef.filePrefix, true, errID);
            description = UMITRigStore.iGetTextField( ...
                resourceInfo, 'description', '', true, errID);

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;

            resourceUUID = UMITRigStore.iGenerateUUID();
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss'));
            compactUUID = strrep(resourceUUID, '-', '');
            fileName = sprintf('%s_%s_%s%s', ...
                resourceDef.filePrefix, timestamp, compactUUID(1:8), lower(extension));

            destinationRel = obj.iBuildResourceRelativePath( ...
                resourceType, obj.Schema.folders.active, fileName);
            destinationPath = obj.iResolveRelativePath(destinationRel);

            transactionPath = obj.iCreateTransactionFolder(['add_', resourceType]);
            stagedFile = fullfile(transactionPath, fileName);
            cleanupStage = onCleanup(@() UMITRigStore.iRemoveFolderIfPresent(transactionPath));

            [ok, message] = copyfile(sourceFile, stagedFile, 'f');
            if ~ok
                error(errID, 'Could not stage resource: %s', message);
            end

            checksum = computeFileChecksum(stagedFile);
            nowTime = datetime('now');
            pointerField = resourceDef.activePointerField;

            if isempty(RigInfo.(pointerField))
                status = 'active';
                RigInfo.(pointerField) = resourceUUID;
            else
                status = 'available';
            end

            record = UMITRigStore.iNewResourceRecord( ...
                resourceUUID, resourceType, fileName, destinationRel, ...
                displayName, description, nowTime, status, checksum, sourceFile);
            RigInfo.resourceRegistry(end+1) = record;
            RigInfo.modifiedOn = nowTime;

            [ok, message] = movefile(stagedFile, destinationPath, 'f');
            if ~ok
                error(errID, 'Could not install resource: %s', message);
            end

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                if isfile(destinationPath)
                    delete(destinationPath);
                end
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog(['add_', resourceType], resourceUUID, 'completed');
        end

        function iSetActiveResource(obj, resourceType, resourceUUID)
            %ISETACTIVERESOURCE Select one non-archived resource as active.

            errID = 'Umitoolbox:UMITRigStore:setActiveResourceFailed';
            resourceUUID = char(string(resourceUUID));

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock(['setActive_', resourceType]); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            recordIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, resourceUUID);

            if isempty(recordIndex)
                error(errID, 'Resource is not owned by this rig: %s', resourceUUID);
            end

            record = RigInfo.resourceRegistry(recordIndex);
            if ~strcmp(record.type, resourceType)
                error(errID, 'Resource "%s" is not of type "%s".', resourceUUID, resourceType);
            end
            if strcmp(record.status, 'archived')
                error(errID, 'Cannot activate an archived resource: %s', resourceUUID);
            end

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            pointerField = resourceDef.activePointerField;
            previousActiveUUID = RigInfo.(pointerField);

            if ~isempty(previousActiveUUID) && ~strcmp(previousActiveUUID, resourceUUID)
                prevIndex = obj.iFindResourceIndex(RigInfo.resourceRegistry, previousActiveUUID);
                if ~isempty(prevIndex)
                    RigInfo.resourceRegistry(prevIndex).status = 'available';
                    RigInfo.resourceRegistry(prevIndex).modifiedOn = datetime('now');
                end
            end

            RigInfo.resourceRegistry(recordIndex).status = 'active';
            RigInfo.resourceRegistry(recordIndex).modifiedOn = datetime('now');
            RigInfo.(pointerField) = resourceUUID;
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog(['setActive_', resourceType], resourceUUID, 'completed');
        end

        function iClearActiveResource(obj, resourceType)
            %ICLEARACTIVERESOURCE Demote the active resource without archiving it.

            errID = 'Umitoolbox:UMITRigStore:clearActiveResourceFailed';

            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock(['clearActive_', resourceType]); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                error(errID, 'Unsupported resource type: %s', resourceType);
            end

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            resourceDef = obj.Schema.resourceTypes.(resourceType);
            pointerField = resourceDef.activePointerField;
            previousActiveUUID = char(string(RigInfo.(pointerField)));

            if isempty(previousActiveUUID)
                return
            end

            recordIndex = obj.iFindResourceIndex( ...
                RigInfo.resourceRegistry, previousActiveUUID);
            if isempty(recordIndex) || ...
                    ~strcmp(RigInfo.resourceRegistry(recordIndex).type, resourceType) || ...
                    ~strcmp(RigInfo.resourceRegistry(recordIndex).status, 'active')
                error(errID, ...
                    'The active pointer for resource type "%s" is inconsistent.', ...
                    resourceType);
            end

            RigInfo.resourceRegistry(recordIndex).status = 'available';
            RigInfo.resourceRegistry(recordIndex).modifiedOn = datetime('now');
            RigInfo.(pointerField) = '';
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog(['clearActive_', resourceType], ...
                previousActiveUUID, 'completed');
        end

        function resource = iGetActiveResource(obj, resourceType)
            %IGETACTIVERESOURCE Resolve and validate the active resource for a type.

            errID = 'Umitoolbox:UMITRigStore:invalidActiveResource';

            if ~isfield(obj.Schema.resourceTypes, resourceType)
                error('Umitoolbox:UMITRigStore:invalidResourceType', ...
                    'Unsupported resource type: %s', resourceType);
            end

            RigInfo = obj.iLoadRigInfo();
            resourceDef = obj.Schema.resourceTypes.(resourceType);
            pointerField = resourceDef.activePointerField;

            if ~isfield(RigInfo, pointerField)
                error(errID, 'Rig metadata is missing active pointer "%s".', pointerField);
            end

            activeUUID = RigInfo.(pointerField);
            if ~(ischar(activeUUID) || (isstring(activeUUID) && isscalar(activeUUID)))
                error(errID, 'Active pointer "%s" must be a text scalar.', pointerField);
            end
            activeUUID = char(string(activeUUID));

            if isempty(activeUUID)
                resource = [];
                return
            end

            matchingIndices = find(strcmp({RigInfo.resourceRegistry.uuid}, activeUUID));
            if isempty(matchingIndices)
                error(errID, ...
                    'Active pointer "%s" does not resolve to a registered resource.', ...
                    pointerField);
            end
            if numel(matchingIndices) > 1
                error('Umitoolbox:UMITRigStore:duplicateResourceUUID', ...
                    'Active resource UUID is registered more than once: %s', activeUUID);
            end

            record = RigInfo.resourceRegistry(matchingIndices);
            if ~strcmp(record.type, resourceType) || ~strcmp(record.status, 'active')
                error(errID, ...
                    ['Active pointer "%s" does not reference an active resource ' ...
                    'of type "%s".'], pointerField, resourceType);
            end

            resourcePath = obj.iResolveRelativePath(record.relativePath);
            if ~isfile(resourcePath)
                error('Umitoolbox:UMITRigStore:missingResourceFile', ...
                    'Active resource file is missing: %s', resourcePath);
            end

            resource = record;
        end

        function index = iFindResourceIndex(~, registry, resourceUUID)
            %IFINDRESOURCEINDEX Return the registry index for one resource UUID.

            index = find(strcmp({registry.uuid}, resourceUUID), 1, 'first');
        end

        function rel = iBuildResourceRelativePath(obj, resourceType, stateFolder, fileName)
            %IBUILDRESOURCERELATIVEPATH Build a canonical resource path.

            resourceDef = obj.Schema.resourceTypes.(resourceType);
            parts = [resourceDef.relativeFolderParts, {stateFolder, fileName}];
            rel = UMITRigStore.iJoinRelative(parts{:});
        end

        function path = iResolveRelativePath(obj, relativePath)
            %IRESOLVERELATIVEPATH Resolve a rig-relative path.

            relativePath = char(string(relativePath));
            relativePath = strrep(relativePath, '/', filesep);
            path = fullfile(obj.RigRoot, relativePath);
        end

        function cleanupObj = iAcquireWriteLock(obj, operation)
            %IACQUIREWRITELOCK Acquire the exclusive rig mutation lock.

            lockFolder = fullfile(obj.RigRoot, obj.Schema.folders.internal, obj.Schema.folders.lock);
            lockMetadataPath = fullfile(lockFolder, obj.Schema.files.lockMetadata);

            lockDirectory = java.io.File(lockFolder);
            if ~lockDirectory.mkdir()
                if ~isfolder(lockFolder)
                    error('Umitoolbox:UMITRigStore:lockCreationFailed', ...
                        'Could not create rig lock directory: %s', lockFolder);
                end

                details = '';
                try
                    existing = obj.getLockInfo();
                    if isfield(existing, 'operation')
                        details = sprintf(' Existing operation: %s.', existing.operation);
                    end
                catch
                end

                error('Umitoolbox:UMITRigStore:rigLocked', ...
                    'Rig is locked by another operation.%s', details);
            end

            try
                LockInfo = struct();
                LockInfo.operation = operation;
                LockInfo.createdOn = datetime('now');
                LockInfo.processID = feature('getpid');
                LockInfo.userName = UMITRigStore.iGetEnvironmentValue({'USERNAME', 'USER'}, 'unknown');
                LockInfo.hostName = UMITRigStore.iGetEnvironmentValue({'COMPUTERNAME', 'HOSTNAME'}, 'unknown');

                saveMatAtomic(lockMetadataPath, 'LockInfo', LockInfo);
            catch ME
                UMITRigStore.iRemoveFolderIfPresent(lockFolder);
                rethrow(ME)
            end

            cleanupObj = onCleanup(@() obj.iReleaseWriteLock(lockFolder));
        end

        function iReleaseWriteLock(~, lockFolder)
            %IRELEASEWRITELOCK Remove the rig mutation lock directory.

            try
                UMITRigStore.iRemoveFolderIfPresent(lockFolder);
            catch ME
                warning('Umitoolbox:UMITRigStore:lockReleaseFailed', ...
                    'Rig operation completed, but the lock directory could not be removed: %s', ...
                    ME.message);
            end
        end

        function iAssertWritable(obj)
            %IASSERTWRITABLE Reject writes when the rig is read-only.

            if obj.IsReadOnly
                error('Umitoolbox:UMITRigStore:readOnlyRig', ...
                    ['Rig is open in read-only mode because validation failed. ' ...
                    'Repair the rig before modifying it.']);
            end
        end

        function iAssertHealthyForMutation(obj)
            %IASSERTHEALTHYFORMUTATION Revalidate before every mutation.

            report = obj.validate('Mode', 'full');
            if ~report.isValid
                obj.IsReadOnly = true;
                error('Umitoolbox:UMITRigStore:invalidRig', ...
                    'Rig validation failed: %s', UMITRigStore.iJoinIssueMessages(report.errors));
            end
        end

        function iAssertValidAfterMutation(obj)
            %IASSERTVALIDAFTERMUTATION Validate the rig after a write.

            report = obj.validate('Mode', 'full');
            if ~report.isValid
                error('Umitoolbox:UMITRigStore:postWriteValidationFailed', ...
                    'Rig failed post-write validation: %s', ...
                    UMITRigStore.iJoinIssueMessages(report.errors));
            end
        end

        function iAssertUpdateStruct(~, updates, errID)
            %IASSERTUPDATESTRUCT Validate an update payload.

            if ~isstruct(updates) || ~isscalar(updates)
                error(errID, '"updates" must be a scalar struct.');
            end
        end

        function target = iApplyEditableUpdates(~, target, updates, allowedFields, errID)
            %IAPPLYEDITABLEUPDATES Apply only explicitly allowed metadata fields.
            %
            %   The 'metadata' field (rig records only) is shallow-merged
            %   into the existing struct rather than replaced, so a caller
            %   adding one new key (e.g. cameraNames) never has to
            %   read-modify-write the whole bucket or clobber keys set by
            %   other callers.

            updateFields = fieldnames(updates);
            invalidFields = setdiff(updateFields, allowedFields);
            if ~isempty(invalidFields)
                error(errID, 'Unsupported editable field(s): %s', strjoin(invalidFields, ', '));
            end

            for iField = 1:numel(updateFields)
                fieldName = updateFields{iField};
                value = updates.(fieldName);

                if ismember(fieldName, {'displayName', 'description'})
                    if ~(ischar(value) || (isstring(value) && isscalar(value)))
                        error(errID, 'Field "%s" must be a text scalar.', fieldName);
                    end
                    value = char(string(value));
                elseif strcmp(fieldName, 'isDefault')
                    if ~islogical(value) || ~isscalar(value)
                        error(errID, 'Field "isDefault" must be a scalar logical.');
                    end
                elseif strcmp(fieldName, 'metadata')
                    if ~isstruct(value) || ~isscalar(value)
                        error(errID, 'Field "metadata" must be a scalar struct.');
                    end
                    existingMetadata = target.metadata;
                    if ~isstruct(existingMetadata) || ~isscalar(existingMetadata)
                        existingMetadata = struct();
                    end
                    newKeys = fieldnames(value);
                    for iKey = 1:numel(newKeys)
                        existingMetadata.(newKeys{iKey}) = value.(newKeys{iKey});
                    end
                    value = existingMetadata;
                end

                target.(fieldName) = value;
            end
        end

        function path = iCreateTransactionFolder(obj, operation)
            %ICREATETRANSACTIONFOLDER Create a unique internal staging folder.

            name = sprintf('%s_%s', operation, UMITRigStore.iGenerateUUID());
            path = fullfile(obj.RigRoot, obj.Schema.folders.internal, ...
                obj.Schema.folders.transactions, name);
            mkdir(path);
        end

        function iAppendLog(obj, operation, targetUUID, result)
            %IAPPENDLOG Append one operation-log entry.

            try
                logFolder = fullfile(obj.RigRoot, obj.Schema.folders.internal, obj.Schema.folders.logs);
                if ~isfolder(logFolder)
                    mkdir(logFolder);
                end

                logFile = fullfile(logFolder, 'operations.log');
                fid = fopen(logFile, 'a');
                if fid == -1
                    return
                end

                fprintf(fid, '%s\t%s\t%s\t%s\n', ...
                    char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
                    operation, targetUUID, result);
                fclose(fid);
            catch
            end
        end

        function report = iNewValidationReport(~, mode)
            %INEWVALIDATIONREPORT Create an empty validation report.

            report = struct();
            report.mode = mode;
            report.isValid = true;
            report.errors = struct('code', {}, 'message', {});
            report.checkedOn = datetime('now');
        end

        function report = iAddIssue(~, report, severity, code, message)
            %IADDISSUE Append one issue to a validation report.

            issue = struct('code', code, 'message', message);
            if strcmp(severity, 'error')
                report.errors(end+1) = issue;
            end
        end
    end

    methods (Static, Access = private)
        function RigInfo = iUpgradeRigInfo(RigInfo)
            %IUPGRADERIGINFO Backfill fields added by later schema versions.
            %
            %   Rig records written under an earlier schema version are
            %   missing fields added since (schema v2 added metadata/
            %   status/archivedOn). Rather than requiring a bulk migration
            %   pass, every load backfills defaults in memory here; the
            %   upgraded shape persists on disk automatically the next
            %   time anything writes this rig's metadata.

            if ~isfield(RigInfo, 'metadata') || ~isstruct(RigInfo.metadata) || ~isscalar(RigInfo.metadata)
                RigInfo.metadata = struct();
            end
            if ~isfield(RigInfo, 'status') || isempty(RigInfo.status)
                RigInfo.status = 'active';
            end
            if ~isfield(RigInfo, 'archivedOn')
                RigInfo.archivedOn = NaT;
            end
            if ~isfield(RigInfo, 'schemaVersion') || RigInfo.schemaVersion < 2
                RigInfo.schemaVersion = 2;
            end
        end

        function id = iNormalizeManagedID(schema, idIn, idLabel)
            %INORMALIZEMANAGEDID Validate one filesystem-backed managed ID.

            if ~(ischar(idIn) || (isstring(idIn) && isscalar(idIn)))
                error('Umitoolbox:UMITRigStore:invalidID', ...
                    '"%s" must be a character vector or string scalar.', idLabel);
            end

            id = char(string(idIn));
            rules = schema.namingRules;

            if isempty(id) || strlength(string(id)) > rules.maxIDLength
                error('Umitoolbox:UMITRigStore:invalidID', ...
                    '"%s" must contain 1 to %d characters.', idLabel, rules.maxIDLength);
            end

            if isempty(regexp(id, rules.idPattern, 'once'))
                error('Umitoolbox:UMITRigStore:invalidID', ...
                    ['"%s" must start with a letter or digit and contain ' ...
                    'only letters, digits, underscores, and hyphens.'], idLabel);
            end

            if any(strcmpi(id, rules.reservedNames))
                error('Umitoolbox:UMITRigStore:invalidID', ...
                    '"%s" uses a reserved filesystem name: %s', idLabel, id);
            end

            for iPrefix = 1:numel(rules.reservedPrefixes)
                if startsWith(id, rules.reservedPrefixes{iPrefix}, 'IgnoreCase', true)
                    error('Umitoolbox:UMITRigStore:invalidID', ...
                        '"%s" uses reserved prefix "%s".', idLabel, rules.reservedPrefixes{iPrefix});
                end
            end
        end

        function iCreateRigFolders(schema, rigPath)
            %ICREATERIGFOLDERS Create canonical rig calibration folders.

            mkdir(rigPath);

            internalBase = fullfile(rigPath, schema.folders.internal);
            mkdir(fullfile(internalBase, schema.folders.transactions));
            mkdir(fullfile(internalBase, schema.folders.recovery));
            mkdir(fullfile(internalBase, schema.folders.logs));

            coregBase = fullfile(rigPath, schema.folders.transforms, schema.folders.cameraCoregistration);
            mkdir(fullfile(coregBase, schema.folders.active));
            mkdir(fullfile(coregBase, schema.folders.archive));

            fileBase = fullfile(rigPath, schema.folders.calibrationFiles);
            mkdir(fullfile(fileBase, schema.folders.active));
            mkdir(fullfile(fileBase, schema.folders.archive));
        end

        function rigID = iMakeDefaultRigID()
            %IMAKEDEFAULTRIGID Build a rig ID for the auto-created default rig.

            hostName = UMITRigStore.iGetEnvironmentValue({'COMPUTERNAME', 'HOSTNAME'}, '');
            safeName = matlab.lang.makeValidName(hostName, 'ReplacementStyle', 'underscore');
            safeName = regexprep(safeName, '_+', '_');
            safeName = regexprep(safeName, '^_+|_+$', '');

            if isempty(safeName)
                rigID = 'DefaultRig';
            else
                rigID = ['Rig_', safeName];
            end
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

        function uuid = iGenerateUUID()
            %IGENERATEUUID Generate a lowercase UUID string.

            uuid = lower(char(java.util.UUID.randomUUID()));
        end

        function uuid = iNormalizeUUIDInput(uuidIn)
            %INORMALIZEUUIDINPUT Normalize one rig/resource UUID input.

            if ~(ischar(uuidIn) || (isstring(uuidIn) && isscalar(uuidIn)))
                error('Umitoolbox:UMITRigStore:invalidUUID', ...
                    'UUID must be a character vector or string scalar.');
            end

            uuid = char(string(uuidIn));
            if isempty(uuid)
                error('Umitoolbox:UMITRigStore:invalidUUID', 'UUID cannot be empty.');
            end
        end

        function path = iAbsolutePath(pathIn)
            %IABSOLUTEPATH Return a canonical absolute filesystem path.

            if ~(ischar(pathIn) || (isstring(pathIn) && isscalar(pathIn)))
                error('Umitoolbox:UMITRigStore:invalidPath', ...
                    'Path must be a character vector or string scalar.');
            end

            pathIn = char(string(pathIn));
            if isempty(pathIn)
                error('Umitoolbox:UMITRigStore:invalidPath', 'Path cannot be empty.');
            end

            path = char(java.io.File(pathIn).getCanonicalPath());
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
                messages = sprintf('%s | ... (%d additional errors)', messages, numel(issues) - count);
            end
        end
    end
end
