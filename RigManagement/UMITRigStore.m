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
    %       getDefaultRig
    %       setDefaultRig
    %       getOrCreateDefaultRig
    %       normalizeIlluminationName
    %       getSpectrum / listSpectra / importSpectrum / removeSpectrum
    %       getFilterSet / listFilterSets / importFilterSet
    %       getRigInfo
    %       duplicate
    %       setCameras
    %       setIlluminations
    %       resolveOpticalConfiguration
    %       updateRigMetadata
    %       renameRigID
    %       archiveRig
    %       restoreRig
    %       addCameraCoregistration
    %       importAndActivateCameraCoregistration
    %       archiveResource
    %       restoreResource
    %       purgeResource
    %       setActiveCameraCoregistration
    %       clearActiveCameraCoregistration
    %       getActiveCameraCoregistration
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
    %         getOrCreateDefaultRig's candidate pool. Default selection is
    %         represented once at store level, never as editable per-Rig
    %         metadata.
    %       - Cameras and canonical illumination definitions are concise,
    %         validated RigInfo fields. Their spectrum IDs reference the
    %         shared human-readable optical repertoire.
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

        function obj = create(rigInfo, varargin)
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
            %       metadata
            %       cameras
            %       illuminations

            errID = 'Umitoolbox:UMITRigStore:createFailed';
            p = inputParser;
            p.FunctionName = 'UMITRigStore.create';
            addParameter(p, 'InternalStoreLockHeld', false, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

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
            if isfield(rigInfo, 'metadata')
                metadata = rigInfo.metadata;
                if ~isstruct(metadata) || ~isscalar(metadata)
                    error(errID, 'Field "metadata" must be a scalar struct.');
                end
            else
                metadata = struct();
            end
            if isfield(rigInfo, 'isDefault')
                error(errID, ...
                    ['"isDefault" is no longer Rig metadata. Create the Rig, ' ...
                    'then call UMITRigStore.setDefaultRig(rigUUID).']);
            end
            cameras = UMITRigStore.iGetStructArrayField(rigInfo, 'cameras');
            illuminations = UMITRigStore.iGetStructArrayField(rigInfo, 'illuminations');
            cameras = UMITRigStore.iValidateCameraRecords(schema, cameras, true);
            illuminations = UMITRigStore.iValidateIlluminationRecords( ...
                schema, illuminations, true);

            if ~p.Results.InternalStoreLockHeld
                storeLock = UMITRigStore.iAcquireStoreLock('createRig'); %#ok<NASGU>
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
            RigInfo.status = 'active';
            RigInfo.archivedOn = NaT;
            RigInfo.metadata = metadata;
            RigInfo.cameras = cameras;
            RigInfo.illuminations = illuminations;
            RigInfo.activeCoregistrationUUID = '';
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
                    IsDefault(end+1, 1) = false; %#ok<AGROW>
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

            defaultUUID = UMITRigStore.iReadDefaultRigUUID();
            if ~isempty(defaultUUID)
                IsDefault = IsReadable & strcmpi(RigUUID, string(defaultUUID));
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

        function definition = getBuiltInRigDefinition(templateID)
            %GETBUILTINRIGDEFINITION Return a hydrated built-in Rig template.
            %
            %   definition = UMITRigStore.getBuiltInRigDefinition('OiS200')
            %
            %   This operation is read-only. It returns a complete rigInfo
            %   structure suitable for review or customization before it is
            %   passed to UMITRigStore.create. It never creates a Rig or
            %   changes the template's canonical Rig ID.

            if nargin < 1 || isempty(templateID)
                templateID = 'OiS200';
            end
            templateID = char(string(templateID));
            if ~strcmpi(templateID, 'OiS200')
                error('Umitoolbox:UMITRigStore:unknownBuiltInRig', ...
                    'Unknown built-in Rig definition: %s', templateID);
            end
            definition = UMITRigStore.iLoadDefaultRigDefinition();
        end

        function canonicalName = normalizeIlluminationName(name)
            %NORMALIZEILLUMINATIONNAME Return a canonical Rig illumination name.
            %
            %   Historical acquisition metadata may use "amber"; it maps
            %   to the canonical semantic name "yellow" without modifying
            %   AcqInfos.mat.

            if ~(ischar(name) || (isstring(name) && isscalar(name)))
                error('Umitoolbox:UMITRigStore:invalidIlluminationName', ...
                    'Illumination name must be a text scalar.');
            end
            canonicalName = lower(strtrim(char(string(name))));
            if strcmp(canonicalName, 'amber')
                canonicalName = 'yellow';
            end
            if ~ismember(canonicalName, getUMITRigSchema().canonicalIlluminations)
                error('Umitoolbox:UMITRigStore:invalidIlluminationName', ...
                    'Unsupported illumination "%s". Expected red, green, or yellow.', ...
                    canonicalName);
            end
        end

        function spectrum = getSpectrum(category, spectrumID)
            %GETSPECTRUM Resolve one canonical 400:700-nm optical spectrum.

            [category, spectrumID] = UMITRigStore.iNormalizeSpectrumIdentity( ...
                category, spectrumID);
            spectrumFile = UMITRigStore.iFindLibraryFile( ...
                category, spectrumID, '.txt');
            if isempty(spectrumFile)
                error('Umitoolbox:UMITRigStore:spectrumNotFound', ...
                    'Spectrum "%s" was not found in category "%s".', ...
                    spectrumID, category);
            end

            values = readmatrix(spectrumFile, 'FileType', 'text', 'CommentStyle', '#');
            values = values(:, 1:min(2, size(values, 2)));
            wavelength = getUMITRigSchema().spectrum.wavelengthNm;
            if size(values, 2) ~= 2 || size(values, 1) ~= numel(wavelength) || ...
                    any(~isfinite(values(:))) || any(values(:, 1) ~= wavelength) || ...
                    any(values(:, 2) < 0) || any(values(:, 2) > 1)
                error('Umitoolbox:UMITRigStore:invalidSpectrum', ...
                    ['Spectrum "%s" must contain exactly 301 finite rows on ' ...
                    '400:700 nm with response in [0,1].'], spectrumID);
            end

            metadata = UMITRigStore.iReadSpectrumMetadata(spectrumFile);
            spectrum = struct( ...
                'id', spectrumID, ...
                'category', category, ...
                'wavelengthNm', values(:, 1), ...
                'response', values(:, 2), ...
                'metadata', metadata, ...
                'file', spectrumFile);
        end

        function spectra = listSpectra(category)
            %LISTSPECTRA List shared spectra without exposing library paths.
            %
            %   spectra = UMITRigStore.listSpectra(category)
            %
            %   CATEGORY is "illumination", "camera", or "filter". The
            %   returned table contains ID, Category, Origin, DisplayName,
            %   Manufacturer, and Model. Origin is "user" or "builtIn".

            [category, ~] = UMITRigStore.iNormalizeSpectrumIdentity( ...
                category, 'LibraryEntry');
            schema = getUMITRigSchema();
            userFolder = UMITRigStore.iUserSpectrumFolder(category);
            builtInFolder = fullfile(fileparts(mfilename('fullpath')), ...
                'resources', schema.spectrum.libraryFolder, ...
                schema.spectrum.spectraFolder, category);
            [files, origins] = UMITRigStore.iListLibraryFiles( ...
                {userFolder, builtInFolder}, ["user", "builtIn"], '.txt');

            ID = strings(0, 1);
            Category = strings(0, 1);
            Origin = strings(0, 1);
            DisplayName = strings(0, 1);
            Manufacturer = strings(0, 1);
            Model = strings(0, 1);
            for iFile = 1:numel(files)
                [~, spectrumID] = fileparts(files{iFile});
                spectrum = UMITRigStore.getSpectrum(category, spectrumID);
                ID(end+1, 1) = string(spectrum.id); %#ok<AGROW>
                Category(end+1, 1) = string(category); %#ok<AGROW>
                Origin(end+1, 1) = origins(iFile); %#ok<AGROW>
                DisplayName(end+1, 1) = UMITRigStore.iMetadataText( ...
                    spectrum.metadata, 'displayName', spectrum.id); %#ok<AGROW>
                Manufacturer(end+1, 1) = UMITRigStore.iMetadataText( ...
                    spectrum.metadata, 'manufacturer', ''); %#ok<AGROW>
                Model(end+1, 1) = UMITRigStore.iMetadataText( ...
                    spectrum.metadata, 'model', ''); %#ok<AGROW>
            end
            spectra = table(ID, Category, Origin, DisplayName, Manufacturer, Model);
            if ~isempty(spectra)
                spectra = sortrows(spectra, 'ID');
            end
        end

        function spectrum = importSpectrum(sourceFile, category, spectrumID, metadata)
            %IMPORTSPECTRUM Normalize and import a shared optical spectrum.
            %
            %   The source must contain wavelength and response columns and
            %   cover the full 400:700-nm interval. Values are interpolated
            %   to 1-nm spacing and scaled by their maximum into [0,1].

            if nargin < 4
                metadata = struct();
            end
            [category, spectrumID] = UMITRigStore.iNormalizeSpectrumIdentity( ...
                category, spectrumID);
            if ~isstruct(metadata) || ~isscalar(metadata)
                error('Umitoolbox:UMITRigStore:invalidSpectrumMetadata', ...
                    'Spectrum metadata must be a scalar struct.');
            end
            sourceFile = UMITRigStore.iAbsolutePath(sourceFile);
            if ~isfile(sourceFile)
                error('Umitoolbox:UMITRigStore:invalidSpectrum', ...
                    'Spectrum source file does not exist: %s', sourceFile);
            end
            source = readmatrix(sourceFile);
            if size(source, 2) < 2
                error('Umitoolbox:UMITRigStore:invalidSpectrum', ...
                    'Spectrum source must contain wavelength and response columns.');
            end
            source = source(:, 1:2);
            source = source(all(isfinite(source), 2), :);
            [source(:, 1), order] = sort(source(:, 1));
            source(:, 2) = source(order, 2);
            if size(source, 1) < 2 || any(diff(source(:, 1)) <= 0) || ...
                    source(1, 1) > 400 || source(end, 1) < 700 || ...
                    any(source(:, 2) < 0) || max(source(:, 2)) <= 0
                error('Umitoolbox:UMITRigStore:invalidSpectrum', ...
                    ['Spectrum must have unique increasing wavelengths covering ' ...
                    '400:700 nm and nonnegative, nonzero responses.']);
            end
            wavelength = getUMITRigSchema().spectrum.wavelengthNm;
            response = interp1(source(:, 1), source(:, 2), wavelength, 'linear');
            response = response ./ max(response);
            if any(~isfinite(response)) || any(response < 0) || any(response > 1)
                error('Umitoolbox:UMITRigStore:invalidSpectrum', ...
                    'Spectrum could not be normalized to the canonical grid.');
            end

            storeLock = UMITRigStore.iAcquireStoreLock('importSpectrum'); %#ok<NASGU>
            destinationFolder = UMITRigStore.iUserSpectrumFolder(category);
            if ~isfolder(destinationFolder)
                mkdir(destinationFolder);
            end
            destination = fullfile(destinationFolder, [spectrumID '.txt']);
            if ~isempty(UMITRigStore.iFindLibraryFile(category, spectrumID, '.txt'))
                error('Umitoolbox:UMITRigStore:spectrumAlreadyExists', ...
                    ['Spectrum "%s" already exists. Create a new physical-profile ' ...
                     'ID instead of replacing a referenced spectrum.'], spectrumID);
            end
            UMITRigStore.iWriteSpectrumFile( ...
                destination, wavelength, response, metadata);
            spectrum = UMITRigStore.getSpectrum(category, spectrumID);
        end

        function removeSpectrum(category, spectrumID)
            %REMOVESPECTRUM Delete an unreferenced user-library spectrum.
            %
            %   Built-in spectra cannot be deleted. A user spectrum is also
            %   protected while referenced by any Rig camera/illumination or
            %   filter-set definition.

            [category, spectrumID] = UMITRigStore.iNormalizeSpectrumIdentity( ...
                category, spectrumID);
            storeLock = UMITRigStore.iAcquireStoreLock('removeSpectrum'); %#ok<NASGU>
            userFile = fullfile(UMITRigStore.iUserSpectrumFolder(category), ...
                [spectrumID '.txt']);
            if ~isfile(userFile)
                error('Umitoolbox:UMITRigStore:spectrumRemovalFailed', ...
                    'Only user-library spectra can be removed: %s.', spectrumID);
            end

            references = strings(0, 1);
            if strcmp(category, 'filter')
                filterFiles = UMITRigStore.iAllFilterSetFiles();
                for iFile = 1:numel(filterFiles)
                    definition = jsondecode(fileread(filterFiles{iFile}));
                    refs = string({definition.excitationSpectrumID, ...
                        definition.emissionSpectrumID});
                    if any(strcmpi(refs, spectrumID))
                        references(end+1, 1) = string(definition.id); %#ok<AGROW>
                    end
                end
            else
                rigs = UMITRigStore.listRigs();
                for iRig = find(rigs.IsReadable)'
                    rigStore = UMITRigStore.open(char(rigs.RigUUID(iRig)));
                    info = rigStore.getRigInfo();
                    if strcmp(category, 'camera')
                        records = info.cameras;
                    else
                        records = info.illuminations;
                    end
                    if ~isempty(records) && any(strcmpi({records.spectrumID}, spectrumID))
                        references(end+1, 1) = string(info.rigID); %#ok<AGROW>
                    end
                end
            end
            if ~isempty(references)
                error('Umitoolbox:UMITRigStore:spectrumInUse', ...
                    'Spectrum "%s" is referenced by: %s.', ...
                    spectrumID, char(strjoin(unique(references), ', ')));
            end
            delete(userFile);
            % Remove sidecars created by the short-lived paired-file format.
            metadataFile = replace(userFile, '.txt', '.json');
            if isfile(metadataFile)
                delete(metadataFile);
            end
        end

        function filterSet = getFilterSet(filterSetID)
            %GETFILTERSET Resolve and validate one filter-set definition.

            schema = getUMITRigSchema();
            filterSetID = UMITRigStore.iNormalizeManagedID(schema, filterSetID, 'filterSetID');
            filterFile = UMITRigStore.iFindFilterSetFile(filterSetID);
            if isempty(filterFile)
                error('Umitoolbox:UMITRigStore:filterSetNotFound', ...
                    'Filter set "%s" was not found.', filterSetID);
            end
            filterSet = jsondecode(fileread(filterFile));
            filterSet = UMITRigStore.iValidateFilterSet(filterSet, filterSetID);
            filterSet.file = filterFile;
        end

        function filterSets = listFilterSets()
            %LISTFILTERSETS List shared filter sets without exposing paths.
            %
            %   The returned table contains ID, Origin, DisplayName,
            %   ExcitationSpectrumID, and EmissionSpectrumID. Origin is
            %   "user" or "builtIn".

            [files, origins] = UMITRigStore.iAllFilterSetFiles();
            ID = strings(0, 1);
            Origin = strings(0, 1);
            DisplayName = strings(0, 1);
            ExcitationSpectrumID = strings(0, 1);
            EmissionSpectrumID = strings(0, 1);
            for iFile = 1:numel(files)
                [~, filterSetID] = fileparts(files{iFile});
                filterSet = UMITRigStore.getFilterSet(filterSetID);
                ID(end+1, 1) = string(filterSet.id); %#ok<AGROW>
                Origin(end+1, 1) = origins(iFile); %#ok<AGROW>
                DisplayName(end+1, 1) = string(filterSet.displayName); %#ok<AGROW>
                ExcitationSpectrumID(end+1, 1) = ...
                    string(filterSet.excitationSpectrumID); %#ok<AGROW>
                EmissionSpectrumID(end+1, 1) = ...
                    string(filterSet.emissionSpectrumID); %#ok<AGROW>
            end
            filterSets = table(ID, Origin, DisplayName, ...
                ExcitationSpectrumID, EmissionSpectrumID);
            if ~isempty(filterSets)
                filterSets = sortrows(filterSets, 'ID');
            end
        end

        function filterSet = importFilterSet(filterSet)
            %IMPORTFILTERSET Validate and save a shared filter-set definition.

            if ~isstruct(filterSet) || ~isscalar(filterSet) || ~isfield(filterSet, 'id')
                error('Umitoolbox:UMITRigStore:invalidFilterSet', ...
                    'Filter-set input must be a scalar struct containing id.');
            end
            filterSet = UMITRigStore.iValidateFilterSet(filterSet, filterSet.id);
            storeLock = UMITRigStore.iAcquireStoreLock('importFilterSet'); %#ok<NASGU>
            libraryRoot = fullfile(UMITRigStore.getRigsRoot(), ...
                getUMITRigSchema().spectrum.libraryFolder, ...
                getUMITRigSchema().spectrum.filterSetsFolder);
            if ~isfolder(libraryRoot)
                mkdir(libraryRoot);
            end
            destination = fullfile(libraryRoot, [filterSet.id '.json']);
            if ~isempty(UMITRigStore.iFindFilterSetFile(filterSet.id))
                error('Umitoolbox:UMITRigStore:filterSetAlreadyExists', ...
                    ['Filter set "%s" already exists. Create a new ID instead ' ...
                     'of replacing a referenced definition.'], filterSet.id);
            end
            serializable = filterSet;
            if isfield(serializable, 'file')
                serializable = rmfield(serializable, 'file');
            end
            UMITRigStore.iWriteJSON(destination, serializable);
        end

        function filterSet = updateFilterSet(filterSetID, filterSet)
            %UPDATEFILTERSET Replace an existing user filter-set definition.
            %
            %   filterSet = UMITRigStore.updateFilterSet(filterSetID, definition)
            %   validates the complete replacement before taking the store
            %   lock. Built-in definitions and filter-set IDs are immutable.

            schema = getUMITRigSchema();
            filterSetID = UMITRigStore.iNormalizeManagedID( ...
                schema, filterSetID, 'filterSetID');
            if ~isstruct(filterSet) || ~isscalar(filterSet)
                error('Umitoolbox:UMITRigStore:invalidFilterSet', ...
                    'Filter-set input must be a scalar struct.');
            end
            filterSet.id = filterSetID;
            filterSet = UMITRigStore.iValidateFilterSet(filterSet, filterSetID);

            userFolder = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.spectrum.libraryFolder, schema.spectrum.filterSetsFolder);
            storeLock = UMITRigStore.iAcquireStoreLock('updateFilterSet'); %#ok<NASGU>
            destination = UMITRigStore.iFindCaseInsensitiveFile( ...
                {userFolder}, [filterSetID '.json']);
            if isempty(destination)
                if isempty(UMITRigStore.iFindFilterSetFile(filterSetID))
                    error('Umitoolbox:UMITRigStore:filterSetNotFound', ...
                        'Filter set "%s" was not found.', filterSetID);
                end
                error('Umitoolbox:UMITRigStore:builtInFilterSetReadOnly', ...
                    'Built-in filter set "%s" cannot be edited.', filterSetID);
            end

            serializable = filterSet;
            if isfield(serializable, 'file')
                serializable = rmfield(serializable, 'file');
            end
            UMITRigStore.iWriteJSON(destination, serializable);
        end

        function obj = getDefaultRig()
            %GETDEFAULTRIG Resolve the active store-level default Rig.
            %
            %   Returns [] when no default pointer exists. A dangling,
            %   unreadable, or archived pointer is a configuration error.

            rigUUID = UMITRigStore.iReadDefaultRigUUID();
            if isempty(rigUUID)
                obj = [];
                return
            end

            try
                obj = UMITRigStore.open(rigUUID);
            catch ME
                error('Umitoolbox:UMITRigStore:invalidDefaultRig', ...
                    'The default Rig pointer does not resolve: %s', ME.message);
            end
            RigInfo = obj.getRigInfo();
            if ~strcmp(RigInfo.status, 'active')
                error('Umitoolbox:UMITRigStore:invalidDefaultRig', ...
                    'The default Rig "%s" is archived.', RigInfo.rigID);
            end
        end

        function obj = setDefaultRig(rigUUID)
            %SETDEFAULTRIG Select exactly one active store-level default Rig.

            storeLock = UMITRigStore.iAcquireStoreLock('setDefaultRig'); %#ok<NASGU>
            obj = UMITRigStore.open(rigUUID);
            rigLock = obj.iAcquireWriteLock('setDefaultRig'); %#ok<NASGU>
            RigInfo = obj.getRigInfo();
            if ~strcmp(RigInfo.status, 'active')
                error('Umitoolbox:UMITRigStore:setDefaultRigFailed', ...
                    'Archived Rig "%s" cannot become the default Rig.', RigInfo.rigID);
            end

            UMITRigStore.iWriteDefaultRigUUID(RigInfo.uuid);
        end

        function [obj, wasCreated, resolution] = getOrCreateDefaultRig()
            %GETORCREATEDEFAULTRIG Resolve "the" rig for this machine.
            %
            %   obj = UMITRigStore.getOrCreateDefaultRig()
            %
            %   Intended for entry points (e.g. raw-data import) that need a
            %   rig but have no way to ask the user which one. Resolution
            %   order:
            %       1) The active Rig named by the store-level pointer.
            %       2) If exactly one active Rig exists, it is promoted to default
            %          (the store pointer is persisted) and returned -- a single rig
            %          per machine is the common case.
            %       3) If no rig exists at all, one is created automatically,
            %          from the built-in OiS200 definition and selected.
            %       4) If multiple active Rigs exist and none is selected,
            %          this is genuinely ambiguous and throws -- mark one as
            %          default via setDefaultRig, or open a specific rig
            %          by rigID/rigUUID instead.
            %
            %   Archived rigs are excluded from every resolution branch --
            %   an archived rig can still be opened explicitly (openByRigID/
            %   open/rigExists), but is never auto-selected for new work.

            errID = 'Umitoolbox:UMITRigStore:noDefaultRig';
            wasCreated = false;
            resolution = 'existingDefault';
            storeLock = UMITRigStore.iAcquireStoreLock('getOrCreateDefaultRig'); %#ok<NASGU>

            obj = UMITRigStore.getDefaultRig();
            if ~isempty(obj)
                return
            end

            legacyDefaultUUID = UMITRigStore.iFindLegacyDefaultUUID();
            if ~isempty(legacyDefaultUUID)
                obj = UMITRigStore.open(legacyDefaultUUID);
                UMITRigStore.iWriteDefaultRigUUID(legacyDefaultUUID);
                resolution = 'promotedLegacyDefault';
                return
            end

            rigs = UMITRigStore.listRigs();
            rigs = rigs(rigs.Status ~= "archived", :);

            readableRigs = rigs(rigs.IsReadable, :);

            if height(readableRigs) == 1
                obj = UMITRigStore.open(char(readableRigs.RigUUID(1)));
                UMITRigStore.iWriteDefaultRigUUID(obj.getRigInfo().uuid);
                resolution = 'promotedExisting';
                return
            end

            if height(readableRigs) == 0
                definition = UMITRigStore.iLoadDefaultRigDefinition();
                if isfolder(fullfile(UMITRigStore.getRigsRoot(), definition.rigID))
                    definition.rigID = [definition.rigID '_' ...
                        strrep(UMITRigStore.iGenerateUUID(), '-', '')];
                    definition.rigID = definition.rigID( ...
                        1:min(64, numel(definition.rigID)));
                end
                obj = UMITRigStore.create(definition, 'InternalStoreLockHeld', true);
                UMITRigStore.iWriteDefaultRigUUID(obj.getRigInfo().uuid);
                wasCreated = true;
                resolution = 'createdDefault';
                return
            end

            error(errID, ...
                ['Multiple active rigs exist and no store-level default is selected. ' ...
                'Call UMITRigStore.setDefaultRig, or open a specific rig by ' ...
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

        function newStore = duplicate(obj, newRigID, varargin)
            %DUPLICATE Create an independent active copy of this Rig.
            %
            %   newStore = sourceStore.duplicate(newRigID)
            %   newStore = sourceStore.duplicate(newRigID, ...
            %       'DisplayName', displayName)
            %
            %   Description, metadata, cameras, and illuminations are
            %   copied. The new Rig receives fresh identity and timestamps.
            %   Archived/default state, active coregistration, the resource
            %   registry, and managed resource files are not copied. By
            %   default the display name is the source display name followed
            %   by " Copy"; callers may provide an explicit DisplayName.

            errID = 'Umitoolbox:UMITRigStore:duplicateFailed';
            p = inputParser;
            p.FunctionName = 'UMITRigStore.duplicate';
            addParameter(p, 'DisplayName', '', ...
                @(x) ischar(x) || (isstring(x) && isscalar(x)));
            parse(p, varargin{:});

            obj.iAssertWritable();
            storeLock = UMITRigStore.iAcquireStoreLock('duplicateRig'); %#ok<NASGU>
            sourceLock = obj.iAcquireWriteLock('duplicateRig'); %#ok<NASGU>
            report = obj.validate('Mode', 'full');
            if ~report.isValid
                error(errID, 'Source Rig validation failed: %s', ...
                    UMITRigStore.iJoinIssueMessages(report.errors));
            end
            sourceInfo = obj.iLoadRigInfo();
            displayName = char(string(p.Results.DisplayName));
            if isempty(displayName)
                if isempty(sourceInfo.displayName)
                    displayName = char(string(newRigID));
                else
                    displayName = [sourceInfo.displayName ' Copy'];
                end
            end
            duplicateInfo = struct( ...
                'rigID', newRigID, ...
                'displayName', displayName, ...
                'description', sourceInfo.description, ...
                'metadata', sourceInfo.metadata, ...
                'cameras', sourceInfo.cameras, ...
                'illuminations', sourceInfo.illuminations);
            newStore = UMITRigStore.create(duplicateInfo, ...
                'InternalStoreLockHeld', true);
        end

        function setCameras(obj, cameras)
            %SETCAMERAS Replace the validated one- or two-camera definition.

            cameras = UMITRigStore.iValidateCameraRecords(obj.Schema, cameras, true);
            obj.iUpdateHardwareField('cameras', cameras);
        end

        function setIlluminations(obj, illuminations)
            %SETILLUMINATIONS Replace canonical Rig illumination definitions.

            illuminations = UMITRigStore.iValidateIlluminationRecords( ...
                obj.Schema, illuminations, true);
            obj.iUpdateHardwareField('illuminations', illuminations);
        end

        function setHardwareConfiguration(obj, cameras, illuminations)
            %SETHARDWARECONFIGURATION Atomically replace all Rig hardware.
            %
            %   store.setHardwareConfiguration(cameras, illuminations)
            %
            %   Both collections are validated before either is persisted.
            %   The update is written under one Rig/store lock and rolls
            %   back to the original metadata if post-write validation
            %   fails. Use this method for staged GUI configuration edits.

            cameras = UMITRigStore.iValidateCameraRecords( ...
                obj.Schema, cameras, true);
            illuminations = UMITRigStore.iValidateIlluminationRecords( ...
                obj.Schema, illuminations, true);

            obj.iAssertWritable();
            storeLock = UMITRigStore.iAcquireStoreLock( ...
                'setHardwareConfiguration'); %#ok<NASGU>
            lockCleanup = obj.iAcquireWriteLock( ...
                'setHardwareConfiguration'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            RigInfo.cameras = cameras;
            RigInfo.illuminations = illuminations;
            RigInfo.modifiedOn = datetime('now');

            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end

            obj.iAppendLog( ...
                'setHardwareConfiguration', RigInfo.uuid, 'completed');
        end

        function optical = resolveOpticalConfiguration(obj, acqInfo, requestedIlluminations, filterSetID, dataFolder)
            %RESOLVEOPTICALCONFIGURATION Preflight acquisition and Rig optics.
            %
            %   Resolves acquisition-specific channel files and CamIdx values
            %   against this Rig's physical illumination/camera spectra and
            %   one shared filter-set definition. Archived Rigs remain valid
            %   for analysis of historical data.

            if nargin < 5
                dataFolder = '';
            end
            if ~isstruct(acqInfo) || ~isscalar(acqInfo) || ...
                    ~isfield(acqInfo, 'ImportedChannels') || ...
                    ~isstruct(acqInfo.ImportedChannels)
                error('Umitoolbox:UMITRigStore:invalidAcqInfos', ...
                    'AcqInfoStream must contain ImportedChannels metadata.');
            end
            if ischar(requestedIlluminations) || isstring(requestedIlluminations)
                requestedIlluminations = cellstr(requestedIlluminations);
            end
            if ~iscell(requestedIlluminations) || isempty(requestedIlluminations)
                error('Umitoolbox:UMITRigStore:invalidIlluminationName', ...
                    'Requested illuminations must be a nonempty text list.');
            end
            requested = cellfun(@UMITRigStore.normalizeIlluminationName, ...
                requestedIlluminations, 'UniformOutput', false);
            if numel(unique(requested)) < 2
                error('Umitoolbox:UMITRigStore:insufficientIlluminations', ...
                    'At least two different illuminations are required.');
            end

            RigInfo = obj.iLoadRigInfo();
            cameras = UMITRigStore.iValidateCameraRecords(obj.Schema, RigInfo.cameras, true);
            rigIlluminations = UMITRigStore.iValidateIlluminationRecords( ...
                obj.Schema, RigInfo.illuminations, true);
            filterSet = UMITRigStore.getFilterSet(filterSetID);
            excitation = UMITRigStore.iResolveFilterSpectrum(filterSet.excitationSpectrumID);
            emission = UMITRigStore.iResolveFilterSpectrum(filterSet.emissionSpectrumID);

            rowNames = obj.Schema.canonicalIlluminations;
            nRows = numel(rowNames);
            nSamples = numel(obj.Schema.spectrum.wavelengthNm);
            optical = struct();
            optical.wavelengthNm = obj.Schema.spectrum.wavelengthNm;
            optical.rowNames = rowNames;
            optical.activeRows = false(nRows, 1);
            optical.illuminationResponse = nan(nRows, nSamples);
            optical.cameraResponse = nan(nRows, nSamples);
            optical.excitationResponse = excitation(:)';
            optical.emissionResponse = emission(:)';
            optical.filterSet = filterSet;
            optical.rigUUID = RigInfo.uuid;
            optical.rigID = RigInfo.rigID;
            optical.channels = struct('name', {}, 'datFile', {}, 'camIdx', {});

            imported = acqInfo.ImportedChannels;
            importedNames = cell(1, numel(imported));
            for iChannel = 1:numel(imported)
                sourceName = '';
                if isfield(imported, 'Color') && ~isempty(imported(iChannel).Color)
                    sourceName = imported(iChannel).Color;
                elseif isfield(imported, 'Tag') && ~isempty(imported(iChannel).Tag)
                    sourceName = imported(iChannel).Tag;
                elseif isfield(imported, 'DatFile')
                    [~, sourceName] = fileparts(imported(iChannel).DatFile);
                end
                try
                    importedNames{iChannel} = UMITRigStore.normalizeIlluminationName(sourceName);
                catch
                    importedNames{iChannel} = '';
                end
            end

            for iRequested = 1:numel(requested)
                name = requested{iRequested};
                channelIdx = find(strcmp(importedNames, name));
                if isempty(channelIdx)
                    error('Umitoolbox:UMITRigStore:illuminationNotAcquired', ...
                        'Requested illumination "%s" was not acquired.', name);
                elseif numel(channelIdx) > 1
                    error('Umitoolbox:UMITRigStore:ambiguousAcquisitionChannel', ...
                        'Requested illumination "%s" maps to multiple imported channels.', name);
                end
                channel = imported(channelIdx);
                if ~isfield(channel, 'DatFile') || isempty(channel.DatFile)
                    error('Umitoolbox:UMITRigStore:invalidAcqInfos', ...
                        'Acquisition channel "%s" has no DatFile.', name);
                end
                if ~isfield(channel, 'CamIdx') || isempty(channel.CamIdx) || ...
                        ~isnumeric(channel.CamIdx) || ~isscalar(channel.CamIdx) || ...
                        ~ismember(channel.CamIdx, [1 2])
                    error('Umitoolbox:UMITRigStore:invalidCameraIndex', ...
                        'Acquisition channel "%s" has no valid CamIdx.', name);
                end
                camIdx = double(channel.CamIdx);
                cameraRecord = cameras([cameras.index] == camIdx);
                if isempty(cameraRecord)
                    error('Umitoolbox:UMITRigStore:cameraNotConfigured', ...
                        'Rig "%s" does not define camera %d.', RigInfo.rigID, camIdx);
                end
                if isempty(cameraRecord.spectrumID)
                    error('Umitoolbox:UMITRigStore:missingCameraSpectrum', ...
                        'Camera %d has no spectral profile configured.', camIdx);
                end
                cameraSpectrum = UMITRigStore.getSpectrum('camera', cameraRecord.spectrumID);

                illuminationRecord = rigIlluminations(strcmp({rigIlluminations.name}, name));
                if isempty(illuminationRecord) || isempty(illuminationRecord.spectrumID)
                    error('Umitoolbox:UMITRigStore:missingIlluminationSpectrum', ...
                        'Rig "%s" does not define a spectrum for illumination "%s".', ...
                        RigInfo.rigID, name);
                end
                illuminationSpectrum = UMITRigStore.getSpectrum( ...
                    'illumination', illuminationRecord.spectrumID);
                row = find(strcmp(rowNames, name), 1);
                optical.activeRows(row) = true;
                optical.illuminationResponse(row, :) = illuminationSpectrum.response(:)';
                optical.cameraResponse(row, :) = cameraSpectrum.response(:)';
                optical.channels(end+1) = struct( ... %#ok<AGROW>
                    'name', name, 'datFile', char(string(channel.DatFile)), 'camIdx', camIdx);
                if ~isempty(dataFolder) && ...
                        ~isfile(fullfile(char(string(dataFolder)), char(string(channel.DatFile))))
                    error('Umitoolbox:UMITRigStore:missingChannelFile', ...
                        'Resolved channel file was not found: %s', ...
                        fullfile(char(string(dataFolder)), char(string(channel.DatFile))));
                end
            end
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
            %   different rig as default first (setDefaultRig).

            errID = 'Umitoolbox:UMITRigStore:archiveRigFailed';
            obj.iAssertWritable();
            storeLock = UMITRigStore.iAcquireStoreLock('archiveRig'); %#ok<NASGU>
            lockCleanup = obj.iAcquireWriteLock('archiveRig'); %#ok<NASGU>
            obj.iAssertHealthyForMutation();

            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;

            if strcmp(RigInfo.status, 'archived')
                error(errID, 'Rig is already archived: %s', RigInfo.rigID);
            end

            defaultUUID = UMITRigStore.iReadDefaultRigUUID();
            if strcmpi(defaultUUID, RigInfo.uuid)
                error(errID, ...
                    ['Cannot archive the default Rig "%s". Select a different ' ...
                    'default with UMITRigStore.setDefaultRig first.'], RigInfo.rigID);
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

        function restoreRig(obj)
            %RESTORERIG Restore an archived Rig to active status.

            errID = 'Umitoolbox:UMITRigStore:restoreRigFailed';
            obj.iAssertWritable();
            lockCleanup = obj.iAcquireWriteLock('restoreRig'); %#ok<NASGU>
            report = obj.validate('Mode', 'full');
            if ~report.isValid
                error(errID, 'Rig validation failed: %s', ...
                    UMITRigStore.iJoinIssueMessages(report.errors));
            end

            RigInfo = obj.iLoadRigInfo();
            if ~strcmp(RigInfo.status, 'archived')
                error(errID, 'Rig is not archived: %s', RigInfo.rigID);
            end
            originalRigInfo = RigInfo;
            RigInfo.status = 'active';
            RigInfo.archivedOn = NaT;
            RigInfo.modifiedOn = datetime('now');
            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end
            obj.iAppendLog('restoreRig', RigInfo.uuid, 'completed');
        end

        function resourceUUID = addCameraCoregistration(obj, sourceFile, resourceInfo)
            %ADDCAMERACOREGISTRATION Import a camera-coregistration transform.
            %
            %   The first imported transform becomes active. Later imports
            %   remain available until explicitly selected.

            if nargin < 3
                resourceInfo = struct();
            end
            resourceUUID = obj.iAddManagedResource( ...
                'cameraCoregistration', sourceFile, resourceInfo, false);
        end

        function resourceUUID = importAndActivateCameraCoregistration(obj, sourceFile, resourceInfo)
            %IMPORTANDACTIVATECAMERACOREGISTRATION Import and activate atomically.

            if nargin < 3
                resourceInfo = struct();
            end
            resourceUUID = obj.iAddManagedResource( ...
                'cameraCoregistration', sourceFile, resourceInfo, true);
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

        function resource = getActiveCameraCoregistration(obj)
            %GETACTIVECAMERACOREGISTRATION Return the active camera transform.

            resource = obj.iGetActiveResource('cameraCoregistration');
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

            try
                UMITRigStore.iValidateRigIdentity(obj.Schema, RigInfo, obj.RigRoot);
            catch ME
                report = obj.iAddIssue(report, 'error', 'invalidRigIdentity', ME.message);
            end

            if isfield(RigInfo, 'status') && isfield(obj.Schema, 'rigStatuses') && ...
                    ~ismember(RigInfo.status, obj.Schema.rigStatuses)
                report = obj.iAddIssue(report, 'error', 'invalidRigStatus', ...
                    sprintf('Rig status "%s" is not one of the allowed values.', RigInfo.status));
            end

            if isfield(RigInfo, 'cameras')
                try
                    UMITRigStore.iValidateCameraRecords( ...
                        obj.Schema, RigInfo.cameras, strcmp(mode, 'full'));
                catch ME
                    report = obj.iAddIssue(report, 'error', 'invalidCameras', ME.message);
                end
            end
            if isfield(RigInfo, 'illuminations')
                try
                    UMITRigStore.iValidateIlluminationRecords( ...
                        obj.Schema, RigInfo.illuminations, strcmp(mode, 'full'));
                catch ME
                    report = obj.iAddIssue(report, 'error', 'invalidIlluminations', ME.message);
                end
            end

            if strcmp(mode, 'full') && isfield(RigInfo, 'resourceRegistry')
                if ~isstruct(RigInfo.resourceRegistry)
                    report = obj.iAddIssue(report, 'error', 'invalidResourceRegistry', ...
                        'resourceRegistry must be a struct array.');
                else
                    resourceUUIDs = strings(0, 1);
                    for iRecord = 1:numel(RigInfo.resourceRegistry)
                        record = RigInfo.resourceRegistry(iRecord);
                        try
                            UMITRigStore.iValidateResourceRecord(obj.Schema, record);
                            resourceUUIDs(end+1, 1) = string(record.uuid); %#ok<AGROW>
                        catch ME
                            report = obj.iAddIssue(report, 'error', 'invalidResourceRecord', ME.message);
                            continue
                        end

                        expectedState = obj.Schema.folders.active;
                        if strcmp(record.status, 'archived')
                            expectedState = obj.Schema.folders.archive;
                        end
                        expectedPath = obj.iBuildResourceRelativePath( ...
                            record.type, expectedState, record.fileName);
                        try
                            if ~strcmpi(strrep(record.relativePath, '\\', '/'), expectedPath)
                                error('Umitoolbox:UMITRigStore:invalidResourcePath', ...
                                    'Resource %s has a non-canonical managed path.', record.uuid);
                            end
                            resourcePath = obj.iResolveRelativePath(record.relativePath);
                        catch ME
                            report = obj.iAddIssue(report, 'error', 'invalidResourcePath', ME.message);
                            continue
                        end
                        if ~isfile(resourcePath)
                            report = obj.iAddIssue(report, 'error', 'missingResourceFile', ...
                                sprintf('Resource %s file is missing: %s', record.uuid, resourcePath));
                            continue
                        end
                        try
                            if ~strcmpi(computeFileChecksum(resourcePath), record.checksum)
                                error('Umitoolbox:UMITRigStore:resourceChecksumMismatch', ...
                                    'Resource %s checksum does not match its registry record.', record.uuid);
                            end
                            if strcmp(record.type, 'cameraCoregistration')
                                UMITRigStore.iValidateCameraCoregistrationFile(resourcePath);
                            end
                        catch ME
                            report = obj.iAddIssue(report, 'error', 'invalidResourcePayload', ME.message);
                        end
                    end
                    if numel(unique(lower(resourceUUIDs))) ~= numel(resourceUUIDs)
                        report = obj.iAddIssue(report, 'error', 'duplicateResourceUUID', ...
                            'resourceRegistry contains duplicate resource UUIDs.');
                    end

                    for iType = fieldnames(obj.Schema.resourceTypes)'
                        resourceType = iType{1};
                        resourceDef = obj.Schema.resourceTypes.(resourceType);
                        pointerField = resourceDef.activePointerField;
                        if ~isfield(RigInfo, pointerField)
                            report = obj.iAddIssue(report, 'error', 'missingActivePointer', ...
                                sprintf('Rig metadata is missing active pointer %s.', pointerField));
                            continue
                        end
                        activeUUID = RigInfo.(pointerField);
                        if ~(ischar(activeUUID) || (isstring(activeUUID) && isscalar(activeUUID)))
                            report = obj.iAddIssue(report, 'error', 'invalidActivePointer', ...
                                sprintf('%s must be a text scalar.', pointerField));
                            continue
                        end
                        if isempty(activeUUID)
                            continue
                        end
                        idx = obj.iFindResourceIndex(RigInfo.resourceRegistry, activeUUID);
                        if isempty(idx)
                            report = obj.iAddIssue(report, 'error', 'danglingActivePointer', ...
                                sprintf('%s points to a resource not in the registry: %s', ...
                                    pointerField, activeUUID));
                        elseif ~strcmp(RigInfo.resourceRegistry(idx).type, resourceType) || ...
                                ~strcmp(RigInfo.resourceRegistry(idx).status, 'active')
                            report = obj.iAddIssue(report, 'error', 'invalidActivePointer', ...
                                sprintf('%s does not point to an active %s resource: %s', ...
                                pointerField, resourceType, activeUUID));
                        end
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
            UMITRigStore.iRemoveObsoleteCalibrationFolder(obj);
        end

        function iUpdateHardwareField(obj, fieldName, value)
            %IUPDATEHARDWAREFIELD Atomically replace one validated hardware field.

            obj.iAssertWritable();
            storeLock = UMITRigStore.iAcquireStoreLock(['set_', fieldName]); %#ok<NASGU>
            lockCleanup = obj.iAcquireWriteLock(['set_', fieldName]); %#ok<NASGU>
            obj.iAssertHealthyForMutation();
            RigInfo = obj.iLoadRigInfo();
            originalRigInfo = RigInfo;
            RigInfo.(fieldName) = value;
            RigInfo.modifiedOn = datetime('now');
            try
                obj.iSaveRigInfo(RigInfo);
                obj.iAssertValidAfterMutation();
            catch ME
                obj.iSaveRigInfo(originalRigInfo);
                rethrow(ME)
            end
            obj.iAppendLog(['set_', fieldName], RigInfo.uuid, 'completed');
        end

        function resourceUUID = iAddManagedResource(obj, resourceType, sourceFile, resourceInfo, makeActive)
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
                sourceProbe = load(sourceFile, '-mat');
            catch ME
                error(errID, 'Source MAT file cannot be loaded: %s', ME.message);
            end
            if strcmp(resourceType, 'cameraCoregistration')
                if ~isfield(sourceProbe, 'tform') || isempty(sourceProbe.tform)
                    error(errID, ...
                        'Camera-coregistration MAT files must contain variable "tform".');
                end
                if isfield(sourceProbe, 'tformInfo') && ...
                        (~isstruct(sourceProbe.tformInfo) || ~isscalar(sourceProbe.tformInfo))
                    error(errID, 'Variable "tformInfo", when present, must be a scalar struct.');
                end
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

            if makeActive || isempty(RigInfo.(pointerField))
                previousActiveUUID = RigInfo.(pointerField);
                if ~isempty(previousActiveUUID)
                    previousIndex = obj.iFindResourceIndex( ...
                        RigInfo.resourceRegistry, previousActiveUUID);
                    if ~isempty(previousIndex)
                        RigInfo.resourceRegistry(previousIndex).status = 'available';
                        RigInfo.resourceRegistry(previousIndex).modifiedOn = datetime('now');
                    end
                end
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

            if ~(ischar(relativePath) || (isstring(relativePath) && isscalar(relativePath)))
                error('Umitoolbox:UMITRigStore:invalidRelativePath', ...
                    'Managed resource paths must be text scalars.');
            end
            relativePath = char(string(relativePath));
            if isempty(relativePath)
                error('Umitoolbox:UMITRigStore:invalidRelativePath', ...
                    'Managed resource paths cannot be empty.');
            end
            relativePath = strrep(relativePath, '/', filesep);
            rootPath = UMITRigStore.iAbsolutePath(obj.RigRoot);
            path = UMITRigStore.iAbsolutePath(fullfile(rootPath, relativePath));
            rootPrefix = [rootPath filesep];
            if ~strcmpi(path, rootPath) && ~startsWith(path, rootPrefix, 'IgnoreCase', true)
                error('Umitoolbox:UMITRigStore:invalidRelativePath', ...
                    'Managed resource path escapes the Rig folder: %s', relativePath);
            end
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
            RigInfo = obj.iLoadRigInfo();
            if strcmp(RigInfo.status, 'archived')
                error('Umitoolbox:UMITRigStore:archivedRigReadOnly', ...
                    'Archived Rig "%s" must be restored before modification.', RigInfo.rigID);
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
        function iValidateRigIdentity(schema, RigInfo, rigRoot)
            %IVALIDATERIGIDENTITY Validate immutable identity and folder binding.

            if ~isstruct(RigInfo) || ~isscalar(RigInfo)
                error('Umitoolbox:UMITRigStore:invalidRigMetadata', ...
                    'Rig metadata must be a scalar struct.');
            end
            UMITRigStore.iAssertUUIDMatchesSchema(schema, RigInfo.uuid, 'Rig UUID');
            rigID = UMITRigStore.iNormalizeManagedID(schema, RigInfo.rigID, 'rigID');
            [~, folderName] = fileparts(UMITRigStore.iAbsolutePath(rigRoot));
            if ~strcmpi(rigID, folderName)
                error('Umitoolbox:UMITRigStore:invalidRigIdentity', ...
                    'Rig metadata rigID does not match its folder name.');
            end
            if ~isnumeric(RigInfo.schemaVersion) || ~isscalar(RigInfo.schemaVersion) || ...
                    RigInfo.schemaVersion ~= schema.version
                error('Umitoolbox:UMITRigStore:invalidRigMetadata', ...
                    'Rig schemaVersion must be %d.', schema.version);
            end
        end

        function iValidateResourceRecord(schema, record)
            %IVALIDATERESOURCERECORD Validate one persisted managed resource.

            if ~isstruct(record) || ~isscalar(record)
                error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                    'Each resourceRegistry entry must be a scalar struct.');
            end
            missing = setdiff(schema.resourceRecordFields, fieldnames(record));
            if ~isempty(missing)
                error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                    'Resource record is missing field(s): %s.', strjoin(missing, ', '));
            end
            UMITRigStore.iAssertUUIDMatchesSchema(schema, record.uuid, 'Resource UUID');
            if ~(ischar(record.type) || (isstring(record.type) && isscalar(record.type))) || ...
                    ~isfield(schema.resourceTypes, char(string(record.type)))
                error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                    'Resource type is invalid.');
            end
            if ~(ischar(record.status) || (isstring(record.status) && isscalar(record.status))) || ...
                    ~ismember(char(string(record.status)), schema.resourceStates)
                error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                    'Resource status is invalid.');
            end
            textFields = {'fileName', 'relativePath', 'displayName', 'description', 'checksum', 'sourceFile'};
            for iField = 1:numel(textFields)
                value = record.(textFields{iField});
                if ~(ischar(value) || (isstring(value) && isscalar(value)))
                    error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                        'Resource field "%s" must be text.', textFields{iField});
                end
            end
            if isempty(record.fileName) || contains(record.fileName, {'/', '\\'}) || ...
                    ~strcmp(record.fileName, char(java.io.File(record.fileName).getName()))
                error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                    'Resource fileName must be a plain file name.');
            end
            if isempty(record.checksum)
                error('Umitoolbox:UMITRigStore:invalidResourceRecord', ...
                    'Resource checksum cannot be empty.');
            end
        end

        function iValidateCameraCoregistrationFile(filePath)
            try
                payload = load(filePath, '-mat');
            catch ME
                error('Umitoolbox:UMITRigStore:invalidCoregistrationPayload', ...
                    'Camera-coregistration MAT file cannot be loaded: %s', ME.message);
            end
            if ~isfield(payload, 'tform') || isempty(payload.tform)
                error('Umitoolbox:UMITRigStore:invalidCoregistrationPayload', ...
                    'Camera-coregistration MAT file must contain tform.');
            end
            if isfield(payload, 'tformInfo') && ...
                    (~isstruct(payload.tformInfo) || ~isscalar(payload.tformInfo))
                error('Umitoolbox:UMITRigStore:invalidCoregistrationPayload', ...
                    'Camera-coregistration tformInfo must be a scalar struct.');
            end
        end

        function iAssertUUIDMatchesSchema(schema, value, label)
            if ~(ischar(value) || (isstring(value) && isscalar(value))) || ...
                    isempty(regexp(char(string(value)), schema.namingRules.uuidPattern, 'once'))
                error('Umitoolbox:UMITRigStore:invalidUUID', ...
                    '%s must be a canonical UUID.', label);
            end
        end

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
            if ~isfield(RigInfo, 'cameras') || ~isstruct(RigInfo.cameras)
                RigInfo.cameras = UMITRigStore.iEmptyCameraRecords();
            end
            if ~isfield(RigInfo, 'illuminations') || ~isstruct(RigInfo.illuminations)
                RigInfo.illuminations = UMITRigStore.iEmptyIlluminationRecords();
            end
            if isfield(RigInfo, 'resourceRegistry') && ~isempty(RigInfo.resourceRegistry) && ...
                    isfield(RigInfo.resourceRegistry, 'type')
                RigInfo.resourceRegistry(strcmp( ...
                    {RigInfo.resourceRegistry.type}, 'calibrationFile')) = [];
            end
            obsoleteFields = intersect(fieldnames(RigInfo), ...
                {'isDefault', 'activeCalibrationFileUUID'});
            if ~isempty(obsoleteFields)
                RigInfo = rmfield(RigInfo, obsoleteFields);
            end
            if ~isfield(RigInfo, 'schemaVersion') || RigInfo.schemaVersion < 3
                RigInfo.schemaVersion = 3;
            end
        end

        function iRemoveObsoleteCalibrationFolder(obj)
            %IREMOVEOBSOLETECALIBRATIONFOLDER Remove retired v1/v2 resources.
            % The folder has no supported schema-v3 meaning. It is removed only
            % after the v3 metadata write succeeds, and its exact path is built
            % from the validated Rig root rather than from persisted metadata.

            legacyFolder = fullfile(obj.RigRoot, 'calibration-files');
            if isfolder(legacyFolder)
                UMITRigStore.iRemoveFolderIfPresent(legacyFolder);
            end
        end

        function records = iEmptyCameraRecords()
            template = struct( ...
                'index', {}, 'displayName', {}, 'manufacturer', {}, ...
                'model', {}, 'serialNumber', {}, 'spectrumID', {});
            records = template;
        end

        function records = iEmptyIlluminationRecords()
            template = struct( ...
                'name', {}, 'displayName', {}, 'manufacturer', {}, ...
                'model', {}, 'spectrumID', {});
            records = template;
        end

        function value = iGetStructArrayField(info, fieldName)
            if ~isfield(info, fieldName) || isempty(info.(fieldName))
                if strcmp(fieldName, 'cameras')
                    value = UMITRigStore.iEmptyCameraRecords();
                else
                    value = UMITRigStore.iEmptyIlluminationRecords();
                end
                return
            end
            value = info.(fieldName);
            if ~isstruct(value)
                error('Umitoolbox:UMITRigStore:createFailed', ...
                    'Field "%s" must be a struct array.', fieldName);
            end
        end

        function cameras = iValidateCameraRecords(schema, cameras, requireReferences)
            if isempty(cameras)
                cameras = UMITRigStore.iEmptyCameraRecords();
                return
            end
            if ~isstruct(cameras) || numel(cameras) > 2
                error('Umitoolbox:UMITRigStore:invalidCameras', ...
                    'Camera definitions must be a struct array with at most two records.');
            end
            missing = setdiff(schema.cameraFields, fieldnames(cameras));
            if ~isempty(missing)
                error('Umitoolbox:UMITRigStore:invalidCameras', ...
                    'Camera definitions are missing field(s): %s.', strjoin(missing, ', '));
            end
            indices = [cameras.index];
            if ~isnumeric(indices) || any(~ismember(indices, [1 2])) || ...
                    numel(unique(indices)) ~= numel(indices)
                error('Umitoolbox:UMITRigStore:invalidCameras', ...
                    'Camera indices must be unique values 1 or 2.');
            end
            textFields = {'displayName', 'manufacturer', 'model', 'serialNumber', 'spectrumID'};
            for iCamera = 1:numel(cameras)
                for iField = 1:numel(textFields)
                    fieldName = textFields{iField};
                    value = cameras(iCamera).(fieldName);
                    if ~(ischar(value) || (isstring(value) && isscalar(value)))
                        error('Umitoolbox:UMITRigStore:invalidCameras', ...
                            'Camera field "%s" must be text.', fieldName);
                    end
                    cameras(iCamera).(fieldName) = char(string(value));
                end
                if requireReferences && ~isempty(cameras(iCamera).spectrumID)
                    UMITRigStore.getSpectrum('camera', cameras(iCamera).spectrumID);
                end
            end
            [~, order] = sort([cameras.index]);
            cameras = cameras(order);
        end

        function illuminations = iValidateIlluminationRecords(schema, illuminations, requireReferences)
            if isempty(illuminations)
                illuminations = UMITRigStore.iEmptyIlluminationRecords();
                return
            end
            if ~isstruct(illuminations)
                error('Umitoolbox:UMITRigStore:invalidIlluminations', ...
                    'Illumination definitions must be a struct array.');
            end
            missing = setdiff(schema.illuminationFields, fieldnames(illuminations));
            if ~isempty(missing)
                error('Umitoolbox:UMITRigStore:invalidIlluminations', ...
                    'Illumination definitions are missing field(s): %s.', strjoin(missing, ', '));
            end
            names = cell(1, numel(illuminations));
            textFields = {'displayName', 'manufacturer', 'model', 'spectrumID'};
            for iIllumination = 1:numel(illuminations)
                names{iIllumination} = UMITRigStore.normalizeIlluminationName( ...
                    illuminations(iIllumination).name);
                illuminations(iIllumination).name = names{iIllumination};
                for iField = 1:numel(textFields)
                    fieldName = textFields{iField};
                    value = illuminations(iIllumination).(fieldName);
                    if ~(ischar(value) || (isstring(value) && isscalar(value)))
                        error('Umitoolbox:UMITRigStore:invalidIlluminations', ...
                            'Illumination field "%s" must be text.', fieldName);
                    end
                    illuminations(iIllumination).(fieldName) = char(string(value));
                end
                if requireReferences && ~isempty(illuminations(iIllumination).spectrumID)
                    UMITRigStore.getSpectrum( ...
                        'illumination', illuminations(iIllumination).spectrumID);
                end
            end
            if numel(unique(names)) ~= numel(names)
                error('Umitoolbox:UMITRigStore:invalidIlluminations', ...
                    'Canonical illumination names must be unique.');
            end
        end

        function [category, spectrumID] = iNormalizeSpectrumIdentity(category, spectrumID)
            schema = getUMITRigSchema();
            if ~(ischar(category) || (isstring(category) && isscalar(category)))
                error('Umitoolbox:UMITRigStore:invalidSpectrumCategory', ...
                    'Spectrum category must be text.');
            end
            category = lower(strtrim(char(string(category))));
            if ~ismember(category, schema.spectrum.categories)
                error('Umitoolbox:UMITRigStore:invalidSpectrumCategory', ...
                    'Spectrum category must be illumination, filter, or camera.');
            end
            spectrumID = UMITRigStore.iNormalizeManagedID( ...
                schema, spectrumID, 'spectrumID');
        end

        function folder = iUserSpectrumFolder(category)
            schema = getUMITRigSchema();
            folder = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.spectrum.libraryFolder, schema.spectrum.spectraFolder, category);
        end

        function filePath = iFindLibraryFile(category, itemID, extension)
            schema = getUMITRigSchema();
            userFolder = UMITRigStore.iUserSpectrumFolder(category);
            builtInFolder = fullfile(fileparts(mfilename('fullpath')), ...
                'resources', schema.spectrum.libraryFolder, ...
                schema.spectrum.spectraFolder, category);
            filePath = UMITRigStore.iFindCaseInsensitiveFile( ...
                {userFolder, builtInFolder}, [itemID extension]);
        end

        function filePath = iFindFilterSetFile(filterSetID)
            schema = getUMITRigSchema();
            userFolder = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.spectrum.libraryFolder, schema.spectrum.filterSetsFolder);
            builtInFolder = fullfile(fileparts(mfilename('fullpath')), ...
                'resources', schema.spectrum.libraryFolder, ...
                schema.spectrum.filterSetsFolder);
            filePath = UMITRigStore.iFindCaseInsensitiveFile( ...
                {userFolder, builtInFolder}, [filterSetID '.json']);
        end

        function [files, origins] = iAllFilterSetFiles()
            schema = getUMITRigSchema();
            folders = { ...
                fullfile(UMITRigStore.getRigsRoot(), ...
                    schema.spectrum.libraryFolder, schema.spectrum.filterSetsFolder), ...
                fullfile(fileparts(mfilename('fullpath')), 'resources', ...
                    schema.spectrum.libraryFolder, schema.spectrum.filterSetsFolder)};
            [files, origins] = UMITRigStore.iListLibraryFiles( ...
                folders, ["user", "builtIn"], '.json');
        end

        function [files, origins] = iListLibraryFiles(folders, folderOrigins, extension)
            %ILISTLIBRARYFILES Enumerate user-first library files by stable ID.

            files = {};
            origins = strings(0, 1);
            seen = strings(0, 1);
            for iFolder = 1:numel(folders)
                if ~isfolder(folders{iFolder})
                    continue
                end
                entries = dir(fullfile(folders{iFolder}, ['*' extension]));
                for iEntry = 1:numel(entries)
                    [~, entryID] = fileparts(entries(iEntry).name);
                    id = lower(string(entryID));
                    if any(seen == id)
                        continue
                    end
                    seen(end+1, 1) = id; %#ok<AGROW>
                    files{end+1} = fullfile(entries(iEntry).folder, entries(iEntry).name); %#ok<AGROW>
                    origins(end+1, 1) = folderOrigins(iFolder); %#ok<AGROW>
                end
            end
        end

        function value = iMetadataText(metadata, fieldName, defaultValue)
            %IMETADATATEXT Read optional display metadata as a string scalar.

            value = string(defaultValue);
            if isfield(metadata, fieldName)
                candidate = metadata.(fieldName);
                if ischar(candidate) || (isstring(candidate) && isscalar(candidate))
                    value = string(candidate);
                end
            end
        end

        function filePath = iFindCaseInsensitiveFile(folders, fileName)
            filePath = '';
            for iFolder = 1:numel(folders)
                folder = folders{iFolder};
                candidate = fullfile(folder, fileName);
                if isfile(candidate)
                    filePath = candidate;
                    return
                end
                if ~isfolder(folder)
                    continue
                end
                entries = dir(folder);
                idx = find(~[entries.isdir] & strcmpi({entries.name}, fileName), 1);
                if ~isempty(idx)
                    filePath = fullfile(folder, entries(idx).name);
                    return
                end
            end
        end

        function filterSet = iValidateFilterSet(filterSet, filterSetID)
            errID = 'Umitoolbox:UMITRigStore:invalidFilterSet';
            if ~isstruct(filterSet) || ~isscalar(filterSet)
                error(errID, 'Filter set "%s" must be a scalar JSON object.', filterSetID);
            end
            required = {'id', 'displayName', 'excitationSpectrumID', 'emissionSpectrumID'};
            missing = setdiff(required, fieldnames(filterSet));
            if ~isempty(missing)
                error(errID, 'Filter set "%s" is missing field(s): %s.', ...
                    filterSetID, strjoin(missing, ', '));
            end
            for iField = 1:numel(required)
                value = filterSet.(required{iField});
                if ~(ischar(value) || (isstring(value) && isscalar(value)))
                    error(errID, 'Filter-set field "%s" must be text.', required{iField});
                end
                filterSet.(required{iField}) = char(string(value));
            end
            filterSet.id = UMITRigStore.iNormalizeManagedID( ...
                getUMITRigSchema(), filterSet.id, 'filterSetID');
            if ~strcmpi(filterSet.id, char(string(filterSetID)))
                error(errID, 'Filter-set file identity does not match its id field.');
            end
            UMITRigStore.iResolveFilterSpectrum(filterSet.excitationSpectrumID);
            UMITRigStore.iResolveFilterSpectrum(filterSet.emissionSpectrumID);
        end

        function response = iResolveFilterSpectrum(spectrumID)
            spectrumID = char(string(spectrumID));
            if isempty(spectrumID) || strcmpi(spectrumID, 'none')
                response = ones(1, 301);
                return
            end
            spectrum = UMITRigStore.getSpectrum('filter', spectrumID);
            response = spectrum.response(:)';
        end

        function definition = iLoadDefaultRigDefinition()
            definitionFile = fullfile(fileparts(mfilename('fullpath')), ...
                'resources', 'defaultRigs', 'OiS200.json');
            if ~isfile(definitionFile)
                error('Umitoolbox:UMITRigStore:missingDefaultDefinition', ...
                    'Built-in OiS200 Rig definition is missing: %s', definitionFile);
            end
            definition = jsondecode(fileread(definitionFile));
            definition.cameras = UMITRigStore.iHydrateDefaultCameras( ...
                definition.cameras);
            definition.illuminations = UMITRigStore.iHydrateDefaultIlluminations( ...
                definition.illuminations);
        end

        function cameras = iHydrateDefaultCameras(cameras)
            %IHYDRATEDEFAULTCAMERAS Fill optional hardware fields from spectra.

            fieldNames = {'displayName', 'manufacturer', 'model'};
            missingFields = ~cellfun(@(name) isfield(cameras, name), fieldNames);
            for iField = find(missingFields)
                [cameras.(fieldNames{iField})] = deal('');
            end
            if ~isfield(cameras, 'serialNumber')
                [cameras.serialNumber] = deal('');
            end
            for iCamera = 1:numel(cameras)
                spectrum = UMITRigStore.getSpectrum( ...
                    'camera', cameras(iCamera).spectrumID);
                for iField = find(missingFields)
                    fieldName = fieldNames{iField};
                    if isfield(spectrum.metadata, fieldName)
                        cameras(iCamera).(fieldName) = char(string( ...
                            spectrum.metadata.(fieldName)));
                    end
                end
            end
        end

        function illuminations = iHydrateDefaultIlluminations(illuminations)
            %IHYDRATEDEFAULTILLUMINATIONS Fill display fields from spectra.

            fieldNames = {'displayName', 'manufacturer', 'model'};
            missingFields = ~cellfun( ...
                @(name) isfield(illuminations, name), fieldNames);
            for iField = find(missingFields)
                [illuminations.(fieldNames{iField})] = deal('');
            end
            for iIllumination = 1:numel(illuminations)
                spectrum = UMITRigStore.getSpectrum( ...
                    'illumination', illuminations(iIllumination).spectrumID);
                for iField = find(missingFields)
                    fieldName = fieldNames{iField};
                    if isfield(spectrum.metadata, fieldName)
                        illuminations(iIllumination).(fieldName) = char(string( ...
                            spectrum.metadata.(fieldName)));
                    end
                end
            end
        end

        function rigUUID = iReadDefaultRigUUID()
            schema = getUMITRigSchema();
            defaultFile = fullfile(UMITRigStore.getRigsRoot(), ...
                schema.store.internalFolder, schema.store.defaultFile);
            rigUUID = '';
            if ~isfile(defaultFile)
                return
            end
            loaded = load(defaultFile, schema.store.defaultVariable, '-mat');
            if ~isfield(loaded, schema.store.defaultVariable)
                error('Umitoolbox:UMITRigStore:invalidDefaultRig', ...
                    'Default Rig metadata variable is missing.');
            end
            DefaultRig = loaded.(schema.store.defaultVariable);
            if ~isstruct(DefaultRig) || ~isscalar(DefaultRig) || ...
                    ~isfield(DefaultRig, 'rigUUID')
                error('Umitoolbox:UMITRigStore:invalidDefaultRig', ...
                    'Default Rig metadata is malformed.');
            end
            rigUUID = UMITRigStore.iNormalizeUUIDInput(DefaultRig.rigUUID);
        end

        function iWriteDefaultRigUUID(rigUUID)
            schema = getUMITRigSchema();
            internalFolder = fullfile(UMITRigStore.getRigsRoot(), schema.store.internalFolder);
            if ~isfolder(internalFolder)
                mkdir(internalFolder);
            end
            DefaultRig = struct( ...
                'rigUUID', UMITRigStore.iNormalizeUUIDInput(rigUUID), ...
                'modifiedOn', datetime('now'));
            saveMatAtomic(fullfile(internalFolder, schema.store.defaultFile), ...
                schema.store.defaultVariable, DefaultRig);
        end

        function rigUUID = iFindLegacyDefaultUUID()
            schema = getUMITRigSchema();
            entries = dir(UMITRigStore.getRigsRoot());
            entries = entries([entries.isdir] & ~ismember({entries.name}, {'.', '..', '.umit', 'library'}));
            matches = {};
            for iEntry = 1:numel(entries)
                rigFile = fullfile(entries(iEntry).folder, entries(iEntry).name, ...
                    schema.files.rigMetadata);
                if ~isfile(rigFile)
                    continue
                end
                try
                    loaded = load(rigFile, schema.metadataVariables.rig, '-mat');
                    info = loaded.(schema.metadataVariables.rig);
                    isActive = ~isfield(info, 'status') || strcmp(info.status, 'active');
                    if isfield(info, 'isDefault') && isequal(info.isDefault, true) && isActive
                        matches{end+1} = info.uuid; %#ok<AGROW>
                    end
                catch
                end
            end
            if numel(matches) > 1
                error('Umitoolbox:UMITRigStore:multipleLegacyDefaults', ...
                    'Multiple legacy Rigs are marked as default. Select one explicitly.');
            elseif isempty(matches)
                rigUUID = '';
            else
                rigUUID = matches{1};
            end
        end

        function cleanupObj = iAcquireStoreLock(operation)
            schema = getUMITRigSchema();
            internalFolder = fullfile(UMITRigStore.getRigsRoot(), schema.store.internalFolder);
            if ~isfolder(internalFolder)
                mkdir(internalFolder);
            end
            lockFolder = fullfile(internalFolder, schema.store.lockFolder);
            if isfolder(lockFolder)
                error('Umitoolbox:UMITRigStore:storeLocked', ...
                    'The Rig store is locked by another operation.');
            end
            [ok, message] = mkdir(lockFolder);
            if ~ok
                error('Umitoolbox:UMITRigStore:storeLockFailed', ...
                    'Could not acquire Rig-store lock: %s', message);
            end
            LockInfo = struct('operation', operation, 'createdOn', datetime('now'));
            save(fullfile(lockFolder, schema.store.lockMetadata), 'LockInfo', '-mat');
            cleanupObj = onCleanup(@() UMITRigStore.iRemoveFolderIfPresent(lockFolder));
        end

        function iWriteJSON(filePath, value)
            jsonText = jsonencode(value, 'PrettyPrint', true);
            folder = fileparts(filePath);
            stagedPath = [tempname(folder) '.json'];
            cleanupStage = onCleanup(@() UMITRigStore.iDeleteFileIfPresent(stagedPath));
            fid = fopen(stagedPath, 'w');
            if fid == -1
                error('Umitoolbox:UMITRigStore:fileWriteFailed', ...
                    'Could not write file: %s', filePath);
            end
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s\n', jsonText);
            clear cleanup
            [ok, message] = movefile(stagedPath, filePath, 'f');
            if ~ok
                error('Umitoolbox:UMITRigStore:fileWriteFailed', ...
                    'Could not install file "%s": %s', filePath, message);
            end
            clear cleanupStage
        end

        function metadata = iReadSpectrumMetadata(filePath)
            %IREADSPECTRUMMETADATA Decode the leading commented JSON header.

            fid = fopen(filePath, 'r');
            if fid == -1
                error('Umitoolbox:UMITRigStore:fileReadFailed', ...
                    'Could not read file: %s', filePath);
            end
            cleanup = onCleanup(@() fclose(fid));
            headerLines = strings(0, 1);
            while ~feof(fid)
                line = fgetl(fid);
                trimmed = strtrim(line);
                if startsWith(trimmed, '#')
                    content = extractAfter(string(trimmed), 1);
                    headerLines(end+1, 1) = strip(content, 'left'); %#ok<AGROW>
                elseif ~isempty(trimmed)
                    break
                end
            end

            if isempty(headerLines)
                % Read old user-library pairs during the format transition.
                legacyFile = replace(filePath, '.txt', '.json');
                if isfile(legacyFile)
                    metadata = jsondecode(fileread(legacyFile));
                else
                    metadata = struct();
                end
                return
            end

            try
                metadata = jsondecode(char(strjoin(headerLines, newline)));
            catch ME
                error('Umitoolbox:UMITRigStore:invalidSpectrumMetadata', ...
                    'Spectrum metadata header is invalid in "%s": %s', ...
                    filePath, ME.message);
            end
            if ~isstruct(metadata) || ~isscalar(metadata)
                error('Umitoolbox:UMITRigStore:invalidSpectrumMetadata', ...
                    'Spectrum metadata header must decode to a scalar structure: %s', ...
                    filePath);
            end
        end

        function iWriteSpectrumFile(filePath, wavelength, response, metadata)
            %IWRITESPECTRUMFILE Write metadata and samples as one text file.

            jsonLines = splitlines(string(jsonencode(metadata, 'PrettyPrint', true)));
            folder = fileparts(filePath);
            stagedPath = [tempname(folder) '.txt'];
            cleanupStage = onCleanup(@() UMITRigStore.iDeleteFileIfPresent(stagedPath));
            fid = fopen(stagedPath, 'w');
            if fid == -1
                error('Umitoolbox:UMITRigStore:fileWriteFailed', ...
                    'Could not write file: %s', filePath);
            end
            cleanup = onCleanup(@() fclose(fid));
            for iLine = 1:numel(jsonLines)
                fprintf(fid, '# %s\n', char(jsonLines(iLine)));
            end
            fprintf(fid, '%d\t%.17g\n', [wavelength(:), response(:)].');
            clear cleanup
            [ok, message] = movefile(stagedPath, filePath, 'f');
            if ~ok
                error('Umitoolbox:UMITRigStore:fileWriteFailed', ...
                    'Could not install spectrum "%s": %s', filePath, message);
            end
            clear cleanupStage
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

        function iDeleteFileIfPresent(path)
            if isfile(path)
                delete(path);
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
