classdef UMITRigLedger
    %UMITRIGLEDGER Folder-bound rig ledger (UMITRigLedger.json).
    %
    %   info   = UMITRigLedger.read(saveFolder)
    %   report = UMITRigLedger.validate(saveFolder)
    %   tf     = UMITRigLedger.exists(saveFolder)
    %
    %   The rig ledger is the single source of rig truth for the data of one
    %   save folder. It lives next to AcqInfos.mat and records, as an
    %   append-only event list, which library rig the folder is bound to and
    %   which content-addressed snapshot parts (optics, coregistration,
    %   filter sets) the folder has committed to. Analysis code reads rig
    %   information ONLY from the ledger: the snapshots embed the spectra
    %   content, so a folder stays self-contained on a machine whose rig
    %   library is empty. Library edits never change a ledger.
    %
    %   Folder states (derived by replaying the events, never stored):
    %       Unknown           - no ledger file
    %       NoRig             - explicit "no rig" (never commits)
    %       BoundUncommitted  - bound; nothing consumed the binding yet
    %       BoundCommitted    - at least one snapshot part was consumed
    %       Conflicted        - a rebind changed a committed part; blocks all
    %                           rig-dependent work until resolved
    %   hadConflict is a permanent mark: true once any conflict occurred.
    %
    %   Read-only operations (never write):
    %       exists, read, validate, previewRebind, filePath
    %       buildOpticsSnapshot, buildCoregistrationSnapshot,
    %       coregistrationSnapshotFromTransform, buildFilterSetSnapshot,
    %       bindingFromRig, makeBinding, sourceRigFromRigInfo,
    %       canonicalJSON, computeHash, shortHash, buildErrorMessage
    %
    %   Writers (each takes a mandatory 'Origin' Name-Value):
    %       writeNoRig, bind, rebind, revertBinding,
    %       dropCoregistrationBinding, commit, recordReimport,
    %       acknowledgeConflict, clearStaleLock
    %
    %   Guards used by rig-dependent consumers (throw, never warn):
    %       assertCanUseOptics, assertCanUseCoregistration,
    %       assertCanUseFilterSet, assertCanReimport
    %
    %   Canonical hash form (computeHash/canonicalJSON):
    %       - SHA-256 over the UTF-8 bytes of a canonical JSON text,
    %         reported as 'sha256:<64 lowercase hex>'.
    %       - The top-level 'hash' field of the hashed struct is excluded.
    %       - Objects: keys sorted by ordinal (code-point) order, written as
    %         {"key":value,...} with no whitespace.
    %       - Lists: struct arrays, cell arrays and numeric vectors are
    %         written as [v1,v2,...] whatever their orientation; numeric
    %         matrices as a list of rows; the fields 'cameras' and
    %         'illuminations' are always lists, even with one element.
    %       - Numbers: 0 is written "0"; any other finite value is rounded
    %         to 6 significant digits and written as <mantissa>e<exponent>,
    %         where the mantissa has no trailing zeros (no "." when it is an
    %         integer) and the exponent is a plain signed integer, e.g.
    %         5.21e2, -1.5e-3, 1e0. NaN and Inf are rejected.
    %       - Text: JSON string with ", \ and control characters escaped
    %         (\" \\ \uXXXX); logical scalars: true/false; empty numeric: [].
    %       Snapshot builders round every stored number to 6 significant
    %       digits, so the canonical form survives the JSON round-trip.
    %
    %   Error identifiers: Umitoolbox:RigLedger:<condition>. See
    %   IMPLEMENTATION_NOTES_rig_ledger.md for the list and messages.
    %
    %   See also UMITRigStore, getUMITRigSchema.

    properties (Constant)
        FileName = 'UMITRigLedger.json'
        LockFolderName = 'UMITRigLedger.lock'
        Format = 'umit-rig-ledger'
        SchemaVersion = 1
        % Mirrors <param.version> in Umit.prj (no runtime version API exists).
        ToolboxVersion = '2.0.0'
        HashPrefix = 'sha256:'
        SignificantDigits = 6
        Parts = {'optics', 'coregistration', 'filterSet'}
    end

    %% Read-only API
    methods (Static)
        function path = filePath(saveFolder)
            %FILEPATH Absolute path of the ledger file of SAVEFOLDER.

            path = fullfile(UMITRigLedger.iFolderText(saveFolder), UMITRigLedger.FileName);
        end

        function tf = exists(saveFolder)
            %EXISTS True when SAVEFOLDER contains a ledger file.

            tf = isfile(UMITRigLedger.filePath(saveFolder));
        end

        function info = read(saveFolder)
            %READ Replay the ledger of SAVEFOLDER. Never writes.
            %
            %   info = UMITRigLedger.read(saveFolder)
            %
            %   Returns a struct with fields folder, ledgerFile, exists,
            %   ledgerUUID, createdOn, state, binding (sourceRig,
            %   opticsHash, coregistrationHash; [] when unbound),
            %   commitments (optics, coregistration, filterSets(id,hash)),
            %   conflictParts, hadConflict, events, snapshots and
            %   hashSources. A folder without a ledger reads as 'Unknown'.
            %   An invalid ledger throws Umitoolbox:RigLedger:invalidLedger.

            folder = UMITRigLedger.iFolderText(saveFolder);
            if ~UMITRigLedger.exists(folder)
                info = UMITRigLedger.iEmptyInfo(folder);
                return
            end
            try
                ledger = UMITRigLedger.iLoad(folder);
            catch ME
                UMITRigLedger.iThrowInvalid('Reading the rig ledger', folder, {ME.message});
            end
            issues = UMITRigLedger.iIntegrityIssues(ledger);
            if ~isempty(issues)
                UMITRigLedger.iThrowInvalid('Reading the rig ledger', folder, issues);
            end
            try
                info = UMITRigLedger.iReplay(ledger, folder);
            catch ME
                UMITRigLedger.iThrowInvalid('Reading the rig ledger', folder, {ME.message});
            end
        end

        function report = validate(saveFolder)
            %VALIDATE Integrity check (schema, hashes, event consistency).
            %
            %   report = UMITRigLedger.validate(saveFolder)
            %
            %   Never throws for an invalid ledger and never writes. REPORT
            %   has fields folder, exists, isValid, issues (cellstr) and
            %   info (the read() result when valid, [] otherwise). A folder
            %   without a ledger is valid (state Unknown).

            folder = UMITRigLedger.iFolderText(saveFolder);
            report = struct('folder', folder, 'exists', false, ...
                'isValid', true, 'issues', {{}}, 'info', []);
            if ~UMITRigLedger.exists(folder)
                report.info = UMITRigLedger.iEmptyInfo(folder);
                return
            end
            report.exists = true;
            try
                ledger = UMITRigLedger.iLoad(folder);
            catch ME
                report.isValid = false;
                report.issues = {ME.message};
                return
            end
            issues = UMITRigLedger.iIntegrityIssues(ledger);
            if isempty(issues)
                try
                    report.info = UMITRigLedger.iReplay(ledger, folder);
                catch ME
                    issues = {ME.message};
                end
            end
            report.isValid = isempty(issues);
            report.issues = issues;
        end

        function conflict = previewRebind(saveFolder, binding)
            %PREVIEWREBIND Report which committed parts a rebind would change.
            %
            %   conflict = UMITRigLedger.previewRebind(saveFolder, binding)
            %
            %   Read-only. CONFLICT has fields hasConflict, parts (cellstr),
            %   committedOpticsHash, committedCoregistrationHash,
            %   requestedOpticsHash, requestedCoregistrationHash and message
            %   (human-readable, empty when there is no conflict).

            info = UMITRigLedger.read(saveFolder);
            binding = UMITRigLedger.iValidateBinding(binding);
            conflict = UMITRigLedger.iDetectConflict(info, binding);
        end

        function snapshot = findSnapshot(info, part, hash)
            %FINDSNAPSHOT Snapshot of PART with HASH from a read() result.
            %
            %   snapshot = UMITRigLedger.findSnapshot(info, part, hash)
            %
            %   Returns [] when the ledger holds no such snapshot.

            snapshot = UMITRigLedger.iFindSnapshot(info, UMITRigLedger.iPartName(part), ...
                UMITRigLedger.iHashText(hash));
        end

        function acknowledgedOn = lastAcknowledgementTime(saveFolder)
            %LASTACKNOWLEDGEMENTTIME Time of the folder's latest acknowledged conflict.
            %
            %   acknowledgedOn = UMITRigLedger.lastAcknowledgementTime(saveFolder)
            %
            %   Read-only. Returns NaT when the folder has no ledger or no
            %   acknowledgeConflict event. The time is converted from the
            %   ledger's UTC timestamp to an unzoned local datetime, the form
            %   PipelineManager uses for dataHistory timestamps. Outputs of
            %   rig-dependent steps computed before this time were computed
            %   with the optics the folder committed to before the conflict.

            acknowledgedOn = NaT;
            info = UMITRigLedger.read(saveFolder);
            for iEvent = numel(info.events):-1:1
                if strcmp(info.events{iEvent}.type, 'acknowledgeConflict')
                    acknowledgedOn = datetime(info.events{iEvent}.timestamp, ...
                        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone', 'UTC');
                    acknowledgedOn.TimeZone = 'local';
                    acknowledgedOn.TimeZone = '';
                    return
                end
            end
        end

        function snapshot = committedFilterSet(saveFolder, filterSetID)
            %COMMITTEDFILTERSET The filter-set snapshot committed for an ID.
            %
            %   snapshot = UMITRigLedger.committedFilterSet(saveFolder, filterSetID)
            %
            %   Read-only. Returns [] when the folder has no commitment for
            %   FILTERSETID (matched case-insensitively, like the library).
            %   Lets consumers analyse a folder whose filter set is no longer
            %   in the rig library.

            info = UMITRigLedger.read(saveFolder);
            snapshot = [];
            entries = info.commitments.filterSets;
            match = find(strcmpi({entries.id}, char(string(filterSetID))), 1);
            if ~isempty(match)
                snapshot = UMITRigLedger.iFindSnapshot(info, 'filterSet', entries(match).hash);
            end
        end

        %% Snapshot builders
        function snapshot = buildOpticsSnapshot(rigSource)
            %BUILDOPTICSSNAPSHOT Optics snapshot from a library rig.
            %
            %   snapshot = UMITRigLedger.buildOpticsSnapshot(rigStore)
            %   snapshot = UMITRigLedger.buildOpticsSnapshot(RigInfo)
            %
            %   Embeds the wavelength grid, each camera (index, spectrumID,
            %   model, spectrumSource, response) sorted by index and each
            %   illumination (name, spectrumID, spectrumSource, response)
            %   sorted by name. Spectra content is resolved from the
            %   library now and copied, so the snapshot no longer depends on
            %   it. Records without a spectrumID keep an empty response.

            RigInfo = UMITRigLedger.iRigInfo(rigSource);
            schema = getUMITRigSchema();

            cameras = struct('index', {}, 'spectrumID', {}, 'model', {}, ...
                'spectrumSource', {}, 'response', {});
            rigCameras = UMITRigLedger.iStructList(RigInfo, 'cameras');
            for iCamera = 1:numel(rigCameras)
                record = rigCameras(iCamera);
                [spectrumID, source, response] = UMITRigLedger.iResolveSpectrum( ...
                    'camera', UMITRigLedger.iTextField(record, 'spectrumID'));
                camera = struct( ...
                    'index', double(record.index), ...
                    'spectrumID', spectrumID, ...
                    'model', UMITRigLedger.iTextField(record, 'model'), ...
                    'spectrumSource', source, ...
                    'response', response);
                cameras(end+1) = camera; %#ok<AGROW>
            end
            if ~isempty(cameras)
                [~, order] = sort([cameras.index]);
                cameras = cameras(order);
            end

            illuminations = struct('name', {}, 'spectrumID', {}, ...
                'spectrumSource', {}, 'response', {});
            rigIlluminations = UMITRigLedger.iStructList(RigInfo, 'illuminations');
            for iIllumination = 1:numel(rigIlluminations)
                record = rigIlluminations(iIllumination);
                [spectrumID, source, response] = UMITRigLedger.iResolveSpectrum( ...
                    'illumination', UMITRigLedger.iTextField(record, 'spectrumID'));
                illumination = struct( ...
                    'name', UMITRigStore.normalizeIlluminationName(record.name), ...
                    'spectrumID', spectrumID, ...
                    'spectrumSource', source, ...
                    'response', response);
                illuminations(end+1) = illumination; %#ok<AGROW>
            end
            if ~isempty(illuminations)
                [~, order] = sort({illuminations.name});
                illuminations = illuminations(order);
            end

            snapshot = struct( ...
                'wavelengthNm', UMITRigLedger.iQuantize(schema.spectrum.wavelengthNm(:)'), ...
                'cameras', cameras, ...
                'illuminations', illuminations);
            % A 1x0 struct field would collapse the struct; keep list fields.
            snapshot = UMITRigLedger.iSetListField(snapshot, 'cameras', cameras);
            snapshot = UMITRigLedger.iSetListField(snapshot, 'illuminations', illuminations);
            snapshot.hash = UMITRigLedger.computeHash(snapshot);
        end

        function snapshot = buildCoregistrationSnapshot(rigStore)
            %BUILDCOREGISTRATIONSNAPSHOT Snapshot of a rig's active transform.
            %
            %   snapshot = UMITRigLedger.buildCoregistrationSnapshot(rigStore)
            %
            %   Returns [] when the rig has no active camera-coregistration
            %   resource. The managed file checksum is verified against the
            %   registry before the transform is copied.

            if ~isa(rigStore, 'UMITRigStore')
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'buildCoregistrationSnapshot requires a UMITRigStore object.');
            end
            resource = rigStore.getActiveCameraCoregistration();
            if isempty(resource)
                snapshot = [];
                return
            end
            resourcePath = rigStore.resolveResourcePath(resource.uuid);
            actualChecksum = computeFileChecksum(resourcePath);
            if ~strcmpi(actualChecksum, resource.checksum)
                error('Umitoolbox:RigLedger:invalidCoregistrationResource', ...
                    ['The active camera-coregistration file of rig "%s" does not match ' ...
                    'its registered checksum: %s'], rigStore.getRigInfo().rigID, resourcePath);
            end
            payload = load(resourcePath, 'tform');
            if ~isfield(payload, 'tform') || isempty(payload.tform)
                error('Umitoolbox:RigLedger:invalidCoregistrationResource', ...
                    'Camera-coregistration file has no tform: %s', resourcePath);
            end
            snapshot = UMITRigLedger.coregistrationSnapshotFromTransform( ...
                resource.uuid, resource.checksum, payload.tform);
        end

        function snapshot = coregistrationSnapshotFromTransform(resourceUUID, sourceChecksum, tform)
            %COREGISTRATIONSNAPSHOTFROMTRANSFORM Coregistration snapshot.
            %
            %   snapshot = UMITRigLedger.coregistrationSnapshotFromTransform( ...
            %       resourceUUID, sourceChecksum, tform)
            %
            %   TFORM is an affine2d/projective2d/rigid2d object (matrix T),
            %   a *tform2d object (matrix A) or a numeric 3x3 matrix. The
            %   matrix is stored as-is in the convention of its class.

            [transformClass, matrix] = UMITRigLedger.iTransformMatrix(tform);
            snapshot = struct( ...
                'resourceUUID', UMITRigLedger.iText(resourceUUID, 'resourceUUID'), ...
                'sourceChecksum', lower(UMITRigLedger.iText(sourceChecksum, 'sourceChecksum')), ...
                'transformClass', transformClass, ...
                'transformMatrix', UMITRigLedger.iQuantize(matrix));
            snapshot.hash = UMITRigLedger.computeHash(snapshot);
        end

        function snapshot = buildFilterSetSnapshot(filterSetID)
            %BUILDFILTERSETSNAPSHOT Filter-set snapshot from the library.
            %
            %   snapshot = UMITRigLedger.buildFilterSetSnapshot(filterSetID)
            %
            %   Embeds id, excitation/emission spectrum IDs and responses.
            %   An empty or 'none' spectrum ID is a flat response of ones.

            filterSet = UMITRigStore.getFilterSet(filterSetID);
            snapshot = struct( ...
                'id', filterSet.id, ...
                'excitationSpectrumID', filterSet.excitationSpectrumID, ...
                'emissionSpectrumID', filterSet.emissionSpectrumID, ...
                'excitationResponse', UMITRigLedger.iFilterResponse(filterSet.excitationSpectrumID), ...
                'emissionResponse', UMITRigLedger.iFilterResponse(filterSet.emissionSpectrumID));
            snapshot.hash = UMITRigLedger.computeHash(snapshot);
        end

        function sourceRig = sourceRigFromRigInfo(rigSource)
            %SOURCERIGFROMRIGINFO Identity of a library rig for bind events.

            RigInfo = UMITRigLedger.iRigInfo(rigSource);
            modifiedOn = '';
            if isfield(RigInfo, 'modifiedOn') && ~isempty(RigInfo.modifiedOn)
                modifiedOn = UMITRigLedger.iDateText(RigInfo.modifiedOn);
            end
            sourceRig = struct( ...
                'uuid', UMITRigLedger.iTextField(RigInfo, 'uuid'), ...
                'rigID', UMITRigLedger.iTextField(RigInfo, 'rigID'), ...
                'displayName', UMITRigLedger.iTextField(RigInfo, 'displayName'), ...
                'libraryModifiedOn', modifiedOn);
        end

        function binding = bindingFromRig(rigStore, varargin)
            %BINDINGFROMRIG Build a binding (identity + snapshots) for a rig.
            %
            %   binding = UMITRigLedger.bindingFromRig(rigStore)
            %   binding = UMITRigLedger.bindingFromRig(rigStore, ...
            %       'IncludeCoregistration', false)
            %
            %   IncludeCoregistration (default true) snapshots the rig's
            %   active coregistration resource; [] when it has none.

            p = inputParser;
            p.FunctionName = 'bindingFromRig';
            addParameter(p, 'IncludeCoregistration', true, @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});
            if ~isa(rigStore, 'UMITRigStore')
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'bindingFromRig requires a UMITRigStore object.');
            end
            coregistration = [];
            if p.Results.IncludeCoregistration
                coregistration = UMITRigLedger.buildCoregistrationSnapshot(rigStore);
            end
            binding = UMITRigLedger.makeBinding( ...
                UMITRigLedger.sourceRigFromRigInfo(rigStore), ...
                UMITRigLedger.buildOpticsSnapshot(rigStore), coregistration);
        end

        function binding = makeBinding(sourceRig, opticsSnapshot, coregistrationSnapshot)
            %MAKEBINDING Assemble a binding from a rig identity and snapshots.
            %
            %   COREGISTRATIONSNAPSHOT may be [] (no coregistration).

            if nargin < 3
                coregistrationSnapshot = [];
            end
            binding = struct('sourceRig', sourceRig, ...
                'optics', opticsSnapshot, ...
                'coregistration', []);
            binding.coregistration = coregistrationSnapshot;
            binding = UMITRigLedger.iValidateBinding(binding);
        end

        %% Canonical serialization and hashing
        function text = canonicalJSON(value)
            %CANONICALJSON Canonical JSON text used for hashing.
            %
            %   See the class help for the exact canonical form. The
            %   top-level 'hash' field of a struct is excluded.

            if isstruct(value) && isscalar(value) && isfield(value, 'hash')
                value = rmfield(value, 'hash');
            end
            text = UMITRigLedger.iCanonical(value, false);
        end

        function hash = computeHash(value)
            %COMPUTEHASH 'sha256:<hex>' of canonicalJSON(VALUE).

            text = UMITRigLedger.canonicalJSON(value);
            digest = java.security.MessageDigest.getInstance('SHA-256');
            bytes = unicode2native(text, 'UTF-8');
            digest.update(typecast(uint8(bytes(:)), 'int8'));
            digestBytes = typecast(digest.digest(), 'uint8');
            hash = [UMITRigLedger.HashPrefix ...
                lower(reshape(dec2hex(digestBytes, 2).', 1, []))];
        end

        function text = shortHash(hash)
            %SHORTHASH 'sha256:' plus the first 12 hex digits, or '(none)'.

            hash = UMITRigLedger.iHashText(hash);
            if isempty(hash)
                text = '(none)';
            else
                text = hash(1:min(numel(hash), numel(UMITRigLedger.HashPrefix) + 12));
            end
        end

        function message = buildErrorMessage(operation, saveFolder, why, fixes)
            %BUILDERRORMESSAGE Shared error text: what, where, why, how to fix.
            %
            %   message = UMITRigLedger.buildErrorMessage(operation, ...
            %       saveFolder, why, fixes)
            %
            %   FIXES is a cellstr; each entry becomes one numbered fix.

            if ischar(fixes) || isstring(fixes)
                fixes = cellstr(fixes);
            end
            lines = cell(1, numel(fixes));
            for iFix = 1:numel(fixes)
                lines{iFix} = sprintf('  (%d) %s', iFix, fixes{iFix});
            end
            message = sprintf('%s was refused for folder "%s".\nWhy: %s\nHow to fix:\n%s', ...
                char(string(operation)), char(string(saveFolder)), ...
                char(string(why)), strjoin(lines, newline));
        end

        %% Guards
        function result = assertCanUseOptics(saveFolder, varargin)
            %ASSERTCANUSEOPTICS Guard for optics consumers (e.g. HemoCompute).
            %
            %   result = UMITRigLedger.assertCanUseOptics(saveFolder, ...
            %       'Operation', 'HemoCompute')
            %
            %   Throws when the folder has no ledger, is NoRig, is
            %   Conflicted or has an invalid ledger. RESULT has fields
            %   snapshot (optics), hash, sourceRig and needsCommit.

            operation = UMITRigLedger.iParseOperation('assertCanUseOptics', varargin{:});
            info = UMITRigLedger.iGuardBase(saveFolder, operation);
            result = struct( ...
                'snapshot', UMITRigLedger.iFindSnapshot(info, 'optics', info.binding.opticsHash), ...
                'hash', info.binding.opticsHash, ...
                'sourceRig', info.binding.sourceRig, ...
                'needsCommit', isempty(info.commitments.optics));
        end

        function result = assertCanUseCoregistration(saveFolder, varargin)
            %ASSERTCANUSECOREGISTRATION Guard for coregistration consumers.
            %
            %   result = UMITRigLedger.assertCanUseCoregistration(saveFolder, ...
            %       'Operation', name, 'RequestedHash', hash)
            %
            %   RequestedHash (optional) is the coregistration snapshot hash
            %   the caller intends to apply. It must equal the bound
            %   coregistration; applying a transform other than the
            %   committed one is a conflict. RESULT has fields snapshot,
            %   hash, sourceRig and needsCommit.

            p = inputParser;
            p.FunctionName = 'assertCanUseCoregistration';
            addParameter(p, 'Operation', 'Using the camera coregistration', @UMITRigLedger.iIsText);
            addParameter(p, 'RequestedHash', '', @(x) isempty(x) || UMITRigLedger.iIsText(x));
            parse(p, varargin{:});
            operation = char(string(p.Results.Operation));
            requested = UMITRigLedger.iHashText(p.Results.RequestedHash);

            info = UMITRigLedger.iGuardBase(saveFolder, operation);
            bound = info.binding.coregistrationHash;
            committed = info.commitments.coregistration;
            if isempty(bound)
                UMITRigLedger.iThrow('Umitoolbox:RigLedger:noCoregistrationBound', ...
                    operation, info.folder, sprintf( ...
                    'the folder is bound to %s without a camera coregistration.', ...
                    UMITRigLedger.iDescribeRig(info.binding.sourceRig)), ...
                    {UMITRigLedger.iFixRebind(), UMITRigLedger.iFixFreshFolder()});
            end
            if ~isempty(requested) && ~strcmp(requested, bound)
                if ~isempty(committed)
                    UMITRigLedger.iThrow('Umitoolbox:RigLedger:coregistrationConflict', ...
                        operation, info.folder, sprintf( ...
                        ['the folder already used the camera coregistration of %s, but a ' ...
                        'different transform was requested. A different transform can only ' ...
                        'be applied through a re-import.'], ...
                        UMITRigLedger.iDescribeHashSource(info, 'coregistration', committed)), ...
                        {UMITRigLedger.iFixUseCommitted('coregistration'), ...
                        UMITRigLedger.iFixFreshFolder()});
                end
                UMITRigLedger.iThrow('Umitoolbox:RigLedger:coregistrationNotBound', ...
                    operation, info.folder, sprintf( ...
                    ['the folder is bound to the camera coregistration of %s, but a ' ...
                    'different transform was requested.'], ...
                    UMITRigLedger.iDescribeHashSource(info, 'coregistration', bound)), ...
                    {['Rebind the folder to the rig whose transform you want to apply ' ...
                    '(rebinding is free while the coregistration is uncommitted).'], ...
                    UMITRigLedger.iFixFreshFolder()});
            end
            result = struct( ...
                'snapshot', UMITRigLedger.iFindSnapshot(info, 'coregistration', bound), ...
                'hash', bound, ...
                'sourceRig', info.binding.sourceRig, ...
                'needsCommit', isempty(committed));
        end

        function assertCanReimport(saveFolder, binding, varargin)
            %ASSERTCANREIMPORT Preflight a raw re-import that binds BINDING.
            %
            %   UMITRigLedger.assertCanReimport(saveFolder, binding, ...
            %       'Operation', name)
            %
            %   Read-only. A re-import resets the coregistration commitment
            %   but keeps the optics commitment (brief section 1.6), so the
            %   re-import is refused when the folder committed to optics
            %   other than BINDING.optics. An empty BINDING (no rig) or a
            %   folder without a ledger always passes. An invalid ledger
            %   throws invalidLedger.

            operation = UMITRigLedger.iParseOperation('assertCanReimport', varargin{:});
            info = UMITRigLedger.read(saveFolder);
            if isempty(binding) || isempty(info.commitments.optics)
                return
            end
            binding = UMITRigLedger.iValidateBinding(binding);
            committed = info.commitments.optics;
            if strcmp(committed, binding.optics.hash)
                return
            end
            UMITRigLedger.iThrow('Umitoolbox:RigLedger:reimportConflict', ...
                operation, info.folder, sprintf( ...
                ['the folder already used the optics of %s, but the re-import would ' ...
                'use %s, whose optics differ. A re-import resets only the coregistration ' ...
                'commitment, never the optics commitment.'], ...
                UMITRigLedger.iDescribeHashSource(info, 'optics', committed), ...
                UMITRigLedger.iDescribeRig(binding.sourceRig)), ...
                {sprintf(['Activate the library rig whose optics the folder used (%s) ' ...
                'and re-import.'], ...
                UMITRigLedger.iDescribeHashSource(info, 'optics', committed)), ...
                UMITRigLedger.iFixFreshFolder(), ...
                ['To accept the new optics deliberately, rebind the folder to the new rig ' ...
                'and acknowledge the conflict (DataViewer: Rig > Bind/Rebind, then ' ...
                'Acknowledge Conflict), then re-import.']});
        end

        function needsCommit = assertCanUseFilterSet(saveFolder, filterSetSnapshot, varargin)
            %ASSERTCANUSEFILTERSET Guard for filter-set consumers.
            %
            %   needsCommit = UMITRigLedger.assertCanUseFilterSet(saveFolder, ...
            %       filterSetSnapshot, 'Operation', name)
            %
            %   Rules: an ID not yet committed is allowed (needsCommit true);
            %   an ID committed with the same hash is allowed (needsCommit
            %   false); an ID committed with a different hash throws
            %   Umitoolbox:RigLedger:filterSetModified.

            operation = UMITRigLedger.iParseOperation('assertCanUseFilterSet', varargin{:});
            snapshot = UMITRigLedger.iValidateSnapshot('filterSet', filterSetSnapshot);
            info = UMITRigLedger.iGuardBase(saveFolder, operation);
            committedHash = UMITRigLedger.iCommittedFilterSetHash(info, snapshot.id);
            needsCommit = isempty(committedHash);
            if ~needsCommit && ~strcmp(committedHash, snapshot.hash)
                UMITRigLedger.iThrowFilterSetModified(operation, info, snapshot, committedHash);
            end
        end
    end

    %% Writers
    methods (Static)
        function wasWritten = writeNoRig(saveFolder, varargin)
            %WRITENORIG Record an explicit "no rig" for a folder.
            %
            %   wasWritten = UMITRigLedger.writeNoRig(saveFolder, 'Origin', origin)
            %
            %   Creates the ledger with a noRig event. A folder that is
            %   already NoRig is left unchanged (returns false). A bound
            %   folder throws Umitoolbox:RigLedger:alreadyBound.

            origin = UMITRigLedger.iParseOrigin('writeNoRig', varargin{:});
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'writeNoRig'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Recording "no rig"');
            switch info.state
                case 'NoRig'
                    wasWritten = false;
                    return
                case 'Unknown'
                otherwise
                    UMITRigLedger.iThrow('Umitoolbox:RigLedger:alreadyBound', ...
                        'Recording "no rig"', folder, sprintf( ...
                        'the folder is already bound to %s; a binding cannot be replaced by "no rig".', ...
                        UMITRigLedger.iDescribeRig(info.binding.sourceRig)), ...
                        {'Keep the existing binding.', UMITRigLedger.iFixFreshFolder()});
            end
            ledger = UMITRigLedger.iAppendEvent(ledger, ...
                UMITRigLedger.iNewEvent('noRig', origin));
            UMITRigLedger.iWrite(folder, ledger);
            wasWritten = true;
        end

        function bind(saveFolder, binding, varargin)
            %BIND Bind an Unknown or NoRig folder to a rig snapshot.
            %
            %   UMITRigLedger.bind(saveFolder, binding, 'Origin', origin)
            %
            %   BINDING comes from bindingFromRig or makeBinding. A folder
            %   that is already bound throws Umitoolbox:RigLedger:alreadyBound
            %   (use rebind).

            origin = UMITRigLedger.iParseOrigin('bind', varargin{:});
            binding = UMITRigLedger.iValidateBinding(binding);
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'bind'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Binding a rig');
            if ~ismember(info.state, {'Unknown', 'NoRig'})
                UMITRigLedger.iThrow('Umitoolbox:RigLedger:alreadyBound', ...
                    'Binding a rig', folder, sprintf('the folder is already bound to %s.', ...
                    UMITRigLedger.iDescribeRig(info.binding.sourceRig)), ...
                    {'Use UMITRigLedger.rebind (DataViewer: Rig > Bind/Rebind) to change the binding.'});
            end
            ledger = UMITRigLedger.iAddBindingSnapshots(ledger, binding);
            event = UMITRigLedger.iNewEvent('bind', origin);
            event.sourceRig = binding.sourceRig;
            event.opticsHash = binding.optics.hash;
            event.coregistrationHash = UMITRigLedger.iNullable(UMITRigLedger.iBindingHash(binding, 'coregistration'));
            ledger = UMITRigLedger.iAppendEvent(ledger, event);
            UMITRigLedger.iWrite(folder, ledger);
        end

        function conflict = rebind(saveFolder, binding, varargin)
            %REBIND Change the binding of a bound folder.
            %
            %   conflict = UMITRigLedger.rebind(saveFolder, binding, 'Origin', origin)
            %
            %   Writes a rebind event. When the new binding changes a
            %   committed part, a conflict event is also written and the
            %   folder becomes Conflicted; CONFLICT (see previewRebind)
            %   describes it. Rebinding to an identical binding writes
            %   nothing.

            origin = UMITRigLedger.iParseOrigin('rebind', varargin{:});
            binding = UMITRigLedger.iValidateBinding(binding);
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'rebind'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Rebinding a rig');
            conflict = UMITRigLedger.iRebindLocked(folder, ledger, info, binding, origin);
        end

        function conflict = revertBinding(saveFolder, varargin)
            %REVERTBINDING Rebind to the committed snapshots (conflict fix 1).
            %
            %   conflict = UMITRigLedger.revertBinding(saveFolder, 'Origin', origin)
            %
            %   Uses the snapshots stored in the ledger, so no library rig
            %   is needed. Committed parts are restored; uncommitted parts
            %   keep their current binding.

            origin = UMITRigLedger.iParseOrigin('revertBinding', varargin{:});
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'revertBinding'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Reverting the rig binding');
            if isempty(info.binding)
                UMITRigLedger.iThrowNotBound('Reverting the rig binding', info);
            end
            opticsHash = info.binding.opticsHash;
            if ~isempty(info.commitments.optics)
                opticsHash = info.commitments.optics;
            end
            coregistrationHash = info.binding.coregistrationHash;
            if ~isempty(info.commitments.coregistration)
                coregistrationHash = info.commitments.coregistration;
            end
            coregistration = [];
            if ~isempty(coregistrationHash)
                coregistration = UMITRigLedger.iFindSnapshot(info, 'coregistration', coregistrationHash);
            end
            sourceRig = UMITRigLedger.iHashSourceRig(info, 'optics', opticsHash);
            binding = UMITRigLedger.makeBinding(sourceRig, ...
                UMITRigLedger.iFindSnapshot(info, 'optics', opticsHash), coregistration);
            conflict = UMITRigLedger.iRebindLocked(folder, ledger, info, binding, origin);
        end

        function conflict = dropCoregistrationBinding(saveFolder, varargin)
            %DROPCOREGISTRATIONBINDING Keep the bound optics, bind no transform.
            %
            %   conflict = UMITRigLedger.dropCoregistrationBinding(saveFolder, ...
            %       'Origin', origin)
            %
            %   Rebinds to the folder's own optics snapshot and source rig
            %   with no coregistration (used when re-imported dual-camera
            %   data could not be coregistered). Writes nothing when no
            %   coregistration is bound. Dropping a committed coregistration
            %   is a conflict, like any rebind.

            origin = UMITRigLedger.iParseOrigin('dropCoregistrationBinding', varargin{:});
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'dropCoregistrationBinding'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Dropping the coregistration binding');
            if isempty(info.binding)
                UMITRigLedger.iThrowNotBound('Dropping the coregistration binding', info);
            end
            binding = UMITRigLedger.makeBinding(info.binding.sourceRig, ...
                UMITRigLedger.iFindSnapshot(info, 'optics', info.binding.opticsHash), []);
            conflict = UMITRigLedger.iDetectConflict(info, binding);
            if isempty(info.binding.coregistrationHash)
                return
            end
            conflict = UMITRigLedger.iRebindLocked(folder, ledger, info, binding, origin);
        end

        function wasWritten = commit(saveFolder, part, value, varargin)
            %COMMIT Record that a snapshot part was consumed by the folder.
            %
            %   wasWritten = UMITRigLedger.commit(saveFolder, 'optics', hash, ...
            %       'Origin', origin, 'FirstUsedBy', name)
            %   wasWritten = UMITRigLedger.commit(saveFolder, 'filterSet', ...
            %       filterSetSnapshot, 'Origin', origin, 'FirstUsedBy', name)
            %
            %   For optics/coregistration VALUE is the bound hash (or its
            %   snapshot). For filterSet VALUE is the snapshot. Commits
            %   already present are not rewritten (returns false). The
            %   folder must be bound and not Conflicted.

            p = inputParser;
            p.FunctionName = 'commit';
            addParameter(p, 'Origin', '', @UMITRigLedger.iIsText);
            addParameter(p, 'FirstUsedBy', '', @UMITRigLedger.iIsText);
            parse(p, varargin{:});
            origin = UMITRigLedger.iRequireOrigin(p.Results.Origin);
            firstUsedBy = char(string(p.Results.FirstUsedBy));
            if isempty(firstUsedBy)
                firstUsedBy = origin;
            end
            part = UMITRigLedger.iPartName(part);
            operation = sprintf('Committing the %s snapshot', part);

            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'commit'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, operation);
            UMITRigLedger.iAssertUsableState(info, operation);

            event = UMITRigLedger.iNewEvent('commit', origin);
            event.part = part;
            if strcmp(part, 'filterSet')
                snapshot = UMITRigLedger.iValidateSnapshot('filterSet', value);
                committedHash = UMITRigLedger.iCommittedFilterSetHash(info, snapshot.id);
                if strcmp(committedHash, snapshot.hash)
                    wasWritten = false;
                    return
                elseif ~isempty(committedHash)
                    UMITRigLedger.iThrowFilterSetModified(operation, info, snapshot, committedHash);
                end
                ledger = UMITRigLedger.iAddSnapshot(ledger, 'filterSet', snapshot);
                event.hash = snapshot.hash;
                event.filterSetID = snapshot.id;
            else
                if isstruct(value)
                    value = UMITRigLedger.iValidateSnapshot(part, value).hash;
                end
                hash = UMITRigLedger.iHashText(value);
                boundHash = info.binding.([part 'Hash']);
                if isempty(boundHash)
                    UMITRigLedger.iThrow('Umitoolbox:RigLedger:noCoregistrationBound', ...
                        operation, folder, 'the folder binding has no camera coregistration.', ...
                        {UMITRigLedger.iFixRebind()});
                end
                if ~strcmp(hash, boundHash)
                    UMITRigLedger.iThrow('Umitoolbox:RigLedger:commitMismatch', ...
                        operation, folder, sprintf( ...
                        'the folder is bound to the %s of %s, but a different %s was used.', ...
                        UMITRigLedger.iPartLabel(part), ...
                        UMITRigLedger.iDescribeHashSource(info, part, boundHash), ...
                        UMITRigLedger.iPartLabel(part)), ...
                        {sprintf('Use the %s snapshot bound in the folder ledger.', part), ...
                        UMITRigLedger.iFixRebind()});
                end
                if strcmp(info.commitments.(part), hash)
                    wasWritten = false;
                    return
                end
                event.hash = hash;
                event.filterSetID = UMITRigLedger.iNullable('');
            end
            event.firstUsedBy = firstUsedBy;
            ledger = UMITRigLedger.iAppendEvent(ledger, event);
            UMITRigLedger.iWrite(folder, ledger);
            wasWritten = true;
        end

        function recordReimport(saveFolder, varargin)
            %RECORDREIMPORT Record a raw re-import; resets coregistration.
            %
            %   UMITRigLedger.recordReimport(saveFolder, 'Origin', origin)
            %
            %   Optics and filter-set commitments are kept. The folder
            %   must already have a ledger.

            origin = UMITRigLedger.iParseOrigin('recordReimport', varargin{:});
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'recordReimport'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Recording a re-import');
            if strcmp(info.state, 'Unknown')
                UMITRigLedger.iThrowNoLedger('Recording a re-import', folder);
            end
            event = UMITRigLedger.iNewEvent('reimport', origin);
            event.previousCoregistrationHash = UMITRigLedger.iNullable(info.commitments.coregistration);
            ledger = UMITRigLedger.iAppendEvent(ledger, event);
            UMITRigLedger.iWrite(folder, ledger);
        end

        function acknowledgeConflict(saveFolder, varargin)
            %ACKNOWLEDGECONFLICT Accept the current binding of a Conflicted folder.
            %
            %   UMITRigLedger.acknowledgeConflict(saveFolder, 'Origin', origin, ...
            %       'Comment', text)
            %
            %   The conflicting parts are committed to the current binding.
            %   Outputs computed before the rebind are NOT recomputed;
            %   hadConflict stays true permanently.

            p = inputParser;
            p.FunctionName = 'acknowledgeConflict';
            addParameter(p, 'Origin', '', @UMITRigLedger.iIsText);
            addParameter(p, 'Comment', '', @UMITRigLedger.iIsText);
            parse(p, varargin{:});
            origin = UMITRigLedger.iRequireOrigin(p.Results.Origin);

            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lock = UMITRigLedger.iAcquireLock(folder, 'acknowledgeConflict'); %#ok<NASGU>
            [ledger, info] = UMITRigLedger.iLoadForWrite(folder, 'Acknowledging a rig conflict');
            if ~strcmp(info.state, 'Conflicted')
                UMITRigLedger.iThrow('Umitoolbox:RigLedger:noConflict', ...
                    'Acknowledging a rig conflict', folder, sprintf( ...
                    'the folder state is %s; there is no open conflict.', info.state), ...
                    {'Nothing to do.'});
            end
            event = UMITRigLedger.iNewEvent('acknowledgeConflict', origin);
            event.comment = char(string(p.Results.Comment));
            event.parts = info.conflictParts;
            event.acceptedOpticsHash = info.binding.opticsHash;
            event.acceptedCoregistrationHash = UMITRigLedger.iNullable(info.binding.coregistrationHash);
            ledger = UMITRigLedger.iAppendEvent(ledger, event);
            UMITRigLedger.iWrite(folder, ledger);
        end

        function wasCleared = clearStaleLock(saveFolder, varargin)
            %CLEARSTALELOCK Remove a leftover ledger lock older than MinimumAge.
            %
            %   wasCleared = UMITRigLedger.clearStaleLock(saveFolder, ...
            %       'MinimumAge', minutes(10))
            %
            %   Returns false when no lock exists. A younger lock throws
            %   Umitoolbox:RigLedger:lockNotStale.

            p = inputParser;
            p.FunctionName = 'clearStaleLock';
            addParameter(p, 'MinimumAge', minutes(10), @(x) isduration(x) && isscalar(x));
            parse(p, varargin{:});
            folder = UMITRigLedger.iExistingFolder(saveFolder);
            lockFolder = fullfile(folder, UMITRigLedger.LockFolderName);
            wasCleared = false;
            if ~isfolder(lockFolder)
                return
            end
            age = UMITRigLedger.iLockAge(lockFolder);
            if age < p.Results.MinimumAge
                error('Umitoolbox:RigLedger:lockNotStale', ...
                    'The rig-ledger lock of "%s" is only %s old; it may belong to a running operation.', ...
                    folder, char(age, 'hh:mm:ss'));
            end
            rmdir(lockFolder, 's');
            wasCleared = true;
        end
    end

    %% Private helpers
    methods (Static, Access = private)
        %% Replay
        function info = iEmptyInfo(folder)
            info = struct( ...
                'folder', folder, ...
                'ledgerFile', fullfile(folder, UMITRigLedger.FileName), ...
                'exists', false, ...
                'ledgerUUID', '', ...
                'createdOn', '', ...
                'state', 'Unknown', ...
                'binding', [], ...
                'commitments', struct('optics', '', 'coregistration', '', ...
                    'filterSets', struct('id', {}, 'hash', {})), ...
                'conflictParts', {{}}, ...
                'hadConflict', false, ...
                'events', {{}}, ...
                'snapshots', struct('optics', {{}}, 'coregistration', {{}}, 'filterSet', {{}}), ...
                'hashSources', struct('part', {}, 'hash', {}, 'sourceRig', {}));
        end

        function info = iReplay(ledger, folder)
            %IREPLAY Derive the folder state from the event list. Throws a
            %plain message on any event-consistency violation.

            info = UMITRigLedger.iEmptyInfo(folder);
            info.exists = true;
            info.ledgerUUID = ledger.ledgerUUID;
            info.createdOn = ledger.createdOn;
            info.events = ledger.events;
            info.snapshots = ledger.snapshots;

            isNoRig = false;
            pendingConflict = {};
            seenUUIDs = {};
            for iEvent = 1:numel(ledger.events)
                event = ledger.events{iEvent};
                where = sprintf('Event %d (%s)', iEvent, event.type);
                if ismember(event.eventUUID, seenUUIDs)
                    error('%s reuses eventUUID %s.', where, event.eventUUID);
                end
                seenUUIDs{end+1} = event.eventUUID; %#ok<AGROW>
                if ~isempty(pendingConflict) && ~strcmp(event.type, 'conflict')
                    error('%s follows a conflicting rebind without its conflict event.', where);
                end

                switch event.type
                    case 'noRig'
                        if ~isempty(info.binding)
                            error('%s appears after a binding.', where);
                        end
                        isNoRig = true;

                    case {'bind', 'rebind'}
                        isRebind = strcmp(event.type, 'rebind');
                        if isRebind && isempty(info.binding)
                            error('%s has no binding to replace.', where);
                        elseif ~isRebind && ~isempty(info.binding)
                            error('%s appears after a binding (expected rebind).', where);
                        end
                        UMITRigLedger.iRequireFields(event, ...
                            {'sourceRig', 'opticsHash', 'coregistrationHash'}, where);
                        opticsHash = UMITRigLedger.iHashText(event.opticsHash);
                        coregistrationHash = UMITRigLedger.iHashText(event.coregistrationHash);
                        UMITRigLedger.iRequireSnapshot(info, 'optics', opticsHash, where);
                        if ~isempty(coregistrationHash)
                            UMITRigLedger.iRequireSnapshot(info, 'coregistration', coregistrationHash, where);
                        end
                        if isRebind
                            UMITRigLedger.iRequireFields(event, ...
                                {'previousOpticsHash', 'previousCoregistrationHash'}, where);
                            if ~strcmp(UMITRigLedger.iHashText(event.previousOpticsHash), info.binding.opticsHash) || ...
                                    ~strcmp(UMITRigLedger.iHashText(event.previousCoregistrationHash), ...
                                    info.binding.coregistrationHash)
                                error('%s does not match the binding it replaces.', where);
                            end
                        end
                        sourceRig = UMITRigLedger.iNormalizeSourceRig(event.sourceRig, where);
                        info.binding = struct('sourceRig', sourceRig, ...
                            'opticsHash', opticsHash, 'coregistrationHash', coregistrationHash);
                        info = UMITRigLedger.iNoteHashSource(info, 'optics', opticsHash, sourceRig);
                        if ~isempty(coregistrationHash)
                            info = UMITRigLedger.iNoteHashSource(info, 'coregistration', ...
                                coregistrationHash, sourceRig);
                        end
                        isNoRig = false;
                        if isRebind
                            pendingConflict = UMITRigLedger.iConflictParts(info);
                        end

                    case 'commit'
                        UMITRigLedger.iRequireFields(event, {'part', 'hash', 'firstUsedBy'}, where);
                        if isempty(info.binding)
                            error('%s appears while the folder is not bound.', where);
                        end
                        if ~isempty(UMITRigLedger.iConflictParts(info))
                            error('%s appears while the folder is Conflicted.', where);
                        end
                        hash = UMITRigLedger.iHashText(event.hash);
                        switch event.part
                            case {'optics', 'coregistration'}
                                if ~strcmp(hash, info.binding.([event.part 'Hash']))
                                    error('%s commits a %s hash that is not bound.', where, event.part);
                                end
                                if ~isempty(info.commitments.(event.part))
                                    error('%s duplicates an existing %s commitment.', where, event.part);
                                end
                                info.commitments.(event.part) = hash;
                            case 'filterSet'
                                UMITRigLedger.iRequireFields(event, {'filterSetID'}, where);
                                snapshot = UMITRigLedger.iRequireSnapshot(info, 'filterSet', hash, where);
                                if ~strcmp(snapshot.id, event.filterSetID)
                                    error('%s filterSetID does not match its snapshot.', where);
                                end
                                if ~isempty(UMITRigLedger.iCommittedFilterSetHash(info, snapshot.id))
                                    error('%s commits filter set "%s" twice.', where, snapshot.id);
                                end
                                info.commitments.filterSets(end+1) = ...
                                    struct('id', snapshot.id, 'hash', hash);
                            otherwise
                                error('%s has unknown part "%s".', where, event.part);
                        end

                    case 'reimport'
                        info.commitments.coregistration = '';

                    case 'conflict'
                        UMITRigLedger.iRequireFields(event, {'parts'}, where);
                        if isempty(pendingConflict) || ...
                                ~isequal(sort(cellstr(event.parts(:)')), sort(pendingConflict))
                            error('%s does not match a preceding conflicting rebind.', where);
                        end
                        pendingConflict = {};
                        info.hadConflict = true;

                    case 'acknowledgeConflict'
                        UMITRigLedger.iRequireFields(event, {'comment'}, where);
                        parts = UMITRigLedger.iConflictParts(info);
                        if isempty(parts)
                            error('%s has no open conflict to acknowledge.', where);
                        end
                        for iPart = 1:numel(parts)
                            info.commitments.(parts{iPart}) = info.binding.([parts{iPart} 'Hash']);
                        end

                    otherwise
                        error('%s has an unknown event type.', where);
                end
            end
            if ~isempty(pendingConflict)
                error('The last rebind changed committed parts but has no conflict event.');
            end

            info.conflictParts = UMITRigLedger.iConflictParts(info);
            if isempty(info.binding)
                if ~isNoRig
                    error('The ledger has no noRig or bind event.');
                end
                info.state = 'NoRig';
            elseif ~isempty(info.conflictParts)
                info.state = 'Conflicted';
            elseif ~isempty(info.commitments.optics) || ...
                    ~isempty(info.commitments.coregistration) || ...
                    ~isempty(info.commitments.filterSets)
                info.state = 'BoundCommitted';
            else
                info.state = 'BoundUncommitted';
            end
        end

        function parts = iConflictParts(info)
            parts = {};
            if isempty(info.binding)
                return
            end
            names = {'optics', 'coregistration'};
            for iPart = 1:numel(names)
                committed = info.commitments.(names{iPart});
                if ~isempty(committed) && ~strcmp(committed, info.binding.([names{iPart} 'Hash']))
                    parts{end+1} = names{iPart}; %#ok<AGROW>
                end
            end
        end

        function conflict = iDetectConflict(info, binding)
            requestedOptics = binding.optics.hash;
            requestedCoreg = UMITRigLedger.iBindingHash(binding, 'coregistration');
            committedOptics = info.commitments.optics;
            committedCoreg = info.commitments.coregistration;
            parts = {};
            if ~isempty(committedOptics) && ~strcmp(committedOptics, requestedOptics)
                parts{end+1} = 'optics';
            end
            if ~isempty(committedCoreg) && ~strcmp(committedCoreg, requestedCoreg)
                parts{end+1} = 'coregistration';
            end
            conflict = struct( ...
                'hasConflict', ~isempty(parts), ...
                'parts', {parts}, ...
                'committedOpticsHash', committedOptics, ...
                'committedCoregistrationHash', committedCoreg, ...
                'requestedOpticsHash', requestedOptics, ...
                'requestedCoregistrationHash', requestedCoreg, ...
                'message', '');
            if conflict.hasConflict
                lines = cell(1, numel(parts));
                for iPart = 1:numel(parts)
                    part = parts{iPart};
                    label = UMITRigLedger.iPartLabel(part);
                    committedRig = UMITRigLedger.iHashSourceRig(info, part, ...
                        conflict.(['committed' UMITRigLedger.iCap(part) 'Hash']));
                    if ~isempty(committedRig) && strcmp(committedRig.uuid, binding.sourceRig.uuid)
                        lines{iPart} = sprintf(['%s: the %s of %s have changed since the ' ...
                            'folder used them.'], label, label, UMITRigLedger.iDescribeRig(committedRig));
                    else
                        lines{iPart} = sprintf(['%s: the folder already used the %s of %s; ' ...
                            'the %s of %s %s.'], label, label, ...
                            UMITRigLedger.iDescribeRig(committedRig), label, ...
                            UMITRigLedger.iDescribeRig(binding.sourceRig), ...
                            UMITRigLedger.iDiffersVerb(part));
                    end
                end
                conflict.message = strjoin(lines, newline);
            end
        end

        function conflict = iRebindLocked(folder, ledger, info, binding, origin)
            operation = 'Rebinding a rig';
            if isempty(info.binding)
                UMITRigLedger.iThrowNotBound(operation, info);
            end
            conflict = UMITRigLedger.iDetectConflict(info, binding);
            newCoreg = UMITRigLedger.iBindingHash(binding, 'coregistration');
            if strcmp(binding.optics.hash, info.binding.opticsHash) && ...
                    strcmp(newCoreg, info.binding.coregistrationHash) && ...
                    strcmp(binding.sourceRig.uuid, info.binding.sourceRig.uuid)
                return
            end
            ledger = UMITRigLedger.iAddBindingSnapshots(ledger, binding);
            event = UMITRigLedger.iNewEvent('rebind', origin);
            event.sourceRig = binding.sourceRig;
            event.opticsHash = binding.optics.hash;
            event.coregistrationHash = UMITRigLedger.iNullable(newCoreg);
            event.previousOpticsHash = info.binding.opticsHash;
            event.previousCoregistrationHash = UMITRigLedger.iNullable(info.binding.coregistrationHash);
            ledger = UMITRigLedger.iAppendEvent(ledger, event);
            if conflict.hasConflict
                conflictEvent = UMITRigLedger.iNewEvent('conflict', origin);
                conflictEvent.parts = conflict.parts;
                conflictEvent.committedOpticsHash = UMITRigLedger.iNullable(conflict.committedOpticsHash);
                conflictEvent.committedCoregistrationHash = UMITRigLedger.iNullable(conflict.committedCoregistrationHash);
                conflictEvent.requestedOpticsHash = conflict.requestedOpticsHash;
                conflictEvent.requestedCoregistrationHash = UMITRigLedger.iNullable(conflict.requestedCoregistrationHash);
                ledger = UMITRigLedger.iAppendEvent(ledger, conflictEvent);
            end
            UMITRigLedger.iWrite(folder, ledger);
        end

        %% Guards
        function info = iGuardBase(saveFolder, operation)
            folder = UMITRigLedger.iFolderText(saveFolder);
            report = UMITRigLedger.validate(folder);
            if ~report.isValid
                UMITRigLedger.iThrowInvalid(operation, folder, report.issues);
            end
            info = report.info;
            UMITRigLedger.iAssertUsableState(info, operation);
        end

        function iAssertUsableState(info, operation)
            switch info.state
                case 'Unknown'
                    UMITRigLedger.iThrowNoLedger(operation, info.folder);
                case 'NoRig'
                    UMITRigLedger.iThrow('Umitoolbox:RigLedger:noRigBound', ...
                        operation, info.folder, ...
                        'the folder is recorded as "no rig" (imported without a rig binding).', ...
                        {UMITRigLedger.iFixBind(), UMITRigLedger.iFixFreshFolder()});
                case 'Conflicted'
                    lines = cell(1, numel(info.conflictParts));
                    for iPart = 1:numel(info.conflictParts)
                        part = info.conflictParts{iPart};
                        lines{iPart} = sprintf(['the %s used so far came from %s, but the ' ...
                            'folder is now bound to %s, whose %s %s'], ...
                            UMITRigLedger.iPartLabel(part), ...
                            UMITRigLedger.iDescribeHashSource(info, part, info.commitments.(part)), ...
                            UMITRigLedger.iDescribeRig(info.binding.sourceRig), ...
                            UMITRigLedger.iPartLabel(part), UMITRigLedger.iDiffersVerb(part));
                    end
                    UMITRigLedger.iThrow('Umitoolbox:RigLedger:conflicted', ...
                        operation, info.folder, sprintf( ...
                        'the folder rig binding is in conflict with its commitments: %s.', ...
                        strjoin(lines, '; ')), UMITRigLedger.iConflictFixes());
            end
        end

        function iThrow(errorID, operation, folder, why, fixes)
            message = UMITRigLedger.buildErrorMessage(operation, folder, why, fixes);
            throwAsCaller(MException(errorID, '%s', message));
        end

        function iThrowInvalid(operation, folder, issues)
            UMITRigLedger.iThrow('Umitoolbox:RigLedger:invalidLedger', operation, folder, ...
                sprintf('the rig ledger "%s" failed integrity validation: %s', ...
                UMITRigLedger.FileName, strjoin(issues, ' | ')), ...
                {sprintf('Restore %s from a backup copy of this folder.', UMITRigLedger.FileName), ...
                UMITRigLedger.iFixFreshFolder()});
        end

        function iThrowNoLedger(operation, folder)
            UMITRigLedger.iThrow('Umitoolbox:RigLedger:noLedger', operation, folder, ...
                sprintf('the folder has no rig ledger (%s), so its rig is unknown.', ...
                UMITRigLedger.FileName), ...
                {UMITRigLedger.iFixBind(), UMITRigLedger.iFixFreshFolder()});
        end

        function iThrowNotBound(operation, info)
            UMITRigLedger.iThrow('Umitoolbox:RigLedger:notBound', operation, info.folder, ...
                sprintf('the folder state is %s; there is no binding to change.', info.state), ...
                {UMITRigLedger.iFixBind()});
        end

        function iThrowFilterSetModified(operation, info, snapshot, ~)
            % The committed hash (4th argument) is compared by the caller;
            % users are only told the filter set differs.
            UMITRigLedger.iThrow('Umitoolbox:RigLedger:filterSetModified', ...
                operation, info.folder, sprintf( ...
                ['filter set "%s" was changed in the rig library after this folder ' ...
                'used it, so the requested filter set differs from the one used.'], ...
                snapshot.id), ...
                {sprintf(['Restore filter set "%s" in the rig library to its committed ' ...
                'content, or save the modified set under a new filter-set ID and use that ID.'], ...
                snapshot.id), UMITRigLedger.iFixFreshFolder()});
        end

        function fixes = iConflictFixes()
            fixes = { ...
                ['Revert the binding to the committed snapshot: ' ...
                'UMITRigLedger.revertBinding(saveFolder, ''Origin'', ...) ' ...
                '(DataViewer: Rig > Bind/Rebind).'], ...
                UMITRigLedger.iFixFreshFolder(), ...
                ['Acknowledge the conflict: UMITRigLedger.acknowledgeConflict(saveFolder, ' ...
                '''Origin'', ..., ''Comment'', ...) (DataViewer: Rig > Acknowledge Conflict). ' ...
                'Outputs computed before the rebind are not recomputed and the folder ' ...
                'stays marked hadConflict.']};
        end

        function text = iFixBind()
            text = ['Bind a library rig to the folder: UMITRigLedger.bind(saveFolder, ' ...
                'UMITRigLedger.bindingFromRig(rigStore), ''Origin'', ...) ' ...
                '(DataViewer: Rig > Bind/Rebind).'];
        end

        function text = iFixRebind()
            text = 'Rebind the folder to a rig that defines the required part (DataViewer: Rig > Bind/Rebind).';
        end

        function text = iFixUseCommitted(part)
            text = sprintf('Use the committed %s snapshot stored in the folder ledger.', part);
        end

        function text = iFixFreshFolder()
            text = 'Start a fresh save folder, or re-import the raw data into a new save folder.';
        end

        function text = iDescribeRig(sourceRig)
            if isempty(sourceRig)
                text = 'an unknown rig';
                return
            end
            name = sourceRig.displayName;
            if isempty(name)
                name = sourceRig.rigID;
            end
            text = sprintf('rig "%s" (%s)', name, sourceRig.rigID);
        end

        function text = iDescribeHashSource(info, part, hash)
            %IDESCRIBEHASHSOURCE The rig a snapshot came from, for user-facing
            %text. Hashes are never shown to users; only the rigs and the
            %fact that their snapshots differ.
            text = UMITRigLedger.iDescribeRig(UMITRigLedger.iHashSourceRig(info, part, hash));
        end

        function verb = iDiffersVerb(part)
            % "optics" is plural.
            if strcmp(part, 'optics')
                verb = 'differ';
            else
                verb = 'differs';
            end
        end

        function text = iPartLabel(part)
            switch part
                case 'coregistration'
                    text = 'camera coregistration';
                case 'filterSet'
                    text = 'filter set';
                otherwise
                    text = part;
            end
        end

        function sourceRig = iHashSourceRig(info, part, hash)
            sourceRig = [];
            for iSource = numel(info.hashSources):-1:1
                entry = info.hashSources(iSource);
                if strcmp(entry.part, part) && strcmp(entry.hash, hash)
                    sourceRig = entry.sourceRig;
                    return
                end
            end
        end

        function info = iNoteHashSource(info, part, hash, sourceRig)
            info.hashSources(end+1) = struct('part', part, 'hash', hash, 'sourceRig', sourceRig);
        end

        function hash = iCommittedFilterSetHash(info, filterSetID)
            hash = '';
            entries = info.commitments.filterSets;
            for iEntry = 1:numel(entries)
                if strcmp(entries(iEntry).id, filterSetID)
                    hash = entries(iEntry).hash;
                    return
                end
            end
        end

        %% Ledger file I/O
        function [ledger, info] = iLoadForWrite(folder, operation)
            if ~UMITRigLedger.exists(folder)
                ledger = struct( ...
                    'format', UMITRigLedger.Format, ...
                    'schemaVersion', UMITRigLedger.SchemaVersion, ...
                    'ledgerUUID', lower(char(java.util.UUID.randomUUID())), ...
                    'createdOn', UMITRigLedger.iNow(), ...
                    'snapshots', struct('optics', {{}}, 'coregistration', {{}}, 'filterSet', {{}}), ...
                    'events', {{}});
                info = UMITRigLedger.iEmptyInfo(folder);
                return
            end
            try
                ledger = UMITRigLedger.iLoad(folder);
            catch ME
                UMITRigLedger.iThrowInvalid(operation, folder, {ME.message});
            end
            issues = UMITRigLedger.iIntegrityIssues(ledger);
            if isempty(issues)
                try
                    info = UMITRigLedger.iReplay(ledger, folder);
                catch ME
                    issues = {ME.message};
                end
            end
            if ~isempty(issues)
                UMITRigLedger.iThrowInvalid(operation, folder, issues);
            end
        end

        function ledger = iLoad(folder)
            path = fullfile(folder, UMITRigLedger.FileName);
            ledger = UMITRigLedger.iDecodeFile(path);
        end

        function ledger = iDecodeFile(path)
            fid = fopen(path, 'r', 'n', 'UTF-8');
            if fid == -1
                error('Umitoolbox:RigLedger:readFailed', 'Could not read rig ledger: %s', path);
            end
            closer = onCleanup(@() fclose(fid));
            text = fread(fid, [1 Inf], '*char');
            clear closer
            try
                raw = jsondecode(text);
            catch ME
                error('Umitoolbox:RigLedger:readFailed', ...
                    'Rig ledger is not valid JSON (%s): %s', path, ME.message);
            end
            ledger = UMITRigLedger.iNormalizeLedger(raw);
        end

        function ledger = iNormalizeLedger(raw)
            %INORMALIZELEDGER Undo jsondecode shape changes (struct arrays,
            %column vectors, null) so the in-memory model is stable.

            if ~isstruct(raw) || ~isscalar(raw)
                error('Rig ledger top level must be a JSON object.');
            end
            required = {'format', 'schemaVersion', 'ledgerUUID', 'createdOn', 'snapshots', 'events'};
            missing = setdiff(required, fieldnames(raw));
            if ~isempty(missing)
                error('Rig ledger is missing field(s): %s.', strjoin(missing, ', '));
            end
            ledger = raw;
            ledger.events = UMITRigLedger.iCellOfStructs(raw.events, 'events');
            if ~isstruct(raw.snapshots) || ~isscalar(raw.snapshots)
                error('Rig ledger "snapshots" must be an object.');
            end
            snapshots = struct('optics', {{}}, 'coregistration', {{}}, 'filterSet', {{}});
            for iPart = 1:numel(UMITRigLedger.Parts)
                part = UMITRigLedger.Parts{iPart};
                if ~isfield(raw.snapshots, part)
                    error('Rig ledger "snapshots" is missing "%s".', part);
                end
                list = UMITRigLedger.iCellOfStructs(raw.snapshots.(part), ['snapshots.' part]);
                for iSnapshot = 1:numel(list)
                    list{iSnapshot} = UMITRigLedger.iNormalizeSnapshotShape(part, list{iSnapshot});
                end
                snapshots.(part) = list;
            end
            ledger.snapshots = snapshots;
        end

        function list = iCellOfStructs(value, label)
            if isempty(value)
                list = {};
            elseif isstruct(value)
                list = num2cell(value(:)');
            elseif iscell(value) && all(cellfun(@(x) isstruct(x) && isscalar(x), value))
                list = value(:)';
            else
                error('Rig ledger "%s" must be a list of objects.', label);
            end
        end

        function snapshot = iNormalizeSnapshotShape(part, snapshot)
            switch part
                case 'optics'
                    if isfield(snapshot, 'wavelengthNm')
                        snapshot.wavelengthNm = snapshot.wavelengthNm(:)';
                    end
                    for listName = {'cameras', 'illuminations'}
                        if isfield(snapshot, listName{1})
                            list = UMITRigLedger.iCellOfStructs(snapshot.(listName{1}), listName{1});
                            for iItem = 1:numel(list)
                                if isfield(list{iItem}, 'response')
                                    list{iItem}.response = list{iItem}.response(:)';
                                end
                            end
                            if isempty(list)
                                records = UMITRigLedger.iEmptyRecords(listName{1});
                            else
                                records = [list{:}];
                            end
                            snapshot = UMITRigLedger.iSetListField(snapshot, listName{1}, records);
                        end
                    end
                case 'filterSet'
                    for name = {'excitationResponse', 'emissionResponse'}
                        if isfield(snapshot, name{1})
                            snapshot.(name{1}) = snapshot.(name{1})(:)';
                        end
                    end
            end
        end

        function records = iEmptyRecords(listName)
            if strcmp(listName, 'cameras')
                records = struct('index', {}, 'spectrumID', {}, 'model', {}, ...
                    'spectrumSource', {}, 'response', {});
            else
                records = struct('name', {}, 'spectrumID', {}, ...
                    'spectrumSource', {}, 'response', {});
            end
        end

        function issues = iIntegrityIssues(ledger)
            %IINTEGRITYISSUES Schema and snapshot-hash checks (not events).

            issues = {};
            if ~strcmp(UMITRigLedger.iMaybeText(ledger.format), UMITRigLedger.Format)
                issues{end+1} = sprintf('format must be "%s".', UMITRigLedger.Format);
            end
            if ~isnumeric(ledger.schemaVersion) || ~isscalar(ledger.schemaVersion) || ...
                    ledger.schemaVersion ~= UMITRigLedger.SchemaVersion
                issues{end+1} = sprintf('schemaVersion %s is not supported (expected %d).', ...
                    mat2str(ledger.schemaVersion), UMITRigLedger.SchemaVersion);
            end
            if isempty(UMITRigLedger.iMaybeText(ledger.ledgerUUID))
                issues{end+1} = 'ledgerUUID is missing.';
            end
            for iPart = 1:numel(UMITRigLedger.Parts)
                part = UMITRigLedger.Parts{iPart};
                list = ledger.snapshots.(part);
                seen = {};
                for iSnapshot = 1:numel(list)
                    try
                        snapshot = UMITRigLedger.iValidateSnapshot(part, list{iSnapshot});
                        if ismember(snapshot.hash, seen)
                            issues{end+1} = sprintf('Duplicate %s snapshot %s.', part, snapshot.hash); %#ok<AGROW>
                        end
                        seen{end+1} = snapshot.hash; %#ok<AGROW>
                    catch ME
                        issues{end+1} = sprintf('%s snapshot %d: %s', part, iSnapshot, ME.message); %#ok<AGROW>
                    end
                end
            end
            if isempty(ledger.events)
                issues{end+1} = 'The ledger has no events.';
            end
            base = {'eventUUID', 'type', 'timestamp', 'toolboxVersion', 'origin'};
            for iEvent = 1:numel(ledger.events)
                missing = setdiff(base, fieldnames(ledger.events{iEvent}));
                if ~isempty(missing)
                    issues{end+1} = sprintf('Event %d is missing field(s): %s.', ...
                        iEvent, strjoin(missing, ', ')); %#ok<AGROW>
                end
            end
        end

        function iWrite(folder, ledger)
            %IWRITE Atomic replace: temp file, verify by re-reading and
            %replaying, back up the old file, rename, restore on failure.

            encoded = ledger;
            for iPart = 1:numel(UMITRigLedger.Parts)
                part = UMITRigLedger.Parts{iPart};
                list = encoded.snapshots.(part);
                for iSnapshot = 1:numel(list)
                    list{iSnapshot} = UMITRigLedger.iEncodableSnapshot(list{iSnapshot});
                end
                encoded.snapshots.(part) = list;
            end
            text = jsonencode(encoded, 'PrettyPrint', true);

            target = fullfile(folder, UMITRigLedger.FileName);
            tempPath = [tempname(folder) '.json'];
            backupPath = [tempname(folder) '.bak.json'];
            tempCleanup = onCleanup(@() UMITRigLedger.iDeleteIfPresent(tempPath));

            fid = fopen(tempPath, 'w', 'n', 'UTF-8');
            if fid == -1
                error('Umitoolbox:RigLedger:writeFailed', ...
                    'Could not write a temporary rig ledger in "%s".', folder);
            end
            closer = onCleanup(@() fclose(fid));
            fprintf(fid, '%s\n', text);
            clear closer

            try
                written = UMITRigLedger.iDecodeFile(tempPath);
                issues = UMITRigLedger.iIntegrityIssues(written);
                if ~isempty(issues)
                    error('%s', strjoin(issues, ' | '));
                end
                UMITRigLedger.iReplay(written, folder);
            catch ME
                error('Umitoolbox:RigLedger:writeFailed', ...
                    'The new rig ledger for "%s" failed verification and was not installed: %s', ...
                    folder, ME.message);
            end

            backupCreated = false;
            if isfile(target)
                [ok, message] = movefile(target, backupPath, 'f');
                if ~ok
                    error('Umitoolbox:RigLedger:writeFailed', ...
                        'Could not back up the rig ledger "%s": %s', target, message);
                end
                backupCreated = true;
            end
            [ok, message] = movefile(tempPath, target, 'f');
            if ~ok
                if backupCreated && isfile(backupPath)
                    movefile(backupPath, target, 'f');
                end
                error('Umitoolbox:RigLedger:writeFailed', ...
                    'Could not install the rig ledger "%s": %s', target, message);
            end
            if backupCreated
                UMITRigLedger.iDeleteIfPresent(backupPath);
            end
            clear tempCleanup
        end

        function snapshot = iEncodableSnapshot(snapshot)
            % Lists of one record must stay JSON arrays.
            for listName = {'cameras', 'illuminations'}
                if isfield(snapshot, listName{1})
                    snapshot.(listName{1}) = num2cell(snapshot.(listName{1})(:)');
                end
            end
        end

        function ledger = iAppendEvent(ledger, event)
            ledger.events{end+1} = event;
        end

        function ledger = iAddBindingSnapshots(ledger, binding)
            ledger = UMITRigLedger.iAddSnapshot(ledger, 'optics', binding.optics);
            if ~isempty(binding.coregistration)
                ledger = UMITRigLedger.iAddSnapshot(ledger, 'coregistration', binding.coregistration);
            end
        end

        function ledger = iAddSnapshot(ledger, part, snapshot)
            list = ledger.snapshots.(part);
            for iSnapshot = 1:numel(list)
                if strcmp(list{iSnapshot}.hash, snapshot.hash)
                    return
                end
            end
            ledger.snapshots.(part){end+1} = snapshot;
        end

        function event = iNewEvent(type, origin)
            event = struct( ...
                'eventUUID', lower(char(java.util.UUID.randomUUID())), ...
                'type', type, ...
                'timestamp', UMITRigLedger.iNow(), ...
                'toolboxVersion', UMITRigLedger.ToolboxVersion, ...
                'origin', origin);
        end

        %% Locking
        function cleanupObj = iAcquireLock(folder, operation)
            lockFolder = fullfile(folder, UMITRigLedger.LockFolderName);
            lockDirectory = java.io.File(lockFolder);
            if ~lockDirectory.mkdir()
                if ~isfolder(lockFolder)
                    error('Umitoolbox:RigLedger:lockCreationFailed', ...
                        'Could not create the rig-ledger lock: %s', lockFolder);
                end
                age = UMITRigLedger.iLockAge(lockFolder);
                error('Umitoolbox:RigLedger:ledgerLocked', ...
                    ['The rig ledger of "%s" is locked by another operation (lock age %s). ' ...
                    'If no other UMIT process is writing to this folder, clear it with ' ...
                    'UMITRigLedger.clearStaleLock(saveFolder).'], folder, char(age, 'hh:mm:ss'));
            end
            try
                LockInfo = struct('operation', operation, ...
                    'createdOn', datetime('now'), 'processID', feature('getpid'));
                saveMatAtomic(fullfile(lockFolder, 'lockInfo.mat'), 'LockInfo', LockInfo);
            catch ME
                UMITRigLedger.iRemoveFolder(lockFolder);
                rethrow(ME)
            end
            cleanupObj = onCleanup(@() UMITRigLedger.iRemoveFolder(lockFolder));
        end

        function age = iLockAge(lockFolder)
            createdOn = NaT;
            try
                loaded = load(fullfile(lockFolder, 'lockInfo.mat'), 'LockInfo');
                createdOn = loaded.LockInfo.createdOn;
            catch
            end
            if isnat(createdOn)
                listing = dir(lockFolder);
                createdOn = datetime(listing(1).datenum, 'ConvertFrom', 'datenum');
            end
            age = datetime('now') - createdOn;
        end

        function iRemoveFolder(path)
            if isfolder(path)
                rmdir(path, 's');
            end
        end

        function iDeleteIfPresent(path)
            if isfile(path)
                delete(path);
            end
        end

        %% Validation of inputs
        function binding = iValidateBinding(binding)
            if ~isstruct(binding) || ~isscalar(binding) || ...
                    ~all(isfield(binding, {'sourceRig', 'optics', 'coregistration'}))
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'A binding must have fields sourceRig, optics and coregistration.');
            end
            binding.sourceRig = UMITRigLedger.iNormalizeSourceRig(binding.sourceRig, 'Binding');
            binding.optics = UMITRigLedger.iValidateSnapshot('optics', binding.optics);
            if ~isempty(binding.coregistration)
                binding.coregistration = UMITRigLedger.iValidateSnapshot( ...
                    'coregistration', binding.coregistration);
            end
        end

        function sourceRig = iNormalizeSourceRig(sourceRig, where)
            fields = {'uuid', 'rigID', 'displayName', 'libraryModifiedOn'};
            if ~isstruct(sourceRig) || ~isscalar(sourceRig) || ~all(isfield(sourceRig, fields))
                error('Umitoolbox:RigLedger:invalidInput', ...
                    '%s sourceRig must have fields %s.', where, strjoin(fields, ', '));
            end
            for iField = 1:numel(fields)
                sourceRig.(fields{iField}) = UMITRigLedger.iMaybeText(sourceRig.(fields{iField}));
            end
            sourceRig = orderfields(rmfield(sourceRig, setdiff(fieldnames(sourceRig), fields)), fields);
            if isempty(sourceRig.uuid)
                error('Umitoolbox:RigLedger:invalidInput', '%s sourceRig.uuid is empty.', where);
            end
        end

        function snapshot = iValidateSnapshot(part, snapshot)
            required.optics = {'wavelengthNm', 'cameras', 'illuminations', 'hash'};
            required.coregistration = {'resourceUUID', 'sourceChecksum', ...
                'transformClass', 'transformMatrix', 'hash'};
            required.filterSet = {'id', 'excitationSpectrumID', 'emissionSpectrumID', ...
                'excitationResponse', 'emissionResponse', 'hash'};
            if ~isstruct(snapshot) || ~isscalar(snapshot)
                error('Umitoolbox:RigLedger:invalidSnapshot', ...
                    'A %s snapshot must be a scalar struct.', part);
            end
            missing = setdiff(required.(part), fieldnames(snapshot));
            if ~isempty(missing)
                error('Umitoolbox:RigLedger:invalidSnapshot', ...
                    'The %s snapshot is missing field(s): %s.', part, strjoin(missing, ', '));
            end
            extra = setdiff(fieldnames(snapshot), required.(part));
            if ~isempty(extra)
                error('Umitoolbox:RigLedger:invalidSnapshot', ...
                    'The %s snapshot has unexpected field(s): %s.', part, strjoin(extra, ', '));
            end
            expected = UMITRigLedger.computeHash(snapshot);
            if ~strcmp(UMITRigLedger.iMaybeText(snapshot.hash), expected)
                error('Umitoolbox:RigLedger:hashMismatch', ...
                    'The %s snapshot hash %s does not match its content (%s).', part, ...
                    UMITRigLedger.shortHash(UMITRigLedger.iMaybeText(snapshot.hash)), ...
                    UMITRigLedger.shortHash(expected));
            end
        end

        function snapshot = iRequireSnapshot(info, part, hash, where)
            snapshot = UMITRigLedger.iFindSnapshot(info, part, hash);
            if isempty(snapshot)
                error('%s references a missing %s snapshot %s.', where, part, hash);
            end
        end

        function snapshot = iFindSnapshot(info, part, hash)
            snapshot = [];
            list = info.snapshots.(part);
            for iSnapshot = 1:numel(list)
                if strcmp(list{iSnapshot}.hash, hash)
                    snapshot = list{iSnapshot};
                    return
                end
            end
        end

        function iRequireFields(event, fields, where)
            missing = setdiff(fields, fieldnames(event));
            if ~isempty(missing)
                error('%s is missing field(s): %s.', where, strjoin(missing, ', '));
            end
        end

        function hash = iBindingHash(binding, part)
            hash = '';
            if ~isempty(binding.(part))
                hash = binding.(part).hash;
            end
        end

        function part = iPartName(part)
            part = char(string(part));
            if ~ismember(part, UMITRigLedger.Parts)
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'Part must be optics, coregistration or filterSet.');
            end
        end

        function origin = iParseOrigin(functionName, varargin)
            p = inputParser;
            p.FunctionName = functionName;
            addParameter(p, 'Origin', '', @UMITRigLedger.iIsText);
            parse(p, varargin{:});
            origin = UMITRigLedger.iRequireOrigin(p.Results.Origin);
        end

        function origin = iRequireOrigin(origin)
            origin = char(string(origin));
            if isempty(origin)
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'Every ledger write requires a nonempty ''Origin'' (the calling tool).');
            end
        end

        function operation = iParseOperation(functionName, varargin)
            p = inputParser;
            p.FunctionName = functionName;
            addParameter(p, 'Operation', 'Rig-dependent processing', @UMITRigLedger.iIsText);
            parse(p, varargin{:});
            operation = char(string(p.Results.Operation));
        end

        function folder = iExistingFolder(saveFolder)
            folder = UMITRigLedger.iFolderText(saveFolder);
            if ~isfolder(folder)
                error('Umitoolbox:RigLedger:invalidInput', 'Save folder does not exist: %s', folder);
            end
        end

        function folder = iFolderText(saveFolder)
            if ~UMITRigLedger.iIsText(saveFolder) || strlength(string(saveFolder)) == 0
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'saveFolder must be a nonempty character vector or string scalar.');
            end
            folder = char(string(saveFolder));
            if endsWith(folder, filesep) && numel(folder) > 1
                folder = folder(1:end-1);
            end
        end

        %% Snapshot content helpers
        function RigInfo = iRigInfo(rigSource)
            if isa(rigSource, 'UMITRigStore')
                RigInfo = rigSource.getRigInfo();
            elseif isstruct(rigSource) && isscalar(rigSource)
                RigInfo = rigSource;
            else
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'Expected a UMITRigStore object or a RigInfo struct.');
            end
        end

        function list = iStructList(info, fieldName)
            list = [];
            if isfield(info, fieldName) && isstruct(info.(fieldName))
                list = info.(fieldName);
            end
        end

        function [spectrumID, source, response] = iResolveSpectrum(category, spectrumID)
            source = '';
            response = [];
            if isempty(spectrumID)
                return
            end
            spectrum = UMITRigStore.getSpectrum(category, spectrumID);
            spectrumID = spectrum.id;
            builtInRoot = fullfile(fileparts(which('UMITRigStore')), 'resources');
            if startsWith(lower(spectrum.file), lower(builtInRoot))
                source = 'builtIn';
            else
                source = 'user';
            end
            response = UMITRigLedger.iQuantize(spectrum.response(:)');
        end

        function response = iFilterResponse(spectrumID)
            spectrumID = char(string(spectrumID));
            if isempty(spectrumID) || strcmpi(spectrumID, 'none')
                response = ones(1, 301);
            else
                spectrum = UMITRigStore.getSpectrum('filter', spectrumID);
                response = UMITRigLedger.iQuantize(spectrum.response(:)');
            end
        end

        function [transformClass, matrix] = iTransformMatrix(tform)
            transformClass = class(tform);
            if isnumeric(tform)
                matrix = double(tform);
            elseif isobject(tform) && isprop(tform, 'A')
                matrix = double(tform.A);
            elseif isobject(tform) && isprop(tform, 'T')
                matrix = double(tform.T);
            else
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'Unsupported camera-coregistration transform class "%s".', transformClass);
            end
            if ~isequal(size(matrix), [3 3]) || any(~isfinite(matrix(:)))
                error('Umitoolbox:RigLedger:invalidInput', ...
                    'The camera-coregistration transform must be a finite 3x3 matrix.');
            end
        end

        function value = iQuantize(value)
            value = round(double(value), UMITRigLedger.SignificantDigits, 'significant');
        end

        function s = iSetListField(s, fieldName, records)
            s.(fieldName) = records(:)';
            if isempty(records)
                s.(fieldName) = records;
            end
        end

        %% Canonical JSON
        function text = iCanonical(value, forceList)
            if isstruct(value)
                if isscalar(value) && ~forceList
                    keys = sort(fieldnames(value));
                    parts = cell(1, numel(keys));
                    for iKey = 1:numel(keys)
                        key = keys{iKey};
                        parts{iKey} = [UMITRigLedger.iCanonicalString(key) ':' ...
                            UMITRigLedger.iCanonical(value.(key), ...
                            ismember(key, {'cameras', 'illuminations'}))];
                    end
                    text = ['{' strjoin(parts, ',') '}'];
                else
                    items = arrayfun(@(x) UMITRigLedger.iCanonical(x, false), value(:)', ...
                        'UniformOutput', false);
                    text = ['[' strjoin(items, ',') ']'];
                end
            elseif iscell(value)
                items = cellfun(@(x) UMITRigLedger.iCanonical(x, false), value(:)', ...
                    'UniformOutput', false);
                text = ['[' strjoin(items, ',') ']'];
            elseif ischar(value) || isstring(value)
                if isstring(value) && ~isscalar(value)
                    text = UMITRigLedger.iCanonical(cellstr(value), false);
                elseif isstring(value) && ismissing(value)
                    text = 'null';
                else
                    text = UMITRigLedger.iCanonicalString(char(value));
                end
            elseif islogical(value) && isscalar(value)
                if value
                    text = 'true';
                else
                    text = 'false';
                end
            elseif isnumeric(value) || islogical(value)
                value = double(value);
                if isempty(value)
                    text = '[]';
                elseif isscalar(value) && ~forceList
                    text = UMITRigLedger.iCanonicalNumber(value);
                elseif isvector(value)
                    items = arrayfun(@UMITRigLedger.iCanonicalNumber, value(:)', ...
                        'UniformOutput', false);
                    text = ['[' strjoin(items, ',') ']'];
                elseif ismatrix(value)
                    rows = cell(1, size(value, 1));
                    for iRow = 1:size(value, 1)
                        rows{iRow} = UMITRigLedger.iCanonical(value(iRow, :), true);
                    end
                    text = ['[' strjoin(rows, ',') ']'];
                else
                    error('Umitoolbox:RigLedger:hashFailed', ...
                        'N-dimensional arrays cannot be hashed.');
                end
            else
                error('Umitoolbox:RigLedger:hashFailed', ...
                    'Values of class "%s" cannot be hashed.', class(value));
            end
        end

        function text = iCanonicalNumber(value)
            if ~isfinite(value)
                error('Umitoolbox:RigLedger:hashFailed', 'NaN and Inf cannot be hashed.');
            end
            if value == 0
                text = '0';
                return
            end
            formatted = sprintf('%.*e', UMITRigLedger.SignificantDigits - 1, value);
            [mantissa, exponentText] = strtok(formatted, 'e');
            mantissa = regexprep(mantissa, '0+$', '');
            mantissa = regexprep(mantissa, '\.$', '');
            text = sprintf('%se%d', mantissa, str2double(exponentText(2:end)));
        end

        function text = iCanonicalString(value)
            value = char(value);
            value = strrep(value, '\', '\\');
            value = strrep(value, '"', '\"');
            control = value < 32;
            if any(control)
                pieces = cell(1, numel(value));
                for iChar = 1:numel(value)
                    if control(iChar)
                        pieces{iChar} = sprintf('\\u%04x', double(value(iChar)));
                    else
                        pieces{iChar} = value(iChar);
                    end
                end
                value = [pieces{:}];
            end
            text = ['"' value '"'];
        end

        %% Small utilities
        function tf = iIsText(value)
            tf = ischar(value) || (isstring(value) && isscalar(value));
        end

        function text = iText(value, label)
            if ~UMITRigLedger.iIsText(value)
                error('Umitoolbox:RigLedger:invalidInput', '%s must be text.', label);
            end
            text = char(string(value));
        end

        function text = iMaybeText(value)
            if ischar(value)
                text = value;
            elseif isstring(value) && isscalar(value) && ~ismissing(value)
                text = char(value);
            else
                text = '';
            end
        end

        function text = iTextField(s, fieldName)
            text = '';
            if isfield(s, fieldName)
                text = UMITRigLedger.iMaybeText(s.(fieldName));
            end
        end

        function hash = iHashText(value)
            hash = UMITRigLedger.iMaybeText(value);
        end

        function value = iNullable(text)
            if isempty(text)
                value = string(missing);
            else
                value = text;
            end
        end

        function text = iCap(text)
            text = [upper(text(1)) text(2:end)];
        end

        function text = iNow()
            text = char(datetime('now', 'TimeZone', 'UTC', ...
                'Format', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
        end

        function text = iDateText(value)
            if isdatetime(value)
                if isempty(value.TimeZone)
                    value.TimeZone = 'local';
                end
                value.TimeZone = 'UTC';
                value.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''';
                text = char(value);
            else
                text = UMITRigLedger.iMaybeText(value);
            end
        end
    end
end
