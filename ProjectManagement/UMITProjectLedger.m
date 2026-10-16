classdef UMITProjectLedger
    %UMITPROJECTLEDGER Folder-bound project ledger (UMITProjectLedger.json).
    %
    %   info   = UMITProjectLedger.read(folder)
    %   report = UMITProjectLedger.validate(folder)
    %   UMITProjectLedger.bind(folder, binding, snapshot, 'Origin', origin)
    %
    %   Replaces UMITProjectBinding.umitlink. One ledger lives in a bound
    %   SaveFolder, plus one in its RawFolder when that is bound and distinct
    %   (folderRole). It is an append-only JSON event ledger built on
    %   UMITLedgerCore (same conventions as UMITRigLedger: format id and
    %   schemaVersion, atomic temp+verify+rename writes, a folder lock).
    %
    %   Authority: the central project store. The project, subject, session
    %   and binding UUIDs (and the folder role) are the ONLY fields UMIT
    %   reads and validates. Each bind and rebind event also carries a
    %   human snapshot (project metadata, subject ID, session ID, time),
    %   frozen at that moment: UMIT never reads, validates or refreshes it,
    %   and any content is tolerated there.
    %
    %   Folder states (derived from the events, never stored):
    %       Unknown - no ledger file
    %       Bound   - the latest binding event binds the folder
    %       Unbound - the ledger exists but its binding was removed (unbind,
    %                 orphan removal); the history is kept
    %   Orphan, mismatch and conflict states are derived by the store
    %   (UMITProjectStore.resolveProjectBinding) from the binding plus the
    %   store, never from the ledger alone.
    %
    %   Events (base: eventUUID, type, timestamp, toolboxVersion, origin):
    %       bind          - binding, snapshot
    %       rebind        - binding, previousBinding, machine, snapshot
    %       unbind        - previousBinding, reason (store removed or moved
    %                       the session's folder binding)
    %       orphanRemoval - previousBinding, machine (the bound project is
    %                       not on this machine)
    %       freshImport   - (none): the raw data were re-imported; the
    %                       binding is unchanged
    %   A binding is {bindingUUID, projectUUID, subjectUUID, sessionUUID,
    %   folderRole}. History is never rewritten. No reference information
    %   (DataParams owns it) and no group ID are stored.
    %
    %   Reading never writes. A folder without a ledger reads as Unknown and
    %   no file is created.
    %
    %   This class depends only on UMITLedgerCore, never on
    %   UMITProjectStore, so import code can record fresh imports.
    %
    %   Error identifiers: Umitoolbox:ProjectLedger:<condition>.
    %
    %   See also UMITLedgerCore, UMITProjectStore, UMITRigLedger.

    properties (Constant)
        FileName = 'UMITProjectLedger.json'
        LockFolderName = 'UMITProjectLedger.lock'
        Format = 'umit-project-ledger'
        SchemaVersion = 1
        FolderRoles = {'processedDataFolder', 'rawDataFolder'}
        BindingFields = {'bindingUUID', 'projectUUID', 'subjectUUID', 'sessionUUID', 'folderRole'}
        EventTypes = {'bind', 'rebind', 'unbind', 'orphanRemoval', 'freshImport'}
    end

    %% Read-only API
    methods (Static)
        function path = filePath(folder)
            %FILEPATH Absolute path of the ledger file of FOLDER.

            path = UMITLedgerCore.filePath(folder, UMITProjectLedger.iSpec());
        end

        function tf = exists(folder)
            %EXISTS True when FOLDER contains a project ledger file.

            tf = isfile(UMITProjectLedger.filePath(folder));
        end

        function info = read(folder)
            %READ Replay the ledger of FOLDER. Never writes.
            %
            %   Returns folder, ledgerFile, exists, ledgerUUID, createdOn,
            %   state ('Unknown', 'Bound', 'Unbound'), binding (struct of
            %   BindingFields plus boundOn and modifiedOn; [] when not
            %   bound), snapshot (the latest bind/rebind snapshot, human
            %   only; [] when none) and events. An invalid ledger throws
            %   Umitoolbox:ProjectLedger:invalidLedger.

            folder = UMITLedgerCore.folderText(folder, UMITProjectLedger.iSpec());
            report = UMITProjectLedger.validate(folder);
            if ~report.isValid
                UMITProjectLedger.iThrowInvalid('Reading the project ledger', folder, report.issues);
            end
            info = report.info;
        end

        function report = validate(folder)
            %VALIDATE Integrity check (format, identity fields, event order).
            %
            %   Never throws for an invalid ledger and never writes. REPORT
            %   has fields folder, exists, isValid, issues and info. A
            %   folder without a ledger is valid (state Unknown). Snapshot
            %   content is never checked.

            spec = UMITProjectLedger.iSpec();
            folder = UMITLedgerCore.folderText(folder, spec);
            report = struct('folder', folder, 'exists', false, 'isValid', true, ...
                'issues', {{}}, 'info', UMITProjectLedger.iEmptyInfo(folder));
            if ~isfile(fullfile(folder, spec.fileName))
                return
            end
            report.exists = true;
            try
                ledger = UMITProjectLedger.iDecode(fullfile(folder, spec.fileName));
                issues = UMITProjectLedger.iIntegrityIssues(ledger);
                if isempty(issues)
                    report.info = UMITProjectLedger.iReplay(ledger, folder);
                end
            catch ME
                issues = {ME.message};
            end
            report.isValid = isempty(issues);
            report.issues = issues;
            if ~report.isValid
                report.info = [];
            end
        end

        function tf = isBound(folder)
            %ISBOUND True when FOLDER has a valid ledger in state Bound.
            %   An invalid ledger counts as not bound here; use read() to
            %   surface the error.

            report = UMITProjectLedger.validate(folder);
            tf = report.isValid && strcmp(report.info.state, 'Bound');
        end

        function snapshot = makeSnapshot(ProjectInfo, subjectID, sessionID)
            %MAKESNAPSHOT Human snapshot for bind/rebind events.
            %
            %   snapshot = UMITProjectLedger.makeSnapshot(ProjectInfo, ...
            %       subjectID, sessionID)
            %
            %   Copies every project metadata field (registries excluded),
            %   plus the subject and session IDs and the time taken, in a
            %   JSON-friendly form (datetimes as ISO-8601 text). Frozen:
            %   never read back, validated or refreshed.

            project = struct();
            if isstruct(ProjectInfo) && isscalar(ProjectInfo)
                names = setdiff(fieldnames(ProjectInfo), {'subjectRegistry'}, 'stable');
                for iName = 1:numel(names)
                    project.(names{iName}) = UMITProjectLedger.iJsonValue(ProjectInfo.(names{iName}));
                end
            end
            snapshot = struct( ...
                'takenOn', UMITLedgerCore.nowText(), ...
                'project', project, ...
                'subjectID', char(string(subjectID)), ...
                'sessionID', char(string(sessionID)));
        end
    end

    %% Writers
    methods (Static)
        function bind(folder, binding, snapshot, varargin)
            %BIND Bind an Unknown or Unbound folder.
            %
            %   UMITProjectLedger.bind(folder, binding, snapshot, 'Origin', origin)
            %
            %   BINDING has exactly the fields bindingUUID, projectUUID,
            %   subjectUUID, sessionUUID and folderRole. A bound folder
            %   throws Umitoolbox:ProjectLedger:alreadyBound (use rebind).

            origin = UMITProjectLedger.iParseOrigin('bind', varargin{:});
            binding = UMITProjectLedger.iCheckBinding(binding, 'binding');
            [folder, ledger, info, lock] = UMITProjectLedger.iBeginWrite(folder, 'bind'); %#ok<ASGLU>
            if strcmp(info.state, 'Bound')
                UMITProjectLedger.iThrow('Umitoolbox:ProjectLedger:alreadyBound', ...
                    'Binding the folder to a project', folder, ...
                    sprintf('the folder is already bound (session %s of project %s).', ...
                    info.binding.sessionUUID, info.binding.projectUUID), ...
                    {'Rebind the folder instead (UMITProjectLedger.rebind).'});
            end
            event = UMITLedgerCore.newEvent('bind', origin);
            event.binding = binding;
            event.snapshot = snapshot;
            UMITProjectLedger.iAppendAndWrite(folder, ledger, event);
        end

        function rebind(folder, binding, snapshot, varargin)
            %REBIND Replace the binding of a bound folder.
            %
            %   UMITProjectLedger.rebind(folder, binding, snapshot, 'Origin', origin)
            %
            %   Records the previous binding and the machine. A folder that
            %   is not bound throws Umitoolbox:ProjectLedger:notBound.

            origin = UMITProjectLedger.iParseOrigin('rebind', varargin{:});
            binding = UMITProjectLedger.iCheckBinding(binding, 'binding');
            [folder, ledger, info, lock] = UMITProjectLedger.iBeginWrite(folder, 'rebind'); %#ok<ASGLU>
            UMITProjectLedger.iRequireBound(info, 'Rebinding the folder');
            event = UMITLedgerCore.newEvent('rebind', origin);
            event.binding = binding;
            event.previousBinding = UMITProjectLedger.iIdentity(info.binding);
            event.machine = UMITProjectLedger.iMachine();
            event.snapshot = snapshot;
            UMITProjectLedger.iAppendAndWrite(folder, ledger, event);
        end

        function unbind(folder, varargin)
            %UNBIND Record that the store removed the folder's binding.
            %
            %   UMITProjectLedger.unbind(folder, 'Reason', reason, 'Origin', origin)
            %
            %   REASON names the store operation (for example
            %   'sessionRemoved', 'relocated', 'rawFolderUnbound'). The
            %   ledger and its history stay; the folder becomes Unbound.

            p = inputParser;
            p.FunctionName = 'unbind';
            addParameter(p, 'Origin', '', @UMITLedgerCore.isText);
            addParameter(p, 'Reason', '', @UMITLedgerCore.isText);
            parse(p, varargin{:});
            origin = UMITLedgerCore.requireOrigin(p.Results.Origin, UMITProjectLedger.iSpec());
            [folder, ledger, info, lock] = UMITProjectLedger.iBeginWrite(folder, 'unbind'); %#ok<ASGLU>
            UMITProjectLedger.iRequireBound(info, 'Unbinding the folder');
            event = UMITLedgerCore.newEvent('unbind', origin);
            event.previousBinding = UMITProjectLedger.iIdentity(info.binding);
            event.reason = char(string(p.Results.Reason));
            UMITProjectLedger.iAppendAndWrite(folder, ledger, event);
        end

        function recordOrphanRemoval(folder, varargin)
            %RECORDORPHANREMOVAL Record removal of a binding to an absent project.
            %
            %   UMITProjectLedger.recordOrphanRemoval(folder, 'Origin', origin)
            %
            %   The caller (the store) decides that the project is absent
            %   on this machine. The folder becomes Unbound.

            origin = UMITProjectLedger.iParseOrigin('recordOrphanRemoval', varargin{:});
            [folder, ledger, info, lock] = UMITProjectLedger.iBeginWrite(folder, 'orphanRemoval'); %#ok<ASGLU>
            UMITProjectLedger.iRequireBound(info, 'Removing an orphan binding');
            event = UMITLedgerCore.newEvent('orphanRemoval', origin);
            event.previousBinding = UMITProjectLedger.iIdentity(info.binding);
            event.machine = UMITProjectLedger.iMachine();
            UMITProjectLedger.iAppendAndWrite(folder, ledger, event);
        end

        function wasWritten = recordFreshImport(folder, varargin)
            %RECORDFRESHIMPORT Record a fresh import into a bound folder.
            %
            %   wasWritten = UMITProjectLedger.recordFreshImport(folder, 'Origin', origin)
            %
            %   The binding is unchanged. A folder without a ledger, or
            %   whose binding was removed, is not bound: nothing is written
            %   (returns false).

            origin = UMITProjectLedger.iParseOrigin('recordFreshImport', varargin{:});
            wasWritten = false;
            if ~UMITProjectLedger.exists(folder)
                return
            end
            [folder, ledger, info, lock] = UMITProjectLedger.iBeginWrite(folder, 'freshImport'); %#ok<ASGLU>
            if ~strcmp(info.state, 'Bound')
                return
            end
            UMITProjectLedger.iAppendAndWrite(folder, ledger, ...
                UMITLedgerCore.newEvent('freshImport', origin));
            wasWritten = true;
        end

        function wasCleared = clearStaleLock(folder, varargin)
            %CLEARSTALELOCK Remove a leftover ledger lock older than MinimumAge.

            p = inputParser;
            p.FunctionName = 'clearStaleLock';
            addParameter(p, 'MinimumAge', minutes(10), @(x) isduration(x) && isscalar(x));
            parse(p, varargin{:});
            wasCleared = UMITLedgerCore.clearStaleLock(folder, UMITProjectLedger.iSpec(), ...
                p.Results.MinimumAge);
        end
    end

    %% Private helpers
    methods (Static, Access = private)
        function spec = iSpec()
            spec = UMITLedgerCore.makeSpec(struct( ...
                'fileName', UMITProjectLedger.FileName, ...
                'lockFolderName', UMITProjectLedger.LockFolderName, ...
                'format', UMITProjectLedger.Format, ...
                'schemaVersion', UMITProjectLedger.SchemaVersion, ...
                'errorPrefix', 'Umitoolbox:ProjectLedger', ...
                'label', 'project ledger', ...
                'clearLockCall', 'UMITProjectLedger.clearStaleLock(folder)', ...
                'folderArgName', 'folder', ...
                'folderLabel', 'Folder'));
        end

        function info = iEmptyInfo(folder)
            info = struct( ...
                'folder', folder, ...
                'ledgerFile', fullfile(folder, UMITProjectLedger.FileName), ...
                'exists', false, ...
                'ledgerUUID', '', ...
                'createdOn', '', ...
                'state', 'Unknown', ...
                'binding', [], ...
                'snapshot', [], ...
                'events', {{}});
        end

        function ledger = iDecode(path)
            spec = UMITProjectLedger.iSpec();
            ledger = UMITLedgerCore.normalizeHeader(UMITLedgerCore.decodeFile(path, spec), spec);
        end

        function issues = iIntegrityIssues(ledger)
            spec = UMITProjectLedger.iSpec();
            issues = [UMITLedgerCore.headerIssues(ledger, spec), UMITLedgerCore.eventIssues(ledger)];
        end

        function info = iReplay(ledger, folder)
            %IREPLAY Derive the state. Throws a plain message on any
            %violation of the identity schema or the event order.

            info = UMITProjectLedger.iEmptyInfo(folder);
            info.exists = true;
            info.ledgerUUID = ledger.ledgerUUID;
            info.createdOn = ledger.createdOn;
            info.events = ledger.events;
            binding = [];
            for iEvent = 1:numel(ledger.events)
                event = ledger.events{iEvent};
                where = sprintf('Event %d (%s)', iEvent, UMITLedgerCore.maybeText(event.type));
                switch event.type
                    case 'bind'
                        if ~isempty(binding)
                            error('%s appears while the folder is bound (expected rebind).', where);
                        end
                        UMITProjectLedger.iRequireFields(event, {'binding', 'snapshot'}, where);
                        binding = UMITProjectLedger.iCheckBinding(event.binding, [where ' binding']);
                        binding.boundOn = event.timestamp;
                        binding.modifiedOn = event.timestamp;
                        info.snapshot = event.snapshot;
                    case 'rebind'
                        UMITProjectLedger.iRequireFields(event, ...
                            {'binding', 'previousBinding', 'machine', 'snapshot'}, where);
                        UMITProjectLedger.iRequirePrevious(binding, event, where);
                        newBinding = UMITProjectLedger.iCheckBinding(event.binding, [where ' binding']);
                        newBinding.boundOn = binding.boundOn;
                        newBinding.modifiedOn = event.timestamp;
                        binding = newBinding;
                        info.snapshot = event.snapshot;
                    case {'unbind', 'orphanRemoval'}
                        UMITProjectLedger.iRequirePrevious(binding, event, where);
                        binding = [];
                    case 'freshImport'
                        if isempty(binding)
                            error('%s appears while the folder is not bound.', where);
                        end
                    otherwise
                        error('%s has an unknown event type.', where);
                end
            end
            info.binding = binding;
            if isempty(binding)
                info.state = 'Unbound';
            else
                info.state = 'Bound';
            end
        end

        function iRequirePrevious(binding, event, where)
            if isempty(binding)
                error('%s appears while the folder is not bound.', where);
            end
            UMITProjectLedger.iRequireFields(event, {'previousBinding'}, where);
            previous = UMITProjectLedger.iCheckBinding(event.previousBinding, [where ' previousBinding']);
            if ~isequal(previous, UMITProjectLedger.iIdentity(binding))
                error('%s does not match the binding it replaces.', where);
            end
        end

        function binding = iCheckBinding(binding, label)
            %ICHECKBINDING Strict identity check (exact fields, canonical UUIDs).

            fields = UMITProjectLedger.BindingFields;
            if ~isstruct(binding) || ~isscalar(binding)
                error('Umitoolbox:ProjectLedger:invalidBinding', '%s must be a scalar object.', label);
            end
            missing = setdiff(fields, fieldnames(binding));
            extra = setdiff(fieldnames(binding), fields);
            if ~isempty(missing) || ~isempty(extra)
                error('Umitoolbox:ProjectLedger:invalidBinding', ...
                    '%s must have exactly the fields %s.', label, strjoin(fields, ', '));
            end
            uuidPattern = '^[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$';
            for name = {'bindingUUID', 'projectUUID', 'subjectUUID', 'sessionUUID'}
                value = UMITLedgerCore.maybeText(binding.(name{1}));
                if isempty(regexp(value, uuidPattern, 'once'))
                    error('Umitoolbox:ProjectLedger:invalidBinding', ...
                        '%s.%s must be a lowercase canonical UUID.', label, name{1});
                end
                binding.(name{1}) = value;
            end
            role = UMITLedgerCore.maybeText(binding.folderRole);
            if ~ismember(role, UMITProjectLedger.FolderRoles)
                error('Umitoolbox:ProjectLedger:invalidBinding', ...
                    '%s.folderRole must be processedDataFolder or rawDataFolder.', label);
            end
            binding.folderRole = role;
            binding = orderfields(binding, fields);
        end

        function identity = iIdentity(binding)
            identity = struct();
            for name = UMITProjectLedger.BindingFields
                identity.(name{1}) = binding.(name{1});
            end
        end

        function iRequireFields(event, fields, where)
            missing = setdiff(fields, fieldnames(event));
            if ~isempty(missing)
                error('%s is missing field(s): %s.', where, strjoin(missing, ', '));
            end
        end

        function iRequireBound(info, operation)
            if ~strcmp(info.state, 'Bound')
                UMITProjectLedger.iThrow('Umitoolbox:ProjectLedger:notBound', operation, ...
                    info.folder, sprintf('the folder is not bound to a project (state %s).', ...
                    info.state), {'Bind the folder first (UMITProjectLedger.bind).'});
            end
        end

        function [folder, ledger, info, lock] = iBeginWrite(folder, operation)
            spec = UMITProjectLedger.iSpec();
            folder = UMITLedgerCore.existingFolder(folder, spec);
            lock = UMITLedgerCore.acquireLock(folder, spec, operation);
            path = fullfile(folder, spec.fileName);
            if ~isfile(path)
                ledger = UMITLedgerCore.newLedger(spec, struct());
                info = UMITProjectLedger.iEmptyInfo(folder);
                return
            end
            report = UMITProjectLedger.validate(folder);
            if ~report.isValid
                UMITProjectLedger.iThrowInvalid('Writing the project ledger', folder, report.issues);
            end
            ledger = UMITProjectLedger.iDecode(path);
            info = report.info;
        end

        function iAppendAndWrite(folder, ledger, event)
            spec = UMITProjectLedger.iSpec();
            ledger = UMITLedgerCore.appendEvent(ledger, event);
            UMITLedgerCore.writeLedger(folder, spec, ledger, ...
                @(tempPath) UMITProjectLedger.iVerifyCandidate(tempPath, folder));
        end

        function iVerifyCandidate(tempPath, folder)
            ledger = UMITProjectLedger.iDecode(tempPath);
            issues = UMITProjectLedger.iIntegrityIssues(ledger);
            if ~isempty(issues)
                error('%s', strjoin(issues, ' | '));
            end
            UMITProjectLedger.iReplay(ledger, folder);
        end

        function origin = iParseOrigin(functionName, varargin)
            p = inputParser;
            p.FunctionName = functionName;
            addParameter(p, 'Origin', '', @UMITLedgerCore.isText);
            parse(p, varargin{:});
            origin = UMITLedgerCore.requireOrigin(p.Results.Origin, UMITProjectLedger.iSpec());
        end

        function iThrow(errorID, operation, folder, why, fixes)
            message = UMITLedgerCore.buildErrorMessage(operation, folder, why, fixes);
            throwAsCaller(MException(errorID, '%s', message));
        end

        function iThrowInvalid(operation, folder, issues)
            UMITProjectLedger.iThrow('Umitoolbox:ProjectLedger:invalidLedger', operation, folder, ...
                sprintf('the project ledger "%s" failed integrity validation: %s', ...
                UMITProjectLedger.FileName, strjoin(issues, ' | ')), ...
                {sprintf('Restore %s from a backup copy of this folder.', UMITProjectLedger.FileName), ...
                'Ask for the folder to be re-bound to its project session.'});
        end

        function machine = iMachine()
            machine = getenv('COMPUTERNAME');
            if isempty(machine)
                machine = getenv('HOSTNAME');
            end
            if isempty(machine)
                try
                    machine = char(java.net.InetAddress.getLocalHost().getHostName());
                catch
                    machine = '';
                end
            end
        end

        function value = iJsonValue(value)
            %IJSONVALUE JSON-friendly copy of a metadata value (snapshot only).
            if isdatetime(value)
                if isscalar(value) && ~isnat(value)
                    value = UMITLedgerCore.dateText(value);
                elseif isscalar(value)
                    value = '';
                else
                    value = arrayfun(@(d) UMITLedgerCore.dateText(d), value, 'UniformOutput', false);
                end
            elseif isstring(value)
                value = cellstr(value);
            elseif isstruct(value)
                names = fieldnames(value);
                for iItem = 1:numel(value)
                    for iName = 1:numel(names)
                        value(iItem).(names{iName}) = UMITProjectLedger.iJsonValue(value(iItem).(names{iName}));
                    end
                end
            elseif iscell(value)
                value = cellfun(@UMITProjectLedger.iJsonValue, value, 'UniformOutput', false);
            elseif ~(ischar(value) || isnumeric(value) || islogical(value))
                value = char(string(class(value)));
            end
        end
    end
end
