classdef UMITLedgerCore
    %UMITLEDGERCORE Shared engine of the folder-bound JSON event ledgers.
    %
    %   spec = UMITLedgerCore.makeSpec(struct(...))
    %   hash = UMITLedgerCore.computeHash(value, hashOptions)
    %   UMITLedgerCore.writeLedger(folder, spec, ledger, verifyFcn)
    %
    %   UMITRigLedger (UMITRigLedger.json) and the project ledger
    %   (UMITProjectLedger.json) are built on this class. It owns everything
    %   that is not specific to one ledger:
    %       - canonical JSON serialization and SHA-256 hashing;
    %       - the ledger header (format, schemaVersion, ledgerUUID,
    %         createdOn, events) and its integrity checks;
    %       - event creation and append;
    %       - atomic temp+verify+rename writes;
    %       - the folder lock (an atomically created lock directory).
    %   Ledger-specific content (snapshots, event semantics, replay, state
    %   derivation) stays in the ledger classes. This class has no state and
    %   never reads or writes anything a caller did not ask for.
    %
    %   Ledger spec (makeSpec): every I/O, lock and input-check function
    %   takes a SPEC struct that names the ledger file, its lock folder, its
    %   format id and schema version, the error-identifier prefix (for
    %   example 'Umitoolbox:RigLedger') and the labels used in messages, so
    %   each ledger keeps its own error identifiers and texts.
    %
    %   Canonical hash form (canonicalJSON/computeHash):
    %       - SHA-256 over the UTF-8 bytes of a canonical JSON text,
    %         reported as 'sha256:<64 lowercase hex>'.
    %       - The top-level 'hash' field of the hashed struct is excluded.
    %       - Objects: keys sorted by ordinal (code-point) order, written as
    %         {"key":value,...} with no whitespace.
    %       - Lists: struct arrays, cell arrays and numeric vectors are
    %         written as [v1,v2,...] whatever their orientation; numeric
    %         matrices as a list of rows. Fields named in
    %         hashOptions.listFields are always lists, even with one element.
    %       - Numbers: 0 is written "0"; any other finite value is rounded
    %         to SignificantDigits (6) significant digits and written as
    %         <mantissa>e<exponent>, where the mantissa has no trailing zeros
    %         (no "." when it is an integer) and the exponent is a plain
    %         signed integer, e.g. 5.21e2, -1.5e-3, 1e0. NaN and Inf are
    %         rejected.
    %       - Text: JSON string with ", \ and control characters escaped
    %         (\" \\ \uXXXX); logical scalars: true/false; empty numeric: [].
    %   Callers that store numbers round them first (quantize), so the
    %   canonical form survives the JSON round trip.
    %
    %   Byte hashing (hashBytes) is a plain SHA-256 over raw bytes with the
    %   same 'sha256:' prefix; it is used for content hashes of arrays.
    %
    %   See also UMITRigLedger.

    properties (Constant)
        HashPrefix = 'sha256:'
        SignificantDigits = 6
        ShortHashDigits = 12
        % Mirrors <param.version> in Umit.prj (no runtime version API exists).
        ToolboxVersion = '2.0.0'
        TimestampFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''
        BaseEventFields = {'eventUUID', 'type', 'timestamp', 'toolboxVersion', 'origin'}
    end

    %% Ledger spec
    methods (Static)
        function spec = makeSpec(options)
            %MAKESPEC Validate and complete a ledger spec.
            %
            %   spec = UMITLedgerCore.makeSpec(options)
            %
            %   Required fields of OPTIONS:
            %       fileName       - ledger file name, e.g. 'UMITRigLedger.json'
            %       lockFolderName - lock directory name
            %       format         - format id written in the file
            %       schemaVersion  - supported schema version (numeric)
            %       errorPrefix    - e.g. 'Umitoolbox:RigLedger'
            %       label          - lowercase noun used in messages,
            %                        e.g. 'rig ledger'
            %       clearLockCall  - how users clear a stale lock, e.g.
            %                        'UMITRigLedger.clearStaleLock(saveFolder)'
            %   Optional fields:
            %       bodyFields     - top-level fields between createdOn and
            %                        events (default {})
            %       folderArgName  - argument name in input errors
            %                        (default 'saveFolder')
            %       folderLabel    - folder noun in input errors
            %                        (default 'Save folder')

            required = {'fileName', 'lockFolderName', 'format', 'schemaVersion', ...
                'errorPrefix', 'label', 'clearLockCall'};
            if ~isstruct(options) || ~isscalar(options) || ~all(isfield(options, required))
                error('Umitoolbox:LedgerCore:invalidSpec', ...
                    'A ledger spec needs fields %s.', strjoin(required, ', '));
            end
            spec = options;
            defaults = struct('bodyFields', {{}}, 'folderArgName', 'saveFolder', ...
                'folderLabel', 'Save folder');
            names = fieldnames(defaults);
            for iName = 1:numel(names)
                if ~isfield(spec, names{iName})
                    spec.(names{iName}) = defaults.(names{iName});
                end
            end
            spec.labelCap = [upper(spec.label(1)) spec.label(2:end)];
            spec.lockLabel = [strrep(spec.label, ' ', '-') ' lock'];
            spec.requiredFields = [{'format', 'schemaVersion', 'ledgerUUID', 'createdOn'}, ...
                spec.bodyFields(:)', {'events'}];
        end
    end

    %% Canonical serialization and hashing
    methods (Static)
        function text = canonicalJSON(value, hashOptions)
            %CANONICALJSON Canonical JSON text used for hashing.
            %
            %   text = UMITLedgerCore.canonicalJSON(value)
            %   text = UMITLedgerCore.canonicalJSON(value, hashOptions)
            %
            %   HASHOPTIONS (optional struct): listFields (cellstr of field
            %   names always written as lists) and errorPrefix (for
            %   <prefix>:hashFailed). See the class help for the exact form.
            %   The top-level 'hash' field of a struct is excluded.

            if nargin < 2
                hashOptions = struct();
            end
            hashOptions = UMITLedgerCore.iHashOptions(hashOptions);
            if isstruct(value) && isscalar(value) && isfield(value, 'hash')
                value = rmfield(value, 'hash');
            end
            text = UMITLedgerCore.iCanonical(value, false, hashOptions);
        end

        function hash = computeHash(value, hashOptions)
            %COMPUTEHASH 'sha256:<hex>' of canonicalJSON(VALUE).

            if nargin < 2
                hashOptions = struct();
            end
            hash = UMITLedgerCore.hashText(UMITLedgerCore.canonicalJSON(value, hashOptions));
        end

        function hash = hashText(text)
            %HASHTEXT 'sha256:<hex>' of the UTF-8 bytes of TEXT.

            bytes = unicode2native(char(string(text)), 'UTF-8');
            hash = UMITLedgerCore.hashBytes(bytes);
        end

        function hash = hashBytes(bytes)
            %HASHBYTES 'sha256:<hex>' of a uint8/int8 byte vector.

            if ~(isa(bytes, 'uint8') || isa(bytes, 'int8'))
                error('Umitoolbox:LedgerCore:invalidInput', ...
                    'hashBytes expects a uint8 or int8 byte array.');
            end
            digest = java.security.MessageDigest.getInstance('SHA-256');
            if ~isempty(bytes)
                if isa(bytes, 'uint8')
                    bytes = typecast(bytes(:), 'int8');
                end
                digest.update(bytes(:));
            end
            digestBytes = typecast(digest.digest(), 'uint8');
            hash = [UMITLedgerCore.HashPrefix ...
                lower(reshape(dec2hex(digestBytes, 2).', 1, []))];
        end

        function text = shortHash(hash)
            %SHORTHASH 'sha256:' plus the first 12 hex digits, or '(none)'.

            hash = UMITLedgerCore.maybeText(hash);
            if isempty(hash)
                text = '(none)';
            else
                text = hash(1:min(numel(hash), ...
                    numel(UMITLedgerCore.HashPrefix) + UMITLedgerCore.ShortHashDigits));
            end
        end

        function value = quantize(value)
            %QUANTIZE Round to SignificantDigits significant digits (double).

            value = round(double(value), UMITLedgerCore.SignificantDigits, 'significant');
        end
    end

    %% Events, time and identifiers
    methods (Static)
        function ledger = newLedger(spec, body)
            %NEWLEDGER Empty ledger: header, BODY fields, then no events.
            %
            %   ledger = UMITLedgerCore.newLedger(spec, body)
            %
            %   BODY is a scalar struct holding spec.bodyFields (or struct()
            %   when the ledger has none). Field order follows the spec.

            ledger = struct( ...
                'format', spec.format, ...
                'schemaVersion', spec.schemaVersion, ...
                'ledgerUUID', UMITLedgerCore.newUUID(), ...
                'createdOn', UMITLedgerCore.nowText());
            for iField = 1:numel(spec.bodyFields)
                ledger.(spec.bodyFields{iField}) = body.(spec.bodyFields{iField});
            end
            ledger.events = {};
        end

        function event = newEvent(type, origin)
            %NEWEVENT Base event: eventUUID, type, timestamp (ISO-8601 UTC),
            %toolboxVersion and origin. No user or host is recorded.

            event = struct( ...
                'eventUUID', UMITLedgerCore.newUUID(), ...
                'type', type, ...
                'timestamp', UMITLedgerCore.nowText(), ...
                'toolboxVersion', UMITLedgerCore.ToolboxVersion, ...
                'origin', origin);
        end

        function ledger = appendEvent(ledger, event)
            %APPENDEVENT Append EVENT; the event list is append-only.

            ledger.events{end+1} = event;
        end

        function uuid = newUUID()
            %NEWUUID Random lowercase v4 UUID.

            uuid = lower(char(java.util.UUID.randomUUID()));
        end

        function text = nowText()
            %NOWTEXT Current time as ISO-8601 UTC text with milliseconds.

            text = char(datetime('now', 'TimeZone', 'UTC', ...
                'Format', UMITLedgerCore.TimestampFormat));
        end

        function text = dateText(value)
            %DATETEXT A datetime as ISO-8601 UTC text (unzoned = local);
            %text is returned unchanged.

            if isdatetime(value)
                if isempty(value.TimeZone)
                    value.TimeZone = 'local';
                end
                value.TimeZone = 'UTC';
                value.Format = UMITLedgerCore.TimestampFormat;
                text = char(value);
            else
                text = UMITLedgerCore.maybeText(value);
            end
        end

        function value = localTime(timestamp)
            %LOCALTIME An event timestamp as an unzoned local datetime.

            value = datetime(timestamp, 'InputFormat', UMITLedgerCore.TimestampFormat, ...
                'TimeZone', 'UTC');
            value.TimeZone = 'local';
            value.TimeZone = '';
        end
    end

    %% Ledger file I/O
    methods (Static)
        function path = filePath(folder, spec)
            %FILEPATH Absolute path of the ledger file in FOLDER.

            path = fullfile(UMITLedgerCore.folderText(folder, spec), spec.fileName);
        end

        function raw = decodeFile(path, spec)
            %DECODEFILE Read and jsondecode a ledger file (no checks).
            %   Throws <prefix>:readFailed when the file cannot be read or
            %   is not JSON.

            fid = fopen(path, 'r', 'n', 'UTF-8');
            if fid == -1
                error([spec.errorPrefix ':readFailed'], 'Could not read %s: %s', spec.label, path);
            end
            closer = onCleanup(@() fclose(fid));
            text = fread(fid, [1 Inf], '*char');
            clear closer
            try
                raw = jsondecode(text);
            catch ME
                error([spec.errorPrefix ':readFailed'], ...
                    '%s is not valid JSON (%s): %s', spec.labelCap, path, ME.message);
            end
        end

        function ledger = normalizeHeader(raw, spec)
            %NORMALIZEHEADER Top-level shape: an object with every required
            %field; events as a 1xN cell of scalar structs. Throws a plain
            %message otherwise. Body fields are left to the caller.

            if ~isstruct(raw) || ~isscalar(raw)
                error('%s top level must be a JSON object.', spec.labelCap);
            end
            missing = setdiff(spec.requiredFields, fieldnames(raw), 'stable');
            if ~isempty(missing)
                error('%s is missing field(s): %s.', spec.labelCap, strjoin(missing, ', '));
            end
            ledger = raw;
            ledger.events = UMITLedgerCore.cellOfStructs(raw.events, 'events', spec);
        end

        function list = cellOfStructs(value, label, spec)
            %CELLOFSTRUCTS A JSON list of objects as a 1xN cell of structs.

            if isempty(value)
                list = {};
            elseif isstruct(value)
                list = num2cell(value(:)');
            elseif iscell(value) && all(cellfun(@(x) isstruct(x) && isscalar(x), value))
                list = value(:)';
            else
                error('%s "%s" must be a list of objects.', spec.labelCap, label);
            end
        end

        function issues = headerIssues(ledger, spec)
            %HEADERISSUES format, schemaVersion and ledgerUUID checks.

            issues = {};
            if ~strcmp(UMITLedgerCore.maybeText(ledger.format), spec.format)
                issues{end+1} = sprintf('format must be "%s".', spec.format);
            end
            if ~isnumeric(ledger.schemaVersion) || ~isscalar(ledger.schemaVersion) || ...
                    ledger.schemaVersion ~= spec.schemaVersion
                issues{end+1} = sprintf('schemaVersion %s is not supported (expected %d).', ...
                    mat2str(ledger.schemaVersion), spec.schemaVersion);
            end
            if isempty(UMITLedgerCore.maybeText(ledger.ledgerUUID))
                issues{end+1} = 'ledgerUUID is missing.';
            end
        end

        function issues = eventIssues(ledger)
            %EVENTISSUES Non-empty event list; base fields; unique UUIDs.

            issues = {};
            if isempty(ledger.events)
                issues{end+1} = 'The ledger has no events.';
            end
            seen = {};
            for iEvent = 1:numel(ledger.events)
                event = ledger.events{iEvent};
                missing = setdiff(UMITLedgerCore.BaseEventFields, fieldnames(event), 'stable');
                if ~isempty(missing)
                    issues{end+1} = sprintf('Event %d is missing field(s): %s.', ...
                        iEvent, strjoin(missing, ', ')); %#ok<AGROW>
                    continue
                end
                uuid = UMITLedgerCore.maybeText(event.eventUUID);
                if ismember(uuid, seen)
                    issues{end+1} = sprintf('Event %d (%s) reuses eventUUID %s.', ...
                        iEvent, UMITLedgerCore.maybeText(event.type), uuid); %#ok<AGROW>
                end
                seen{end+1} = uuid; %#ok<AGROW>
            end
        end

        function writeLedger(folder, spec, ledger, verifyFcn)
            %WRITELEDGER Atomic replace of the ledger file in FOLDER.
            %
            %   UMITLedgerCore.writeLedger(folder, spec, ledger, verifyFcn)
            %
            %   LEDGER is jsonencoded (pretty-printed) to a temp file in
            %   FOLDER. VERIFYFCN(tempPath) re-reads the candidate and must
            %   throw if it is not a valid ledger. The current file is moved
            %   to a backup, the candidate renamed into place (the backup is
            %   restored on failure) and the backup deleted. Any failure
            %   throws <prefix>:writeFailed and leaves the previous file.

            text = jsonencode(ledger, 'PrettyPrint', true);
            target = fullfile(folder, spec.fileName);
            tempPath = [tempname(folder) '.json'];
            backupPath = [tempname(folder) '.bak.json'];
            tempCleanup = onCleanup(@() UMITLedgerCore.deleteIfPresent(tempPath));

            fid = fopen(tempPath, 'w', 'n', 'UTF-8');
            if fid == -1
                error([spec.errorPrefix ':writeFailed'], ...
                    'Could not write a temporary %s in "%s".', spec.label, folder);
            end
            closer = onCleanup(@() fclose(fid));
            fprintf(fid, '%s\n', text);
            clear closer

            try
                verifyFcn(tempPath);
            catch ME
                error([spec.errorPrefix ':writeFailed'], ...
                    'The new %s for "%s" failed verification and was not installed: %s', ...
                    spec.label, folder, ME.message);
            end

            backupCreated = false;
            if isfile(target)
                [ok, message] = movefile(target, backupPath, 'f');
                if ~ok
                    error([spec.errorPrefix ':writeFailed'], ...
                        'Could not back up the %s "%s": %s', spec.label, target, message);
                end
                backupCreated = true;
            end
            [ok, message] = movefile(tempPath, target, 'f');
            if ~ok
                if backupCreated && isfile(backupPath)
                    movefile(backupPath, target, 'f');
                end
                error([spec.errorPrefix ':writeFailed'], ...
                    'Could not install the %s "%s": %s', spec.label, target, message);
            end
            if backupCreated
                UMITLedgerCore.deleteIfPresent(backupPath);
            end
            clear tempCleanup
        end
    end

    %% Locking
    methods (Static)
        function cleanupObj = acquireLock(folder, spec, operation)
            %ACQUIRELOCK Create the folder lock atomically; release on cleanup.
            %
            %   lock = UMITLedgerCore.acquireLock(folder, spec, operation)
            %
            %   The lock is the directory spec.lockFolderName, created with
            %   java.io.File.mkdir (atomic). It holds lockInfo.mat
            %   (operation, createdOn, processID; no user or host). A held
            %   lock throws <prefix>:ledgerLocked with the lock age. The lock
            %   is removed when the returned onCleanup object is destroyed.

            lockFolder = fullfile(folder, spec.lockFolderName);
            lockDirectory = java.io.File(lockFolder);
            if ~lockDirectory.mkdir()
                if ~isfolder(lockFolder)
                    error([spec.errorPrefix ':lockCreationFailed'], ...
                        'Could not create the %s: %s', spec.lockLabel, lockFolder);
                end
                age = UMITLedgerCore.lockAge(lockFolder);
                error([spec.errorPrefix ':ledgerLocked'], ...
                    ['The %s of "%s" is locked by another operation (lock age %s). ' ...
                    'If no other UMIT process is writing to this folder, clear it with ' ...
                    '%s.'], spec.label, folder, char(age, 'hh:mm:ss'), spec.clearLockCall);
            end
            try
                LockInfo = struct('operation', operation, ...
                    'createdOn', datetime('now'), 'processID', feature('getpid'));
                saveMatAtomic(fullfile(lockFolder, 'lockInfo.mat'), 'LockInfo', LockInfo);
            catch ME
                UMITLedgerCore.removeFolder(lockFolder);
                rethrow(ME)
            end
            cleanupObj = onCleanup(@() UMITLedgerCore.removeFolder(lockFolder));
        end

        function wasCleared = clearStaleLock(folder, spec, minimumAge)
            %CLEARSTALELOCK Remove a lock older than MINIMUMAGE (duration).
            %   Returns false when no lock exists; a younger lock throws
            %   <prefix>:lockNotStale.

            folder = UMITLedgerCore.existingFolder(folder, spec);
            lockFolder = fullfile(folder, spec.lockFolderName);
            wasCleared = false;
            if ~isfolder(lockFolder)
                return
            end
            age = UMITLedgerCore.lockAge(lockFolder);
            if age < minimumAge
                error([spec.errorPrefix ':lockNotStale'], ...
                    'The %s of "%s" is only %s old; it may belong to a running operation.', ...
                    spec.lockLabel, folder, char(age, 'hh:mm:ss'));
            end
            rmdir(lockFolder, 's');
            wasCleared = true;
        end

        function age = lockAge(lockFolder)
            %LOCKAGE Age of a lock directory (lockInfo.mat, else its date).

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
    end

    %% Messages, input checks and small utilities
    methods (Static)
        function message = buildErrorMessage(operation, folder, why, fixes)
            %BUILDERRORMESSAGE Shared error text: what, where, why, how to fix.
            %
            %   message = UMITLedgerCore.buildErrorMessage(operation, ...
            %       folder, why, fixes)
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
                char(string(operation)), char(string(folder)), ...
                char(string(why)), strjoin(lines, newline));
        end

        function folder = folderText(folder, spec)
            %FOLDERTEXT Folder argument as char without a trailing separator.

            if ~UMITLedgerCore.isText(folder) || strlength(string(folder)) == 0
                error([spec.errorPrefix ':invalidInput'], ...
                    '%s must be a nonempty character vector or string scalar.', spec.folderArgName);
            end
            folder = char(string(folder));
            if endsWith(folder, filesep) && numel(folder) > 1
                folder = folder(1:end-1);
            end
        end

        function folder = existingFolder(folder, spec)
            %EXISTINGFOLDER folderText, and the folder must exist.

            folder = UMITLedgerCore.folderText(folder, spec);
            if ~isfolder(folder)
                error([spec.errorPrefix ':invalidInput'], '%s does not exist: %s', ...
                    spec.folderLabel, folder);
            end
        end

        function origin = requireOrigin(origin, spec)
            %REQUIREORIGIN Every ledger write names its calling tool.

            origin = char(string(origin));
            if isempty(origin)
                error([spec.errorPrefix ':invalidInput'], ...
                    'Every ledger write requires a nonempty ''Origin'' (the calling tool).');
            end
        end

        function tf = isText(value)
            %ISTEXT char row or string scalar.

            tf = ischar(value) || (isstring(value) && isscalar(value));
        end

        function text = maybeText(value)
            %MAYBETEXT char for text values; '' for anything else.

            if ischar(value)
                text = value;
            elseif isstring(value) && isscalar(value) && ~ismissing(value)
                text = char(value);
            else
                text = '';
            end
        end

        function value = nullable(text)
            %NULLABLE Empty text becomes a missing string (JSON null).

            if isempty(text)
                value = string(missing);
            else
                value = text;
            end
        end

        function deleteIfPresent(path)
            %DELETEIFPRESENT Delete a file if it exists.

            if isfile(path)
                delete(path);
            end
        end

        function removeFolder(path)
            %REMOVEFOLDER Remove a folder tree if it exists.

            if isfolder(path)
                rmdir(path, 's');
            end
        end
    end

    %% Private helpers
    methods (Static, Access = private)
        function hashOptions = iHashOptions(hashOptions)
            if ~isfield(hashOptions, 'listFields')
                hashOptions.listFields = {};
            end
            if ~isfield(hashOptions, 'errorPrefix')
                hashOptions.errorPrefix = 'Umitoolbox:LedgerCore';
            end
        end

        function text = iCanonical(value, forceList, hashOptions)
            if isstruct(value)
                if isscalar(value) && ~forceList
                    keys = sort(fieldnames(value));
                    parts = cell(1, numel(keys));
                    for iKey = 1:numel(keys)
                        key = keys{iKey};
                        parts{iKey} = [UMITLedgerCore.iCanonicalString(key) ':' ...
                            UMITLedgerCore.iCanonical(value.(key), ...
                            ismember(key, hashOptions.listFields), hashOptions)];
                    end
                    text = ['{' strjoin(parts, ',') '}'];
                else
                    items = arrayfun(@(x) UMITLedgerCore.iCanonical(x, false, hashOptions), ...
                        value(:)', 'UniformOutput', false);
                    text = ['[' strjoin(items, ',') ']'];
                end
            elseif iscell(value)
                items = cellfun(@(x) UMITLedgerCore.iCanonical(x, false, hashOptions), ...
                    value(:)', 'UniformOutput', false);
                text = ['[' strjoin(items, ',') ']'];
            elseif ischar(value) || isstring(value)
                if isstring(value) && ~isscalar(value)
                    text = UMITLedgerCore.iCanonical(cellstr(value), false, hashOptions);
                elseif isstring(value) && ismissing(value)
                    text = 'null';
                else
                    text = UMITLedgerCore.iCanonicalString(char(value));
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
                    text = UMITLedgerCore.iCanonicalNumber(value, hashOptions);
                elseif isvector(value)
                    items = arrayfun(@(x) UMITLedgerCore.iCanonicalNumber(x, hashOptions), ...
                        value(:)', 'UniformOutput', false);
                    text = ['[' strjoin(items, ',') ']'];
                elseif ismatrix(value)
                    rows = cell(1, size(value, 1));
                    for iRow = 1:size(value, 1)
                        rows{iRow} = UMITLedgerCore.iCanonical(value(iRow, :), true, hashOptions);
                    end
                    text = ['[' strjoin(rows, ',') ']'];
                else
                    error([hashOptions.errorPrefix ':hashFailed'], ...
                        'N-dimensional arrays cannot be hashed.');
                end
            else
                error([hashOptions.errorPrefix ':hashFailed'], ...
                    'Values of class "%s" cannot be hashed.', class(value));
            end
        end

        function text = iCanonicalNumber(value, hashOptions)
            if ~isfinite(value)
                error([hashOptions.errorPrefix ':hashFailed'], 'NaN and Inf cannot be hashed.');
            end
            if value == 0
                text = '0';
                return
            end
            formatted = sprintf('%.*e', UMITLedgerCore.SignificantDigits - 1, value);
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
    end
end
