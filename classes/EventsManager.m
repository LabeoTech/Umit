classdef EventsManager < handle
    % EVENTSMANAGER Manage trigger detection and event metadata for umIT.
    %
    % This class detects triggers from LabeoTech analog input recordings,
    % manages event IDs and names, and builds the event information needed
    % to create an "events.mat" file for downstream trial-based analysis.
    %
    % The class also supports a simplified external-signal workflow through
    % GETTRIGGERSFROMSIGNAL, where a synthetic or user-provided trigger
    % vector is processed without relying on acquisition channel metadata.

    properties
        trigThr = 'auto'; % Trigger detection threshold in volts. When set to 'auto', the threshold is placed at 80% of the signal span above the minimum value.
        trigType char = 'EdgeSet' % Trigger interpretation mode.
        trigPolarity char = 'positive' % Trigger polarity. {'positive','negative'}.
        trigChanName = {''}; % Name of AI channel(s) containing the triggers.
        minInterStim single {mustBeNonnegative} = 2 % Minimal inter-stimulus interval in seconds used to merge burst stimuli.
        RawFolder char % Folder containing the ai_xxxx.bin files and optionally the event list file.
    end

    properties (SetAccess = private)
        dictAIChan cell = {'CameraTrig', 'CameraTrig2', 'InputTrigger',...
            'StimAna1','StimAna2', 'StimDig', 'Optogen', 'Unused',...
            'AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}; % List of valid Analog IN channel names.
        AIChanList cell % List of Analog IN channels obtained from AcqInfo.
        SaveFolder char % Folder containing "AcqInfos.mat" and default location for "events.mat".
        AcqInfo struct % Acquisition metadata loaded from AcqInfos.mat or info.txt.
        AnalogIN single % Analog input data array.
        sr single {mustBePositive} = 10000; % Sample rate of the current AnalogIN data in Hz.
        timestamps single {mustBeNonnegative} % Trigger timestamps in seconds.
        state logical {mustBeNonnegative} % Trigger state at each timestamp (1 = ON, 0 = OFF).
        eventID uint16 {mustBePositive} % Condition index for each event transition.
        repetitionID uint16 {mustBePositive} % Repetition index for each event transition.
        eventNameList cell % List of event/condition names.
        baselinePeriod single {mustBePositive} % Trial baseline duration in seconds used for event-based trial splitting.
        selectedEvents logical % Logical event-selection vector with the same size as eventID (true = keep, false = ignore).
        b_isDigital logical = false; % TRUE for digital stimulation when OiS200 is the master trigger source.
        minTrigAmp single {mustBePositive} = .15; % Minimal signal amplitude in volts used when trigThr is 'auto'. Ignored otherwise.
        EventFileName char % Name of the event file containing event information (.csv, .txt, .vpixx, ...).
        b_hasExternalSignal logical = false; % TRUE when current AnalogIN/event data came from getTriggersFromSignal.
    end

    properties (Access = private)
        warnOrigState
        ParseMethods = {'none','csv','vpixx'};
        privateEventFileParseMethod;
        b_LP_applied = false
        triggerTrimInfo struct = struct( ...
            'nLeadDropped', 0, ...
            'nTrailDropped', 0, ...
            'bTrailingRiseKept', false, ...
            'trimApplied', false);
    end

    properties (Dependent)
        EventFileParseMethod char % Name of the method used to read the event file. {'none','csv','vpixx'}.
    end

    methods
        function obj = EventsManager(varargin)
            %EVENTSMANAGER Construct an EventsManager object.
            %
            % This constructor initializes the object from acquisition metadata and/or
            % previously saved event information.
            %
            % Syntax:
            %   obj = EventsManager()
            %   obj = EventsManager(SaveFolder)
            %   obj = EventsManager(SaveFolder, RawFolder)
            %   obj = EventsManager(SaveFolder, RawFolder, ParseMethod)
            %
            % Inputs:
            %   SaveFolder  - Folder containing "AcqInfos.mat" and the default
            %                 location for saving/loading "events.mat".
            %   RawFolder   - Folder containing "ai_*.bin" files and optionally the
            %                 external event file.
            %   ParseMethod - Event file parsing method:
            %                 {'none','csv','vpixx'}.
            %
            % Notes:
            %   - If "events.mat" already exists in SaveFolder, it is loaded first.
            %   - If no "events.mat" is found, acquisition metadata are loaded and the
            %     object proceeds with normal initialization.
            %   - If RawFolder is provided and no saved event file exists, AnalogIN is
            %     read and trigger detection is attempted automatically when possible.

            % Input validation:
            p = inputParser;
            addOptional(p, 'SaveFolder', pwd, @(x) isfolder(convertStringsToChars(x)))
            addOptional(p, 'RawFolder', '', @(x) ischar(x) || isStringScalar(x))
            addOptional(p, 'ParseMethod', 'csv', @(x) ismember(lower(convertStringsToChars(x)), obj.ParseMethods));
            parse(p, varargin{:});

            % Disable backtrace warnings:
            obj.warnOrigState = warning;
            warning('off', 'backtrace');

            % Set main properties:
            obj.SaveFolder = convertStringsToChars(p.Results.SaveFolder);
            obj.RawFolder = convertStringsToChars(p.Results.RawFolder);
            obj.EventFileParseMethod = lower(convertStringsToChars(p.Results.ParseMethod));

            % Load saved events first, when available:
            if isfile(fullfile(obj.SaveFolder, 'events.mat'))
                obj.loadEvents;
                return
            end

            % Otherwise initialize from acquisition metadata:
            try
                obj.setInfo;
            catch
                warning(['The folder "%s" does not contain valid acquisition metadata. ' ...
                    'Set a valid SaveFolder/RawFolder or load an existing "events.mat" file.'], ...
                    obj.SaveFolder);
                return
            end

            % Read AnalogIN if RawFolder was provided:
            if ~isempty(obj.RawFolder)
                obj.setAnalogIN;

                % Detect internally-generated triggers:
                if ~isempty([obj.trigChanName{:}])
                    obj.getTriggers;
                end
            end

            % Force EventFileParseMethod to "none" if the stimulus is digital.
            if obj.b_isDigital && ~strcmpi(obj.EventFileParseMethod, 'none')
                warning('Event file parsing method set to "none" because the stimulus is digital.')
                obj.EventFileParseMethod = 'none';
            end
        end

        function set.EventFileParseMethod(obj, parseMethod)
            %SET.EVENTFILEPARSEMETHOD Validate and set the event-file parser.
            %
            % Valid parsing methods are stored in the private property
            % "ParseMethods".

            parseMethod = lower(convertStringsToChars(parseMethod));
            validateattributes(parseMethod, {'char'}, {'scalartext'}, 'set.EventFileParseMethod');
            msg = [sprintf('Invalid parse method "%s". It must be one of the following:', parseMethod), ...
                sprintf('\n%s', obj.ParseMethods{:})];
            assert(ismember(parseMethod, obj.ParseMethods), msg);
            obj.privateEventFileParseMethod = parseMethod;
        end

        function set.trigChanName(obj, trigName)
            %SET.TRIGCHANNAME Validate and set trigger channel name(s).
            %
            % This setter normalizes the input to a cell array of character vectors and
            % validates each channel name against "AIChanList" using case-insensitive
            % matching.
            %
            % Input:
            %   trigName : char | string scalar | cell array of char/string
            %
            % Notes:
            %   - Empty input is accepted and stored as {''}.
            %   - Valid channel names are canonicalized to the capitalization used in
            %     "AIChanList".

            if ischar(trigName) || isStringScalar(trigName)
                trigName = {convertStringsToChars(trigName)};
            elseif isstring(trigName)
                trigName = cellstr(trigName(:)');
            elseif iscell(trigName)
                assert(all(cellfun(@(x) ischar(x) || isStringScalar(x), trigName)), ...
                    'All trigger channel names must be char vectors or string scalars.');
                trigName = cellfun(@convertStringsToChars, trigName, 'UniformOutput', false);
            else
                error('Format "%s" not supported. Use char, string, or cell array of text.', class(trigName));
            end

            trigName = reshape(trigName, 1, []);

            if isempty(trigName) || (numel(trigName) == 1 && isempty(trigName{1}))
                obj.trigChanName = {''};
                return
            end

            assert(~isempty(obj.AIChanList), ...
                'AIChanList is empty. Acquisition info must be loaded before setting trigChanName.');

            trigNameCanonical = cell(size(trigName));
            for ii = 1:numel(trigName)
                idx = find(strcmpi(trigName{ii}, obj.AIChanList), 1, 'first');
                if isempty(idx)
                    error('Invalid trigger channel "%s". It must match one of the entries in AIChanList.', trigName{ii});
                end
                trigNameCanonical{ii} = obj.AIChanList{idx};
            end

            obj.trigChanName = trigNameCanonical;
        end

        function set.trigType(obj, trigType)
            %SET.TRIGTYPE Validate and set trigger interpretation mode.
            %
            % Valid trigger types are:
            %   'EdgeSet'    : rising edge = ON, falling edge = OFF
            %   'EdgeToggle' : consecutive rising edges alternate ON/OFF
            %
            % Matching is case-insensitive. The stored value is canonicalized for
            % consistency.

            trigType = convertStringsToChars(trigType);
            validateattributes(trigType, {'char'}, {'scalartext'}, 'set.trigType');

            if strcmpi(trigType, 'EdgeSet')
                obj.trigType = 'EdgeSet';
            elseif strcmpi(trigType, 'EdgeToggle')
                obj.trigType = 'EdgeToggle';
            else
                error('Invalid trigger type "%s". It must be either "EdgeSet" or "EdgeToggle".', trigType);
            end
        end

        function set.trigPolarity(obj, trigPolarity)
            %SET.TRIGPOLARITY Validate and set trigger polarity.
            %
            % Valid trigger polarities are:
            %   'positive' : event ON is detected on an upward threshold crossing
            %   'negative' : event ON is detected on a downward threshold crossing
            %
            % Matching is case-insensitive. The stored value is canonicalized for
            % consistency.

            trigPolarity = convertStringsToChars(trigPolarity);
            validateattributes(trigPolarity, {'char'}, {'scalartext'}, 'set.trigPolarity');

            if strcmpi(trigPolarity, 'positive')
                obj.trigPolarity = 'positive';
            elseif strcmpi(trigPolarity, 'negative')
                obj.trigPolarity = 'negative';
            else
                error('Invalid trigger polarity "%s". It must be either "positive" or "negative".', trigPolarity);
            end
        end

        function setBaselinePeriod(obj, baselinePeriod)
            %SETBASELINEPERIOD Set the pre-event baseline duration used for trial splitting.
            %
            % Syntax:
            %   obj.setBaselinePeriod()
            %   obj.setBaselinePeriod(baselinePeriod)
            %
            % Input:
            %   baselinePeriod : scalar baseline duration in seconds
            %
            % Notes:
            %   - If omitted and 2 or more triggers exist, the baseline is set to 20%
            %     of the shortest inter-trigger interval.
            %   - If omitted and only a single trigger exists, the baseline is set to
            %     20% of the remaining analog acquisition duration after trigger onset.
            %   - For a single trigger, long baselines are allowed and any out-of-range
            %     frames are later padded with NaNs by "getFrameMatrix".

            assert(~isempty(obj.state), ...
                'Failed to set trial interval. No triggers detected yet.');

            tm_on = obj.timestamps(obj.state);
            assert(~isempty(tm_on), ...
                'Failed to set baseline period. No trigger onset timestamps available.');

            framePeriod = 1 / obj.AcqInfo.FrameRateHz;

            if nargin < 2
                if numel(tm_on) >= 2
                    trialLength = diff(tm_on);
                    baselinePeriod = 0.2 * min(trialLength);
                else
                    if ~isempty(obj.AnalogIN)
                        acqDuration = size(obj.AnalogIN, 1) / obj.sr;
                        remainingDuration = max(acqDuration - tm_on(1), framePeriod);
                    else
                        remainingDuration = framePeriod;
                    end
                    baselinePeriod = 0.2 * remainingDuration;
                end
            else
                validateattributes(baselinePeriod, {'numeric'}, ...
                    {'scalar', 'real', 'finite', 'positive'}, ...
                    'setBaselinePeriod', 'baselinePeriod');

                if numel(tm_on) >= 2
                    maxBaselineAllowed = min(diff(tm_on) - framePeriod);
                    assert(baselinePeriod < maxBaselineAllowed, ...
                        'Baseline time period is too long. It must be shorter than %0.2f seconds.', ...
                        maxBaselineAllowed);
                end

                assert(baselinePeriod >= framePeriod, ...
                    'Baseline time period is too short. It must be larger than %0.2f seconds.', ...
                    framePeriod);
            end

            obj.baselinePeriod = baselinePeriod;
            fprintf('Baseline period set to %0.2f seconds.\n', obj.baselinePeriod)
        end

        function setAnalogIN(obj)
            %SETANALOGIN Read and validate ai_*.bin files into the "AnalogIN" property.
            %
            % This method:
            %   1) validates the analog input filenames
            %   2) checks for missing zero-indexed files
            %   3) validates file payload length
            %   4) concatenates all files into "AnalogIN"
            %   5) crops the signal to the first/last camera trigger interval
            %
            % File naming convention:
            %   ai_00000.bin, ai_00001.bin, ai_00002.bin, ...
            %
            % Notes:
            %   - File numbering must be zero-indexed and consecutive.
            %   - File payload must be compatible with:
            %         header offset = 5 * 4 bytes
            %         data type      = double
            %         block length   = 1e4 samples
            %         channels       = obj.AcqInfo.AINChannels

            aiFilesList = dir(fullfile(obj.RawFolder, 'ai_*.bin'));
            if isempty(aiFilesList)
                warning('Analog Input files (ai_xxxxx.bin) not found in "%s".', obj.RawFolder)
                return
            end

            if isempty(obj.AcqInfo)
                obj.setInfo;
            end

            expr = '^ai_(\d{5})\.bin$';
            fileNames = {aiFilesList.name};
            tok = regexp(fileNames, expr, 'tokens', 'once');

            assert(all(~cellfun(@isempty, tok)), ...
                ['Invalid analog input file name detected. File names must follow ' ...
                'the convention "ai_0000N.bin" (zero-indexed, 5 digits).']);

            fileIdx = cellfun(@(x) str2double(x{1}), tok);

            assert(numel(unique(fileIdx)) == numel(fileIdx), ...
                'Duplicate analog input file indices detected in RawFolder.');

            [fileIdx, sortIdx] = sort(fileIdx);
            aiFilesList = aiFilesList(sortIdx);

            expectedIdx = 0:(numel(fileIdx)-1);
            if ~isequal(fileIdx(:)', expectedIdx)
                missingIdx = setdiff(expectedIdx, fileIdx);
                if isempty(missingIdx)
                    error('Analog input files must be zero-indexed and consecutive starting at ai_00000.bin.');
                else
                    missingFiles = arrayfun(@(x) sprintf('ai_%05d.bin', x), missingIdx, 'UniformOutput', false);
                    error('Missing analog input file(s): %s', strjoin(missingFiles, ', '));
                end
            end

            disp('Reading analog inputs...')

            obj.AnalogIN = [];
            nChan = obj.AcqInfo.AINChannels;
            headerBytes = 5 * 4;
            bytesPerSample = 8;
            blockLen = 1e4;

            for ind = 1:numel(aiFilesList)
                filePath = fullfile(aiFilesList(ind).folder, aiFilesList(ind).name);

                payloadBytes = aiFilesList(ind).bytes - headerBytes;
                assert(payloadBytes >= 0, ...
                    'File "%s" is smaller than the expected header size.', aiFilesList(ind).name);
                assert(mod(payloadBytes, bytesPerSample) == 0, ...
                    'File "%s" has an invalid payload length for double-precision data.', aiFilesList(ind).name);

                nSamples = payloadBytes / bytesPerSample;
                assert(mod(nSamples, blockLen * nChan) == 0, ...
                    ['File "%s" contains %d data samples after the header. This is not ' ...
                    'compatible with block length %d and %d analog input channels.'], ...
                    aiFilesList(ind).name, nSamples, blockLen, nChan);

                data = memmapfile(filePath, ...
                    'Offset', headerBytes, ...
                    'Format', 'double', ...
                    'Repeat', inf);

                tmp = data.Data;
                tmp = reshape(tmp, blockLen, nChan, []);
                tmp = permute(tmp, [1 3 2]);
                tmp = reshape(tmp, [], nChan);

                obj.AnalogIN = [obj.AnalogIN; tmp]; %#ok<AGROW>
            end

            camT = diff(obj.AnalogIN(:,1) > 2.5);
            camT = [camT; NaN];
            camTOn = find(camT == 1, 1, 'first');
            camTOff = find(camT == -1, 1, 'last');

            if isempty(camTOn) || isempty(camTOff) || camTOff < camTOn
                warning('Failed to locate valid camera trigger bounds. AnalogIN was loaded without cropping.');
            else
                obj.AnalogIN = obj.AnalogIN(camTOn:camTOff, :);
            end

            obj.AnalogIN = single(obj.AnalogIN);
            obj.b_hasExternalSignal = false;
            disp('Done')
        end

        function out = get.EventFileParseMethod(obj)
            %GET.EVENTFILEPARSEMETHOD Return the current event-file parser name.
            out = obj.privateEventFileParseMethod;
        end

        function out = get.repetitionID(obj)
            %GET.REPETITIONID Generate repetition indices from the current event stream.
            %
            % Output:
            %   out - Repetition index vector with the same size as eventID.
            %
            % Notes:
            %   - Repetitions are counted independently for each condition.
            %   - The output is duplicated for ON/OFF transitions.
            %   - This getter preserves the current assumption that condition IDs are
            %     contiguous from 1 to numel(eventNameList).

            out = [];
            if isempty(obj.eventID); return; end

            IDlist = obj.eventID(obj.state);
            out = zeros(size(IDlist), class(IDlist));
            for ii = 1:length(obj.eventNameList)
                idx = IDlist == ii;
                out(idx) = 1:sum(idx);
            end
            out = repelem(out, 2, 1);
        end

        function getTriggers(obj, varargin)
            %GETTRIGGERS Detect triggers from one or more analog input channels.
            %
            % This method detects triggers from the channel(s) stored in
            % "obj.trigChanName" or provided as input. It updates:
            %   - timestamps
            %   - state
            %   - eventID
            %   - eventNameList
            %   - selectedEvents
            %
            % Inputs:
            %   ChannelName (optional): char | string | cellstr
            %       Trigger channel name(s). If omitted, uses obj.trigChanName.
            %   b_verbose (optional, logical, default = true):
            %       If true, displays basic trigger detection statistics.
            %
            % Name-Value Pairs:
            %   'FilterFreq' (scalar, default = 0):
            %       Low-pass cutoff frequency in Hz. Use 0 to disable filtering.
            %
            % Notes:
            %   - Filtering is applied to a working copy of the signal only. The
            %     stored AnalogIN data is preserved.
            %   - If trigThr = 'auto', the threshold is recomputed independently for
            %     each channel.
            %   - If boundary-partial pulses are discarded, the trimming information is
            %     stored in "triggerTrimInfo" and may later be used to align external
            %     event metadata.

            p = inputParser;
            addRequired(p,'obj');
            addOptional(p,'ChannelName','',@(x) ischar(x) || isStringScalar(x) || isstring(x) || iscellstr(x));
            addOptional(p,'b_verbose',true,@(x) islogical(x) && isscalar(x))
            addParameter(p,'FilterFreq',0,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0)
            parse(p,obj,varargin{:})

            chanName = p.Results.ChannelName;
            b_verbose = p.Results.b_verbose;
            FilterFreq = p.Results.FilterFreq;

            if ~isempty(chanName)
                obj.trigChanName = chanName;
            elseif isempty([obj.trigChanName{:}])
                warning('Trigger channel name not set. Trigger detection aborted.')
                return
            end

            if obj.b_hasExternalSignal
                warning(['Current event data came from an external signal. ' ...
                    'Reloading AnalogIN from RawFolder before running metadata-aware trigger detection.']);
                if isempty(obj.RawFolder) || ~isfolder(obj.RawFolder)
                    warning('RawFolder is not available. Trigger detection aborted.')
                    return
                end
                obj.setAnalogIN;
            elseif isempty(obj.AnalogIN)
                obj.setAnalogIN;
            end

            if isempty(obj.AnalogIN)
                warning('AnalogIN is empty. Trigger detection aborted.')
                return
            end

            if FilterFreq > 0
                assert(FilterFreq < obj.sr/2, ...
                    'FilterFreq must be smaller than the Nyquist frequency (%0.3f Hz).', obj.sr/2);
                f = fdesign.lowpass('N,F3dB', 4, FilterFreq, obj.sr);
                lpass = design(f,'butter');
            end

            analogIN_orig = obj.AnalogIN;
            trigThr_orig = obj.trigThr;

            obj.timestamps = [];
            obj.state = [];
            obj.eventID = [];
            obj.eventNameList = {};
            obj.selectedEvents = [];
            obj.triggerTrimInfo = struct( ...
                'nLeadDropped', 0, ...
                'nTrailDropped', 0, ...
                'bTrailingRiseKept', false, ...
                'trimApplied', false);
            obj.b_hasExternalSignal = false;

            nDetected = 0;
            lastTrigThr = trigThr_orig;

            for ii = 1:numel(obj.trigChanName)
                idxCh = strcmpi(obj.trigChanName{ii}, obj.AIChanList);
                if ~any(idxCh)
                    error('Trigger channel "%s" not found in AIChanList.', obj.trigChanName{ii});
                end

                chanData = analogIN_orig(:, idxCh);

                if FilterFreq > 0
                    chanData = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(chanData(:)))');
                else
                    chanData = chanData(:)';
                end

                obj.trigThr = trigThr_orig;

                [tmstmp, chanState, trimInfo] = obj.detectTrig(chanData);

                if isempty(tmstmp)
                    warning('Failed to detect triggers in channel "%s".', obj.trigChanName{ii});
                    continue
                end

                disp(['Triggers detected in channel "' obj.trigChanName{ii} '".']);
                nDetected = nDetected + 1;

                obj.timestamps = [obj.timestamps; tmstmp]; %#ok<AGROW>
                obj.state = [obj.state; logical(chanState)]; %#ok<AGROW>
                obj.eventID = [obj.eventID; repmat(uint16(nDetected), numel(tmstmp), 1)]; %#ok<AGROW>

                if trimInfo.trimApplied || trimInfo.bTrailingRiseKept
                    if nDetected > 1
                        warning(['Boundary-partial trigger handling occurred in more than one trigger channel. ' ...
                            'Automatic alignment of a single external event file may be ambiguous and should be reviewed manually.']);
                    end
                    obj.triggerTrimInfo = trimInfo;
                end

                num = erase(obj.trigChanName{ii}, 'StimAna');
                if startsWith(obj.trigChanName{ii}, 'StimAna', 'IgnoreCase', true) && ...
                        isfield(obj.AcqInfo, ['Stimulation' num '_Name'])
                    obj.eventNameList{nDetected} = obj.AcqInfo.(['Stimulation' num '_Name']);
                else
                    obj.eventNameList{nDetected} = obj.trigChanName{ii};
                end

                lastTrigThr = obj.trigThr;
            end

            obj.AnalogIN = analogIN_orig;

            if isempty(obj.timestamps)
                obj.trigThr = trigThr_orig;
                obj.selectedEvents = [];
                return
            end

            obj.trigThr = lastTrigThr;

            [obj.timestamps, indx] = sort(obj.timestamps);
            obj.state = obj.state(indx);
            obj.eventID = obj.eventID(indx);

            obj.selectedEvents = true(size(obj.eventID));

            if obj.b_isDigital
                eventOrder = obj.AcqInfo.Events_Order(:);

                if obj.triggerTrimInfo.nLeadDropped > 0
                    assert(numel(eventOrder) > obj.triggerTrimInfo.nLeadDropped, ...
                        'Events_Order is too short after trimming leading partial trigger(s).');
                    warning(['Leading partial digital trigger(s) were discarded during detection. ' ...
                        'Trimming the first %d event(s) from AcqInfo.Events_Order to preserve alignment.'], ...
                        obj.triggerTrimInfo.nLeadDropped);
                    eventOrder = eventOrder(obj.triggerTrimInfo.nLeadDropped+1:end);
                end

                if obj.triggerTrimInfo.nTrailDropped > 0
                    assert(numel(eventOrder) > obj.triggerTrimInfo.nTrailDropped, ...
                        'Events_Order is too short after trimming trailing partial trigger(s).');
                    warning(['Trailing partial digital trigger(s) were discarded during detection. ' ...
                        'Trimming the last %d event(s) from AcqInfo.Events_Order to preserve alignment.'], ...
                        obj.triggerTrimInfo.nTrailDropped);
                    eventOrder = eventOrder(1:end-obj.triggerTrimInfo.nTrailDropped);
                end

                assert(numel(eventOrder) == nnz(obj.state), ...
                    ['The number of detected digital trigger onsets does not match ' ...
                    'the number of entries in AcqInfo.Events_Order after boundary trimming.']);

                eventIDtmp = zeros(size(obj.state), 'uint16');
                eventIDtmp(obj.state == 1) = uint16(eventOrder);
                eventIDtmp(obj.state == 0) = uint16(eventOrder);
                obj.eventID = eventIDtmp;

                fn = regexp(fieldnames(obj.AcqInfo),'Stim\d+','match','once');
                fn(cellfun(@isempty,fn)) = [];
                IDs = cellfun(@(x) obj.AcqInfo.(x).ID, fn);
                Names = cellfun(@(x) obj.AcqInfo.(x).Name, fn, 'UniformOutput', false);
                [~,idx] = sort(IDs);
                obj.eventNameList = Names(idx);

                disp('Trigger timestamps generated.');
            end

            obj.selectedEvents = true(size(obj.eventID));

            errID = 'Umitoolbox:EventsManager:IncompatibleArraySizes';
            msg = 'IncompatibleArraySizes: eventID, state, and timestamps must have the same length.';
            assert(isequal(length(obj.eventID), length(obj.state), length(obj.timestamps)), errID, msg)

            if ~isempty(obj.eventNameList)
                msg = 'The unique values of eventID do not match the number of elements in eventNameList.';
                assert(isequal(numel(unique(obj.eventID)), numel(obj.eventNameList)), errID, msg);
            end

            if b_verbose
                disp('Trigger detection completed.')
                disp('---------- Trigger info ----------')
                deltaT = [diff(obj.timestamps); nan];

                if isnumeric(obj.trigThr)
                    thrStr = sprintf('%0.2f', obj.trigThr);
                else
                    thrStr = char(obj.trigThr);
                end

                if strcmpi(obj.trigPolarity, 'negative')
                    trialStateStr = 'LOW';
                    interTrialStateStr = 'HIGH';
                else
                    trialStateStr = 'HIGH';
                    interTrialStateStr = 'LOW';
                end

                fprintf(['Total number of triggers: %d\n' ...
                    'Total number of conditions: %d\n' ...
                    'Average trial length (%s state): %0.3f s\n' ...
                    'Average inter-trial length (%s state): %0.3f s\n' ...
                    'Trigger detection threshold: %s V\n' ...
                    'Trigger type: %s\n' ...
                    'Trigger polarity: %s\n'], ...
                    sum(obj.state), length(obj.eventNameList), ...
                    trialStateStr, mean(deltaT(obj.state == 1), 'omitnan'), ...
                    interTrialStateStr, mean(deltaT(obj.state == 0), 'omitnan'), ...
                    thrStr, obj.trigType, obj.trigPolarity)
                disp('--------------------------------')
            end

            if isempty(obj.baselinePeriod)
                obj.setBaselinePeriod;
            end
        end

        function out = getTriggersFromSignal(obj, signal, sr, varargin)
            %GETTRIGGERSFROMSIGNAL Detect triggers from an external 1-D signal and
            % store the generated event data in the object.
            %
            % This method applies the same optional low-pass filtering and trigger
            % detection logic used internally by the class, but it does not use
            % acquisition metadata to build conditions. The external signal is treated
            % as a single generic condition named "extSignal".
            %
            % Syntax:
            %   out = obj.getTriggersFromSignal(signal, sr)
            %   out = obj.getTriggersFromSignal(signal, sr, b_verbose)
            %   out = obj.getTriggersFromSignal(signal, sr, b_verbose, 'FilterFreq', Fc)
            %
            % Inputs:
            %   signal    - Numeric vector containing the external trigger signal.
            %   sr        - Sampling rate of "signal" in Hz.
            %   b_verbose - (Optional, logical, default = true) If true, display basic
            %               trigger detection statistics.
            %
            % Name-Value Pairs:
            %   'FilterFreq' - Low-pass cutoff frequency in Hz. Use 0 to disable
            %                  filtering. Default = 0.
            %
            % Output:
            %   out - Structure containing the detected trigger/event information:
            %       .timestamps       - Trigger timestamps in seconds.
            %       .state            - Trigger state at each timestamp (true = ON,
            %                           false = OFF).
            %       .eventID          - Event IDs. Always 1 for detected events.
            %       .repetitionID     - Consecutive repetition indices.
            %       .eventNameList    - {'extSignal'}
            %       .selectedEvents   - Logical array with the same size as eventID.
            %       .trigThr          - Trigger threshold used for detection.
            %       .trigType         - Trigger type used for detection.
            %       .sr               - Input sampling rate.
            %       .triggerTrimInfo  - Boundary trimming info returned by detectTrig.
            %       .baselinePeriod   - Baseline period stored in the object after
            %                           trigger detection. Empty when no triggers were
            %                           detected.
            %
            % Notes:
            %   - This method is limited to a single external signal vector.
            %   - Digital-stimulation logic is disabled for this method.
            %   - The generated event data are stored in the object and the object is
            %     switched to external-signal mode.

            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'signal', @(x) isnumeric(x) && isvector(x) && ...
                ~isempty(x) && isreal(x) && all(isfinite(x(:))));
            addRequired(p, 'sr', @(x) isnumeric(x) && isscalar(x) && ...
                isfinite(x) && x > 0);
            addOptional(p, 'b_verbose', true, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FilterFreq', 0, @(x) isnumeric(x) && isscalar(x) && ...
                isfinite(x) && x >= 0);
            parse(p, obj, signal, sr, varargin{:});

            b_verbose = p.Results.b_verbose;
            FilterFreq = p.Results.FilterFreq;

            % Preserve current object state used by detectTrig and for rollback:
            srOrig = obj.sr;
            trigThrOrig = obj.trigThr;
            bIsDigitalOrig = obj.b_isDigital;
            AnalogINOrig = obj.AnalogIN;
            timestampsOrig = obj.timestamps;
            stateOrig = obj.state;
            eventIDOrig = obj.eventID;
            eventNameListOrig = obj.eventNameList;
            selectedEventsOrig = obj.selectedEvents;
            baselinePeriodOrig = obj.baselinePeriod;
            triggerTrimInfoOrig = obj.triggerTrimInfo;
            bHasExternalSignalOrig = obj.b_hasExternalSignal;

            signal = single(signal(:));
            lastTrigThr = trigThrOrig;

            try
                % External signal processing is always generic (non-digital).
                obj.sr = sr;
                obj.b_isDigital = false;

                if FilterFreq > 0
                    assert(FilterFreq < sr/2, ...
                        'FilterFreq must be smaller than the Nyquist frequency (%0.3f Hz).', sr/2);
                    f = fdesign.lowpass('N,F3dB', 4, FilterFreq, sr);
                    lpass = design(f, 'butter');
                    signalProc = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(signal(:)))');
                else
                    signalProc = signal(:)';
                end

                obj.trigThr = trigThrOrig;
                [timestamps, state, trimInfo] = obj.detectTrig(signalProc);
                lastTrigThr = obj.trigThr;

                if isempty(timestamps)
                    eventID = zeros(0,1,'uint16');
                    repetitionID = zeros(0,1,'uint16');
                    selectedEvents = [];
                    baselinePeriod = [];
                else
                    eventID = ones(numel(timestamps), 1, 'uint16');

                    nOn = nnz(state);
                    repetitionID = repelem(uint16((1:nOn)'), 2, 1);
                    repetitionID = repetitionID(1:numel(eventID));

                    selectedEvents = true(size(eventID));
                    baselinePeriod = [];
                end

                out = struct();
                out.timestamps = timestamps;
                out.state = logical(state);
                out.eventID = eventID;
                out.repetitionID = repetitionID;
                out.eventNameList = {'extSignal'};
                out.selectedEvents = selectedEvents;
                out.trigThr = lastTrigThr;
                out.trigType = obj.trigType;
                out.sr = sr;
                out.triggerTrimInfo = trimInfo;
                out.baselinePeriod = baselinePeriod;

                % Store generated event-related info in the object:
                obj.AnalogIN = signal;
                obj.sr = sr;
                obj.timestamps = out.timestamps;
                obj.state = out.state;
                obj.eventID = out.eventID;
                obj.eventNameList = out.eventNameList;
                obj.selectedEvents = out.selectedEvents;
                obj.trigThr = out.trigThr;
                obj.triggerTrimInfo = out.triggerTrimInfo;
                obj.baselinePeriod = out.baselinePeriod;
                obj.b_hasExternalSignal = true;

                % Compute baseline only when triggers exist.
                if ~isempty(obj.timestamps) && isempty(obj.baselinePeriod)
                    obj.setBaselinePeriod;
                    out.baselinePeriod = obj.baselinePeriod;
                end

                if b_verbose
                    disp('Trigger detection completed.')
                    disp('---------- Trigger info ----------')
                    deltaT = [diff(obj.timestamps); nan];

                    if isnumeric(obj.trigThr)
                        thrStr = sprintf('%0.2f', obj.trigThr);
                    else
                        thrStr = char(obj.trigThr);
                    end

                    if strcmpi(obj.trigPolarity, 'negative')
                        trialStateStr = 'LOW';
                        interTrialStateStr = 'HIGH';
                    else
                        trialStateStr = 'HIGH';
                        interTrialStateStr = 'LOW';
                    end

                    fprintf(['Total number of triggers: %d\n' ...
                        'Total number of conditions: %d\n' ...
                        'Average trial length (%s state): %0.3f s\n' ...
                        'Average inter-trial length (%s state): %0.3f s\n' ...
                        'Trigger detection threshold: %s V\n' ...
                        'Trigger type: %s\n' ...
                        'Trigger polarity: %s\n'], ...
                        sum(obj.state), length(obj.eventNameList), ...
                        trialStateStr, mean(deltaT(obj.state == 1), 'omitnan'), ...
                        interTrialStateStr, mean(deltaT(obj.state == 0), 'omitnan'), ...
                        thrStr, obj.trigType, obj.trigPolarity)
                    disp('--------------------------------')
                end


            catch ME
                % Full rollback on failure.
                obj.sr = srOrig;
                obj.trigThr = trigThrOrig;
                obj.b_isDigital = bIsDigitalOrig;
                obj.AnalogIN = AnalogINOrig;
                obj.timestamps = timestampsOrig;
                obj.state = stateOrig;
                obj.eventID = eventIDOrig;
                obj.eventNameList = eventNameListOrig;
                obj.selectedEvents = selectedEventsOrig;
                obj.baselinePeriod = baselinePeriodOrig;
                obj.triggerTrimInfo = triggerTrimInfoOrig;
                obj.b_hasExternalSignal = bHasExternalSignalOrig;
                rethrow(ME)
            end

            % Restore only the temporary digital flag. Keep the external signal
            % data stored in the object.
            obj.b_isDigital = bIsDigitalOrig;
        end

        function f = plot(obj, chanName, varargin)
            %PLOT Plot analog input channels and overlay detected triggers.
            %
            % Supports a special external-signal mode created by
            % GETTRIGGERSFROMSIGNAL. In that case, the current AnalogIN is treated as a
            % single generic signal and channel selection is ignored with a warning.

            if isempty(obj.AnalogIN)
                warning('No signal to plot!')
                f = [];
                return
            end

            if isempty(obj.eventID)
                obj.selectedEvents = [];
            elseif isempty(obj.selectedEvents) || ...
                    ~islogical(obj.selectedEvents) || ...
                    ~isequal(size(obj.selectedEvents), size(obj.eventID))
                warning(['selectedEvents was missing or invalid. ' ...
                    'Resetting it to all TRUE with the same size as eventID.']);
                obj.selectedEvents = true(size(obj.eventID));
            end

            f = [];
            bDownsample = true;
            figTag = 'EventsManager_AnalogINPlot';

            if nargin >= 2 && ~isempty(chanName) && isscalar(chanName) && ishghandle(chanName)
                f = chanName;
                chanName = [];
            elseif nargin < 2
                chanName = [];
            end

            if ~isempty(varargin)
                if isscalar(varargin{1}) && ishghandle(varargin{1})
                    f = varargin{1};
                    if numel(varargin) >= 2
                        bDownsample = varargin{2};
                    end
                else
                    bDownsample = varargin{1};
                end
            end

            validateattributes(bDownsample, {'logical'}, {'scalar'}, 'plot', 'bDownsample');

            if isempty(f)
                oldFig = findall(groot, 'Type', 'figure', 'Tag', figTag);
                if ~isempty(oldFig)
                    delete(oldFig(ishghandle(oldFig)));
                end

                f = figure( ...
                    'Name', 'Analog Inputs', ...
                    'NumberTitle', 'off', ...
                    'Tag', figTag, ...
                    'CreateFcn', {@movegui,'northwest'});
            end

            traceColor = [0 0 0];
            bgColor = [];
            if ishghandle(f)
                figH = ancestor(f, 'figure');
                if isempty(figH) && strcmpi(get(f, 'Type'), 'figure')
                    figH = f;
                end
                if ~isempty(figH) && isprop(figH, 'Color')
                    try
                        bgColor = figH.Color;
                    catch
                    end
                end
                if isempty(bgColor) && isprop(f, 'Color')
                    try
                        bgColor = f.Color;
                    catch
                    end
                end
            end
            if isnumeric(bgColor) && numel(bgColor) == 3 && all(isfinite(bgColor))
                luminance = 0.2126 * bgColor(1) + 0.7152 * bgColor(2) + 0.0722 * bgColor(3);
                if luminance < 0.5
                    traceColor = [1 1 1];
                else
                    traceColor = [0 0 0];
                end
            end

            dsFactor = 1;
            if bDownsample
                dsFactor = 10;
            end

            xVec = (0:size(obj.AnalogIN,1)-1) ./ obj.sr;
            axYSize = [min(obj.AnalogIN(:)), max(obj.AnalogIN(:))];

            if obj.b_hasExternalSignal
                if nargin >= 2 && ~isempty(chanName) && ~(isscalar(chanName) && ishghandle(chanName))
                    warning('Channel selection is ignored when plotting an external signal.')
                end

                ax = subplot(1,1,1, 'Parent', f, 'PlotBoxAspectRatio', [1, .35, 1]);
                ax.YLabel.String = 'amp. (V)';
                ax.XLabel.String = 'time (s)';
                title(ax, 'extSignal', 'Interpreter', 'none');
                hold(ax, 'on');

                ptc = gobjects(0);
                evNames = {};

                if ~isempty(obj.timestamps) && ~isempty(obj.eventID) && any(obj.selectedEvents)
                    idxEv = unique(obj.eventID(obj.selectedEvents), 'stable');
                    colorArr = jet(64);
                    colorArr = colorArr(round(linspace(1, 64, numel(idxEv))), :);
                    evNames = obj.eventNameList(idxEv);

                    for kk = 1:numel(idxEv)
                        xOn = obj.timestamps(obj.eventID == idxEv(kk) & obj.state == 1 & obj.selectedEvents);
                        xOff = obj.timestamps(obj.eventID == idxEv(kk) & obj.state == 0 & obj.selectedEvents);

                        nPatch = min(numel(xOn), numel(xOff));
                        if nPatch == 0
                            continue
                        end

                        xOn = xOn(1:nPatch);
                        xOff = xOff(1:nPatch);

                        x = [xOn xOff xOff xOn];
                        y = repmat([axYSize(1) axYSize(1) axYSize(2) axYSize(2)], nPatch, 1);

                        ptc(kk) = patch( ...
                            ax, ...
                            x', y', colorArr(kk,:), ...
                            'FaceAlpha', 0.25, ...
                            'EdgeColor', 'none', ...
                            'Tag', 'TrigPatch');
                    end
                end

                line( ...
                    xVec(1:dsFactor:end), ...
                    obj.AnalogIN(1:dsFactor:end), ...
                    'LineStyle', '-', ...
                    'Color', traceColor, ...
                    'Parent', ax);

                if isnumeric(obj.trigThr) && isscalar(obj.trigThr) && isfinite(obj.trigThr)
                    ln = line(ax, [xVec(1) xVec(end)], [obj.trigThr obj.trigThr], 'Color', 'r');
                    ln.Tag = 'thrLn';
                end

                if ~isempty(ptc) && ~isempty(evNames)
                    validPatch = isgraphics(ptc);
                    if any(validPatch)
                        legend(ax, ptc(validPatch), evNames(validPatch), ...
                            'Location', 'northeast', ...
                            'Interpreter', 'none');
                    end
                end

                hold(ax, 'off');
                return
            end

            if isempty(chanName)
                chanName = obj.AIChanList;
            end

            if ischar(chanName) || isStringScalar(chanName)
                chanName = {convertStringsToChars(chanName)};
            elseif isstring(chanName)
                chanName = cellstr(chanName(:)');
            elseif iscell(chanName)
                assert(all(cellfun(@(x) ischar(x) || isStringScalar(x), chanName)), ...
                    'Channel names must be text.');
                chanName = cellfun(@convertStringsToChars, chanName, 'UniformOutput', false);
            else
                error('Invalid chanName input. Use char, string, or cell array of text.');
            end

            chanName = reshape(chanName, 1, []);
            chanName = unique(chanName, 'stable');

            b_chanExists = false(size(chanName));
            chanNameCanonical = cell(size(chanName));
            for ii = 1:numel(chanName)
                idx = find(strcmpi(chanName{ii}, obj.AIChanList), 1, 'first');
                if ~isempty(idx)
                    b_chanExists(ii) = true;
                    chanNameCanonical{ii} = obj.AIChanList{idx};
                end
            end

            if all(~b_chanExists)
                error(['Invalid channel name(s). Valid names are:', ...
                    sprintf('\n"%s"', obj.AIChanList{:})]);
            elseif any(~b_chanExists)
                warning(['The following channel(s) do not exist and will be ignored:', ...
                    sprintf('\n"%s"', chanName{~b_chanExists})])
            end

            chanName = chanNameCanonical(b_chanExists);
            [~, chanIndx] = ismember(chanName, obj.AIChanList);

            nRows = ceil(sqrt(numel(chanName)));
            nCols = ceil(numel(chanName) / nRows);

            evNames = {};
            if ~isempty(obj.eventID) && any(obj.selectedEvents)
                evNames = obj.eventNameList(unique(obj.eventID(obj.selectedEvents), 'stable'));
            end

            s = gobjects(1, numel(chanName));

            for ii = 1:numel(chanName)
                s(ii) = subplot(nRows, nCols, ii, 'Parent', f, 'PlotBoxAspectRatio', [1, .35, 1]);
                s(ii).YLabel.String = 'amp. (V)';
                title(s(ii), chanName{ii}, 'Interpreter', 'none');
                hold(s(ii), 'on');

                ptc = gobjects(0);

                if ~isempty(obj.timestamps) && ~isempty(obj.eventID) && any(obj.selectedEvents)
                    idxEv = unique(obj.eventID(obj.selectedEvents), 'stable');

                    colorArr = jet(64);
                    colorArr = colorArr(round(linspace(1, 64, numel(idxEv))), :);

                    for kk = 1:numel(idxEv)
                        xOn = obj.timestamps(obj.eventID == idxEv(kk) & obj.state == 1 & obj.selectedEvents);
                        xOff = obj.timestamps(obj.eventID == idxEv(kk) & obj.state == 0 & obj.selectedEvents);

                        nPatch = min(numel(xOn), numel(xOff));
                        if nPatch == 0
                            continue
                        end

                        xOn = xOn(1:nPatch);
                        xOff = xOff(1:nPatch);

                        x = [xOn xOff xOff xOn];
                        y = repmat([axYSize(1) axYSize(1) axYSize(2) axYSize(2)], nPatch, 1);

                        ptc(kk) = patch( ...
                            s(ii), ...
                            x', y', colorArr(kk,:), ...
                            'FaceAlpha', 0.25, ...
                            'EdgeColor', 'none', ...
                            'Tag', 'TrigPatch');
                    end
                end

                line( ...
                    xVec(1:dsFactor:end), ...
                    obj.AnalogIN(1:dsFactor:end, chanIndx(ii)), ...
                    'LineStyle', '-', ...
                    'Color', traceColor, ...
                    'Parent', s(ii));

                if isnumeric(obj.trigThr) && isscalar(obj.trigThr) && isfinite(obj.trigThr)
                    ln = line(s(ii), [xVec(1) xVec(end)], [obj.trigThr obj.trigThr], 'Color', 'r');
                    ln.Tag = 'thrLn';
                end

                if ii == 1 && ~isempty(ptc) && ~isempty(evNames)
                    validPatch = isgraphics(ptc);
                    if any(validPatch)
                        legend(s(ii), ptc(validPatch), evNames(validPatch), ...
                            'Location', 'northeast', ...
                            'Interpreter', 'none');
                    end
                end

                hold(s(ii), 'off');
            end

            s(end).XLabel.String = 'time (s)';

            if numel(s) > 1
                linkaxes(s, 'xy')
            end
        end

        function status = readConditionFile(obj, varargin)
            %READCONDITIONFILE Read event metadata from an external event file.
            %
            % The file is parsed using the method selected in EventFileParseMethod.
            % If boundary-partial trigger events were discarded during trigger
            % detection, this method trims the parsed event sequence accordingly when
            % the trim is unambiguous.
            %
            % Syntax:
            %   status = obj.readConditionFile()
            %   status = obj.readConditionFile(EventFileName)
            %   status = obj.readConditionFile(..., 'CSVcols', {'Col1','Col2'})
            %
            % Inputs:
            %   EventFileName : char | string scalar, optional
            %       Full path or file name of the event file to read.
            %
            % Name-Value Pairs:
            %   CSVcols : cell array of char
            %       Selected CSV columns to merge into event labels.
            %
            % Output:
            %   status : logical
            %       TRUE if the file was successfully parsed.
            %
            % Notes:
            %   - When an explicit file is provided by the user, parsing errors are not
            %     silently suppressed.
            %   - During automatic file discovery, incompatible files are skipped.

            p = inputParser;
            addRequired(p,'obj')
            addOptional(p,'EventFileName','',@(x) ischar(x) || isStringScalar(x) || isempty(x));
            addParameter(p,'CSVcols',{''},@(x) iscell(x) && ischar([x{:}]));
            parse(p,obj,varargin{:});

            evFile = convertStringsToChars(p.Results.EventFileName);
            colNames = p.Results.CSVcols;
            b_userSpecifiedFile = ~isempty(evFile);

            status = false;

            if strcmpi(obj.EventFileParseMethod, 'none')
                if ~isempty(evFile) || ~isempty(colNames{:})
                    warning('Operation aborted. EventFileParseMethod is not set yet.')
                end
                return
            end

            if isempty(obj.eventID)
                warning(['Operation aborted. Triggers missing. Detect triggers first, then read the ' ...
                    'event file to update the event IDs and names.']);
                return
            end

            if ~isempty(obj.EventFileName) && isempty(evFile)
                evFile = obj.EventFileName;
            end

            if ~isempty(evFile)
                [folder, name, ext] = fileparts(evFile);
                if isempty(folder)
                    folder = obj.RawFolder;
                end

                if isempty(ext)
                    switch lower(obj.EventFileParseMethod)
                        case 'csv'
                            evFile = fullfile(folder, [name '.csv']);

                        case 'vpixx'
                            if exist(fullfile(folder, [name '.vpixx']), 'file')
                                evFile = fullfile(folder, [name '.vpixx']);
                            else
                                evFile = fullfile(folder, [name '.txt']);
                            end

                        otherwise
                            error('Unknown file parsing method.')
                    end
                end

                if ~exist(evFile, 'file')
                    warning('Operation aborted. Events file "%s" does not exist.', evFile);
                    return
                end

                evFile = {evFile};

            else
                fList = [dir(fullfile(obj.RawFolder,'*.csv')); ...
                    dir(fullfile(obj.RawFolder,'*.vpixx')); ...
                    dir(fullfile(obj.RawFolder,'*.txt'))];

                evFile = fullfile({fList.folder}, {fList.name});

                if isempty([evFile{:}])
                    warning('Event list not found. Event IDs and names will not be updated.')
                    return
                end
            end

            evID = [];
            evNames = {};
            evSeqNames = {};

            for ii = 1:numel(evFile)
                obj.EventFileName = evFile{ii};

                try
                    if strcmpi(obj.EventFileParseMethod,'csv') && endsWith(lower(obj.EventFileName), '.csv')
                        [evID, evNames, evSeqNames] = obj.readCSVfile(colNames);

                    elseif strcmpi(obj.EventFileParseMethod,'vpixx') && ...
                            (endsWith(lower(obj.EventFileName), '.vpixx') || endsWith(lower(obj.EventFileName), '.txt'))
                        [evID, evNames, evSeqNames] = obj.readVpixxFile;

                    else
                        continue
                    end

                catch ME
                    if b_userSpecifiedFile
                        rethrow(ME)
                    else
                        continue
                    end
                end

                if ~isempty(evID)
                    disp(['Event list updated from file "' obj.EventFileName '"'])
                    break
                end
            end

            if ~b_userSpecifiedFile && isempty(evID)
                warning('Event file parsing aborted. Failed to parse "%s" file in "%s".', ...
                    upper(obj.EventFileParseMethod), obj.RawFolder)
                return
            end

            assert(~isempty(evID), ...
                'Failed to parse "%s". Is this a valid event file?', obj.EventFileName)

            if obj.triggerTrimInfo.nLeadDropped > 0 || obj.triggerTrimInfo.nTrailDropped > 0
                assert(numel(unique(obj.eventID(obj.state))) <= 1, ...
                    ['Automatic trimming of the external event file is ambiguous because ' ...
                    'triggers were detected from more than one channel. Review manually.']);

                nEvTotal = numel(evSeqNames);
                nLead = obj.triggerTrimInfo.nLeadDropped;
                nTrail = obj.triggerTrimInfo.nTrailDropped;

                assert((nLead + nTrail) < nEvTotal, ...
                    'Boundary trimming would remove all events from "%s".', obj.EventFileName);

                warning(['Boundary-partial trigger event(s) were discarded during trigger detection. ' ...
                    'Trimming %d leading and %d trailing event(s) from "%s" to preserve alignment.'], ...
                    nLead, nTrail, obj.EventFileName);

                evSeqNames = evSeqNames((1+nLead):(end-nTrail));

                evNames = unique(evSeqNames, 'stable');
                evID_on = zeros(numel(evSeqNames), 1, 'uint16');
                for jj = 1:numel(evNames)
                    evID_on(strcmp(evSeqNames, evNames{jj})) = jj;
                end
                evID = repelem(evID_on, 2, 1);
            end

            assert(isequaln(length(obj.timestamps), length(evID)), ...
                ['Failed to update events from file "%s". There is a mismatch between ' ...
                'the number of trigger transitions and the number of event labels after trimming.'], ...
                obj.EventFileName)

            obj.eventID = evID;
            obj.eventNameList = evNames;
            obj.selectedEvents = true(size(obj.eventID));

            status = true;
        end

        function removeCondition(obj, conditionName)
            %REMOVECONDITION Ignore all repetitions of a condition.
            %
            % This method marks all events belonging to the selected condition as
            % ignored by setting the corresponding entries in "selectedEvents" to
            % false.
            %
            % Input:
            %   conditionName : char | string scalar
            %       Name of the condition to ignore. Matching is case-insensitive.
            %
            % Notes:
            %   - Both ON and OFF transitions of the selected condition are ignored.
            %   - If all conditions become ignored, a warning is raised.

            conditionName = convertStringsToChars(conditionName);

            if ~obj.validateCondition(conditionName)
                return
            end

            condID = find(strcmpi(conditionName, obj.eventNameList), 1, 'first');
            obj.selectedEvents(obj.eventID == condID) = false;

            fprintf('Condition "%s" will be ignored!\n', obj.eventNameList{condID})

            if ~any(obj.selectedEvents)
                warning('All conditions were ignored. Data splitting by events will be impossible!')
            end
        end

        function removeRepetition(obj, conditionName, repetitionIndex)
            %REMOVEREPETITION Ignore one or more repetitions from a condition.
            %
            % This method marks selected repetitions from a given condition as ignored
            % by setting the corresponding entries in "selectedEvents" to false.
            %
            % Inputs:
            %   conditionName    : char | string scalar
            %       Name of the condition. Matching is case-insensitive.
            %   repetitionIndex  : positive integer scalar or vector
            %       Repetition index or indices to ignore.
            %
            % Notes:
            %   - Both ON and OFF transitions of each selected repetition are ignored.
            %   - If all repetitions from the selected condition become ignored, a
            %     warning is raised.

            conditionName = convertStringsToChars(conditionName);
            validateattributes(repetitionIndex, {'numeric'}, ...
                {'vector', 'real', 'finite', 'positive', 'integer'}, ...
                'removeRepetition', 'repetitionIndex');

            repetitionIndex = unique(repetitionIndex(:)', 'stable');

            if ~obj.validateCondition(conditionName)
                return
            end

            for ii = 1:numel(repetitionIndex)
                if ~obj.validateRepetition(conditionName, repetitionIndex(ii))
                    return
                end
            end

            condID = find(strcmpi(conditionName, obj.eventNameList), 1, 'first');
            idxCond = (obj.eventID == condID);
            idxRep = ismember(obj.repetitionID, repetitionIndex);

            obj.selectedEvents(idxCond & idxRep) = false;

            if isscalar(repetitionIndex)
                fprintf('Repetition "%d" from condition "%s" will be ignored!\n', ...
                    repetitionIndex, obj.eventNameList{condID});
            else
                fprintf('Repetitions [%s] from condition "%s" will be ignored!\n', ...
                    num2str(repetitionIndex), obj.eventNameList{condID});
            end

            if all(~obj.selectedEvents(idxCond))
                warning('All repetitions from condition "%s" were ignored. The whole condition will be ignored!', ...
                    obj.eventNameList{condID})
            end
        end

        function clearIgnoredEvents(obj)
            %CLEARIGNOREDEVENTS Reset selectedEvents to the default include-all state.
            %
            % Policy:
            %   - If eventID is empty, selectedEvents is set to [].
            %   - Otherwise, selectedEvents is a logical array of TRUE with the same
            %     size as eventID.

            if isempty(obj.eventID)
                obj.selectedEvents = [];
                return
            end

            obj.selectedEvents = true(size(obj.eventID));
            fprintf('Condition and repetition lists successfully reset!\n');
        end

        function saveEvents(obj, saveFolder)
            %SAVEEVENTS Save trigger and event information to "events.mat".
            %
            % Syntax:
            %   obj.saveEvents()
            %   obj.saveEvents(saveFolder)
            %
            % Input:
            %   saveFolder : char | string scalar
            %       Destination folder for "events.mat". If omitted, uses obj.SaveFolder.
            %
            % Notes:
            %   - Saving is aborted if no triggers were detected.
            %   - This method also saves sr, trigThr, and trigPolarity so
            %     external-signal round-trips preserve the signal sampling rate and
            %     trigger detection settings.

            if nargin < 2 || isempty(saveFolder)
                saveFolder = obj.SaveFolder;
            else
                saveFolder = convertStringsToChars(saveFolder);
                validateattributes(saveFolder, {'char'}, {'scalartext'}, 'saveEvents', 'saveFolder');
                assert(isfolder(saveFolder), 'The save folder "%s" does not exist.', saveFolder);
                obj.SaveFolder = saveFolder;
            end

            if isempty(obj.timestamps)
                warning('Unable to create events.mat file. No triggers found!')
                return
            end

            RawFolder = obj.RawFolder; %#ok<NASGU>
            AnalogIN = obj.AnalogIN; %#ok<NASGU>
            sr = obj.sr; %#ok<NASGU>
            trigThr = obj.trigThr; %#ok<NASGU>
            trigPolarity = obj.trigPolarity; %#ok<NASGU>
            eventID = obj.eventID; %#ok<NASGU>
            timestamps = obj.timestamps; %#ok<NASGU>
            state = obj.state; %#ok<NASGU>
            eventNameList = obj.eventNameList; %#ok<NASGU>
            baselinePeriod = obj.baselinePeriod; %#ok<NASGU>
            trigChanName = obj.trigChanName; %#ok<NASGU>
            selectedEvents = obj.selectedEvents; %#ok<NASGU>
            EventFileName = obj.EventFileName; %#ok<NASGU>
            EventFileParseMethod = obj.EventFileParseMethod; %#ok<NASGU>
            b_hasExternalSignal = obj.b_hasExternalSignal; %#ok<NASGU>

            save(fullfile(saveFolder, 'events.mat'), ...
                'RawFolder', ...
                'AnalogIN', ...
                'sr', ...
                'trigThr', ...
                'trigPolarity', ...
                'eventID', ...
                'state', ...
                'timestamps', ...
                'eventNameList', ...
                'baselinePeriod', ...
                'trigChanName', ...
                'selectedEvents', ...
                'EventFileName', ...
                'EventFileParseMethod', ...
                'b_hasExternalSignal');

            disp(['Events MAT file saved in ' saveFolder]);
        end

        function loadEvents(obj, Folder)
            %LOADEVENTS Load variables from events.mat and enforce class invariants.
            %
            % Syntax:
            %   obj.loadEvents()
            %   obj.loadEvents(Folder)
            %
            % Input:
            %   Folder : char | string scalar
            %       Folder containing "events.mat". If omitted, uses obj.SaveFolder
            %       when available, otherwise the current folder.
            %
            % Notes:
            %   - Metadata-aware fields are initialized before loading trigChanName.
            %   - External-signal events can be loaded without acquisition metadata.
            %   - selectedEvents is normalized after loading.

            if nargin < 2 || isempty(Folder)
                if isfile(fullfile(obj.SaveFolder, 'events.mat'))
                    Folder = obj.SaveFolder;
                else
                    Folder = pwd;
                end
            end

            Folder = convertStringsToChars(Folder);
            validateattributes(Folder, {'char'}, {'scalartext'}, 'loadEvents', 'Folder');

            assert(isfile(fullfile(Folder, 'events.mat')), ...
                'Failed to load event info. "events.mat" not found in "%s".', Folder);

            obj.SaveFolder = Folder;
            evInfo = load(fullfile(Folder, 'events.mat'));

            % ---------------------------------------------------------------------
            % Load fields that do not depend on AIChanList / setInfo
            % ---------------------------------------------------------------------
            directFields = { ...
                'RawFolder', ...
                'AnalogIN', ...
                'sr', ...
                'trigThr', ...
                'trigPolarity', ...
                'timestamps', ...
                'state', ...
                'eventID', ...
                'eventNameList', ...
                'baselinePeriod', ...
                'selectedEvents', ...
                'EventFileName', ...
                'EventFileParseMethod'};

            for ii = 1:numel(directFields)
                thisField = directFields{ii};
                if isfield(evInfo, thisField)
                    obj.(thisField) = evInfo.(thisField);
                end
            end

            if isfield(evInfo, 'b_hasExternalSignal')
                obj.b_hasExternalSignal = logical(evInfo.b_hasExternalSignal);
            else
                obj.b_hasExternalSignal = false;
            end

            % ---------------------------------------------------------------------
            % Initialize metadata before loading trigChanName
            % ---------------------------------------------------------------------
            bCanLoadInfo = isfile(fullfile(Folder, 'AcqInfos.mat')) || ...
                (~isempty(obj.RawFolder) && isfolder(obj.RawFolder));

            if bCanLoadInfo
                obj.setInfo;
            elseif obj.b_hasExternalSignal
                obj.AIChanList = {};
            else
                error(['Failed to load acquisition metadata for saved events. ' ...
                    'Neither "AcqInfos.mat" nor a valid RawFolder are available.']);
            end

            % ---------------------------------------------------------------------
            % Load trigChanName only after AIChanList is available
            % ---------------------------------------------------------------------
            if isfield(evInfo, 'trigChanName')
                if ~isempty(obj.AIChanList)
                    obj.trigChanName = evInfo.trigChanName;
                else
                    obj.trigChanName = {''};
                end
            end

            % ---------------------------------------------------------------------
            % Enforce selectedEvents policy
            % ---------------------------------------------------------------------
            if isempty(obj.eventID)
                obj.selectedEvents = [];
            else
                if isempty(obj.selectedEvents) || ...
                        ~islogical(obj.selectedEvents) || ...
                        ~isequal(size(obj.selectedEvents), size(obj.eventID))
                    warning(['selectedEvents is missing or invalid in "events.mat". ' ...
                        'Resetting it to all TRUE with the same size as eventID.']);
                    obj.selectedEvents = true(size(obj.eventID));
                else
                    obj.selectedEvents = logical(obj.selectedEvents);
                end
            end

            % Rebuild baseline only for metadata-aware mode.
            if isempty(obj.baselinePeriod) && ~obj.b_hasExternalSignal && ~isempty(obj.timestamps)
                obj.setBaselinePeriod;
            end

            fprintf('Events info loaded from folder %s\n', Folder);
        end

        function [frMat, conditionIDlist, repetitionList] = getFrameMatrix(obj, datLen, varargin)
            %GETFRAMEMATRIX Generate a repetition-by-frame index matrix for trial splitting.
            %
            % This method outputs frame indices for each selected trial. It supports:
            %   1) all conditions and repetitions
            %   2) one selected condition
            %   3) one selected condition and one or more repetitions
            %
            % Inputs:
            %   datLen          : positive integer scalar
            %       Number of frames in the imaging time series.
            %   conditionName   : optional condition name
            %   repetitionIndex : optional repetition index or indices
            %
            % Outputs:
            %   frMat           : repetition-by-frame matrix of frame indices
            %   conditionIDlist   : condition ID for each returned trial
            %   repetitionList  : repetition index for each returned trial
            %
            % Notes:
            %   - Out-of-bounds frames are set to NaN.
            %   - For a single-trigger acquisition, the trial extends from the computed
            %     start frame to the end of the dataset.
            %   - baselinePeriod must be defined before trial indices are built.

            p = inputParser();
            addRequired(p, 'obj');
            addRequired(p, 'datLen', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x == round(x));
            addOptional(p, 'conditionName', '', @(x) ischar(x) || isStringScalar(x))
            addOptional(p, 'repetitionIndex', [], @(x) (isnumeric(x) && all(x > 0)) || isempty(x))
            parse(p, obj, datLen, varargin{:});

            conditionName = convertStringsToChars(p.Results.conditionName);
            repetitionIndex = p.Results.repetitionIndex;

            if isempty(conditionName) && ~isempty(repetitionIndex)
                error('Repetition index cannot be set without the condition name.')
            end

            if isempty(obj.baselinePeriod)
                obj.setBaselinePeriod;
            end
            assert(~isempty(obj.baselinePeriod), ...
                'baselinePeriod is not set. Detect triggers and set a valid baseline before calling getFrameMatrix.');

            [evIdx, conditionIDlist, repetitionList] = obj.getEventIndex(conditionName, repetitionIndex);
            if isempty(evIdx)
                frMat = [];
                return
            end

            evOnsetFrame = round(obj.baselinePeriod * obj.AcqInfo.FrameRateHz);
            frOn = ceil(obj.timestamps(obj.state) * obj.AcqInfo.FrameRateHz);
            frStart = frOn - evOnsetFrame + 1;

            if numel(frOn) >= 2
                trialLen = [diff(frOn); max(diff(frOn))];
            else
                trialLen = datLen - frStart + 1;
            end

            trialLen = max(trialLen, 1);
            frMat = nan(sum(obj.state), max(trialLen));

            for ii = 1:numel(frStart)
                frVec = frStart(ii):frStart(ii) + trialLen(ii) - 1;
                b_outbound = frVec < 1 | frVec > datLen;
                frVec(b_outbound) = nan;
                frMat(ii, 1:numel(frVec)) = frVec;
            end

            frMat = frMat(evIdx, :);
        end

        function [dataByEv, conditionIDlist, repetitionList] = splitDataByEvents(obj, data, varargin)
            %SPLITDATABYEVENTS Split a 3-D image time series into event-locked trials.
            %
            % Input:
            %   data : numeric 3-D array
            %       Imaging time series with dimensions Y x X x T.
            %
            % Name-Value Pairs:
            %   'condition'  : char | string scalar
            %       Condition name used to filter trials.
            %   'repetition' : numeric vector
            %       Repetition indices to include.
            %
            % Outputs:
            %   dataByEv      : 4-D single array
            %       Trial-split data with dimensions E x Y x X x Ttrial.
            %   conditionIDlist : vector
            %       Condition ID for each returned trial.
            %   repetitionList: vector
            %       Repetition index for each returned trial.
            %
            % Notes:
            %   - Frames outside the valid range are padded with NaNs.
            %   - If trials have unequal lengths, the output is cropped to the shortest
            %     valid trial length.

            p = inputParser();
            addRequired(p, 'obj');
            addRequired(p, 'data', @(x) isnumeric(x) && ndims(x) == 3);
            addParameter(p, 'condition', '', @(x) ischar(x) || isStringScalar(x))
            addParameter(p, 'repetition', [], @(x) (isnumeric(x) && all(x > 0)) || isempty(x))
            parse(p, obj, data, varargin{:});

            [frMat, conditionIDlist, repetitionList] = obj.getFrameMatrix( ...
                size(data,3), p.Results.condition, p.Results.repetition);

            dataByEv = nan(size(frMat,1), size(data,1), size(data,2), size(frMat,2), 'single');

            if isempty(frMat)
                return
            end

            b_validFrames = ~isnan(frMat);
            for ii = 1:size(frMat,1)
                dataByEv(ii,:,:,b_validFrames(ii,:)) = data(:,:,frMat(ii,b_validFrames(ii,:)));
            end

            if any(isnan(frMat(:)))
                disp('Cropping trials to shortest length...');
                firstNaNCol = find(any(isnan(frMat),1), 1, 'first');
                if ~isempty(firstNaNCol)
                    dataByEv(:,:,:,firstNaNCol:end) = [];
                end
            end

            disp('Finished splitting data by events');
        end

        function evInfo = exportEventInfo(obj)
            %EXPORTEVENTINFO Package event-related information into a structure.
            %
            % Output:
            %   evInfo : struct
            %       Contains:
            %           - eventNameList
            %           - baselinePeriod
            %           - FrameRateHz
            %           - eventID        (ON events only)
            %           - selectedEvents (ON events only)

            if isempty(obj.eventID)
                obj.selectedEvents = [];
            elseif isempty(obj.selectedEvents) || ...
                    ~islogical(obj.selectedEvents) || ...
                    ~isequal(size(obj.selectedEvents), size(obj.eventID))
                warning(['selectedEvents was missing or invalid. ' ...
                    'Resetting it to all TRUE with the same size as eventID.']);
                obj.selectedEvents = true(size(obj.eventID));
            end

            fieldsToExport = {'eventNameList','baselinePeriod'};
            evInfo = struct();
            for ii = 1:length(fieldsToExport)
                evInfo.(fieldsToExport{ii}) = obj.(fieldsToExport{ii});
            end

            evInfo.FrameRateHz = obj.AcqInfo.FrameRateHz;
            evInfo.eventID = obj.eventID(obj.state);
            evInfo.selectedEvents = obj.selectedEvents(obj.state);
        end

        function [tmstmp, state] = getConditionTimestamps(obj, conditionName, varargin)
            %GETCONDITIONTIMESTAMPS Return timestamps and states for a condition.
            %
            % Syntax:
            %   [tmstmp, state] = obj.getConditionTimestamps(conditionName)
            %   [tmstmp, state] = obj.getConditionTimestamps(conditionName, repetitionIndex)
            %
            % Inputs:
            %   conditionName   : char | string scalar
            %       Condition name. Matching is case-insensitive.
            %   repetitionIndex : numeric scalar or vector, optional
            %       Repetition index or indices to include.
            %
            % Outputs:
            %   tmstmp : vector
            %       Timestamps for the selected condition/repetitions.
            %   state  : logical vector
            %       Trigger states corresponding to tmstmp.
            %
            % Notes:
            %   - This method returns both ON and OFF transitions.
            %   - It does not assume ON/OFF pairs are adjacent in time, so it is safe
            %     when events from different conditions or channels are interleaved.

            tmstmp = [];
            state = [];

            conditionName = convertStringsToChars(conditionName);
            if ~obj.validateCondition(conditionName)
                return
            end

            repetitionIndex = [];
            if nargin > 2
                repetitionIndex = varargin{1};
            end

            condID = find(strcmpi(conditionName, obj.eventNameList), 1, 'first');
            fullMask = (obj.eventID == condID) & obj.selectedEvents;

            if ~isempty(repetitionIndex)
                validateattributes(repetitionIndex, {'numeric'}, ...
                    {'vector', 'real', 'finite', 'positive', 'integer'}, ...
                    'getConditionTimestamps', 'repetitionIndex');

                repetitionIndex = unique(repetitionIndex(:)', 'stable');

                validReps = unique(obj.repetitionID(obj.state & obj.eventID == condID & obj.selectedEvents), 'stable');
                assert(all(ismember(repetitionIndex, validReps)), ...
                    'One or more requested repetition indices do not exist for condition "%s".', ...
                    obj.eventNameList{condID});

                fullMask = fullMask & ismember(obj.repetitionID, repetitionIndex);
            end

            tmstmp = obj.timestamps(fullMask);
            state = obj.state(fullMask);
        end

        function [evIdx, conditionIDlist, repetitionList, eventNameList] = getEventIndex(obj, conditionName, repetitionIndex)
            %GETEVENTINDEX Retrieve ON-event indices based on condition and repetition.
            %
            %   [evIdx, conditionIDlist, repetitionList, eventNameList] = ...
            %       getEventIndex(obj)
            %
            %   [evIdx, conditionIDlist, repetitionList, eventNameList] = ...
            %       getEventIndex(obj, conditionName)
            %
            %   [evIdx, conditionIDlist, repetitionList, eventNameList] = ...
            %       getEventIndex(obj, conditionName, repetitionIndex)
            %
            %   This helper method returns a logical index over ON events only, together
            %   with the corresponding condition IDs, repetition indices, and event
            %   names for the selected events.
            %
            %   Inputs:
            %       conditionName   : char | string scalar | cellstr
            %           Condition name(s) used to filter events. If empty, all
            %           conditions are used.
            %
            %       repetitionIndex : numeric vector
            %           Repetition indices to keep. If empty, all repetitions are used.
            %
            %   Outputs:
            %       evIdx           : logical vector
            %           Logical index over ON events only.
            %
            %       conditionIDlist : numeric vector
            %           Condition IDs for the selected ON events.
            %
            %       repetitionList  : numeric vector
            %           Repetition indices for the selected ON events.
            %
            %       eventNameList   : cell array of character vectors
            %           Event names for the selected ON events. The output has the same
            %           number of elements and ordering as "conditionIDlist".
            %
            %   Notes:
            %       - Filtering is performed on ON events only.
            %       - Ignored events are excluded automatically.
            %       - Repetition filtering is validated against the currently selected
            %         condition subset.
            %       - eventNameList is generated from obj.eventNameList using the
            %         selected condition IDs.

            if nargin < 2
                conditionName = '';
            end
            if nargin < 3
                repetitionIndex = [];
            end

            evIdx = [];
            conditionIDlist = [];
            repetitionList = [];
            eventNameList = {};

            if isempty(obj.eventID)
                return
            end

            if isempty(obj.selectedEvents) || ...
                    ~islogical(obj.selectedEvents) || ...
                    ~isequal(size(obj.selectedEvents), size(obj.eventID))
                warning(['selectedEvents was missing or invalid. ' ...
                    'Resetting it to all TRUE with the same size as eventID.']);
                obj.selectedEvents = true(size(obj.eventID));
            end

            if isempty(conditionName)
                condID = unique(obj.eventID(obj.state), 'stable')';
            else
                if ischar(conditionName) || isStringScalar(conditionName)
                    conditionName = {convertStringsToChars(conditionName)};
                elseif isstring(conditionName)
                    conditionName = cellstr(conditionName(:)');
                elseif iscell(conditionName)
                    assert(all(cellfun(@(x) ischar(x) || isStringScalar(x), conditionName)), ...
                        'conditionName must contain text values.');
                    conditionName = cellfun(@convertStringsToChars, conditionName, 'UniformOutput', false);
                else
                    error('Invalid conditionName input.');
                end

                conditionName = reshape(conditionName, 1, []);
                conditionName = unique(conditionName, 'stable');

                condID = zeros(1, numel(conditionName), 'like', obj.eventID);
                for ii = 1:numel(conditionName)
                    if ~obj.validateCondition(conditionName{ii})
                        return
                    end
                    condID(ii) = find(strcmpi(conditionName{ii}, obj.eventNameList), 1, 'first');
                end
            end

            condIdxFull = ismember(obj.eventID, condID);
            onMaskBase = obj.state & condIdxFull & obj.selectedEvents;

            if isempty(repetitionIndex)
                evIdxFull = onMaskBase;
            else
                validateattributes(repetitionIndex, {'numeric'}, ...
                    {'vector', 'real', 'finite', 'positive', 'integer'}, ...
                    'getEventIndex', 'repetitionIndex');

                repetitionIndex = unique(repetitionIndex(:)', 'stable');

                validReps = unique(obj.repetitionID(onMaskBase), 'stable');
                assert(all(ismember(repetitionIndex, validReps)), ...
                    'One or more requested repetition indices are not available in the selected condition subset.');

                evIdxFull = onMaskBase & ismember(obj.repetitionID, repetitionIndex);
            end

            conditionIDlist = obj.eventID(evIdxFull);
            repetitionList = obj.repetitionID(evIdxFull);
            evIdx = evIdxFull(obj.state);

            % Return event names with the same size/order as conditionIDlist.
            if isempty(conditionIDlist)
                eventNameList = {};
            else
                eventNameList = obj.eventNameList(conditionIDlist);
            end
        end

    end
    methods (Access = private)
        function setInfo(obj)
            %SETINFO Read acquisition metadata and update dependent class properties.
            %
            % This method loads the acquisition info structure from "AcqInfos.mat" in
            % SaveFolder when available. If that file is not present, it falls back to
            % reading the raw acquisition metadata from RawFolder using ReadInfoFile.
            %
            % It also:
            %   - builds AIChanList
            %   - validates channel names
            %   - updates the analog input sample rate
            %   - infers trigger channels for analog/digital stimulation
            %   - updates minInterStim from stimulation timing metadata
            %
            % Notes:
            %   - AICh fields are sorted numerically (AICh1, AICh2, ..., AICh10) to
            %     avoid lexical ordering issues.
            %   - Older metadata layouts without explicit AICh fields are handled for
            %     retrocompatibility.

            acqInfoFile = fullfile(obj.SaveFolder, 'AcqInfos.mat');

            if isfile(acqInfoFile)
                a = load(acqInfoFile);
                assert(isfield(a, 'AcqInfoStream'), ...
                    '"%s" does not contain the variable "AcqInfoStream".', acqInfoFile);
                obj.AcqInfo = a.AcqInfoStream;
            else
                assert(~isempty(obj.RawFolder) && isfolder(obj.RawFolder), ...
                    ['Failed to load acquisition info. "AcqInfos.mat" was not found in SaveFolder ' ...
                    'and RawFolder is empty or invalid.']);
                obj.AcqInfo = ReadInfoFile(obj.RawFolder);
            end

            fn = fieldnames(obj.AcqInfo);

            idxChan = startsWith(fn, 'AICh', 'IgnoreCase', true);

            if ~any(idxChan)
                obj.AIChanList = cell(obj.AcqInfo.AINChannels, 1);

                switch obj.AcqInfo.AINChannels
                    case 1
                        obj.AIChanList(1) = {'CameraTrig'};
                        idx = 2;

                    case {2,6,10}
                        obj.AIChanList(1:2) = {'CameraTrig', 'StimAna1'};
                        idx = 3;

                    case {3,7,11}
                        obj.AIChanList(1:3) = {'CameraTrig', 'StimAna1', 'StimAna2'};
                        idx = 4;

                    otherwise
                        error('Failed to infer channel names from acquisition metadata.')
                end

                for ii = 0:obj.AcqInfo.AINChannels-idx
                    obj.AIChanList{ii+idx} = ['AI', num2str(ii+1)];
                end

                for ii = 1:numel(obj.AIChanList)
                    obj.AcqInfo.(['AICh' num2str(ii)]) = obj.AIChanList{ii};
                end

            else
                chanFields = fn(idxChan);

                chanNum = nan(size(chanFields));
                for ii = 1:numel(chanFields)
                    tok = regexp(chanFields{ii}, '^AICh(\d+)$', 'tokens', 'once', 'ignorecase');
                    assert(~isempty(tok), ...
                        'Unexpected AI channel field name "%s" found in acquisition metadata.', chanFields{ii});
                    chanNum(ii) = str2double(tok{1});
                end
                [~, ord] = sort(chanNum);
                chanFields = chanFields(ord);

                obj.AIChanList = cellfun(@(x) obj.AcqInfo.(x), chanFields, 'UniformOutput', false);
                obj.AIChanList = cellfun(@convertStringsToChars, obj.AIChanList, 'UniformOutput', false);
            end

            assert(isequaln(numel(obj.AIChanList), obj.AcqInfo.AINChannels), ...
                'Mismatch found between the AI channel list and the total number of channels.');

            assert(all(cellfun(@(x) ischar(x) && isrow(x), obj.AIChanList)), ...
                'All AI channel names must be character vectors.');

            assert(all(ismember(lower(obj.AIChanList), lower(obj.dictAIChan))), ...
                'Unexpected AnalogIN channel name found in acquisition metadata.');

            obj.sr = obj.AcqInfo.AISampleRate;

            if any(startsWith(fn, 'event', 'IgnoreCase', true))
                obj.b_isDigital = true;
                obj.trigChanName = 'StimDig';
                return
            end

            obj.b_isDigital = false;

            b_hasStim1 = any(cellfun(@any, regexpi(fn, 'stimulation1_')));
            b_hasStim2 = any(cellfun(@any, regexpi(fn, 'stimulation2_')));

            if ~any([b_hasStim1, b_hasStim2])
                return
            end

            names = {'StimAna1', 'StimAna2'};
            obj.trigChanName = names([b_hasStim1, b_hasStim2]);

            idxBurst = contains(fn, 'burst_delay', 'IgnoreCase', true);
            idxPeriod = contains(fn, '_period', 'IgnoreCase', true);
            idxFreq = contains(fn, '_frequency', 'IgnoreCase', true);

            if any(idxBurst) && obj.AcqInfo.(fn{find(idxBurst,1,'first')}) > 0
                burst_delay = obj.AcqInfo.(fn{find(idxBurst,1,'first')});
                stim_period = obj.AcqInfo.(fn{find(idxPeriod,1,'first')});
                obj.minInterStim = 1.15 * (burst_delay + stim_period) / 1000;

            elseif any(idxPeriod)
                obj.minInterStim = 1.15 * obj.AcqInfo.(fn{find(idxPeriod,1,'first')}) / 1000;

            elseif any(idxFreq)
                obj.minInterStim = 1.15 * (1 / obj.AcqInfo.(fn{find(idxFreq,1,'first')}));
            end
        end


        function [timestamps, state, trimInfo] = detectTrig(obj, data)
            %DETECTTRIG Detect trigger transitions from a 1-D signal.
            %
            % Input:
            %   data : numeric vector
            %       Trigger signal.
            %
            % Outputs:
            %   timestamps : single column vector
            %       Trigger transition timestamps in seconds.
            %   state : logical column vector
            %       Trigger state at each timestamp:
            %           true  = ON / onset
            %           false = OFF / offset
            %   trimInfo : struct
            %       Information about boundary-partial triggers handled during
            %       detection.
            %
            % Notes:
            %   - If trigThr = 'auto', threshold is set to 80% of the signal span above
            %     the minimum value.
            %   - trigPolarity controls whether onset is detected on an upward
            %     ('positive') or downward ('negative') threshold crossing.
            %   - Biphasic signals are rectified with abs().
            %   - Leading/trailing boundary-partial pulses generate warnings.
            %   - Ambiguous mid-stream mismatches are rejected with an error.

            timestamps = single.empty(0,1);
            state = false(0,1);
            trimInfo = struct( ...
                'nLeadDropped', 0, ...
                'nTrailDropped', 0, ...
                'bTrailingRiseKept', false, ...
                'trimApplied', false);

            validateattributes(data, {'single','double'}, {'vector','nonempty','real','finite'}, ...
                'detectTrig', 'data');

            data = double(data(:));

            % -------------------------------------------------------------
            % Robust biphasic classification
            % -------------------------------------------------------------
            % Do not classify a signal as biphasic only because a few noisy samples
            % go below zero. Instead, require meaningful positive and negative
            % excursions relative to the signal center.
            dataCenter = median(data);
            dataCentered = data - dataCenter;
            dataRange = max(dataCentered) - min(dataCentered);

            ampThr = max(0.02, 0.1 * max(dataRange, obj.minTrigAmp));
            minSideSamples = max(3, ceil(0.001 * numel(data))); % ~0.1% of samples, at least 3

            bHasPositivePhase = nnz(dataCentered > ampThr) >= minSideSamples;
            bHasNegativePhase = nnz(dataCentered < -ampThr) >= minSideSamples;

            bRectifiedBiphasic = bHasPositivePhase && bHasNegativePhase;

            if bRectifiedBiphasic
                disp('Biphasic signal detected.')
                data = abs(data);
            end

            % -------------------------------------------------------------
            % Resolve threshold
            % -------------------------------------------------------------
            if ischar(obj.trigThr) || isStringScalar(obj.trigThr)
                assert(strcmpi(obj.trigThr, 'auto'), ...
                    'Invalid trigger threshold mode. Use a numeric scalar or ''auto''.');
                trigAmp = 0.8 * (max(data) - min(data));
                if trigAmp < obj.minTrigAmp
                    warning(['Operation aborted. The trigger amplitude is too low for ' ...
                        'automatic detection. Manually set a threshold and try again.']);
                    return
                end
                obj.trigThr = trigAmp + min(data);
            else
                validateattributes(obj.trigThr, {'numeric'}, {'scalar','real','finite'}, ...
                    'detectTrig', 'obj.trigThr');
            end

            % -------------------------------------------------------------
            % Determine onset/offset crossings
            % -------------------------------------------------------------
            % For biphasic signals, rectification turns the problem into a
            % positive-polarity one.
            if bRectifiedBiphasic || strcmpi(obj.trigPolarity, 'positive')
                tmRise = find(data(1:end-1) < obj.trigThr & data(2:end) > obj.trigThr); % onset
                tmFall = find(data(1:end-1) > obj.trigThr & data(2:end) < obj.trigThr); % offset
            else
                tmRise = find(data(1:end-1) > obj.trigThr & data(2:end) < obj.trigThr); % onset
                tmFall = find(data(1:end-1) < obj.trigThr & data(2:end) > obj.trigThr); % offset
            end

            if isempty(tmRise) && isempty(tmFall)
                return
            end

            % Leading unmatched offset
            if ~isempty(tmFall) && (isempty(tmRise) || tmFall(1) < tmRise(1))
                trimInfo.nLeadDropped = 1;
                trimInfo.trimApplied = true;
                tmFall(1) = [];
                warning(['Signal appears to start while the trigger is already ON. ' ...
                    'Discarding the leading partial trigger event to preserve alignment.']);
            end

            % Trailing unmatched onset
            if ~isempty(tmRise) && (isempty(tmFall) || tmRise(end) > tmFall(end))
                if obj.b_isDigital || strcmpi(obj.trigType, 'EdgeToggle')
                    trimInfo.bTrailingRiseKept = true;
                    warning(['Signal appears to end before the final trigger returns to the idle state. ' ...
                        'Keeping the trailing onset. OFF time will be reconstructed from known durations ' ...
                        'for digital triggers or interpreted by toggle semantics.']);
                else
                    trimInfo.nTrailDropped = 1;
                    trimInfo.trimApplied = true;
                    tmRise(end) = [];
                    warning(['Signal appears to end before the final trigger returns to the idle state. ' ...
                        'Discarding the trailing partial trigger event to preserve alignment.']);
                end
            end

            % Validate mismatch after boundary handling
            if obj.b_isDigital || strcmpi(obj.trigType, 'EdgeToggle')
                bOk = (numel(tmRise) == numel(tmFall)) || (numel(tmRise) == numel(tmFall) + 1);
                assert(bOk, ...
                    ['Trigger edge mismatch is not limited to a simple boundary case. ' ...
                    'Detection aborted because event alignment would be ambiguous.']);
            else
                assert(numel(tmRise) == numel(tmFall), ...
                    ['Trigger edge mismatch is not limited to a simple boundary case. ' ...
                    'Detection aborted because event alignment would be ambiguous.']);
            end

            % -------------------------------------------------------------
            % Build raw transition list
            % -------------------------------------------------------------
            timestamps = single([tmRise; tmFall] ./ obj.sr);
            state = [true(numel(tmRise),1); false(numel(tmFall),1)];
            [timestamps, idx] = sort(timestamps);
            state = state(idx);

            % Toggle mode: use consecutive onset edges as ON/OFF states
            if strcmpi(obj.trigType, 'EdgeToggle')
                riseTimes = timestamps(state == 1);
                timestamps = riseTimes;
                state = true(numel(riseTimes),1);
                state(2:2:end) = false;
                return
            end

            % -------------------------------------------------------------
            % Burst consolidation for non-digital EdgeSet triggers only
            % -------------------------------------------------------------
            if ~obj.b_isDigital
                onTimes = timestamps(state == 1);
                stimLim = find(diff(onTimes) >= obj.minInterStim);
                nStim = numel(stimLim) + 1;

                if nStim < nnz(state)
                    disp('Burst stim detected.')

                    burstOff = [stimLim; numel(tmFall)];
                    burstOn = [1; stimLim + 1];

                    tmRise = tmRise(burstOn);
                    tmFall = tmFall(burstOff);

                    timestamps = single([tmRise; tmFall] ./ obj.sr);
                    state = [true(numel(tmRise),1); false(numel(tmFall),1)];

                    [timestamps, idx] = sort(timestamps);
                    state = state(idx);
                end
                return
            end

            % -------------------------------------------------------------
            % Digital stimulation: derive OFF timestamps from known durations
            % -------------------------------------------------------------
            fn = regexp(fieldnames(obj.AcqInfo),'Stim\d+','match','once');
            fn(cellfun(@isempty,fn)) = [];
            stimInfo = cellfun(@(x) obj.AcqInfo.(x), fn, 'UniformOutput', false);
            stimInfo = [stimInfo{:}];

            durationMap = containers.Map(num2cell([stimInfo.ID]), num2cell([stimInfo.Duration]));
            eventOrder = obj.AcqInfo.Events_Order(:);

            if trimInfo.nLeadDropped > 0
                eventOrder = eventOrder(trimInfo.nLeadDropped+1:end);
            end
            if trimInfo.nTrailDropped > 0
                eventOrder = eventOrder(1:end-trimInfo.nTrailDropped);
            end

            assert(numel(eventOrder) == nnz(state), ...
                ['The number of detected digital trigger onsets does not match ' ...
                'AcqInfo.Events_Order after boundary trimming.']);

            onTimes = timestamps(state == 1);
            durationOrder = zeros(numel(eventOrder), 1, 'single');
            for ii = 1:numel(eventOrder)
                assert(isKey(durationMap, eventOrder(ii)), ...
                    'Stimulus ID %d from Events_Order was not found in AcqInfo Stim definitions.', eventOrder(ii));
                durationOrder(ii) = single(durationMap(eventOrder(ii)));
            end

            timestamps = zeros(2*numel(onTimes), 1, 'single');
            state = false(2*numel(onTimes), 1);
            timestamps(1:2:end) = onTimes;
            timestamps(2:2:end) = onTimes + durationOrder;
            state(1:2:end) = true;

            [timestamps, idx] = sort(timestamps);
            state = state(idx);
        end

        function [eventID, eventNameList, eventSeqNameList] = readCSVfile(obj, colName)
            %READCSVFILE Read a CSV event file.
            %
            % If no column is provided, the first column is read. If multiple columns
            % are selected, their values are concatenated with a hyphen.
            %
            % Outputs:
            %   eventID          : duplicated event IDs (ON/OFF)
            %   eventNameList    : unique event names in stable order
            %   eventSeqNameList : event name sequence, one entry per trial/event

            opts = detectImportOptions(obj.EventFileName);

            if ~exist('colName','var') || isempty(colName)
                colName = {''};
            elseif ischar(colName)
                colName = {colName};
            end

            if ~isempty([colName{:}])
                opts.SelectedVariableNames = colName;
            else
                colName = {''};
                opts.DataLines(1) = 2;
            end

            opts.VariableTypes = repmat({'char'}, size(opts.VariableTypes));
            opts.MissingRule = 'omitrow';
            opts.ExtraColumnsRule = 'ignore';

            out = table2cell(readtable(obj.EventFileName, opts));

            eventSeqNameList = cell(size(out,1),1);
            for ii = 1:size(out,1)
                if ~isempty([colName{:}])
                    eventSeqNameList{ii} = strjoin(reshape(vertcat(colName, out(ii,:)), 1, []), '-');
                else
                    eventSeqNameList{ii} = strjoin(out(ii,:), '-');
                end
            end

            eventNameList = unique(eventSeqNameList, 'stable');
            nameMap = containers.Map(eventNameList, num2cell(1:numel(eventNameList)));

            eventID_on = zeros(numel(eventSeqNameList),1,'uint16');
            for ii = 1:numel(eventSeqNameList)
                eventID_on(ii) = uint16(nameMap(eventSeqNameList{ii}));
            end

            eventID = repelem(eventID_on, 2, 1);
        end

        function [eventID, eventNameList, eventSeqNameList] = readVpixxFile(obj)
            %READVPIXXFILE Extract event sequence from a .vpixx or .txt file.
            %
            % Outputs:
            %   eventID          : duplicated event IDs (ON/OFF)
            %   eventNameList    : unique event names in stable order
            %   eventSeqNameList : event name sequence, one entry per event/trial
            %
            % Notes:
            %   - The parser looks for the block between "RAW DATA" and "SORTED".
            %   - If those markers are missing or malformed, the function returns empty
            %     outputs instead of indexing invalid text ranges.

            eventID = [];
            eventNameList = {};
            eventSeqNameList = {};

            filetext = fileread(obj.EventFileName);

            idxStart = regexp(filetext, 'RAW DATA', 'end', 'once');
            idxStop = regexp(filetext, 'SORTED', 'start', 'once');

            if isempty(idxStart) || isempty(idxStop) || idxStop <= idxStart
                return
            end

            filetext = strtrim(filetext(idxStart+1:idxStop-1));
            if isempty(filetext)
                return
            end

            tab = strsplit(filetext, '\n')';

            outID = [];
            eventSeqNameList = {};

            myCols = [2 3];

            for i = 2:numel(tab)
                str = strsplit(tab{i}, '\t');

                if numel(str) < max(myCols)
                    continue
                end

                if isempty(str{1})
                    continue
                end

                condID = str2double(str{myCols(1)});
                stimName = str{myCols(2)};

                if isnan(condID) || isempty(stimName)
                    continue
                end

                outID(end+1,1) = condID; %#ok<AGROW>
                eventSeqNameList(end+1,1) = {stimName}; %#ok<AGROW>
            end

            if isempty(eventSeqNameList)
                return
            end

            eventNameList = unique(eventSeqNameList, 'stable');
            eventID_on = zeros(numel(eventSeqNameList), 1, 'uint16');

            for ii = 1:numel(eventNameList)
                eventID_on(strcmp(eventSeqNameList, eventNameList{ii})) = ii;
            end

            eventID = repelem(eventID_on, 2, 1);
        end

        function bOk = validateCondition(obj, conditionName)
            %VALIDATECONDITION Validate whether a condition exists and is not ignored.
            %
            % Input:
            %   conditionName : char | string scalar
            %       Condition name. Matching is case-insensitive.
            %
            % Output:
            %   bOk : logical
            %       True when the condition exists and at least one event from that
            %       condition is still selected.

            conditionName = convertStringsToChars(conditionName);
            validateattributes(conditionName, {'char'}, {'scalartext'}, 'validateCondition', 'conditionName');

            if isempty(obj.eventID)
                error('No events are currently available.')
            end

            if isempty(obj.selectedEvents) || ...
                    ~islogical(obj.selectedEvents) || ...
                    ~isequal(size(obj.selectedEvents), size(obj.eventID))
                warning(['selectedEvents was missing or invalid. ' ...
                    'Resetting it to all TRUE with the same size as eventID.']);
                obj.selectedEvents = true(size(obj.eventID));
            end

            bOk = false;
            indxCond = find(strcmpi(conditionName, obj.eventNameList), 1, 'first');

            if isempty(indxCond)
                error('The condition with name "%s" does not exist in the "eventNameList".', conditionName);
            end

            if all(~obj.selectedEvents(obj.eventID == indxCond))
                warning('The condition with name "%s" is already listed as ignored!', obj.eventNameList{indxCond})
                return
            end

            bOk = true;
        end

        function bOk = validateRepetition(obj, conditionName, repetitionIndex)
            %VALIDATEREPETITION Validate whether a repetition exists and is not ignored.
            %
            % Inputs:
            %   conditionName   : char | string scalar
            %       Condition name. Matching is case-insensitive.
            %   repetitionIndex : positive integer scalar
            %       Repetition index to validate.
            %
            % Output:
            %   bOk : logical
            %       True when the repetition exists and is still selected.

            conditionName = convertStringsToChars(conditionName);
            validateattributes(repetitionIndex, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
                'validateRepetition', 'repetitionIndex');

            bOk = false;

            if ~obj.validateCondition(conditionName)
                return
            end

            condID = find(strcmpi(conditionName, obj.eventNameList), 1, 'first');
            validReps = unique(obj.repetitionID(obj.state & obj.eventID == condID), 'stable');

            assert(ismember(repetitionIndex, validReps), ...
                'The repetition index %d is out of bounds for condition "%s".', ...
                repetitionIndex, obj.eventNameList{condID});

            repIdx = (obj.eventID == condID) & (obj.repetitionID == repetitionIndex);

            if all(~obj.selectedEvents(repIdx))
                warning('The repetition "%d" from condition "%s" is already listed as ignored!', ...
                    repetitionIndex, obj.eventNameList{condID});
                return
            end

            bOk = true;
        end

        function delete(obj)
            %DELETE Restore the original MATLAB warning state on object destruction.
            warning(obj.warnOrigState)
        end
    end
end
