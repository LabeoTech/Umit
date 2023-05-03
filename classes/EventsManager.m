classdef EventsManager < handle
    % EVENTSMANAGER. A class to manage triggers and meta data (event IDs)
    % in order to create an "events.mat" file to be used by umIT.
    properties
        trigThr = 'auto'; % Trigger detection Threshold in volts (Default = 'auto'). The 'auto' means that the threshold will be at 70% of the signal amplitude.
        trigType char {mustBeMember(trigType,{'EdgeSet','EdgeToggle'})} = 'EdgeSet' % Trigger type.
        trigChanName = 'Internal-main'; % Default name of AI channel containing the triggers.
        minInterStim single {mustBeNonnegative} = 2 % Minimal inter stim time (in seconds). This param. is used to identify Bursts in analogIN. !This value should be higher than the inter-burst value to work!
    end
    properties (SetAccess = private)
        dictAIChan cell = {'CameraTrig','Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}; %  List of all available channels from Imaging system:
        DataFolder char % Folder containing the ai_xxxx.bin and (optionally) the meta data file.
        EventFileName char % Name of the meta data file containing events information (.csv, .txt, .vpixx ...).
        EventFileParseMethod char % Name of the method to read the Event File. {'default','vpixx'};
        AcqInfo struct % Content of the info.txt file as a structure.
        AnalogIN single % Array of analog input data.
        sr single {mustBePositive} = 10000; % Sample rate of analog input channels in Hz.
        timestamps single {mustBeNonnegative} % Time stamps in seconds of triggers.
        state uint8 {mustBeNonnegative} % State of triggers (1= ON, 0 = OFF)
        eventID uint16 {mustBePositive} % Index of each condition (1, 2, 3 ...)
        eventNameList cell % Name of each condition.
        b_isDigital logical = false; % TRUE for digital stimulation using OiS200 as master.
        minTrigAmp single {mustBePositive} = .15; % Minimal signal amplitude (in Volts). Used when the trigger threshold is "auto". Ignored otherwise.
    end
    
    methods
        function obj = EventsManager(DataFolder,varargin)
            % This is the constructor method.
            % It instantiates the object and sets the values for the
            % DataFolder and EventFileName properties.
            % It reads the analog IN channels and the info.txt as well.
            
            % Input validation:
            p = inputParser;
            addRequired(p,'DataFolder',@isfolder)
            addOptional(p,'EventFileName','',@(x) isfile(fullfile(DataFolder,x)));
            addParameter(p, 'ParseMethod','default',@(x) ismember(lower(x),{'default','vpixx'}));
            parse(p,DataFolder,varargin{:});
            % Set main properties:
            obj.DataFolder = p.Results.DataFolder;
            obj.EventFileName = p.Results.EventFileName;
            obj.EventFileParseMethod = p.Results.ParseMethod;
            % Set data from info.txt file:
            obj.setInfo;
            % Set Analog inputs:
            obj.setAnalogIN;   
            % Get event info if the stim is digital:
            if obj.b_isDigital
                obj.getTriggers;
            end
        end
        
        function out = get.trigChanName(obj)
            if ischar(obj.trigChanName)
                out = {obj.trigChanName};
            else
                out = obj.trigChanName;
            end
        end
                    
        function plotAnalogIN(obj, varargin)
            % PLOTANALOGIN plots the Analog input signals and overlays the
            % threshold as well as the detected triggers, if existent.
            p = inputParser;
            addRequired(p,'obj');
            addOptional(p, 'chanName','',@(x) isempty(x) || all(ismember(lower(x),lower(obj.dictAIChan))));
            parse(p,obj,varargin{:});
            chanName = p.Results.chanName;            
            if isempty(obj.AnalogIN)
                return
            end
            if isempty(chanName)
                chanIndx = 1:size(obj.AnalogIN,2);
            else
                [~,chanIndx] = ismember(chanName,obj.dictAIChan);
            end                         
            
            xVec = [0:size(obj.AnalogIN,1)-1]./obj.sr;% Use X-axis in seconds.
            axYSize = [min(obj.AnalogIN(:)), max(obj.AnalogIN(:))];
            % Get list of titles for each channel:
            fn = fieldnames(obj.AcqInfo);            
            fn = fn(startsWith(fn,'AICh', 'IgnoreCase',true));
            fn = fn(chanIndx);
            % Show 4 plot per figure:
            nFigs = ceil(mod(length(fn)/4,4));
            cnt = 1;   
            b_trigsPlotted= false;
            for ii = 1:nFigs
                f(ii) = figure('Name',sprintf('Analog Inputs %d/%d',ii,nFigs),...
                    'Visible','off','NumberTitle','off', 'Position',[0 0 560 720],...
                    'CreateFcn',{@movegui,'center'},'CloseRequestFcn', @closeAllFigs);                
                for jj = 1:4
                    if cnt > length(fn)
                        break
                    end
                    s(cnt) = subplot(4,1,jj,'Parent',f(ii));
                    % Plot analogIN traces:
                    plot(s(cnt),xVec,obj.AnalogIN(:,chanIndx(cnt)),'ko-','MarkerSize',2, 'Color',[.8 .8 .8], 'MarkerEdgeColor','k');
                    if jj == 1
                        % Set axes labels:
                        s(cnt).XLabel.String = 'time (s)';
                        s(cnt).YLabel.String = 'amp.(V)';
                    end
                    % Plot detected triggers and threshold lines,
                    if ~b_trigsPlotted
                        hold(s(cnt),'on');
                        % Plot threshold line:
                        if ~ischar(obj.trigThr)
                            ln = line([xVec(1) xVec(end)],[obj.trigThr obj.trigThr],'Color','r');
                            ln.Tag = 'thrLn';
                        end
                        % Plot Trigger patches:
                        if ~isempty(obj.timestamps)
                            idx = unique(obj.eventID);
                            % Create semi-transparent patches to represent HIGH state of triggers:                          
                            % Trigger color code
                            colorArr = jet(64);
                            colorArr = colorArr(round(linspace(1,64,numel(obj.eventNameList))),:);
                            for kk = 1:length(idx)
                                xOn = obj.timestamps(obj.eventID == idx(kk) & obj.state == 1);
                                xOff = obj.timestamps(obj.eventID == idx(kk) & obj.state == 0);
                                x = [xOn xOff xOff xOn];
                                y = repmat([axYSize(1) axYSize(1) axYSize(2) axYSize(2)], size(xOn));
                                ptc(kk) = patch(s(cnt), x',y',colorArr(kk,:), 'FaceAlpha', .5, 'EdgeColor', 'none', 'Tag','TrigPatch');
                            end
                            % Put legend on the first plot:
                            if jj == 1
                                legend(s(cnt),ptc,obj.eventNameList,'Location','best')
                            end
                        end
                        hold(s(cnt),'off');
                    else
                        % copy the content of the first axis
                        if exist('ln','var')
                            copyobj(ln,s(cnt));
                        end
                        if exist('ptc','var')
                            copyobj(ptc,s(cnt));
                        end
                    end
                    title(s(cnt), obj.AcqInfo.(fn{cnt}));
                    cnt = cnt+1;                    
                end                
            end
                        
            for ii = 1:length(f)
                f(ii).UserData = f;
                f(ii).Visible = 'on';
            end
            % Link all axes together
            linkaxes(s,'xy');
            set(s(1),'YLim',axYSize);
            
            % CloseFig callback:
            function closeAllFigs(src,~)
                h = src.UserData;
                delete(h)
            end            
        end
        
        function getTriggers(obj,varargin)
            % GETTRIGGERS detects the triggers from one or more analog IN channels.
            % It records the timestamps and state of each event.
            
            p = inputParser;
            addRequired(p,'obj');
            addOptional(p,'TriggerChannelName','Internal-main',@(x) ischar(x) || (iscell(x) && ischar(x{1})));            
            parse(p,obj,varargin{:})
            obj.trigChanName = p.Results.TriggerChannelName;
            
            % Reset event info:
            obj.timestamps = [];
            obj.state = [];
            obj.eventID = [];
            obj.eventNameList = {};            
            for ii = 1:length(obj.trigChanName)
                idxCh = strcmpi(obj.trigChanName{ii},obj.dictAIChan);
                [tmstmp,chanState] = obj.detectTrig(obj.AnalogIN(:,idxCh));               
                obj.timestamps = [obj.timestamps;tmstmp];
                obj.state = [obj.state;chanState];
                % Create dummy event ID:
                obj.eventID = [obj.eventID; ones(length(tmstmp),1).*ii];
                % Control for failed detections:
                if isempty(obj.timestamps)
                    warning(['Failed to detect triggers in channel ' obj.trigChanName{ii}])
                else
                    disp(['Triggers detected in channel ' obj.trigChanName{ii}]);
                end
            end 
            
            % Create dummy list of event names:
            obj.eventNameList = arrayfun(@num2str,unique(obj.eventID),'UniformOutput',false);
            
            [obj.timestamps,indx] = sort(obj.timestamps);
            obj.state = obj.state(indx);
            % Populate event ID and event Name lists:
            
            % For digital stimulation:
            if obj.b_isDigital
                % Overwrite dummy eventIDs with real ones from the info.txt.
                obj.eventID = ones(size(obj.state));   
                obj.eventID(obj.state == 1) = obj.AcqInfo.Events_Order;
                obj.eventID(obj.state == 0) = obj.AcqInfo.Events_Order;
                fn = regexp(fieldnames(obj.AcqInfo),'Stim\d+','match','once');
                fn(cellfun(@isempty,fn)) = [];                
                IDs = cellfun(@(x) obj.AcqInfo.(x).ID,fn);
                Names = cellfun(@(x) obj.AcqInfo.(x).Name,fn, 'UniformOutput',false);
                [~,idx] = sort(IDs);
                obj.eventNameList = Names(idx);
                disp('Trigger timestamps generated!');
                return
            end
            
            % For Analog stimulation:
            
            % Case # 1: Multiple channels:
            if length(obj.trigChanName) > 1
                % Keep dummy event IDs and create eventNameList with
                % channel names:
                [~,idx] = ismember(lower(obj.trigChanName), lower(obj.dictAIChan));
                obj.eventNameList = obj.dictAIChan(idx);                                                
                disp('Trigger timestamps generated!')
                return
            end                                         
            % Case # 2: Single channel with event list from text file:
            if ~isempty(obj.EventFileName)
                switch lower(obj.EventFileParseMethod)
                    case 'default'
                        [obj.eventID, obj.eventNameList] = obj.readDefaultFile;
                    case 'vpixx'
                        [obj.eventID, obj.eventNameList] = obj.readVpixxFile;
                    otherwise
                        error(['Unknown event file parsing method ' obj.EventFileParseMethod]);
                end
                % Repeat each event ID item to account for the offset (state == 0)
                if length(obj.eventID) < length(obj.state)
                    obj.eventID = repelem(obj.eventID,2);
                end
                disp(['Event ID list read from file "' obj.EventFileName '"']);
            end                  
        end               
    end
    
    methods (Access = private)
        
        function setInfo(obj)
            % SETINFO calls the "ReadInfoFile.m" function or simply loads
            % the AcqInfo.m file from the DataFolder and stores the
            % structure in the "AcqInfo" property.
            
            % Load existing info:
            if exist(fullfile(obj.DataFolder, 'AcqInfos.mat'), 'file')
                a = load(fullfile(obj.DataFolder, 'AcqInfos.mat'));
                obj.AcqInfo = a.AcqInfoStream;
            else
                % Read Info file:
                obj.AcqInfo = ReadInfoFile(obj.DataFolder);
            end
            % For retrocompatibility:
            if ~any(startsWith(fieldnames(obj.AcqInfo),'AICh', 'IgnoreCase',true))
                obj.AcqInfo.AICh1 = 'CameraTrig';
                obj.AcqInfo.AICh2 = 'StimAna1';
                for ii = 3:obj.AcqInfo.AINChannels
                    obj.AcqInfo.(['AICh' num2str(ii)]) = ['AI' num2str(ii-2)];
                end
            end
            
            % Update Analog input sample rate:
            obj.sr = obj.AcqInfo.AISampleRate;
            % Check if the stimulation is digital:
            if any(startsWith(fieldnames(obj.AcqInfo),'event','IgnoreCase',true))
                obj.b_isDigital = true;                
            end            
        end
        
        function setAnalogIN(obj)
            % SETANALOGIN reads the ai_xxxx.bin files and saves it in "AnalogIN" property.
            % List of analog files containing raw data:
            
            aiFilesList = dir(fullfile(obj.DataFolder,'ai_*.bin'));
            if isempty(aiFilesList)
                warning(['Analog Input files (ai_xxxx.bin) not found in "' obj.DataFolder '"'])
                return
            end            
            disp('Reading analog inputs...')
            % Opening of the files:
            obj.AnalogIN = [];
            for ind = 1:size(aiFilesList,1)
                data = memmapfile(fullfile(obj.DataFolder,aiFilesList(ind).name),...
                    'Offset', 5*4, 'Format', 'double', 'repeat', inf);
                tmp = data.Data;
                tmp = reshape(tmp, 1e4, obj.AcqInfo.AINChannels, []);
                tmp = permute(tmp,[1 3 2]);
                tmp = reshape(tmp,[],obj.AcqInfo.AINChannels);
                obj.AnalogIN = [obj.AnalogIN; tmp];
            end
            % Crop to first and last camera triggers:
            camT = diff(obj.AnalogIN(:,1) > 2.5); camT = [camT;NaN];
            camTOn = find(camT == 1,1,'first');
            camTOff = find(camT == -1,1,'last');
            obj.AnalogIN = obj.AnalogIN(camTOn:camTOff,:);
            disp('Done')
        end
        
        % Trigger detection methods        
        function [timestamps, state] = detectTrig(obj,data)
            % DETECTTRIGGERS detects the triggers from a given signal and
            % outputs the timestamps and state. THis function is called by
            % the method "getTriggers".
            
            % Input:
            %    data(1xN num. vector): vector containing the triggers.
            timestamps = single.empty();
            state = uint8.empty();
            % Check if "data" is a vector;
            validateattributes(data,{'single','double'},{'vector'});
            %
            if size(data,1) > size(data,2)
                data = data';
            end
            if strcmpi(obj.trigThr, 'auto')
                trigAmp = (.7*(max(data(:)) - min(data(:))));
                % Control for the minimal amplitude for trigger detection
                if trigAmp < obj.minTrigAmp
                    warning('Operation Aborted! Trigger Amplitude threshold too low! Manually set a threshold and try again.');
                    return
                end
                % Update trigger threshold:
                obj.trigThr = (.7*(max(data(:)) - min(data(:)))) + min(data(:));
            end
            
            % Find samples that cross the threshold (rising and falling):            
            tmRise = find(data < obj.trigThr & [data(2:end) nan] > obj.trigThr);
            tmFall = find(data > obj.trigThr & [data(2:end) nan] < obj.trigThr);            
            if isempty(tmRise)
                disp('Triggers not found!')
                return
            elseif numel(tmRise) ~= numel(tmFall)
                disp('Number of rising and falling edges are not equal!')
                return
            end
            timestamps = single([tmRise tmFall]./obj.sr);
            state = [ones(1,numel(tmRise), 'uint8') zeros(1,numel(tmFall), 'uint8')];
            % Sort arrays by time and flip:
            [timestamps,idx] = sort(timestamps);
            timestamps = timestamps';
            state = state(idx)';
            % For Toggle type triggers:
            if strcmpi(obj.trigType, 'edgetoggle')                
                % Use 2nd signal onset as the trial "OFF" state:
                timestamps = timestamps(state==1);
                % Overwite state:
                state = ones(sum(state),1,'uint8');
                state(2:2:end) = 0;
            end
            
            % Deal with bursts stimuli:
            StimLim = find(diff(timestamps(state ==1)) >= obj.minInterStim); %             
            NbStim = length(StimLim)+1;
            if NbStim < sum(state) && ~obj.b_isDigital
                disp('Burst stim detected!')                                                  
                BurstOff = [StimLim; length(tmFall)]; 
                BurstOn = [1; StimLim+1];
                % Update timestamps and state:                
                tmRise = tmRise(BurstOn);
                tmFall = tmFall(BurstOff); 
                % Recalculate the timestamps and state:
                timestamps = single([tmRise tmFall]./obj.sr);
                state = [ones(1,numel(tmRise), 'uint8') zeros(1,numel(tmFall), 'uint8')];
                % Sort arrays by time and flip:
                [timestamps,idx] = sort(timestamps);
                timestamps = timestamps';
                state = state(idx)';
            end                 
            % Update timestamps and states for digital stim. Digital stim
            % will only have a toggle marking the stim onset and the
            % duration (in seconds) is recorded in the "info.txt" file.
            if obj.b_isDigital
                disp('Digital stim found!')
                fn = regexp(fieldnames(obj.AcqInfo),'Stim\d+','match','once');
                fn(cellfun(@isempty,fn)) = [];  
                stimInfo = cellfun(@(x) obj.AcqInfo.(x), fn);
                durationMap = containers.Map([stimInfo.ID],[stimInfo.Duration]);
                durationOrder = arrayfun(@(x) durationMap(x), obj.AcqInfo.Events_Order);                
                timestamps = reshape([timestamps(state == 1)'; timestamps(state ==1)' + durationOrder],numel(state),[]); % Create timestamps to mark the end of the stim (trigger + duration_sec);               
            end                       
        end
        
        % Event File parsers:                
        function [eventID, eventNameList] = readDefaultFile(obj)
            % READDEFAULTFILE reads a .CSV file with a single column containing
            % the event names in the presentation order. (NO HEADER!)
            
            opts = delimitedTextImportOptions('Delimiter',',','VariableTypes','char', 'VariableNames','NAME');
            out = readtable(fullfile(obj.DataFolder ,obj.EventFileName),opts);                        
            eventNameList = unique(out.NAME); % Get list of event names in ascending order
            nameMap = containers.Map(eventNameList,1:length(eventNameList)); 
            eventID = zeros(height(out),1);
            for ii = 1:length(out.NAME)
                eventID(ii) = nameMap(out.NAME{ii});
            end 
        end
        
        function [eventID, eventNameList] = readVpixxFile(obj)
            % READVPIXXFILE extracts the condition ID, name and timestamps (if applicable)
            % from the RAW DATA section of a .txt or .vpixx file.
            
            % Extract only the part of the text containing the real sequence of
            % conditions (RAW DATA section):
            filetext = fileread(fullfile(obj.DataFolder,obj.EventFileName));
            [~,idx_start] = regexp(filetext, 'RAW DATA');
            idx_stop = regexp(filetext, 'SORTED');
            filetext = strip(filetext(idx_start+1:idx_stop-1));
            tab = strsplit(filetext, '\n')';
            % Get stimulus ID and order. Ignore Event and Time columns for now!
            out = {};
            myCols = [2 3]; % Keep just Condition and Stimulus.
            for i = 2:length(tab)
                str = strsplit(tab{i},'\t');
                if isempty(str{1})
                    % Skip Event rows.
                    continue
                end
                out = [out;str(myCols)];
            end            
            out(:,1) = cellfun(@str2double, out(:,1), 'UniformOutput',false);
            eventID = [out{:,1}]';
            IDs = unique(eventID);
            eventNameList = cell(length(IDs),1);            
            for ii = 1:length(IDs)
                idx = find(eventID == IDs(ii),1,'first');
                eventNameList(ii) = out(idx,2);
            end
        end
                
    end
end


