classdef EventsManager < handle
    % EVENTSMANAGER. A class to manage triggers and meta data (event IDs)
    % in order to create an "events.mat" file to be used by umIT.
    properties
        trigThr single {mustBeNumeric} = 2.5; % Trigger detection Threshold in volts (Default = 2.5 V).
        trigType char {mustBeMember(trigType,{'EdgeSet','EdgeToggle'})} = 'EdgeSet' % Trigger type.
        trigChanName char  = 'Internal-main'; % Default name of AI channel containing the triggers.
    end
    properties (SetAccess = private)
        dictAIChan cell = {'Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}; %  List of all available channels from Imaging system:
        DataFolder char % Folder containing the ai_xxxx.bin and (optionally) the meta data file.
        EventFileName char % Name of the meta data file containing events information (.csv, .txt, .vpixx ...).
        AcqInfo struct % Content of the info.txt file as a structure.
        AnalogIN single % Array of analog input data.
        sr = 10000; % Sample rate of analog input channels in Hz.
        timestamps single {mustBeNonnegative} % Time stamps in seconds of triggers.
        state uint8 {mustBeNonnegative}% State of triggers (1= ON, 0 = OFF)
        eventID uint16 {mustBePositive} % Index of each condition (1, 2, 3 ...)
        eventNameList cell% Name of each condition.
        minInterStim = 1.5 % Minimal inter stim time (in seconds). This param. is used to identify Bursts in analogIN. !This value should be higher than the inter-burst value to work!
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
            addOptional(p,'EventFileName','',@isfile);
            parse(p,DataFolder,varargin{:});
            % Set main properties:
            obj.DataFolder = p.Results.DataFolder;
            obj.EventFileName = p.Results.EventFileName;
            % Set data from info.txt file:
            obj.setInfo;
            % Set Analog inputs:
            obj.setAnalogIN;
            %
        end
        
        function plotAnalogIN(obj)
            % PLOTANALOGIN plots the Analog input signals and overlays the
            % threshold as well as the detected triggers, if existent.
            
            if isempty(obj.AnalogIN)
                return
            end
            % Show 4 plot per figure:
            nFigs = mod(size(obj.AnalogIN,2),4);
            cnt = 1;
            xVec = [0:size(obj.AnalogIN,1)-1]./obj.sr;% Use X-axis in seconds.
            axYSize = [min(obj.AnalogIN(:)), max(obj.AnalogIN(:))];
            % Get list of titles for each channel:
            fn = fieldnames(obj.AcqInfo);
            fn = fn(startsWith(fn,'AICh'));
            for ii = 1:nFigs
                f(ii) = figure('Name',sprintf('Analog Inputs %d/%d',ii,nFigs),...
                    'Visible','off','NumberTitle','off', 'Position',[0 0 560 720],...
                    'CreateFcn',{@movegui,'center'},'CloseRequestFcn', @closeAllFigs);
                for jj = 1:4
                    if cnt > size(obj.AnalogIN,2)
                        break
                    end
                    s(cnt) = subplot(4,1,jj,'Parent',f(ii));
                    % Plot analogIN traces:
                    plot(s(cnt),xVec,obj.AnalogIN(:,cnt),'ko-','MarkerSize',2, 'Color',[.8 .8 .8], 'MarkerEdgeColor','k');
                    if jj == 1
                        % Set axes labels:
                        s(cnt).XLabel.String = 'time (s)';
                        s(cnt).YLabel.String = 'amp.(V)';
                    end
                    % Plot detected triggers and threshold lines,
                    if cnt == 2
                        hold(s(cnt),'on');
                        % Plot threshold line:
                        ln = line([xVec(1) xVec(end)],[obj.trigThr obj.trigThr],'Color','r');
                        ln.Tag = 'thrLn';
                        % Plot Trigger patches:
                        if ~isempty(obj.timestamps)
                            idx = unique(obj.eventID);
                            % Create semi-transparent patches to represent HIGH state of triggers:                          
                            % Trigger color code
                            colorArr = jet(64);
                            colorArr = colorArr(round(linspace(1,32,numel(obj.eventNameList))),:);
                            for kk = 1:length(idx)
                                xOn = obj.timestamps(obj.eventID == idx(kk) & obj.state == 1);
                                xOff = obj.timestamps(obj.eventID == idx(kk) & obj.state == 0);
                                x = [xOn xOff xOff xOn];
                                y = repmat([axYSize(1) axYSize(1) axYSize(2) axYSize(2)], size(xOn));
                                ptc = patch(s(cnt), x',y',colorArr(kk,:), 'FaceAlpha', .5, 'EdgeColor', 'none', 'Tag','TrigPatch');
                            end
                        end
                        hold(s(cnt),'off');
                    elseif cnt > 2
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
            % GETTRIGGERS detects the triggers from the analog IN channels.
            % It records the timestamps and state of each event.
            
            
            [obj.timestamps, obj.state] = obj.detectTrig(obj.AnalogIN(:,2));
            obj.eventID = ones(size(obj.state));
            obj.eventNameList = {'1'};
                       
            
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
                return
            end
            % Read Info file:
            obj.AcqInfo = ReadInfoFile(obj.DataFolder);
            % Update Analog input sample rate:
            obj.sr = obj.AcqInfo.AISampleRate;
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
            camT = diff(AnalogIN(:,1) > 2.5); camT = [camT;NaN];
            camTOn = find(camT == 1,1,'first');
            camTOff = find(camT == -1,1,'last');
            AnalogIN = AnalogIN(camTOn:camTOff,:);
            disp('Done')
        end
        
        % Trigger detection methods
        
        function [timestamps, state] = detectTrig(obj,data)
            % DETECTTRIGGERS detects the triggers from a given signal and
            % outputs the timestamps and state. THis function is called by
            % the method "getTriggers".
            
            % Input:
            %    data(1xN num. vector): vector containing the triggers.
            
            p = inputParser;
            addRequired(p,'obj');
            addRequired(p, 'data', @(x) isnumeric(x) && isvector(x));
            parse(p,obj, data);
            % Initialize variables:
            data = p.Results.data;
            thr = obj.trigThr;
            clear p
            %
            if size(data,1) > size(data,2)
                data = data';
            end
            % Find samples that cross the threshold (rising and falling):            
            tmRise = find(data < thr & [data(2:end) nan] > thr);
            tmFall = find(data > thr & [data(2:end) nan] < thr);            
            timestamps = single([tmRise tmFall]./obj.sr);
            state = [ones(1,numel(tmRise), 'uint8') zeros(1,numel(tmFall), 'uint8')];
            % Sort arrays by time and flip:
            [timestamps,idx] = sort(timestamps);
            timestamps = timestamps';
            state = state(idx)';
            % For Toggle type triggers:
            if strcmpi(obj.trigType, 'edgetoggle')
                disp('Setting toggle...');
                % Use 2nd signal onset as the trial "OFF" state:
                timestamps = timestamps(state==1);
                % Overwite state:
                state = ones(sum(state),1,'uint8');
                state(2:2:end) = 0;
            end
            % Deal with bursts stimuli:
            StimLim = find(diff(timestamps(state ==1)) >= obj.minInterStim); %             
            NbStim = length(StimLim)+1;
            if NbStim < sum(state)
                disp('Burst stim detected!')
                BurstStim = find(diff(timestamps(state == 1)) < obj.minInterStim);
                
                plot(data);
                hold on
                plot(tmRise(StimLim), data(StimLim),'gx')
                plot(tmFall(StimLim),data(StimLim), 'go')
                plot(tmRise(BurstStim), data(BurstStim),'rx')
                plot(tmFall(BurstStim), data(BurstStim),'ro')
                
                
                
            end
            
            
            
            
            
        end
        
    end
end


