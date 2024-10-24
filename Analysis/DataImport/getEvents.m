function getEvents(RawFolder,SaveFolder, varargin)
% GETEVENTS detects events from LabeoTech Imaging system's analog channels
% (e.g., ai_0000x.bin) and saves the event information to "events.mat" file
% in the SaveFolder to be used by other umIT functions.

% Defaults:
default_opts = struct('StimChannel','Internal-main', 'Threshold','auto','TriggerType','EdgeSet', 'minInterStimTime', 2,'ConditionFileType','none','ConditionFileName','auto','CSVColNames','all','baselinePeriod','auto');
opts_values = struct('StimChannel',{{'Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}'},'Threshold',{{'auto',Inf}}, 'TriggerType', {{'EdgeSet', 'EdgeToggle'}},'minInterStimTime',[0.5, Inf],'ConditionFileType',{{'none','CSV','Vpixx'}}, 'ConditionFileName',{{'auto'}},'CSVColNames',{{'all'}},'baselinePeriod', {{'auto',Inf}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
% Arguments validation:
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p, RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
%%%%
if ischar(opts.StimChannel)
    opts.StimChannel = {opts.StimChannel};
end
if all(cellfun(@isempty,opts.StimChannel))
    warning('No stim channels selected. Event file creation aborted.')
    return
end
% Check and prepare parameters for EventsManager object:
opts.ConditionFileName = strip(opts.ConditionFileName);
opts.CSVColNames = strip(opts.CSVColNames);
if strcmpi(opts.ConditionFileName, 'auto') && ~strcmpi(opts.CSVColNames,'all') && ~strcmpi(opts.ConditionFileType,'CSV')
    % Force file type to 'CSV' if column names were set
    opts.ConditionFileType = 'CSV';
    warning('Condition file type set to "CSV" given that specific columns were stated in parameters!');    
end

if strcmpi(opts.CSVColNames, 'all')
    opts.CSVColNames = {''};
else    
    opts.CSVColNames = strsplit(opts.CSVColNames,',');
end

if strcmpi(opts.ConditionFileName, 'auto')
    opts.ConditionFileName = '';
end

% Instantiate EventsManager class:
evObj = EventsManager(SaveFolder,RawFolder,opts.ConditionFileType);

for ii = 1:length(opts.StimChannel)
    % Update internal channel names. These may change depending on the OiS
    % acquisition software version.
    if strcmpi(opts.StimChannel{ii}, 'internal-main')
        warning(['Translated ''Internal-Main'' channel to ''' evObj.AIChanList{2} '''']);
        opts.StimChannel{ii} = evObj.AIChanList{2};
    elseif strcmpi(opts.StimChannel{ii}, 'internal-aux')
        warning(['Translated ''Internal-Aux'' channel to ''' evObj.AIChanList{3} '''']);
        opts.StimChannel{ii} = evObj.AIChanList{3};
    end
end
% Update EventsManager object properties:
evObj.trigThr = opts.Threshold;
evObj.trigType = opts.TriggerType;
if ~isempty(evObj.trigChanName{:})
    evObj.trigChanName = opts.StimChannel;
end
evObj.minInterStim = opts.minInterStimTime;
% Detect triggers:
evObj.getTriggers(evObj.trigChanName,true); % Verbose
% Update event names:
evObj.readConditionFile(opts.ConditionFileName,'CSVcols', opts.CSVColNames);
% Set pre and post event times:
if ~strcmpi(opts.baselinePeriod,'auto')
    evObj.setBaselinePeriod(opts.baselinePeriod);
end    
% Save events to file:
evObj.saveEvents(SaveFolder);
end