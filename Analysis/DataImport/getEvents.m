function getEvents(RawFolder,SaveFolder, varargin)
% GETEVENTS detects events from LabeoTech Imaging system's analog channels
% (e.g., ai_00000.bin) and saves the event information to "events.mat" file
% in the SaveFolder to be used by other umIT functions.

% Defaults:
default_opts = struct('StimChannel','Internal-main', 'Threshold','auto','TriggerType','EdgeSet', 'minInterStimTime', 2,'ConditionFileType','none','ConditionFileName','auto','CSVColNames','all');
opts_values = struct('StimChannel',{{'Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}'},'Threshold',{{'auto',Inf}}, 'TriggerType', {{'EdgeSet', 'EdgeToggle'}},'minInterStimTime',[0.5, Inf],'ConditionFileType',{{'none','CSV','Vpixx'}}, 'ConditionFileName',{{'auto'}},'CSVColNames',{{'all'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
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
if isempty(opts.StimChannel)
    warning('No stim channels selected. Event file creation aborted.')
    return
end
% Check and prepare parameters for EventsManager object:
opts.ConditionFileName = strip(opts.ConditionFileName);
opts.CSVVolNames = strip(opts.CSVColNames);
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
evObj = EventsManager(RawFolder,opts.ConditionFileType);
% Update EventsManager object properties:
evObj.trigThr = opts.Threshold;
evObj.trigChanName = opts.StimChannel;
evObj.minInterStim = opts.minInterStimTime;
% Detect triggers:
evObj.getTriggers(true); % Verbose;
% Update event names:
evObj.readEventFile(opts.ConditionFileName,'CSVcols', opts.CSVColNames);
% Save events to file:
evObj.saveEvents(SaveFolder);
end