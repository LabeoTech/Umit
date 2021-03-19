function out = get_timestampsFromSingleChannel(RawFolder, SaveFolder, varargin)
% GET_TIMESTAMPSFROMSINGLECHANNEL creates a TTL_events .MAT file in SAVEFOLDER
% from Analog signals recorded from ONE CHANNEL of LabeoTech Imaging systems.
% Inputs:
%   RawFolder: directory containing ai_xxxx.bin files.
%   SaveFolder: directory to save .MAT eventsfile.
%       optional parameters:
%       opts.Channel = Analog channel index with the TTL signal. If not
%       provided, the function will try to find the one with the largest Standard
%       Deviation value.
%       opts.threshold = Threshold  to be used in the signal detection.
%       Default = 2.5 (v).
% Output:
%   TTL_events.mat file containing channel ID, state and timestamps. 
%   For details, see function CREATE_TTL_EVENTSFILE.m.
%
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'RawFolder', @isfolder)% For Raw Folder as input
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('channel', -1, 'threshold', 2.5);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File: 
default_Output = ''; % There is no Output file for this function.
% Parse inputs:
parse(p,RawFolder, SaveFolder, varargin{:});
% Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
%%%%
cd(RawFolder)
txt = fileread(fullfile(RawFolder, 'info.txt'));
sr = regexp(txt, '(?<=AISampleRate:\s*)\d+', 'match', 'once');sr = str2double(sr);
tAIChan = regexp(txt, '(?<=AINChannels:\s*)\d+', 'match', 'once'); tAIChan = str2double(tAIChan);
aiFilesList = dir('ai_*.bin'); 
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile(aiFilesList(ind).name, 'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, sr, tAIChan, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],tAIChan);
    AnalogIN = [AnalogIN; tmp];
end
% Crop to first and last camera triggers:
camT = diff(AnalogIN(:,1) > opts.threshold); camT = [camT;NaN];
camTOn = find(camT == 1,1,'first');
camTOff = find(camT == -1,1,'last');
AnalogIN = AnalogIN(camTOn:camTOff,:);

if opts.channel == -1
    STDev = std(AnalogIN(:,2:end), 0, 1);% exclude Cam triggers from search.
    sigChan = find(STDev == max(STDev)) + 1;
    sigChan = sigChan(1); % if find more than one, pick the first (arbitrary choice here...)
else
    sigChan = opts.channel;
end
% create TTL_events file:
create_TTL_eventsFile(SaveFolder, AnalogIN(:,sigChan),sr,opts.threshold)
end



