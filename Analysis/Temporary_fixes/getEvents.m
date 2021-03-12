function outFile = getEvents(RawFolder, SaveFolder, opts, Output)
% GETEVENTS retrieves the timestamps of trigger signals in Analog inputs of
% LabeoTech systems.
%
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'RawFolder', @isfolder)% For Raw Folder as input
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure
default_opts = struct('EdgeDetection', 'rising', 'OptionN', 'value_OptionN'); % EdgeDetection can be {"rising", "falling", "toggle"}.
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File: 
default_Output = {'AnalogTimeStamps.dat'};
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
% Parse inputs:
parse(p,RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder; 
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%
cd(RawFolder);
aiFiles = dir('ai*.bin');
AnalogIN = [];
for ind = 1:size(aiFiles,1)
    data = memmapfile(aiFiles(ind).name, 'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.data;
    tmp = reshape(tmp,1e4, 11, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],11);
    AnalogIN = [AnalogIN; tmp];
end
% Try to find channel that most likely has the triggers in it:
STDev = std(AnalogIN(:,2:end), 0, 1);% exclude Cam triggers from search.
sigChan = find(STDev == max(STDev)) + 1;
txt = fileread('info.txt');
str = regexp(txt, 'Illumination\d+:');
nChan = numel(str);
thr = 2.0; % Detection threshold
idx_camT = (AnalogIN(:,1) < thr & [AnalogIN(2:end,1); NaN] > thr);
sig = AnalogIN(idx_camT, sigChan);
sig = sig(4:nChan:end);
% find edge;
switch opts.EdgeDetection
    case 'rising'
        event_stmp  = (sig < thr & [sig(2:end) ; NaN] > thr);
    case 'falling'
        event_stmp  = (sig > thr & [sig(2:end) ; NaN] < thr);
    case 'toggle'
        event_stmp  = (sig > thr);
end
save(fullfile(SaveFolder, 'events.mat'), 'event_stmp', 'opts.EdgeDetection');

% Output file names
outFile = Output;
end



