function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
% RUN_IMAGESCLASSIFICATION calls the function
% IMAGESCLASSIFICATION from the IOI library (LabeoTech).

% Defaults:
default_Output = {'fluo_475.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; % This is here only as a reference for PIPELINEMANAGER.m. The real outputs will be stored in OUTFILE.
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_IgnoreStim', false, 'StimChannel','Internal-main');
opts_values = struct('BinningSpatial', 2.^[0:5], 'BinningTemp',2.^[0:5],'b_IgnoreStim',[false, true], 'StimChannel', {{'Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}'});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
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
% outFile = {};
clear p
%%%%
% Get existing ImagesClassification files in directory:
existing_ChanList  = dir(fullfile(SaveFolder,'*.dat'));
idxName = ismember({existing_ChanList.name}, default_Output);
existing_ChanList = existing_ChanList(idxName);
% Calls function from IOI library. Temporary for now.
ImagesClassification(RawFolder, SaveFolder, opts.BinningSpatial, opts.BinningTemp,...
    opts.b_IgnoreStim, 0, opts.StimChannel);
% Get only new files created during ImagesClassification:
chanList = dir(fullfile(SaveFolder,'*.dat'));
idx = ismember({chanList.name}, default_Output);
chanList = chanList(idx);
idxName = ismember({chanList.name}, {existing_ChanList.name});
idxDate = ismember([chanList.datenum], [existing_ChanList.datenum]);
idxNew = ~all([idxName; idxDate],1);
chanList = {chanList(idxNew).name};
% If there is Stimulation, add "eventID" and "eventNameList" to the output
% files of ImagesClassification.
if ~opts.b_IgnoreStim
    % Here the first channel fom "chanList" is chosen to retrieve the
    % "Stim" data:
    try
        stimInfo = load(fullfile(SaveFolder, 'StimParameters.mat'));
        if isempty(stimInfo) || sum(stimInfo.Stim) == 0
            warning('Stim signal not found! Skipped Event file creation.')
        else
            % This works for one channel for now:
            genEventsFromStimParameters(SaveFolder);
        end
    catch
        warning('Stim Parameters file not found! Skipped event file creation.')
    end
    
end
outFile = fullfile(SaveFolder, chanList);
end