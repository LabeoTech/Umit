function varargout = ReadAnalogsIn(FolderPath, SaveFolder, Infos, chanName, trigPolarity, b_LPF)
% ReadAnalogsIn  Read and process analog input recordings to detect triggers.
%
%   out = ReadAnalogsIn(FolderPath, SaveFolder, Infos, chanName, trigPolarity, b_LPF)
%   reads analog binary files (ai_*.bin) from FolderPath, detects camera and
%   stimulation triggers, and saves stimulation parameters to
%   SaveFolder/StimParameters.mat. Infos is a struct with recording metadata
%   (e.g., AISampleRate, AINChannels, FrameRateHz). chanName specifies the
%   stimulation channel (e.g., 'Internal-Main', 'Internal-Aux' or a custom
%   channel label). trigPolarity can be 'positive' or 'negative'. b_LPF is
%   an optional boolean to apply a low-pass filter to external stim
%   channels.
%
%   [out] = ReadAnalogsIn(...) returns a struct with detected stimulation
%   information and frame rate. If no stimulations are detected, a default
%   StimParameters.mat is saved and out is empty.
%
%   Notes:
%     - The function expects analog files named ai_*.bin with a 5*4 byte
%       header and doubles stored thereafter, organized per Infos.AINChannels.
%     - If chanName is empty or 'Internal-Main', the internal main channel is
%       used. If the provided chanName is not found, defaults to internal main.
out = [];
if( ~strcmp(FolderPath, filesep) )
    FolderPath = strcat(FolderPath, filesep);
end
if ~exist("b_LPF",'var')
    b_LPF = false;
end
% List of analog files containing raw data:
aiFilesList = dir([FolderPath 'ai_*.bin']);
% Check if all files exist:
aiFileNames = sort({aiFilesList.name})';
aiFileIndx = erase(aiFileNames,'ai_');aiFileIndx = erase(aiFileIndx,'.bin');
aiFileIndx = str2double(aiFileIndx);
if ~strcmpi(aiFileNames{1}, 'ai_00000.bin') || any(diff(aiFileIndx)~=1)
    error('Analog IN binary files missing! Impossible to read Analog In. Please ensure that all files are available and try again.')
end
% Opening of the files:
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderPath aiFilesList(ind).name],...
        'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, Infos.AINChannels, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],Infos.AINChannels);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data aiFilesList;

% Get list of channels from info.txt file:
fn = fieldnames(Infos);fn = fn(startsWith(fn,'AICh','IgnoreCase',true));
if ~isempty(fn)
    chanNameList = cellfun(@(x) Infos.(x),fn,'UniformOutput',false);
    % Translate the generic names of the internal Analog channels to the
    % actual ones
    if any(startsWith(chanName,'Internal','IgnoreCase',true))
        nameTranslator = containers.Map({'internal-main','internal-aux'},chanNameList([2 3]));
        for ii = 1:length(chanName)
            if isKey(nameTranslator, lower(chanName{ii}))
                chanName{ii} = nameTranslator(lower(chanName{ii}));
            end
        end
    end

else
    % FALLBACK - Consider the following channel organization when there are no channel names in the info.txt file:
    chanNameList = {'CameraTrig','Internal-Main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8','StimDig'}; % List of existing Analog channel names.
end

% Get channel indices list
[~,stimChan] = ismember(upper(chanName), upper(chanNameList));
for ii = 1:length(stimChan)
    if stimChan(ii) == 0
        warning('Invalid channel name "%s"! The "Internal-Main" channel will be read instead.', chanName{ii});
        stimChan(ii) = 2;
    end
end
% Detect Triggers in each channel:
% save current warning modes
prev = warning;

% disable stack trace and verbose suppression message
warning off backtrace
warning off verbose

Stim = {};
for i = 1:length(stimChan)
    Stim{i} = detectTriggers(stimChan(i), Infos, AnalogIN, trigPolarity,b_LPF);
end


% Be sure to delete existing StmParameters.mat file
if isfile(fullfile(SaveFolder,'StimParameters.mat'))
    delete(fullfile(SaveFolder,'StimParameters.mat'));
end
disp('Checking stim info...')
idxMiss = cellfun(@(x) isequaln(sum(x),0), Stim);
if all(idxMiss)
    % If no stim is detected, save StimParameters.mat file with default
    % values:
    disp('No Stimulations detected. Resting State experiment?');
    Stim = 0;
    save([SaveFolder filesep 'StimParameters.mat'], 'Stim');
    if nargout
        varargout{:} = out;
    end
    % restore previous warning modes
    warning(prev)
    return
end
% Remove data with missing triggers:
Stim(idxMiss) = [];%#ok; It is used in an eval fcn below.
stimChan(idxMiss) = [];
out = struct();
for i = 1:length(stimChan)
    if isfield(Infos,['AICh' num2str(stimChan(i))])
        chanName = Infos.(['AICh' num2str(stimChan(i))]);
    else
        chanName = num2str(stimChan(i));
    end
    v = genvarname(['Stim_' chanName]);
    eval(['out.' v ' = Stim{i};']);
end
out.Stim = 1; % Indicates that stim triggers were found.
out.FrameRateHz = Infos.FrameRateHz; % Add recording frame rate to StimParameters file.
% Save Stim parameters:
save([SaveFolder filesep 'StimParameters.mat'], '-struct', 'out');
disp('StimParameters saved!')
if nargout
    varargout{:} = out;
end
% restore previous warning modes
warning(prev)
end

% Local functions:
function Stim = detectTriggers(stimChan, Infos, AnalogIN, trigPolarity,b_LPF)
Stim = 0;
% CamTrig is on the first channel:
CamTrig = find((AnalogIN(1:(end-1),1) < 1.25) & (AnalogIN(2:end,1) >= 1.25))+1;
if strcmpi(trigPolarity,'negative')
    % Flip signal if it's negative:
    sigMin = min(AnalogIN(:,stimChan));
    sigMax = max(AnalogIN(:,stimChan));
    AnalogIN(:,stimChan) = (1 - (AnalogIN(:,stimChan) - sigMin)./(sigMax - sigMin)).*...
        (sigMax - sigMin) + sigMin;
end
% Detect Stimulation triggers in channel 2:
% StimTrig is on the second channel (except if slave):
if stimChan > 3

    if b_LPF
        % Also, we filter the signal to remove high-frequency noise. This is
        % common with photodiodes, for instance:

        f = fdesign.lowpass('N,F3dB', 4, 200, 20000); % Apply low-pass filter @200Hz to remove high-frequency noise.
        lpass = design(f,'butter');

        AnalogIN(:,stimChan) = filtfilt(lpass.sosMatrix, lpass.ScaleValues, AnalogIN(:,stimChan)')';
    end

    % If the stim channel is external, set the amplitude as the half of the
    % signal amplitude:
    minThr = 0.15; % Minimal threshold value for detection.
    sigAmp = max(AnalogIN(:,stimChan)) - min(AnalogIN(:,stimChan));
    if sigAmp > minThr
        thr = min(AnalogIN(:,stimChan)) + (sigAmp/2);
    else
        % This will force the detection to fail (not ideal).
        thr = minThr;
    end

elseif ( ~isfield(Infos, 'Stimulation1_Amplitude') )
    % Set threshold amplitude for internal channels to 2.5V when the amplitude
    % value is not available (retrocompatibility issue)
    %     Infos.Stimulation1_Amplitude = 5;
    thr = 2.5;

else
    %     thr = Infos.Stimulation1_Amplitude/2;
    thr = 1;
end

% Detect trigger rising edges:
StimTrig = find((AnalogIN(1:(end-1), stimChan) < thr) &...
    (AnalogIN(2:end, stimChan) >= thr))+1;
StimTrigOff = find((AnalogIN(1:(end-1), stimChan) >thr) &...
    (AnalogIN(2:end, stimChan) <= thr))+1;

if isfield(Infos,['AICh' num2str(stimChan)])
    chanName = Infos.(['AICh' num2str(stimChan)]);
else
    chanName = num2str(stimChan);
end

if numel(StimTrig) ~= numel(StimTrigOff)
    % Raise a warning if the number of rising and falling edges are
    % not equal:
    warning('Failed to detect triggers in channel "%s"! The number of rising and falling edges are not equal!',chanName)
    return
end
if isempty(StimTrig)
    warning('Missing triggers in channel "%s"!',chanName)
    return
end
% Add Stimulation field for retrocompatibility:
if( ~isfield(Infos, 'Stimulation') || Infos.Stimulation == 0 )
    Infos.Stimulation = 1;
end

if Infos.Stimulation == 1
    Period = median(StimTrig(2:end)-StimTrig(1:(end-1)))/Infos.AISampleRate;
    StimLim = find(diff(StimTrig)>20000); % Force minimum 2s interstim.
    NbStim = length(StimLim)+1;
    if( NbStim == length(StimTrig) ) %Single Pulse trigged Stims
        StimLim = StimTrigOff; % Update StimLim
        if mean(StimLim - StimTrig)./Infos.AISampleRate < (CamTrig(2) - CamTrig(1))/Infos.AISampleRate
            StimLim = StimLim + 3*(CamTrig(2) - CamTrig(1));
        end
        Stim = zeros(length(AnalogIN(:,stimChan)),1);
        for indS = 1:NbStim
            Stim(StimTrig(indS):StimLim(indS)) = 1;
        end
    else
        % Pulses train Stim
        Stim = zeros(length(AnalogIN(:,stimChan)),1);
        if( NbStim > 1 )
            Stim(StimTrig(1):StimTrig(StimLim(1))) = 1;
            for indS = 2:length(StimLim)
                Stim(StimTrig(StimLim(indS-1)+1):StimTrig((StimLim(indS)))) = 1;
            end
            Stim(StimTrig(StimLim(end)+1):StimTrig(end)) = 1;
        else
            Stim(StimTrig(1):StimTrig(end)) = 1;
        end
    end
    Stim = Stim(CamTrig);
elseif Infos.Stimulation == 2
    if b_LPF
        warning('Low-Pass Filter not applied to Digital Filter');
    end
    NbStimAI    = length(StimTrig);
    NbStimCycle = Infos.Stimulation_Repeat;
    NbStim      = sum(~cellfun(@isempty, regexpi(fieldnames(Infos), 'stim\d{1}')));
    NbColIll    = sum(startsWith(fieldnames(Infos), 'Illumination'));

    % -------------------------------------------------------------
    % Pulse-onset detection
    % -------------------------------------------------------------

    % 1) Threshold AI signal (preserve pulse width)
    StimAI = AnalogIN(:,stimChan) > thr;

    % Optional smoothing if noisy (uncomment if needed)
    % StimAI = movmean(AnalogIN(:,stimChan),5) > thr;

    % 2) Detect rising edges at AI resolution
    StimTrigAI = find(diff(StimAI) == 1) + 1;

    if isempty(StimTrigAI)
        warning('No stimulation onsets detected after thresholding.')
        Stim = zeros(length(CamTrig),1,'single');
        return
    end

    % 3) Map AI indices to imaging frame indices
    % Each frame k spans [CamTrig(k), CamTrig(k+1))
    frameIdx = discretize(StimTrigAI, [CamTrig; Inf]);

    % Remove invalid mappings
    frameIdx = frameIdx(~isnan(frameIdx));

    % Frame-level onset vector
    StimFrameOnsets = zeros(length(CamTrig),1,'single');

    % If multiple pulses fall inside same frame, keep first only
    uniqueFrames = unique(frameIdx,'stable');
    StimFrameOnsets(uniqueFrames) = 1;

    % -------------------------------------------------------------
    % Sanity check
    % -------------------------------------------------------------
    if NbStimAI ~= NbStimCycle*NbStim
        disp('Acquisition might have been stopped before the end. Not all stimulations were acquired!');
    end

    % -------------------------------------------------------------
    % Encode Stim "Codes" into frame timestamps
    % -------------------------------------------------------------

    StimTrig = find(StimFrameOnsets > 0);

    StimIDs = [];
    StimDurations = [];
    for indS = 1:NbStim
        eval(['StimIDs = cat(1, StimIDs, Infos.Stim' int2str(indS) '.ID);']);
        eval(['StimDurations = cat(1, StimDurations, Infos.Stim' int2str(indS) '.Duration);']);
    end

    % Reset Stim vector (frame resolution)
    Stim = zeros(length(CamTrig),1,'single');

    for ind = 1:length(StimTrig)

        if isfield(Infos,'Events_Order')
            ID = Infos.Events_Order(ind);
        else
            ID = mod(ind-1, NbStim) + 1;
        end

        St = StimTrig(ind);
        En = round(StimDurations(ID) * Infos.FrameRateHz + St);

        % Clip to valid frame range
        En = min(En, length(Stim));

        Stim(St:En) = StimIDs(ID);
    end

end
end