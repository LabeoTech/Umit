function varargout = ReadAnalogsIn(FolderPath, SaveFolder, Infos, stimChan)
out = [];
if( ~strcmp(FolderPath, filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

% List of analog files containing raw data:
aiFilesList = dir([FolderPath 'ai_*.bin']);

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

% Detect Triggers in each channel:
Stim = {};
for i = 1:length(stimChan)
    Stim{i} = detectTriggers(stimChan(i), Infos, AnalogIN);
end
disp('Checking stim info...')
idxMiss = cellfun(@(x) isequaln(sum(x),0), Stim);
if all(idxMiss)
    % If no stim is detected, save StimParameters.mat file with default
    % values:
    disp('No Stimulations detected. Resting State experiment?');    
    Stim = 0;
    save([SaveFolder filesep 'StimParameters.mat'], 'Stim');
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
end

% Local functions:
function Stim = detectTriggers(stimChan, Infos, AnalogIN)
Stim = 0;
% CamTrig is on the first channel:
CamTrig = find((AnalogIN(1:(end-1),1) < 1.25) & (AnalogIN(2:end,1) >= 1.25))+1;
% Detect Stimulation triggers in channel 2:
% StimTrig is on the second channel (except if slave):
if stimChan > 3
    % If the stim channel is external, set the amplitude as the half of the
    % signal amplitude:
    thr = min(AnalogIN(:,stimChan)) + ((max(AnalogIN(:,stimChan)) - min(AnalogIN(:,stimChan)))/2);
    % Also, we filter the signal to remove high-frequency noise. This is
    % common with photodiodes, for instance:
    f = fdesign.lowpass('N,F3dB', 4, 200, 10000); % Apply low-pass filter @200Hz to remove high-frequency noise.
    lpass = design(f,'butter');
    AnalogIN(:,stimChan) = filtfilt(lpass.sosMatrix, lpass.ScaleValues, AnalogIN(:,stimChan)')';    
elseif ( ~isfield(Infos, 'Stimulation1_Amplitude') )
    % Set threshold amplitude for internal channels to 2.5V when the amplitude
    % value is not available (retrocompatibility issue)
%     Infos.Stimulation1_Amplitude = 5;
    thr = 2.5;
else
    thr = Infos.Stimulation1_Amplitude/2;
end
  
StimTrig = find((AnalogIN(1:(end-1), stimChan) < thr) &...
    (AnalogIN(2:end, stimChan) >= thr))+1;
if isempty(StimTrig)
    if isfield(Infos,['AICh' num2str(stimChan)])
        str = Infos.(['AICh' num2str(stimChan)]);
    else
        str = num2str(stimChan);
    end
    disp(['Missing triggers in channel ' str '!'])
    return
end
% Add Stimulation field for retrocompatibility:
if( ~isfield(Infos, 'Stimulation') )
    Infos.Stimulation = 1;
end

if Infos.Stimulation == 1
    Period = median(StimTrig(2:end)-StimTrig(1:(end-1)))/Infos.AISampleRate;
    Freq = 1/Period;
    StimLim = find(diff(StimTrig)>20000);
    NbStim = length(StimLim)+1;
    if( NbStim == length(StimTrig) ) %Single Pulse trigged Stims
        StimLim = find((AnalogIN(1:(end-1), stimChan) >thr) &...
            (AnalogIN(2:end, stimChan) <= thr))+1;
        StimLength = mean(StimLim - StimTrig)./Infos.AISampleRate;
        if StimLength < (CamTrig(2) - CamTrig(1))/Infos.AISampleRate
            StimLength = 3*(CamTrig(2) - CamTrig(1))/Infos.AISampleRate;
            StimLim = StimLim + 3*(CamTrig(2) - CamTrig(1));
        end
        Stim = zeros(length(AnalogIN(:,stimChan)),1);
        for indS = 1:NbStim
            Stim(StimTrig(indS):StimLim(indS)) = 1;
        end
    else %Pulses train Stim
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
    NbStimAI = length(StimTrig);
    NbStimCycle = Infos.Stimulation_Repeat;
    NbStim = sum(~cellfun(@isempty, regexpi(fieldnames(Infos), 'stim\d{1}'))); % Search for field "Stim" + one digit;
    NbColIll = sum(startsWith(fieldnames(Infos), 'Illumination'));
    InterFrame = mean(diff(CamTrig));
    if (NbColIll == 1)
        expander = [zeros(1,ceil(NbColIll*InterFrame*1.1)) ones(1,ceil(NbColIll*InterFrame*1.1))]./ceil(NbColIll*InterFrame*1.1);
    else
        expander = [zeros(1,ceil((NbColIll-1)*InterFrame*1.1)) ones(1,ceil((NbColIll-1)*InterFrame*1.1))];
    end
    
    Stim = diff(AnalogIN(:,stimChan),1,1)>2.5;
    Stim = conv(Stim,expander,'same')>0.1;
    Stim = Stim(CamTrig);
    
    if( NbStimAI ~= NbStimCycle*NbStim )
        disp('Acquisition might have been stopped before the end. Not all stimulations were acquired!');
    end
    
    StimTrig = find(diff(Stim)>0.5)+1;
    StimIDs = [];
    StimDurations = [];
    for indS = 1:NbStim
        eval(['StimIDs = cat(1, StimIDs, Infos.Stim' int2str(indS) '.Code);']);
        eval(['StimDurations = cat(1, StimDurations, Infos.Stim' int2str(indS) '.Duration);']);
    end
    % Encode Stim "Codes" into Trigger timestamps:
    Stim = zeros(size(CamTrig),'single'); % Reset "Stim"
    for ind = 1:length(StimTrig)
        if isfield(Infos,'Events_Order')
            ID = Infos.Events_Order(ind);
        else
            ID = mod(ind-1, NbStim) + 1;
        end
        St = StimTrig(ind);
        En = round(StimDurations(ID).*Infos.FrameRateHz + St);
        Stim(St:En) = StimIDs(ID);
    end
end
end