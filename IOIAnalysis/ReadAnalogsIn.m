function out = ReadAnalogsIn(FolderPath, SaveFolder, Infos, stimChan)

% Instantiate some variableS:
Stim = 0;
StimLength = 0;
NbStim = 0;
InterStim_min = 0;
InterStim_max = 0;

%
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

% CamTrig is on the first channel:
CamTrig = find((AnalogIN(1:(end-1),1) < 1.25) & (AnalogIN(2:end,1) >= 1.25))+1;
% Detect Stimulation triggers in channel 2:
% StimTrig is on the second channel (except if slave):
if( ~isfield(Infos, 'Stimulation1_Amplitude') )
    Infos.Stimulation1_Amplitude = 5;
end
StimTrig = find((AnalogIN(1:(end-1), stimChan) < Infos.Stimulation1_Amplitude/2) &...
    (AnalogIN(2:end, stimChan) >= Infos.Stimulation1_Amplitude/2))+1;
% If no stim is detected, save StimParameters.mat file with default
% values:
if isempty(StimTrig)
    disp('No Stimulations detected. Resting State experiment?');    
    save([SaveFolder filesep 'StimParameters.mat'], 'CamTrig', 'Stim',...
        'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
    return
end

if Infos.Stimulation == 1
    Period = median(StimTrig(2:end)-StimTrig(1:(end-1)))/Infos.AISampleRate;
    Freq = 1/Period;
    Width = sum(AnalogIN(StimTrig(1):StimTrig(2),stimChan) > 2.5)...
        /(Period*Infos.AISampleRate);    
    StimLim = find(diff(StimTrig)>20000);
    NbStim = length(StimLim)+1;
    if( NbStim == length(StimTrig) ) %Single Pulse trigged Stims
        StimLim = find((AnalogIN(1:(end-1), stimChan) >Infos.Stimulation1_Amplitude/2) &...
            (AnalogIN(2:end, stimChan) <= Infos.Stimulation1_Amplitude/2))+1;
        StimLength = mean(StimLim - StimTrig)./Infos.AISampleRate;
        if StimLength < (CamTrig(2) - CamTrig(1))/Infos.AISampleRate
            StimLength = 3*(CamTrig(2) - CamTrig(1))/Infos.AISampleRate;
            StimLim = StimLim + 3*(CamTrig(2) - CamTrig(1));
        end
        InterStim_min = min((StimTrig(2:end) - StimLim(1:(end-1)))./10000);
        InterStim_max = max((StimTrig(2:end) - StimLim(1:(end-1)))./10000);
        Stim = zeros(length(AnalogIN(:,stimChan)),1);
        for indS = 1:NbStim
           Stim(StimTrig(indS):StimLim(indS)) = 1; 
        end
    else %Pulses train Stim
        StimLength = round(length(StimTrig)/(NbStim*Freq));
        InterStim_min = min((StimTrig(StimLim + 1) - StimTrig(StimLim))./10000);
        InterStim_max = max((StimTrig(StimLim + 1) - StimTrig(StimLim))./10000);
        InterStim_min = InterStim_min - StimLength;
        InterStim_max = InterStim_max - StimLength;
    
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
    
    StimLength = 2*InterFrame;
    NbStim = NbStimAI;
    StStart = find(diff(AnalogIN(:,stimChan)) > Infos.Stimulation1_Amplitude/2);
    InterStim_min = mean(StStart(2:end) - StStart(1:(end-1)))/1e4;
    InterStim_max = InterStim_min;    
end
% Save Stim parameters:
save([SaveFolder filesep 'StimParameters.mat'],'CamTrig', 'Stim', 'StimLength',...
    'NbStim', 'InterStim_min', 'InterStim_max');
end

