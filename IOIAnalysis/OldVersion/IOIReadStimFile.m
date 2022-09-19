function IOIReadStimFile(ExpeFolder)
AnalogIN = load([ExpeFolder filesep 'IOI_aux.mat']);
CamTrig = find(diff(AnalogIN.aux(:,2))>10000);

Stim = zeros(size(CamTrig));
StimLength = 0;
NbStim = 0;

StimON = find(diff(AnalogIN.aux(:,1))>=7500);
if( ~isempty(StimON) )
    Period = (StimON(2)-StimON(1))/10000;
    Freq = 1/Period;
    Width = sum(AnalogIN.aux(StimON(1):StimON(2),1)>1000)/1000;
    
    StimLim = find(diff(StimON)>50000);
    NbStim = length(StimLim)+1;
    StimLength = round(length(StimON)/(NbStim*Freq));
    InterStim_min = min((StimON(StimLim + 1) - StimON(StimLim))./10000);
    InterStim_max = max((StimON(StimLim + 1) - StimON(StimLim))./10000);
    InterStim_min = InterStim_min - StimLength;
    InterStim_max = InterStim_max - StimLength;
    
    Stim = zeros(length(AnalogIN.aux(:,1)),1);
    Stim(StimON(1):StimON(StimLim(1))) = 1;
    for indS = 2:length(StimLim)
        Stim(StimON(StimLim(indS-1)+1):StimON((StimLim(indS)))) = 1;
    end
    Stim(StimON(StimLim(end)+1):StimON(end)) = 1;
    
    Stim = Stim(CamTrig);
    
end
save([ExpeFolder filesep 'StimParameters.mat'], 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
end