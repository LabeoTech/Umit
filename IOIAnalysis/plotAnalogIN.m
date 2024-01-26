function varargout = plotAnalogIN(AnalogIN,Infos,chanIndx)
% PLOTANALOGIN plots the Analog input signals
% Input:
%   AnalogIN(matrix): AnalogIN data (Time x Channel).
%   Infos (struct): Acquisition Info obtained from the info.txt file.
%   chanIndx(optional, array): Index or list of indices of the channel to plot. If not
%   provided, all channels will be plotted.
% Output (optional):
%   f (handle): array of figure handle(s).


if ~exist('chanIndx','var')
    chanIndx = 1:size(AnalogIN,2);
else
    % Avoid repetitions:
    chanIndx = unique(chanIndx,'stable');
end
% Get list of channel names from the Info file:
fn = fieldnames(Infos);
chanList = cellfun(@(x) Infos.(x), fn(startsWith(fn,'AICh')),'UniformOutput',false);
if isempty(chanList)
    chanList = arrayfun(@(x) sprintf('AICh%d',x), 1:size(AnalogIN,2),'UniformOutput',false);
end
chanList = chanList(chanIndx);
% Calculate the number of figures with a maximum of 4 subplots per figure:
nFigs = ceil(length(chanList)/4);
nAxPerFig = ones(1,nFigs)*4;
remAx = mod(length(chanList),4);
if remAx > 0
    nAxPerFig(end) = remAx;
end
cnt = 1;
% Set axes properties:
xVec = [0:size(AnalogIN,1)-1]./Infos.AISampleRate;% Use X-axis in seconds
axYSize = [min(AnalogIN(:)), max(AnalogIN(:))]; % Get min-max for Yscale

for ii = 1:nFigs
    f(ii) = figure('Name',sprintf('Analog Inputs %d/%d (downsampled to %d Hz)',ii,nFigs, Infos.AISampleRate/10),...
        'Visible','off','NumberTitle','off', 'Position',[0 0 560 nAxPerFig(ii)*200],...
        'CreateFcn',{@movegui,'center'},'CloseRequestFcn', @closeAllFigs);     
    for jj = 1:nAxPerFig(ii)
        s(cnt) = subplot(nAxPerFig(ii),1,jj, 'Parent',f(ii));
        s(cnt).XLabel.String = 'time (s)';
        s(cnt).YLabel.String = 'amp.(V)';
        % Plot analogIN traces (1/10 downsampled to save space):
        line(xVec(1:10:end),AnalogIN(1:10:end,chanIndx(cnt)),'LineStyle','-', 'Color',[.3 .3 .3],'Parent',s(cnt));
        if jj == 1
            % Set axes labels:
            s(cnt).XLabel.String = 'time (s)';
            s(cnt).YLabel.String = 'amp.(V)';
        end
        % copy the content of the first axis
        if exist('ln','var')
            copyobj(ln,s(cnt));
        end
        title(s(cnt), chanList{cnt});
        cnt = cnt+1;
    end
end
% Distribute figures on screen:
scrsz = get(0,'screensize');
totalFigWidth = 560*nFigs + 10*(nFigs-1);
f(1).Position(1) = (scrsz(3) - totalFigWidth)/2;
for ii = 2:nFigs
    f(ii).Position(1) = sum([f(ii-1).Position([1 3]),10]);    
end  
arrayfun(@(x) set(x, 'Visible','on','UserData',f),f);
% Link all axes together
linkprop(s,{'XLim','YLim'});
set(s(1),'YLim',axYSize);

if nargout
    varargout{:} = f;
end
% CloseFig callback:
    function closeAllFigs(src,~)
        % CLOSEALLFIGS closes all figures when the selected figure is closed.        
        delete(src.UserData)
    end
end