function out = DataNavigator(varargin)
out = '0';
%%%%%%%
% Vars init
%%%%%%%
h.paths.FolderName = '';
h.flags.saveROIS = false; %flag to know if any changes were made to ROIs
h.flags.saveEvnts = false;
h.flags.IsThereHbO = false;
h.flags.IsThereHbR = false;
h.flags.IsThereHbT = false;
h.flags.IsThereFlow = false;
h.flags.IsThereFluo = false;
h.flags.IsThereGreen = false;
h.flags.IsThereYellow = false;
h.flags.IsThereRed = false;
h.flags.VideoPlaying = false;
%Map
h.data.Map = [];
h.data.EventBuf = 0;
h.data.Stim.PreStimLength = 5;
%%%%%%%
% GUI init
%%%%%%%

% Interface Generation
h.ui.fig = figure('Name', 'DataExplorer', 'Position', [100 25 1000 550], 'CloseRequestFcn', @my_closereq);

%%% ROI 
% ROIs view and management
h.ui.ROIs = uipanel('Parent', h.ui.fig, 'Title','ROI','FontSize',12,...
             'Position',[.01 .25 .495 .74]);
h.ui.AddButton = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.01 0.925 0.1 0.075],...
    'String','+', 'Callback', @AddNewROI);
h.ui.RemButton = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.11 0.925 0.1 0.075],...
    'String','-', 'Callback', @RemoveROI);
h.ui.ROIsPanel = uipanel('Parent', h.ui.ROIs,...
    'Position',[0.01 0.09 0.21 0.83]);
h.ui.ROIsListax = axes('Parent', h.ui.ROIsPanel, 'Position',[0 0 1 1], 'XLim', [0 1], 'YLim', [0 1]);
axis(h.ui.ROIsListax,'off');
h.ui.ROIsMap = axes('Parent', h.ui.ROIs, 'Position',[0.23 0.05 0.75 0.92]);

h.ui.SaveROIpb = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.01 0.01 0.1 0.075],...
    'String','save', 'Callback', @SaveROIs);
h.ui.LoadROIpb = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.11 0.01 0.1 0.075],...
    'String', 'load', 'Callback', @LoadROIs);

set(h.ui.AddButton, 'Enable', 'off'); 
set(h.ui.RemButton, 'Enable', 'off');
set(h.ui.SaveROIpb, 'Enable', 'off');
set(h.ui.LoadROIpb, 'Enable', 'off');

%%% Events 
% Events view and management
h.ui.EventsPan = uipanel('Parent', h.ui.fig, 'Title','Events','FontSize',12,...
             'Position',[.51 .25 .485 .74]);
h.ui.EventsDispPan.Container = uipanel('Parent', h.ui.EventsPan, 'Title', 'Selection', 'FontSize', 12,...
              'Position', [0.0 0.5 1.0 0.4]);
h.ui.EventsDispPan.Ax = axes('Parent', h.ui.EventsDispPan.Container, 'Position', ...
                    [0.175 0.15 0.72 0.85]);
h.ui.EventsDispPan.Cbox = uicontrol('Style','checkbox','Parent', h.ui.EventsDispPan.Container,...
    'Units', 'normalized', 'Position', [0.01 0.4 0.1 0.2],...
    'String','1', 'Callback', @OnEditEventsClicked);
h.ui.EventsDispPan.Slider = uicontrol('Parent', h.ui.EventsDispPan.Container,...
             'Style', 'slider',  'Min', 1, 'Max', 10, 'Value', 1,...
             'Units', 'normalized', 'Position', [0.95 0.0 0.05 1], 'SliderStep', [0.01 0.3],...
             'Callback', @MoveEventsDisp); 
h.ui.ROIsSelector = uicontrol('Style','popupmenu','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.20 0.01 0.15 0.075],...
    'String', 'Empty', 'Callback', @SrcChange);
h.ui.ChannelSelector = uicontrol('Style','popupmenu','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.37 0.01 0.15 0.075],...
    'String', 'Empty', 'Callback', @SrcChange);
h.ui.EventsDispRegen = uicontrol('Style','pushbutton','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.01 0.01 0.15 0.075],...
    'String','Gen', 'Callback', @PopulateEvntsDisplay);
h.ui.FilteringOpt = uicontrol('Style','checkbox','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.575 0.01 0.25 0.075],...
    'String','Filter data', 'Callback', @FiltOptChange);
h.ui.GlobalOpt = uicontrol('Style','checkbox','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.775 0.01 0.25 0.075],...
    'String','Global signal', 'Callback', @GlobalSigChange);
h.ui.EventsMeanPan.Ax = axes('Parent', h.ui.EventsPan, 'Position', ...
                    [0.15 0.15 0.825 0.325]);
h.ui.EventsMeanPan.Tag = axes('Parent', h.ui.EventsPan, 'Position', ...
                    [0.0 0.15 0.10 0.18]); axis(h.ui.EventsMeanPan.Tag, 'off');
text(h.ui.EventsMeanPan.Tag, 0.225, 0.5, 'Mean', 'Rotation', 90); 

set(h.ui.EventsDispPan.Slider, 'Enable', 'off'); 
set(h.ui.ROIsSelector, 'Enable', 'off'); 
set(h.ui.ChannelSelector, 'Enable', 'off'); 
set(h.ui.FilteringOpt, 'Enable', 'off'); 
set(h.ui.GlobalOpt, 'Enable', 'off'); 
set(h.ui.EventsDispRegen, 'Enable', 'off'); 
%%% Videos 
% Display animated sequences 
h.ui.AnimatedPan = uipanel('Parent', h.ui.fig, 'Title','Sequence','FontSize',12,...
             'Position',[.01 .01 .125 .105]);
h.ui.StartVideo = uicontrol('Style','pushbutton','Parent', h.ui.AnimatedPan,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Start', 'Callback', @StartVideo);
h.ui.Overlay = 0;

%%% Data Loading:
h.ui.NewData = uipanel('Parent', h.ui.fig, 'Title','Load DataSet','FontSize',12,...
             'Position',[.01 .125 .125 .105]);
h.ui.LoadData = uicontrol('Style','pushbutton','Parent', h.ui.NewData,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Load', 'Callback', @OpenFolder);

%%% Data graphs:
h.ui.Graph = uipanel('Parent', h.ui.fig, 'Title','Figures','FontSize',12,...
             'Position',[.150 .125 .125 .105]);
h.ui.GGraph = uicontrol('Style','pushbutton','Parent', h.ui.Graph,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Generate', 'Callback', @GenerateGraphs);

%%% Data Export:
h.ui.XLS = uipanel('Parent', h.ui.fig, 'Title','Spreadsheet','FontSize',12,...
             'Position',[.150 .01 .125 .105]);
h.ui.Eport = uicontrol('Style','pushbutton','Parent', h.ui.XLS,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Export', 'Callback', @exportXLS);

%%% Data Export:
h.ui.PStimL = uipanel('Parent', h.ui.fig, 'Title','PreStim','FontSize',12,...
             'Position',[.290 .125 .125 .105]);
h.ui.SetPS = uicontrol('Style','popupmenu','Parent', h.ui.PStimL,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String',{'0s', '2s', '5s', '10s'}, 'Value', 1, 'Callback', @setPreStimLength);

%%% Intensity checkup:
h.ui.Icheck = uipanel('Parent', h.ui.fig, 'Title','Intensity','FontSize',12,...
             'Position',[.290 .01 .125 .105]);
h.ui.IChckButton = uicontrol('Style','pushbutton','Parent', h.ui.Icheck,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String', 'Check', 'Callback', @validateIntensity);

%%%%%%%%%%%%%%%
%Fonctions & Callbacks:
%%%%%%%%%%%%%%%
    function validateIntensity(~,~,~)
        cmap = gray(4096);
        cmap(1:512,1) = 1;
        cmap(1:512,2:3) = 0;
        cmap((end-511):end,1:2) = 0;
        cmap((end-511):end,3) = 1;
        if( h.flags.IsThereGreen )
            Datptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
            d = memmapfile(Datptr.datFile, 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Green Channel');
            colormap(cmap);
            axis image; axis off; colorbar;
        end
        if( h.flags.IsThereYellow )
            Datptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
            d = memmapfile(Datptr.datFile, 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Yellow Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end
        if( h.flags.IsThereRed )
            Datptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
            d = memmapfile(Datptr.datFile, 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Red Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end
        if( h.flags.IsThereFlow )
            d = memmapfile([h.paths.FolderName filesep 'sChan.dat'], 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Flow Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end 
        if( h.flags.IsThereFluo )
            d = memmapfile([h.paths.FolderName filesep 'fChan.dat'], 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Fluo Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end 
    end

    function setPreStimLength(Src,~,~)
        sID = get(h.ui.SetPS, 'Value');
        if( sID == 1 )
            h.data.Stim.PreStimLength = 0.11; 
        elseif( sID == 2 )
            h.data.Stim.PreStimLength = 2; 
        elseif( sID == 3 )
            h.data.Stim.PreStimLength = 5; 
        elseif( sID == 4 )
            h.data.Stim.PreStimLength = 10;
        end
        OpenFolder(Src);
    end

    function GenerateGraphs(~,~,~)
        %Waiting Dlg...                
        GraphsDlg = dialog('Position',[500 500 250 150],'Name','Graphs');
        GraphsStr = uicontrol('Parent', GraphsDlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Saving Map...');
        pause(0.1);
        Map = zeros(size(h.data.Map,1), size(h.data.Map,2), 3);
        G = Map(:,:,1);
        if( h.flags.IsThereGreen )
            G = reshape(h.data.gDatPtr.Data(1:length(h.data.Map(:))), size(h.data.Map));
        end
        Y = Map(:,:,1);
        if( h.flags.IsThereYellow )
            Y = reshape(h.data.yDatPtr.Data(1:length(h.data.Map(:))), size(h.data.Map));
        end
        R = Map(:,:,1);
        if( h.flags.IsThereRed )
            R = reshape(h.data.rDatPtr.Data(1:length(h.data.Map(:))), size(h.data.Map));
        end
        Map(:,:,1) = (84.46*R + 77.5326*Y + 27.6758*G);
        Map(:,:,2) = (25.07*R + 52.1235*Y + 143.7679*G);
        Map(:,:,3) = (17.95*R + 21.7241*Y + 66.4394*G);
        
        Map = bsxfun(@rdivide, Map, permute(max(reshape(Map,[],3), [], 1), [3 1 2]));
        Map = bsxfun(@rdivide, Map, permute([1 1.275 1.35], [3 1 2]));
        fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
        imshow(Map);
        axis image; axis off;
        title('Imaged Area');
        figName = 'RawMap';
        saveas(fig, [h.paths.Graphs figName], 'png');
        close(fig);
        
       
        eLen = floor(h.data.AcqFreq*...
                (h.data.Stim.StimLength + h.data.Stim.InterStim_min));
        T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.StimLength + h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
        %for each ROI
        Map = zeros(size(h.data.Map));
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            Map = Map + mask;
            set(GraphsStr, 'String', ['ROI #' int2str(indR) ' Map saving ...']);
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            image(Map); hold 'on';
            image(~repmat(mask,1,1,3), 'AlphaData', ...
                mask.*0.10 + (~imerode(mask,strel('diamond',1))&mask)*0.9);
            title([h.data.ROIs{indR}.name ' Map']);
            axis image; axis off;
            figName = [h.data.ROIs{indR}.name ' Map'];
            saveas(fig, [h.paths.Graphs figName], 'png');
            close(fig);
            
            AccumGreen = zeros(eLen, length(h.data.EvntList),'single');
            AccumYellow = zeros(eLen, length(h.data.EvntList),'single');
            AccumRed = zeros(eLen, length(h.data.EvntList),'single');
            AccumHbO = zeros(eLen, length(h.data.EvntList),'single');
            AccumHbR = zeros(eLen, length(h.data.EvntList),'single');
            AccumFlow = zeros(eLen, length(h.data.EvntList),'single');
            for indE = 1:length(h.data.EvntList)
                set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Colours...']);
                
                fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                ax = axes('Parent', fig);
                hold(ax,'on');
                maxi = 1.01;
                mini = 0.99;

                if( h.flags.IsThereGreen )
                    %Open
                    Datptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
                    d = memmapfile(Datptr.datFile, 'Format', 'single');
                    d = d.Data((length(h.data.Map(:))*(h.data.G_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.G_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                     Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     d = d./L;
                
                     if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     %Plot
                     plot(ax, T, d, 'Color', [0.0 0.75 0.0], 'LineWidth', 2);
                     AccumGreen(:,indE) = d;
                end
                if( h.flags.IsThereYellow )
                    %Open
                    Datptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
                    d = memmapfile(Datptr.datFile, 'Format', 'single');
                    d = d.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                     Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     d = d./L;
                
                      if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     
                     %Plot
                     plot(ax, T, d, 'Color', [0.75 0.75 0.0], 'LineWidth', 2);
                     AccumYellow(:,indE) = d;
                end
                if( h.flags.IsThereRed )
                    %Open
                    Datptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
                    d = memmapfile(Datptr.datFile, 'Format', 'single');
                    d = d.Data((length(h.data.Map(:))*(h.data.R_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.R_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                     Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     d = d./L;
                
                      if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     
                     %Plot
                     plot(ax, T, d, 'Color', [0.75 0.0 0.0], 'LineWidth', 2);
                     AccumRed(:,indE) = d;
                end
                clear d;
                box(ax,'on');
                set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                title(['{\Delta}Reflectance over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                ylabel('Normalized Reflectance');
                xlabel('Time (sec)');
                xlim([T(1), T(end)]);
                ylim([mini, maxi]);
                line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                line(ax, [0 0], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
                line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
                %Save figure for colours:
                figName = ['Colours_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                saveas(fig, [h.paths.Graphs figName], 'png');
                close(fig);
                
                %Figure for Flow:
                if( h.flags.IsThereFlow )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Flow...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 1.25;
                    mini = 0.75;
                    %Open
                    dF = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
                    dF = dF.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dF(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dF((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     dF = dF./L;
                                    
                     if( max(dF) > maxi) 
                         maxi = max(dF);
                     end
                     if( min(dF) < mini) 
                         mini = min(dF);
                     end
                     
                     %Plot
                     plot(ax, T, dF, 'Color', [0.0 0.0 0.0], 'LineWidth', 2);
                     box(ax,'on');
                     set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                     title(['{\Delta}Flow over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                     ylabel('{\Delta}Flow');
                     xlabel('Time (sec)');
                     
                     line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                     line(ax, [0 0], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Flow_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                     AccumFlow(:,indE) = dF;
                     clear dF;
                end
                
                if( h.flags.IsThereFluo )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Fluo...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 1.25;
                    mini = 0.75;
                    %Open
                    dF = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
                    dF = dF.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dF(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dF((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     dF = dF./L;
                                    
                     if( max(dF) > maxi) 
                         maxi = max(dF);
                     end
                     if( min(dF) < mini) 
                         mini = min(dF);
                     end
                     
                     %Plot
                     plot(ax, T, dF, 'Color', [0.0 0.0 0.0], 'LineWidth', 2);
                     box(ax,'on');
                     set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                     title(['{\Delta}Flow over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                     ylabel('{\Delta}Flow');
                     xlabel('Time (sec)');
                     
                     line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                     line(ax, [0 0], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Flow_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                     AccumFlow(:,indE) = dF;
                     clear dF;
                end
                
                %Figure for Hbs:
                if( h.flags.IsThereHbO && h.flags.IsThereHbR )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ' Hbs...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 5;
                    mini = -5;
                    %Open
                    dO = memmapfile(h.data.HBinfos.datFileHbO, 'Format', 'single');
                    dO = dO.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dR = memmapfile(h.data.HBinfos.datFileHbR, 'Format', 'single');
                    dR = dR.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dO = reshape(dO, [], eLen);
                    dO = mean(dO(mask(:) == 1, :), 1);
                    dR = reshape(dR, [], eLen);
                    dR = mean(dR(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dO(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dO((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     dO = dO - L;
                     Pstart = median(dR(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dR((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     dR = dR - L;
                
                     if( max(dO) > maxi) 
                         maxi = max(dO);
                     end
                     if( min(dO) < mini) 
                         mini = min(dO);
                     end
                     if( max(dR) > maxi)
                         maxi = max(dR);
                     end
                     if( min(dR) < mini)
                         mini = min(dR);
                     end
                     
                     %Plot
                     plot(ax, T, dR, 'Color', [0.0 0.0 1.0], 'LineWidth', 2);
                     plot(ax, T, dO, 'Color', [1.0 0.0 0.0], 'LineWidth', 2);
                     plot(ax, T, dR + dO, 'Color', [0.0 1.0 0.0], 'LineWidth', 2);
                     box(ax,'on');
                     set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                     title(['{\Delta}Hb over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                     ylabel('{\Delta}Hb Concentration');
                     xlabel('Time (sec)');
                     
                     line(ax, [T(1) T(end)], [0 0], 'Color', 'k', 'LineStyle',':');
                     line(ax, [0 0], [-5 5], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [-5 5], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Hb_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                     AccumHbO(:,indE) = dO;
                     AccumHbR(:,indE) = dR;
                     clear dO dR;
                end
                
            end
            set(GraphsStr, 'String', ['ROI #' int2str(indR) ' average graphs...']);
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            hold(ax,'on');
            if(h.flags.IsThereGreen)
                errorbar(ax, T, mean(AccumGreen,2), std(AccumGreen,1,2)/sqrt(indE),...
                    'Color', [0.0 0.75 0.0], 'LineWidth', 2);
            end
            if(h.flags.IsThereYellow)
                errorbar(ax, T, mean(AccumYellow,2), std(AccumYellow,1,2)/sqrt(indE),...
                    'Color', [0.75 0.75 0.], 'LineWidth', 2);
            end
            if(h.flags.IsThereRed)
                errorbar(ax, T, mean(AccumRed,2), std(AccumRed,1,2)/sqrt(indE),...
                    'Color', [0.75 0.0 0.0], 'LineWidth', 2);
            end
            box(ax,'on');
            set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
            title(['Mean reflectance over ' h.data.ROIs{indR}.name]);
            ylabel('Reflectance');
            xlabel('Time (sec)');
            line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
            line(ax, [0 0], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
            line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
            axis tight;
            %Save figure for colours:
            figName = ['MeanColours_' h.data.ROIs{indR}.name];
            saveas(fig, [h.paths.Graphs figName], 'png');
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            hold(ax,'on');
            if(h.flags.IsThereHbO)
                errorbar(ax, T, mean(AccumHbO,2), std(AccumHbO,1,2)/sqrt(indE),...
                    'Color', [1. 0.0 0.0], 'LineWidth', 2);
                errorbar(ax, T, mean(AccumHbR,2), std(AccumHbR,1,2)/sqrt(indE),...
                    'Color', [0.0 0.0 1.], 'LineWidth', 2);
                errorbar(ax, T, mean(AccumHbR+AccumHbO,2), std(AccumHbR+AccumHbO,1,2)/sqrt(indE),...
                    'Color', [0.0 1. 0.0], 'LineWidth', 2);
            end
            box(ax,'on');
            set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
            title(['{\Delta}Hb over ' h.data.ROIs{indR}.name]);
            ylabel('{\Delta}Hb Concentration');
            xlabel('Time (sec)');
            line(ax, [T(1) T(end)], [0 0], 'Color', 'k', 'LineStyle',':');
            line(ax, [0 0], [-5 5], 'Color', 'k', 'LineStyle','--');
            line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [-5 5], 'Color', 'k', 'LineStyle','--');
            axis tight;
            %Save figure for colours:
            figName = ['MeanHbs_' h.data.ROIs{indR}.name];
            saveas(fig, [h.paths.Graphs figName], 'png');
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            hold(ax,'on');
            if( h.flags.IsThereFlow )
                errorbar(ax, T, mean(AccumFlow,2), std(AccumFlow,1,2)/sqrt(indE),...
                    'Color', 'k', 'LineWidth', 2);
                title(['Mean {\Delta}Flow over ' h.data.ROIs{indR}.name]);
                ylabel('{\Delta}Flow');
                xlabel('Time (sec)');
                
                line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                line(ax, [0 0], [0.9 1.1], 'Color', 'k', 'LineStyle','--');
                line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [0.9 1.1], 'Color', 'k', 'LineStyle','--');
                axis tight;
                %Save figure for colours:
                figName = ['MeanFlow_' h.data.ROIs{indR}.name];
                saveas(fig, [h.paths.Graphs figName], 'png');
                close(fig);
            end
        end
        set(GraphsStr, 'String', 'Video Sequences...');
        
        Map(Map > 1) = 1;
        Map = Map ==1;
        %Video sequence of each Colour
        if( h.flags.IsThereGreen )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                Datptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
                dat = memmapfile(Datptr.datFile, 'Format', 'single');
                dat = dat.Data((length(h.data.Map(:))*(h.data.G_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.G_eflag(indE) +eLen - 1)) );
                dat = reshape(dat, size(Map,1), size(Map,2), []);
                Accum = Accum + dat;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = bsxfun(@times, Accum, Map);
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'AbsGreen.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Green Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = bsxfun(@rdivide, Accum, mean(Accum,3));
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelGreen.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Green Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
        end
        
        if( h.flags.IsThereYellow )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                Datptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
                dat = memmapfile(Datptr.datFile, 'Format', 'single');
                dat = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                dat = reshape(dat, size(Map,1), size(Map,2), []);
                Accum = Accum + dat;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = bsxfun(@times, Accum, Map);
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'AbsYellow.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Yellow Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = bsxfun(@rdivide, Accum, mean(Accum,3));
            Accum = reshape(Accum,[],size(Accum,3));
            maxi = max(max(Accum(Map(:),:),[],1),[],2);
            mini = min(min(Accum(Map(:),:),[],1),[],2);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelYellow.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [mini maxi]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Yellow Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
        end
        
        if( h.flags.IsThereRed )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                Datptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
                dat = memmapfile(Datptr.datFile, 'Format', 'single');
                dat = dat.Data((length(h.data.Map(:))*(h.data.R_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.R_eflag(indE) +eLen - 1)) );
                dat = reshape(dat, size(Map,1), size(Map,2), []);
                Accum = Accum + dat;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = bsxfun(@times, Accum, Map);
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'AbsRed.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Red Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = bsxfun(@rdivide, Accum, mean(Accum,3));
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelRed.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Red Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
        end
        
        %Video sequence of HbO, HbR & HbT
        if( h.flags.IsThereHbO )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            AccumO = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                Datptr = matfile([h.paths.FolderName filesep 'Data_Hbs.mat']);
                dat = memmapfile(Datptr.datFileHbO, 'Format', 'single');
                dat = dat.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                dat = reshape(dat, size(Map,1), size(Map,2), []);
                AccumO = AccumO + dat;
            end
            AccumO = AccumO./sum(h.data.EvntList);
            AccumO = bsxfun(@times, AccumO, Map);
            AccumO = reshape(AccumO,[],size(AccumO,3));
            P = prctile(reshape(AccumO(Map(:),:),[],1),[5 95]);
            AccumO = reshape(AccumO, size(Map,1),[],size(AccumO,2));
            v = VideoWriter([h.paths.Graphs filesep 'HbO.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(AccumO,3)
                imagesc(ax, squeeze(AccumO(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['HbO Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            AccumR = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                Datptr = matfile([h.paths.FolderName filesep 'Data_Hbs.mat']);
                dat = memmapfile(Datptr.datFileHbR, 'Format', 'single');
                dat = dat.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                dat = reshape(dat, size(Map,1), size(Map,2), []);
                AccumR = AccumR + dat;
            end
            AccumR = AccumR./sum(h.data.EvntList);
            AccumR = bsxfun(@times, AccumR, Map);
            AccumR = reshape(AccumR,[],size(AccumR,3));
            P = prctile(reshape(AccumR(Map(:),:),[],1),[5 95]);
            AccumR = reshape(AccumR, size(Map,1),[],size(AccumR,2));
            v = VideoWriter([h.paths.Graphs filesep 'HbR.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(AccumR,3)
                imagesc(ax, squeeze(AccumR(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['HbO Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = AccumO + AccumR;
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'HbT.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['HbT Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
        end
        
        %Video sequence of flow
        if( h.flags.IsThereFlow )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                Datptr = matfile([h.paths.FolderName filesep 'Flow_infos.mat']);
                dat = memmapfile(Datptr.datFile, 'Format', 'single');
                dat = dat.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                dat = reshape(dat, size(Map,1), size(Map,2), []);
                Accum = Accum + dat;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = bsxfun(@times, Accum, Map);
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'Flow.avi']);
            v.FrameRate = 1.25;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap jet;
                colorbar;
                title(['Flow Intensity at: ' num2str(T(indF)) ' sec']);
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
        end
        delete(GraphsDlg);
    end

    function exportXLS(~,~,~)
        %Waiting Dlg...                
        ExportDlg = dialog('Position',[500 500 250 150],'Name','Export');
        uicontrol('Parent', ExportDlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Exporting data...');
        pause(0.1);
        
        eLen = floor(h.data.AcqFreq*...
                (h.data.Stim.StimLength + h.data.Stim.InterStim_min));
        T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.StimLength + h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
        array = zeros(eLen, size(h.data.ROIs,2)*length(h.data.EvntList)+1, 'single');
        filename = [h.paths.FolderName filesep 'DataExport.xls'];
        array(:, 1) = T; names = {'T'};
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            for indE = 1:length(h.data.EvntList)
                dO =   h.data.hoDatPtr.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                dO = reshape(dO, [], eLen);
                dO = mean(dO(mask(:) == 1, :), 1);
                array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dO;
                names = cat(1, names, {['R' int2str(indR) 'E' int2str(indE)]});
            end
        end
        Tabl = array2table(array, 'VariableNames', names);
        writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','HbO','Range','A1');
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            for indE = 1:length(h.data.EvntList)
                dR =   h.data.hrDatPtr.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                dR = reshape(dR, [], eLen);
                dR = mean(dR(mask(:) == 1, :), 1);
               array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dR;
            end
        end
        Tabl = array2table(array, 'VariableNames', names);
        writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','HbR','Range','A1');
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            for indE = 1:length(h.data.EvntList)
                dF =   h.data.fDatPtr.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                dF = reshape(dF, [], eLen);
                dF = mean(dF(mask(:) == 1, :), 1);
               array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
            end
        end
        Tabl = array2table(array, 'VariableNames', names);
        writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Flow','Range','A1');
        
        %Green
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            for indE = 1:length(h.data.EvntList)
                dF =   h.data.gDatPtr.Data((length(h.data.Map(:))*(h.data.G_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.G_eflag(indE) +eLen - 1)) );
                dF = reshape(dF, [], eLen);
                dF = mean(dF(mask(:) == 1, :), 1);
               array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
            end
        end
        Tabl = array2table(array, 'VariableNames', names);
        writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Green','Range','A1');
        %Yellow
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            for indE = 1:length(h.data.EvntList)
                dF =   h.data.yDatPtr.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                dF = reshape(dF, [], eLen);
                dF = mean(dF(mask(:) == 1, :), 1);
               array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
            end
        end
        Tabl = array2table(array, 'VariableNames', names);
        writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Yellow','Range','A1');
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            for indE = 1:length(h.data.EvntList)
                dF =   h.data.rDatPtr.Data((length(h.data.Map(:))*(h.data.R_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.R_eflag(indE) +eLen - 1)) );
                dF = reshape(dF, [], eLen);
                dF = mean(dF(mask(:) == 1, :), 1);
               array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
            end
        end
        Tabl = array2table(array, 'VariableNames', names);
        writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Red','Range','A1');
        
        delete(ExportDlg);
    end

    function OpenFolder(Src, ~, ~)
        if( strcmp(Src.String, 'Load') )
            h.paths.FolderName = uigetdir();
        elseif( strcmp(h.paths.FolderName, '') )
            return;
        end
        
        h.flags.bsaveROIS = false; %flag to know if any changes were made to ROIs
        h.flags.bsaveEvnts = false;
        h.flags.Stim = false;
        h.flags.IsThereHbO = false;
        h.flags.IsThereHbR = false;
        h.flags.IsThereHbT = false;
        h.flags.IsThereFlow = false;
        h.flags.IsThereGreen = false;
        h.flags.IsThereRed = false;
        h.flags.IsThereYellow = false;
               
        % Files Path
        h.paths.HbFile = [h.paths.FolderName filesep 'Data_Hbs.mat'];
        h.paths.ROIsFile = [h.paths.FolderName filesep 'ROIs.mat'];
        h.paths.EVNTsFile = [h.paths.FolderName filesep 'Events.mat'];
        h.paths.StimProto = [h.paths.FolderName filesep 'StimParameters.mat'];
        h.paths.Graphs = [h.paths.FolderName filesep 'Graphs' filesep];
        h.paths.Flow = [h.paths.FolderName filesep 'Flow_infos.mat'];
        h.paths.Fluo = [h.paths.FolderName filesep 'Data_Fluo.mat'];
        
        if( strcmp(Src.String, 'Load') )
            if( exist(h.paths.Graphs , 'dir') )
                ButtonName = questdlg('A folder containing figures already exist. Do you want to overwrite it?', ...
                    'Figures folder', ...
                    'Yes', 'Change', 'Yes');
                switch ButtonName
                    case 'Yes'
                        disp('Erasing old figures...');
                        delete([h.paths.Graphs '*.*']);
                    case 'Change'
                        dname = uigetdir(h.paths.FolderName);
                        h.paths.Graphs = dname;
                end % switch
                
            else
                mkdir(h.paths.Graphs);
            end
        end
        
        %Waiting Dlg...                
        opendlg = dialog('Position',[500 500 250 150],'Name','Loading...');
        uicontrol('Parent', opendlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Loading data. Please wait...');
        pause(0.1);
        
        %Load stimulation parameters
        if(  exist(h.paths.StimProto, 'file') )
            h.data.Stim = load(h.paths.StimProto);
            sID = get(h.ui.SetPS, 'Value');
            if( sID == 1 )
                h.data.Stim.PreStimLength = 0.11;
            elseif( sID == 2 )
                h.data.Stim.PreStimLength = 2;
            elseif( sID == 3 )
                h.data.Stim.PreStimLength = 5;
            elseif( sID == 4 )
                h.data.Stim.PreStimLength = 10;
            end
            h.flags.Stim = true;
        else
            h.flags.Stim = false;
            h.data.Stim.NbStim = 1;
            h.data.Stim.PreStimLength = 0;
            h.data.Stim.StimLength = 0;
            h.data.Stim.InterStim_min = 0;
        end
        
        h.data.AcqFreq = 0;
        Ts = 1e6;
                
        if( exist(h.paths.HbFile, 'file') )
            h.flags.IsThereHbO = true;
            h.flags.IsThereHbR = true;
            h.flags.IsThereHbT = true;
            
            h.data.HBinfos = matfile(h.paths.HbFile,'Writable', true);
            h.data.AcqFreq = h.data.HBinfos.Freq;
            if( ~exist(h.data.HBinfos.datFileHbO,'file') )
                h.data.HBinfos.datFileHbO = [h.paths.FolderName filesep 'HbO.dat'];
                h.data.HBinfos.datFileHbR = [h.paths.FolderName filesep 'HbR.dat'];
            end
            h.data.hoDatPtr = memmapfile(h.data.HBinfos.datFileHbO, 'Format', 'single');
            h.data.hrDatPtr = memmapfile(h.data.HBinfos.datFileHbR, 'Format', 'single');
            nframes = h.data.HBinfos.datLength;
            Ts = min(nframes, Ts);
             
            stim = h.data.HBinfos.Stim;
            if( ~h.flags.Stim )
                h.data.Stim.StimLength = length(stim)/h.data.HBinfos.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.Stim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
             if( isempty(Start) )
                 h.data.H_eflag = 1;
             else
                h.data.H_eflag = Start;
             end
        else            
            disp('No Hb concentrations were computed for this experiment!');
            h.flags.IsThereHbO = false;
            h.flags.IsThereHbR = false;
            h.flags.IsThereHbT = false;
        end
  
        Map = [];
        RawDatFiles = dir([h.paths.FolderName filesep 'Data_*.mat']);
        if( isempty(RawDatFiles) )
            %No compatible files were found
            h.flags.IsThereHbO = false;
            h.flags.IsThereHbR = false;
            h.flags.IsThereHbT = false;
            h.flags.Stim = false;
            
            disp(['No data files found in ' FolderName ' Folder.']);
            disp('There is nothing to show you... sorry!');
            return;
        end
        %Green channel:
        if( ~isempty(strfind([RawDatFiles.name],'green')) ) %#ok<*STREMP>
            h.flags.IsThereGreen = true;
            Dat_Gptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
            nrows = Dat_Gptr.datSize(1,1);
            ncols = Dat_Gptr.datSize(1,2);
            nframes = Dat_Gptr.datLength;
            Freq =  Dat_Gptr.Freq;
            h.data.AcqFreq = Freq;
            Ws = ncols;
            Hs = nrows;
            Ts = min(nframes, Ts);
            stim = Dat_Gptr.Stim;
            if( ~h.flags.Stim )
                h.data.Stim.StimLength = length(stim)/Dat_Gptr.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.Stim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            
            if( isempty(Start) )
                 h.data.G_eflag = 1;
             else
                h.data.G_eflag = Start;
            end
             if( ~exist(Dat_Gptr.datFile,'file') )
                Dat_Gptr.datFile = [h.paths.FolderName filesep 'gChan.dat'];
            end
            h.data.gDatPtr = memmapfile(Dat_Gptr.datFile,...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.gDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        
        %Yellow channel:
        if( ~isempty(strfind([RawDatFiles.name],'yellow')) )
            h.flags.IsThereYellow = true;
            Dat_Yptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
            nrows = Dat_Yptr.datSize(1,1);
            ncols = Dat_Yptr.datSize(1,2);
            nframes = Dat_Yptr.datLength;
            Freq =  Dat_Yptr.Freq;
            h.data.AcqFreq = Freq;
            Ws = ncols;
            Hs = nrows;
            Ts = min(Ts, nframes);
            stim = Dat_Yptr.Stim;
            if( ~h.flags.Stim )
                h.data.Stim.StimLength = length(stim)/Dat_Yptr.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.Stim.PreStimLength);
            if( idxS < 1 )
                idxS = 1;
            end
            Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            if( ~exist(Dat_Yptr.datFile,'file') )
                Dat_Yptr.datFile = [h.paths.FolderName filesep 'gChan.dat'];
            end
            h.data.yDatPtr = memmapfile(Dat_Yptr.datFile,...
                'Format', 'single');
            
              if( isempty(Start) )
                 h.data.Y_eflag = 1;
             else
                h.data.Y_eflag = Start;
             end
            if( isempty(Map) )
                Map = reshape(h.data.yDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        
        %Red channel:
        if( ~isempty(strfind([RawDatFiles.name],'red')) )
            h.flags.IsThereRed = true;
            Dat_Rptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
            nrows = Dat_Rptr.datSize(1,1);
            ncols = Dat_Rptr.datSize(1,2);
            nframes = Dat_Rptr.datLength;
            Freq =  Dat_Rptr.Freq;
            h.data.AcqFreq = Freq;
            Ws = ncols;
            Hs = nrows;
            Ts = min(Ts, nframes);
            stim = Dat_Rptr.Stim;
            if( ~h.flags.Stim )
                h.data.Stim.StimLength = length(stim)/Dat_Rptr.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.Stim.PreStimLength);
            if( idxS < 1 )
                idxS = 1;
            end
            Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
           
               if( isempty(Start) )
                 h.data.R_eflag = 1;
             else
                h.data.R_eflag = Start;
             end
            if( ~exist(Dat_Rptr.datFile,'file') )
                Dat_Rptr.datFile = [h.paths.FolderName filesep 'gChan.dat'];
            end
            h.data.rDatPtr = memmapfile(Dat_Rptr.datFile,...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.rDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        if( exist(h.paths.Flow, 'file') )
             h.flags.IsThereFlow = true;
             
             h.data.fInfo = matfile(h.paths.Flow,'Writable',true);
             h.data.AcqFreq = h.data.fInfo.Freq;
             nrows = h.data.fInfo.datSize(1,1);
             ncols = h.data.fInfo.datSize(1,2);
             Ws = ncols;
             Hs = nrows;
             if( ~exist(h.data.fInfo.datFile,'file') )
                h.data.fInfo.datFile = [h.paths.FolderName filesep 'sChan.dat'];
             end
             if( exist([h.paths.FolderName filesep 'Flow.dat'], 'file') )
                h.data.fInfo.datFile = [h.paths.FolderName filesep 'Flow.dat'];
             end
             h.data.fDatPtr = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
             nframes = h.data.fInfo.datLength;
             Ts = min(nframes, Ts);
             
             stim = h.data.fInfo.Stim;
             if( size(stim,2) > size(stim,1) )
                 stim = stim';
             end
             idxS = floor(h.data.AcqFreq*h.data.Stim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            
             if( isempty(Start) )
                 h.data.F_eflag = 1;
             else
                h.data.F_eflag = Start;
             end
             
             if( isempty(Map) )
                Map = reshape(h.data.fDatPtr.Data(1:(ncols*nrows)),nrows,[]);
             end
        else
             disp('No flow measures for this experiment!');
             h.flags.IsThereFlow = false;
        end
        if( exist(h.paths.Fluo, 'file') )
             h.flags.IsThereFluo = true;
             
             h.data.fInfo = matfile(h.paths.Fluo,'Writable', true);
             h.data.AcqFreq = h.data.fInfo.Freq;
             nrows = h.data.fInfo.datSize(1,1);
             ncols = h.data.fInfo.datSize(1,2);
             Ws = ncols;
             Hs = nrows;
             if( ~exist(h.data.fInfo.datFile,'file') )
                h.data.fInfo.datFile = [h.paths.FolderName filesep 'fChan.dat'];
             end
             h.data.fDatPtr = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
             nframes = h.data.fInfo.datLength;
             Ts = min(nframes, Ts);
             
             stim = h.data.fInfo.Stim;
             if( ~h.flags.Stim )
                 h.data.Stim.StimLength = length(stim)/h.data.fInfo.Freq;
             end
             if( size(stim,2) > size(stim,1) )
                 stim = stim';
             end
            idxS = floor(h.data.AcqFreq*h.data.Stim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            
             if( isempty(Start) )
                 h.data.F_eflag = 1;
             else
                h.data.F_eflag = Start;
             end
             
             if( isempty(Map) )
                Map = reshape(h.data.fDatPtr.Data(1:(ncols*nrows)),nrows,[]);
             end
        else
             disp('No fluorescence measures for this experiment!');
             h.flags.IsThereFluo = false;
        end
        h.data.NCols = Ws;
        h.data.NRows = Hs;
        h.data.NFrames = Ts;
        
        if( round(h.data.AcqFreq*h.data.Stim.StimLength) > Ts )
            h.data.Stim.StimLength = Ts/h.data.AcqFreq;
        end
        
        %Map
        h.data.Map =  double(Map)./max(double(Map(:)));
        imshow(Map,[],'Parent',h.ui.ROIsMap);
        axis(h.ui.ROIsMap,'off', 'image');
        
        %ROIs file:
        ROIs = {};
        if( exist(h.paths.ROIsFile, 'file') )
            load(h.paths.ROIsFile);
            h.data.ROIs = ROIs;
            clear ROIs;
        else
            h.data.ROIs = ROIs;
        end
        
        %Events file:
        E = ones(1, h.data.Stim.NbStim);
        if(  exist(h.paths.EVNTsFile, 'file')  )
            load([h.paths.FolderName filesep 'Events.mat']);
        end
        h.data.EvntList = E;
        h.ui.EventsDispPan.Slider.Max = h.data.Stim.NbStim;
        clear E;
        
        Str = {};
        if( h.flags.IsThereGreen )
            Str{end+1} = 'Green';
        end
        if( h.flags.IsThereRed )
            Str{end+1} = 'Red';          
        end
        if( h.flags.IsThereYellow )    
            Str{end+1} = 'Yellow'; 
        end
        if( h.flags.IsThereHbO )    
            Str{end+1} = 'HbO';        
        end
        if( h.flags.IsThereHbR )    
            Str{end+1} = 'HbR';         
        end
        if( h.flags.IsThereHbT )
            Str{end+1} = 'HbT';         
        end
        if( h.flags.IsThereFlow )
            Str{end+1} = 'Flow';         
        end
        if( h.flags.IsThereFluo )
            Str{end+1} = 'Fluo';         
        end
        set(h.ui.ChannelSelector,'String', Str);
        
        set(h.ui.AddButton, 'Enable', 'on');
        set(h.ui.RemButton, 'Enable', 'on');
        set(h.ui.SaveROIpb, 'Enable', 'on');
        set(h.ui.LoadROIpb, 'Enable', 'on');
        set(h.ui.EventsDispPan.Slider, 'Enable', 'on');
        set(h.ui.ROIsSelector, 'Enable', 'on');
        set(h.ui.ChannelSelector, 'Enable', 'on');
        set(h.ui.FilteringOpt, 'Enable', 'on');
        set(h.ui.GlobalOpt, 'Enable', 'on');
        set(h.ui.StartVideo, 'Enable', 'on');

        RefreshLoop('All');
        
        delete(opendlg);
    end

    function StartVideo(~,~,~)
        
        if( ~h.flags.VideoPlaying )
            h.ui.VScreen = figure('Name', 'Video', 'Position', [200 200 500 500], 'Visible', 'off');
            colormap(h.ui.VScreen, 'jet');
            h.ui.ScreenAx = axes('Parent', h.ui.VScreen);
            
            Str = {};
            if( h.flags.IsThereGreen )    
                Str{end+1} = 'Green';   
            end
            if( h.flags.IsThereRed )    
                Str{end+1} = 'Red';          
            end            
            if( h.flags.IsThereYellow )
                Str{end+1} = 'Yellow'; 
            end
            if( h.flags.IsThereHbO )    
                Str{end+1} = 'HbO';        
            end
            if( h.flags.IsThereHbR )    
                Str{end+1} = 'HbR';         
            end
            if( h.flags.IsThereHbR && h.flags.IsThereHbO)  
                Str{end+1} = 'HbT'; 
            end
            if( h.flags.IsThereFlow )  
                Str{end+1} = 'Flow'; 
            end
            if( h.flags.IsThereFluo )
                Str{end+1} = 'Fluo';         
            end
            [selchan, valid] = listdlg('PromptString', 'Select channel:',...
                'SelectionMode', 'single',...
                'ListString', Str);
            if( ~valid )
                return;
            end
            SelectedSrc = Str{selchan};
            
            if( size(h.data.ROIs,2) > 0 )
                Str = cell(size(h.data.ROIs,2)+1,1);
                Str{1} = 'All';
                for indR = 1:size(h.data.ROIs,2)
                    Str{indR+1} = h.data.ROIs{indR}.name;
                end
            else
                Str = 'All';
            end
            [selroi, valid] = listdlg('PromptString', 'Select Roi:',...
                'SelectionMode', 'single',...
                'ListString', Str);
            if( ~valid )
                return;
            end
                      
            %Waiting Dlg...
            opendlg = dialog('Position',[500 500 250 150],'Name','Loading...');
            uicontrol('Parent', opendlg, 'Style','text',...
                'Position',[20 80 210 40], 'String', 'Loading data. Please wait...');
            pause(0.1);
            
            data = 0;
            isZeroCentered = 0;
            isHbT = 0;
            if( ~any(contains(SelectedSrc, {'HbT', 'Flow', 'Fluo'})) )
                isHbT = 0;
                isZeroCentered = 0;
                eval(['StartPts = h.data.' SelectedSrc(1) '_eflag;']);
                eval(['data = h.data.' lower(SelectedSrc(1)) 'DatPtr;']);
                h.data.vidClim = [0.99 1.01];
            elseif( contains(SelectedSrc, 'Flow') )
                isHbT = 0;
                isZeroCentered = 0;
                StartPts = h.data.F_eflag;
                data = h.data.fDatPtr;
                h.data.vidClim = [0 4e5];
            elseif( contains(SelectedSrc,'Fluo') )
                isHbT = 0;
                isZeroCentered = 1;
                StartPts = h.data.F_eflag;
                data = h.data.fDatPtr;
                h.data.vidClim = [-0.05 0.05];
            elseif( contains(SelectedSrc,'HbT') )
                isHbT = 1;
                isZeroCentered = 1;
                StartPts = h.data.H_eflag;
                dO = h.data.hoDatPtr;
                dR = h.data.hoDatPtr;
                h.data.vidClim = [-5 5];
            else
                isHbT = 0;
                isZeroCentered = 1;
                StartPts = h.data.H_eflag;
                eval(['data = h.data.h' lower(SelectedSrc(3)) 'DatPtr;']);
                h.data.vidClim = [-5 5];
            end
            prompt = {'Colormap minimum:','Colormap maximum:'};
            dlg_title = 'Colormap Limits';
            num_lines = 1;
            defaultans = {num2str(h.data.vidClim(1)),num2str(h.data.vidClim(2))};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            h.data.vidClim(1) = str2double(answer{1});
            h.data.vidClim(2) = str2double(answer{2});
            
            if(selroi == 1)
                h.data.vidMap = ones(size(h.data.Map));
            else
                h.data.vidMap =  h.data.ROIs{selroi-1}.mask;
            end
            if( h.data.Stim.NbStim > 1 )
                eLen = floor(h.data.AcqFreq*...
                    (h.data.Stim.StimLength + h.data.Stim.InterStim_min));
                Accum = zeros(size(h.data.Map,1), size(h.data.Map,2), eLen);
                T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.StimLength + h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
                for indE = 1:h.data.Stim.NbStim
                    if( ~isHbT )
                        d = data.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    else
                        d = dO.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                        d = d + dR.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    end
                    d = reshape(d, size(Accum));
                    
                    Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                    Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                    m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                    L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                        (m*T(round(end - h.data.Stim.PreStimLength/2))));
                    if( ~isZeroCentered )
                        d = d./L;
                    else
                        d = d - L;
                    end
                    
                    Accum = Accum + d;
                end
                Accum = Accum./indE;
            else
                eLen = h.data.NFrames;
                T = linspace(0, eLen/h.data.AcqFreq, eLen);
                d = reshape(data.Data, h.data.NRows, h.data.NCols, []);
                Accum = d;
            end
           
            
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            h.data.Accumulator = Accum;
            h.data.vidInd = 1;
            h.data.vidLimit = eLen;
            h.data.vidTimeVect = T;
            h.data.vidChan = SelectedSrc;
            set(h.ui.StartVideo, 'String', 'Stop');
            set(h.ui.VScreen, 'Visible', 'on');
            %timer...
            h.data.Timer = timer;
            h.data.Timer.TimerFcn = @ShowNextFrame;
            h.data.Timer.Period = 0.2;
            h.data.Timer.ExecutionMode = 'fixedRate';
            start(h.data.Timer);
            h.flags.VideoPlaying = true;
            delete(opendlg);
        else
            if( isvalid( h.data.Timer ) )
                stop(h.data.Timer);
                delete(h.data.Timer);
            end
            set(h.ui.StartVideo, 'String', 'Start');
            if( isvalid(h.ui.VScreen) )
                delete(h.ui.VScreen);
            end
            h.flags.VideoPlaying = false;
        end
    end

    function ShowNextFrame(~,~,~)
        
        if( ~isvalid(h.ui.VScreen) )
            stop(h.data.Timer);
            delete(h.data.Timer);
            set(h.ui.StartVideo, 'String', 'Start');
            h.flags.VideoPlaying = false;
            return;
        end
        figure(h.ui.VScreen);
        if( h.ui.Overlay )
            image(h.ui.ScreenAx, repmat(h.data.Map,1,1,3));
            hold on;
            i = imagesc(h.ui.ScreenAx, squeeze(h.data.Accumulator(:, :, h.data.vidInd)));
            hold off;
            title([h.data.vidChan ' channel at:' num2str(h.data.vidTimeVect(h.data.vidInd)) 's']);
            colorbar;
            
            Center = mean(h.data.vidClim);
            tMax = h.data.vidClim(2) - Center;
            tMin = h.data.vidClim(1) - Center;
            aMap = squeeze(h.data.Accumulator(:, :, h.data.vidInd));
            %  aMap = imfilter(aMap, fspecial('gaussian',32,8),'same','symmetric');
            aMap((aMap > Center + 0.25*tMin) & (aMap < Center + 0.25*tMax)) = 0;
            aMap(aMap > 0 ) = 1;
            
            
            alpha(i, aMap.*h.data.vidMap);
            set(h.ui.ScreenAx, 'CLim', h.data.vidClim);
            axis(h.ui.ScreenAx, 'image');
            h.data.vidInd = h.data.vidInd + 1;
            if( h.data.vidInd > h.data.vidLimit )
                h.data.vidInd = 1;
            end
        else
            i = imagesc(h.ui.ScreenAx, squeeze(h.data.Accumulator(:, :, h.data.vidInd)));
            set(h.ui.ScreenAx, 'CLim', h.data.vidClim);
            axis(h.ui.ScreenAx, 'image');
            title([h.data.vidChan ' channel at:' num2str(h.data.vidTimeVect(h.data.vidInd)) 's']);
            colorbar;
            h.data.vidInd = h.data.vidInd + 1;
            if( h.data.vidInd > h.data.vidLimit )
                h.data.vidInd = 1;
            end
        end
        
    end

    function GlobalSigChange(~, ~, ~)
        set(h.ui.EventsDispRegen, 'Enable', 'on'); 
    end

    function SrcChange(~,~,~)
        set(h.ui.EventsDispRegen, 'Enable', 'on'); 
    end

    function FiltOptChange(~, ~, ~)
        set(h.ui.EventsDispRegen, 'Enable', 'on'); 
    end

    function ret = FilterData(data, type)
        if( strcmp(type, 'IOI') )
            f = fdesign.lowpass('N,F3dB', 4, 0.4, h.data.AcqFreq);
            hpass = design(f,'butter');
            f = fdesign.lowpass('N,F3dB', 4, 1/60, h.data.AcqFreq);
            lpass = design(f,'butter');
            
            Hd= filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(data));
            Ld = filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(data));
            ret = single(Hd./Ld);
        elseif( strcmp(type, 'Fluo') )
            
        end
    end

    function MoveEventsDisp(~,~,~)
        
        idx = round(get(h.ui.EventsDispPan.Slider, 'Value'));
        T = h.data.EventBuf(1, :);
        d = h.data.EventBuf(idx+1, :);
        plot(h.ui.EventsDispPan.Ax, T, d, 'k');
        line(h.ui.EventsDispPan.Ax, [0 0], [min(d) max(d)],'Color', 'g', 'LineStyle','--');
        line(h.ui.EventsDispPan.Ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [min(d) max(d)],...
            'Color', 'r', 'LineStyle','--');
        if( mean(d) > 0.5 )
            line(h.ui.EventsDispPan.Ax, [T(1) T(end)], [1 1],...
                'Color', 'k', 'LineStyle',':');
        else
            line(h.ui.EventsDispPan.Ax, [T(1) T(end)], [0 0],...
                'Color', 'k', 'LineStyle',':');
        end
        xlim(h.ui.EventsDispPan.Ax, [T(1), T(end)]);
        
        %Creer le checkbox associe
        set(h.ui.EventsDispPan.Cbox,'String',  int2str(idx) );        
        set(h.ui.EventsDispPan.Cbox,'Value',  h.data.EvntList(idx) );        
    end

    function PopulateEvntsDisplay(~, ~, ~)
         if( isempty(h.data.ROIs) )
             return;
         end
              
        sID = get(h.ui.ChannelSelector, 'Value');
        sStr = get(h.ui.ChannelSelector, 'String');
        SelectedSrc = sStr{sID};
        
        rID = get(h.ui.ROIsSelector, 'Value');
        rStr = get(h.ui.ROIsSelector, 'String');
        SelectedROI = rStr{rID};
           
        if( ValidateEvntSrc(SelectedSrc, SelectedROI) )
            
            waitdialog = dialog('Position',[500 500 250 150],'Name','Computing...');
            uicontrol('Parent', waitdialog, 'Style','text',...
                'Position',[20 80 210 40], 'String', 'Computing Events curve. Please wait...');
            pause(0.1);
            
            data = 0;
            StartPts = 1;
            isHbT = 0; isHb = 0;
            if( ~contains(SelectedSrc, {'Hb', 'Flow', 'Fluo'}) )
                eval(['StartPts = h.data.' SelectedSrc(1) '_eflag;']);
                eval(['data = h.data.' lower(SelectedSrc(1)) 'DatPtr;']);
                isDivide = 1;
                isHb = 0;
            elseif( contains(SelectedSrc, 'Flow') )
                StartPts = h.data.F_eflag;
                data = h.data.fDatPtr;
                isDivide = 1;
                isHb = 0;
            elseif( contains(SelectedSrc, 'Fluo') )
                StartPts = h.data.F_eflag;
                data = h.data.fDatPtr;
                isDivide = 0;
                isHb = 0;
            else
                StartPts = h.data.H_eflag;
                isDivide = 0;
                isHb = 1;
                if( SelectedSrc == 'HbR' )
                    eval('data = h.data.hrDatPtr;');
                elseif( SelectedSrc == 'HbO' )
                    eval('data = h.data.hoDatPtr;');
                else
                    isHbT = 1;
                end
            end
            eLen = floor(h.data.AcqFreq*...
                (h.data.Stim.StimLength + h.data.Stim.InterStim_min));
            
            if( eLen == 0 )
               eLen = h.data.NFrames;
            end
            
            if( strcmp(SelectedROI, 'AllPixels') )
                mask = ones(size(h.data.Map));
            else
                idx = arrayfun(@(a) strcmp(h.data.ROIs{a}.name, SelectedROI), 1:size(h.data.ROIs,2));
                mask = h.data.ROIs{idx == 1}.mask;                                                                       
            end
            
            if( h.data.Stim.NbStim == 1 )
                h.ui.EventsDispPan.Visible = false;
            else
                h.ui.EventsDispPan.Visible = true;
                h.ui.EventsDispPan.min = 1 - (h.data.Stim.NbStim)*0.6;
                if( h.ui.EventsDispPan.min > 0 )
                    h.ui.EventsDispPan.min = 0;
                end
            end
            
            if( h.data.Stim.NbStim > 1 )
                T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.StimLength + h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
                h.data.EventBuf = zeros(h.data.Stim.NbStim + 1, eLen, 'single');
                h.data.EventBuf(1,:) = T;
            else
                T = linspace(0, eLen/h.data.AcqFreq, eLen);
                h.data.EventBuf = zeros(h.data.Stim.NbStim + 1, eLen, 'single');
                h.data.EventBuf(1,:) = T;
            end
            for indE = 1:h.data.Stim.NbStim
                
                if( isHbT == 0 )
                    d = data.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                else
                    d1 = h.data.hoDatPtr.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    d2 = h.data.hrDatPtr.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    d = d1+d2;
                end
                d = reshape(d, [], eLen);
                Glob = mean(d, 1);
                d = mean(d(mask(:)==1,:),1);
                
                if( (isHb == 0) && get(h.ui.GlobalOpt, 'Value') )
                    d = d./Glob;
                end
                if( (isHb == 0) && get(h.ui.FilteringOpt, 'Value') )
                    d = FilterData( d, 'IOI');
                end
                
                %Detrend
                Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                if( isDivide )
                    d = d./L;
                else
                    d = d - L;
                end
                
                h.data.EventBuf(indE + 1, :) = d;

            end
            MoveEventsDisp();
            
        end
        delete(waitdialog);
        MeanRecalculation();
        set(h.ui.EventsDispRegen, 'Enable', 'off'); 
    end

    function OnEditEventsClicked(~, src, ~)
        V = get(src.Source,'Value');
        h.data.EvntList(str2double(get(src.Source, 'String'))) = V;
        MeanRecalculation();
        h.flags.saveEvnts = true;
    end

    function MeanRecalculation()
        if( h.data.Stim.NbStim > 1 )
            d = zeros(1, floor(h.data.AcqFreq*...
                (h.data.Stim.StimLength + h.data.Stim.InterStim_min)));
            for indE = 1:length(h.data.EvntList)
                if(h.data.EvntList(indE))
                    d = d + h.data.EventBuf(indE + 1,:);
                end
            end
            d = d/sum(h.data.EvntList);
            T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.StimLength + h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, length(d));            
        else
            d = h.data.EventBuf(2,:);
            T = h.data.EventBuf(1,:);
        end        
        
        plot(h.ui.EventsMeanPan.Ax, T, d, 'k');
        line(h.ui.EventsMeanPan.Ax, [0 0], [min(d) max(d)],'Color', 'g', 'LineStyle','--');
        line(h.ui.EventsMeanPan.Ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [min(d) max(d)],...
            'Color', 'r', 'LineStyle','--');
        sID = get(h.ui.ChannelSelector, 'Value');
        sStr = get(h.ui.ChannelSelector, 'String');
        SelectedSrc = sStr{sID};
        if( isempty(strfind(SelectedSrc, 'Hb')) && isempty(strfind(SelectedSrc, 'Fluo')))
            line(h.ui.EventsMeanPan.Ax, [T(1) T(end)], [1 1],...
                'Color', 'k', 'LineStyle',':');
        else
            line(h.ui.EventsMeanPan.Ax, [T(1) T(end)], [0 0],...
                'Color', 'k', 'LineStyle',':');
        end
        xlim(h.ui.EventsMeanPan.Ax,[T(1), T(end)]);
    end

    function bRet = ValidateEvntSrc(cSrc, rSrc)
        if( cSrc(1) == 'G' )
            bRet = h.flags.IsThereGreen;
        elseif( cSrc(1) == 'R' )
            bRet = h.flags.IsThereRed;
        elseif( cSrc(1) == 'Y' )
            bRet = h.flags.IsThereYellow;
        elseif( cSrc(1) == 'F' )
            bRet = h.flags.IsThereFlow || h.flags.IsThereFluo;
        elseif( cSrc(3) == 'O' )
            bRet = h.flags.IsThereHbO;
        elseif( cSrc(3) == 'R' )
            bRet = h.flags.IsThereHbR;
        elseif( cSrc(3) == 'T' )
            bRet = h.flags.IsThereHbT;
        else
            bRet = false;
        end
        
        if( ~strcmp(rSrc, 'AllPixels') && ~any(arrayfun(@(r) strcmp(rSrc, h.data.ROIs{r}.name), 1:size(h.data.ROIs,2))) )
            bRet = false;
        end
    end

    function RefreshLoop(Option)
        if( strcmp(Option, 'ROIs') )
            RefreshROIsList();
            RefreshMapDisplay();
        elseif( strcmp(Option, 'Events') )
            if( h.data.Stim.NbStim > 1 )
                PopulateEvntsDisplay([],[],[]);
            end
        elseif( strcmp(Option, 'All') )
            RefreshROIsList();
            RefreshMapDisplay();
            if( h.data.Stim.NbStim > 1 )
                PopulateEvntsDisplay([],[],[]);
            end
        end
    end

    function RefreshMapDisplay(~,~,~)
        lst = get(h.ui.ROIsMap, 'Children');
        while( ~isempty(lst) )
           delete(lst(1));
           lst = get(h.ui.ROIsMap, 'Children');
        end
        
        if( isfield(h.ui, 'MskAxes') )
            for indA = 1:size(h.ui.MskAxes,2)
                cla(h.ui.MskAxes(indA));
            end
            h.ui = rmfield(h.ui, 'MskAxes');
        end
        
        cmap = gray(64);
        image(h.ui.ROIsMap, round(h.data.Map*64), 'AlphaData', 1);
        colormap(h.ui.ROIsMap, cmap);
        axis(h.ui.ROIsMap, 'off');
        
        NbRois = size(h.data.ROIs,2);
        for indR = 1:NbRois
            m = h.data.ROIs{indR}.mask;
            h.ui.MskAxes(indR) = axes('Parent', h.ui.ROIs, 'Position',[0.23 0.05 0.75 0.92]);
            image(h.ui.MskAxes(indR), m, 'AlphaData', m.*0.25);
            colormap(h.ui.MskAxes(indR), cat(1,[0 0 0], h.data.ROIs{indR}.color));
            axis(h.ui.MskAxes(indR), 'off');
        end 
    end

    function AddNewROI(~,~,~)
        h.flags.saveROIS = true;
        
        lst = get(h.ui.ROIsMap, 'Children');
        while( ~isempty(lst) )
            delete(lst(1));
            lst = get(h.ui.ROIsMap, 'Children');
        end
           
        if( isfield(h.ui, 'MskAxes') )
            for indA = 1:size(h.ui.MskAxes,2)
                cla(h.ui.MskAxes(indA));
            end
            h.ui = rmfield(h.ui, 'MskAxes');
        end
        
        image(h.ui.ROIsMap, round(h.data.Map*64));
        colormap(h.ui.ROIsMap, gray(64));
        axis(h.ui.ROIsMap, 'off');
        prompt = {'Name:'};
        dlg_title = 'New ROI';
        num_lines = 1;
        Tmp = {};
        if( ~isempty(h.data.ROIs) )
            Tmp = cell(size(h.data.ROIs,2),1);
            for indR = 1:size(h.data.ROIs,2)
                Tmp{indR} = h.data.ROIs{indR}.name;
            end
        end
        NewROI_ID = 1;
        while( sum(ismember(Tmp,['ROI_' int2str(NewROI_ID)])) )
            NewROI_ID = NewROI_ID + 1;
        end
        def = {['ROI_' int2str(NewROI_ID)]};
        
        answer = inputdlg(prompt, dlg_title,num_lines,def);
        if( ~isempty(answer) )
            Type = questdlg('ROI selection method:','Method selection',...
                'Circle', 'Polygon', 'Surround', 'Circle');
            
            h_im = get(h.ui.ROIsMap,'Children');
            
            switch Type
                case 'Circle'
                    e = imellipse(h.ui.ROIsMap);
                    mask = createMask(e, h_im(end));
                    delete(e);
                    clear e;
                case 'Rectangle'
                    r = imrect(h.ui.ROIsMap);
                    mask = createMask(r, h_im(end));
                    delete(r);
                    clear r;
                case 'Free Hand'
                    f = imfreehand(h.ui.ROIsMap);
                    mask = createMask(f, h_im(end));
                    delete(f);
                    clear f;
                case 'Polygon'
                    p = impoly(h.ui.ROIsMap);
                    mask = createMask(p, h_im(end));
                    delete(p);
                    clear p;
                 case 'Surround'
                     [orig, valid] = listdlg('PromptString', 'From wich ROI?',...
                         'SelectionMode', 'single',...
                         'ListString', Tmp);
                     if( ~valid )
                         return;
                     end
                     width = inputdlg( 'Width of the surround area:', 'Width', [1, 50], {'25'}); 
                     width = str2double(width{1});
                     
                     m = imfill(h.data.ROIs{orig}.mask,'holes');
                     mask = imdilate(m, strel('disk',width)) & ~m;
                     
            end
            h.data.ROIs{end+1} = struct('name',answer, 'mask', mask, 'color',  [0 0 1]);           
        end
        
        RefreshLoop('ROIs');
    end

    function RemoveROI(~,~,~)
        
        if( isempty(h.data.ROIs) )
            return;
        end
       
        Tmp = cell(size(h.data.ROIs,2),1);
        for indR = 1:size(h.data.ROIs,2)
            Tmp{indR} = h.data.ROIs{indR}.name;
        end
       
        [sel, valid] = listdlg('PromptString', 'Select the ROI to be removed:',...
            'SelectionMode', 'single',...
            'ListString', Tmp);
                     
        if( ~valid )
            return;
        end
        
        h.data.ROIs(sel) = [];  
        h.flags.saveROIS = true;
        RefreshLoop('ROIs');    
    end
        
    function OnChangeColorClicked(~, src, ~)
        h.flags.saveROIS = true;
        OldC = get(src.Source,'BackgroundColor');
        NewC = uisetcolor(OldC);
        ID = get(src.Source, 'UserData');
        h.data.ROIs{ID}.color = NewC;
        set(src.Source,'BackgroundColor',NewC);
        RefreshLoop('ROIs');
    end

    function RefreshROIsList(~ ,~, ~)
        lst = get(h.ui.ROIsPanel, 'Children');
        while( ~isempty(lst) )
            delete(lst(1));
            lst = get(h.ui.ROIsPanel, 'Children');
        end
        
        if( ~isempty(h.data.ROIs) )
            Str = cell(size(h.data.ROIs,2) + 1,1);
            Str{1} = 'AllPixels';
            for indR = 1:size(h.data.ROIs,2)
                uicontrol('Style','text','Parent',h.ui.ROIsPanel,...
                    'Units', 'normalized', 'Position',[0.01 0.99-indR*0.05 0.75 0.05],...
                    'String',h.data.ROIs{indR}.name, 'HorizontalAlignment','left');
                h.ui.ROIsColor(indR) = uicontrol('Style','pushbutton',...
                    'Parent',h.ui.ROIsPanel,'Units', 'Normalized',...
                    'Position',[0.80 0.99-indR*0.05 0.175 0.05], ...
                    'BackGroundColor', h.data.ROIs{indR}.color,...
                    'Callback',@OnChangeColorClicked,...
                    'UserData', indR);
                Str{indR+1} = h.data.ROIs{indR}.name;
            end
            set( h.ui.ROIsSelector, 'String', Str );
        end
    end

    function SaveROIs(~,~,~)
        if( h.flags.saveROIS )
            ROIs = h.data.ROIs; %#ok<*NASGU>
            save(h.paths.ROIsFile, 'ROIs');
        end
        h.flags.saveROIS = false;
        SaveEvnts();
    end

    function SaveEvnts(~,~,~)
        if( h.flags.saveEvnts )
            E = h.data.EvntList;
            save(h.paths.EVNTsFile, 'E');
        end
        h.flags.saveEvnts = false;
    end

    function LoadROIs(~,~,~)
%         [FileName,PathName,FilterIndex] = uigetfile('ROIs.mat');
%         if( any(FileName ~= 'ROIs.mat') || FilterIndex == 0 )
%             return;
%         end
%         Tmp = load([PathName FileName]);
%         if( any(size(Tmp.ROIs{1}.mask) ~= size(h.data.Map)) )
%             msgbox('ROIs file selected does not fit data dimension.','Load ROIs');
%             return;
%         end
%         h.data.ROIs = Tmp.ROIs;
%         h.flags.saveROIS = true;
%         RefreshLoop('ROIs');
    end

    function my_closereq(~,~)
        if( h.flags.saveROIS )
            selection = questdlg('Do you want to save modified ROIs list?',...
                'Before closing...',...
                'Yes','No','Yes');
            if( strcmp(selection, 'Yes') )
                SaveROIs();
                SaveEvnts();
            end
        end
        
        delete(gcf);
    end

end
