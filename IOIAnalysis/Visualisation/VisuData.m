function Ret = VisuData()

%%% Todo List:
% 1- Option to set spacial filter size
% 2- Hemodynamique Correction
% 3- Correction artefact de stim?

dParams.Folder = '';
dParams.Chan = '';
OldChan = '';
Data = [];
cData = [];
eData = [];
Ims = [];
Infos = [];
currentPixel = [1 1];
Conditions = [];
CondSequence = 0;
StimTrig = 1;
AcqInfoStream = [];

% Figures:
%Principale
hParams.figP = uifigure('Name', 'Parametres', 'NumberTitle','off',...
    'Position', [20 300 250 700], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'CloseRequestFcn', @FermeTout);
%Raw data:
hParams.figR = uifigure('Name', 'Images', 'NumberTitle','off',...
    'Position', [285 510 500 500], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off',...
    'WindowButtonDownFcn', @ChangePtPos, 'CloseRequestFcn', @NeFermePas);
%Corr data:
hParams.figC = uifigure('Name', 'Correlation', 'NumberTitle','off',...
    'Position', [285 200 250 250], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off', ...
    'CloseRequestFcn', @NeFermePas);
%Decours Temporel:
hParams.figT = uifigure('Name', 'Signal Temporel', 'NumberTitle','off',...
    'Position', [600 200 750 250], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off',...
    'CloseRequestFcn', @NeFermePas);
%Par Condition:
hParams.figE = uifigure('Name', 'Condition', 'NumberTitle','off',...
    'Position', [800 510 500 500], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off',...
    'CloseRequestFcn', @NeFermePas,'WindowButtonDownFcn', @ChangePtPos);

% GUI
%Pour le chemin d'acces vers le data a visualiser:
hParams.ExpLabel = uilabel(hParams.figP, 'Text','Experience:',...
    'Position',[5, 655, 100, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left');
hParams.ExpEdit = uieditfield(hParams.figP, 'text',...
    'Value', 'Choisir un dossier', 'Position',[5, 640, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left');
hParams.ExpPb = uibutton(hParams.figP, 'push',...
    'Text', '', 'Position',[200, 640, 35, 29], 'BackgroundColor', 'w',...
    'Icon', 'FolderIcon.png', 'ButtonPushedFcn', @ChangeFolder);

%Pre-Analyse:
hParams.PreAPB = uibutton(hParams.figP, 'push', ...
    'Text','Pre-Analyse', 'Position',[45, 590, 150, 35],...
    'BackgroundColor','w', 'FontName', 'Calibri', 'FontSize', 12,...
    'ButtonPushedFcn', @RunPreAna, 'visible', 'off');
hParams.PreALabel = uilabel(hParams.figP, 'Text','Pre-Analyse en cours. Patientez svp...',...
    'Position',[20, 550, 200, 35], 'BackgroundColor','w', ...
    'FontName', 'Calibri', 'FontSize', 12, 'visible', 'off');

%Cannal a utiliser:
hParams.ChanLabel = uilabel(hParams.figP, 'Text','Canal d''imagerie:',...
    'Position',[5, 590, 200, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.ChanPopMenu = uidropdown(hParams.figP, 'Items', {'Choisir'},...
    'Position',[5, 575, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10, 'visible', 'off',...
    'ValueChangedFcn', @OuvrirData);

%Type d'experience:
hParams.TypeLabel = uilabel(hParams.figP, 'Text','Type d''enregistrement',...
    'Position',[5, 525, 200, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.TypePopMenu = uidropdown(hParams.figP, 'Items',{'Choisir', 'RestingState', 'Episodique'},...
    'Position',[5, 505, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,...
    'visible', 'off', 'ValueChangedFcn', @ChangeType);

%Interaction Communes:
hParams.dFsFPB = uibutton(hParams.figP, 'push',...
    'Text','DF/F', 'Position',[20, 450, 75, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @DFsF);
hParams.GSRPB = uibutton(hParams.figP, 'push',...
    'Text','GSR', 'Position',[20, 400, 75, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @GSR);
hParams.Hemodyn = uibutton(hParams.figP, 'push',...
    'Text','Hemo', 'Position',[135, 450, 75, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12, 'Enable', 'on',...
    'visible', 'off', 'ButtonPushedFcn', @HemodynCorr);
hParams.SpatFilt = uibutton(hParams.figP, 'push',...
    'Text','F-Spacial', 'Position',[135, 400, 75, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @FiltreSpat);
hParams.Print = uibutton(hParams.figP, 'push',...
    'Text','Sauvegarde Figs', 'Position',[45, 15, 150, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @Print);

% Pour l'episodique:
% Fichier Vpixx reference:
hParams.VpixxLabel = uilabel(hParams.figP, 'Text','Fichier Stimulation:',...
    'Position',[5, 360, 100, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.VpixxEdit = uieditfield(hParams.figP, 'text',...
    'Value', 'Choisir un fichier', 'Position',[5, 345, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left','visible', 'off');
hParams.VpixxPb = uibutton(hParams.figP, 'push',...
    'Text', '', 'Position',[200, 345, 35, 29], 'BackgroundColor', 'w',...
    'Icon', 'FolderIcon.png', 'ButtonPushedFcn', @SelectFichier,'visible', 'off');
%Entree Analogique:
hParams.StimChanLabel = uilabel(hParams.figP, 'Text','Entree Analogique:',...
    'Position',[5, 305, 100, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.StimChanPopMenu = uidropdown(hParams.figP, 'Items',{'Choisir'},...
    'Position',[5, 290, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,...
    'visible', 'off', 'ValueChangedFcn', @ChangeStimSignal);
%Timing des decoupages:
hParams.PreStimLabel = uilabel(hParams.figP, 'Text','PreStim (s):',...
    'Position',[5, 250, 75, 25], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.PreStimEdit = uieditfield(hParams.figP, 'numeric',...
    'Value', 1, 'Position',[150, 250, 75, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'right',...
    'visible', 'off', 'ValueChangedFcn', @TimingValidation);
hParams.StimLabel = uilabel(hParams.figP, 'Text','Stim (s):',...
    'Position',[5, 225, 75, 25], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.StimEdit = uieditfield(hParams.figP, 'numeric',...
    'Value', 3, 'Position',[150, 225, 75, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'right',...
    'visible', 'off', 'ValueChangedFcn', @TimingValidation);
hParams.PostStimLabel = uilabel(hParams.figP, 'Text','PostStim (s):',...
    'Position',[5, 200, 75, 25], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.PostStimEdit = uieditfield(hParams.figP, 'numeric',...
    'Value', 5, 'Position',[150, 200, 75, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'right',...
    'visible', 'off', 'ValueChangedFcn', @TimingValidation, ...
    'Enable', 'off');
%Preparation du data:
hParams.SegEvnt = uibutton(hParams.figP, 'push',...
    'Text','Decoupage', 'Position',[45, 150, 150, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @StimDecoupe);


% Visualisation des images brutes:
% Graph:
hParams.axR1 = uiaxes(hParams.figR, 'Position', [67.5 100 375 375]);
% Boutons:
hParams.CurrentImageSl = uislider( hParams.figR,... 
    'Value', 1, 'Limits', [1 2], 'MajorTicks', [1 2],...
    'Position',[50, 90, 400, 3], 'ValueChangedFcn', @ChangeImage);
hParams.CI_MinLabel = uilabel( hParams.figR, 'Text', 'Minimum:', ...
    'Position',[100, 15, 75, 25],'FontName', 'Calibri', 'FontSize', 12,...
    'HorizontalAlignment', 'left');
hParams.CI_MaxLabel = uilabel( hParams.figR, 'Text', 'Maximum:', ...
    'Position',[300, 15, 75, 25],'FontName', 'Calibri', 'FontSize', 12,...
    'HorizontalAlignment', 'left');
hParams.CI_Min_Edit = uieditfield( hParams.figR, 'numeric',...
    'Value', 1, 'Position',[175, 15, 50, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left',...
    'ValueChangedFcn', @AdjustImage);
hParams.CI_Max_Edit = uieditfield( hParams.figR, 'numeric',...
    'Value', 4096, 'Position',[375, 15, 50, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left',...
    'ValueChangedFcn', @AdjustImage);

% Visualisation Corrélation
%Graph:
hParams.axC1 = uiaxes(hParams.figC, 'Position', [5 5 240 240]);

% Visualisation Décours Temporel:
%Graph:
hParams.axT1 = uiaxes(hParams.figT, 'Position', [5 5 740 240]);

%Visuatisation des Conditions:
hParams.axE1 = uiaxes(hParams.figE, 'Position', [67.5 100 375 375]);
hParams.CondSl = uislider( hParams.figE,... 
    'Value', 1, 'Limits', [1 2], 'MajorTicks', [1 2],...
    'Position',[50, 90, 400, 3], 'ValueChangedFcn', @RefreshImageCond);
hParams.Cond_MinLabel = uilabel( hParams.figE, 'Text', 'Min:', ...
    'Position',[315, 15, 25, 25],'FontName', 'Calibri', 'FontSize', 12,...
    'HorizontalAlignment', 'left');
hParams.Cond_MaxLabel = uilabel( hParams.figE, 'Text', 'Max:', ...
    'Position',[400, 15, 25, 25],'FontName', 'Calibri', 'FontSize', 12,...
    'HorizontalAlignment', 'left');
hParams.Cond_Min_Edit = uieditfield( hParams.figE, 'numeric',...
    'Value', 1, 'Position',[340, 15, 50, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left',...
    'ValueChangedFcn', @RefreshImageCond);
hParams.Cond_Max_Edit = uieditfield( hParams.figE, 'numeric',...
    'Value', 4096, 'Position',[425, 15, 50, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left',...
    'ValueChangedFcn', @RefreshImageCond);
hParams.Cond_PlayStop = uibutton(hParams.figE, 'push',...
    'Text','Play', 'Position',[15, 15, 45, 25], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'on', 'ButtonPushedFcn', @PlayStopCond);
hParams.Cond_Sel = uidropdown(hParams.figE, 'Items',{'Choisir'},...
    'Position',[75, 15, 105, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,...
    'visible', 'on', 'ValueChangedFcn', @NewImageCond);
hParams.Cond_Reps = uidropdown(hParams.figE, 'Items',{'Choisir'},...
    'Position',[200, 15, 85, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,...
    'visible', 'on', 'ValueChangedFcn', @NewImageCond);

% Initialisation de l'interface:
ChangeMode('Ouverture');

% Fonctions et Callbacks:
    function ChangeFolder(~,~,~)
        selpath = uigetdir(path);
        if( selpath == 0 )
            return;
        end
        
        if( ~strcmp(selpath(end), filesep) )
            selpath = strcat(selpath, filesep);
        end
        dParams.sFolder = selpath;
        hParams.ExpEdit.Value = selpath;
        ChargerDossier();
    end

    function ChargerDossier()
         %Validation du dossier
        list = dir([dParams.sFolder '*.dat']);
        
        Channels{1} = 'Choisir un Canal';
        for ind = 1:size(list,1) 
            Channels{ind+1} = list(ind).name;
        end
        
        hParams.ChanPopMenu.Items = Channels;
        hParams.ChanPopMenu.Value = Channels{1};
        if( size(list,1) == 0 )
            ChangeMode('PreAna');
        else
            ChangeMode('SelectParams');
        end
       % CheckDefaultParams();
       FigsOnTop();
       AcqInfoStream = ReadInfoFile(dParams.sFolder);
       hParams.StimChanPopMenu.Items = {'Choisir', 'CameraTrig', 'StimInterne'};
       for ind = 1:(AcqInfoStream.AINChannels - 2)
           hParams.StimChanPopMenu.Items{end+1} = ['Analog In #' int2str(ind)];
       end
    end

    function OuvrirData(~,~,~)
        
        if( iscell(hParams.ChanPopMenu.Items) )
            dParams.Chan = hParams.ChanPopMenu.Value;
            if( contains(hParams.ChanPopMenu.Items{1}, 'Choisir') )
                hParams.ChanPopMenu.Items = hParams.ChanPopMenu.Items(2:end);
            end
        else
            dParams.Chan = hParams.ChanPopMenu.Value;
        end
        
        if( ~strcmp(dParams.Chan, OldChan) )
            fid = fopen([dParams.sFolder dParams.Chan]);
            Data = fread(fid,inf, 'single=>single');
            Tmp = dir([dParams.sFolder 'Data_*.mat']);
            Infos = matfile([dParams.sFolder Tmp(1).name]);
            Data = reshape(Data, Infos.datSize(1,1), Infos.datSize(1,2),[]);
            Data = imresize3(Data, [256 256 size(Data,3)]);
            fclose(fid);
            
            hParams.GSRPB.Enable = 'on';
            hParams.dFsFPB.Enable = 'on';
            hParams.HemodynCorr.Enable = 'on';
        end    
        OldChan = dParams.Chan;
        hParams.dFsFPD.Enable = 'on'; 
        
        ChangeMode('SelectParams')
       
        hParams.CurrentImageSl.Limits = [1 size(Data,3)];
        hParams.CurrentImageSl.Value = 1;
        hParams.CurrentImageSl.MajorTicks = 1:1000:size(Data,3);
        
        ResetScale();
        ChangeImage();
        
        if( strcmp(AcqInfoStream.Camera_Model, 'D1024') )
            hParams.BNoise.Enable = 'on';
        end
    end

    function ResetScale()
        P = prctile(Data(:), [2.5 97.5]);
        Mini = P(1);
        Maxi = P(2);
        hParams.CI_Min_Edit.Value = double(Mini);
        hParams.CI_Max_Edit.Value = double(Maxi);
    end

    function ChangeImage(~,~,~)
        Id = round(hParams.CurrentImageSl.Value);
        Im = imresize(squeeze(Data(:,:,Id)),[256 256]);
        imagesc(hParams.axR1, Im);
        colorbar(hParams.axR1);
        caxis(hParams.axR1, [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value]);
        title(hParams.axR1,['Image #: ' int2str(Id)]);
        axis(hParams.axR1, 'off', 'image');
        hold(hParams.axR1, 'on');
        plot(hParams.axR1, currentPixel(1), currentPixel(2), 'or');
        hold(hParams.axR1, 'off');
        DecoursTemp();    
    end

    function AdjustImage(~,~,~)
        
        caxis(hParams.axR1, [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value]);
        DecoursTemp();
    end

    function RunPreAna(~,~,~)
        
        prompt = {'Binning Spatial (1 si aucun binning):',...
            'Binning Temporel(1 si aucun binning):',...
            'Redifinir une Region d''interet? (0:non; 1:oui)',...
            'Ignorer le signal de stimulation interne du systeme? (0:non; 1:oui)'};
        dlgtitle = 'Pre-Analyse';
        dims = [1 50];
        definput = {'1','1', '0', '0'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        
        if( isempty(answer) )
            return;
        end
        hParams.PreALabel.Visible = 'on';
        pause(0.01);
        try
        ImagesClassification(dParams.sFolder, str2double(answer{1}),...
            str2double(answer{2}), str2double(answer{3}), str2double(answer{4}));
        catch e
            disp(e);
            hParams.PreALabel.Value = 'Une erreur est survenue durant la pre-analyse.';
        end
        hParams.PreALabel.Visible = 'off';
        ChargerDossier();
    end

    function ChangeType(~,~,~)
        
        dParams.sExpType = hParams.TypePopMenu.Value;
        if( contains(hParams.TypePopMenu.Items{1}, 'Choisir') )
           hParams.TypePopMenu.Items = hParams.TypePopMenu.Items(2:end);
        end
        if( any(contains(hParams.ChanPopMenu.Items, 'Choisir')) )
            ChangeMode('SelectParams')
        elseif( strcmp(dParams.sExpType, 'RestingState') )
            ChangeMode('RestingState');
            CorrMap();
        else
            ChangeMode('Episodique Prepa');
        end
    end
        
    function ChangeMode(NewMode)
       
        switch( NewMode )
            case 'Ouverture'
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                hParams.figE.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.TypeLabel.Visible = 'off';
                hParams.TypePopMenu.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'off';
                hParams.ChanPopMenu.Visible = 'off';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.Hemodyn.Visible = 'off';
                hParams.SpatFilt.Visible = 'off';
                hParams.VpixxLabel.Visible = 'off';
                hParams.VpixxEdit.Visible = 'off';
                hParams.VpixxPb.Visible = 'off';
                hParams.StimChanPopMenu.Visible = 'off'; 
                hParams.PreStimLabel.Visible = 'off'; 
                hParams.PreStimEdit.Visible = 'off'; 
                hParams.StimLabel.Visible = 'off'; 
                hParams.StimEdit.Visible = 'off'; 
                hParams.PostStimLabel.Visible = 'off'; 
                hParams.PostStimEdit.Visible = 'off'; 
                hParams.SegEvnt.Visible = 'off'; 
                hParams.Print.Visible = 'off'; 
                
            case 'PreAna'
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                hParams.figE.Visible = 'off';
                
                hParams.PreAPB.Visible = 'on';
                hParams.TypeLabel.Visible = 'off';
                hParams.TypePopMenu.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'off';
                hParams.ChanPopMenu.Visible = 'off';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.Hemodyn.Visible = 'off';
                hParams.BNoise.Visible = 'off';
                hParams.SpatFilt.Visible = 'off';
                hParams.VpixxLabel.Visible = 'off';
                hParams.VpixxEdit.Visible = 'off';
                hParams.VpixxPb.Visible = 'off';
                hParams.StimChanPopMenu.Visible = 'off';
                hParams.StimChanLabel.Visible = 'off'; 
                hParams.PreStimLabel.Visible = 'off'; 
                hParams.PreStimEdit.Visible = 'off'; 
                hParams.StimLabel.Visible = 'off'; 
                hParams.StimEdit.Visible = 'off'; 
                hParams.PostStimLabel.Visible = 'off'; 
                hParams.PostStimEdit.Visible = 'off'; 
                hParams.SegEvnt.Visible = 'off'; 
                hParams.Print.Visible = 'off'; 
                
            case 'SelectParams'
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                hParams.figE.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.BNoise.Visible = 'off';
                hParams.Hemodyn.Visible = 'off';
                hParams.SpatFilt.Visible = 'off';
                hParams.VpixxLabel.Visible = 'off';
                hParams.VpixxEdit.Visible = 'off';
                hParams.VpixxPb.Visible = 'off';
                hParams.StimChanPopMenu.Visible = 'off';
                hParams.StimChanLabel.Visible = 'off';
                hParams.PreStimLabel.Visible = 'off'; 
                hParams.PreStimEdit.Visible = 'off'; 
                hParams.StimLabel.Visible = 'off'; 
                hParams.StimEdit.Visible = 'off'; 
                hParams.PostStimLabel.Visible = 'off'; 
                hParams.PostStimEdit.Visible = 'off'; 
                hParams.SegEvnt.Visible = 'off'; 
                hParams.Print.Visible = 'off'; 
                
            case 'RestingState'
                hParams.figR.Visible = 'on';
                hParams.figC.Visible = 'on';
                hParams.figT.Visible = 'on';
                hParams.figE.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'on';
                hParams.GSRPB.Visible = 'on';
                hParams.BNoise.Visible = 'on';
                hParams.Hemodyn.Visible = 'on';
                hParams.SpatFilt.Visible = 'on';
                hParams.VpixxLabel.Visible = 'off';
                hParams.VpixxEdit.Visible = 'off';
                hParams.VpixxPb.Visible = 'off';
                hParams.StimChanPopMenu.Visible = 'off';
                hParams.StimChanLabel.Visible = 'off';
                hParams.PreStimLabel.Visible = 'off'; 
                hParams.PreStimEdit.Visible = 'off'; 
                hParams.StimLabel.Visible = 'off'; 
                hParams.StimEdit.Visible = 'off'; 
                hParams.PostStimLabel.Visible = 'off'; 
                hParams.PostStimEdit.Visible = 'off'; 
                hParams.SegEvnt.Visible = 'off'; 
                hParams.Print.Visible = 'on'; 
            
            case 'Episodique Prepa'    
                hParams.figR.Visible = 'on';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                hParams.figE.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'on';
                hParams.GSRPB.Visible = 'on';
                hParams.BNoise.Visible = 'on';
                hParams.Hemodyn.Visible = 'on';
                hParams.SpatFilt.Visible = 'on';
                hParams.VpixxLabel.Visible = 'on';
                hParams.VpixxEdit.Visible = 'on';
                hParams.VpixxPb.Visible = 'on';
                hParams.StimChanPopMenu.Visible = 'on';
                hParams.StimChanLabel.Visible = 'on';
                hParams.PreStimLabel.Visible = 'off'; 
                hParams.PreStimEdit.Visible = 'off'; 
                hParams.StimLabel.Visible = 'off'; 
                hParams.StimEdit.Visible = 'off'; 
                hParams.PostStimLabel.Visible = 'off'; 
                hParams.PostStimEdit.Visible = 'off'; 
                hParams.SegEvnt.Visible = 'off'; 
                hParams.Print.Visible = 'off'; 
            
            case 'Episodique Decoup'    
                hParams.figR.Visible = 'on';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                hParams.figE.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'on';
                hParams.GSRPB.Visible = 'on';
                hParams.BNoise.Visible = 'on';
                hParams.Hemodyn.Visible = 'on';
                hParams.SpatFilt.Visible = 'on';
                hParams.VpixxLabel.Visible = 'on';
                hParams.VpixxEdit.Visible = 'on';
                hParams.VpixxPb.Visible = 'on';
                hParams.StimChanPopMenu.Visible = 'on';
                hParams.StimChanLabel.Visible = 'on';
                hParams.PreStimLabel.Visible = 'on'; 
                hParams.PreStimEdit.Visible = 'on'; 
                hParams.StimLabel.Visible = 'on'; 
                hParams.StimEdit.Visible = 'on'; 
                hParams.PostStimLabel.Visible = 'on'; 
                hParams.PostStimEdit.Visible = 'on'; 
                hParams.SegEvnt.Visible = 'on'; 
                hParams.Print.Visible = 'off'; 
                
            case 'Episodique'      
                hParams.figR.Visible = 'on';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'on';
                hParams.figE.Visible = 'on';
                
                hParams.PreAPB.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'on';
                hParams.GSRPB.Visible = 'on';
                hParams.BNoise.Visible = 'on';
                hParams.Hemodyn.Visible = 'on';
                hParams.SpatFilt.Visible = 'on';
                hParams.VpixxLabel.Visible = 'on';
                hParams.VpixxEdit.Visible = 'on';
                hParams.VpixxPb.Visible = 'on';
                hParams.StimChanPopMenu.Visible = 'on';
                hParams.StimChanLabel.Visible = 'on';
                hParams.PreStimLabel.Visible = 'on'; 
                hParams.PreStimEdit.Visible = 'on'; 
                hParams.StimLabel.Visible = 'on'; 
                hParams.StimEdit.Visible = 'on'; 
                hParams.PostStimLabel.Visible = 'on'; 
                hParams.PostStimEdit.Visible = 'on'; 
                hParams.SegEvnt.Visible = 'on'; 
                hParams.Print.Visible = 'on'; 
                         
            otherwise
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.TypeLabel.Visible = 'off';
                hParams.TypePopMenu.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'off';
                hParams.ChanPopMenu.Visible = 'off';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.BNoise.Visible = 'off';
                hParams.Hemodyn.Visible = 'off';
                hParams.SpatFilt.Visible = 'off';
                hParams.VpixxLabel.Visible = 'off';
                hParams.VpixxEdit.Visible = 'off';
                hParams.VpixxPb.Visible = 'off';
                hParams.StimChanPopMenu.Visible = 'off';
                hParams.StimChanLabel.Visible = 'off';
                hParams.PreStimLabel.Visible = 'off'; 
                hParams.PreStimEdit.Visible = 'off'; 
                hParams.StimLabel.Visible = 'off'; 
                hParams.StimEdit.Visible = 'off'; 
                hParams.PostStimLabel.Visible = 'off'; 
                hParams.PostStimEdit.Visible = 'off'; 
                hParams.SegEvnt.Visible = 'off'; 
                hParams.Print.Visible = 'off'; 
        end
        dParams.Mode = NewMode;
        
    end

    function DFsF(~,~,~)
        hParams.dFsFPB.Text = '';
        hParams.dFsFPB.Icon = 'Sablier.png';
        drawnow;
        if( mean(reshape(Data,[],size(Data,3)),1) < 0.5 )
            Data = Data + 1;
        end
        
        dims = size(Data);
        Data = reshape(Data, [], dims(3));
        
        if( contains(hParams.ChanPopMenu.Value, 'f') )
            lp_cutoff = 1/5;
            hp_cutoff = Infos.Freq/3;
        else 
            lp_cutoff = 1/120;
            hp_cutoff = 1;
        end
            
        f = fdesign.lowpass('N,F3dB', 4, lp_cutoff, Infos.Freq);
        lpass = design(f,'butter');
        lpData = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(Data)'))';
        Data = Data./lpData;
        clear lpData;
        f = fdesign.lowpass('N,F3dB', 4, hp_cutoff, Infos.Freq);
        lpass = design(f,'butter');
        Data = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(Data)'))';        
        
        Data = reshape(Data,dims);
        ResetScale();
        ChangeImage();
                
        if( strcmp(dParams.sExpType, 'RestingState') )
            CorrMap();
            DecoursTemp();
        elseif( strcmp(dParams.Mode, 'Episodique') ) 
            StimDecoupe();
            NewImageCond();
        end
        hParams.dFsFPB.Icon = '';
        hParams.dFsFPB.Text = 'DF/F';
        hParams.dFsFPB.Enable = 'off';
        drawnow;
    end

    function GSR(~,~,~)
        hParams.GSRPB.Text = '';
        hParams.GSRPB.Icon = 'Sablier.png';
        drawnow;
        dims = size(Data);
        Data = reshape(Data, [], dims(3));
        Signal = mean(Data,1);
        Signal = Signal./mean(Signal);
        X = [ones(1,dims(3)); Signal];
        B = X'\Data';
        A = (X'*B)';
        Data = Data - A;
        Data = reshape(Data, dims);
        ResetScale();
        ChangeImage();
        hParams.GSRPB.Enable = 'off';
        
        if( strcmp(dParams.sExpType, 'RestingState') )
            CorrMap();
            DecoursTemp();
        elseif( strcmp(dParams.Mode, 'Episodique') )
            StimDecoupe();
            NewImageCond();
        end
        hParams.GSRPB.Icon = '';
        hParams.GSRPB.Text = 'GSR';
        drawnow;
    end

    function FiltreSpat(~,~,~)
        answer = inputdlg('Largeur du filtre a appliquer:', 'Filtre Gaussien',...
            1, {'1.5'});
        if( isempty(answer) )
            return;
        end
        answer = str2num(answer{:});
        
        Data = imgaussfilt(Data, answer);
        ResetScale();
        ChangeImage();
        if( strcmp(dParams.sExpType, 'RestingState') )
            CorrMap();
            DecoursTemp();
        elseif( strcmp(dParams.Mode, 'Episodique') )
            StimDecoupe();
            NewImageCond();
        end
    end

    function BNoise(~,~,~)
        
        h = waitbar(0, 'Filtre Banding Noise');
        for indF = 1:size(Data,3)
           Data(:,:,indF) = single(xRemoveStripesVertical(squeeze(Data(:,:,indF)), 4 + nextpow2(256), 'db4', 2));
            waitbar(indF/size(Data,3),h);
        end
        close(h);
        
        ChangeImage();
       if( strcmp(dParams.sExpType, 'RestingState') )
            CorrMap();
            DecoursTemp();
        elseif( strcmp(dParams.Mode, 'Episodique') )
            StimDecoupe();
            NewImageCond();
        end
    end    

    function HemodynCorr(~,~,~)
        
        if( ~exist([dParams.sFolder 'fluoCorrected.dat'],'file') )
            try
                Data = HemoCorrection(dParams.sFolder, {'Red', 'Amber'});
            catch
                return;
            end
            fid = fopen([dParams.sFolder 'fluoCorrected.dat'],'w');
            fwrite(fid, Data, 'single');
            fclose(fid);
        else
            fid = fopen([dParams.sFolder 'fluoCorrected.dat']);
            Data = fread(fid, inf, '*single');
            Data = reshape(Data, Infos.datSize(1,1), Infos.datSize(1,2),[]);
            fclose(fid);
        end
        hParams.GSRPB.Enable = 'on';
        hParams.dFsFPB.Enable = 'on';
        hParams.HemodynCorr.Enable = 'off';
        ResetScale();
        ChangeImage();
        if( strcmp(dParams.sExpType, 'RestingState') )
            CorrMap();
            DecoursTemp();
        elseif( strcmp(dParams.Mode, 'Episodique') )
            StimDecoupe();
            NewImageCond();
        end
        
    end
       
    function CorrMap(~,~,~)
       cData = imresize(Data,[64 64]);
       cData = cData - mean(cData,3);
       cData = corr(reshape(cData,[],size(cData,3))');
       
       RefreshCorrMap();
    end

    function RefreshCorrMap(~,~,~)
        if( isempty(cData) )
            return;
        end
        
        Facteur = 64/256;
        Id = floor((currentPixel(1)-1)*Facteur)*64 + floor(currentPixel(2)*Facteur);
        if( Id < 1 ) 
            Id = 1;
        end
        if( Id > 4096)
            Id = 4096;
        end
        imagesc(hParams.axC1, reshape(cData(Id,:),64,64),[0 1]);
        colorbar(hParams.axC1);
        axis(hParams.axC1, 'off', 'image');
    end

    function ChangePtPos(Obj, Evnt)
       if( strcmp(Obj.Name, 'Images') )
            Pos = round(Obj.Children(6).CurrentPoint);
       elseif (strcmp(Obj.Name, 'Condition') )
           Pos = round(Obj.Children(9).CurrentPoint);
       end
       Pos = Pos(1,1:2);
       if( any(Pos < 1) | any(Pos > 255) )
           return;
       end
       currentPixel = Pos;
       ChangeImage();
       if( strcmp(dParams.sExpType, 'RestingState') )
           RefreshCorrMap();
           DecoursTemp();
       elseif( strcmp(dParams.Mode, 'Episodique') )
           RefreshImageCond();
       end
    end

    function DecoursTemp()
       % disp('DecoursTemp')
        if( strcmp(dParams.Mode, 'Episodique') )
            Cond_id = find(cellfun(@(x) strcmp(x, hParams.Cond_Sel.Value), hParams.Cond_Sel.Items));
            if( strcmp(hParams.Cond_Reps.Value, 'Moyenne') )
                temporal_signal = squeeze(mean(eData(currentPixel(2),currentPixel(1),:,Cond_id,:),5));
            else
                Reps_id = find(cellfun(@(x) strcmp(x, hParams.Cond_Reps.Value), hParams.Cond_Reps.Items)) - 1;
                temporal_signal = squeeze(eData(currentPixel(2),currentPixel(1),:,Cond_id,Reps_id));
            end
        else
            temporal_signal = squeeze(Data(currentPixel(2),currentPixel(1),:));
        end
        hold(hParams.axT1, 'off');
        plot(hParams.axT1, temporal_signal);
        hold(hParams.axT1, 'on');
        if( strcmp(dParams.Mode, 'Episodique') )
            X = hParams.CondSl.Value;
            Y = [hParams.Cond_Min_Edit.Value, hParams.Cond_Max_Edit.Value];
        else
            X = hParams.CurrentImageSl.Value;
            Y = [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value];
        end
        line(hParams.axT1,[X, X],...
            [Y(1), Y(2)],'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        ylim(hParams.axT1, [Y(1), Y(2)]);
    end

    function SelectFichier(~,~,~)
        [filename, pathname] = uigetfile('*.txt', 'Choisir le fichier Vpixx');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        end
                
        %Lire le fichier txt:
        filetext = fileread([pathname filename]);
        
        expr = '[^\n]*Start Trial';
        
        file_lim = strfind(filetext, 'SORTED');
        fileread_info = regexp(filetext(1:file_lim), expr, 'match');
        CondSequence = zeros(length(fileread_info), 1);
        for ind = 1:length(fileread_info)
            a = sscanf(fileread_info{ind}(1:strfind(fileread_info{ind}, '[')), '%f');
            if isempty(a) ==1
                error('No condition listed in Vpixx file')
            end
            CondSequence(ind) = a(2);
        end
        
        file_lim = strfind(filetext, 'SUMMARY');
        expr = '\[[^\]]*\]';
        Conditions = regexp(filetext(file_lim:end), expr, 'match');        
        hParams.VpixxEdit.Value = filename;
        
        hParams.Cond_Sel.Items = Conditions;
        hParams.Cond_Reps.Items = {'Moyenne'};
        for ind = 1:round(length(CondSequence)/size(Conditions,2))
            hParams.Cond_Reps.Items{end+1} = ['Rep #' int2str(ind)];
        end       
                
        FigsOnTop();
    end

    function StimDecoupe(~,~,~)
        
        TotalLength = floor((hParams.PostStimEdit.Value + ...
            hParams.PreStimEdit.Value + hParams.StimEdit.Value)*Infos.Freq); 
        
        eData = [];
        Cond = size(Conditions,2);
        Reps = size(CondSequence,1)/size(Conditions,2);
        rCntr = ones(1,Cond);
        eData = zeros(256, 256, TotalLength, Cond, Reps); 
        Debut = StimTrig - round(hParams.PreStimEdit.Value*Infos.Freq);
        Fin = Debut + TotalLength - 1;
        for indE = 1:length(CondSequence)
            eData(:,:,:,CondSequence(indE), rCntr(CondSequence(indE))) = ...
                Data(:,:, Debut(indE):Fin(indE));
            rCntr(CondSequence(indE)) = rCntr(CondSequence(indE)) + 1;           
        end
        
        hParams.CondSl.Limits = [1 TotalLength];
        hParams.CondSl.Value = 1;
        hParams.CondSl.MajorTicks = 1:10:TotalLength;
        
        ChangeMode('Episodique');
        NewImageCond();
    end

    function TimingValidation(~,~,~)
        TotalLength =  mean(diff(StimTrig,1,1))/Infos.Freq;
        
        hParams.PostStimEdit.Value = floor(TotalLength -...
            hParams.PreStimEdit.Value - hParams.StimEdit.Value); 
      
    end
      
    function ChangeStimSignal(~,~,~)
        dParams.StimChan = hParams.StimChanPopMenu.Value;
        if( contains(hParams.StimChanPopMenu.Items{1}, 'Choisir') )
            hParams.StimChanPopMenu.Items = hParams.StimChanPopMenu.Items(2:end);
        end
        
        idx = find(cellfun(@(x) strcmp(x, hParams.StimChanPopMenu.Value), hParams.StimChanPopMenu.Items));
        %Lire entrees analogiques:
        AnalogIN = [];
        aiFilesList = dir([dParams.sFolder 'ai*.bin']);
        for ind = 1:size(aiFilesList,1) %for each "ai_.bin" file:
            data = memmapfile([dParams.sFolder aiFilesList(ind).name], ...
                'Offset', 5*4, 'Format', 'double', 'repeat', inf);
            tmp = data.Data; %Read data
            tmp = reshape(tmp, AcqInfoStream.AISampleRate, AcqInfoStream.AINChannels, []); %reshape data based on the number of AIs
            tmp = permute(tmp,[1 3 2]); %Permutation of dimension because AIs are interlaced when saved.
            tmp = reshape(tmp,[], AcqInfoStream.AINChannels);
            AnalogIN = [AnalogIN; tmp]; %Concatenation with previous ai_.bin files
        end
        clear tmp ind data;
        Signal = AnalogIN(:,idx);
        Cam = AnalogIN(:,1);
        [~, Cam] = ischange(Cam);
        Cam = find(diff(Cam,1,1)>0);
        clear AnalogIN aiFilesList;
        
        NbColors = sum(contains(fieldnames(AcqInfoStream), 'Illumination'));
        Colors = {};
        for ind = 1:NbColors
            eval(['Colors{' int2str(ind) '} = AcqInfoStream.Illumination' int2str(ind) '.Color;']);
        end
        tag = dParams.Chan(1);
        switch(tag)
            case 'r'
                idx = find(contains(Colors, 'Red'));
            case 'y'
                idx = find(contains(Colors, 'Amber'));
            case 'g'
                idx = find(contains(Colors, 'Green'));
            case 'f'
                idx = find(contains(Colors, 'Fluo'));
        end
        
        [~, Signal] = ischange(Signal,'Threshold', 2);
        Signal = Signal(Cam);
        Signal = Signal(idx:NbColors:end);
        Signal = Signal(1:size(Data,3));
        StimTrig = find(diff(Signal,1,1)>0);
        
        hParams.PreStimEdit.Limits = [0, (StimTrig(1) - 1)/Infos.Freq];
        hParams.StimEdit.Limits = [1, mean(diff(StimTrig,1,1))/Infos.Freq];
        hParams.PostStimEdit.Limits = [0, mean(diff(StimTrig,1,1))/Infos.Freq];
        
        ChangeMode('Episodique Decoup');
    end
    
    function NewImageCond(~,~,~)
        %disp('NewImageCond');
        Ims = [];
        Cond_id = find(cellfun(@(x) strcmp(x, hParams.Cond_Sel.Value), hParams.Cond_Sel.Items));
        
        if( strcmp(hParams.Cond_Reps.Value, 'Moyenne') )
            Ims = squeeze(mean(eData(:,:,:,Cond_id,:),5));            
        else
            Reps_id = find(cellfun(@(x) strcmp(x, hParams.Cond_Reps.Value), hParams.Cond_Reps.Items)) - 1;
            Ims = squeeze(eData(:,:,:,Cond_id,Reps_id));            
        end
        
        P = prctile(Ims(:),[1 99]);
        hParams.Cond_Min_Edit.Value = P(1);
        hParams.Cond_Max_Edit.Value = P(2);
        RefreshImageCond();
    end

    function RefreshImageCond(~,~,~)                     
      %  disp('RefreshImageCond');
        Id = round(hParams.CondSl.Value);
        Im = imresize(squeeze(Ims(:,:,Id)),[256 256]);
        imagesc(hParams.axE1, Im);
        colorbar(hParams.axE1);
        caxis(hParams.axE1, [hParams.Cond_Min_Edit.Value, hParams.Cond_Max_Edit.Value]);
        title(hParams.axE1,['Image #: ' int2str(Id)]);
        axis(hParams.axE1, 'off', 'image');
        hold(hParams.axE1, 'on');
        plot(hParams.axE1, currentPixel(1), currentPixel(2), 'or');
        hold(hParams.axE1, 'off');
        % Update line on plot:
        DecoursTemp();
    end

    function PlayStopCond(~,~,~)
        
    end

    function Print(~,~,~)
           % Sauvegarde les images affichees dans des fichiers en format .png 
       filename = inputdlg('Nom du fichier de sauvegarde:' , 'Sauvegarde de figures');
      % disp('Sauvegarde de figures en cours...');
       fields = fieldnames(hParams);
       fields = regexp(fields, 'fig\w*[^P]', 'match');fields = [fields{:}];
       for i = 1:length(fields)
           if strcmp(hParams.(fields{i}).Visible,'on')
               name = [filename{:} '_' hParams.(fields{i}).Name '.png'];
               idx = arrayfun(@(x) isa(x, 'matlab.ui.control.UIAxes'), hParams.(fields{i}).Children);
               handle = hParams.(fields{i}).Children(idx);
               fig = figure('Visible', 'off', 'Position', hParams.(fields{i}).Position);
               ax = axes(fig);
               copyobj(handle.Children, ax);
               % Save all parameters of the UIAxes
               uiAxParams = get(handle);
               uiAxParamNames = fieldnames(handle);
               % Get list of editable params in new axis
               editableParams = fieldnames(set(ax));
               % Remove the UIAxes params that aren't editable in the new axes (add others you don't want)
               badFields = uiAxParamNames(~ismember(uiAxParamNames, editableParams));
               badFields = [badFields; 'Parent'; 'Children'; 'XAxis'; 'YAxis'; 'ZAxis';'Position';'OuterPosition'];
               uiAxGoodParams = rmfield(uiAxParams,badFields);
               % set editable params on new axes
               ax.Title.String = handle.Title.String;
               set(ax, uiAxGoodParams);
               
               if strcmp(fields{i}, 'figE')
                   condName = hParams.Cond_Sel.Value;
                   repName = hParams.Cond_Reps.Value;
                   title(strjoin({condName, repName}, '--'));
               end
               saveas(fig, fullfile(hParams.ExpEdit.Value, name), 'png')               
           end
       end
     %  disp('Fini!')
       uiwait(msgbox(['Figures sauvegardes dans ' hParams.ExpEdit.Value], 'Sauvegarde reussite'));
       close all
       ChangeImage();
       RefreshImageCond();
    end

    function NeFermePas(~,~,~)
        
    end

    function FermeTout(~,~,~)
        delete(hParams.figT);
        delete(hParams.figR);
        delete(hParams.figC);
        delete(hParams.figE);
        delete(hParams.figP);        
    end

    function FigsOnTop()
        if( strcmp(hParams.figR.Visible, 'on') )
            figure(hParams.figR);
        end
        if( strcmp(hParams.figT.Visible, 'on') )
            figure(hParams.figT);
        end
        if( strcmp(hParams.figC.Visible, 'on') )
            figure(hParams.figC);
        end
        if( strcmp(hParams.figE.Visible, 'on') )    
            figure(hParams.figE);
        end
        figure(hParams.figP);
    end
end