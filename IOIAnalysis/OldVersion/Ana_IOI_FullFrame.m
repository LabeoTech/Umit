function out = Ana_IOI_FullFrame(FolderName, verbose, b_tFilter, b_fitFilter, OStream)

%%%%%%%%%%%
% Opening %
%%%%%%%%%%%
%List of files to be included in the analysis
FileList = dir([FolderName filesep 'Data_*.mat']);
Ws = []; Hs = []; Ts = []; Fs = [];
if( isempty(FileList) )
    %No compatible files were found
    disp(['No data files found in ' FolderName ' Folder.']);
    disp('Analysis will not run');
    return;
end

%Green channel detected
IsThereGreen = false;
if( ~isempty(strfind([FileList.name],'green')) )
    IsThereGreen = true;
    Dat_Gptr = matfile([FolderName filesep 'Data_green.mat'], 'Writable', true);
    nrows = Dat_Gptr.datSize(1,1);
    ncols = Dat_Gptr.datSize(1,2);
    nframes = Dat_Gptr.datLength;
    Ws = ncols;
    Hs = nrows;
    Ts = nframes - 1;
    Fs = Dat_Gptr.Freq;
    if( ~exist(Dat_Gptr.datFile,'file') )
        Dat_Gptr.datFile = [FolderName filesep 'gChan.dat'];
    end
    gDatPtr = memmapfile(Dat_Gptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes
end
%Yellow channel detected
IsThereYellow = false;
if( ~isempty(strfind([FileList.name],'yellow')) )
    IsThereYellow = true;
    Dat_Yptr = matfile([FolderName filesep 'Data_yellow.mat'], 'Writable', true);
    nrows = Dat_Yptr.datSize(1,1);
    ncols = Dat_Yptr.datSize(1,2);
    nframes = Dat_Yptr.datLength;
    Ws = [Ws, ncols];
    Hs = [Hs, nrows];
    Ts = [Ts, nframes-1];
    Fs = [Fs, Dat_Yptr.Freq];
    if( ~exist(Dat_Yptr.datFile,'file') )
        Dat_Yptr.datFile = [FolderName filesep 'yChan.dat'];
    end
    yDatPtr = memmapfile(Dat_Yptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes
end
%Red channel detected
IsThereRed = false;
if( ~isempty(strfind([FileList.name],'red')) )
    IsThereRed = true;
    Dat_Rptr = matfile([FolderName filesep 'Data_red.mat'], 'Writable', true);
    nrows = Dat_Rptr.datSize(1,1);
    ncols = Dat_Rptr.datSize(1,2);
    nframes = Dat_Rptr.datLength;
    Ws = [Ws, ncols];
    Hs = [Hs, nrows];
    Ts = [Ts, nframes-1];
    Fs = [Fs, Dat_Rptr.Freq];
    if( ~exist(Dat_Rptr.datFile,'file') )
        Dat_Rptr.datFile = [FolderName filesep 'rChan.dat'];
    end
    rDatPtr = memmapfile(Dat_Rptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes
end

%Is all required colors available for HB calculation?
if( IsThereRed + IsThereYellow + IsThereGreen < 2 )
    disp('*** Impossible to compute Hb concentrations. More color channels needed.');
    fprintf('\n');
    return;
end

%Confirm data dimensions
if( length(unique(Ws)) > 1 )
    disp('Channels have unmatching dimensions');
    disp('Analysis will stop here.');
    return;
end
iWidth = double(Ws(1));
if( length(unique(Hs)) > 1 )
    disp('Channels have unmatching dimensions');
    disp('Analysis will stop here.');
    return;
end
iHeight = double(Hs(1));
NbFrames = double(min(Ts));
FreqHb = min(Fs);
clear Ws Hs Ts Fs;

%%%%%%%%%%%%%%%%%
% Output Config %
%%%%%%%%%%%%%%%%%
%Creation of output file
if( exist([FolderName filesep 'Data_Hbs.mat'], 'file') )
    delete([FolderName filesep 'Data_Hbs.mat']);
    delete([FolderName filesep 'HbO.dat']);
    delete([FolderName filesep 'HbR.dat']);
end
OutputFile = matfile([FolderName filesep 'Data_Hbs.mat']);
OutputFile.datFileHbO = [FolderName filesep 'HbO.dat'];
OutputFile.datFileHbR = [FolderName filesep 'HbR.dat'];
OutputFile.datLength = NbFrames;
OutputFile.datSize = [iWidth, iHeight];
OutputFile.Freq = FreqHb;
fHbO = fopen([FolderName filesep 'HbO.dat'], 'W');
fHbR = fopen([FolderName filesep 'HbR.dat'], 'W');

%Filtering Parameters
fbase = ceil(60*FreqHb);
fioi = ceil(2*FreqHb);
clear FileList nframes
ExpFun = @(P,x) P(1).*exp(-P(2).*x) + P(3).*exp(-P(4).*x) + P(5);
Opt = optimset(@fminsearch);
Opt.Display = 'off';
warning('OFF', 'MATLAB:rankDeficientMatrix');

%%%%%%%%%%%%%%%%%%%%%%%%%
% Hb Concentration Calc %
%%%%%%%%%%%%%%%%%%%%%%%%%
%HbO and HbR computation parameters
whichCurve = 'Silico';
rescaling_factor = 1e6;

baseline_hbt = 100;
baseline_hbo = 60;
baseline_hbr = 40;

eps_pathlength = ioi_epsilon_pathlength(whichCurve,baseline_hbt,baseline_hbo,baseline_hbr, 0);
%load('eps_path_MC.mat');
if( IsThereGreen && IsThereYellow && IsThereRed )
    A = eps_pathlength;
elseif( IsThereYellow && IsThereRed )
    A = [eps_pathlength(1,:); eps_pathlength(3,:)];
elseif( IsThereGreen && IsThereRed )
    A = [eps_pathlength(1,:); eps_pathlength(2,:)];
elseif( IsThereGreen && IsThereYellow )
    A = [eps_pathlength(2,:); eps_pathlength(3,:)];
end
Ainv=single(rescaling_factor*pinv(A)); % A Inv devrait etre 3x2 ou 2x3 (3=couleurs, 2=hbo,hbr)
clear which* rescaling_factor lambda* npoints baseline_* eps_pathlength A

%For each row, filter each channels and compute HbO, HbR...
PrcTags = linspace(1, double(iHeight), 11); indPr = 2;
if( ~isempty(OStream) )
    StaticStr = OStream.String;
end
for ind = 1:iHeight
    %Read, filter and detrend data
    if( IsThereRed )
        pR = rDatPtr.Data(double(ind):iHeight:(iHeight*iWidth*NbFrames));
        pR = reshape(pR, iWidth, [])';
        
        if( b_fitFilter )
            S = mean(pR,2);
            B = fminsearch(@(P) norm(double(S) - ExpFun(P,(1:size(pR,1))')),...
                [30 0.0025 20 0.015 double(mean(S))],Opt);
            Approx = ExpFun([B(1:4) 0],1:size(pR,1));
            Pred =[ones(1, size(pR,1)); linspace(0,1,size(pR,1)); linspace(0,1,size(pR,1)).^2; Approx]'; 
            b = Pred\pR;
            Rbase = (Pred*b);
        else
            Rbase = medfilt1(pR,fbase,[],1,'truncate');
        end
        if( b_tFilter )
            Rioi= medfilt1(pR,fioi,[],1,'truncate');
        else
            Rioi = pR;
        end
        Rnorm = Rioi./Rbase;
        clear Rioi Rbase pR;
    end
    if( IsThereYellow )
        pY = yDatPtr.Data(double(ind):iHeight:(iHeight*iWidth*NbFrames));
        pY = reshape(pY, iWidth, [])';
        
        if( b_fitFilter )
            S = mean(pY,2);
            B = fminsearch(@(P) norm(double(S) - ExpFun(P,(1:size(pY,1))')),...
                [30 0.0025 20 0.015 double(mean(S))],Opt);
            Approx = ExpFun([B(1:4) 0],1:size(pY,1));
            Pred =[ones(1, size(pY,1)); linspace(0,1,size(pY,1)); linspace(0,1,size(pY,1)).^2; Approx]';
            b = Pred\pY;
            Ybase = (Pred*b);
        else
            Ybase = medfilt1(pY,fbase,[],1,'truncate');
        end
        if( b_tFilter )
            Yioi = medfilt1(pY,fioi,[],1,'truncate');
        else
            Yioi = pY;
        end
        Ynorm = Yioi./Ybase;
        clear Yioi Ybase pY;
    end
    if( IsThereGreen )
        pG = gDatPtr.Data(ind:iHeight:(iHeight*iWidth*NbFrames));
        pG = reshape(pG, iWidth, [])';
        if( b_fitFilter )
            S = mean(pG,2);
            B = fminsearch(@(P) norm(double(S) - ExpFun(P,(1:size(pG,1))')),...
                [30 0.0025 20 0.015 double(mean(S))],Opt);
            Approx = ExpFun([B(1:4) 0],1:size(pG,1));
            Pred =[ones(1, size(pG,1)); linspace(0,1,size(pG,1)); linspace(0,1,size(pG,1)).^2; Approx]';
            b = Pred\pG;
            Gbase = (Pred*b);
        else
            Gbase = medfilt1(pG,fbase,[],1,'truncate');
        end
        if( b_tFilter )
            Gioi = medfilt1(pG,fioi,[],1,'truncate');
        else
            Gioi = pG;
        end
        Gnorm = Gioi./Gbase;
        clear Gioi Gbase pG;
    end
    
    %Compute HbO & HbR
    if( IsThereGreen && IsThereYellow && IsThereRed )
        Cchan = cat(2, Rnorm(:), Gnorm(:), Ynorm(:));
    elseif( IsThereYellow && IsThereRed )
        Cchan = cat(2, Rnorm(:), Ynorm(:));
    elseif( IsThereGreen && IsThereRed )
        Cchan = cat(2, Rnorm(:), Gnorm(:));
    elseif( IsThereGreen && IsThereYellow )
        Cchan = cat(2, Gnorm(:), Ynorm(:));
    end
    
    LogCchan = -log(Cchan);
    Hbs = Ainv*LogCchan';
    
    %Save
    fwrite(fHbO,reshape(Hbs(1,:), [], iWidth), 'single');
    %HbO(ind,:,:) = reshape(Hbs(1,:), [], iWidth)';
    fwrite(fHbR,reshape(Hbs(2,:), [], iWidth), 'single');
    %HbR(ind,:,:) = reshape(Hbs(2,:), [], iWidth)';
    
    if( ind >= PrcTags(indPr) )
        if( isempty(OStream) )
            fprintf('%d%% .. ', 10*(indPr-1));
        else
            OStream.String = sprintf('%s\r%s',...
                ['Completion: ' int2str(10*(indPr-1)) '%'],...
                StaticStr);
            drawnow;
        end
        indPr = indPr + 1;
    end
end
if( isempty(OStream) )
    fprintf('\n');
else
    OStream.String = StaticStr;
    OStream.String = sprintf('%s\r%s',...
        'Done.',...
        OStream.String);
    drawnow;
end
if( IsThereRed )
   OutputFile.Stim = Dat_Rptr.Stim;
elseif( IsThereGreen )
   OutputFile.Stim = Dat_Gptr.Stim;
else
   OutputFile.Stim = Dat_Yptr.Stim;
end
warning('ON', 'MATLAB:rankDeficientMatrix');

fclose(fHbO);
fclose(fHbR);
fHbO = fopen([FolderName filesep 'HbO.dat'], 'r+');
dat = fread(fHbO,inf,'single');
frewind(fHbO);
dat = reshape(dat, NbFrames, iWidth, iHeight);
dat = permute(dat, [3 2 1]);
fwrite(fHbO, dat,'single');
fclose(fHbO);
fHbR = fopen([FolderName filesep 'HbR.dat'], 'r+');
dat = fread(fHbR,inf,'single');
frewind(fHbR);
dat = reshape(dat, NbFrames, iWidth, iHeight);
dat = permute(dat, [3 2 1]);
fwrite(fHbR, dat,'single');
fclose(fHbR);
end
