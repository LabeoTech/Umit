function Data = MappingAnalysis(FolderPath, DetectionThreshold, HysteresisThreshold)
%Infos files opening
if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end
InfoStim = Read_OptoGenParams_File(FolderPath);
%InfoAcq = ReadInfoFile(FolderPath);
InfoAcq.AISampleRate = 10000;
InfoAcq.AINChannels = 11;

%Analog Inputs reading:
aiFiles = dir([FolderPath 'ai_*.bin']);
AnalogIn = [];
for indF = 1:length(aiFiles)
    fid = fopen([FolderPath aiFiles(indF).name]);
    fseek(fid,5*4,'bof');
    dat = fread(fid,inf,'double');
    dat = reshape(dat, InfoAcq.AISampleRate, InfoAcq.AINChannels, []);
    dat = permute(dat,[2 1 3]);
    dat = reshape(dat,11,[]);
    AnalogIn = [AnalogIn, dat];
end
clear dat indF aiFiles fid

%Stimulation detection: 
Stim = AnalogIn(2,:);
iStim = find(diff(Stim,1,2)>2.50);
Xpos = InfoStim.Positions(:,2);
Ypos = InfoStim.Positions(:,1);
sX = length(unique(Xpos));
sY = length(unique(Ypos));
iX = round((sX-1)*(Xpos - min(Xpos(:)))/(max(Xpos(:)) - min(Xpos(:))) + 1);
iY = round((sY-1)*(Ypos - min(Ypos(:)))/(max(Ypos(:)) - min(Ypos(:))) + 1);
Xpos = unique(Xpos);
Ypos = unique(Ypos);
if( size(InfoStim.Positions,2) == 3 )
    Lpow = unique(InfoStim.Positions(:,3));
    iPow = ones(size(iY));
    for indP = 1:length(Lpow)
        iPow(InfoStim.Positions(:,3) == Lpow(indP)) = indP;
    end
else
    Lpow = InfoStim.LaserPower;
    iPow = ones(size(iY));
end

%Analog Inputs reshaping based on stimulations:
%Channels for Paw 1: AI #1(4 in Matlab), AI #2(5 in Matlab), AI #3(6 in Matlab)
%Channels for Paw 2: AI #5(8 in Matlab), AI 6#(9 in Matlab), AI #7(10 in Matlab)
f = fdesign.lowpass('N,F3dB', 4, 250, 10000);
lpass = design(f,'butter');
clear f
AnalogIn = filtfilt(lpass.sosMatrix, lpass.ScaleValues, AnalogIn');
clear lpass

NbRep = mean(accumarray(iX,1))/sY;
NbRep = NbRep./length(Lpow);
Reps = ones(sY, sX, length(Lpow));
Mp1 = zeros(sY, sX, 4000, NbRep, length(Lpow));
Mp2 = zeros(sY, sX, 4000, NbRep, length(Lpow));
for indS = 1:length(iStim)
    idx = iStim(indS) + (-500:3499);
    indR = Reps(iY(indS), iX(indS), iPow(indS));
    Reps(iY(indS), iX(indS), iPow(indS)) = Reps(iY(indS), iX(indS), iPow(indS)) + 1;
    Mp1(iY(indS), iX(indS), :, indR, iPow(indS)) = sqrt(AnalogIn(idx,4).^2 +...
        AnalogIn(idx,5).^2 + AnalogIn(idx,6).^2);
    Mp2(iY(indS), iX(indS), :, indR, iPow(indS)) = sqrt(AnalogIn(idx,8).^2 +...
        AnalogIn(idx,9).^2 + AnalogIn(idx,10).^2);
end
clear aIn dat idx ind iStim iX iY Stim sX sY iPow Reps

rawMp1 = flipud(rot90(Mp1));
rawMp2 = flipud(rot90(Mp2));
dims = size(rawMp1);

Mp1 = squeeze(mean(rawMp1,4));
Mp2 = squeeze(mean(rawMp2,4));

Mp1 = bsxfun(@minus, Mp1, mean(Mp1,3));
Mp2 = bsxfun(@minus, Mp2, mean(Mp2,3));

mMp1 = mean(Mp1(:,:,1:500,:),3);
sMp1 = std(Mp1(:,:,1:500,:),0,3);
zMp1 = bsxfun(@rdivide, bsxfun(@minus, Mp1, mMp1), sMp1);

mMp2 = mean(Mp2(:,:,1:500),3);
sMp2 = std(Mp2(:,:,1:500),0,3);
zMp2 = bsxfun(@rdivide, bsxfun(@minus, Mp2, mMp2), sMp2);
clear mMp* sMp*

[vMaxZ_1, tMax1] = max(zMp1,[],3);
vMaxZ_1 = squeeze(vMaxZ_1); tMax1 = squeeze(tMax1);
MapActiv1 = vMaxZ_1 >= DetectionThreshold;
MapActiv1(tMax1<500) = false;
MapActiv1(tMax1>1000) = false;
tMax1 = tMax1 - 500;
tMax1(~MapActiv1) = 0;
tMax1 = tMax1/10;

[vMaxZ_2, tMax2] = max(zMp2,[],3);
vMaxZ_2 = squeeze(vMaxZ_2); tMax2 = squeeze(tMax2);
MapActiv2 = vMaxZ_2 >= DetectionThreshold;
MapActiv2(tMax2<500) = false;
MapActiv2(tMax2>1000) = false;
tMax2 = tMax2 - 500;
tMax2(~MapActiv2) = 0;
tMax2 = tMax2/10;

if( length(dims) < 5 )
    tThresh_1 = zeros(dims(1), dims(2));
    tThresh_2 = zeros(dims(1), dims(2));
    D5 = 1;
else
    tThresh_1 = zeros(dims(1), dims(2), dims(5));
    tThresh_2 = zeros(dims(1), dims(2), dims(5));
    D5 = dims(5);
end
clear ind*
for indR = 1:D5
    for indX = 1:dims(2)
        for indY = 1:dims(1)
            tmp1 = find(zMp1(indY, indX,500:end, indR) >= DetectionThreshold, 1,'first');
            tmp1 = 500 + tmp1;
            tmp1 = tmp1 - find(zMp1(indY, indX,tmp1:-1:500, indR) <= HysteresisThreshold, 1,'first');
            
            tmp2 = find(zMp2(indY, indX,500:end, indR) >= DetectionThreshold, 1,'first');
            tmp2 = 500 + tmp2;
            tmp2 = tmp2 - find(zMp2(indY, indX,tmp2:-1:500, indR) <= HysteresisThreshold, 1,'first');
               
            if( isempty(tmp1) )
                tThresh_1(indY, indX, indR) = 0;
            else
                tThresh_1(indY, indX, indR) = tmp1;
            end
            if( isempty(tmp2) )
                tThresh_2(indY, indX, indR) = 0;
            else
                tThresh_2(indY, indX, indR) = tmp2;
            end
        end
    end
end
tThresh_1 = tThresh_1 - 500;
tThresh_1 = tThresh_1/10;
tThresh_1(~MapActiv1) = 0;

tThresh_2 = tThresh_2 - 500;
tThresh_2 = tThresh_2/10;
tThresh_2(~MapActiv2) = 0;
clear tmp* ind*

for indR = 1:length(Lpow)
    eval(['Data.tOnset_Paw1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR).*tThresh_1(:,:,indR);']);
    eval(['Data.tOnset_Paw2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR).*tThresh_2(:,:,indR);']);
    eval(['Data.tMax_Paw1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR).*tMax1(:,:,indR);']);
    eval(['Data.tMax_Paw2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR).*tMax2(:,:,indR);']);
    eval(['Data.vMax_Paw1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR).*squeeze(max(Mp1(:,:,:,indR),[],3));']);
    eval(['Data.vMax_Paw2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR).*squeeze(max(Mp2(:,:,:,indR),[],3));']);
    eval(['Data.vMin_Paw1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR).*squeeze(min(Mp1(:,:,:,indR),[],3));']);
    eval(['Data.vMin_Paw2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR).*squeeze(min(Mp2(:,:,:,indR),[],3));']);
    eval(['Data.zMax_Paw1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR).*squeeze(max(zMp1(:,:,:,indR),[],3));']);
    eval(['Data.zMax_Paw2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR).*squeeze(max(zMp2(:,:,:,indR),[],3));']);
    eval(['Data.zMin_Paw1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR).*squeeze(min(zMp1(:,:,:,indR),[],3));']);
    eval(['Data.zMin_Paw2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR).*squeeze(min(zMp2(:,:,:,indR),[],3));']);
    eval(['Data.MapActiv1_P' int2str(Lpow(indR)) ' = MapActiv1(:,:,indR);']);
    eval(['Data.MapActiv2_P' int2str(Lpow(indR)) ' = MapActiv2(:,:,indR);']);
end    
Data.sagittal_axis = Xpos;
Data.coronal_axis = Ypos;
Data.Infos = InfoStim;

save([FolderPath 'DataMapping.mat'], 'Data');
%Generate Figures:
Xpix = round(InfoStim.RefX + Xpos/InfoStim.MMpPix);
Ypix = round(InfoStim.RefY + Ypos/InfoStim.MMpPix);
%Max: 
hfig = figure;
axRef = axes('Parent', hfig);
axOvr = axes('Parent', hfig);
axis(axOvr, 'off');

IList = dir([FolderPath '*.png']);
if( length(IList) > 1 )
    str = {IList.name};
    [v,c] = listdlg('PromptString','Select a file for Reference:',...
                'SelectionMode','single',...
                'ListString',str);
    if( c < 1 )
        IRef = zeros(1024,1024,3);
    else
        IList = IList(v);
        IRef = imread([FolderPath IList.name]);
    end
elseif( isempty(IList) )
    IRef = zeros(1024,1024,3);
else
    IRef = imread([FolderPath IList.name]);
end

Xaxmm = ((1:1024) - InfoStim.RefX)*InfoStim.MMpPix;
Yaxmm = ((1:1024) - InfoStim.RefY)*InfoStim.MMpPix;
imagesc(axRef, Xaxmm, Yaxmm, IRef);
hold(axRef,'on');
plot(axRef, 0, 0, 'or');
text(axRef, 0.1, -0.1, '{\beta}','FontSize',16,'FontWeight', 'bold', 'Color', 'r')


for indP = 1:length(Lpow)
    %MAX 3-4-5
    imagesc(axOvr, Ypos, Xpos, MapActiv1(:,:,indP).*squeeze(max(Mp1(:,:,:,indP),[],3)),...
        'AlphaData', MapActiv1(:,:,indP)*0.5);
    colormap(axRef,'gray');
    title(['Max Amplitude Channels Paw#1 Power ' int2str(Lpow(indP))])
    xlabel(axRef, 'Coronal axis (mm)');
    ylabel(axRef, 'Sagital axis (mm)');
    axis(axRef, 'image');
    axis(axOvr, 'image');
    axis(axOvr,'off');
    colorbar('AxisLocation','in');
    linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
        'CameraPosition', 'XLim', 'YLim'});
    saveas(hfig, [FolderPath 'MaxAmpPaw1_' int2str(Lpow(indP)) '.png']);
    
    %MAX 7-8-9
    imagesc(axOvr, Ypos, Xpos, MapActiv2(:,:,indP).*squeeze(max(Mp2(:,:,:,indP),[],3)),...
        'AlphaData', MapActiv2(:,:,indP)*0.5);
    colormap(axRef,'gray');
    title(['Max Amplitude Channels Paw#2 Power ' int2str(Lpow(indP))])
    xlabel(axRef, 'Coronal axis (mm)');
    ylabel(axRef, 'Sagital axis (mm)');
    axis(axRef, 'image');
    axis(axOvr, 'image');
    axis(axOvr,'off');
    colorbar('AxisLocation','in');
    linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
        'CameraPosition', 'XLim', 'YLim'});
    saveas(hfig, [FolderPath 'MaxAmpPaw2_' int2str(Lpow(indP)) '.png']);
    
    %Tmax 3-4-5
    imagesc(axOvr, Ypos, Xpos, MapActiv1(:,:,indR).*tMax1(:,:,indR),'AlphaData', MapActiv1(:,:,indP)*0.5);
    colormap(axRef,'gray');
    title(['Rising Time to Maximum Channels Paw#1 Power ' int2str(Lpow(indP))])
    xlabel(axRef, 'Coronal axis (mm)');
    ylabel(axRef, 'Sagital axis (mm)');
    axis(axRef, 'image');
    axis(axOvr, 'image');
    axis(axOvr,'off');
    colorbar('AxisLocation','in');
    linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
        'CameraPosition', 'XLim', 'YLim'});
    saveas(hfig, [FolderPath 'TmaxPaw1_' int2str(Lpow(indP)) '.png']);
    
    %Tmax 7-8-9
    imagesc(axOvr, Ypos, Xpos, MapActiv2(:,:,indR).*tMax2(:,:,indR),'AlphaData', MapActiv2(:,:,indP)*0.5);
    colormap(axRef,'gray');
    title(['Rising Time to Maximum Channels Paw#2 Power ' int2str(Lpow(indP))])
    xlabel(axRef, 'Coronal axis (mm)');
    ylabel(axRef, 'Sagital axis (mm)');
    axis(axRef, 'image');
    axis(axOvr, 'image');
    axis(axOvr,'off');
    colorbar('AxisLocation','in');
    linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
        'CameraPosition', 'XLim', 'YLim'});
    saveas(hfig, [FolderPath 'TmaxPaw2_' int2str(Lpow(indP)) '.png']);
    
    %TOnset 3-4-5
    imagesc(axOvr, Ypos, Xpos, MapActiv1(:,:,indR).*tThresh_1(:,:,indR), 'AlphaData', MapActiv1(:,:,indP)*0.5);
    colormap(axRef,'gray');
    title(['Onset Time to Maximum Channels Paw#1 Power ' int2str(Lpow(indP))])
    xlabel(axRef, 'Coronal axis (mm)');
    ylabel(axRef, 'Sagital axis (mm)');
    axis(axRef, 'image');
    axis(axOvr, 'image');
    axis(axOvr,'off');
    caxis(axOvr,[0 50]);
    colorbar('AxisLocation','in');
    linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
        'CameraPosition', 'XLim', 'YLim'});
    saveas(hfig, [FolderPath 'TonChanPaw1_' int2str(Lpow(indP)) '.png']);
    
    %TOnset 7-8-9
    imagesc(axOvr, Ypos, Xpos, MapActiv2(:,:,indR).*tThresh_2(:,:,indR), 'AlphaData', MapActiv2(:,:,indP)*0.5);
    colormap(axRef,'gray');
    title(['Onset Time to Maximum Channels Paw#2 Power ' int2str(Lpow(indP))])
    xlabel(axRef, 'Coronal axis (mm)');
    ylabel(axRef, 'Sagital axis (mm)');
    axis(axRef, 'image');
    axis(axOvr, 'image');
    axis(axOvr,'off');
    caxis(axOvr,[0 50]);
    colorbar('AxisLocation','in');
    linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
        'CameraPosition', 'XLim', 'YLim'});
    saveas(hfig, [FolderPath 'TonChanPaw2_' int2str(Lpow(indP)) '.png']);
end
close(hfig);

end