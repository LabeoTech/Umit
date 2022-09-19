function out = IClass_BFly(FolderName, Binning)

%%%%DEFINES -> THESE CONSTANTS ARE HARDCODED!!!! DO NOT CHANGE THEM.
NOFPF = 256;
DEF_VISUEL = 0;

%%%%%%%%%%%%%%%%%%%%%
% Acq. Info file:
%%%%%%%%%%%%%%%%%%%%%
AcqInfoStream = ReadInfoFile(FolderName);


disp('Recovering stimulation parameters')
disp('**************************');

tAIChan = AcqInfoStream.AINChannels;

%%%%%%%%%%%%%%%%%%%%%
% Stimulation detected
%%%%%%%%%%%%%%%%%%%%%
if( AcqInfoStream.Stimulation > 0 )
    ReadAnalogsIn(FolderName, FolderName, AcqInfoStream);
else
    fprintf('No stimulation detected. \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images sequence validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgFilesList = dir([FolderName filesep 'img_0*.bin']);
aiFilesList = dir([FolderName filesep 'ai_*.bin']);


hWima = 5;
hWai = 5;
header = memmapfile([FolderName filesep imgFilesList(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);

nx=header.Data.header(2);
ny=header.Data.header(3);
frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx*ny*2 + 3*8;


NombreImage = 0;
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    NombreImage = NombreImage+size(data.Data,1);
end

idImg = zeros(NombreImage, 3);

for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    idImg((NOFPF*(ind-1)+1):(NOFPF*(ind-1)+size(data.Data,1)),:) = cell2mat(arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1),'UniformOutput',false))';
end
clear nx ny data header ind;

% Verbose
disp(['Opening of: ' FolderName]);
disp(['Number of Frames acquired: ' int2str(NombreImage)]);
disp(['Frames'' resolution: ' int2str(ImRes_XY(1)) ' pix X ' int2str(ImRes_XY(2)) ' pix']);
% end of Verbose

%%%%%%%%%%%%%%%%%%%%%
% Binning
%%%%%%%%%%%%%%%%%%%%%
if( Binning )
    %Verbose
    disp('Binning option is ON');
    %end of Verbose
    
    Rx = round(ImRes_XY(1)/Binning);
    Ry = round(ImRes_XY(2)/Binning);
    
else
    Rx = ImRes_XY(1);
    Ry = ImRes_XY(2);
end


%%%%%%%%%%%%%%%%%%%%%
% Analog inputs
%%%%%%%%%%%%%%%%%%%%%
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', hWai*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, AcqInfoStream.AISampleRate, tAIChan, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],tAIChan);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data;

%%%%
%Stimulation Params
%%%%
Str = [];
if( exist([FolderName filesep 'StimParameters.mat'], 'file') )
    load([FolderName filesep 'StimParameters.mat']);
    if( NbStim > 0 )
        Str = sprintf('%s\r%s%s\r%s%s%s',...
            'Stim detected: yes',...
            'Number of events: ', int2str(NbStim),...
            'Length of each event: ', int2str(StimLength), ' sec');
        bStim = 1;
    else
        Str = sprintf('%s',...
            'Stim detected: no');
        bStim = 0;
    end
    fprintf(Str);
    fprintf('\n');
    
else
    bStim = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Illumination Sequence
%%%%%%%%%%%%%%%%%%%%%%%
% Red           = 0001;
% Yellow/Amber  = 0010;
% Green         = 0100;
% Other         = 1000;
nbColors = sum(cellfun(@(X) contains(X,'Illumination'),fieldnames(AcqInfoStream)));

bFluo = 0; nFluo = 0;
bGreen = 0; nGreen = 0;
bRed = 0; nRed = 0;
bYellow = 0; nYellow = 0;
bSpeckle = 0;nSpeckle = 0;
for indC = 1:nbColors
   Tag = eval(['AcqInfoStream.Illumination' int2str(indC) ';']);
   switch(Tag)
       case 'Fluo'
           bFluo = 1;
           nFluo = indC;
       case 'Green'
           bGreen = 1;
           nGreen = indC;
       case 'Amber'
           bYellow = 1;
           nYellow = indC;
       case 'Red'
           bRed = 1;
           nRed = indC;
       case 'Speckle'
           bSpeckle = 1;
           nSpeckle = indC;
   end
end

Freq = AcqInfoStream.FrameRateHz;

if( bFluo )
       if( exist([FolderName filesep 'Data_Fluo.mat'],'file') )
            delete([FolderName filesep 'Data_Fluo.mat']);
        end
        fFluo = matfile([FolderName filesep 'Data_Fluo.mat'],'Writable',true);
        fFluo.datFile = [FolderName filesep 'fChan.dat'];
        fFluo.datSize = [Rx, Ry];
        fFluo.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
        fFluo.Freq = Freq/nbColors;
        cFluo = 1;
        fidF = fopen([FolderName filesep 'fChan.dat'],'w');
end    
if( bSpeckle )
        if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
            delete([FolderName filesep 'Data_speckle.mat']);
        end
        fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);
        fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
        fSpeckle.datSize = [Rx, Ry];
        fSpeckle.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
        fSpeckle.Freq = Freq/nbColors;
        cSpeckle = 1;
        fidS = fopen([FolderName filesep 'sChan.dat'],'w');

end
if( bRed )
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(floor(NombreImage/nbColors), 1, 'single');
    fRed.Freq = Freq/nbColors;
    cRed = 1;
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
end
if( bYellow )
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
    fYellow.Freq = Freq/nbColors;
    cYellow = 1;
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
end
if( bGreen )    
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(floor(NombreImage/nbColors), 1, 'single');
    fGreen.Freq = Freq/nbColors;
    cGreen = 1;
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
end

%Interpolation for bad or missing frames
%SkipNFirst = sum(idImg(:,1) == 0);
SkipNFirst = 0;
MissingOffset = idImg(:,2);
idImg(:,1) = idImg(:,1) + 1;
goodFrames = find(accumarray(idImg((SkipNFirst+1):end,1),1)==1)';
ConseqFromLeft = [1 diff(goodFrames,1,2)==1];
ConseqFromRight = fliplr([true diff(fliplr(goodFrames),1,2)==-1]);
goodFrames = goodFrames(ConseqFromLeft|ConseqFromRight);
badFrames = 1:max(goodFrames(:));
badFrames = badFrames(~ismember(badFrames, goodFrames));
%%% Lookup Table For missing frames
InterpLUT = zeros(8,1);
if( ~isempty(badFrames) )
    InterpLUT = zeros(8,size(badFrames,2));
    % 1: Frame before
    % 2: File Number where to find this frame
    % 3: Image ID in this file
    % 4: Frame After
    % 5: File Number where to find this frame
    % 6: Image ID in this file
    % 7: Ratio between these frame and the one missing
    % 8: Frame tag id
    for ind = 1:size(badFrames,2)
        tmpID = badFrames(ind);
        tmpBefore = tmpID - (nbColors:nbColors:(tmpID-1));
        idx = find(ismember(tmpBefore,idImg(:,1))&~ismember(tmpBefore,badFrames),1,'first');
        tmpBefore = tmpBefore(idx);
        InterpLUT(1,ind) = tmpBefore;
        idx = find(tmpBefore == idImg(:,1),1,'first');
        InterpLUT(2,ind) = floor((idx-1)/256) + 1;
        InterpLUT(3,ind) = rem((idx-1),256) + 1;
        
        tmpAfter = tmpID + (nbColors:nbColors:(NombreImage));
        idx = find(ismember(tmpAfter,idImg(:,1))&~ismember(tmpAfter,badFrames),1,'first');
        tmpAfter = tmpAfter(idx);
        InterpLUT(4,ind) = tmpAfter;
        idx = find(tmpAfter == idImg(:,1), 1, 'first');
        InterpLUT(5,ind) = floor((idx-1)/256) + 1;
        InterpLUT(6,ind) = rem((idx-1),256) + 1;
        
        tmpRatio =  (tmpID - InterpLUT(1,ind))./...
            (InterpLUT(4,ind) - InterpLUT(1,ind));
        InterpLUT(7,ind) = tmpRatio;
        InterpLUT(8,ind) = badFrames(ind);
    end
    clear tmpRatio tmpAfter tmpBefore;
    %%% Interpolation of missing frames
    TmpFrames = struct('framej',[], 'imgj',[]);
    for ind = 1:size(InterpLUT,2)
        dBefore = memmapfile([FolderName filesep...
            imgFilesList(InterpLUT(2,ind)).name],...
            'Offset', hWima*4 + (InterpLUT(3,ind) - 1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        dAfter = memmapfile([FolderName filesep...
            imgFilesList(InterpLUT(5,ind)).name],...
            'Offset', hWima*4 + (InterpLUT(6,ind) - 1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        TmpFrames(ind).imgj = uint16(round(InterpLUT(7,ind)*double(dAfter.Data.imgj - dBefore.Data.imgj))) + dBefore.Data.imgj;
        TmpFrames(ind).framej = uint64([InterpLUT(8,ind), 1, 1]);
    end
    
    fid = fopen([FolderName filesep 'img_interp.bin'],'w');
    for ind = 1:size(InterpLUT,2)
        fwrite(fid, TmpFrames(ind).framej, 'uint64');
        fwrite(fid, TmpFrames(ind).imgj, 'uint16');
    end
    fclose(fid);
end

%Rebuilding addresses for each frames...
NombreImage = (NombreImage - SkipNFirst);
ImAddressBook = zeros(NombreImage,2);
for ind = 1:NombreImage
    if( ismember(ind, badFrames) )
        fidx = find( ind == InterpLUT(8,:), 1, 'first');
        ImAddressBook(ind,1) = size(imgFilesList,1) + 1;
        ImAddressBook(ind,2) = fidx;
    elseif( ismember(ind, goodFrames) )
        fidx = find( ind == idImg(:,1), 1, 'first');
        ImAddressBook(ind,1) = floor((fidx-1)/256) + 1;
        ImAddressBook(ind,2) = rem(fidx-1, 256) + 1;
    end
end

%Saving infos...
if( ~strcmp(FolderName(end), filesep) )
    FolderName = [FolderName filesep];
end
save([FolderName 'ImagesLUT.mat'], 'ImAddressBook');

%%%%
% Images Classification and filtering
%%%%
if( bFluo )
   
    disp('Fluorescence channel classification:');
   
    tags = nFluo:nbColors:NombreImage;
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        
        if( Binning )
            img = imresize(single(dat.Data.imgj),1/Binning);
        else
            img = dat.Data.imgj;
        end
        Images = single(img);
        fwrite(fidF, Images, 'single');
        
        if( bStim )
            fFluo.Stim(cFluo,1) = single(Stim(indF));
        else
            fFluo.Stim(cFluo,1) = 0;
        end
        cFluo = cFluo + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
          
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            indT = indT + 1;
        end
    end
    ind = ind + 1;
    fFluo.datLength = cFluo - 1;
    fFluo.FirstDim = 'y';
    fclose(fidF);
    disp('done');
end
if( bSpeckle )
   
    disp('Speckle channel classification:');
   
    tags = nSpeckle:nbColors:NombreImage;
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        
        if( Binning )
            img = imresize(single(dat.Data.imgj),1/Binning);
        else
            img = dat.Data.imgj;
        end
        Images = single(img);
        fwrite(fidS, Images, 'single');
        
        if( bStim )
            fSpeckle.Stim(cSpeckle,1) = single(Stim(indF));
        else
            fSpeckle.Stim(cSpeckle,1) = 0;
        end
        cSpeckle = cSpeckle + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
          
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            indT = indT + 1;
        end
    end
    ind = ind + 1;
    fSpeckle.datLength = cSpeckle - 1;
    fSpeckle.FirstDim = 'y';
    fclose(fidS);
    disp('done');
end
if( bRed )
    disp('Red channel classification:');
    
    tags = nRed:nbColors:NombreImage;
    
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        
        if( Binning )
            img = imresize(single(dat.Data.imgj),1/Binning);
        else
            img = dat.Data.imgj;
        end
        Images = single(img);
        fwrite(fidR, Images, 'single');
        
        if( bStim )
            fRed.Stim(cRed,1) = single(Stim(indF));
        else
            fRed.Stim(cRed,1) = 0;
        end
        cRed = cRed + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
           
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
         
            indT = indT + 1;
        end
    end
    
    ind = ind + 1;
    fRed.datLength = cRed - 1;
    fRed.FirstDim = 'y';
    fclose(fidR);

        disp('Done');

end
if( bYellow )
        disp('Yellow channel classification:');

      
    tags = nYellow:nbColors:NombreImage;
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        
        if( Binning )
            img = imresize(single(dat.Data.imgj),1/Binning);
        else
            img = dat.Data.imgj;
        end
        Images = single(img);
        fwrite(fidY, Images, 'single');
        
        if( bStim )
            fYellow.Stim(cYellow,1) = single(Stim(indF));
        else
            fYellow.Stim(cYellow,1) = 0;
        end
        cYellow = cYellow + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
           
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            indT = indT + 1;
        end
    end
    
    ind = ind + 1;
    fYellow.datLength = cYellow - 1;
    fYellow.FirstDim = 'y';
    fclose(fidY);

        disp('Done.');

end
if( bGreen )
   
        disp('Green channel classification:');
   
    
    tags = nGreen:nbColors:NombreImage;
    
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        
        if( Binning )
            img = imresize(single(dat.Data.imgj),1/Binning);
        else
            img = dat.Data.imgj;
        end
        Images = single(img);
        fwrite(fidG, Images, 'single');
        
        if( bStim )
            fGreen.Stim(cGreen,1) = single(Stim(indF));
        else
            fGreen.Stim(cGreen,1) = 0;
        end
        cGreen = cGreen + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
           
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            indT = indT + 1;
        end
    end
    
    ind = ind + 1;
    fGreen.datLength = cGreen - 1;
    fGreen.FirstDim = 'y';
    fclose(fidG);
        disp('done');
 
end

fprintf('\n');
%Verbose

    fprintf(['Done!']);
    fprintf('\n');

%end of Verbose

end
