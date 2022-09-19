function out = OpenIOI_OldSyst(FolderName, Binning, OStream)

disp('Computing stimulation parameters')
disp('**************************');
IOIReadStimFile(FolderName);
            
%%%%
%Images memory maping
%%%%
NomSeq = [FolderName filesep 'IOI_scan.seq'];
SizeImage = memmapfile(NomSeq,'Offset',580,'Format','uint32','Repeat',1);
NombreImage = memmapfile(NomSeq,'Offset',572,'Format','uint32','Repeat',1);
ImRes_XY = memmapfile(NomSeq,'Offset',548,'Format','uint32','Repeat',2);
SizeImage = double(SizeImage.Data)/2;
ImRes_XY = double(ImRes_XY.Data);
NombreImage = double(NombreImage.Data);
data = memmapfile(NomSeq,'Offset',1024,'Format',{'uint16', [ImRes_XY(1) ImRes_XY(2)], 'framej';'int32', 1, 'datej';'uint16', 2, 'timej';'uint16', SizeImage-4-ImRes_XY(1)*ImRes_XY(2), 'junkj'},'repeat',inf);

% Verbose
disp(['Opening of: ' FolderName]);
disp(['Number of Frames acquired: ' int2str(NombreImage)]);
disp(['Frames'' resolution: ' int2str(ImRes_XY(1)) ' pix X ' int2str(ImRes_XY(2)) ' pix']);
% end of Verbose

%Empty frame detection at the beginning of the file
FrameToSkip_Start = 0;
while( sum(std(single(data.Data(FrameToSkip_Start + 1).framej),0,2)) == 0 )
    FrameToSkip_Start = FrameToSkip_Start + 1;
end
%Verbose
disp(['Before acquisition false images: ' int2str(FrameToSkip_Start)]);
%end of Verbose

%Empty frame detection at the end of the file
FrameToSkip_End = 0;
while( sum(std(single(data.Data(NombreImage - FrameToSkip_End).framej),0,2)) == 0 )
    FrameToSkip_End = FrameToSkip_End + 1;
end
%Verbose
disp(['After acquisition false images: ' int2str(FrameToSkip_End)]);
%end of Verbose

NombreImage = NombreImage - FrameToSkip_End - FrameToSkip_Start;

%If binning...
if( Binning )
    %Verbose
    disp('Binning option is ON');
    %end of Verbose
end

%%%%
%Stimulation Params
%%%%
load([FolderName filesep 'StimParameters.mat']);
AnalogIN = load([FolderName filesep 'IOI_aux.mat']);
CamTrig = find(diff(AnalogIN.aux(:,2))>10000);
StartDelay = round(CamTrig(1)/10);
EndDelay = round((length(AnalogIN.aux(:,2)) - CamTrig(end))/10);

% Verbose
disp(['Camera Trigs detected: ' int2str(length(CamTrig))]);
disp(['Recording of analog inputs starts ' int2str(StartDelay) ' ms before the first trigger.']) 
disp(['Recording of analog inputs ends ' int2str(EndDelay) ' ms after the last trigger.']) 
% end of Verbose

if( (StartDelay < 100)&&(EndDelay < 100) )
    disp('Both start and end delays of the analog input record are inside the camera trigger period.')
    disp('Analysis of this dataset is not possible.')
    disp('Impossible to set any time refs between illumination and camera acquisition.');
    fprintf('This situation can occurs when the saving of images or/and analog inputs \nare started and terminated while the IOI acquisition is already running.\n');
    fprintf(['The recording sequence should be:\n' ...
        ' - click save button in StreamPix\n'...
        ' - click save button in IOI analog inputs window\n'...
        ' - click start button in IOI software window\n'...
        ' - wait until acquisition is done\n'...
        ' - click stop button in StreamPix\n'...
        ' - click stop button in IOI analog inputs window\n'...
        ' - click ''add experiment'' button in IOI software window\n']);
    out = 'Error';
    return
end

if( length(CamTrig) > NombreImage )
    Timing = zeros(1,NombreImage);
    FrameTime = datevec(double(data.Data(1).datej)/86400 + datenum(1970,1,1));
    Tz = (double(FrameTime(4))*60*60+double(FrameTime(5))*60+double(FrameTime(6)))*1000;
    Tz = Tz + double(data.Data(1).timej(1)) + double(data.Data(1).timej(2))/1000.0;
    for indF = 1:NombreImage
        FrameTime = datevec(double(data.Data(indF).datej)/86400 + datenum(1970,1,1));
        Timing(indF) = (double(FrameTime(4))*60*60 + double(FrameTime(5))*60+double(FrameTime(6)))*1000;
        Timing(indF) = Timing(indF) + double(data.Data(indF).timej(1)) + double(data.Data(indF).timej(2))/1000.0 - Tz;
    end
    
    IdxSaut = find(diff(Timing) > mean(diff(Timing))+0.5*mean(diff(Timing)));
    LtSaut = round((Timing(IdxSaut+1) - Timing(IdxSaut))/mean(diff(Timing)))-1;
    NoFrameSkip = (IdxSaut+1):(IdxSaut+LtSaut);
    NombreImage = NombreImage + LtSaut;
    
   disp(['Frames missing in seq file: ' int2str(NoFrameSkip) ' at position ' int2str(IdxSaut)]);
end

if( length(CamTrig) < (NombreImage - FrameToSkip_End - FrameToSkip_Start)  ) 
    disp('IOI Error: Camera was in "Free Running" mode. Acquisition is then not usable. Please set the camera to "IO Trigger"');
    out = 'Error';
    return
end

%%%%
% Color Sequence
%%%%
FrameSeq = load([FolderName filesep 'IOI_scaninfo.mat']);
Tmp = abs(fft(FrameSeq.Signaux,[],1));
[~, Freq] = max(Tmp(2:100,:),[],1);

ColorsString = {};
NbColors = sum(Freq(2:5)>1);
rPos = find(FrameSeq.Signaux(:,2),1,'first');
gPos = find(FrameSeq.Signaux(:,3),1,'first');
yPos = find(FrameSeq.Signaux(:,4),1,'first');
sPos = find(FrameSeq.Signaux(:,5),1,'first');
Seq = sort([rPos, gPos, yPos, sPos]);
if( Binning )
    Rx = round(ImRes_XY(1)/2);
    Ry = round(ImRes_XY(2)/2);
else
    Rx = ImRes_XY(1);
    Ry = ImRes_XY(2);
end
if(Freq(2) > 1)
    ColorsString{end+1} = 'Red';
    disp('Red illumination detected');
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(1, 2, 'single');
    fRed.Freq = Freq(2);
    cRed = 1;
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
end
if(Freq(3) > 1)
    ColorsString{end+1} =  'Green';
    disp('Green illumination detected');
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(1, 2, 'single');
    fGreen.Freq = Freq(3);
    cGreen = 1;
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
end
if(Freq(4) > 1)
    ColorsString{end+1} =  'Yellow';
    disp('Yellow illumination detected');
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(1, 2, 'single');
    fYellow.Freq = Freq(4);
    cYellow = 1;
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
end
if(Freq(5) > 1)
    ColorsString{end+1} =  'Speckle';
    disp('Speckle illumination detected');
    if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
        delete([FolderName filesep 'Data_speckle.mat']);
    end
    fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);
    fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
    fSpeckle.datSize = [Rx, Ry];
    fSpeckle.Stim = zeros(1, 2, 'single');
    fSpeckle.Freq = Freq(5);
    cSpeckle = 1;
    fidS = fopen([FolderName filesep 'sChan.dat'],'w');
end

% Verbose
if( sum(Stim(:)) > 0 )
    disp('Stim detected: yes');
    disp(['Number of events: ' int2str(NbStim)]);
    disp(['Length of each event: ' int2str(StimLength) 'sec']);
end
% end of Verbose

%%%%
% Images Classification and filtering
%%%%
NbCycles = floor(NombreImage/size(ColorsString,2));
Marks = round(linspace(0, NbCycles, 10));
Marks = Marks*size(ColorsString,2);
fprintf('Progress: ');
data = memmapfile(NomSeq,'Offset',1024 + 2*FrameToSkip_Start*SizeImage,...
    'Format',{'uint16', [ImRes_XY(1) ImRes_XY(2)], 'framej';...
    'int32', 1, 'datej';'uint16', 2, 'timej';...
    'uint16', SizeImage-4-ImRes_XY(1)*ImRes_XY(2), 'junkj'},...
    'repeat', inf);
for ind = 2:length(Marks)
    
    Frame = data.Data((Marks(ind-1)+1):Marks(ind));
    Frame = reshape([Frame(:).framej],ImRes_XY(1),ImRes_XY(2),[]);
        
    if( Binning )
        Frame = imresize(single(Frame), 0.5);
    else
        Frame = single(Frame);
    end
    
    tColor = blanks(size(Frame,3));
    for indF = 1:size(Frame,3)
        cTag = ColorsString{mod(indF-1,size(ColorsString,2)) + 1};
        tColor(indF) = cTag(1);
    end
    
    for indC = 1:size(ColorsString,2)
        cTag = ColorsString{indC};
        idx = strfind(tColor, cTag(1));
        nbI = length(idx);
        eval(['fwrite(fid' cTag(1) ', Frame(:,:, idx),''single'');']);
        eval(['f' cTag '.Stim(1, c' cTag ':(c' cTag ' + nbI - 1)) = Stim(idx + Marks(ind-1))'';']);
        eval(['c' cTag ' = c' cTag ' + nbI;']);
    end
    fprintf('%d%%...', 11*(ind-1));
end
for indC = 1:size(ColorsString,2)
    cTag = ColorsString{indC};
    eval(['f' cTag '.datLength = c' cTag ';']);
end
fprintf('\n');

%Verbose
disp(['Done with file ' FolderName]);
str = ['************* ' sprintf('\r')];
disp(str);
%end of Verbose
end
