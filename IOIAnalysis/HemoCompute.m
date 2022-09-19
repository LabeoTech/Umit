function [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination)

%Inputs Validation
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = strcat(DataFolder, filesep);
end
bSave = true;
if isempty(SaveFolder)
    bSave = false; % Disable saving if SaveFolder was not provided.
elseif( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = strcat(SaveFolder, filesep);
end
    
if( ~contains(lower(FilterSet), {'gcamp', 'jrgeco', 'none'}) )
    disp('Invalide Filter set name');
    return;
end
idx = contains(lower(Illumination), 'amber');
if( any(idx) )
    Illumination{idx} = 'yellow';
end
if( sum(contains({'red', 'yellow', 'green'}, lower(Illumination))) < 2 )
    disp('At least two different illumination wavelengths are needed for Hb computation'); 
    return;
end

%Files Opening:
NbFrames = inf;
fidR = 0;
fidY = 0;
fidG = 0;
for indC = 1:size(Illumination,2)
    switch lower(Illumination{indC})
        case 'red'
            fidR = fopen([DataFolder 'red.dat']);
            iRed = matfile([DataFolder 'red.mat']);
            NbFrames = min([NbFrames, iRed.datLength]);
            NbPix = iRed.datSize;
            Freq = iRed.Freq;
        case 'green'
            fidG = fopen([DataFolder 'green.dat']);
            iGreen = matfile([DataFolder 'green.mat']);
            NbFrames = min([NbFrames, iGreen.datLength]);
            NbPix = iGreen.datSize;
            Freq = iGreen.Freq;
        case 'yellow'
            fidY = fopen([DataFolder 'yellow.dat']);
            iYellow = matfile([DataFolder 'yellow.mat']);
            NbFrames = min([NbFrames, iYellow.datLength]);
            NbPix = iYellow.datSize;
            Freq = iYellow.Freq;
        otherwise
            disp('Unknown colour');
    end
end
NbPix = double(NbPix);
% clear iRed iGreen iYellow indC;


% Filter setting
switch( lower(FilterSet) )
    case 'gcamp'
        Filters.Excitation = 'GCaMP';
        Filters.Emission = 'GCaMP';
    case 'jrgeco'
        Filters.Excitation = 'none';
        Filters.Emission = 'jRGECO';
    otherwise
        Filters.Excitation = 'none';
        Filters.Emission = 'none';
end
Infos = load([DataFolder 'AcqInfos.mat']);
Filters.Camera = Infos.AcqInfoStream.Camera_Model;
clear Infos;

%Computation itself:
A = ioi_epsilon_pathlength('Hillman', 100, 60, 40, Filters);

MemFact = 16;
f = fdesign.lowpass('N,F3dB', 4, 1, Freq); %Low Pass
lpass_high = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq); %Low Pass
lpass_low = design(f,'butter');
NbPts = floor(NbFrames/100);
Precision = [int2str(NbPix(1)*MemFact) '*single'];
HbO = zeros(NbPix(1), NbPix(2), NbFrames, 'single');
HbR = zeros(NbPix(1), NbPix(2), NbFrames, 'single');

% Computation loop
h = waitbar(0,'Computing');
nIter = NbPix(2)/MemFact;
for indP = 1:nIter
    if( fidR )
        fseek(fidR, (indP-1)*NbPix(1)*MemFact*4,'bof');
        Red = fread(fidR,[NbPix(1)*MemFact, NbFrames],Precision,(NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4);
        Red = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Red)'))';
        tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Red)'))';
        tmp(tmp<min(Red(:))) = min(Red(:));
        Red = (Red)./(tmp);
        Red = -log10(Red);
    end
    if( fidG )
        fseek(fidG, (indP-1)*NbPix(1)*MemFact*4,'bof');
        Green = fread(fidG,[NbPix(1)*MemFact, NbFrames],Precision,(NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4);
        Green = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Green)'))';
        tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Green)'))';
        tmp(tmp<min(Green(:))) = min(Green(:));
        Green = (Green)./(tmp);
        Green = -log10(Green);
    end
    if( fidY )
        fseek(fidY, (indP-1)*NbPix(1)*MemFact*4,'bof');
        Yel = fread(fidY,[NbPix(1)*MemFact, NbFrames],Precision,(NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4);
        Yel = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Yel)'))';
        tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Yel)'))';
        tmp(tmp<min(Yel(:))) = min(Yel(:));
        Yel = (Yel)./(tmp);
        Yel = -log10(Yel);
    end
    clear tmp;
    if(  fidR*fidG*fidY > 0)
        Ainv = pinv(A);
        Hbs = Ainv*([Red(:), Green(:), Yel(:)]') .* 1e6;
        clear Red Green Yel;
    elseif( fidR*fidG > 0)
        Ainv = pinv(A(1:2,:));
        Hbs = Ainv*([Red(:), Green(:)]') .* 1e6;
        clear Red Green;
    elseif( fidG*fidY > 0)
        Ainv = pinv(A(2:3,:));
        Hbs = Ainv*([Green(:), Yel(:)]') .* 1e6;
        clear Green Yel;
    else
        Ainv = pinv(A([1 3],:));
        Hbs = Ainv*([Red(:), Yel(:)]') .* 1e6;
        clear Red Yel;
    end
    
    Hbs = reshape(Hbs, 2, NbPix(1), MemFact, []);
    Hbs = real(Hbs);
    HbO(:,(indP-1)*MemFact + (1:MemFact),:) = squeeze(Hbs(1,:,:,:));
    HbR(:,(indP-1)*MemFact + (1:MemFact),:) = squeeze(Hbs(2,:,:,:));
    waitbar(indP/nIter,h);
end
close(h);
% Save File management:
if( bSave )
    % Save HbO:
    % Save .DAT file:
    fidHbO = fopen([SaveFolder 'HbO.dat'],'W');
    fwrite(fidHbO, HbO, '*single');
    fclose(fidHbO);
    % Get handles from one of the color channels (data from
    % ImagesClassification):
    cTags = {'iRed', 'iGreen', 'iYellow'};
    for i = 1:3
        if exist(cTags{i},'var')
            break
        end
    end
    % Save .MAT file:
    fHbO = matfile([SaveFolder 'HbO.mat'], 'Writable', true);
    fHbO.datFile = 'HbO.dat';
    fHbO.datSize = [size(HbO,1), size(HbO,2)];    
    fHbO.datLength = size(HbO,3);
    fHbO.FirstDim = 'y';
    fHbO.Datatype = 'single';
    fHbO.datName = 'data';
    fHbO.dim_names = {'Y', 'X', 'T'};
    fHbO.Freq = Freq;
    % Copy "Stim" and "tExposure" from input color channel:
    eval(['fHbO.Stim = ' cTags{i} '.Stim;']); 
    eval(['fHbO.tExposure = ' cTags{i} '.tExposure;']); 
    
    % Save HbR:
    fidHbR = fopen([SaveFolder 'HbR.dat'],'W');    
    fwrite(fidHbR, HbR, '*single');
    fclose(fidHbR);
    
     % Save .MAT file:
    fHbR = matfile([SaveFolder 'HbR.mat'], 'Writable', true);
    fHbR.datFile = 'HbR.dat';
    fHbR.datSize = [size(HbR,1), size(HbR,2)];    
    fHbR.datLength = size(HbR,3);
    fHbR.FirstDim = 'y';
    fHbR.Datatype = 'single';
    fHbR.datName = 'data';
    fHbR.dim_names = {'Y', 'X', 'T'};
    fHbR.Freq = Freq;
    % Copy "Stim" and "tExposure" from one of the input color channels (data from ImagesClassification):
    eval(['fHbR.Stim = ' cTags{i} '.Stim;']); 
    eval(['fHbR.tExposure = ' cTags{i} '.tExposure;']); 
end
