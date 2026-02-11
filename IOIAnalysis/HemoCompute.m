 function [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize, b_RAMsafeMode)
% HEMOCOMPUTE approximates concentration variation of oxygenated (HbO) and
% de-oxygenated (HbR) hemoglobin from two or three illumination wavelengths 
% of intrinsic signals. 
% Inputs:
%   DataFolder (char): path to folder where the input files "red",
%   "green","yellow" are stored.
%   SaveFolder (char): path where to save the HbO and HbR files. If empty,
%   the data will not be saved to a .dat file.
%   FilterSet (char): type of excitation/emission set of filters used. It is one
%    of the folowing:
%           - gCaMP 
%           - jrGECO
%           - none
%   Illumination (cell): names of illumination colors used ("red", "green",
%    "yellow"). A minimal of two must be provided. 
%   b_normalize (bool): Set to TRUE, to normalize the data (when the input
%   data is not normalized already).
% Outputs:
%   HbO (numerical array): data with the same dimensions of the reflectance
%   channels containing the approximate variations of the oxygenated
%   hemoglobin.
%   HbR (numerical array): data with the same dimensions of the reflectance
%   channels containing the approximate variations of the deoxygenated
%   hemoglobin.

% Inputs Validation
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
    disp('Invalid Filter set name');
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

if ~exist('b_RAMsafeMode','var')
    b_RAMsafeMode = false;
end
% Tags:
fTags = {'fidR', 'fidG', 'fidY'};
cTags = {'iRed', 'iGreen', 'iYellow'};
% colors = {'Red', 'Green', 'Yellow'};
%Files Opening:
NbFrames = inf;
fidR = 0;
fidY = 0;
fidG = 0;
for indC = 1:size(Illumination,2)
    switch lower(Illumination{indC})
        case 'red'
            fidR = fopen([DataFolder 'red.dat']);
            c_r = onCleanup(@() safeFclose(fidR));
            iRed = matfile([DataFolder 'red.mat']);
            NbFrames = min([NbFrames, iRed.datLength(1,end)]);
            NbPix = iRed.datSize;
            Freq = iRed.Freq;
        case 'green'
            fidG = fopen([DataFolder 'green.dat']);
            c_g = onCleanup(@() safeFclose(fidG));
            iGreen = matfile([DataFolder 'green.mat']);
            NbFrames = min([NbFrames, iGreen.datLength(1,end)]);
            NbPix = iGreen.datSize;
            Freq = iGreen.Freq;
        case 'yellow'
            fidY = fopen([DataFolder 'yellow.dat']);
            c_y = onCleanup(@() safeFclose(fidY));
            iYellow = matfile([DataFolder 'yellow.mat']);
            NbFrames = min([NbFrames, iYellow.datLength(1,end)]);
            NbPix = iYellow.datSize;
            Freq = iYellow.Freq;
        otherwise
            disp('Unknown colour');
    end
end
% Get data size:
for i = 1:3
    if exist(cTags{i},'var')
        break
    end
end
eval(['iFile = ' cTags{i} ';']);
datsz = [iFile.datSize, iFile.datLength];
datsz(end) = NbFrames; % Update trial/movie length.    
NbPix = double(NbPix);
% Check if the data is normalized:
indxNorm = [-2 -2 -2];

disp('Checking channel data...')
for i = 1:3
    if eval([fTags{i} '== 0'])
        continue
    end
    % Subset 10 frames to assess if the data was normalized or not.
    frIdx = unique(floor(linspace(1,datsz(3),10)));
    tmp = zeros(datsz(1),datsz(2),'single');          
    for jj = 1:length(frIdx)
        
        eval(['fseek(' fTags{i} ',(frIdx(jj) - 1)*prod(datsz([1 2]))*4,''bof'');']);
        eval(['frame = fread(' fTags{i} ',datsz([1 2]),''*single'');']);        
        tmp = tmp + frame';
    end
    tmp = tmp./length(frIdx);
    Mdat = mean(tmp,'all','omitnan');
    
    clear mapFile
    if Mdat >.75 && Mdat <1.25
        % If the data is centered at one.
        indxNorm(i) = 1;
    elseif Mdat > -.25 && Mdat <.25
        % If the data is centered at zero
        indxNorm(i) = 0;
    else
        % Data not normalized!
        indxNorm(i) = -1;
    end
end
% Check if all channels have the same profile:
if ~all(indxNorm)
    error('The input data is heterogeneous! All channels must be preprocessed in the same way.')
end
indxNorm = indxNorm(1);
if indxNorm == -1 && ~b_normalize    
    error('Operation aborted! The channels must be normalized or set "b_normalize" input to TRUE.')
end
if indxNorm == 0   
    warning('The channels are centered at zero. They will be shifted to be centered at one.')
end
if b_normalize && (indxNorm == 1 || indxNorm == 0)
    b_normalize = false;
    warning('The input data is already normalized. Normalization will be skipped!')
end
disp('Data checked!')
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

% Computation itself:
A = ioi_epsilon_pathlength('Hillman', 100, 60, 40, Filters);

f = fdesign.lowpass('N,F3dB', 4, 1, Freq); %Low Pass
lpass_high = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq); %Low Pass
lpass_low = design(f,'butter');

% Calculate the number of chunks to use
if numel(datsz) == 4    
    % For 4D data with dimensions {'E','Y','X','T}, chunk over the "Events"
    % dimension
    hasEvents = true;
    szYXTE = datsz;
    szYXTE = szYXTE([2 3 4 1]);    
else
    % For 3D data with dimensions {'Y','X','T}:
    hasEvents = false;
    % Calculate number of chunks over X.
    nChunks = calculateMaxChunkSize(prod(datsz)*4,12,.1);          
    chunkX = ceil(NbPix(2) / nChunks);  
end
  
if ~b_RAMsafeMode
    HbO = zeros(datsz, 'single');
    HbR = zeros(datsz, 'single');
else
    metaData = load(iFile.Properties.Source);
    
    metaData.datFile = 'HbO.dat';
    preallocateDatFile(fullfile(SaveFolder,'HbO.dat'),metaData);
    fid_hbo = fopen(fullfile(SaveFolder,'HbO.dat'),'r+');
    c_hbo = onCleanup(@() safeFclose(fid_hbo));
    
    metaData.datFile = 'HbR.dat';    
    preallocateDatFile(fullfile(SaveFolder,'HbR.dat'),metaData);
    fid_hbr = fopen(fullfile(SaveFolder,'HbR.dat'),'r+');
    c_hbr = onCleanup(@() safeFclose(fid_hbr));
end   
    

% Computation loop
h = waitbar(0,'Computing');

for indP = 1:nChunks
    
    % Set indices for chunking
    if ~hasEvents            
        xStart  = (indP-1)*chunkX + 1;
        xEnd    = min(xStart + chunkX -1, NbPix(2));
        xIdx    = xStart:xEnd;
    end
    
    if b_RAMsafeMode
        h.Name = ['HemoCompute (chunk ' num2str(indP) '/' num2str(nChunks) ')'];drawnow()
    end
    if( fidR )
        waitbar(indP/nChunks,h,'Red channel [Reading file...]')
        % Read file
        if hasEvents
            % Read one full trial          
            Red = readTrial(fidR, indP, datsz, 'single');           
        else
            Red = spatialSlabIO('read',fidR,NbPix(1),NbPix(2),NbFrames,xIdx,'single');
        end
        Red = reshape(Red,[],NbFrames);   
        
        if b_normalize 
            waitbar(indP/nChunks,h,'Red channel [Normalizing data...]')
            Red = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Red)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Red)'))';
            tmp(tmp<min(Red(:))) = min(Red(:));
            Red = (Red)./(tmp);
        end        
        if indxNorm == 0 % data centered at zero
            Red = Red + 1;
        end
        Red = -log(Red);
    end
    if( fidG )
        waitbar(indP/nChunks,h,'Green channel [Reading file...]')
        % Read file
        if hasEvents
            % Read one full trial            
            Green = readTrial(fidG, indP, datsz, 'single'); 
        else
            Green = spatialSlabIO('read',fidG,NbPix(1),NbPix(2), NbFrames,xIdx,'single');
        end
        Green = reshape(Green,[],NbFrames);

        if b_normalize
            waitbar(indP/nChunks,h,'Green channel [Normalizing data...]')
            Green = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Green)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Green)'))';
            tmp(tmp<min(Green(:))) = min(Green(:));
            Green = (Green)./(tmp);
        end
        if indxNorm == 0 % data centered at zero
            Green = Green + 1;
        end
        Green = -log(Green);
    end
    if( fidY )        
        waitbar(indP/nChunks,h,'Yellow channel [Reading file...]')
        % Read file
        if hasEvents
            % Read one full trial            
            Yel = readTrial(fidY, indP, datsz, 'single'); 
        else
            Yel = spatialSlabIO('read',fidY,NbPix(1),NbPix(2), NbFrames,xIdx,'single');
        end
        
        Yel = reshape(Yel,[],NbFrames);
        if b_normalize
            waitbar(indP/nChunks,h,'Yellow channel [Normalizing data...]')
            Yel = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Yel)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Yel)'))';
            tmp(tmp<min(Yel(:))) = min(Yel(:));
            Yel = (Yel)./(tmp);
        end
        if indxNorm == 0 % data centered at zero
            Yel = Yel + 1;
        end
        Yel = -log(Yel);
    end
    clear tmp;   
    
    waitbar(indP/nChunks,h,'Computing [HbO] and [HbR]...')
    
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
        
    if numel(datsz) == 4
        Hbs = reshape(Hbs, [2 1, datsz(2:end)]);
        Hbs = real(Hbs);
        if ~b_RAMsafeMode
            HbO(indP,:,:,:) = squeeze(Hbs(1,:,:,:,:));
            HbR(indP,:,:,:) = squeeze(Hbs(2,:,:,:,:));
        end            
                  
    else
        Hbs = reshape(Hbs, 2, NbPix(1), numel(xIdx), []);
        Hbs = real(Hbs);
        if ~b_RAMsafeMode
            HbO(:,xIdx,:) = squeeze(Hbs(1,:,:,:));
            HbR(:,xIdx,:) = squeeze(Hbs(2,:,:,:));
        end
    end
    
    if b_RAMsafeMode
        waitbar(indP/nChunks,h,'Writing to files...')
        if hasEvents
            writeTrial_YXTE(fid_hbo,indP,squeeze(Hbs(1,:,:,:,:)),szYXTE,'single'); % HbO
            writeTrial_YXTE(fid_hbr,indP,squeeze(Hbs(2,:,:,:,:)),szYXTE,'single'); % HbR                 
        else
            spatialSlabIO('write',fid_hbo,NbPix(1),NbPix(2), NbFrames,xIdx,'single',squeeze(Hbs(1,:,:,:))); % HbO
            spatialSlabIO('write',fid_hbr,NbPix(1),NbPix(2), NbFrames,xIdx,'single',squeeze(Hbs(2,:,:,:))); % HbR
        end
    end
    
end

    
close(h);
if b_RAMsafeMode
    fclose(fid_hbr);
    fclose(fid_hbo);
    HbO = 'HbO.dat';
    HbR = 'HbR.dat';
    if hasEvents
        % Permute back to EYXT
        permuteDat_YXTE_to_EYXT_inplace(fullfile(SaveFolder,'HbO.dat'),szYXTE,'single');
        permuteDat_YXTE_to_EYXT_inplace(fullfile(SaveFolder,'HbR.dat'),szYXTE,'single');
    end
    
    return
end
% Save File management:
if bSave
    % Delete existing HbO and HbR files in folder before creating new ones:
    warning('off')
    delete([SaveFolder 'HbO.dat']);
    delete([SaveFolder 'HbO.mat']);
   
    delete([SaveFolder 'HbR.dat']);
    delete([SaveFolder 'HbR.mat']);
    warning('on')
    
    % Save HbO:    
    % Save .DAT file:
    fidHbO = fopen([SaveFolder 'HbO.dat'],'W');
    fwrite(fidHbO, HbO, '*single');
    fclose(fidHbO);    
    fn = setdiff(fieldnames(iFile), {'Properties', 'datFile'});
    
    % Save .MAT file:
    fHbO = matfile([SaveFolder 'HbO.mat'], 'Writable', true);
    fHbO.datFile = 'HbO.dat';
    for i = 1:numel(fn)
        fHbO.(fn{i}) = iFile.(fn{i});
    end
    
    % Save HbR:
    % Save .DAT file:
    fidHbR = fopen([SaveFolder 'HbR.dat'],'W');
    fwrite(fidHbR, HbR, '*single');
    fclose(fidHbR);
    
    % Save .MAT file:
    fHbR = matfile([SaveFolder 'HbR.mat'], 'Writable', true);
    fHbR.datFile = 'HbR.dat';
    for i = 1:numel(fn)
        fHbR.(fn{i}) = iFile.(fn{i});
    end 
end

end