function [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize,varargin)
% HEMOCOMPUTE approximates concentration variation of oxygenated (HbO) and
% de-oxygenated (HbR) hemoglobin from two or three illumination wavelengths 
% of intrinsic signals. 
% Inputs:
%   DataFolder (char): path to folder where the input files "red",
%   "green", "yellow" are stored.
%   SaveFolder (char): path where to save the HbO and HbR files. If empty,
%   the data will not be saved to a .dat file.
%   FilterSet (char): type of excitation/emission set of filters used. It is one
%    of the following:
%           - gCaMP 
%           - jrGECO
%           - none
%   Illumination (cell): names of illumination colors used ("red", "green",
%    "yellow"). A minimum of two must be provided. 
%   b_normalize (bool): Set to TRUE to normalize the data (when the input
%   data is not normalized already).
%
%   HbT_uM (optional, numerical): Total hemoglobin concentration in micromolar. 
%   Default is 100 uM if not provided.
%
%   O2_sat (optional, numerical): Oxygen saturation percentage. Default is 60% 
%   if not provided.
%
% Outputs:
%   HbO (numerical array): data with the same dimensions as the reflectance
%   channels containing the approximate variations of the oxygenated
%   hemoglobin.
%   HbR (numerical array): data with the same dimensions as the reflectance
%   channels containing the approximate variations of the deoxygenated
%   hemoglobin.



% Check for extra inputs:
if nargin == 6
    HbT_uM = varargin{1};
elseif nargin == 7
    HbT_uM = varargin{1};
    O2_sat = varargin{2};
else
    HbT_uM = 100;
    O2_sat = 60;
end



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

% if( ~contains(lower(FilterSet), {'gcamp', 'jrgeco', 'none'}) )
%     disp('Invalid Filter set name');
%     return;
% end
idx = contains(lower(Illumination), 'amber');
if( any(idx) )
    Illumination{idx} = 'yellow';
end
if( sum(contains({'red', 'yellow', 'green'}, lower(Illumination))) < 2 )
    disp('At least two different illumination wavelengths are needed for Hb computation');
    return;
end
% Tags and native-channel metadata:
fTags = {'fidR', 'fidG', 'fidY'};
channelNames = {'red', 'green', 'yellow'};
nativeLength = nan(1,3);
nativeFreq = nan(1,3);
nativeDatSize = cell(1,3);
selectedChannels = false(1,3);
% Files Opening:
fidR = 0;
fidY = 0;
fidG = 0;
for indC = 1:size(Illumination,2)
    switch lower(Illumination{indC})
        case 'red'
            fidR = fopen([DataFolder 'red.dat']);
            iRed = load([DataFolder 'red.mat']);
            nativeLength(1) = iRed.datLength(1,end);
            nativeFreq(1) = iRed.Freq;
            nativeDatSize{1} = iRed.datSize;
            selectedChannels(1) = true;
        case 'green'
            fidG = fopen([DataFolder 'green.dat']);
            iGreen = load([DataFolder 'green.mat']);
            nativeLength(2) = iGreen.datLength(1,end);
            nativeFreq(2) = iGreen.Freq;
            nativeDatSize{2} = iGreen.datSize;
            selectedChannels(2) = true;
        case 'yellow'
            fidY = fopen([DataFolder 'yellow.dat']);
            iYellow = load([DataFolder 'yellow.mat']);
            nativeLength(3) = iYellow.datLength(1,end);
            nativeFreq(3) = iYellow.Freq;
            nativeDatSize{3} = iYellow.datSize;
            selectedChannels(3) = true;
        otherwise
            disp('Unknown colour');
    end
end

if ~any(selectedChannels)
    error('No valid illumination channels were selected.');
end

% Validate dimensions and choose the lowest-frequency timeline.
idxSelected = find(selectedChannels);
targetFreq = min(nativeFreq(idxSelected));
freqTol = max(1e-6, abs(targetFreq) * 1e-6);
idxTargetCandidates = idxSelected(abs(nativeFreq(idxSelected) - targetFreq) <= freqTol);
[~, idxMinLength] = min(nativeLength(idxTargetCandidates));
idxTarget = idxTargetCandidates(idxMinLength);
NbFrames = nativeLength(idxTarget);
Freq = nativeFreq(idxTarget);
NbPix = nativeDatSize{idxTarget};

for i = idxSelected
    if ~isequal(nativeDatSize{i}, NbPix)
        error(['All selected illumination channels must have matching spatial/event ' ...
            'dimensions before Hb computation. Channel "%s" does not match the ' ...
            'lowest-frequency reference channel.'], channelNames{i});
    end
end

% Get target metadata object.
switch idxTarget
    case 1
        iFile = iRed;
    case 2
        iFile = iGreen;
    case 3
        iFile = iYellow;
end

datsz = [iFile.datSize, NbFrames];
NbPix = double(NbPix);
% Check if the data is normalized:
indxNorm = [-2 -2 -2];

disp('Checking channel data...')
for i = 1:3
    if eval([fTags{i} '== 0'])
        continue
    end
    thisDatsz = [nativeDatSize{i}, nativeLength(i)];
    frIdx = unique(floor(linspace(1,thisDatsz(end),10)));
    tmp = zeros(thisDatsz(1),thisDatsz(2),'single');
    for jj = 1:length(frIdx)
        eval(['fseek(' fTags{i} ',(frIdx(jj) - 1)*prod(thisDatsz([1 2]))*4,''bof'');']);
        eval(['frame = fread(' fTags{i} ',thisDatsz([1 2]),''*single'');']);
        tmp = tmp + frame';
    end
    tmp = tmp./length(frIdx);
    Mdat = mean(tmp,'all','omitnan');

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
% Check if all selected channels have the same profile:
selectedNorm = indxNorm(indxNorm ~= -2);
if isempty(selectedNorm) || any(selectedNorm ~= selectedNorm(1))
    error('The input data is heterogeneous! All channels must be preprocessed in the same way.')
end
indxNorm = selectedNorm(1);
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
Infos = load([DataFolder 'AcqInfos.mat']);
%Computation itself:
A = ioi_epsilon_pathlength('Hillman', HbT_uM, O2_sat, 100 - O2_sat, FilterSet, Infos.AcqInfoStream.Camera_Model);
clear Infos;
f = fdesign.lowpass('N,F3dB', 4, 1, Freq); %Low Pass
lpass_high = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq); %Low Pass
lpass_low = design(f,'butter');
% NbPts = floor(NbFrames/100);
if numel(datsz) == 4  
    % For 4D data with dimensions {'E','Y','X','T}:
    nIter = double(datsz(1));
    offset = 4;
    Size = [prod(datsz(2:3)), datsz(4)];
    Precision = '*single';
    Skip = (datsz(1)-1)*4;
else
    % For 3D data with dimensions {'Y','X','T}:
    MemFact = 16;
    Precision = [int2str(NbPix(1)*MemFact) '*single'];    
    Skip = (NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4;
    nIter = NbPix(2)/MemFact;
    offset = NbPix(1)*MemFact*4;
    Size = [NbPix(1)*MemFact, NbFrames];
end
HbO = zeros(datsz, 'single');
HbR = zeros(datsz, 'single');

% Computation loop
h = waitbar(0,'Computing');
for indP = 1:nIter
    if( fidR )
%         fseek(fidR, (indP-1)*NbPix(1)*MemFact*4,'bof');
%         Red = fread(fidR,[NbPix(1)*MemFact, NbFrames],Precision,(NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4);
        fseek(fidR,(indP-1)*offset,'bof');
        redSize = Size;
        redSize(end) = nativeLength(1);
        Red = fread(fidR,redSize,Precision,Skip);
        Red = iResampleChannelToLowestFrequency(Red, iRed, NbFrames, Freq);
        if b_normalize 
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
%         fseek(fidG, (indP-1)*NbPix(1)*MemFact*4,'bof');
%         Green = fread(fidG,[NbPix(1)*MemFact, NbFrames],Precision,(NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4);
        fseek(fidG,(indP-1)*offset,'bof');
        greenSize = Size;
        greenSize(end) = nativeLength(2);
        Green = fread(fidG,greenSize,Precision,Skip);
        Green = iResampleChannelToLowestFrequency(Green, iGreen, NbFrames, Freq);
        if b_normalize
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
%         fseek(fidY, (indP-1)*NbPix(1)*MemFact*4,'bof');
%         Yel = fread(fidY,[NbPix(1)*MemFact, NbFrames],Precision,(NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4);
        fseek(fidY,(indP-1)*offset,'bof');
        yellowSize = Size;
        yellowSize(end) = nativeLength(3);
        Yel = fread(fidY,yellowSize,Precision,Skip);
        Yel = iResampleChannelToLowestFrequency(Yel, iYellow, NbFrames, Freq);
        if b_normalize
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
        HbO(indP,:,:,:) = squeeze(Hbs(1,:,:,:,:));
        HbR(indP,:,:,:) = squeeze(Hbs(2,:,:,:,:));
    else
        Hbs = reshape(Hbs, 2, NbPix(1), MemFact, []);
        Hbs = real(Hbs);
        HbO(:,(indP-1)*MemFact + (1:MemFact),:) = squeeze(Hbs(1,:,:,:));
        HbR(:,(indP-1)*MemFact + (1:MemFact),:) = squeeze(Hbs(2,:,:,:));
    end
    
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
    fn = setdiff(fieldnames(iFile), {'Properties', 'datFile'});
    
    % Save .MAT file:
    fHbO = matfile([SaveFolder 'HbO.mat'], 'Writable', true);
    fHbO.datFile = 'HbO.dat';
    for i = 1:numel(fn)
        fHbO.(fn{i}) = iFile.(fn{i});
    end
    fHbO.datLength = NbFrames;
    fHbO.Freq = Freq;
    
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
    fHbR.datLength = NbFrames;
    fHbR.Freq = Freq;
end

end

function data = iResampleChannelToLowestFrequency(data, sourceMeta, targetNt, targetFreq)
%IRESAMPLECHANNELTOLOWESTFREQUENCY Match one channel to the target timeline.
%
%   data = iResampleChannelToLowestFrequency(data, sourceMeta, targetNt,
%   targetFreq) resamples data along the last dimension so it has targetNt
%   samples. If the source frequency is higher than the target frequency, an
%   anti-aliasing low-pass filter is applied before interpolation. If the
%   source frequency matches the target but the lengths differ, interpolation
%   is applied without anti-aliasing. Lower-frequency sources are not
%   upsampled.

sourceNt = sourceMeta.datLength(1,end);
sourceFreq = sourceMeta.Freq;

if sourceNt == targetNt
    return
end

freqTol = max(1e-6, abs(targetFreq) * 1e-6);
inputSize = size(data);
if inputSize(end) ~= sourceNt
    error(['Temporal resampling failed. The input data length does not match ' ...
        'the source metadata datLength.']);
end

data = reshape(data, [], sourceNt);

if sourceFreq > targetFreq + freqTol
    aaCutoff = 0.45 * targetFreq;
    if aaCutoff > 0 && aaCutoff < sourceFreq/2
        f = fdesign.lowpass('N,F3dB', 4, aaCutoff, sourceFreq);
        lpass = design(f,'butter');
        data = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(data)'))';
    end
elseif sourceFreq < targetFreq - freqTol
    error(['Cannot upsample a lower-frequency channel for Hb computation. ' ...
        'Select channels with compatible frequencies or downsample all ' ...
        'channels to the lowest-frequency timeline.']);
end

xSource = linspace(0, 1, sourceNt);
xTarget = linspace(0, 1, targetNt);
data = interp1(xSource, single(data)', xTarget, 'linear', 'extrap')';

inputSize(end) = targetNt;
data = reshape(single(data), inputSize);

end
