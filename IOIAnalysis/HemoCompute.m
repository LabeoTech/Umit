function [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, Illumination, b_normalize, b_RAMsafeMode)
% HEMOCOMPUTE Approximate HbO/HbR from two or three reflectance channels.
%
%   [HbO, HbR] = HemoCompute(DataFolder, SaveFolder, FilterSet, ...
%       Illumination, b_normalize, b_RAMsafeMode)
%
%   Estimates concentration variation of oxygenated (HbO) and
%   deoxygenated (HbR) hemoglobin from two or three intrinsic-signal
%   illumination wavelengths. Supported channels are red, green, and yellow
%   (amber is mapped to yellow).
%
%   If selected channels have different temporal sampling rates or lengths,
%   higher-frequency channels are anti-aliased and resampled to match the
%   lowest-frequency selected channel. Same-frequency channels with small
%   length differences are resampled without anti-aliasing. Lower-frequency
%   channels are not upsampled.
%
%   Inputs:
%       DataFolder      - Folder containing input .dat/.mat channel files.
%       SaveFolder      - Folder where HbO/HbR files are saved. If empty,
%                         files are not saved.
%       FilterSet       - 'gCaMP', 'jrGECO', or 'none'.
%       Illumination    - Cell array with at least two of: red, green,
%                         yellow/amber.
%       b_normalize     - If true, normalize channels before Hb computation.
%       b_RAMsafeMode   - If true, stream spatial chunks to disk.
%
%   Outputs:
%       HbO             - HbO array, or 'HbO.dat' in RAM-safe mode.
%       HbR             - HbR array, or 'HbR.dat' in RAM-safe mode.

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

% Tags and native-channel metadata.
fTags = {'fidR', 'fidG', 'fidY'};
cTags = {'iRed', 'iGreen', 'iYellow'};
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
            c_r = onCleanup(@() safeFclose(fidR)); %#ok<NASGU>
            iRed = load([DataFolder 'red.mat']);
            nativeLength(1) = iRed.datLength(1,end);
            nativeFreq(1) = iRed.Freq;
            nativeDatSize{1} = iRed.datSize;
            selectedChannels(1) = true;
        case 'green'
            fidG = fopen([DataFolder 'green.dat']);
            c_g = onCleanup(@() safeFclose(fidG)); %#ok<NASGU>
            iGreen = load([DataFolder 'green.mat']);
            nativeLength(2) = iGreen.datLength(1,end);
            nativeFreq(2) = iGreen.Freq;
            nativeDatSize{2} = iGreen.datSize;
            selectedChannels(2) = true;
        case 'yellow'
            fidY = fopen([DataFolder 'yellow.dat']);
            c_y = onCleanup(@() safeFclose(fidY)); %#ok<NASGU>
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

% Validate spatial/event dimensions and choose the lowest-frequency timeline.
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

% Check if the data is normalized.
indxNorm = [-2 -2 -2];

disp('Checking channel data...')
for i = 1:3
    if eval([fTags{i} '== 0'])
        continue
    end

    % Subset 10 frames from each channel native timeline to assess whether
    % the data were normalized.
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

    clear mapFile
    if Mdat >.75 && Mdat <1.25
        % If the data is centered at one.
        indxNorm(i) = 1;
    elseif Mdat > -.25 && Mdat <.25
        % If the data is centered at zero.
        indxNorm(i) = 0;
    else
        % Data not normalized.
        indxNorm(i) = -1;
    end
end

% Check if all selected channels have the same preprocessing profile.
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

% Filter setting.
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

% Computation itself.
A = ioi_epsilon_pathlength('Hillman', 100, 60, 40, Filters);

f = fdesign.lowpass('N,F3dB', 4, 1, Freq); % Low Pass
lpass_high = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq); % Low Pass
lpass_low = design(f,'butter');

% Calculate the number of chunks to use.
if numel(datsz) == 4
    % For 4D data with dimensions {'E','Y','X','T}, chunk over the Events
    % dimension.
    hasEvents = true;
    nChunks = datsz(1);
    szYXTE = datsz;
    szYXTE = szYXTE([2 3 4 1]);
else
    % For 3D data with dimensions {'Y','X','T}, chunk over X.
    hasEvents = false;
    nChunks = calculateMaxChunkSize(prod(datsz)*4,12,.1);
    chunkX = ceil(NbPix(2) / nChunks);
end

if ~b_RAMsafeMode
    HbO = zeros(datsz, 'single');
    HbR = zeros(datsz, 'single');
else
    metaData = load(iFile.Properties.Source);
    metaData.datLength = NbFrames;
    metaData.Freq = Freq;

    metaData.datFile = 'HbO.dat';
    preallocateDatFile(fullfile(SaveFolder,'HbO.dat'),metaData);
    fid_hbo = fopen(fullfile(SaveFolder,'HbO.dat'),'r+');
    c_hbo = onCleanup(@() safeFclose(fid_hbo)); 

    metaData.datFile = 'HbR.dat';
    preallocateDatFile(fullfile(SaveFolder,'HbR.dat'),metaData);
    fid_hbr = fopen(fullfile(SaveFolder,'HbR.dat'),'r+');
    c_hbr = onCleanup(@() safeFclose(fid_hbr)); 
end

% Computation loop.
h = waitbar(0,'Computing');

for indP = 1:nChunks

    % Set indices for chunking.
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
        if hasEvents
            redDatsz = [nativeDatSize{1}, nativeLength(1)];
            Red = readTrial(fidR, indP, redDatsz, 'single');
        else
            Red = spatialSlabIO('read',fidR,NbPix(1),NbPix(2),nativeLength(1),xIdx,'single');
        end
        Red = iResampleChannelToLowestFrequency(Red, iRed, NbFrames, Freq);
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
        if hasEvents
            greenDatsz = [nativeDatSize{2}, nativeLength(2)];
            Green = readTrial(fidG, indP, greenDatsz, 'single');
        else
            Green = spatialSlabIO('read',fidG,NbPix(1),NbPix(2),nativeLength(2),xIdx,'single');
        end
        Green = iResampleChannelToLowestFrequency(Green, iGreen, NbFrames, Freq);
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
        if hasEvents
            yellowDatsz = [nativeDatSize{3}, nativeLength(3)];
            Yel = readTrial(fidY, indP, yellowDatsz, 'single');
        else
            Yel = spatialSlabIO('read',fidY,NbPix(1),NbPix(2),nativeLength(3),xIdx,'single');
        end
        Yel = iResampleChannelToLowestFrequency(Yel, iYellow, NbFrames, Freq);
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
            spatialSlabIO('write',fid_hbo,NbPix(1),NbPix(2),NbFrames,xIdx,'single',squeeze(Hbs(1,:,:,:))); % HbO
            spatialSlabIO('write',fid_hbr,NbPix(1),NbPix(2),NbFrames,xIdx,'single',squeeze(Hbs(2,:,:,:))); % HbR
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
        % Permute back to EYXT.
        permuteDat_YXTE_to_EYXT_inplace(fullfile(SaveFolder,'HbO.dat'),szYXTE,'single');
        permuteDat_YXTE_to_EYXT_inplace(fullfile(SaveFolder,'HbR.dat'),szYXTE,'single');
    end
    return
end

% Save File management.
if bSave
    % Delete existing HbO and HbR files in folder before creating new ones.
    warning('off')
    delete([SaveFolder 'HbO.dat']);
    delete([SaveFolder 'HbO.mat']);

    delete([SaveFolder 'HbR.dat']);
    delete([SaveFolder 'HbR.mat']);
    warning('on')

    % Save HbO .dat file.
    fidHbO = fopen([SaveFolder 'HbO.dat'],'W');
    fwrite(fidHbO, HbO, '*single');
    fclose(fidHbO);
    fn = setdiff(fieldnames(iFile), {'Properties', 'datFile'});

    % Save HbO .mat file.
    fHbO = matfile([SaveFolder 'HbO.mat'], 'Writable', true);
    fHbO.datFile = 'HbO.dat';
    for i = 1:numel(fn)
        fHbO.(fn{i}) = iFile.(fn{i});
    end
    fHbO.datLength = NbFrames;
    fHbO.Freq = Freq;

    % Save HbR .dat file.
    fidHbR = fopen([SaveFolder 'HbR.dat'],'W');
    fwrite(fidHbR, HbR, '*single');
    fclose(fidHbR);

    % Save HbR .mat file.
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
