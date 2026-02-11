function varargout = NormalisationFiltering(FolderData, FileData, lowFreq,...
    highFreq, bDivide, bExpfit, varargin)
%%%%% Data Normalisation by low-pass filtering %%%%%
% 
% General Infos:
%
% This function can be used to normalise channels (delta F/F or delta R/R),
% or to do a low-pass filtering.
%
% Inputs:
%
% Option A: data to be normalised will be opened within this function
%
%   1- FolderData:  Folder containing the data to be oppened
%   2- FileData:    Channel to open('red', 'green', 'fluo_475', etc.)
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore
%   5- bDivide:     if 1, the data returned (below highFreq) is normalised
%                   by the low freq signal (below lowFreq)
%                   if 0, the low freq signal (below lowFreq) is
%                   substracted from the data returned (below highFreq)
%   6- bExpFit:     if 1, a double exponential curve is fit on each pixels
%                   to correct for the illumination decay that can occurs when temperature
%                   is not well managed
%
%   Ex: Dat = NormalisationFiltering(pwd, 'red', 1/120, 1, 1);
%   
%   This call would return the red channel with a low-pass at 1 Hz,
%   normalized (through a division) by a low-pass at 1/120 Hz.
%
% Option B: data to be normalised is given by one of the argument
%
%   1- FolderData:  Folder containing the data to be oppened
%   2- FileData:    Data, as a 3D matrix (Y, X, Time).
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore 
%   5- bDivide:     if 1, the data returned (below highFreq) is normalised
%                   by the low freq signal (below lowFreq)
%                   if 0, the low freq signal (below lowFreq) is
%                   substracted from the data returned (below highFreq)
%   6- bExpFit:     if 1, a double exponential curve is fit on each pixels
%                   to correct for the illumination decay that can occurs when temperature
%                   is not well managed
%   7- Freq (Optional): Data sample rate.
%
%   Ex: Dat = NormalisationFiltering(pwd, dat, 0, 1, 1);
%   
%   This call would return the data contained in the variable dat with a
%   low-pass at 1 Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'FolderData',@isfolder)
addRequired(p,'FileData', @(x) ischar(x) || isnumeric(x))
addRequired(p,'lowFreq',@(x) isscalar(x) & isnumeric(x))
addRequired(p,'highFreq', @(x) isscalar(x) & isnumeric(x))
addRequired(p,'bDivide', @(x) isscalar(x) & (isnumeric(x) | islogical(x)))
addRequired(p,'bExpFit', @(x) isscalar(x) & (isnumeric(x) | islogical(x)))
addOptional(p,'Freq', [], @(x) (isscalar(x) & isnumeric(x)) | isempty(x))
addOptional(p,'saveFilename','',@ischar);

% Parse inputs:
parse(p, FolderData, FileData, lowFreq, highFreq, bDivide, bExpfit, varargin{:});
% Select function based on "FileData" variable type:
if( isa(p.Results.FileData, 'char') )
    if nargout
        OutData = NormFiltFromFile(p.Results.FolderData, p.Results.FileData,...
            p.Results.lowFreq, p.Results.highFreq, p.Results.bDivide, p.Results.bExpFit,true,p.Results.saveFilename);
        varargout{1} = OutData;
    else
        NormFiltFromFile(p.Results.FolderData, p.Results.FileData,...
            p.Results.lowFreq, p.Results.highFreq, p.Results.bDivide, p.Results.bExpFit,false, p.Results.saveFilename);
    end
    
else
    OutData = NormFiltDirect(p.Results.FolderData, p.Results.FileData,...
        p.Results.lowFreq, p.Results.highFreq, p.Results.bDivide,...
        p.Results.bExpFit, p.Results.Freq);
    varargout{1} = OutData;
end

end

function OutData = NormFiltDirect(FolderData, OutData, lowFreq, highFreq, bDivide, bExpFit, Freq)

ExpFun = @(P,x) abs(P(1)).*exp(-abs(P(2)).*x) + abs(P(3)).*exp(-abs(P(4)).*x)...
    - abs(P(5))*x + P(6);
Opt = optimset(@fminsearch);
Opt.Display = 'off';

if isempty(Freq)
    FileData = dir(fullfile(FolderData, '*.mat'));
    Tags = {'red','green','yellow','fluo_'};
    idx = arrayfun(@(x) contains(FileData(x).name, Tags), 1:size(FileData,1));
    idx = find(idx,1,'first');
    Infos = matfile(fullfile(FolderData, FileData(idx).name));
    Freq = Infos.Freq;
end

% Temporal filtering
if( lowFreq > 0 )
    UseLPFilt = 1;
    if( (1/lowFreq) > (size(OutData,3)/Freq) )
        lowFreq = 1/((size(OutData,3)/Freq));
    end
    f = fdesign.lowpass('N,F3dB', 4, lowFreq, Freq); %Fluo lower Freq
    lpass = design(f,'butter');
else
    UseLPFilt = 0;
end
if( highFreq > 0 )
    UseHPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, highFreq, Freq);   %Fluo Higher Freq
    hpass = design(f,'butter');
else
    UseHPFilt = 0;
end

dims = size(OutData);
if( bExpFit )
    rng('shuffle');
    S = mean(reshape(OutData,[], dims(3)),1);
    B = fminsearch(@(P) sum((double(S) - ExpFun(P,(1:size(S,2)))).^2),...
        rand(1,6).*[30 1 20 1 1 double(mean(S))],Opt);
    Approx = ExpFun([B(1:4) 0 0],1:size(S,2));
    Pred =[ones(1, size(S,2)); linspace(0,1,size(S,2)); Approx]';
end
PrcLims = round(linspace(1, dims(1), 11));
fprintf('Progress: ');
for ind = 1:dims(1)
    Signal = double(squeeze(OutData(ind,:,:)));
    % Exponential Fit
    if( bExpFit )
       B = Pred\Signal';
       Approx = (Pred*B)';
       Signal = Signal./Approx;
    end
    if( UseLPFilt )
        LP_lowCutOff = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
    else
        LP_lowCutOff = ones(size(Signal));
    end
    if( UseHPFilt )
        LP_highCutOff = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
    else
        LP_highCutOff = Signal;
    end
    
    if( bDivide )
        OutData(ind,:,:) = single(LP_highCutOff./LP_lowCutOff);
    else
        OutData(ind,:,:) = single(LP_highCutOff-LP_lowCutOff);
    end
    
    if( any(ind == PrcLims) )
        idx = find(ind == PrcLims);
        fprintf('%d%%..', 10*(idx-1));
    end
    
end
fprintf('\n');   
end

function OutData = NormFiltFromFile(FolderData, FileName, lowFreq, highFreq, bDivide, bExpFit, bReturn, outFile)

if ~strcmp(FolderData(end), filesep)
    FolderData = [FolderData filesep];
end

% --- Load metadata
fMetaData = load([FolderData strrep(FileName,'.dat','.mat')]);
dimNames  = fMetaData.dim_names;
datSize   = [fMetaData.datSize fMetaData.datLength];
Freq      = fMetaData.Freq;

hasEvents = any(strcmpi(dimNames,'E'));

% --- Output filenames
if isempty(outFile)
    outName = [erase(FileName,'.dat') '_NormFilt.dat'];
else
    [~,outName,~] = fileparts(outFile);
end
outDat  = [FolderData outName '.dat'];
outMat  = [FolderData outName '.mat'];

% --- Preallocate output file
preallocateDatFile(outDat, fMetaData);

fidIn  = fopen([FolderData FileName], 'r');
cIn = onCleanup(@() safeFclose(fidIn));
fidOut = fopen(outDat, 'r+');
cOut = onCleanup(@() safeFclose(fidOut));

% ---------------- Filter design (unchanged) ----------------
if lowFreq > 0
    f = fdesign.lowpass('N,F3dB',4,lowFreq,Freq);
    lpass = design(f,'butter');
    UseLPFilt = true;
else
    UseLPFilt = false;
end

if highFreq > 0
    f = fdesign.lowpass('N,F3dB',4,highFreq,Freq);
    hpass = design(f,'butter');
    UseHPFilt = true;
else
    UseHPFilt = false;
end

% ---------------- Exponential fit helpers ----------------
ExpFun = @(P,x) abs(P(1)).*exp(-abs(P(2)).*x) + ...
                abs(P(3)).*exp(-abs(P(4)).*x) ...
                - abs(P(5))*x + P(6);
Opt = optimset('Display','off');

fprintf('NormalisationFiltering (file mode)\n');

%% =========================================================
% EVENT MODE: E,Y,X,T  
%% =========================================================
if hasEvents
    
    Ne = datSize(1);
    Ny = datSize(2);
    Nx = datSize(3);
    Nt = datSize(4);

    elemsPerTrial = Ny * Nx * Nt;
    bytesPerTrial = elemsPerTrial * getByteSize(fMetaData.Datatype); % single precision

    for e = 1:Ne
        fprintf('Trial %d / %d\n', e, Ne);

        % --- Read one full trial
        fseek(fidIn, (e-1)*bytesPerTrial, 'bof');
        slab = fread(fidIn, elemsPerTrial, '*single');
        slab = reshape(slab, Ny, Nx, Nt);

        % --- Exp fit (trial-wise)
        if bExpFit
            S = mean(reshape(slab,[],Nt),1);
            B = fminsearch(@(P) sum((double(S)-ExpFun(P,1:Nt)).^2), ...
                rand(1,6).*[30 1 20 1 1 mean(S)], Opt);
            Approx = ExpFun([B(1:4) 0 0],1:Nt);
            Pred = [ones(1,Nt); linspace(0,1,Nt); Approx]';
        end

        % --- Process Y lines
        for y = 1:Ny
            Signal = double(squeeze(slab(y,:,:)));

            if bExpFit
                B = Pred \ Signal';
                Signal = Signal ./ (Pred*B)';
                clear B
            end
            
            if UseLPFilt
                LP_low = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
            else
                LP_low = ones(size(Signal));
            end

            if UseHPFilt
                LP_high = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
            else
                LP_high = Signal;
            end
            clear Signal
            if bDivide
                slab(y,:,:) = single(LP_high ./ LP_low);
            else
                slab(y,:,:) = single(LP_high - LP_low);
            end
            clear LP_high LP_low
        end

        % --- Write trial back
        fseek(fidOut, (e-1)*bytesPerTrial, 'bof');
        fwrite(fidOut, slab, 'single');
        clear slab
        
    end

%% =========================================================
% NO EVENTS: Y,X,T  ? chunked over time
% =========================================================
else
    Ny = datSize(1);
    Nx = datSize(2);
    Nt = datSize(3);
    if bExpFit
        nChunks = calculateMaxChunkSize(prod(datSize)*4,2,.3);
    else
        nChunks = calculateMaxChunkSize(prod(datSize)*4,1,.3);
    end    
    
    chunkX  = ceil(Nx / nChunks);
    fprintf('Filtering data (%i chunk(s))...\n',nChunks)
    lastPct = -1;   % ensures 0% prints
    for c = 1:nChunks
                
        xStart = (c-1)*chunkX + 1;
        xEnd   = min(xStart + chunkX - 1, Nx);
        xIdx   = xStart:xEnd;
        fprintf('Chunk #%i [Reading data from file...]\n',c)
        slab = spatialSlabIO( ...
            'read', fidIn, Ny, Nx, Nt, xIdx, 'single');
        % --- Exp fit (slab-wise)
        if bExpFit
            S = mean(reshape(slab,[],Nt),1);
            B = fminsearch(@(P) sum((double(S)-ExpFun(P,1:Nt)).^2), ...
                rand(1,6).*[30 1 20 1 1 mean(S)], Opt);
            Approx = ExpFun([B(1:4) 0 0],1:Nt);
            Pred = [ones(1,Nt); linspace(0,1,Nt); Approx]';
        end
        
        fprintf('Chunk #%i [Temporal filtering...]\n',c)
        fprintf('Progress:')
        for y = 1:Ny
            
            Signal = double(squeeze(slab(y,:,:)));
            
            if bExpFit
                B = Pred \ Signal';
                Signal = Signal ./ (Pred*B)';
                clear B
            end
            
            if UseLPFilt
                LP_low = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
            else
                LP_low = ones(size(Signal));
            end

            if UseHPFilt
                LP_high = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
            else
                LP_high = Signal;
            end
            clear Signal
            if bDivide
                slab(y,:,:) = single(LP_high ./ LP_low);
            else
                slab(y,:,:) = single(LP_high - LP_low);
            end
            clear LP_high LP_low
            
            % Print progress
            pct = floor(100 * y / Ny);
            if pct ~= lastPct && mod(pct,10) == 0
                fprintf('%d%% ', pct);
                lastPct = pct;
            end
        end
        fprintf('\nChunk #%i [Writing to file...]\n',c)

        spatialSlabIO( ...
            'write', fidOut, Ny, Nx, Nt, xIdx, 'single', slab);
        
        fprintf('Chunk #%i [Completed]\n',c)
        clear slab
        
    end
end

fclose(fidIn);
fclose(fidOut);


% Output
if bReturn
    % Treat the outDat file as temporary, read it and delete
    fid = fopen(outDat,'r');
    OutData = fread(fid, inf, '*single');
    OutData = reshape(OutData, datSize);
    fclose(fid);
    % Delete 
    delete(outDat)
    delete(outMat)
else
    OutData = [];
end

end
