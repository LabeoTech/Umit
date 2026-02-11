function outFile = genRetinotopyMaps(data, metaData, SaveFolder, varargin)
% GENRETINOTOPYMAPS create amplitude and phase maps for each cardinal
% direction as well as the azimuth and elevation maps when at least two
% perpendicular directions are present in the input data.
%
% The timestamps marking the individual directions must be encoded in the "events.mat" file.

% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T").
%   metaData: .mat file with meta data associated with "data".
%   opts (optional): structure containing extra parameters. See "default_opts" variable below for details!
%
% Outputs:
%   outFile: filenames of the generated data.

% Defaults:
default_Output = {'AzimuthMap.dat' 'ElevationMap.dat'};  %#ok. This line is here just for Pipeline management.
% Parse inputs:
default_opts = struct('nSweeps', 1, 'b_useAverageMovie', false, 'ViewingDist_cm', 0,'ScreenXsize_cm',0,'ScreenYsize_cm',0);
opts_values = struct('nSweeps', [1,Inf], 'b_useAverageMovie', [true,false], 'ViewingDist_cm', [0,Inf],'ScreenXsize_cm',[0,Inf],'ScreenYsize_cm',[0,Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
p = inputParser;
addRequired(p,'data');
addRequired(p,'metaData');
addRequired(p,'SaveFolder', @isfolder);
addOptional(p,'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p,data, metaData, SaveFolder, varargin{:});
opts = p.Results.opts;
SaveFolder = p.Results.SaveFolder;
dataIn = p.Results.data;
clear p

% Validate inputs
errID1 = 'umIToolbox:genRetinotopyMaps:WrongInput';
errID2 = 'umIToolbox:genRetinotopyMaps:MissingInput';
assert(all(ismember({'Y','X','T'},metaData.dim_names)), errID1,'Input data must have dimensions "Y","X","T".')
assert(isfile(fullfile(SaveFolder,'events.mat')), errID2,'"events.mat" file not found!');
evntInfo = load(fullfile(SaveFolder,'events.mat'));
assert(all(ismember(evntInfo.eventNameList, {'0','90','180','270'})), ...
    errID2,'Invalid directions found in "events.mat" file! Must be 0,90,180,270.')

outFile = {};

% Decide mode
if ischar(dataIn) % Low-RAM mode
    outFile = RAMsafeMode(dataIn, metaData, SaveFolder, evntInfo, opts, outFile);
else % Standard in-memory mode
    outFile = standardMode(dataIn, metaData, SaveFolder, evntInfo, opts, outFile);
end

%% ==================== NESTED STANDARD MODE ====================
    function outFile = standardMode(data, metaData, SaveFolder, evntInfo, opts, outFile)
        ampMaps = cell(size(evntInfo.eventNameList));
        phiMaps = ampMaps;
        freqFFT = opts.b_useAverageMovie * 0 + round(opts.nSweeps)+1;
        
        framestamps = round(evntInfo.timestamps*metaData.Freq);
        w = waitbar(0,'Calculating FFT ...','Name','genRetinotopyMaps');
        w.Children.Title.Interpreter = 'none';
        
        for ind = 1:numel(evntInfo.eventNameList)
            w.Children.Title.String = ['Calculating FFT for direction ' evntInfo.eventNameList{ind}]; drawnow
            indxOn = find(evntInfo.eventID == ind & evntInfo.state == 1);
            indxOff = find(evntInfo.eventID == ind & evntInfo.state == 0);
            
            if opts.b_useAverageMovie
                % Average DeltaR movie
                bsln_len = round(mean(framestamps(indxOn(2:end)) - framestamps(indxOff(1:end-1))));
                trial_len = round(mean(framestamps(indxOff) - framestamps(indxOn)));
                avg_mov = zeros([metaData.datSize, trial_len + bsln_len],'single');
                for ii = 1:length(indxOn)
                    DeltaR = data(:,:,framestamps(indxOn(ii))-bsln_len : framestamps(indxOn(ii))+trial_len-1) - ...
                        median(data(:,:,framestamps(indxOn(ii))-bsln_len : framestamps(indxOn(ii))-1), 3,'omitnan');
                    avg_mov = avg_mov + DeltaR;
                end
                avg_mov = avg_mov/length(indxOn);
                fDat = fft(avg_mov,[],3);
            else
                frOn = [];
                for ii = 1:length(indxOn)
                    frOn = [frOn, framestamps(indxOn(ii)):framestamps(indxOff(ii))-1];
                end
                fDat = fft(data(:,:,frOn),[],3);
            end
            ampMaps{ind} = (abs(fDat(:,:,freqFFT)) * 2) / size(fDat,3);
            phiMaps{ind} = mod(-angle(fDat(:,:,freqFFT)),2*pi);
            waitbar(ind/numel(evntInfo.eventNameList),w);
        end
        close(w);
        
        outFile = saveMaps(ampMaps, phiMaps, metaData, opts, SaveFolder, evntInfo, outFile);
    end

%% ==================== NESTED RAM-SAFE MODE ====================
    function outFile = RAMsafeMode(datFile, metaData, SaveFolder, evntInfo, opts, outFile)
        nY = metaData.datSize(1);
        nX = metaData.datSize(2);
        nT = metaData.datLength;
        
        
        ampMaps = cell(size(evntInfo.eventNameList));
        phiMaps = ampMaps;
        framestamps = round(evntInfo.timestamps*metaData.Freq);
        freqFFT = opts.b_useAverageMovie * 0 + round(opts.nSweeps)+1;
        
        w = waitbar(0,'Calculating FFT (Low RAM usage) ...','Name','genRetinotopyMaps');
        w.Children.Title.Interpreter = 'none';
        fidIn = fopen(datFile,'r');
        cIn = onCleanup(@() safeFclose(fidIn));
        for ind = 1:numel(evntInfo.eventNameList)
            w.Children.Title.String = ['Calculating FFT for direction ' evntInfo.eventNameList{ind}]; drawnow
            indxOn = find(evntInfo.eventID == ind & evntInfo.state==1);
            indxOff = find(evntInfo.eventID == ind & evntInfo.state==0);
            
            % Prepare aggregated movie
            ampMaps{ind} = zeros(nY,nX,'single');
            phiMaps{ind} = zeros(nY,nX,'single');
            
            if opts.b_useAverageMovie
                % Compute baseline and trial length
                bsln_len = round(mean(framestamps(indxOn(2:end)) - framestamps(indxOff(1:end-1))));
                trial_len = round(mean(framestamps(indxOff) - framestamps(indxOn)));
                total_len = trial_len + bsln_len;
                
                % Preallocate averaged movie
                avg_mov = zeros(nY, nX, total_len, 'single');
                
                % Calculate number of elements per frame
                bytesPerElem = 4; % single
                elemsPerFrame = nY * nX;
                
                for ii = 1:length(indxOn)
                    % Determine frame indices for this trial
                    tStart = max(framestamps(indxOn(ii)) - bsln_len,1);
                    tEnd   = min(framestamps(indxOn(ii)) + trial_len - 1,nT);
                    nFrames = tEnd - tStart + 1;
                    assert(nFrames>0, 'Failed to split data by trials')
                    
                    % Preallocate trial data
                    trialData = zeros(nY, nX, nFrames, 'single');
                    
                    % Read frames sequentially from file
                    for f = 1:nFrames
                        % Compute byte offset: (frameIndex-1) * bytesPerFrame
                        fseek(fidIn, (tStart+f-2) * elemsPerFrame * bytesPerElem, 'bof');
                        % Read a single frame
                        frame = fread(fidIn, elemsPerFrame, '*single');
                        trialData(:,:,f) = reshape(frame, [nY, nX]);
                    end
                    
                    % Compute baseline from pre-stimulus period
                    baseline = median(trialData(:,:,1:bsln_len), 3, 'omitnan');
                    
                    % Compute DeltaR and sum into avg_mov
                    avg_mov = avg_mov + (trialData - baseline);
                end
                
                % Divide by number of trials to get average
                avg_mov = avg_mov / length(indxOn);
                
                % Compute FFT along the 3rd dimension (time)
                fDat = fft(avg_mov, [], 3);
                ampMaps{ind} = (abs(fDat(:,:,freqFFT))*2)/size(fDat,3);
                phiMaps{ind} = mod(-angle(fDat(:,:,freqFFT)),2*pi);
                waitbar(ind/numel(evntInfo.eventNameList),w);
            
            else
                % Collect all frames for this direction
                frOn = [];
                for ii = 1:length(indxOn)
                    frOn = [frOn, framestamps(indxOn(ii)):framestamps(indxOff(ii))-1];
                end
                nFrames = length(frOn);
                ampMap = ampMaps{ind};
                phiMap = phiMaps{ind};
                % Chunk along X
                nChunks = calculateMaxChunkSize(nY * nX * nT * 4,1,.1);  % number of chunks
                chunkX  = ceil(nX / nChunks);
                
                for c = 1:nChunks
                    fprintf('Processing spatial slab [%i/%i] for direction %s\n',c,nChunks,evntInfo.eventNameList{ind})
                    xStart = (c-1)*chunkX + 1;
                    xEnd   = min(xStart + chunkX - 1, nX);
                    xIdx   = xStart:xEnd;
                    
                    slabData = spatialSlabIO('read', fidIn, nY, nX, nT, xIdx, 'single');
                    slabData = slabData(:,:,frOn);
                    
                    fSlab = fft(slabData,[],3);
                    ampSlab = (abs(fSlab(:,:,freqFFT))*2)/size(fSlab,3);
                    phiSlab = mod(-angle(fSlab(:,:,freqFFT)),2*pi);
                    ampMap(:,xIdx) = ampSlab;
                    phiMap(:,xIdx) = phiSlab;
                    clear fSlab
                end
                ampMaps{ind} = ampMap;
                phiMaps{ind} = phiMap;
                waitbar(ind/numel(evntInfo.eventNameList),w);
            end
        end
        close(w);
        fclose(fidIn);
        % Save maps
        outFile = saveMaps(ampMaps, phiMaps, metaData, opts, SaveFolder, evntInfo, outFile);
    end

%% ==================== HELPER FUNCTION TO SAVE MAPS ====================
    function outFile = saveMaps(ampMaps, phiMaps, metaData, opts, SaveFolder, evntInfo, outFile)
        AzimMap = zeros([metaData.datSize 2],'single');
        ElevMap = AzimMap;
        [idxAz,indxAz] = ismember({'0','180'}, evntInfo.eventNameList);
        [idxEl,indxEl] = ismember({'90','270'}, evntInfo.eventNameList);
        
        if all(idxAz)
            AzimMap(:,:,1) = mean(cat(3,ampMaps{indxAz}),3);
            AzimMap(:,:,2) = pi + ((phiMaps{indxAz(1)} - phiMaps{indxAz(2)})/2);
        end
        if all(idxEl)
            ElevMap(:,:,1) = mean(cat(3,ampMaps{indxEl}),3);
            ElevMap(:,:,2) = pi + ((phiMaps{indxEl(1)} - phiMaps{indxEl(2)})/2);
        end
        
        if all([opts.ViewingDist_cm, opts.ScreenXsize_cm, opts.ScreenYsize_cm])>0
            va_az = atand(opts.ScreenXsize_cm/(2*opts.ViewingDist_cm)); va_az = [-va_az va_az];
            va_el = atand(opts.ScreenYsize_cm/(2*opts.ViewingDist_cm)); va_el = [-va_el va_el];
            AzimMap(:,:,2) = rescale(AzimMap(:,:,2), va_az(1), va_az(2));
            ElevMap(:,:,2) = rescale(ElevMap(:,:,2), va_el(1), va_el(2));
        end
        
        if sum(AzimMap(:)) ~= 0
            md = genMetaData(AzimMap,{'Y','X','F'}, metaData);
            filename = fullfile(SaveFolder, 'AzimuthMap.dat');
            save2Dat(filename,AzimMap,md);
            outFile = [outFile;{filename}];
        end
        if sum(ElevMap(:)) ~= 0
            md = genMetaData(ElevMap,{'Y','X','F'}, metaData);
            filename = fullfile(SaveFolder, 'ElevationMap.dat');
            save2Dat(filename,ElevMap,md);
            outFile = [outFile;{filename}];
        end
    end

end
