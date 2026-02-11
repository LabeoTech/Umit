function [outData, metaData] = split_data_by_event(data, metaData, SaveFolder, varargin)
% SPLIT_DATA_BY_EVENT Reshapes an image time series dataset into a 4D matrix
% organized by event.
%
%   outData dimensions: {E, Y, X, T}
%
% Inputs:
%   data        - 3D numeric matrix {Y,X,T} or char path to .dat file
%   metaData    - Metadata structure or matlab.io.MatFile
%   SaveFolder  - Folder containing events.mat and output location
%
% Optional inputs (name-value via struct):
%   preEventTime_sec  - Seconds before event ('auto' or numeric)
%   postEventTime_sec - Seconds after event ('auto' or numeric)
%   PadWith           - Padding value for out-of-bound snippets
%                       ('mean', 'NaN', or numeric)
%
% Outputs:
%   outData  - 4D matrix {E,Y,X,T} (empty if low-RAM mode)
%   metaData - Updated metadata for output data

% Defaults (used by PipelineManager)
default_Output = 'data_splitByEvent.dat'; %#ok<NASGU>
default_opts = struct('preEventTime_sec','auto', 'postEventTime_sec','auto', 'PadWith','mean');
opts_values = struct( 'preEventTime_sec',{{'auto',Inf}}, 'postEventTime_sec',{{'auto',Inf}}, 'PadWith',{{'mean','NaN',Inf}});%#ok For PipelineManager

%% Input parsing
p = inputParser;
addRequired(p,'data', @(x) (isnumeric(x) && ndims(x)==3) || ischar(x));
addRequired(p,'metaData', @(x) isstruct(x) || isa(x,'matlab.io.MatFile'));
addRequired(p,'SaveFolder', @isfolder);
addOptional(p,'opts', default_opts, @(x) (isstruct(x) && ~isempty(x) && ...
     isfield(x,'PadWith') && ismember(x.PadWith,{'mean','NaN'})) || isnumeric(x));
parse(p, data, metaData, SaveFolder, varargin{:});

data     = p.Results.data;
metaData = p.Results.metaData;
folder   = p.Results.SaveFolder;
opts     = p.Results.opts;

clear p

%% Load events file
evFile = fullfile(folder,'events.mat');
if ~isfile(evFile)
    error('umIToolbox:split_data_by_event:FileNotFound', ...
        'Event file ("events.mat") not found in %s', folder);
end
evDat = load(evFile);

%% Determine pre/post event times
if strcmp(opts.preEventTime_sec,'auto') || strcmp(opts.postEventTime_sec,'auto')
    if isfield(metaData,'preEventTime_sec')
        opts.preEventTime_sec  = metaData.preEventTime_sec;
        opts.postEventTime_sec = metaData.postEventTime_sec;
    else
        tmTrial = round(mean(diff(evDat.timestamps(evDat.state==1)),'omitnan'));
        opts.preEventTime_sec  = round(0.2*tmTrial,2);
        opts.postEventTime_sec = tmTrial - opts.preEventTime_sec;
    end
end

%% Derived parameters
szdat     = [metaData.datSize metaData.datLength]; % [Y X T]
sr        = metaData.Freq;

preFr     = round(sr * opts.preEventTime_sec);
postFr    = round(sr * opts.postEventTime_sec);
centralFr = preFr + 1;

len_trial = preFr + postFr;
timestamps = evDat.timestamps(evDat.state==1);
n_trial    = numel(timestamps);

new_dims = {'E','Y','X','T'};

%% Low-RAM vs in-memory mode
if ischar(data)
    b_RAMSafeMode = true;

    newMD              = metaData;
    newMD.dim_names    = new_dims;
    newMD.datSize      = [n_trial szdat(1)];
    newMD.datLength    = [szdat(2) len_trial];
    szYXTE = [newMD.datSize, newMD.datLength];
    szYXTE = szYXTE([2 3 4 1]);
    outFile = fullfile(SaveFolder,'DATABYEVENTS.dat');
    preallocateDatFile(outFile, newMD);

    fidOut = fopen(outFile,'r+');
    cOut = onCleanup(@() safeFclose(fidOut));

    [~,fname,ext] = fileparts(data);
    fidIn = fopen(fullfile(SaveFolder,[fname ext]),'r');
    cIn = onCleanup(@() safeFclose(fidIn));
else
    b_RAMSafeMode = false;
    outData = nan([n_trial szdat(1) szdat(2) len_trial],'single');
end

disp('Splitting data by events...')

%% Constants for I/O
bytesPerElem = getByteSize('single');
frameBytes   = szdat(1) * szdat(2) * bytesPerElem;

%% Main loop
for i = 1:n_trial
    fix_snippet = false;

    trialFr = round(sr * timestamps(i));
    startFrIn = trialFr - preFr;
    stopFrIn  = trialFr + postFr - 1;

    if startFrIn < 1
        startFrIn = 1;
        fix_snippet = true;
    end
    if stopFrIn > szdat(3)
        stopFrIn = szdat(3);
        fix_snippet = true;
    end

    % Read snippet
    if b_RAMSafeMode
        fseek(fidIn, (startFrIn-1)*frameBytes, 'bof');
        nFrames = stopFrIn - startFrIn + 1;
        snippet = fread(fidIn, szdat(1)*szdat(2)*nFrames, '*single');
        snippet = reshape(snippet, szdat(1), szdat(2), []);
    else
        snippet = data(:,:,startFrIn:stopFrIn);
    end

    % Placement indices
    startFrOut = centralFr - (trialFr - startFrIn);
    stopFrOut  = startFrOut + size(snippet,3) - 1;

    trialData = nan(szdat(1), szdat(2), len_trial, 'single');

    if fix_snippet
        switch opts.PadWith
            case 'mean'
                avg = mean(snippet,3,'omitnan');
                trialData = repmat(avg,1,1,len_trial);
            case 'NaN'
                % already NaN
            otherwise
                trialData(:) = opts.PadWith;
        end
    end
    
    trialData(:,:,startFrOut:stopFrOut) = snippet;
    
    if b_RAMSafeMode
        % Write trial to output file (permuted to YXTE for faster writing)
        fprintf('Writing trial #%i to file...\n',i)
        writeTrial_YXTE(fidOut,i,trialData,szYXTE,'single');        
    else
        outData(i,:,:,:) = trialData;
    end
end




disp('Done!')

%% Update metadata

extraParams = metaData;
extraParams.preEventTime_sec  = opts.preEventTime_sec;
extraParams.postEventTime_sec = opts.postEventTime_sec;
extraParams.eventID           = evDat.eventID(evDat.state==1);
extraParams.eventNameList     = evDat.eventNameList;

if ~b_RAMSafeMode
    metaData = genMetaData(outData, new_dims, extraParams);
else
    metaData = newMD;
    metaData.preEventTime_sec  = opts.preEventTime_sec;
    metaData.postEventTime_sec = opts.postEventTime_sec;
    metaData.eventID           = evDat.eventID(evDat.state==1);
    metaData.eventNameList     = evDat.eventNameList;
    
    % Close files
    fclose(fidIn);
    fclose(fidOut);
    % Permute back the output file
    
    permuteDat_YXTE_to_EYXT_inplace(outFile,szYXTE,'single');
    outData = outFile;
end

end
