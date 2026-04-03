function [outData, metaData] = apply_aggregate_function(data, metaData, varargin)
% APPLY_AGGREGATE_FUNCTION applies an aggregate function along Time (T) or Event (E)
% dimensions in low RAM mode.
%
% Usage:
%   outData = apply_aggregate_function(data, metaData, opts)
%
% Inputs:
%   data     : numerical matrix or filename of .dat file
%   metaData : structure or matfile associated with data
%   opts     : struct with fields:
%              aggregateFcn  : 'mean','median','mode','std','max','min','sum'
%              dimensionName : 'T' or 'E'
%
% Outputs:
%   outData  : aggregated data (or filename in low RAM mode)
%   metaData : updated metadata

% Default options
default_Output = 'aggFcn_applied.dat'; %#ok This is here for PIPELINEMANAGER.M.
default_opts = struct('aggregateFcn','mean','dimensionName','T');
opts_values = struct('aggregateFcn',{{'mean','std','max','min','sum'}},'dimensionName',{{'T','E'}});

p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) || ischar(x));
addRequired(p,'metaData', @(x) isstruct(x) || isa(x,'matlab.io.MatFile'));
addOptional(p,'opts',default_opts,@(x) isstruct(x) && ...
    ismember(x.aggregateFcn, opts_values.aggregateFcn) && ...
    ismember(x.dimensionName, opts_values.dimensionName));
parse(p,data,metaData,varargin{:});
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p

% Determine if low RAM mode
bLowRAM = ischar(data);

% Map dimension
dimName = upper(opts.dimensionName);
data_dim_names = metaData.dim_names;

if ~ismember(dimName, data_dim_names)
    error('apply_aggregate_function:InvalidDim',...
        'Dimension %s not found in metadata.', dimName);
end
dimIdx = find(strcmp(dimName,data_dim_names));


if bLowRAM
    % Preallocate output file
    outFile = fullfile(pwd,'AGGFCN_DATA.dat'); % default name
    outMeta = metaData;
    sz = [metaData.datSize metaData.datLength];
    % Remove aggregated dimension
    newSz = sz;
    if strcmpi(opts.dimensionName,'t')
        newSz(dimIdx) = [];
        outMeta.dim_names(dimIdx) = [];
    end
    outMeta.datSize = newSz([1:2]);
    outMeta.datLength = newSz([3:end]);
    
    
    fidIn = fopen(data,'r');
    c_in = onCleanup(@() safeFclose(fidIn));
    fidOut = fopen(outFile,'w');
    c_out = onCleanup(@() safeFclose(fidOut));
    % Determine chunk size
    
    dims = sz;
    % Only T or E aggregation supported
    switch dimName
        case 'T'
            % Aggregate along time for data with or without events                        
            
            if numel(dims) == 4
                % Data has E,Y,X,T
                nE  = dims(1);
                nY  = dims(2);
                nX  = dims(3);
                nT  = dims(4);
                
                % Preallocate output: aggregate along time => size [E,Y,X]
                aggData = zeros(nE, nY, nX, 'single');
                
                for e = 1:nE
                    trialData = readTrial(fidIn,e,[nE,nY,nX,nT],'single');
                    trialData = permute(trialData,[3 1 2]);                   
                    % Aggregate along T (3rd dimension)
                    aggData(e,:,:) = calcAgg(trialData, opts.aggregateFcn);
                end
                
            elseif numel(dims) == 3
                % Data has Y,X,T
                nY = dims(1);
                nX = dims(2);
                nT = dims(3);
                
                % Preallocate output: aggregate along time => size [Y,X]
                aggData = zeros(nY, nX, 'single');
                
                nElem = nY*nX*nT;
                nIter = calculateMaxChunkSize(nY*nX*nT*4,2);
                if nIter >1
                    % Determine X chunking
                    xPerIter = ceil(nX / nIter);
                    
                    xStart = 1;
                    
                    for ii = 1:nIter
                        
                        % X indices for this chunk
                        xEnd   = min(xStart + xPerIter - 1, nX);                                                
                        
                        % Read slab: [Y, xCount, T]
                        slab = spatialSlabIO('read',fidIn, nY,nX,nT, xStart:xEnd, 'single');
                        
                        % Rearrange to [T, Y, Xchunk]
                        slab = permute(slab, [3 1 2]);
                        
                        % Aggregate along T ? [Y, Xchunk]
                        aggChunk = calcAgg(slab, opts.aggregateFcn);
                        
                        % Store result into output
                        aggData(:, xStart:xEnd) = aggChunk;
                        
                        % Advance X index
                        xStart = xEnd + 1;
                    end

                else
                    % Read the entire block
                    trialData = fread(fidIn, nElem, ['*' metaData.Datatype]);
                    trialData = reshape(trialData, nY, nX, nT);
                    trialData = permute(trialData,[3 1 2]);
                    % Aggregate along T (3rd dimension)
                    aggData(:,:) = calcAgg(trialData, opts.aggregateFcn);
                end
                
            else
                error('Unsupported data dimensions for case ''T''.');
            end
            
            % Write output
            fwrite(fidOut, cast(aggData,'single'), 'single');
            
            
        case 'E'
            % Aggregate along Event dimension by condition for interleaved data [E,Y,X,T]
            nE  = dims(1);
            nY  = dims(2);
            nX  = dims(3);
            nT  = dims(4);                       
            
            % Unique conditions
            evIdx = unique(metaData.eventID);
            nCond = numel(evIdx);
            
            % Preallocate output depending on aggregation
            switch lower(opts.aggregateFcn)
                case {'mean','sum','std'}
                    aggData = zeros(nCond, nY, nX, nT, 'double'); % double for stability
                case {'min','max'}
                    aggData = zeros(nCond, nY, nX, nT, 'single');
                otherwise
                    error('Unsupported aggregation function');
            end
                        
            for c = 1:nCond
                condEvents = find(metaData.eventID == evIdx(c));
                nEventsCond = numel(condEvents);
                
                % Reset temporary aggregation for this condition
                switch lower(opts.aggregateFcn)
                    case {'mean','sum'}
                        tempAgg = zeros(nY, nX, nT, 'double');
                    case {'min','max'}
                        tempAgg = zeros(nY, nX, nT, 'single');
                    case 'std'
                        meanVal = zeros(nY, nX, nT, 'double');
                        M2 = zeros(nY, nX, nT, 'double');
                end
                
                % Loop over all events of this condition
                for ee = 1:nEventsCond
                    e = condEvents(ee);
                    
                    % Read event e from interleaved file
                    trialData = readTrial(fidIn,e,[nE,nY,nX,nT],'single');
                                        
                    % Incremental aggregation
                    switch lower(opts.aggregateFcn)
                        case {'mean','sum'}
                            tempAgg = tempAgg + double(trialData);
                        case 'max'
                            if ee == 1
                                tempAgg = trialData;
                            else
                                tempAgg = max(tempAgg, trialData);
                            end
                        case 'min'
                            if ee == 1
                                tempAgg = trialData;
                            else
                                tempAgg = min(tempAgg, trialData);
                            end
                        case 'std'
                            if ee == 1
                                meanVal = double(trialData);
                                M2 = zeros(nY,nX,nT,'double');
                            else
                                delta = double(trialData) - meanVal;
                                meanVal = meanVal + delta / ee;
                                M2 = M2 + delta .* (double(trialData) - meanVal);
                            end
                    end
                end
                
                % Finalize aggregation for this condition
                switch lower(opts.aggregateFcn)
                    case 'mean'
                        aggData(c,:,:,:) = tempAgg / nEventsCond;
                    case 'sum'
                        aggData(c,:,:,:) = tempAgg;
                    case 'std'
                        aggData(c,:,:,:) = sqrt(M2 / (nEventsCond - 1));
                    otherwise
                        aggData(c,:,:,:) = tempAgg; % min/max already finalized
                end
            end
            
            % Write output
            fwrite(fidOut, cast(aggData, 'single'), 'single');
    end
    metaData = genMetaData(aggData,outMeta.dim_names,outMeta);
    save(strrep(outFile,'.dat','.mat'),'-struct','metaData');
    % Save aggregate data to file
    fclose(fidIn);
    fclose(fidOut);
    outData = outFile;
        
else
    newDims = metaData.dim_names;
    % Standard: full RAM
    if strcmpi(opts.dimensionName, 't')
        outData = applyAggFcnFull(data, opts.aggregateFcn, dimIdx);
        newDims(dimIdx) = [];
    else
        % For events, perform aggregate function by condition
        disp('working on events')
        sz = [metaData.datSize metaData.datLength];
        evIdx = unique(metaData.eventID);
        outData = zeros(numel(evIdx),sz(2),sz(3),sz(4),'single');
        for ii = 1:length(evIdx)
            idx = metaData.eventID == evIdx(ii);
            outData(ii,:,:,:) = applyAggFcnFull(data(idx,:,:,:),opts.aggregateFcn,dimIdx);
        end
        % Update reduced event ID list in meta data
        metaData.eventID = evIdx;
    end
    
    metaData = genMetaData(outData,newDims,metaData);
    disp('Done')
end
end

%% Local functions
function out = calcAgg(vals, aggFcn)
switch aggFcn
    case 'mean',   out = mean(vals,1,'omitnan');
    case 'median', out = median(vals,1,'omitnan');
    case 'mode',   out = mode(vals,1);
    case 'std',    out = std(vals,0,1,'omitnan');
    case 'max',    out = max(vals,[],1,'omitnan');
    case 'min',    out = min(vals,[],1,'omitnan');
    case 'sum',    out = sum(vals,1,'omitnan');
    otherwise,     out = vals;
end
end

function outData = applyAggFcnFull(data, aggFcn, dimIdx)
% Permute to bring dimIdx first
permOrder = [dimIdx setdiff(1:ndims(data), dimIdx)];
dataP = permute(data, permOrder);
szP = size(dataP);
dataP = reshape(dataP, szP(1), []);
out = calcAgg(dataP, aggFcn);
outData = reshape(out, [1 szP(2:end)]);
% Permute back
outData = ipermute(outData, permOrder);
end
