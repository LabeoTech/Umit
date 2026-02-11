function outData = genCorrelationMatrix(data, metaData, varargin)
% GENCORRELATIONMATRIX creates a correlation matrix from an Image time
% series with dimensions Y,X,T using Regions of Interest (ROIs) stored in a
% "ROImasks_xxxx.mat" file (created by the ROImanager App).
% The function offers 3 distinct ways of calculating the correlation
% between ROIs:
%   1- Centroid vs centroid: correlates the centroid pixels between ROIs.
%   2- Centroid vs Aggregation: correlates the centroid pixel of a source
%       ROI agains the aggregation of pixels in the target area. Here the
%       aggregation function applied comprises: max, min, mean and median.
%   3- Average vs Average: correlates the spatial average of pixels from
%       each ROI.
%
% Inputs:
%   data (3D numerical matrix): Image time series with dimensions {'Y','X','T}.
%   metaData (struct): structure containing smeta data associated with "data".
% Output:
%   outData: structure containing the correlation values of each ROI.

% Defaults: IMPORTANT, keep all default statements in one line each so the
% Pipeline Managers will be able to read it!
default_Output = 'corrMatrix.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('ROImasks_filename', 'ROImasks_data.mat', 'CorrAlgorithm', 'centroid_vs_centroid', 'SpatialAggFcn', 'mean','b_FisherZ_transform', false, 'b_genSPCMaps', false);
opts_values = struct('ROImasks_filename', {{'ROImasks_data.mat'}}, 'CorrAlgorithm',{{'centroid_vs_centroid','centroid_vs_agg', 'avg_vs_avg'}}, 'SpatialAggFcn', {{'mean', 'max', 'min', 'median'}},'b_FisherZ_transform',[true,false],'b_genSPCMaps',[true,false]);%  % This is here only as a reference for PIPELINEMANAGER.m.

default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) | ischar(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Optional Parameters:
% Validation criteria for optinal params:
valid_spatAgg = @(x) ismember(x.SpatialAggFcn, opts_values.SpatialAggFcn);
valid_corrAlg = @(x) ismember(x.CorrAlgorithm, opts_values.CorrAlgorithm);
% Add optional params:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && valid_spatAgg(x) && valid_corrAlg(x));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality'));
% Initialize Variables:
parse(p, data, metaData, varargin{:});
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
object = p.Results.object;
clear p
%%%%
% Find ROI file:
opts.ROImasks_filename = findMyROIfile(opts.ROImasks_filename,object);

% Validate data dimensions
errID = 'Umitoolbox:genCorrelationMatrix:WrongFormat';
errMsg = 'Invalid data format. The input Data must be an Image time series with dimensions {"Y","X","T"}';
assert(isequal(metaData.dim_names, {'Y','X','T'}),errID,errMsg);
% Load ROI data
roi_data = load(opts.ROImasks_filename,'-mat');
% Dispatch
if isnumeric(data)
    disp('Calculating Correlation matrix...')
    outData = genCorr_standard(data, metaData, roi_data, opts, object);
else
    disp('Calculating Correlation matrix (RAM-safe mode)...')
    outData = genCorr_RAMsafe(data, metaData, roi_data, opts, object);
end
disp('Done')

end

%% Local functions --------------------------------------------------------
function outData = genCorr_standard(data, metaData, roi_data, opts, object)

data_sz = size(data);
data = reshape(data, [], data_sz(3));

roi_names = {roi_data.ROI_info.Name}';
roi_corrVals = cell(size(roi_names));

centroid_list = arrayfun(@(x) ...
    find(bwmorph(x.Stats.ROI_binary_mask,'shrink',Inf)), ...
    roi_data.ROI_info, 'UniformOutput', false);

centroid_list = cellfun(@(x) x(1),centroid_list);

% ----- SPC Maps -----
if opts.b_genSPCMaps
    fprintf('Calculating Seed Pixel correlation maps...')
    SPCMaps = cell(length(centroid_list),1);
    
    for i = 1:length(centroid_list)
        fprintf('\nROI %i/%i=>',i,length(centroid_list));
        lastPct = -1;
        tmp_out = zeros(1,size(data,1),'single');
        for j = 1:size(tmp_out,2)
            rho = corrcoef(data(centroid_list(i),:)', data(j,:)');
            tmp_out(j) = rho(1,2);
            % Print progress
            pct = floor(100 *j/ size(tmp_out,2));
            if pct ~= lastPct && mod(pct,10) == 0
                fprintf('%d%% ', pct);
                lastPct = pct;
            end
        end
        SPCMaps{i} = reshape(tmp_out, data_sz(1), data_sz(2));
        if opts.b_FisherZ_transform
            SPCMaps{i} = ZFisher_truncated(SPCMaps{i});
        end
    end
end

% ----- Correlation matrix -----
if startsWith(opts.CorrAlgorithm,'centroid','IgnoreCase',true)
    
    switch opts.CorrAlgorithm
        
        case 'centroid_vs_centroid'
            roiVals = data(centroid_list,:);
            B = corrcoef(roiVals');
            
        case 'centroid_vs_agg'
            B = zeros(length(centroid_list),'single');
            target = arrayfun(@(x) data(x.Stats.ROI_binary_mask(:),:), ...
                roi_data.ROI_info,'UniformOutput',false);
            
            for i = 1:length(centroid_list)
                source = data(centroid_list(i),:);
                rho_vals = cell(size(target));
                for j = 1:numel(target)
                    tmp = target{j};
                    vals = zeros(1,size(tmp,1),'single');
                    for k = 1:size(tmp,1)
                        rho_tmp = corrcoef(source,tmp(k,:));
                        vals(k) = rho_tmp(1,2);
                    end
                    rho_vals{j} = vals;
                end
                
                switch opts.SpatialAggFcn
                    case 'mean'
                        B(i,:) = cellfun(@(x) mean(x,'omitnan'), rho_vals);
                    case 'median'
                        B(i,:) = cellfun(@(x) median(x,'omitnan'), rho_vals);
                    case 'min'
                        B(i,:) = cellfun(@(x) min(x,[],'omitnan'), rho_vals);
                    case 'max'
                        B(i,:) = cellfun(@(x) max(x,[],'omitnan'), rho_vals);
                end
            end
    end
    
else
    roiVals = zeros(numel(roi_corrVals),data_sz(3),'single');
    for i = 1:numel(roi_corrVals)
        roiVals(i,:) = mean( ...
            data(roi_data.ROI_info(i).Stats.ROI_binary_mask(:),:), ...
            1,'omitnan');
    end
    B = corrcoef(roiVals');
end

if opts.b_FisherZ_transform
    B = ZFisher_truncated(B);
end

out = arrayfun(@(i) B(i,:),1:size(B,1),'UniformOutput',false)';
dim_names = {'O','O'};

outData = save2Mat('', out ,roi_names, dim_names, ...
    'label',roi_names,'appendMetaData', metaData,'genFile', false);

if exist('SPCMaps','var')
    dim_names = {'Y','X','O'};
    if isempty(object)
        save2Mat(fullfile(pwd,'corrMat_SPCMaps.mat'),SPCMaps,roi_names,...
            dim_names,'appendMetaData',metaData,'genFile',true);
    else
        save2Mat(fullfile(object.SaveFolder,'corrMat_SPCMaps.mat'),SPCMaps,...
            roi_names,dim_names,'appendMetaData',metaData,'genFile',true);
    end
end

end

function outData = genCorr_RAMsafe(datFile, metaData, roi_data, opts, object)

% RAM-safe correlation matrix computation from .dat file
% Uses block-based streaming for SPCMaps

datatype = metaData.Datatype;
bytes    = getByteSize(datatype);

Ny = metaData.datSize(1);
Nx = metaData.datSize(2);
Nt = metaData.datLength;

roi_names = {roi_data.ROI_info.Name}';
nROI = numel(roi_names);

roi_masks = arrayfun(@(x) x.Stats.ROI_binary_mask(:), ...
    roi_data.ROI_info,'UniformOutput',false);

centroid_list = arrayfun(@(x) ...
    find(bwmorph(x.Stats.ROI_binary_mask,'shrink',Inf)), ...
    roi_data.ROI_info,'UniformOutput',false);

centroid_list = cellfun(@(x) x(1),centroid_list);

fid = fopen(datFile,'r');
c = onCleanup(@() safeFclose(fid));

bytesPerFrame = Ny * Nx * bytes;

%% =====================================================
% LOAD ROI TIME SERIES (centroids or averages)
%% =====================================================

roiVals = zeros(nROI, Nt, 'single');

for t = 1:Nt
    
    fseek(fid,(t-1)*bytesPerFrame,'bof');
    frame = fread(fid,Ny*Nx,['*' datatype]);
    
    for r = 1:nROI
        
        switch opts.CorrAlgorithm
            
            case 'centroid_vs_centroid'
                roiVals(r,t) = frame(centroid_list(r));
                
            case 'avg_vs_avg'
                roiVals(r,t) = mean(frame(roi_masks{r}),'omitnan');
                
            case 'centroid_vs_agg'
                roiVals(r,t) = frame(centroid_list(r));
        end
    end
end

%% =====================================================
% BLOCK-BASED SPC MAPS
%% =====================================================

if opts.b_genSPCMaps
    
    totalBytes = bytesPerFrame * Nt;
    nChunks = calculateMaxChunkSize(totalBytes,2,.2);
    framesPerChunk = ceil(Nt / nChunks);
    
    SPCMaps = cell(nROI,1);
    w = waitbar(0,'Creating SPCMaps...');
    
    for r = 1:nROI
        w.Name = ['ROI ' num2str(r) ' / ' num2str(nROI)];
        seed_ts = double(roiVals(r,:)');  % Nt × 1
        
        sum_x  = zeros(Ny*Nx,1,'double');
        sum_x2 = zeros(Ny*Nx,1,'double');
        sum_xy = zeros(Ny*Nx,1,'double');
        
        sum_y  = 0;
        sum_y2 = 0;
        Ntotal = 0;
        
        for cChunk = 1:nChunks
            
            tStart = (cChunk-1)*framesPerChunk + 1;
            tEnd   = min(cChunk*framesPerChunk, Nt);
            B = tEnd - tStart + 1;
            
            fseek(fid,(tStart-1)*bytesPerFrame,'bof');
            block = fread(fid, Ny*Nx*B, ['*' datatype]);
            block = reshape(block, Ny*Nx, B);
            block = double(block);
            
            seed_block = seed_ts(tStart:tEnd)';
            
            sum_x  = sum_x  + sum(block,2);
            sum_x2 = sum_x2 + sum(block.^2,2);
            sum_xy = sum_xy + block * seed_block';
            
            sum_y  = sum_y  + sum(seed_block);
            sum_y2 = sum_y2 + sum(seed_block.^2);
            
            Ntotal = Ntotal + B;
            waitbar(cChunk/nChunks,w);
        end
        
        numerator = Ntotal*sum_xy - sum_x*sum_y;
        denom_x   = Ntotal*sum_x2 - sum_x.^2;
        denom_y   = Ntotal*sum_y2 - sum_y^2;
        
        rho = numerator ./ sqrt(denom_x .* denom_y);
        rho(isnan(rho)) = 0;
        
        SPCMaps{r} = reshape(single(rho), Ny, Nx);
        
        if opts.b_FisherZ_transform
            SPCMaps{r} = ZFisher_truncated(SPCMaps{r});
        end
                
    end
    
    delete(w);
    disp('Seed Pixel Correlation Maps created!');
end

%% =====================================================
% CORRELATION MATRIX
%% =====================================================

if strcmp(opts.CorrAlgorithm,'centroid_vs_centroid') || ...
        strcmp(opts.CorrAlgorithm,'avg_vs_avg')
    
    B = corrcoef(roiVals');
    
else
    % centroid_vs_agg
    
    B = zeros(nROI,'single');
    
    for i = 1:nROI
        
        source = roiVals(i,:);
        
        for j = 1:nROI
            
            target_ts = zeros(1,Nt,'single');
            
            for t = 1:Nt
                fseek(fid,(t-1)*bytesPerFrame,'bof');
                frame = fread(fid,Ny*Nx,['*' datatype]);
                target_ts(t) = mean(frame(roi_masks{j}),'omitnan');
            end
            
            rho_tmp = corrcoef(source,target_ts);
            B(i,j) = rho_tmp(1,2);
        end
    end
end

if opts.b_FisherZ_transform
    disp('Applying Z-Fisher transform to correlation matrix...')
    B = ZFisher_truncated(B);
end

%% =====================================================
% SAVE CORRELATION MATRIX
%% =====================================================

out = arrayfun(@(i) B(i,:),1:size(B,1),'UniformOutput',false)';
dim_names = {'O','O'};

outData = save2Mat('', out ,roi_names, dim_names, ...
    'label',roi_names,'appendMetaData', metaData,'genFile', false);

%% =====================================================
% SAVE SPCMaps
%% =====================================================

if exist('SPCMaps','var')
    
    dim_names = {'Y','X','O'};
    
    if isempty(object)
        save2Mat(fullfile(pwd,'corrMat_SPCMaps.mat'),SPCMaps,roi_names,...
            dim_names,'appendMetaData',metaData,'genFile',true);
    else
        save2Mat(fullfile(object.SaveFolder,'corrMat_SPCMaps.mat'),...
            SPCMaps,roi_names,dim_names,...
            'appendMetaData',metaData,'genFile',true);
    end
end

end


function out = ZFisher_truncated(data)
% This function truncates the ZFisher transformed data between -0.998 and
% 0.998 Pearson's rho values. This avoids the creation of Infinite values
% at Pearsons' correlation values of -1 and 1.

lim = 0.998;
data(data < -lim) = -lim;
data(data > lim) = lim;
out = atanh(data);
end


