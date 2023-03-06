function outData = genCorrelationMatrix(data, metaData,SaveFolder, varargin)
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
default_opts = struct('ROImasks_filename', 'ROImasks_data.mat', 'CorrAlgorithm', 'centroid_vs_centroid', 'SpatialAggFcn', 'mean','b_FisherZ_transform', false, 'SPCMapFileName', 'none');
opts_values = struct('ROImasks_filename', {{'ROImasks_data.mat'}}, 'CorrAlgorithm',{{'centroid_vs_centroid','centroid_vs_agg', 'avg_vs_avg'}}, 'SpatialAggFcn', {{'mean', 'max', 'min', 'median'}},'b_FisherZ_transform',[true,false],'SPCMapFileName',{{'none'}});%  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% Validation criteria for optinal params:
valid_spatAgg = @(x) ismember(x.SpatialAggFcn, opts_values.SpatialAggFcn);
valid_corrAlg = @(x) ismember(x.CorrAlgorithm, opts_values.CorrAlgorithm);
% Add optional params:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && valid_spatAgg(x) && valid_corrAlg(x));
% Initialize Variables:
parse(p, data, metaData,SaveFolder, varargin{:});
data = p.Results.data;
metaData = p.Results.metaData;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
%%%%
% Check if ROImasks file exist:
[~,ROIfilename,~] = fileparts(opts.ROImasks_filename);
opts.ROImasks_filename = fullfile(SaveFolder, [ROIfilename,'.mat']);
folder = strrep(SaveFolder, '\', '\\');
assert(isfile(opts.ROImasks_filename),'Umitoolbox:genCorrelationMatrix:FileNotFound',['ROI file not found in ' folder]);

% Validate data dimensions
errID = 'Umitoolbox:genCorrelationMatrix:WrongFormat';
errMsg = 'Invalid data format. The input Data must be an Image time series with dimensions {"Y","X","T"}';
assert(isequal(metaData.dim_names, {'Y','X','T'}),errID,errMsg);
%
% Load ROI file:
roi_data = load(opts.ROImasks_filename);
% reshape matrix:
data_sz = size(data);
data = reshape(data, [],data_sz(3));
% Calculate the correlation for each ROI:
roi_names = {roi_data.ROI_info.Name}';
roi_corrVals = cell(size(roi_names));
% Create list of centroids:
centroid_list = arrayfun(@(x) find(bwmorph(x.Stats.ROI_binary_mask,'shrink',Inf)),...
    roi_data.ROI_info, 'UniformOutput', false);
% In cases where there are non-contiguous ROIs, get first centroid (Arbitrary decision here!):
centroid_list = cellfun(@(x) x(1),centroid_list);
% Create SPCMaps for each centroid:
if ~strcmpi(opts.SPCMapFileName,'none')
    SPCMaps = cell(length(centroid_list),1);
    w = waitbar(0,'Creating SPCMaps...');
    for i = 1:length(centroid_list)
        tmp_out = zeros(1,size(data,1),'single');
        for j = 1:size(tmp_out,2)
            rho = corrcoef(data(centroid_list(i),:)', data(j,:)');
            tmp_out(j) = rho(1,2);
        end
        SPCMaps{i} = reshape(tmp_out, data_sz(1), data_sz(2));
        if opts.b_FisherZ_transform
            % Apply Z-Fisher transform to SPCMaps
            SPCMaps{i} = ZFisher_truncated(SPCMaps{i});
        end
        waitbar(i/length(centroid_list),w);
    end
    delete(w);
    disp('Seed Pixel Correlation Maps created!');
end

% Get centroids, if applicable:
if startsWith(opts.CorrAlgorithm, 'centroid', 'IgnoreCase', true)
    switch opts.CorrAlgorithm
        case 'centroid_vs_centroid'
            roiVals = data(centroid_list,:);
            clear data;
            % Calculate Pearson's correlation:
            B = corrcoef(roiVals');
        case 'centroid_vs_agg'
            B = zeros(length(centroid_list),length(centroid_list), 'single');
            target = arrayfun(@(x) data(x.Stats.ROI_binary_mask(:),:), roi_data.ROI_info,...
                'UniformOutput',false);
            for i = 1:length(centroid_list)
                source = data(centroid_list(i),:);
                rho_vals = cell(size(target));
                for j = 1:numel(target)
                    tmp = target{j};
                    vals = zeros(1,size(target{j},1), 'single');
                    for k = 1:size(tmp,1)
                        rho_tmp = corrcoef(source,tmp(k,:));
                        vals(k) = rho_tmp(1,2);
                    end
                    rho_vals{j} = vals;
                end
                switch opts.SpatialAggFcn
                    %'mean', 'max', 'min', 'median'
                    case 'mean'
                        B(i,:) = cellfun(@(x) mean(x, 'omitnan'), rho_vals);
                    case 'median'
                        B(i,:) = cellfun(@(x) median(x, 'omitnan'), rho_vals);
                    case 'min'
                        B(i,:) = cellfun(@(x) min(x,[],'omitnan'), rho_vals);
                    case 'max'
                        B(i,:) = cellfun(@(x) max(x,[],'omitnan'), rho_vals);
                    otherwise
                        % This shouldn't be reached. Only in cases where
                        % we want to make a huge matrix such as the fig. 4
                        % in Bauer et al., 2017??
                        disp('Huge Matrix not available yet!')
                        error('Centroid versus No aggregation option is unavailable for now!')
                end
            end
    end
else
    % In case where avg_vs_avg
    roiVals = zeros(numel(roi_corrVals),data_sz(3), 'single');
    for i = 1:numel(roi_corrVals)
        roiVals(i,:) = mean(data(roi_data.ROI_info(i).Stats.ROI_binary_mask(:),:),1,'omitnan');
    end
    clear data;
    % Calculate Pearson's correlation:
    B = corrcoef(roiVals');
end
% Apply Z-Fisher transformation to the correlation matrix:
if opts.b_FisherZ_transform
    disp('Applying Z-Fisher transform to correlation matrix...')
    B = ZFisher_truncated(B);
end
%%% Save Data to .mat file:
% Create cell array per observation:
out = cell(size(B,1),1);
for i = 1:size(B,1)
    out{i} = B(i,:);
end
% Create new dimension names as {'O', 'O'}:
dim_names = {'O','O'};
outData = genDataMetaStructure(out, roi_names, dim_names, metaData, 'label',roi_names);
% Create .MAT files with SPCMaps:
if exist('SPCMaps', 'var')
    dim_names = {'Y', 'X','O'};
    % When the function is called from DataViewer
    out = genDataMetaStructure(SPCMaps, roi_names, dim_names, metaData);
    [~,outName,~] = fileparts(opts.SPCMapFileName);
    save(fullfile(SaveFolder,[outName '.mat']),'-struct','out','-v7.3');
end
end

% Local function
function out = ZFisher_truncated(data)
% This function truncates the ZFisher transformed data between -0.998 and
% 0.998 Pearson's rho values. This avoids the creation of Infinite values
% at Pearsons' correlation values of -1 and 1.

lim = 0.998;
data(data < -lim) = -lim;
data(data > lim) = lim;
out = atanh(data);
end


