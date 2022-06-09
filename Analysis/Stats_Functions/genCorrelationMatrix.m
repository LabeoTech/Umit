function outDataStat = genCorrelationMatrix(data, metaData, varargin)
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
%   outDataStat: structure containing the correlation values of each ROI.

% Defaults: IMPORTANT, keep all default statements in one line each so the
% Pipeline Managers will be able to read it!
default_Output = 'corrMatrix.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('ROImasks_filename', 'ROImasks_data.mat', 'CorrAlgorithm', 'centroid_vs_centroid', 'SpatialAggFcn', 'mean','b_FisherZ_transform', false);
opts_values = struct('ROImasks_filename', {{'ROImasks_data.mat'}}, 'CorrAlgorithm',{{'centroid_vs_centroid','centroid_vs_agg', 'avg_vs_avg'}}, 'SpatialAggFcn', {{'none','mean', 'max', 'min', 'median'}},'b_FisherZ_transform',[true,false]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.

default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
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
%
% Load ROI file:
roi_data = load(opts.ROImasks_filename);
% reshape matrix:
data_sz = size(data);
data = reshape(data, [],data_sz(3));
% Calculate the correlation for each ROI:
roi_names = {roi_data.ROI_info.Name}';
roi_corrVals = cell(size(roi_names));

% Get centroids, if applicable:
if startsWith(opts.CorrAlgorithm, 'centroid', 'IgnoreCase', true)
    centroid_list = arrayfun(@(x) find(bwmorph(x.Stats.ROI_binary_mask,'shrink',Inf)),...
        roi_data.ROI_info, 'UniformOutput', false); 
    % In cases where there are non-contiguous ROIs, get first centroid (Arbitrary decision here!):
    centroid_list = cellfun(@(x) x(1),centroid_list); 
    switch opts.CorrAlgorithm
        case 'centroid_vs_centroid'
            roiVals = data(centroid_list,:);            
            clear data;
            % Calculate Pearson's correlation:
            B = corr(roiVals');                    
        case 'centroid_vs_agg'
            B = zeros(length(centroid_list),length(centroid_list), 'single');
            target = arrayfun(@(x) data(x.Stats.ROI_binary_mask(:),:), roi_data.ROI_info,...
                'UniformOutput',false);
            for i = 1: length(centroid_list)
                source = data(centroid_list(i),:);                
                rho_vals = cellfun(@(x) corr(source',x'), target, 'UniformOutput',false);
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
                        error('Centroic versus No aggregation option is unavailable for now!')                        
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
    B = corr(roiVals');    
end
% Apply Z-Fisher transformation to the correlation matrix:
if opts.b_FisherZ_transform
    disp('Applying Z-Fisher transform to correlation matrix...')
   B = atanh(B);
end                     
%%% Save Data to .mat file:
% Create cell array per observation:
out = cell(size(B,1),1);
for i = 1:size(B,1)
    out{i} = B(i,:);
end
% Create new dimension names as {'O', 'O'}:
dim_names = {'O','O'};
outDataStat = save2Mat([], out ,roi_names, dim_names, 'label',roi_names ,...
    'appendMetaData', metaData,'genFile', false);
end
