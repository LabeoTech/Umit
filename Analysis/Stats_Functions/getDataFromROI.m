function outData = getDataFromROI(data, metaData, varargin)
% GETDATAFROMROI extracts and aggregates data from regions of interest
% (ROIs) in imaging data using an "ROI_xxxxxx.mat" file located in
% subject's folder.

% Inputs:
%   data: numerical matrix containing imaging data.
%   metaData: .mat file with meta data associated with "data".
% Output:
%   outData: structure containing stats-ready data extracted from ROIs.

% Defaults:
default_Output = 'ROI_data.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('ROImasks_filename', 'ROImasks_data.mat', 'SpatialAggFcn', 'mean');
opts_values = struct('ROImasks_filename', {{'ROImasks_data.mat'}}, 'SpatialAggFcn',{{'none','mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}});% This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) | ischar(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ...
    ismember(x.SpatialAggFcn, opts_values.SpatialAggFcn));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality'));

% Parse inputs:
parse(p,data, metaData, varargin{:});
% Initialize Variables:
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
object = p.Results.object;
clear p
%%%%%%%%%%%%%%%%

% Parse File path to find subject folder:
opts.ROImasks_filename = findMyROIfile(opts.ROImasks_filename,object);

% Load ROI file:
roi_data = load(opts.ROImasks_filename);
if isnumeric(data)
    % Execute standard mode
    outData = getROIdata_standardMode(data,metaData,roi_data,opts, object);
else
    % Execute in RAM safe mode
    outData = getROIdata_RAMsafeMode(data, metaData, roi_data, opts, object);
end

end

%% Local functions --------------------------------------------------------

function outData = getROIdata_standardMode(data,metaData,roi_data,opts,object)
% STANDARD MODE FOR ROI DATA EXTRACTION


% locate "X" and "Y" dimensions in metaData and in ROI info:
dim_names = metaData.dim_names;
[~,yxLoc] = ismember({'Y','X'}, dim_names);

% Check if frame size is the same as the one in ROI file:
data_sz = size(data);
errID = 'Umitoolbox:getDataFromROI:IncompatibleSizes';
errMsg = 'Data file frame size is different from the one in ROI file.';
assert(isequal(data_sz(yxLoc), size(roi_data.img_info.imageData)), errID, errMsg)
% permute matrix:
orig_dim = 1:ndims(data);
new_dim = [yxLoc setdiff(orig_dim, yxLoc)];
data = permute(data, new_dim);
dim_names = dim_names(new_dim);
% reshape matrix:
data_sz = data_sz(new_dim);
data = reshape(data, prod(data_sz([1 2])),[]);
% extract ROI pixel values from data:
roi_names = {roi_data.ROI_info.Name}';
roi_pixVals = cell(size(roi_names));
disp(['Extracting ' opts.SpatialAggFcn ' values from ROIs...'])
for i = 1:length(roi_pixVals)
    roi_msk = roi_data.ROI_info(i).Stats.ROI_binary_mask;
    pixVals = data(roi_msk(:),:);
    % Perform aggregate operation:
    pixVals = applyAggFcn(pixVals, opts.SpatialAggFcn);
    % reshape data back to retrieve dimensions:
    if length(data_sz) > 2
        pixVals = reshape(pixVals, [size(pixVals,1) data_sz(3:end)]);
    end
    roi_pixVals{i} = pixVals;
end
disp('Done');
% Create new dimension names:
if strcmp(opts.SpatialAggFcn, 'none')
    % This is an special case where the User wants to keep all individual
    % pixels from the ROI. In this case, the dimension "P" (for pixel) is
    % created after the "O" (observation).
    new_dim_names = {'O', 'P', dim_names{3:end}};
else
    new_dim_names ={'O', dim_names{3:end}};
end
if isa(metaData, 'matlab.io.MatFile')
    metaData.Properties.Writable = true;
end
metaData.ROIfile = opts.ROImasks_filename;
if isempty(object)
    outData = save2Mat('', roi_pixVals, roi_names,...
        new_dim_names, 'appendMetaData', metaData, 'genFile', false);
else
    outData = save2Mat('', roi_pixVals, roi_names,...
        new_dim_names, 'appendMetaData', metaData, 'genFile', false,...
        'appendObjectInfo',object);
end
end

function outData = getROIdata_RAMsafeMode(datFile, metaData, roi_data, opts, object)
% RAM-safe ROI extraction from .dat file
%
% EYXT  -> processed trial by trial
% YXT   -> processed frame by frame
%
% Output is identical to standard mode.

dim_names = metaData.dim_names;
datatype  = metaData.Datatype;
bytes     = getByteSize(datatype);

hasEvents = any(strcmpi(dim_names,'E'));

roi_names   = {roi_data.ROI_info.Name}';
roi_pixVals = cell(size(roi_names));

% Locate Y and X
[~,yxLoc] = ismember({'Y','X'}, dim_names);
if strcmp(opts.SpatialAggFcn, 'none')
    new_dim_names = [{'P'},setdiff(dim_names,{'Y','X'})];
else
    new_dim_names = setdiff(dim_names,{'Y','X'});
end

datSize = [metaData.datSize metaData.datLength];

fprintf('Extracting %s values from ROIs (RAM-safe mode)...\n',opts.SpatialAggFcn)

fid = fopen(datFile,'r');
c = onCleanup(@() safeFclose(fid));

%% =========================================================
% EVENT MODE: E,Y,X,T
%% =========================================================
if hasEvents
    
    Ne = datSize(1);
    Ny = datSize(2);
    Nx = datSize(3);
    Nt = datSize(4);
    
    elemsPerTrial = Ny * Nx * Nt;
    bytesPerTrial = elemsPerTrial * bytes;
    
    % Preallocate per ROI
    for r = 1:numel(roi_pixVals)
        if strcmp(opts.SpatialAggFcn,'none')
            nPix = nnz(roi_data.ROI_info(r).Stats.ROI_binary_mask);
            roi_pixVals{r} = zeros(Ne, nPix, Nt, datatype);
        else
            roi_pixVals{r} = zeros(Ne, Nt, datatype);
        end
    end
    
    for e = 1:Ne
        
        fseek(fid,(e-1)*bytesPerTrial,'bof');
        slab = fread(fid, elemsPerTrial, ['*' datatype]);
        slab = reshape(slab, Ny, Nx, Nt);
        
        for r = 1:numel(roi_pixVals)
            
            roi_msk = roi_data.ROI_info(r).Stats.ROI_binary_mask;
            linIdx  = find(roi_msk(:));
            
            tmp = reshape(slab, Ny*Nx, Nt);
            pixVals = tmp(linIdx,:);
            
            pixVals = applyAggFcn(pixVals, opts.SpatialAggFcn);
            
            roi_pixVals{r}(e,:,:) = pixVals;
        end
        
        clear slab tmp
    end
    
    
    %% =========================================================
    % NO EVENTS: Y,X,T  -> frame-by-frame
    %% =========================================================
else
    
    Ny = datSize(1);
    Nx = datSize(2);
    Nt = datSize(3);
    
    elemsPerFrame = Ny * Nx;
    bytesPerFrame = elemsPerFrame * bytes;
    
    % Preallocate
    for r = 1:numel(roi_pixVals)
        if strcmp(opts.SpatialAggFcn,'none')
            nPix = nnz(roi_data.ROI_info(r).Stats.ROI_binary_mask);
            roi_pixVals{r} = zeros(nPix, Nt, datatype);
        else
            roi_pixVals{r} = zeros(1, Nt, datatype);
        end
    end
    w = waitbar(0,'Extracting ROI data ...');
    w.Name  = 'getDataFromROI';
    for t = 1:Nt
        
        fseek(fid,(t-1)*bytesPerFrame,'bof');
        frame = fread(fid, elemsPerFrame, ['*' datatype]);
        frame = reshape(frame, Ny, Nx);
        
        frame2D = reshape(frame, Ny*Nx, 1);
        
        for r = 1:numel(roi_pixVals)
            
            roi_msk = roi_data.ROI_info(r).Stats.ROI_binary_mask;
            pixVals = frame2D(roi_msk(:));
            pixVals = applyAggFcn(pixVals, opts.SpatialAggFcn);
            
            roi_pixVals{r}(:,t) = pixVals;
        end
        if mod(t,50) == 0 || t == 1 || t == Nt
            waitbar(t/Nt,w);
            clear frame frame2D
        end
        
        
    end
    new_dim_names = [{'O'}, new_dim_names];  % O,T
    fprintf('Done\n');
    
    % Metadata handling (unchanged)
    if isa(metaData,'matlab.io.MatFile')
        metaData.Properties.Writable = true;
    end
    metaData.ROIfile = opts.ROImasks_filename;
    
    if isempty(object)
        outData = save2Mat('', roi_pixVals, roi_names,...
            new_dim_names, 'appendMetaData', metaData, 'genFile', false);
    else
        outData = save2Mat('', roi_pixVals, roi_names,...
            new_dim_names, 'appendMetaData', metaData, 'genFile', false,...
            'appendObjectInfo', object);
    end
    
end
close(w);
end


function out = applyAggFcn(vals, fcn_name)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the 1st
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!

switch fcn_name
    case 'mean'
        out = mean(vals, 1, 'omitnan');
    case 'median'
        out = median(vals, 1, 'omitnan');
    case 'mode'
        out = mode(vals, 1);
    case 'std'
        out = std(vals, 0, 1, 'omitnan');
    case 'max'
        out = max(vals, [], 1, 'omitnan');
    case 'min'
        out = min(vals, [], 1, 'omitnan');
    case 'sum'
        out = sum(vals, 1, 'omitnan');
    otherwise
        out = vals;
end
end


