function [outData, metaData] = apply_aggregate_function(data, metaData, varargin)
% APPLY_AGGREGATE_FUNCTION applies an aggregate function to one
% dimensions of a .DAT file. 
% !! Removed option to perform aggregation over multiple dimensions. BrunoO
% (09-06-2022).

% Inputs:
%   data: numerical matrix containing imaging data.
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing the function's parameters:
%       aggregateFcn (default = "mean") : name of the aggregate function.
%       dimensionName (default = "T") : name of the dimension(s) to perform
%       the calculation.
% Output:
%   outData: numerical matrix containing aggregated imaging data.   
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_Output = 'aggFcn_applied.dat'; %ok This is here for PIPELINEMANAGER.M.
default_opts = struct('aggregateFcn', 'mean', 'dimensionName', 'T');
opts_values = struct('aggregateFcn', {{'mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}}, 'dimensionName',{{'X','Y','Z','T','E'}});% This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
% Parse inputs:
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ismember(x.aggregateFcn, opts_values.aggregateFcn) && ...
    ismember(x.dimensionName, opts_values.dimensionName));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
data = p.Results.data; 
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%%%%%%%%%%%%%%%%%%%%

% Parse dimension names from opts struct:
str = opts.dimensionName;
str = split(str, ',');
dim_names = cellfun(@(x) upper(strip(x)), str, 'UniformOutput', false); 
% Validate if dimension name(s) is(are) in File meta data:
errID = 'Umitoolbox:apply_aggregate_function:InvalidName';
errMsg = 'Input dimension name(s) was not found in file meta data';
[idx_dims, dimVec] = ismember(dim_names, metaData.dim_names);
assert(all(idx_dims), errID, errMsg);
data_dim_names = metaData.dim_names;
% Look for "E"vent dimension:
idx_evnt = strcmp('E',data_dim_names(dimVec));

% Permute data to bring "dimVec" to firsts dimensions:
orig_sz = size(data);
data = permute(data,[dimVec, setdiff(1:ndims(data), dimVec)]);
perm_dim_names = metaData.dim_names([dimVec, setdiff(1:ndims(data), dimVec)]); 
perm_sz = size(data);
% Apply aggregate function in each dimension:
for i = 1:length(dimVec)
    if idx_evnt(i)
        [data, new_eventID] = applyAggFcn(data, opts.aggregateFcn, i, metaData.eventID);
    else
        data = applyAggFcn(data, opts.aggregateFcn, i, []);
    end
end
% Find singleton dimensions:
singDims = ( size(data) == 1 );
% Permute data back and remove singleton dimensions:
[~,locB] = ismember(metaData.dim_names,perm_dim_names);
outData = squeeze(permute(data, locB));
% Update dim_names:
singDims = singDims(locB);
% new_dim_names = metaData.dim_names;
data_dim_names (singDims) = [];

% Create metaData structure based on aggregated data:
extraParams = metaData;
if exist('new_eventID', 'var')
    % IF "E"vent dimension was processed, update eventID variable in metaData file:
    extraParams.eventID = new_eventID;    
end
metaData = genMetaData(outData, data_dim_names,extraParams);
end

% Local function:
function [out, evntID_out] = applyAggFcn(vals, aggfcn, idxDim, evntList)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!
% Inputs:
% vals = multi-D numeric data.
% idxDim = index of the dimension of vals to be aggregated.
% evntList (int array OR empty) = array of indices to group by (used in
% "Event" dimension).

if ~isempty(evntList)
    % If the dimension is an Event, group by event and aggregate each group:
    % Permute "vals" to have "idxDim" as first dimension and reshape "vals" to 
    % force to have 2 dimensions:
    orig_sz = size(vals);
    % permute:
    dimorder = [idxDim,setdiff(1:ndims(vals), idxDim)];
    vals = permute(vals,dimorder);
    perm_sz = size(vals);
    % reshape:
    vals = reshape(vals, size(vals,1), []);
    % Create empty matrix with size = unique(evntList) x size(vals,2):
    ID_list = unique(evntList);
    out = zeros(numel(ID_list), size(vals,2), class(vals));
    % Calculate aggregate function across events :
    for i = 1:numel(ID_list)
        idx = (evntList == ID_list(i));
        out(i,:) = calcAgg(vals(idx,:), 1);
    end
    % Reshape data to match permuted "vals":
    out = reshape(out,[size(out,1), perm_sz(2:end)]);
    % Permute data to original size of vals:
    out = ipermute(out,dimorder);    
    evntID_out = ID_list;
else
    out = calcAgg(vals,idxDim);
end

    function out = calcAgg(vals,idxDim)
        switch aggfcn
            case 'mean'
                out = mean(vals, idxDim,'omitnan');
            case 'median'
                out = median(vals, idxDim, 'omitnan');
            case 'mode'
                out = mode(vals, idxDim);
            case 'std'
                out = std(vals, 0, idxDim, 'omitnan');
            case 'max'
                out = max(vals, [], idxDim, 'omitnan');
            case 'min'
                out = min(vals, [], idxDim, 'omitnan');
            case 'sum'
                out = sum(vals, idxDim, 'omitnan');
            otherwise
                out = vals;
        end
    end
end
