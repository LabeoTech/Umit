function outFile = apply_aggregate_function(File, SaveFolder, varargin)
% APPLY_AGGREGATE_FUNCTION applies an aggregate function to one or more
% dimensions of a .DAT file. 
% Inputs:
%   File : fullpath of functional imaging .DAT file.
%   SaveFolder : path to save the output file.
%   Output (optional) : Name of outFile.
%   opts (optional) : structure containing the function's parameters:
%       aggregateFcn (default = "mean") : name of the aggregate function.
%       dimensionName (default = "T") : name of the dimension(s) to perform
%       the calculation.
% Output:
%   outFile : name of Output file.

% Defaults:
default_opts = struct('aggregateFcn', 'mean', 'dimensionName', 'T');
default_Output = 'aggFcn_applied.dat'; 
%%% Arguments parsing and validation %%%
% Parse inputs:
p = inputParser;
addRequired(p,'File',@(x) isfile(x) & endsWith(x,'.dat'))
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ismember(x.aggregateFcn, {'mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}) && ...
    ~isempty(x.dimensionName) && (iscell(x.dimensionName) && ischar(x.dimenqsionName{:}) || ...
    ischar(x.dimensionName)));
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;

% Map .DAT file to memory:
[mData, metaData] = mapDatFile(File);

% Parse dimension names from opts struct:
str = opts.dimensionName;
str = split(str, ',');
dim_names = cellfun(@(x) upper(strip(x)), str, 'UniformOutput', false); 
% Validate if dimension name(s) is(are) in File meta data:
errID = 'Umitoolbox:apply_aggregate_function:InvalidName';
errMsg = 'Input dimension name(s) was not found in file meta data';
[idx_dims, dimVec] = ismember(dim_names, metaData.dim_names);
assert(all(idx_dims), errID, errMsg);
% Look for "E"vent dimension:
idx_evnt = strcmp('E', dim_names(dimVec));

% Load Data:
data = mData.Data.(metaData.datName);
% Permute data to bring "dimVec" to firsts dimensions:
orig_sz = size(data);
data = permute(data,[dimVec, setdiff(1:ndims(data), dimVec)]);
perm_sz = size(data);
% Apply aggregate function in each dimension:
for i = 1:length(dimVec)
    if idx_evnt(i)
        [data, new_eventID] = applyAggFcn(data, opts.aggregateFcn, dimVec(i), metaData.eventID);
    else
        data = applyAggFcn(data, opts.aggregateFcn, dimVec(i), []);
    end
end
% Find singleton dimensions
singDims = ( size(data) == 1 );
% Permute data:
[~,locB] = ismember(orig_sz, perm_sz);
data = permute(data, locB);
% Remove singleton dimensions:
data = squeeze(data);
% Update dim_names:
singDims = singDims(locB);
new_dim_names = metaData.dim_names;
new_dim_names (singDims) = [];
% Save DATA, METADATA and DIMENSION_NAMES to DATFILE:
[~,filename, ~] = fileparts(mData.Filename);
filename = [opts.aggregateFcn '_over_' [dim_names{:}] '_' filename '.dat'];
datFile = fullfile(SaveFolder, filename);
save2Dat(datFile, data, new_dim_names);
% IF "E"vent dimension was processed, update eventID variable in metaData file:
if exist('new_eventID', 'var')
    [~, metaData] = mapDatFile(datFile);
    metaData.Properties.Writable = true;
    metaData.eventID = new_eventID;
end
% Output filename:
outFile = filename;
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
    vals = permute(vals, [idxDim,setdiff(1:ndims(vals), idxDim)]);
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
    [~,locB] = ismember(orig_sz, perm_sz);
    out = permute(out, locB);
    evntID_out = ID_list;
else
    out = calcAgg(vals,idxDim);
end

    function out = calcAgg(vals,idxDim)
        switch aggfcn
            case 'mean'
                out = nanmean(vals, idxDim);
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
