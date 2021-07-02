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
addRequired(p,'File',@isfile & endsWith('.dat'))
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x) && ...
    ismember(x.aggregateFcn, {'mean', 'max', 'min', 'median', 'mode', 'sum', 'std'}));
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%

% Try to find dimension names in file metadata. If not, throw error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save DATA, METADATA and DIMENSION_NAMES to DATFILE:
datFile = fullfile(SaveFolder, Output);
save2Dat(datFile, Amp, new_dim_names);
% Output file names
outFile = Output;
end

% Local function:
function out = applyAggFcn(vals, aggfcn)
% APPLYAGGFCN performs the aggregate function of name "fcn_name" on the 1st
% dimension of the data "vals". All aggregate functions EXCLUDE NaNs!

switch aggfcn
    case 'mean'
        out = nanmean(vals, 1);
    case 'median'
        out = median(vals, 1, 'omitnan');
%     case 'mode'
%         out = mode(vals, 1);
%     case 'std'
%         out = std(vals, 0, 1, 'omitnan');
    case 'max'
        out = max(vals, [], 1, 'omitnan');
    case 'min'
        out = min(vals, [], 1, 'omitnan');
%     case 'sum'
%         out = sum(vals, 1, 'omitnan');
    otherwise
        out = repmat(aggfcn,size(vals));
end
end
