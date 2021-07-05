function outFile = calculate_response_amplitude(File, SaveFolder, varargin)
% CALCULATE_RESPONSE_AMPLITUDE in data split into events ("E" dimension in
% "dim_names" variable). It calculates the signal amplitude in time ("T").
% The calculation consists on a value "postEventTime" calculated in a time window
% between the event trigger and a given time in seconds ("timeWindow") minus a value
% "postEventTime". Default values for post and preEventTimes are 'max' and
% 'median'.
% (Note to self: improve function description.)
% 
% Inputs:
%   File: fullpath of functional imaging .DAT file.
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
%   opts (optional) : structure containing extra parameters.
% Output:
%   outFile: name of Output file.

% Defaults:
default_opts = struct('preEvent_value', 'median', 'postEvent_value', 'max', 'timeWindow', -1);
default_Output = 'amplitude_Map.dat'; 
%%% Arguments parsing and validation %%%
% Parse inputs:
p = inputParser;
addRequired(p,'File',@(x) isfile(x) & endsWith(x,'.dat'))
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
addOptional(p, 'Output', default_Output)
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%

% Optional arguments validation:
errID = 'Umitoolbox:calculate_respose_amplitude:InvalidDataType';
errMsg = 'Invalid Input. Input must be member of {"mean", "median", "min","max"} or a numeric scalar';
valid_Opts1 = @(x) (ismember(char(x), {'mean', 'median', 'min','max'}) || (isscalar(x) && isnumeric(x)));
valid_Opts2 = @(x) (isscalar(x) && (x == -1 || x>0));
assert(valid_Opts1(opts.preEvent_value), errID, errMsg);
assert(valid_Opts1(opts.postEvent_value), errID, errMsg);
assert(valid_Opts2(opts.timeWindow), errID, 'Input must be a scalar positive numeric value or "-1".');

% Map Data and metadata:
[mData, metaData] = mapDatFile(File);

% Get dimension names:
dims = metaData.dim_names;

% Validate if data has the following dimension names "E" and "T":
errID = 'Umitoolbox:calculate_respose_amplitude:WrongInput';
errMsg = 'Input Data is invalid. It must have an "Event" and a "Time" dimensions.';
assert(all(ismember({'E', 'T'}, dims)), errID, errMsg);

% load data:
data = mData.Data.(metaData.datName);

% Identify "T" dimension and permute data so Time is the first dimension:
idxT = find(strcmp('T', dims));
orig_dim_indx = 1:numel(dims);
new_dim_indx = [idxT setdiff(orig_dim_indx, idxT)];
data = permute(data, new_dim_indx);

% Store data size:
data_sz = size(data);

% Reshape data:
data = reshape(data,data_sz(1), []);

% Perform amplitude calculation:
trigFrame = round(metaData.preEventTime_sec * metaData.Freq);
switch opts.timeWindow
    case -1
        endFrame  = size(data,1);
    otherwise
        endFrame  = round(metaData.Freq * opts.timeWindow);
        if endFrame > size(data,1)
            endFrame = size(data,1);
        end
end

% Get baseline (pre_trigger) and postTrigger data:
bsln = data(1:metaData.preEventTime_sec,:);
postTrig = data(trigFrame:endFrame,:);

% Use aggregate function OR value defined by User:
bsln = applyAggFcn(bsln, opts.preEvent_value);
postTrig = applyAggFcn(postTrig, opts.postEvent_value);

% calculate amplitude:
Amp = postTrig - bsln;

% Reshape Amp to match data size:
new_sz = data_sz;
new_sz(1) = 1;
Amp = reshape(Amp, new_sz);

% Find singleton dimensions:
singletonDims = size(Amp) == 1;

% Permute Amp to original data size and remove singleton dimensions:
Amp = squeeze(permute(Amp,new_dim_indx)); 

% Do the same in dimension names:
singletonDims = singletonDims(new_dim_indx);
new_dim_names = dims(~singletonDims); 

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
