function [outData, metaData] = trim_movie(data, metaData, varargin)
% TRIM_MOVIE crops a Image time series by N Frames at the beginning and/OR end of the movie.

% Inputs:
%   data: numerical matrix containing image time series (dimensions {"Y", "X", "T"}).
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing extra parameters.

% Outputs: 
%   outData: trimmed numerical matrix.
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_opts = struct('crop_start_sec', 0, 'crop_end_sec',0);
default_Output = 'cropped_mov.dat'; %#ok
%%% Arguments parsing and validation %%%
% Parse inputs:
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
outData = p.Results.data; 
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p
%%%%
dim_names = metaData.dim_names;
% Input File validation:
errID = 'Umitoolbox:trim_movie:InvalidDataType';
errMsg = 'Invalid Input. Input must be a 3D matrix with third dimension  = "T"';
assert(numel(dim_names) == 3 && strcmp(dim_names{3}, 'T'), errID, errMsg);
errMsg = 'must be a positive number';
valid_Opts = @(x) x>=0;
assert(valid_Opts(opts.crop_start_sec), errID, ['Start crop time ' errMsg]);
assert(valid_Opts(opts.crop_end_sec), errID, ['End crop time ' errMsg]);


% Identify "T" dimension and permute data so Time is the first dimension:
idxT = find(strcmp('T', dim_names));
% Calculate cropping frames:
fr_start = round(metaData.Freq*opts.crop_start_sec);
if fr_start == 0
    fr_start = 1;
end
fr_stop = round(size(outData,idxT) - metaData.Freq*opts.crop_end_sec);

orig_dim_indx = 1:numel(dim_names);
new_dim_indx = [idxT setdiff(orig_dim_indx, idxT)];
outData = permute(outData, new_dim_indx);
% Store data size:
data_sz = size(outData);
% Reshape data:
outData = reshape(outData,data_sz(1), []);
% Validate if cropped movie has a minimal length of 1 frame:
errMsg = ['Failed to crop from frame ' num2str(fr_start) ' to frame ' num2str(fr_stop)];
assert(~isempty(fr_start:fr_stop), errID, errMsg);
% Crop Movie:
disp('Cropping movie...')
% Recover data dimensions:
new_sz = data_sz;
new_sz(1) = length(fr_start:fr_stop);
outData = outData(fr_start:fr_stop,:);
outData = reshape(outData,new_sz);

outData = permute(outData,[2:numel(dim_names) 1]);
% Create metaData: 
metaData = genMetaData(outData, dim_names, metaData);

end