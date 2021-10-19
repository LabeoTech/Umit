function outFile = trim_movie(File, SaveFolder, varargin)
% TRIM_MOVIE crops a Image time series by N Frames at the beginning and/OR end of the movie.
% 
% Inputs:
%   File: fullpath of functional imaging .DAT file.
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
%   opts (optional) : structure containing extra parameters.
% Output:
%   outFile: name of Output file.

% Defaults:
default_opts = struct('crop_start_sec', 0, 'crop_end_sec',0);
default_Output = 'cropped_mov.dat'; 
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
% Map Data and metadata:
[mData, metaData] = mapDatFile(File);
dim_names = metaData.dim_names;
% Input File validation:
errID = 'Umitoolbox:trim_movie:InvalidDataType';
errMsg = 'Invalid Input. Input must be a 3D matrix with third dimension  = "T"';
assert(numel(dim_names) == 3 && strcmp(dim_names{3}, 'T'), errID, errMsg);
errMsg = 'must be a positive number';
valid_Opts = @(x) x>=0;
assert(valid_Opts(opts.crop_start_sec), errID, ['Start crop time ' errMsg]);
assert(valid_Opts(opts.crop_end_sec), errID, ['End crop time ' errMsg]);

% load data:
data = mData.Data.(metaData.datName);

% Identify "T" dimension and permute data so Time is the first dimension:
idxT = find(strcmp('T', dim_names));
% Calculate cropping frames:
fr_start = round(metaData.Freq*opts.crop_start_sec);
if fr_start == 0
    fr_start = 1;
end
fr_stop = round(size(data,idxT) - metaData.Freq*opts.crop_end_sec);

orig_dim_indx = 1:numel(dim_names);
new_dim_indx = [idxT setdiff(orig_dim_indx, idxT)];
data = permute(data, new_dim_indx);
% Store data size:
data_sz = size(data);
% Reshape data:
data = reshape(data,data_sz(1), []);
% Validate if cropped movie has a minimal length of 1 frame:
errMsg = ['Failed to crop from frame ' num2str(fr_start) ' to frame ' num2str(fr_stop)];
assert(~isempty(fr_start:fr_stop), errID, errMsg);
% Crop Movie:
disp('Cropping movie...')
% Recover data dimensions:
new_sz = data_sz;
new_sz(1) = length(fr_start:fr_stop);
data = data(fr_start:fr_stop,:);
data = reshape(data,new_sz);
data = permute(data,[2:numel(dim_names) 1]);

% Save DATA, METADATA and DIMENSION_NAMES to DATFILE:
[~,filename,ext] = fileparts(File);
outFile = ['cropped_' filename ext];
datFile = fullfile(SaveFolder, outFile);
save2Dat(datFile, data, dim_names);
end