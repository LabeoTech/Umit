function outFile = SeedPixCorr(File, SaveFolder, varargin)
% This function reduces the matrix DATA in DATFILENAME to a 64 x 64, calculates
% a pixel-wise temporal correlation.

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('imageResizeTo', 64, 'FisherZ_transform', false);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File:
default_Output = 'SeedPixCorr.dat';
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
%%%%
% Open memMapfile:
mmData = mapDatFile(File);
metaData = matfile(strrep(File, '.dat', '_info.mat'));
% Load data:
data = mmData.Data.data;
% Calculate SeedPixel Correlation:
% Preserve data Aspect Ratio:
data_size = size(data);
if opts.imageResizeTo == -1
    dataAspectRatio = 1;
    xy_size = data_size([1 2]);
    A = data;
else
    dataAspectRatio = opts.imageResizeTo/max(data_size([1 2]));
    xy_size = round(data_size([1 2]).*dataAspectRatio);
    A = imresize3(data, [xy_size(1), xy_size(2), size(data,3)], 'nearest');
end
B = reshape(A, [], size(A,3))';
[CM, P] = corr(B);
% Apply Z Fisher transformation to corr Data:
if opts.FisherZ_transform
    CM = atanh(CM);
end
CM = reshape(CM, [xy_size(1) xy_size(2) xy_size(1)*xy_size(2)]);
P = single(reshape(P, [xy_size(1) xy_size(2) xy_size(1)*xy_size(2)]));
% SAVING DATA :
% Generate .DAT and .MAT file Paths:
datFile = fullfile(SaveFolder, Output);
% Create MetaData structure:
szCM = size(CM);
metaDat = struct('datName', {'CM', 'P'}, 'datSize', {szCM([1 2]), szCM([1 2])},...
    'datLength', {szCM(3) szCM(3)}, 'Datatype', {'single', 'single'}, 'datFile', datFile);
% Save CM and METADAT to DATFILE:
save2Dat(datFile, CM,'-w', metaDat)
% Append "P" to DATFILE:
save2Dat(datFile, P, '-a')
outFile = Output;
if ~isprop(metaData, 'BregmaXY')
    warning('MATLAB:UMIToolbox:MissingInfo', 'Anatomical landmarks not found in input File metaData!')
else
    metaData_out = matfile(strrep(datFile, '.dat', '_info.mat'));
    metaData_out.Properties.Writable = true;
    metaData_out.BregmaXY = metaData.BregmaXY*dataAspectRatio;
    metaData_out.LambdaXY = metaData.LambdaXY*dataAspectRatio;
end
end