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
default_opts = struct('imageResizeTo', 64);
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
imsz = opts.imageResizeTo;

% Open memMapfile:
mmData = mapDatFile(File);
% Load data:
data = mmData.Data.data;
% Calculate SeedPixel Correlation:
A = imresize3(data, [imsz, imsz, size(data,3)]);
B = reshape(A, [], size(A,3))';
[CM, P] = corr(B);
CM = reshape(CM, [imsz imsz imsz^2]);
P = single(reshape(P, [imsz imsz imsz^2]));
clear A B data

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
end


