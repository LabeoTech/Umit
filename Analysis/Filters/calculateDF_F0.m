function outFile = calculateDF_F0(File, SaveFolder, varargin)
% CALCULATEDF_F0 generates DeltaF/F0 movie from FILE.
% Inputs:
%   File: fullpath of functional imaging .DAT file.
%   SaveFolder: path to save the output file.
%   Output (optional) : Name of outFile.
% Output:
%   outFile: name of Output file.
% Defaults:
default_Output = 'deltaF_F0.dat'; % This is here only as a reference for PIPELINEMANAGER.m. 
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Output file:
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x));
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%Initialize Variables:
File = p.Results.File; 
SaveFolder = p.Results.SaveFolder;
%%%%

% Load Data:
[mData, metaData] = mapDatFile(File);
data = mData.Data.data;
% Find "T"ime dimension:
indxT = find(strcmp('T', metaData.dim_names));
% Calculate baseline over time
bsln = mean(data,indxT);
% Calculate DeltaF/F0:
data = (data - bsln)./ bsln;
% % Replace NaNs with zeros:
% idx = isnan(data);
% data(idx) = 0;
% SAVING DATA :
% Generate .DAT and .MAT file Paths:
[~,filename,ext] = fileparts(File);
outFile = ['deltaF_F0_' filename ext];
save2Dat(fullfile(SaveFolder, outFile), data, metaData.dim_names);
end



