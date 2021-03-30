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
Output = p.Results.Output;
%%%%
% Load Data:
mData = mapDatFile(File);
data = mData.Data.data;
% Calculate baseline
bsln = mean(data,3);
% Calculate DeltaF/F0:
data = (data - bsln)./ bsln;
%Save data using save2dat.m function
datFile = fullfile(SaveFolder, Output);
save2Dat(datFile, data);
% Output file names
outFile = Output;
end



