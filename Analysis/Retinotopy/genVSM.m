function [data, metaData] = genVSM(SaveFolder, varargin)
% GENVSM creates a Visual Sign Map (Sereno et al. 1994, Sereno et al. 1995)
% from the phase component of the Azimuth and Elevation maps created with
% the function "genRetinotopyMaps.m".
% Inputs:
%   SaveFolder (char): full path of the folder containing the files
%    "AzimuthMap.dat" and "ElevationMap.dat".
%   opts (optional): structure containing extra parameters. See "default_opts" variable below for details!
%
% Outputs:
%   data (2D matrix): matrix containing the Visual Sign Map (VSM).
%   metaData (struct): structure containing the meta data associated with
%   "data".

% Defaults:
default_Output = 'VSM.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('SpatialFilter_Sigma', 0);
opts_values = struct('SpatialFilter_Sigma', [0, Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
% Notes on Spatial filter:
% - A value of zero indicates that no filters will apply!
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables and remove inputParser object:
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
%%%%
% Further input validation:
% Check if the save folder contains the files "Azimuth.dat" and
% "Elevation.dat"
AzFile = fullfile(SaveFolder, 'AzimuthMap.dat');
ElFile = fullfile(SaveFolder, 'ElevationMap.dat');
errID = 'umIToolbox:genVSM:MissingInput';
errMsg = 'The files "AzimuthMap.dat" and "ElevationMap.dat" were not found in the Save Folder!';
assert(isfile(AzFile) & isfile(ElFile),errID, errMsg);
% Open Maps
azMap= loadDatFile(AzFile);
elMap = loadDatFile(ElFile);
% Calculate Visual Sign Map:
phaseAz = azMap(:,:,2); 
phaseEl = elMap(:,:,2); 
% Remove NaNs from phase maps:
phaseAz(isnan(phaseAz)) = 1000;
phaseEl(isnan(phaseAz)) = 1000;
disp('Calculating visual sign map...');
[gradAzx, gradAzy] = gradient(phaseAz);
[gradElx, gradEly] = gradient(phaseEl);
gradDirAz = atan2(gradAzy, gradAzx);
gradDirEl = atan2(gradEly, gradElx);

data = sin(angle(exp(1i.*gradDirAz).*exp(-1i.*gradDirEl)));
data(isnan(data)) = 0;
% Filter map
if opts.SpatialFilter_Sigma > 0
    data = imgaussfilt(data, opts.SpatialFilter_Sigma);
end
% Create meta data:
metaData = genMetaData(data,{'Y','X'});
disp('Done');
end
