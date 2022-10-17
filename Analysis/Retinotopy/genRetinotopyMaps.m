function outFile = genRetinotopyMaps(data, metaData,SaveFolder, varargin)
% GENRETINOTOPYMAPS create amplitude and phase maps for each cardinal
% direction as well as the azimuth and elevation maps when at least two
% perpendicular directions are present in the input data. 
%
% Necessary info:
%   The metaData file must contain the variable "Stim_trialMarker" containing
%   the markers that identify the frames of each direction. In addition,
%   the file must contain the variables "eventID" and "eventNameList" with the 
%   indices and direction identifiers, respectively. The identifiers encoding one
%   or more of the cardinal angles are: 0,90,180 and 270.
%
% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T").
%   metaData: .mat file with meta data associated with "data".
%   opts (optional): structure containing extra parameters. See "default_opts" variable below for details!
%
% Outputs:
%   outFile: filenames of the generated data.

% Defaults:
default_Output = {'AzimuthMap.dat' 'ElevationMap.dat'};  %#ok. This line is here just for Pipeline management.
default_opts = struct('nSweeps', 1);
opts_values = struct('nSweeps', [1,Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is numerical
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, metaData,SaveFolder, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables and remove inputParser object:
data = p.Results.data;
metaData = p.Results.metaData;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
%%%%
% Further input validation:
% Check if the data is an image time series:
errID1 = 'umIToolbox:genRetinotopyMaps:WrongInput';
errID2 = 'umIToolbox:genRetinotopyMaps:MissingInput';
assert(all(ismember({'Y','X','T'},metaData.dim_names)), errID1,...
    'Input data must be an image time series with dimensions "Y","X","T".')
% Check if the metaData has the necessary variables:
assert(all(ismember({'Stim_trialMarker','eventID','eventNameList'},fieldnames(metaData))),errID2,...
    'The one or more of the variables "Stim_trialMarker", "eventID" and "eventNameList" are missing in meta data file.')
%%%%
outFile = {};
% Prepare data for FFT calculation:
% Find NaNs and replace them with zeros:
% idx_nan = isnan(data);
% data(idx_nan) = 0;

% Split data and calculate fft over time for each direction:
ampMaps = cell(size(metaData.eventNameList));
phiMaps = ampMaps;

fr_rise = find(metaData.Stim_trialMarker < .5 & [metaData.Stim_trialMarker(:,2:end) nan] > .5);
fr_fall = find(metaData.Stim_trialMarker > .5 & [metaData.Stim_trialMarker(:,2:end) nan] < .5);
freqFFT = round(opts.nSweeps)+1; 
w = waitbar(0,'Calculating FFT ...', 'Name', 'genRetinotopyMaps');
w.Children.Title.Interpreter = 'none';
for i = 1:numel(metaData.eventNameList)
    w.Children.Title.String = ['Calculating FFT for direction ' metaData.eventNameList{i}];
    % Calculate FFT:
    fDat = fft(data(:,:,fr_rise(i):fr_fall(i)),[],3);      
    % Create Amplitude/Phase maps for the input frequency (from Zhuang et al., 2017)  
    ampMaps{i} = (abs(fDat(:,:,freqFFT)) * 2) / size(fDat,3);
    phiMaps{i} = mod(-1.*angle(fDat(:,:,freqFFT)),2*pi);
    waitbar(i/numel(metaData.eventNameList),w);
end
clear fDat md
close(w);
%%% Calculate Azimuth and Elevation maps
% Check if all directions exist:
AzimMap = zeros([metaData.datSize 1],'single');
ElevMap = AzimMap;
[idxAz,indxAz] = ismember({'0','180'}, metaData.eventNameList);
[idxEl,indxEl] = ismember({'90','270'}, metaData.eventNameList);
% Azimuth:
if all(idxAz)
    disp('Calculating Azimuth map...')
    pow = mean(cat(3,ampMaps{indxAz}),3); % Average amplitude;
    phi = (phiMaps{indxAz(1)} - phiMaps{indxAz(2)})/2; % From Kalatsky and Stryker, 2003
    % Package amplitude and phase into complex numbers:
    AzimMap(:,:) = pow.*(cos(phi) + i.*sin(phi));     
end  
% Elevation:
if all(idxEl)  
    disp('Calculating Elevation map...')
    pow = mean(cat(3,ampMaps{indxEl}),3); % Average amplitude;
    phi = (phiMaps{indxEl(1)} - phiMaps{indxEl(2)})/2; % From Kalatsky and Stryker, 2003
    % Package amplitude and phase into complex numbers:
    ElevMap(:,:) = pow.*(cos(phi) + i.*sin(phi));     
end       
% Save raw FFT data:
%%% Save Maps:
% Save Azimuth map:
if sum(AzimMap(:)) ~= 0    
    md = genMetaData(AzimMap,metaData.dim_names([1 2]));
    filename = fullfile(SaveFolder, 'AzimuthMap.dat');
    save2Dat(filename,AzimMap,md);
    outFile = [outFile;{filename}];
end
% Save Elevation map:
if sum(ElevMap(:)) ~= 0    
    md = genMetaData(ElevMap,metaData.dim_names([1 2]));
    filename = fullfile(SaveFolder, 'ElevationMap.dat');
    save2Dat(filename,ElevMap,md);
    outFile = [outFile;{filename}];
end
disp('Done!');

end