function outFile = genRetinotopyMaps(data, metaData, varargin)
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
default_Output = {'FFT_0deg.dat', 'FFT_90deg.dat', 'FFT_180deg.dat', 'FFT_270deg.dat', 'AzimuthMap.dat' 'ElevationMap.dat'};  %#ok. This line is here just for Pipeline management.
default_opts = struct('nSweeps', 1, 'sigSource','Fluorophore', 'b_UseMask', 1, 'MaskFile', 'ImagingReferenceFrame.mat');
opts_values = struct('nSweeps', [1,Inf], 'sigSource',{{'Intrinsic', 'Fluorophore'}},'b_UseMask',[false,true],'MaskFile',{{'ImagingReferenceFrame.mat']});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is numerical
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Modality'));
% Parse inputs:
parse(p,data, metaData, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables and remove inputParser object:
data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
object = p.Results.object;
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
% Check if the logical mask exists in the mask file:
if opts.b_UseMask
    opts.MaskFile = findMyROIfile(opts.MaskFile,object);
    assert(isfile(opts.MaskFile), errID2,'Mask file not found!')
    a = load(opts.MaskFile);
    assert(isfield(a,'logical_mask'),errID2, 'Logical mask not found in Mask file!')
    assert(isequal(size(a.logical_mask), metaData.datSize), errID1,...
        'Logicak mask as a different size from the data!');
    logical_mask = a.logical_mask;
else
    logical_mask = true(metaData.datSize);
end
%%%%
% Prepare data for FFT calculation:
% Find NaNs and replace them with zeros:
idx_nan = isnan(data);
data(idx_nan) = 0;

% Split data and calculate fft over time for each direction:
fftArr = cell(size(metaData.eventNameList));
ampMaps = fftArr; phiMaps = fftArr;

fr_rise = find(data < thr & [data(:,2:end) nan] > thr);
fr_fall = find(data > thr & [data(:,2:end) nan] < thr);
nReps = sum(metaData.eventID == metaData.eventID(1)); % Get the number of triggers.
fr_rise = fr_rise(1:nReps:end);
fr_fall = fr_fall(1:nReps:end);
freqFFT= round(opts.nSweeps)+1; 
for i = unique(metaData.eventID)
    % Calculate FFT:
    fDat = fft(data(:,:,fr_rise(i):fr_fall(i)),[],3);
    fftArr{i} = fDat;
    % Create Amplitude/Phase maps for the input frequency (from Zhuang et al., 2017)  
    ampMaps{i} = (abs(fDat(:,:,freqFFT)) * 2) / size(fDat,3);
    phiMaps{i} = mod(-1*angle(fDat(:,:,freqFFT)),2*pi);
end
%%% Calculate Azimuth and Elevation maps
% Check if all directions exist:
AzimMap = [];
ElevMap = [];
[idxAz,indxAz] = ismember({'0','180'}, metaData.eventNameList);
[idxEl,indxEl] = ismember({'90','270'}, metaData.eventNameList);
% Azimuth:
if all(idxAz)
    AzimMap = zeros(metaData.datSize,2,'single');
    AzimMap(:,:,1) = mean(cat(3,ampMaps{indxAz}),3); % Average amplitude;
    if strcmpi(opts.sigSource, 'intrinsic')
        AzimMap(:,:,2) = phiMaps{indxAz(1)} - phiMaps{indxAz(2)}; % From Kalatsky and Stryker, 2003
    else
        AzimMap(:,:,2) = mean(cat(3,ampMaps{indxAz}),3); % From Zhuang et al., 2017
    end
end  
% Elevation:
if all(idxEl)
    ElevMap = zeros(metaData.datSize,2,'single');
    ElevMap(:,:,1) = mean(cat(3,ampMaps{indxEl}),3); % Average amplitude;
    if strcmpi(opts.sigSource, 'intrinsic')
        ElevMap(:,:,2) = phiMaps{indxEl(1)} - phiMaps{indxEl(2)}; % From Kalatsky and Stryker, 2003
    else
        ElevMap(:,:,2) = mean(cat(3,ampMaps{indxEl}),3); % From Zhuang et al., 2017
    end
end  
    
    




































































end