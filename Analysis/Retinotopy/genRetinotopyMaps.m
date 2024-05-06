function outFile = genRetinotopyMaps(data, metaData,SaveFolder, varargin)
% GENRETINOTOPYMAPS create amplitude and phase maps for each cardinal
% direction as well as the azimuth and elevation maps when at least two
% perpendicular directions are present in the input data.
%
% The timestamps marking the individual directions must be encoded in the "events.mat" file.

% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T").
%   metaData: .mat file with meta data associated with "data".
%   opts (optional): structure containing extra parameters. See "default_opts" variable below for details!
%
% Outputs:
%   outFile: filenames of the generated data.

% Defaults:
default_Output = {'AzimuthMap.dat' 'ElevationMap.dat'};  %#ok. This line is here just for Pipeline management.
default_opts = struct('nSweeps', 1, 'b_useAverageMovie', false, 'ViewingDist_cm', 0,'ScreenXsize_cm',0,'ScreenYsize_cm',0);
opts_values = struct('nSweeps', [1,Inf], 'b_useAverageMovie', [true,false], 'ViewingDist_cm', [0,Inf],'ScreenXsize_cm',[0,Inf],'ScreenYsize_cm',[0,Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
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
% Check if the "events.mat" file exists and if the event Names correspond to the directions 0, 90, 180 and 270:
assert(isfile(fullfile(SaveFolder, 'events.mat')), errID2, '"events.mat" file not found!')
evntInfo = load(fullfile(SaveFolder, 'events.mat'));
assert(all(ismember(evntInfo.eventNameList, {'0','90','180','270'})),...
    errID2, 'Invalid directions found in "events.mat" file! They must be one of the following: 0,90,180,270.')
%%%%
outFile = {};
% Prepare data for FFT calculation:
% Find NaNs and replace them with zeros:
% idx_nan = isnan(data);
% data(idx_nan) = 0;

% Split data and calculate fft over time for each direction:
ampMaps = cell(size(evntInfo.eventNameList));
phiMaps = ampMaps;
% Select Frequency:
if opts.b_useAverageMovie
    freqFFT = 2;
else
    freqFFT = round(opts.nSweeps)+1;
end

w = waitbar(0,'Calculating FFT ...', 'Name', 'genRetinotopyMaps');
w.Children.Title.Interpreter = 'none';
framestamps = round(evntInfo.timestamps*metaData.Freq);
for ind = 1:numel(evntInfo.eventNameList)
    w.Children.Title.String = ['Calculating FFT for direction ' evntInfo.eventNameList{ind}]; drawnow
    indxOn = find(evntInfo.eventID == ind & evntInfo.state == 1);
    indxOff = find(evntInfo.eventID == ind & evntInfo.state == 0);
    if opts.b_useAverageMovie        
        % Create average DeltaR movie for each direction:
        bsln_len = round(mean(framestamps(indxOn(2:end)) - framestamps(indxOff(1:end-1))));
        if isempty(bsln_len) || bsln_len <= 0
            error('Could not compute a valid baseline period! Does the recording have interstimulus time? If not, set the option "b_useAverageMovie" to FALSE.');
        end
        trial_len = round(mean(framestamps(indxOff) - framestamps(indxOn)));
        avg_mov = zeros([metaData.datSize, trial_len],'single');
        disp(['Creating average Delta R movie for direction ' evntInfo.eventNameList{ind}]);
        for ii = 1:length(indxOn)
            DeltaR = data(:,:,framestamps(indxOn(ii)): framestamps(indxOn(ii)) + trial_len -1) - ...
                median(data(:,:,framestamps(indxOn(ii)) - bsln_len: framestamps(indxOn(ii)) - 1), 3,'omitnan'); % trial period minus the intertrial period.
            avg_mov = [avg_mov + DeltaR];
        end
        avg_mov = avg_mov/length(indxOn); % Average DeltaR movie.
        % Calculate FFT of average movie:
        fDat = fft(avg_mov,[],3);
    else        
        % Calculate FFT of whole stimulus:
        fDat = fft(data(:,:,framestamps(indxOn(1)):framestamps(indxOff(end))),[],3);        
    end
    % Create Amplitude/Phase maps for the input frequency (from Zhuang et al., 2017)
    ampMaps{ind} = (abs(fDat(:,:,freqFFT)) * 2) / size(fDat,3);
    phiMaps{ind} = mod(-1.*angle(fDat(:,:,freqFFT)),2*pi);
    waitbar(ind/numel(evntInfo.eventNameList),w);
end
clear fDat md avg_mov 
close(w);
%%% Calculate Azimuth and Elevation maps
% Check if all directions exist:
AzimMap = zeros([metaData.datSize 2],'single');
ElevMap = AzimMap;
[idxAz,indxAz] = ismember({'0','180'}, evntInfo.eventNameList);
[idxEl,indxEl] = ismember({'90','270'}, evntInfo.eventNameList);
% Azimuth:
if all(idxAz)
    disp('Calculating Azimuth map...')
    AzimMap(:,:,1) = mean(cat(3,ampMaps{indxAz}),3); % Average amplitude;
    AzimMap(:,:,2) = pi + ((phiMaps{indxAz(1)} - phiMaps{indxAz(2)})/2); % From Kalatsky and Stryker, 2003; Shift phase to be between 0 - 2pi  
end
% Elevation:
if all(idxEl)
    disp('Calculating Elevation map...')
    ElevMap(:,:,1) = mean(cat(3,ampMaps{indxEl}),3); % Average amplitude;
    ElevMap(:,:,2) = pi + ((phiMaps{indxEl(1)} - phiMaps{indxEl(2)})/2); % From Kalatsky and Stryker, 2003; Shift phase to be between 0 - 2pi      
end
% Rescale phase maps:
if all([opts.ViewingDist_cm, opts.ScreenXsize_cm, opts.ScreenYsize_cm])>0
    disp('Rescaling phase maps to visual angle...')
    va_az = atand(opts.ScreenXsize_cm/(2*opts.ViewingDist_cm)); va_az = [-va_az va_az];
    va_el = atand(opts.ScreenYsize_cm/(2*opts.ViewingDist_cm)); va_el = [-va_el va_el];
    AzimMap(:,:,2) = rescale(AzimMap(:,:,2), va_az(1), va_az(2));
    ElevMap(:,:,2) = rescale(ElevMap(:,:,2), va_el(1), va_el(2));    
elseif any([opts.ViewingDist_cm, opts.ScreenXsize_cm, opts.ScreenYsize_cm])>0
    warning('Cannot calculate visual angle! Phase values will be shown as radians!')
end
%%% Save Maps:
% Save Azimuth map:
if sum(AzimMap(:)) ~= 0
    md = genMetaData(AzimMap,{'Y','X','T'}, metaData);
    filename = fullfile(SaveFolder, 'AzimuthMap.dat');
    save2Dat(filename,AzimMap,md);
    outFile = [outFile;{filename}];
end
% Save Elevation map:
if sum(ElevMap(:)) ~= 0
    md = genMetaData(ElevMap,{'Y','X','T'}, metaData);
    filename = fullfile(SaveFolder, 'ElevationMap.dat');
    save2Dat(filename,ElevMap,md);
    outFile = [outFile;{filename}];
end
disp('Done!');

end