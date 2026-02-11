function [outData, metaData] = calculate_response_amplitude(data, metaData, varargin)
% CALCULATE_RESPONSE_AMPLITUDE calculates signal amplitude for event-split data.
% Works on data with "E" dimension.
%
% Inputs:
%   data: 4D numeric matrix (E,Y,X,T) or filename of .dat file
%   metaData: struct or matfile with associated metadata
%   opts (optional): structure with fields:
%       preEvent_value: 'mean','median','min','max','AUC' (default: 'median')
%       postEvent_value: 'mean','median','min','max','AUC' (default: 'max')
%       TimeWindow_sec: [start,end] relative to trigger or 'all' (default: 'all')
%
% Outputs:
%   outData: 3D numerical matrix (E,Y,X) with response amplitudes
%   metaData: updated metadata structure

% Defaults
default_Output = 'amplitude_Map.dat'; %#ok
default_opts = struct('preEvent_value', 'median', 'postEvent_value', 'max', 'TimeWindow_sec', 'all');
opts_values = struct('preEvent_value', {{'mean','median','min','max','AUC'}},'postEvent_value', {{'mean','median','min','max','AUC'}}, 'TimeWindow_sec', {{'all',Inf}});

% Parse inputs
p = inputParser;
addRequired(p,'data', @(x) isnumeric(x) || ischar(x));
addRequired(p,'metaData', @(x) isstruct(x) || isa(x,'matlab.io.MatFile'));
addOptional(p,'opts', default_opts, @(x) isstruct(x));
parse(p, data, metaData, varargin{:});

data = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
clear p

% Validate options
validFcn = @(x) ismember(x, opts_values.preEvent_value);
assert(validFcn(opts.preEvent_value), 'Invalid preEvent_value');
assert(validFcn(opts.postEvent_value), 'Invalid postEvent_value');
if isnumeric(opts.TimeWindow_sec)
    assert(numel(opts.TimeWindow_sec)==2 && diff(opts.TimeWindow_sec)>0, ...
        'TimeWindow_sec must be a 2-element increasing vector');
else
    assert(strcmpi(opts.TimeWindow_sec,'all'), 'TimeWindow_sec must be "all" or numeric');
end

% Validate that data has "E" and "T"
dims = metaData.dim_names;
assert(all(ismember({'E','T'}, dims)), 'Data must have "E" and "T" dimensions');

% Determine data dimensions
bIsFile = ischar(data);
if bIsFile
    inFile = data;
    sz = [metaData.datSize, metaData.datLength];
    fidIn = fopen(inFile,'r');
    cIn = onCleanup(@() safeFclose(fidIn));
    Ne = metaData.datSize(1); Ny = metaData.datSize(2); Nx = metaData.datLength(1); Nt = metaData.datLength(2);
    outFile = fullfile(fileparts(inFile), 'CALCRESPONSEAMP.dat');
    % Set output file's meta data
    md = metaData;
    md.datSize = [Ne,Ny];
    md.datLength = Nx;
    md.dim_names = {'E','Y','X'};
    
    
else
    [Ne, Ny, Nx, Nt] = size(data);
end
outData = zeros(Ne, Ny, Nx, 'single');

% Determine trigger frame
trigFrame = round(metaData.preEventTime_sec*metaData.Freq);

% Determine post-event frames
if isnumeric(opts.TimeWindow_sec)
    frOn = round(opts.TimeWindow_sec(1)*metaData.Freq) + trigFrame;
    frOff = round(opts.TimeWindow_sec(2)*metaData.Freq) + trigFrame;
    if frOn > Nt || frOff > Nt || frOn>frOff
        warning('TimeWindow_sec out of range, using full post-event window');
        frOn = trigFrame+1; frOff = Nt;
    end
else
    frOn = trigFrame+1; frOff = Nt;
end

disp('Calculating response amplitude...');

% Loop over trials
for e = 1:Ne
    if bIsFile
        % Read trial from file                
        trialData = readTrial(fidIn,e,sz,'single');                
    else
        trialData = squeeze(data(e,:,:,:));
    end
    
    % Reshape to 2D (pixels x time)
    trial2D = reshape(trialData, Ny*Nx, Nt);
    
    % Baseline and post-trigger
    bslnData = trial2D(:,1:trigFrame);
    postData = trial2D(:,frOn:frOff);
    
    % Apply aggregation functions
    if strcmpi(opts.preEvent_value,'AUC') || strcmpi(opts.postEvent_value,'AUC')
        opts.preEvent_value = 'AUC';
        opts.postEvent_value = 'AUC';
    end
    
    bslnVal = applyAggFcn(bslnData, opts.preEvent_value);
    postVal = applyAggFcn(postData, opts.postEvent_value);
    
    % Calculate amplitude
    amp = postVal - bslnVal;
    amp = reshape(amp, Ny, Nx);
    outData(e,:,:) = amp;
    
end

if bIsFile
    % Write the whole Amp array to output file
    fidOut = fopen(outFile,'w');
    fwrite(fidOut,outData,'single');
    fclose(fidOut);
    fclose(fidIn);
    outData = outFile;
    metaData = md;
    save(strrep(outFile,'.dat','.mat'),'-struct','metaData');
    return
end

% Update metadata
new_dim_names = {'E','Y','X'};
metaData = genMetaData(outData, new_dim_names, metaData);
disp('Done!');

end

% -----------------------
function out = applyAggFcn(vals, fcn)
switch fcn
    case 'mean'
        out = mean(vals,2,'omitnan');
    case 'median'
        out = median(vals,2,'omitnan');
    case 'max'
        out = max(vals,[],2,'omitnan');
    case 'min'
        out = min(vals,[],2,'omitnan');
    case 'AUC'
        out = trapz(vals,2);
end
end
