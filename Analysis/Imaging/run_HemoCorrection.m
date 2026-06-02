function outData = run_HemoCorrection(SaveFolder,data, metaData, varargin)
% RUN_HEMOCORRECTION calls the function
% HEMOCORRECTION from the IOI library (LabeoTech).
% In brief, this function applies a pixelwise linear regression 
% of the fluorescence signal onto the reflectance signals (Valley et al,
% 2020).
%
% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T").
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing the Reflectance Channels to be
%   used in the correction.
%
% Note: 
% The algorithms used here are:
% Linear Regresssion: This one corresponds to the "Regression model" described
% in: 
%   Valley, Matthew & Moore, Michael & Zhuang, Jun & Mesa, Natalia & Castelli, Dan & Sullivan,
%   David & Reimers, Mark & Waters, Jack. (2020). Separation of hemodynamic signals from
%   GCaMP fluorescence measured with widefield imaging. Journal of Neurophysiology. 123.
%   10.1152/jn.00304.2019.
% Ratiometric:
%   Uses one reflectance/reference channel. If the reference channel has a
%   different number of frames than the fluorescence channel, it is temporally
%   resampled to the fluorescence timeline before subtraction. If it has a
%   higher sampling rate, an anti-aliasing low-pass filter is applied first.
%   Described in:
%   Wekselblatt, Joseph B., Erik D. Flister, Denise M. Piscopo, and Cristopher M. Niell. 2016. 
%   “Large-Scale Imaging of Cortical Dynamics during Sensory Perception and Behavior.”
%   Journal of Neurophysiology 115 (6): 2852–66. https://doi.org/10.1152/jn.01056.2015.


% Defaults:
default_Output = 'hemoCorr_fluo.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('Algorithm','LinearRegression','Red', true, 'Green', true, 'Amber', true,'Other','');
opts_values = struct('Algorithm',{{'LinearRegression','Ratiometric'}},'Red',[false, true], 'Green',[false, true],'Amber',[false, true],'Other',{{''}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, data, metaData, varargin{:});

fluoMetaData = p.Results.metaData;
if isa(fluoMetaData, 'matlab.io.MatFile')
    fluoMetaData = load(fluoMetaData.Properties.Source);
end

if( ~strcmp(SaveFolder(end),filesep) )
    SaveFolder = strcat(SaveFolder, filesep);
end

% Translate opts to char cell array:
list = {'Red','Green','Amber','Other'};
idx = cellfun(@(x) any(p.Results.opts.(x)), list);
list(~idx) = [];

% Run HemoCorrection function from IOI library:
fprintf('Performing hemodynamic correction in fluo channel using %s algorithm...\n',p.Results.opts.Algorithm)
if strcmpi(p.Results.opts.Algorithm,'LinearRegression')
    if ~isempty(p.Results.opts.Other)
        warning('Other files will be ignored with Linear Regression algorithm. Only RED, GREEN and AMBER are allowed.')
    end
    outData = HemoCorrection(SaveFolder,data,false,list,'fMetaData',fluoMetaData);
else 
    assert(length(list)==1, 'Failed to perform Ratiometric hemodynamic correction! Only a single channel is allowed!')
    % Control for single file in list or OTHER.
    fprintf('Using channel "%s" in hemodynamic correction...\n',list{1});
    if strcmpi(list{1}, 'Other')
        refFile = p.Results.opts.Other;
        [~,refName,refExt] = fileparts(refFile);
        if isempty(refExt)
            refFile = [refName '.dat'];
        else
            refFile = [refName refExt];
        end
    else
        switch lower(list{1})
            case 'red'
                if( exist([SaveFolder 'rChan.dat'], 'file') )
                    refFile = 'rChan.dat';
                else
                    refFile = 'red.dat';
                end
            case {'amber', 'yellow'}
                if( exist([SaveFolder 'yChan.dat'], 'file') )
                    refFile = 'yChan.dat';
                else
                    refFile = 'yellow.dat';
                end
            case 'green'
                if( exist([SaveFolder 'gChan.dat'], 'file') )
                    refFile = 'gChan.dat';
                else
                    refFile = 'green.dat';
                end
        end
    end

    [~,refName,~] = fileparts(refFile);
    refMatFile = fullfile(SaveFolder, [refName '.mat']);
    assert(exist(refMatFile, 'file') == 2, ...
        'Failed to perform Ratiometric hemodynamic correction! Missing metadata file "%s".', refMatFile);
    refMetaData = load(refMatFile);

    % Load reference channel and resample it to the fluorescence timeline
    % before normalization/subtraction when frame counts differ.
    tmp = loadDatFile(fullfile(SaveFolder,refFile));
    assert(isequal(size(tmp,1), size(data,1)) && isequal(size(tmp,2), size(data,2)), ...
        'Reference and fluorescence channels must have the same spatial size')
    tmp = iResampleReferenceToFluoTimeline(tmp, refMetaData, fluoMetaData);
    assert(all(size(data) == size(tmp)), ...
        'Reference and fluorescence channels must have the same size after temporal resampling')

    % Perform normalization
    datSize = size(data);
    data = reshape(data,[],datSize(end));
    m_data = mean(data,2);
    % Calculate dF/Fmean for fluo channel
    data = (data - m_data)./m_data;    
    
    tmpSize = size(tmp);
    % Calculate dF/Fmean for reference channel
    tmp = reshape(tmp,[],tmpSize(end));
    tmp = (tmp- mean(tmp,2))./mean(tmp,2);
    % Subtract fluo by reference to perform hemodynamic correction
    data = data - tmp;
    % Put original average fluorescence back to the normalized data:
    data = (data.*m_data) + m_data;
    clear tmp;
    % Reshape the corrected fluo channel to its original shape
    outData = reshape(data,datSize);
end
disp('Finished hemodynamic correction.')
end

function tmp = iResampleReferenceToFluoTimeline(tmp, refMetaData, fMetaData)
%IRESAMPLEREFERENCETOFLUOTIMELINE Match one reference channel to fluo T.
%
%   tmp = iResampleReferenceToFluoTimeline(tmp, refMetaData, fMetaData)
%   applies anti-aliasing before downsampling when the reference channel has
%   a higher sampling rate than the fluorescence channel, then resamples tmp
%   along the third dimension to fMetaData.datLength.

NtRef = refMetaData.datLength(1,end);
NtFluo = fMetaData.datLength(1,end);
freqRef = refMetaData.Freq;
freqFluo = fMetaData.Freq;

if( size(tmp,3) ~= NtRef )
    error('Input reference channel length does not match its metadata.');
end

% Resampling is valid only when both channels span the same recording
% duration. A mismatched duration usually indicates that one channel was
% cropped, truncated, or no longer matches its metadata.
durationTolSec = 1e-3;
refDurationSec = double(NtRef) ./ double(freqRef);
fluoDurationSec = double(NtFluo) ./ double(freqFluo);
assert(abs(refDurationSec - fluoDurationSec) <= durationTolSec, ...
    'Umitoolbox:run_HemoCorrection:DurationMismatch', ...
    ['Reference channel does not span the same recording duration as the ' ...
     'fluorescence channel. Reference: Length=%d, FrameRateHz=%g, ' ...
     'Duration=%0.6f s. Fluorescence: Length=%d, FrameRateHz=%g, ' ...
     'Duration=%0.6f s.'], ...
    NtRef, freqRef, refDurationSec, NtFluo, freqFluo, fluoDurationSec);

if( NtRef == NtFluo )
    tmp = single(tmp);
    return
end

if( freqRef > freqFluo && NtRef > NtFluo )
    aaCutoff = 0.45 * freqFluo;
    if( aaCutoff > 0 && aaCutoff < freqRef/2 )
        sz = size(tmp);
        tmp = reshape(tmp, [], sz(3));
        f = fdesign.lowpass('N,F3dB', 4, aaCutoff, freqRef);
        lpass = design(f, 'butter');
        tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
        tmp = reshape(tmp, sz);
    end
end

sz = size(tmp);
xRef = linspace(0, 1, NtRef);
xFluo = linspace(0, 1, NtFluo);
tmp = reshape(tmp, [], NtRef);
tmp = interp1(xRef, single(tmp)', xFluo, 'linear', 'extrap')';
tmp = reshape(single(tmp), sz(1), sz(2), NtFluo);

end
