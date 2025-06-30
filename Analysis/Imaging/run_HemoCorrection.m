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
    outData = HemoCorrection(SaveFolder,data,false,list,'fMetaData',p.Results.metaData);
else 
    assert(length(list)==1, 'Failed to perform Ratiometric hemodynamic correction! Only a single channel is allowed!')
    % Control for single file in list or OTHER 
    fprintf('Using channel "%s" in hemodynamic correction...\n',list{1});
    refFile = list{1};
    [~,refFile,~] = fileparts(refFile); refFile = [refFile '.dat'];
    % Load Reference channel
    tmp = loadDatFile(fullfile(SaveFolder,refFile));
    % Control for equal dimensions
    assert(all(size(data) == size(tmp)),'Reference and fluorescence channels must have the same size')    
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