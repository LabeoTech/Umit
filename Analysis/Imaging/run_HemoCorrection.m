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
% The algorithm used here corresponds to the "Regression model" described
% in: 
% Valley, Matthew & Moore, Michael & Zhuang, Jun & Mesa, Natalia & Castelli, Dan & Sullivan,
% David & Reimers, Mark & Waters, Jack. (2020). Separation of hemodynamic signals from
% GCaMP fluorescence measured with widefield imaging. Journal of Neurophysiology. 123.
% 10.1152/jn.00304.2019.


% Defaults:
default_Output = 'hemoCorr_fluo.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('Red', true, 'Green', true, 'Amber', true);
opts_values = struct('Red',[false, true], 'Green',[false, true],'Amber',[false, true]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, data, metaData, varargin{:});

% Translate opts to char cell array:
fields = fieldnames(p.Results.opts);
idx = cellfun(@(x) p.Results.opts.(x), fields);
list = fields(idx)';

% Run HemoCorrection function from IOI library:
disp('Performing hemodynamic correction in fluo channel...')
outData = HemoCorrection(p.Results.SaveFolder, p.Results.data,false,list,'fMetaData',p.Results.metaData);

disp('Finished hemodynamic correction.')
end