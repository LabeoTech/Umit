function outData = spatialGaussFilt(data, varargin)
% SPATIALGAUSSFILT performs a spatial gaussian filter in an image time
% series.


% Limitations:
% The data must be an Image time series OR a Correlation map with
% dimensions: {Y,X,T}.

% Defaults:
default_Output = 'spatFilt.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('Sigma',1);
opts_values = struct('Sigma',[0,Inf]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is a 3-D numerical matrix:
addOptional(p,'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, varargin{:});
%Initialize Variables:
outData = p.Results.data;
opts = p.Results.opts;
clear p
%%%%
% Apply gaussian filter:
disp('Filtering frames...')
outData = imgaussfilt(outData,opts.Sigma, 'FilterDomain','spatial');% Forced FilterDomain to "spatial" to avoid problems when the data has Infs or NaNs.);
disp('Finished with Spatial Gaussian filter.')      
end


