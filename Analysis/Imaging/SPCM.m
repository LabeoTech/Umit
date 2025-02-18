function outData = SPCM(data, varargin)
% This function creates a Seed-Pixel Correlation Map (SPCM) by calculating
% Pearson's zero-lag temporal correlation between pixels.
% Note that this function calculates the correlation for all pixels,
% which may use a significant amount of RAM if the data is large.
%
% Note:
% If the goal is to create a correlation matrix for specific regions, consider
% using the "genCorrelationMatrix" function instead.
%
% Inputs:
%   data (3D numerical matrix): Image time series with dimensions {'Y', 'X', 'T'}.
%   FisherZ_transform (bool, optional): If true, Fisher's Z transformation is
%       applied to the data. Otherwise, the data is expressed as Pearson's
%       correlation values.
%
% Outputs:
%   outData (struct): Structure containing Seed-Pixel correlation maps.


% Defaults:
default_Output = 'SPCMap.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('b_FisherZ_transform', false);
opts_values = struct('b_FisherZ_transform',[false, true]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, varargin{:});
% Initialize Variables:
opts = p.Results.opts;
clear p
%%%%

% Calculate SeedPixel Correlation:
sz_dat = size(data);
data = reshape(data, [], size(data,3))';
disp('Calculating Pearson''s correlation...');
% Calculate Pearson's correlation:
spcmap = corrcoef(data);
clear data
% Apply Z Fisher transformation to corr Data:
if opts.b_FisherZ_transform
    disp('Applying Fisher''s Z transform...');    
    % Z-Fisher Transform:
    spcmap = atanh(spcmap);   
end
spcmap = reshape(spcmap, [sz_dat(1) sz_dat(2) sz_dat(1)*sz_dat(2)]);
% Packlage the data in structure
outData = genDataStructure(spcmap);
disp('Finished SPCM');
end