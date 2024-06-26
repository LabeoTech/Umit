function outData = normalizeZScore(data)
% NORMALIZEZSCORE normalizes an image time series with zero mean and unit standard 
% deviation.

% Limitations:
% The data must be an Image time series with dimensions
% {Y,X,T}.

% Defaults:
default_Output = 'normZ.dat'; %#ok. This line is here just for Pipeline management.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
% Parse inputs:
parse(p,data);
%Initialize Variables:
outData = p.Results.data;
clear p
%%%%

disp('Calculating zscore...')
orig_sz = size(outData);
outData = reshape(outData,[],orig_sz(3));
outData = (outData - mean(outData,2,'omitnan'))./std(outData,0,2,'omitnan');
outData = reshape(outData,orig_sz);
end
    