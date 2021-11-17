function outData = GSR(data, metaData)
% GSR performs global signal regression to data with path specified as
% ARGS.INPUT in order to remove global fluctuations from signal.


% Limitations:
% The data must be an Image time series with dimensions
% {Y,X,T}.

% Defaults:
default_Output = 'GSR.dat'; %#ok. This line is here just for Pipeline management.

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
% Parse inputs:
parse(p,data, metaData);
%Initialize Variables:
outData = p.Results.data;
metaData = p.Results.metaData;
clear p
%%%%
% Validate if "data" is an Image Time Series:
errID = 'umIToolbox:GSR:InvalidInput';
errMsg = 'Wrong Input Data type. Data must be an Image time series with dimensions "X", "Y" and "T".';
assert(all(ismember(metaData.dim_names,{'Y', 'X', 'T'})), errID, errMsg);

% Find NaNs and replace them with zeros:
idx_nan = isnan(outData);
outData(idx_nan) = 0;
% Reshape data:
szData = size(outData);
outData = reshape(outData, [], szData(3), 1);
% Calculate GSR:
disp('Calculating Global signal regression...')
mData = mean(outData,3);
Sig = mean(outData);
Sig = Sig / mean(Sig);
X = [ones(szData(3),1), Sig'];
B = X\outData';
A = X*B;
clear X B Sig
outData = outData - A';%Center at Zero
outData = outData + mData; %Center over constant mean value.
outData = reshape(outData,szData);
% Put NaNs back to data:
outData(idx_nan) = NaN;
disp('Finished GSR.')
end


