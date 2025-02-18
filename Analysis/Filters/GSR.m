function outData = GSR(data, varargin)
% GSR performs global signal regression to data in order to remove global fluctuations from signal.
% Limitations:
% The data must be an Image time series with dimensions Y,X,T.
%
% Inputs:
%   data: numerical matrix containing image time series.
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing extra parameters
%       - MaskFile (char): Name of the ".mat" file containing the
%       "logical_mask" variable. If "none", no mask will be applied.
%
% Outputs: 
%   outData: numerical matrix with dimensions {Y,X,T}.   

% Defaults:
default_Output = 'GSR.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('MaskFile','none');
opts_values = struct('MaskFile',{{'none'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data, varargin{:});
%Initialize Variables:
outData = data;
opts = p.Results.opts;
clear p data
%%%%

if ~strcmpi(opts.MaskFile,'none')
 %%%%% TBD %%%%%
else
    logical_mask = true([size(data,1),size(data,2)]);
end
% Find NaNs and replace them with zeros:
idx_nan = isnan(outData(:,:,1));
% Reshape data:
szData = size(outData);
outData = reshape(outData, [], szData(3));
outData(idx_nan(:),:) = 0;
mData = mean(outData,'all');
% Calculate GSR:
disp('Calculating Global signal regression...');
Sig = mean(outData(logical_mask(:),:),1);
Sig = Sig / mean(Sig);
X = [ones(szData(3),1), Sig'];
clear Sig
% 
A = zeros(length(X),size(outData,1),'single');
nChunks = calculateMaxChunkSize(outData,7);
indxChk = round(linspace(0,size(outData,1),nChunks));
if nChunks > 1    
    fprintf('\n%i%% ...\n', 0);
    for ii = 1:length(indxChk)-1
        fprintf('\n%11.0f%% ...\n', 100*ii/nChunks);
        A(:,indxChk(ii)+1:indxChk(ii+1)) = X*(X\outData(indxChk(ii)+1:indxChk(ii+1),:)');        
    end
    fprintf('\n%i%%\n', 100);
else
    A = X*(X\outData');
end
clear X 
outData = outData - A';% Center over constant mean value.
outData = reshape(outData,szData);
outData = outData + mData; % Center over constant mean value.
%
disp('Finished GSR.')
end