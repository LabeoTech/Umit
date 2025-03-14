function outData = GSR(data, metaData, varargin)
% GSR performs global signal regression to data in order to remove global fluctuations from signal.
% Limitations:
% The data must be an Image time series with dimensions
% {Y,X,T}.
%
% Inputs:
%   data: numerical matrix containing image time series (with dimensions "Y", "X", "T").
%   metaData: .mat file with meta data associated with "data".
%   opts (optional) : structure containing extra parameters
%       - b_UseMask (bool): If TRUE, the function loads a logical mask from
%           the " MaskFile to calculate the Global Signal only from pixels inside the mask.
%       - MaskFile (char): Name of the ".mat" file containing the "logical_mask" variable.
%   object (used by PipelineManager class ONLY): Protocol objecf of type "Modality".
%       If provided, the function will look for the "MaskFile" in the "Subject" folder.
%
% Outputs: 
%   outData: numerical matrix with dimensions {Y,X,T}.   
%   metaData: .mat file with meta data associated with "outData".

% Defaults:
default_Output = 'GSR.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('b_UseMask', false, 'MaskFile','ImagingReferenceFrame.mat');
opts_values = struct('b_UseMask', [true,false], 'MaskFile',{{'ImagingReferenceFrame.mat'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p,'data',@(x) isnumeric(x) & ndims(x) == 3); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality'));

% Parse inputs:
parse(p,data, metaData, varargin{:});
%Initialize Variables:
outData = p.Results.data;
metaData = p.Results.metaData;
opts = p.Results.opts;
object = p.Results.object;
clear p data
%%%%
% Validate if "data" is an Image Time Series:
errID = 'umIToolbox:GSR:InvalidInput';
errMsg = 'Wrong Input Data type. Data must be an Image time series with dimensions "X", "Y" and "T".';
assert(all(ismember(metaData.dim_names,{'Y', 'X', 'T'})), errID, errMsg);
if opts.b_UseMask
    % Retrieve file containing logical mask
    opts.MaskFile = findMyROIfile(opts.MaskFile,object);
    a = load(opts.MaskFile, 'logical_mask');
    if isempty(fieldnames(a))
        msg = 'Variable "logical_mask" not found in mask file!';
        msgID = 'umIToolbox:GSR:MissingInput';
        error(msgID,msg);
    end
    disp(['Using logical mask from file: ' opts.MaskFile ]);
    logical_mask = a.logical_mask;
    % Check if the logical mask has the same frame size than the data:
    assert(isequal(size(logical_mask), metaData.datSize), 'umIToolbox:GSR:InvalidInput',...
        'The logical mask size is different from the frame size in data');
else
    logical_mask = true(metaData.datSize);
end


% Reshape data:
szData = size(outData);
% Find NaNs:
idx_nan = isnan(outData(:,:,1));idx_nan = idx_nan(:);
outData = reshape(outData, [], szData(3));
% Calculate GSR:
disp('Calculating Global signal regression...');
Sig = mean(outData(logical_mask(:),:),1,'omitnan');
Sig = Sig / mean(Sig);
X = [ones(szData(3),1), Sig'];
clear Sig
% 
A = zeros(length(X),size(outData,1),'single');
nChunks = calculateMaxChunkSize(outData,7);
indxChk = round(linspace(0,size(outData,1),nChunks));
mData = 0;
if nChunks > 1    
    fprintf('\n%i%% ...\n', 0);
    for ii = 1:length(indxChk)-1
        fprintf('%11.0f%% ...\n', 100*ii/nChunks);
        nan_msk = idx_nan(indxChk(ii)+1:indxChk(ii+1));
        dataChunk = outData(indxChk(ii)+1:indxChk(ii+1),:);
        dataChunk(nan_msk,:) = 0;
%         A(:,indxChk(ii)+1:indxChk(ii+1)) = X*(X\outData(indxChk(ii)+1:indxChk(ii+1),:)');        
        A(:,indxChk(ii)+1:indxChk(ii+1)) = X*(X\dataChunk');        
        mData = mData + sum(dataChunk,'all');
    end
    mData = mData/numel(outData);
    clear dataChunk
    fprintf('%i%%\n', 100);
else
    mData = mean(outData,'all');
    A = X*(X\outData');
end
clear X 
outData = outData - A';% Center over constant mean value.
clear A
% Put NaNs back to data:
outData(idx_nan(:),:) = NaN;
outData = reshape(outData,szData);
outData = outData + mData; % Center over constant mean value.
disp('Finished GSR.')
end