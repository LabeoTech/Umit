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
% Find NaNs and replace them with zeros:
idx_nan = isnan(outData);
outData(idx_nan) = 0;
% Reshape data:
szData = size(outData);
outData = reshape(outData, [], szData(3));
% Calculate GSR:
disp('Calculating Global signal regression...');
mData = mean(outData,3);
Sig = mean(outData(logical_mask(:),:),1);
Sig = Sig / mean(Sig);
X = [ones(szData(3),1), Sig'];
clear Sig
A = X*(X\outData');
clear X 
outData = outData - A';%Center at Zero
outData = outData + mData; %Center over constant mean value.
outData = reshape(outData,szData);
% Put NaNs back to data:
outData(idx_nan) = NaN;
disp('Finished GSR.')
end


