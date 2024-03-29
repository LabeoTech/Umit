function [outData, metaData] = genVSM(SaveFolder, varargin) 
% GENVSM creates a Visual Sign Map (Sereno et al. 1994,1995; Zhuang et al. 2017)
% from the phase component of the Azimuth and Elevation maps created with
% the function "genRetinotopyMaps.m".
% Inputs:
%   SaveFolder (char): full path of the folder containing the files
%    "AzimuthMap.dat" and "ElevationMap.dat".
%   opts (optional): structure containing extra parameters. See "default_opts" variable below for details!
%
% Outputs:
%   data (2D matrix): matrix containing the Visual Sign Map (VSM).
%   metaData (struct): structure containing the meta data associated with
%   "data".

% Defaults:
default_Output = 'visualSignMap.dat'; %#ok. This line is here just for Pipeline management.
default_opts = struct('PhaseMapFilter_Sigma', 0,'VSMFilter_Sigma',0, 'b_CreatePatches', false, 'PatchFileName', 'VisualCtxAreas.mat', 'b_UseMask', false, 'MaskFile','ImagingReferenceFrame.mat');
opts_values = struct('PhaseMapFilter_Sigma', [0, Inf], 'VSMFilter_Sigma', [0, Inf], 'b_CreatePatches', [true; false], 'PatchFileName', {{'VisualCtxAreas.mat'}},'b_UseMask', [true,false], 'MaskFile',{{'ImagingReferenceFrame.mat'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
% Notes on Spatial filter:
% - A value of zero indicates that no filters will apply!
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Modality'));
% Parse inputs:
parse(p,SaveFolder, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables and remove inputParser object:
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
object = p.Results.object;
clear p
%%%%
% Further input validation:
% Check if the save folder contains the files "Azimuth.dat" and
% "Elevation.dat"
AzFile = fullfile(SaveFolder, 'AzimuthMap.dat');
ElFile = fullfile(SaveFolder, 'ElevationMap.dat');
errID = 'umIToolbox:genVSM:MissingInput';
errMsg = 'The files "AzimuthMap.dat" and "ElevationMap.dat" were not found in the Save Folder!';
assert(isfile(AzFile) & isfile(ElFile),errID, errMsg);
% Open Maps
[azMap, metaData]= loadDatFile(AzFile);
elMap = loadDatFile(ElFile);
% Calculate Visual Sign Map:
phaseAz = azMap(:,:,2); 
phaseEl = elMap(:,:,2); 
% Remove NaNs from phase maps:
phaseAz(isnan(phaseAz)) = 1000;
phaseEl(isnan(phaseAz)) = 1000;
% Filter Phase masps
if opts.PhaseMapFilter_Sigma > 0
    disp('Filtering phase maps...')
    phaseAz = imgaussfilt(phaseAz, opts.PhaseMapFilter_Sigma);
    phaseEl = imgaussfilt(phaseEl, opts.PhaseMapFilter_Sigma);
end
disp('Calculating visual sign map...');
[gradAzx, gradAzy] = gradient(phaseAz);
[gradElx, gradEly] = gradient(phaseEl);
gradDirAz = atan2(gradAzy, gradAzx);
gradDirEl = atan2(gradEly, gradElx);

outData = sin(angle(exp(1i.*gradDirAz).*exp(-1i.*gradDirEl)));
outData(isnan(outData)) = 0;
% Filter Visual Sign Map
if opts.VSMFilter_Sigma > 0
    disp('Filtering phase maps...')
    outData = imgaussfilt(outData, opts.VSMFilter_Sigma);    
end
% Create meta data:
metaData = genMetaData(outData,{'Y','X'}, metaData);
disp('Finished creating VSM');
% Segment the Visual Sign Map and create a labeled matrix:
if opts.b_CreatePatches
    genVAmask(outData, metaData, SaveFolder, opts, object)
end
end

% Local function

function genVAmask(vsm, metaData, SaveFolder, opts, object)
% GENVAMASK generates a .mat file with the ROIs from the segmented 
% Visual Sign Map created using the function "genVSM.m".
% Inputs:
%   data: numerical matrix containing the visual sign map (with dimensions "Y", "X").
%   metaData: .mat file with meta data associated with "data".
%   SaveFolder (char): full path to the folder containing the "VSM.dat"
%       file. The ROImasks file will be saved here as well.
%   opts (optional) : structure containing extra parameters:
%       - fileName (char): Name of the ".mat" file to save the masks.
%       - b_UseMask (bool): If TRUE, the function loads a logical mask from
%           the " MaskFile to create the ROIs.
%       - MaskFile (char): Name of the ".mat" file containing the
%          "logical_mask" variable.
%   object (used by PipelineManager class ONLY): Protocol objecf of type "Modality".
%       If provided, the function will look for the "MaskFile" in the "Subject" folder.
% Output: 
%   No output. The function will save a "ROImasks" file in the SaveFolder.
%   This file is used by the ROImanager app.
            
%%% Arguments parsing and validation %%%

if ~endsWith(opts.PatchFileName, '.mat')
    opts.PatchFileName = [opts.PatchFileName, '.mat'];
end
    
% Validate if "data" is an Image Time Series:
errID = 'umIToolbox:genVAmask:InvalidInput';
errMsg = 'Wrong Input Data type. Data must be an Image with dimensions "Y", "X".';
assert(all(ismember(metaData.dim_names,{'Y', 'X'})), errID, errMsg);
if opts.b_UseMask
    if isempty(fileparts(opts.MaskFile))
        % Update path for mask file:
        opts.MaskFile = fullfile(SaveFolder, opts.MaskFile);
    end
    % Retrieve file containing logical mask
    opts.MaskFile = findMyROIfile(opts.MaskFile,object);
    a = load(opts.MaskFile, 'logical_mask');
    if isempty(fieldnames(a))
        msg = 'Variable "logical_mask" not found in mask file!';
        msgID = 'umIToolbox:genVAmask:MissingInput';
        error(msgID,msg);
    end
    disp(['Using logical mask from file: ' opts.MaskFile ]);
    logical_mask = a.logical_mask;
    % Check if the logical mask has the same frame size than the data:
    assert(isequal(size(logical_mask), metaData.datSize), 'umIToolbox:genVAmask:InvalidInput',...
        'The logical mask size is different from the frame size in data');
else
    logical_mask = true(metaData.datSize);
end    
% Segment the Visual Sign Map:
% Set threshold as +-1.5 STD of the VSM ( Garret et al. 2014):
thr = 1.5*std(vsm(:));
% First, binarize the VSM and intersect with logical mask:
rawMap = (imbinarize(abs(vsm), thr) & logical_mask);
% Perform opening of Binary Image
rawMap = bwmorph(rawMap,'open',inf);
[L,nPatch] = bwlabel(rawMap);
% patchSign = zeros(size(nPatch));
% Close patches:
patchMap = zeros(size(rawMap),'single');
for i = 1:nPatch 
%     if mean(vsm(L == i),'all') > 0
%         patchSign(i) = 1;
%     elseif mean(vsm(L ==i),'all') < 0
%         patchSign(i) = -1;
%     end
    patchMap = patchMap + bwmorph(imfill(L == i, 'holes'), 'close', Inf);
end
total_area = bwmorph(imfill(~imerode(~patchMap, strel('diamond',15)),...
    'holes'), 'majority', Inf);
% Estimate patch borders:
patchBorder = total_area - patchMap;
patchBorder = bwmorph(patchBorder, 'skel', Inf);
patchBorder = bwmorph(patchBorder, 'spur',Inf); % Remove small edges from the border.
patches = imfill(patchBorder,'holes') & ~patchBorder;
patches = bwlabel(patches,4);
% Save patches to ".mat" file:
filename = fullfile(SaveFolder, opts.PatchFileName);
save(filename, 'patches');
disp(['Patches saved as ' filename]);
end