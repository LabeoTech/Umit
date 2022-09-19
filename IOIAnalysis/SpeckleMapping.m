function DatOut = SpeckleMapping(folderPath, sType, channel, bSaveMap, bLogScale)
%%%%%%%%%%%%%%%%%%%% Speckle Mapping function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the standard deviation (spatialy or temporaly) of speckle
% acquisition. This measure is proportional to the strength of blood flow
% in vessels.
%
% INPUTS:
%
% 1- folderPath: Folder containing the speckle data (called speckle.dat)
%
% 2- sType: how the stdev should be computed. Two options:
%       - Spatial: stdev will be computed on a 5x5 area in the XY plane
%       - Temporal: stdev will be computed on a 5x1 vector in the time dimension 
%
% 3- channel (optional): Channel to analyse, for example 'green', 'red',
% etc. (speckle by default)
%
% 4- bSaveMap: Save a map of averaged stddev over the acquisition (.tiff
% file)
%
% 5- bLogScale: boolean flag to put data on a -log10 scale
%           - true: ouput data is equal to -log10(data)
%           - false: data = data;
%
% OUTPUT:
%
% 1- DatOut: StDev variation averaged over time.

if(nargin < 3)
    channel = 'speckle';
%     bSaveStack = 1;
    bSaveMap = 1;
    bLogScale = 1;
end

if( ~strcmp(folderPath(end), filesep) )
    folderPath = strcat(folderPath, filesep);
end

channel = lower(channel);
if(~exist([folderPath channel '.dat'],'file') )
    disp([channel '.dat file is missing. Did you run ImagesClassificiation?']);
    return;
end

disp(['Opening ' channel '.dat']);

try
    Infos = matfile([folderPath channel '.mat']);
    fid = fopen([folderPath channel '.dat']);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat, Infos.datSize(1,1), Infos.datSize(1,2),[]);
    dat = dat./mean(dat,3);
catch 
    disp(['Failed to open ' channel ' files'])
    return
end

disp('Mapping Computation');
switch lower(sType)
    case 'spatial'
        Kernel = zeros(5,5,1,'single');
        Kernel(:,:,1) = single(fspecial('disk',2)>0);        
    case 'temporal'
        Kernel = ones(1,1,5,'single');        
end

DatOut = stdfilt(dat,Kernel);
DatOut = single(DatOut);

%Remove outliers
pOutlier = prctile(DatOut(:), 99);
DatOut(DatOut>pOutlier) = pOutlier;
% Get the average speckle contrast map:
DatOut = mean(DatOut./mean(dat,3),3);

if( bLogScale )
    DatOut = -log10(DatOut);
end

%Generate output
% copyfile([folderPath channel '.mat'], [folderPath flow '.mat'])
% if( bSaveStack )
%     disp('Saving');
%     mFileOut = matfile([folderPath 'flow.mat'], 'Writable', true);
%     mFileOut.FirstDim = Infos.FirstDim;
%     mFileOut.Freq = Infos.Freq;
%     mFileOut.Stim = Infos.Stim;
%     mFileOut.datLength = Infos.datLength;
%     mFileOut.datSize = Infos.datSize;
%     mFileOut.datFile = 'flow.dat';
% 
%     fid = fopen([folderPath 'flow.dat'],'w'); 
%     fwrite(fid, single(DatOut), 'single');
%     fclose(fid);
% end
if( bSaveMap )
%     Map = mean(DatOut,3);
    obj = Tiff(fullfile(folderPath, 'std_speckle.tiff'), 'w');    
    setTag(obj, 'ImageWidth', size(DatOut,2));
    setTag(obj, 'ImageLength', size(DatOut,1));
    setTag(obj, 'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(obj, 'SampleFormat',Tiff.SampleFormat.IEEEFP);
    setTag(obj, 'BitsPerSample', 32);
    setTag(obj, 'SamplesPerPixel', 1);
    setTag(obj, 'Compression',Tiff.Compression.None);
    setTag(obj, 'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    write(obj, DatOut);
end
disp('Done');
end
