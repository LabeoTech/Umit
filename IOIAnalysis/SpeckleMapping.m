function DatOut = SpeckleMapping(dataSource, sType, channel, bSaveMap, bLogScale)
%%%%%%%%%%%%%%%%%%%% Speckle Mapping function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the standard deviation (spatialy or temporaly) of speckle
% acquisition. This measure is proportional to the strength of blood flow
% in vessels.
%
% INPUTS:
%
% 1- dataSource (char | 3D num. array): Folder containing the speckle data
% (called speckle.dat) OR the speckle data itself with dimensions (Y,X,T);
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

if ischar(dataSource)
    % If the data source is a path to a folder, open the speckle file:
    if( ~strcmp(dataSource(end), filesep) )
        dataSource = strcat(dataSource, filesep);
    end
    
    channel = lower(channel);
    if(~exist([dataSource channel '.dat'],'file') )
        disp([channel '.dat file is missing. Did you run ImagesClassificiation?']);
        return;
    end
    
    disp(['Opening ' channel '.dat']);
    
    try
        Infos = matfile([dataSource channel '.mat']);
        fid = fopen([dataSource channel '.dat']);
        DatOut = fread(fid, inf, '*single');
        DatOut = reshape(DatOut, Infos.datSize(1,1), Infos.datSize(1,2),[]);
        %
    catch
        disp(['Failed to open ' channel ' files'])
        return
    end
else
    %     % If the data source is already the speckle data:
    DatOut = dataSource + 0; % This operation is just to force Matlab to load the data in RAM so the availabe memory estimation works properly.    
    clear dataSource
end
% Normalize data by temporal mean:
DatOut = DatOut./mean(DatOut,3,'omitnan');

disp('Mapping Computation...');
switch lower(sType)
    case 'spatial'
        Kernel = zeros(5,5,1,'double');
        Kernel(:,:,1) = double(fspecial('disk',2)>0);
    case 'temporal' 
        datSz = size(DatOut);
        Kernel = ones(1,1,5,'double');
end
% Check if the data needs to be chunked to be processed by "stdfilt":
nChunks = calculateMaxChunkSize(numel(DatOut),10);
if nChunks > 1
    % If there is not enough RAM to run the std filter in the whole data,
    % split it in smaller sizes
    
    if strcmpi(sType, 'temporal')
        indxChk = round(linspace(0,datSz(1),nChunks));
        for ii = 1:length(indxChk)-1            
%             DatOut(indxChk(ii)+1:indxChk(ii+1),:,:)= single(stdfilt(DatOut(indxChk(ii)+1:indxChk(ii+1),:,:),Kernel)./...
%                 imfilter(DatOut(indxChk(ii)+1:indxChk(ii+1),:,:), Kernel/(5*5),'conv'));
            DatOut(indxChk(ii)+1:indxChk(ii+1),:,:)= single(stdfilt(DatOut(indxChk(ii)+1:indxChk(ii+1),:,:),Kernel));
            fprintf('[%1.0f%%]\n', 100*ii/(length(indxChk)-1));
        end
    else
        indxChk = round(linspace(0,size(DatOut,3),nChunks));
        for ii = 1:length(indxChk)-1
%             DatOut(:,:,indxChk(ii)+1:indxChk(ii+1)) = single(stdfilt(DatOut(:,:,indxChk(ii)+1:indxChk(ii+1)), Kernel)./...
%                 imfilter(DatOut(:,:,indxChk(ii)+1:indxChk(ii+1)),Kernel/(5*5),'conv'));
            DatOut(:,:,indxChk(ii)+1:indxChk(ii+1)) = single(stdfilt(DatOut(:,:,indxChk(ii)+1:indxChk(ii+1)), Kernel));
            fprintf('[%1.0f%%]\n', 100*ii/(length(indxChk)-1));
        end
    end
    fprintf('[Completed]\n')
else    
%     DatOut = single(stdfilt(DatOut,Kernel)./imfilter(DatOut, Kernel./(5*5),'conv'));
      DatOut = single(stdfilt(DatOut,Kernel));

end

% Remove outliers
% remOutlier; % This was taking too much RAM for large data.

% Get the average speckle contrast map:
DatOut = mean(DatOut,3);

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
if( bSaveMap && exist('dataSource','var') )
    %     Map = mean(DatOut,3);
    obj = Tiff(fullfile(dataSource, 'std_speckle.tiff'), 'w');
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
% Local function
% ------------------------------------------------------------------------
%     function remOutlier
%         % Custom function to remove outliers from data by imputing the value of the
%         % 99th percentile to higher values.
%         dataSz = size(data);
%         data = data(:);
%         nChunk = calculateMaxChunkSize(data,4);
%         val99th = zeros(1,nChunk,'single');
%         indx = round(linspace(1,length(data),nChunk+1));
%         disp('Removing outliers...')
%         for jj = 1:length(indx)-1
%             snippet = sort(data(indx(jj):indx(jj+1)));
%             idxNaN = isnan(snippet);
%             snippet(idxNaN) = [];
%             idx99 = .99*(length(snippet));
%             if round(idx99) == idx99
%                 val99th(jj) = snippet(idx99+1);
%             else
%                 val99th(jj) = (snippet(floor(idx99+1))+snippet(ceil(idx99+1)))/2;
%             end
%             fprintf('[%1.0f%%]\n', 100*jj/(length(indx)-1));
%         end
%         fprint('[Completed]\n')
%         clear snippet idxNaN indx
%         % Estimate the outliers as the average of the 99th perc. from
%         % each chunk:
%         valOutLier = mean(val99th);
%         data(data>valOutLier) = valOutLier;
%         % rebuild data:
%         data = reshape(data,dataSz);
%     end
end


