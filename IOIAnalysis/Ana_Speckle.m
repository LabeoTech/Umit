function varargout = Ana_Speckle(Folder, bNormalize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Infos:
%
% This function is used to calculate the blood flow from Laser speckle
% data.
%
% Inputs:
% 1. Folder: Folder contaning the dataset to work with.
% Ouput:
% - If an output is set, the 3D matrix containing the blood flow data it's meta data
%    will be given back through this output.
% - If no output is specified, the "Flow.dat" and "Flow.mat" files will  be saved
%    in Folder.
% 
% Exemples:
%
% 1- Ana_Speckle(Folder); 
% The blood flow will be calculated from the "speckle.dat" in  "Folder". 
% The "Flow.dat" and "Flow.mat" files will be created in same folder.
% 2- data = Ana_Speckle(Folder);
% The 3D matrix containing the blood flow data is given as output.
% 3- [data, metaData] = Ana_Speckle(Folder);
% A structure containing the data's meta data is also given as output.

disp('Running Ana speckle...')
FileList = dir([Folder filesep 'speckle.mat']);
if( isempty(FileList) )
    disp(['No speckle data files found in ' Folder ' Folder.']);
    disp('Speckle Analysis will not run');
    return;
else
    Iptr = matfile([Folder filesep 'speckle.mat']);
    nx = Iptr.datSize(1,2);
    ny = Iptr.datSize(1,1);
    nt = Iptr.datLength;
    tFreq = Iptr.Freq;
    speckle_int_time = Iptr.tExposure/1000.;    
    Dptr = memmapfile(fullfile(Folder,Iptr.datFile), 'Format', 'single');
end

% Parameters
% speckle_window_size = 5;

fprintf('Opening files.\n');

fid = fopen([Folder filesep 'speckle.dat']);
dat = fread(fid,inf,'single');
fclose(fid);
dat = reshape(dat, ny, nx,[]);
MeanMap = mean(dat,3);
% Need to convert to contrast
fprintf('Flow Conversion:\n');
speckle_window = fspecial('disk',2)>0;
OPTIONS.GPU = 0;
OPTIONS.Power2Flag = 0;
OPTIONS.Brep = 0;
dat = zeros(nt - 1, ny, nx, 'single');
prcflg = linspace(1, nt-1, 11); indP = 2;
for i3 = 1:nt-1
    if( i3 >= prcflg(indP) )
         fprintf('%d%%...', 10*(indP-1));
         indP = indP+1;
     end
    tmp_laser = Dptr.Data(((i3-1)*nx*ny + 1):(i3*nx*ny));
    tmp_laser = reshape(tmp_laser, ny, []);  
    tmp_laser = tmp_laser./MeanMap;
    std_laser = imgaussfilt(stdfilt(tmp_laser,speckle_window),1);
    mean_laser = imgaussfilt(convnfft(tmp_laser,speckle_window,'same',1:2,OPTIONS)/sum(speckle_window(:)),1);
    contrast=std_laser./mean_laser;
    dat(i3, :, :) = single(private_flow_from_contrast(contrast,speckle_int_time));
end
clear tmp_laser std_laser contrast mean_laser;

fprintf('\nFiltering:\n');
fW = ceil(0.5*tFreq);
dat = medfilt1(dat, fW, [], 1, 'truncate');
fprintf('100%%.');
fprintf('\nSaving...\n');
% Save/Output data and meta data:
dat = permute(dat, [2 3 1]);

%Normalization:
if( bNormalize )
    dat = dat./mean(dat,3);
end

% Create meta data structure:
metaData = struct();
metaData.datName = 'data';
metaData.Stim = Iptr.Stim;
metaData.datLength = Iptr.datLength-1;
metaData.datSize = Iptr.datSize;
metaData.Freq = Iptr.Freq;
metaData.Datatype = class(dat);
metaData.dim_names = Iptr.dim_names;
metaData.datFile = [Folder filesep 'Flow.dat'];

out_vars = {'dat', 'metaData'};
if nargout > 0
    for i = 1:nargout
        eval(['varargout{' num2str(i) '} = ' out_vars{i}, ';']);
    end    
else
    disp('Saving data to file : "Flow.dat"...')
    fFlow = fopen([Folder filesep 'Flow.dat'], 'w');
    fwrite(fFlow, dat, 'single');
    fclose(fFlow);
    save(fullfile(Folder, 'Flow.mat'),  '-struct', 'metaData') 
end

fprintf('Done!\n');
end

function speed = private_flow_from_contrast(contrast,T)
% trouve la vitesse � partir des images de contraste
%[nx ny] = size(contrast);

% Correct for points that cannot be...
contrast(isnan(contrast)|contrast<0)=0;

% Not sure about this, verify
contrast2=contrast(3:end-2,3:end-2);
mmean=mean(contrast2(1:end));
sstd=std(contrast2(1:end));

% Build non-linear curve between contrast and correlation time (tau)
tau=(logspace(-15,0,60).^.5); % Correlation time
K  = ((tau/(2*T)).*(1-exp(-2*T*ones(size(tau))./tau))).^(1/2);
% Find values for which the mean contrast is in the middle
[~, index1]=find(K>(mmean-3*sstd),1);
[~, index2]=find(K>(mmean+3*sstd),1);
if isempty(index1), index1=1; end
if isempty(index2)||index2==index1, index2=60; end

% For these values, build a log-linear vector on which contrast is computed
Tau2=(logspace(log10(tau(index1)),log10(tau(index2)),40));
K  = ((Tau2/(2*T)).*(1-exp(-2*T*ones(size(Tau2))./Tau2))).^(1/2);

% Add limit points for interpolation
Tau2=[Tau2(1) Tau2 Tau2(end)];
K= [ 0 K 1e30];
% Interpolate contrast image on these values to obtain correlation time
Tau3=interp1(K,Tau2,contrast); %This is SLOW
% Get speed from correlation time
speed=1./Tau3;
end