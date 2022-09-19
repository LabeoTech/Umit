function MapAveraging(FolderPath)

if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

Flist = dir(FolderPath);
Flist = Flist(3:end);
idx = arrayfun(@(x) Flist(x).isdir, 1:size(Flist,1));
Flist = Flist(idx);

idx = zeros(size(Flist));
for ind = 1:size(idx,1)
   SubList = dir([FolderPath Flist(ind).name filesep]);
   idx(ind) = contains([SubList(:).name],'DataMapping.mat');
end
Flist = Flist(idx>0);
clear idx SubList ind

disp('Experiments to be averaged:')
for ind = 1:size(Flist,1)
    disp(Flist(ind).name);
end

disp('Loading...')
Dat = {};
for ind = 1:size(Flist,1)
   load([FolderPath Flist(ind).name filesep 'DataMapping.mat']); 
   Dat{ind} = Data;
end
clear Data ind

disp('CoRegistration...');
Yax = linspace(-5,5,50);
Xax = linspace(-5,5,50);
[X, Y] = meshgrid(Yax, Xax);
Max_123 = zeros(size(Dat,2),length(Yax),length(Xax));
Max_567 = zeros(size(Dat,2),length(Yax),length(Xax));
Tmax_123 = zeros(size(Dat,2),length(Yax),length(Xax));
Tmax_567 = zeros(size(Dat,2),length(Yax),length(Xax));
Ton_123 = zeros(size(Dat,2),length(Yax),length(Xax));
Ton_567 = zeros(size(Dat,2),length(Yax),length(Xax));
for ind = 1:size(Dat,2)
    F = scatteredInterpolant(repmat(Dat{ind}.coronal_axis, size(Dat{ind}.sagittal_axis,1),1),...
        repelem(Dat{ind}.sagittal_axis, size(Dat{ind}.coronal_axis,1),1),...
        reshape(flipud(rot90(Dat{ind}.valueMax_123)),[],1), 'linear', 'none');
    Max_123(ind,:,:) = F(X,Y);
    
    F = scatteredInterpolant(repmat(Dat{ind}.coronal_axis, size(Dat{ind}.sagittal_axis,1),1),...
        repelem(Dat{ind}.sagittal_axis, size(Dat{ind}.coronal_axis,1),1),...
        reshape(flipud(rot90(Dat{ind}.valueMax_567)),[],1), 'linear', 'none');
    Max_567(ind,:,:) = F(X,Y);
    
    F = scatteredInterpolant(repmat(Dat{ind}.coronal_axis, size(Dat{ind}.sagittal_axis,1),1),...
        repelem(Dat{ind}.sagittal_axis, size(Dat{ind}.coronal_axis,1),1),...
        reshape(flipud(rot90(Dat{ind}.timeMax_123)),[],1), 'linear', 'none');
    Tmax_123(ind,:,:) = F(X,Y);
    
    F = scatteredInterpolant(repmat(Dat{ind}.coronal_axis, size(Dat{ind}.sagittal_axis,1),1),...
        repelem(Dat{ind}.sagittal_axis, size(Dat{ind}.coronal_axis,1),1),...
        reshape(flipud(rot90(Dat{ind}.timeMax_567)),[],1), 'linear', 'none');
    Tmax_567(ind,:,:) = F(X,Y);
   
    F = scatteredInterpolant(repmat(Dat{ind}.coronal_axis, size(Dat{ind}.sagittal_axis,1),1),...
        repelem(Dat{ind}.sagittal_axis, size(Dat{ind}.coronal_axis,1),1),...
        reshape(flipud(rot90(Dat{ind}.timeOnset_123)),[],1), 'linear', 'none');
    Ton_123(ind,:,:) = F(X,Y);
    
    F = scatteredInterpolant(repmat(Dat{ind}.coronal_axis, size(Dat{ind}.sagittal_axis,1),1),...
        repelem(Dat{ind}.sagittal_axis, size(Dat{ind}.coronal_axis,1),1),...
        reshape(flipud(rot90(Dat{ind}.timeOnset_567)),[],1), 'linear', 'none');
    Ton_567(ind,:,:) = F(X,Y);
end
clear F;
Max_123(isnan(Max_123)) = 0;
Max_567(isnan(Max_567)) = 0;
Tmax_123(isnan(Tmax_123)) = 0;
Tmax_567(isnan(Tmax_567)) = 0;
Ton_123(isnan(Ton_123)) = 0;
Ton_567(isnan(Ton_567)) = 0;

disp('Generate figures...')
hFig = figure('Position', [2661 331 600 600]);
imagesc(Yax, Xax, squeeze(mean(Max_123,1)));
axis image;
title('Max, ch. 123');
saveas(hFig, [FolderPath 'MaxAmpChan123.png']);

imagesc(Yax, Xax, squeeze(mean(Max_567,1)));
axis image;
title('Max, ch. 567');
saveas(hFig, [FolderPath 'MaxAmpChan567.png']);

imagesc(Yax, Xax, squeeze(mean(Tmax_123,1)));
axis image;
title('Time to Max, ch. 123');
saveas(hFig, [FolderPath 'TmaxChan123.png']);

imagesc(Yax, Xax, squeeze(mean(Tmax_567,1)));
axis image;
title('Time to Max, ch. 567');
saveas(hFig, [FolderPath 'TmaxChan567.png']);

Tmp = Ton_123;
Tmp(Tmp <= 0) = nan;
imagesc(Yax, Xax, squeeze(mean(Tmp,1,'omitnan')));
axis image;
title('Onset Time, ch. 123');
saveas(hFig, [FolderPath 'TonChan123.png']);

Tmp = Ton_567;
Tmp(Tmp <= 0) = nan;
imagesc(Yax, Xax, squeeze(mean(Tmp,1,'omitnan')));
axis image;
title('Onset Time, ch. 567');
saveas(hFig, [FolderPath 'TonChan567.png']);

disp('Saving...')
save([FolderPath 'AveragedData.mat'],'Max_123','Max_567','Tmax_123','Tmax_567','Ton_123','Ton_567');

disp('Done.');
end