function Coreg2Cams(DataFolder, bSave)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end
load([DataFolder 'AcqInfos.mat']);

disp('Opening');
NbChan = sum(contains(fieldnames(AcqInfoStream), 'Illumination'));
indC1 = 1;
indC2 = 1;
for ind = 1:NbChan
    eval(['Color = AcqInfoStream.Illumination' int2str(ind) '.Color;']);
    if( contains(Color, '#') )
        idx = strfind(Color,' nm');
        Color = ['fluo_' Color(idx - 3 + (0:2))];
    end
    Infos = matfile([DataFolder lower(Color) '.mat']);
    fid = fopen([DataFolder lower(Color) '.dat']);
    tmp = fread(fid, Infos.datSize(1,1)*Infos.datSize(1,2), '*single');
    tmp = reshape(tmp,Infos.datSize(1,1),Infos.datSize(1,2));
    eval(['cIdx = AcqInfoStream.Illumination' int2str(ind) '.CamIdx;']);
    if( cIdx == 1 )
        Cam1(:,:,indC1) = tmp;
        indC1 = indC1 + 1;
    else
        Cam2(:,:,indC2) = tmp;
        indC2 = indC2 + 1;
    end
end
disp('Coregistering');
Cam1 = sum(Cam1,3);
Cam2 = sum(Cam2,3);
Cam2 = (Cam2 - min(Cam2(:)));
Cam1 = (Cam1 - min(Cam1(:)));

T1 = affine2d();
T1.T = [1 0 0; 0 1 0; 0 0 1];
[opt1, metric] = imregconfig('monomodal');
T2 = imregtform(stdfilt(Cam2), stdfilt(Cam1), 'rigid', opt1, metric);

save([DataFolder 'CoregTransf.mat'], 'T1', 'T2');
if( bSave )
    disp('Saving');
    for ind = 1:NbChan
        eval(['Color = AcqInfoStream.Illumination' int2str(ind) '.Color;']);
        if( contains(Color, '#') )
            idx = strfind(Color,' nm');
            Color = ['fluo_' Color(idx - 3 + (0:2))];
        end
        Infos = matfile([DataFolder lower(Color) '.mat']);
        
        eval(['cIdx = AcqInfoStream.Illumination' int2str(ind) '.CamIdx;']);
        if( cIdx == 2 )
            fid = fopen([DataFolder lower(Color) '.dat']);
            tmp = fread(fid, inf, '*single');
            tmp = reshape(tmp, Infos.datSize(1,1), Infos.datSize(1,2), []);
            fclose(fid);
            for indF = 1:size(tmp,3)
               tmp(:,:,indF) = imwarp(squeeze(tmp(:,:,indF)), T2,...
                   'OutputView', imref2d(Infos.datSize)); 
            end
            fid = fopen([DataFolder lower(Color) '.dat'],'w');
            fwrite(fid, tmp,'single');
            fclose(fid);
        end
    end
end
end