function out = FluoCorrGen(FolderName, Wlength)

if( ~strcmp(FolderName(end), filesep) )
    FolderName = strcat(FolderName, filesep);
end

disp('Opening HB Files...');
% Open HB Files
%%%%%%%%%%%%%%%
if( ~exist([FolderName 'HbO.dat'], 'file') )
    disp('No HbO file found. Run Ana_IOI_FullFrame.');
    return;
end

HBInfo = matfile([FolderName filesep 'Data_Hbs.mat']);
fid = fopen([FolderName 'HbO.dat']);
dat_HbO = fread(fid,'single');
fclose(fid);
dat_HbO = single(reshape(dat_HbO, HBInfo.datSize(1,1), HBInfo.datSize(1,2),[]));

if( ~exist([FolderName 'HbR.dat'], 'file') )
    disp('No HbR file found. Run Ana_IOI_FullFrame.');
    return;
end

fid = fopen([FolderName 'HbR.dat']);
dat_HbR = fread(fid,'single');
fclose(fid);
dat_HbR = single(reshape(dat_HbR, HBInfo.datSize(1,1), HBInfo.datSize(1,2),[]));

%%TODO
dat_HbR = dat_HbR + 40;
dat_HbR(dat_HbR > 120) = 120;
dat_HbR(dat_HbR < 0) = 0;
dat_HbO = dat_HbO + 60;
dat_HbO(dat_HbO > 120) = 120;
dat_HbO(dat_HbO < 0) = 0;
dat_HbR = round(dat_HbR*4) + 1;
dat_HbO = round(dat_HbO*4) + 1;
%%Fin

disp('Loading HB lookup table...');
disp(['Computing Correction for ' Wlength 'nm.']);
    
eval(['File = ''FactCorrHB_' Wlength '.mat'';']);
            
if( ~exist(File,'file') )
    HBLookupTableGen(Wlength);
end

FC = load(File);

Corr = FC.Corr(sub2ind(size(FC.Corr),dat_HbR(:),dat_HbO(:)));
Corr = reshape(Corr, size(dat_HbO));

disp('Saving Correction...');
fid = fopen([FolderName 'hCorr' Wlength '.dat'],'w');
fwrite(fid, Corr, 'single');
fclose(fid);

disp('Done.');

end