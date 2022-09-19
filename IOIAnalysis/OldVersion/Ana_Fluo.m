function out = Ana_Fluo(FolderName, b_HbCorr, b_HilbertF)
%%%%%%%%%%%%%%%%%%%%%%
% Validation & opening
%%%%%%%%%%%%%%%%%%%%%%

if( ~strcmp(FolderName(end), filesep) )
    FolderName = strcat(FolderName, filesep);
end

FileList = dir([FolderName '*Fluo*.mat']);
if( isempty(FileList) )
    disp(['No Fluorescence data files found in ' FolderName ' Folder.']);
    disp('Fluorescence Analysis will not run');
    return;
end

for indF = 1:size(FileList,1)
    FInfo =  matfile([FileList(indF).folder filesep FileList(indF).name], 'Writable', true);
    fprintf(['Opening file: ' cell2mat(FInfo.Wavelength(1,1)) 'nm.\n']);
    fid = fopen([FolderName 'fChan_' cell2mat(FInfo.Wavelength(1,1)) '.dat']);
    Fluo = fread(fid, inf, 'single=>single');
    fclose(fid);

    if( exist([FolderName 'Data_Hbs.mat'], 'file') && b_HbCorr )        
        fprintf('Applying HB Correction.\n');
        
        FluoCorrGen(FolderName, cell2mat(FInfo.Wavelength(1,1)));
    
        HBInfo = matfile([FolderName 'Data_Hbs.mat']);
        fid = fopen([FolderName 'hCorr' cell2mat(FInfo.Wavelength(1,1)) '.dat']);
        Corr = fread(fid, inf, 'single');
        Corr = single(Corr);
        fclose(fid);
        
        Fluo = Fluo(1:length(Corr));
        Fluo = Fluo.*Corr;
    end
    
    Fluo = reshape(Fluo, FInfo.datSize(1,1), FInfo.datSize(1,2),[]);
    FInfo.datLength = size(Fluo,3);
    
    fprintf('Filtering.\n');
    if( ~b_HilbertF )
        dims = size(Fluo);
        
        %Filter parameters
        f = fdesign.lowpass('N,F3dB', 4, FInfo.Freq/2, FInfo.Freq);
        hpass = design(f,'butter');
        f = fdesign.lowpass('N,F3dB', 4, 1/10, FInfo.Freq);
        lpass = design(f,'butter');
        
        %Memory management
        tSiz = (prod(dims)*8)/1e9;
        tSiz = (16 + 8*tSiz);
        nbStep = ceil(tSiz / 32);
        
        lims = round(linspace(1, dims(1)*dims(2) + 1, nbStep+1));
        Fluo = reshape(Fluo, [], dims(3));
        for indL = 2:length(lims)
            imin = lims(indL - 1);
            imax = lims(indL) - 1;
            dFh = single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(Fluo(imin:imax,:)')))';
            dFl = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(Fluo(imin:imax,:)')))';
            Fluo(imin:imax,:) = (dFh - dFl)./dFl;
            clear dFh dFl;
        end
        
        Fluo  = reshape(Fluo, dims(1), dims(2), dims(3));
    else
        dims = size(Fluo);
        [~, ylower] = envelope(reshape(Fluo,[],dims(3))',4*FInfo.Freq,'peak');
        ylower = reshape(ylower', dims);
        Fluo = (Fluo - ylower)./ylower;
        for ind = 1:dims(3)
            Fluo(:,:,ind) = medfilt2(squeeze(Fluo(:,:,ind)),[5 5],'symmetric');
        end
    end
    
    fprintf('Saving.\n');
    fid = fopen([FolderName 'fChan_' cell2mat(FInfo.Wavelength(1,1)) '.dat'],'w');
    fwrite(fid,Fluo,'single');
    fclose(fid);

end

