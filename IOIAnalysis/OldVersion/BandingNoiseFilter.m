function out = BandingNoiseFilter(FolderName, OStream)

%%%%%%%%%%% 
% Opening %
%%%%%%%%%%%
%List of files to be included in the analysis
FileList = dir([FolderName filesep 'Data_*.mat']);
if( isempty(FileList) )
    %No compatible files were found
    disp(['No data files found in ' FolderName ' Folder.']);
    disp('Analysis will not run');
    return;
end

AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
    'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);

%Green channel detected
if( ~isempty(strfind([FileList.name],'green')) )
    if( isempty(OStream) )
        fprintf('Green channel banding noise filtering:\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Green channel banding noise filtering:',...
            OStream.String);
        drawnow;
    end
    
    Dat_ptr = matfile([FolderName filesep 'Data_green.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    nfram = Dat_ptr.datLength;
    BN_Filter(Dat_ptr.datFile, nrows, ncols, nfram);
    
    if( isempty(OStream) )
         fprintf('Done.\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Done.',...
            OStream.String);
        drawnow;
    end    
end
%Yellow channel detected
if( ~isempty(strfind([FileList.name],'yellow')) )
     if( isempty(OStream) )
        fprintf('Yellow channel banding noise filtering:\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Yellow channel banding noise filtering:',...
            OStream.String);
        drawnow;
     end
     
    Dat_ptr = matfile([FolderName filesep 'Data_yellow.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    nfram = Dat_ptr.datLength;
    BN_Filter(Dat_ptr.datFile, nrows, ncols, nfram);
    
    if( isempty(OStream) )
         fprintf('Done.\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Done.',...
            OStream.String);
        drawnow;
    end    
end
%Red channel detected
if( ~isempty(strfind([FileList.name],'red')) )
     if( isempty(OStream) )
        fprintf('Red channel banding noise filtering:\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Red channel banding noise filtering:',...
            OStream.String);
        drawnow;
    end
    Dat_ptr = matfile([FolderName filesep 'Data_red.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    nfram = Dat_ptr.datLength;
    
    BN_Filter(Dat_ptr.datFile, nrows, ncols, nfram);
    if( isempty(OStream) )
         fprintf('Done.\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Done.',...
            OStream.String);
        drawnow;
    end    
end
%Speckle channel detected
if( ~isempty(strfind([FileList.name],'speckle')) )
    if( isempty(OStream) )
        fprintf('Speckle channel banding noise filtering:\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Speckle channel banding noise filtering:',...
            OStream.String);
        drawnow;
    end
    Dat_ptr = matfile([FolderName filesep 'Data_speckle.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    nfram = Dat_ptr.datLength;
    
    BN_Filter(Dat_ptr.datFile, nrows, ncols, nfram);
    if( isempty(OStream) )
         fprintf('Done.\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Done.',...
            OStream.String);
        drawnow;
    end    
end
%Fluo channel detected
if( ~isempty(strfind([FileList.name],'Fluo')) )
    if( isempty(OStream) )
        fprintf('Fluo channel banding noise filtering:\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Fluo channel banding noise filtering:',...
            OStream.String);
        drawnow;
    end
    Dat_ptr = matfile([FolderName filesep 'Data_Fluo.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    nfram = Dat_ptr.datLength;
    
    BN_Filter(Dat_ptr.datFile, nrows, ncols, nfram);
    if( isempty(OStream) )
         fprintf('Done.\n')
    else
        OStream.String = sprintf('%s\r%s',...
            'Done.',...
            OStream.String);
        drawnow;
    end    
end

function BN_Filter(DatFile, nrows, ncols, nFrames)

SizeTot = 4*double(nrows)*double(ncols)*double(nFrames)/1.024e9; %in GByte;
if( ispc )
    SizeLim = memory;
    SizeLim = SizeLim.MaxPossibleArrayBytes/1.024e9;
    SizeLim = SizeLim*0.125;
else
    SizeLim = 4;
end
NbInc = round(SizeTot/SizeLim);
Lims = round(linspace(1, nFrames, NbInc+1));

minD = min(ncols,nrows);
Tags = round(linspace(0, nFrames, 20));
indT = 1;
indFt = 0;
if( ~isempty(OStream) )
    StaticStr = OStream.String;
else
    StaticStr = [];
end
if( NbInc > 2 )
    for indI = 2:NbInc
        Frame_ptr = memmapfile(DatFile,'Writable', true,...
            'Offset', 4*(Lims(indI-1)-1)*double(ncols)*double(nrows), ...
            'Format', 'single', 'repeat', (Lims(indI) - Lims(indI-1))*double(ncols)*double(nrows));
        Frame = reshape(Frame_ptr.Data, nrows, ncols, []);
        mFrame = mean(Frame,3);
        
        for indF = 1:(Lims(indI) - Lims(indI-1))
            indFt = indFt + 1;
            if( indF == 1 )
                Tmp = squeeze(Frame(:,:,indF)-median(Frame(:,:,2:10),3));
            else
                Tmp = squeeze(Frame(:,:,indF) - Frame(:,:,indF-1));
            end
            Frame(:,:,indF) = Frame(:,:,indF) - ...
                (Tmp - ...
                single(xRemoveStripesVertical(Tmp, nextpow2(minD)-4, 'db4', 2)));
            if( indFt >= Tags(indT) )
                P = round((100*Tags(indT))/nFrames);
                if( isempty(OStream) )
                    fprintf('%d%% .. ', P);
                    if( indT == 10 )
                        fprintf('\n');
                    end
                    
                else
                    OStream.String = sprintf('%s\r%s',...
                        ['Completion: ' int2str(P) '%'],...
                        StaticStr);
                    drawnow;
                end
                indT = indT + 1;
            end
        end
    end
    OStream.String = StaticStr;
    Frame_ptr.Data = Frame;
else
     Frame_ptr = memmapfile(DatFile,'Writable', true, 'Offset', 0, ...
            'Format', 'single');
     Frame = reshape(Frame_ptr.Data, nrows, ncols, []);
     mFrame = mean(Frame,3);
        
     for indF = 1:size(Frame,3)
        Tmp = Frame(:,:,indF)-mFrame;
        Frame(:,:,indF) = Frame(:,:,indF) - ...
             (Tmp - ...
             single(xRemoveStripesVertical(Tmp, 4 + nextpow2(minD), 'db4', 2)));
            
         if( indF >= Tags(indT) )
             P = round((100*Tags(indT))/nFrames);
             if( isempty(OStream) )
                 fprintf('%d%% .. ', P);
                 if( indT == 10 )
                     fprintf('\n');
                 end
                 
             else
                 OStream.String = sprintf('%s\r%s',...
                     ['Completion: ' int2str(P) '%'],...
                     StaticStr);
                 drawnow;
             end
             indT = indT + 1;
         end
     end
     if( ~isempty(OStream) )
        OStream.String = StaticStr;
     end
     Frame_ptr.Data = Frame;
end
end
end