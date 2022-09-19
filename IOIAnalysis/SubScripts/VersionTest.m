function V = VersionTest(FolderName)


tOld = dir( [FolderName filesep 'IOI_scan.seq']);
tNew = dir( [FolderName filesep 'img_*.bin']);
if( numel(tOld) > 0 )
%    disp('System version detected: IOI 1.0');
    V = '1.0';
    %OpenIOI_OldSyst(FolderName, Binning);
elseif( numel(tNew) > 0 )
    header = memmapfile([FolderName filesep 'img_00000.bin'], ...
    'Offset', 0, 'Format', 'int32', 'repeat', 1);
    %FLAG for header version here!
    headerVer = header.Data;
    if( headerVer == 3 )
        AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
                'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);
            
        if( strcmp(AcqInfoStream(5,1).Var1,'AISampleRate') )
%        disp('System version detected: IOI 2.2');
            V = '2.2';
            AcqInfoStream2 = ReadInfoFile(FolderName);
            if( isfield(AcqInfoStream2, 'Camera_Model') )
               if (strcmp(AcqInfoStream2.Camera_Model, 'D1024') || strcmp(AcqInfoStream2.Camera_Model, 'D1312')) 
                   V = '2.3';
               end
            end
            
        else
            V = '2.3';
        end
    elseif( headerVer == 2 )
%        disp('System version detected: IOI 2.1');
        V = '2.1';
    elseif( headerVer == 1 )
        V = '2.0';
%        disp('System version detected: IOI 2.0');
    end
else
    V = '0.0';
%      disp('No Valid data set found.');
end

end