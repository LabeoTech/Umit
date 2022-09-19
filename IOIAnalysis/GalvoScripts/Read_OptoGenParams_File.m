function out = Read_OptoGenParams_File(FolderPath)

if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

fid = fopen([FolderPath 'OptoGenParams.txt']);
tline = fgetl(fid);

if( contains(tline, 'Single') )
    tline = fgetl(fid);
    while( ~feof(fid) )
       tline = fgetl(fid);
       tag = tline(1:(strfind(tline,':')-1));
       switch(tag)
           case 'PreCond Power'
               tline = regexprep(tline, '\t', ' ');
               indices = strfind(tline((strfind(tline,':')+1):end), ' ') + (strfind(tline,':') + 1);
               Pwrs = str2num(tline(indices(1):end));
               out.PreCondPwr = Pwrs;             
           case 'PreCond Ton'
               out.PreCondTON = str2num(tline((strfind(tline,':')+1):end));
           case 'PreCond Xpos (pix)' 
               out.PreCondXpix = str2num(tline((strfind(tline,':')+1):end));
           case 'PreCond Ypos (pix)'
                out.PreCondYpix = str2num(tline((strfind(tline,':')+1):end));
           case 'PreCond Xpos (mm)'
                out.PreCondXmm = str2num(tline((strfind(tline,':')+1):end));
           case 'PreCond Ypos (mm)'
               out.PreCondYmm = str2num(tline((strfind(tline,':')+1):end));
           case 'InterStim Delay'
               out.IStimDelay = str2num(tline((strfind(tline,':')+1):end));
           case 'OptoGen Power'
               out.OptoGenPwr = str2num(tline((strfind(tline,':')+1):end));
           case 'OptoGen Ton'
               out.OptoGenTON = str2num(tline((strfind(tline,':')+1):end));
           case 'OptoGen Xpos (pix)'
               out.OptoGenXpix = str2num(tline((strfind(tline,':')+1):end));
           case 'OptoGen Ypos (pix)'
               out.OptoGenYpix = str2num(tline((strfind(tline,':')+1):end));
           case 'OptoGen Xpos (mm)'
               out.OptoGenXmm = str2num(tline((strfind(tline,':')+1):end));
           case 'OptoGen Ypos (mm)'
               out.OptoGenYmm = str2num(tline((strfind(tline,':')+1):end));
           case 'NbReps'
               out.NbReps = str2num(tline((strfind(tline,':')+1):end));
           case 'Trig Offset (msec)'
               out.TOffset = str2num(tline((strfind(tline,':')+1):end));
           case 'NbEvents'
               out.NbEvnts = str2num(tline((strfind(tline,':')+1):end));
           case 'Events Order'
               out.EOrder = str2num(tline((strfind(tline,':')+1):end));
           case 'Events Description'
               tline = fgetl(fid);
               nbParam = length(regexp(tline,'\t'));
               out.EDesc = zeros(out.NbEvnts, nbParam);
               for ind = 1:out.NbEvnts
                   tline = fgetl(fid);
                   Tmp = str2num(tline);
                   out.EDesc(ind,:) = Tmp(2 + (0:(nbParam - 1)));
               end
           case 'Use Null Condition?'
               answer = tline((strfind(tline,':')+1):end);
               if( strcmp(answer, 'Yes') )
                   out.NullCond = 1;
               else
                   out.NullCond = 0;
               end
           case 'Reference (X; 1 to 512)'
               answer = tline((strfind(tline,':')+1):end);
               answer = str2double(answer);
               out.RefX = answer;
           case 'Reference (Y; 1 to 512)'
               answer = tline((strfind(tline,':')+1):end);
               answer = str2double(answer);
               out.RefY = answer;
           case 'MM per pix'
               answer = tline((strfind(tline,':')+1):end);
               answer = str2double(answer);
               out.MMxPix = answer;
           case 'Pulsed'
               answer = tline((strfind(tline,':')+1):end);
               answer = str2double(answer);
               out.isPulsed = answer;
           case 'Pulse Freq'
               answer = tline((strfind(tline,':')+1):end);
               answer = str2double(answer);
               out.PulseFreq = answer;
           case 'Pulse Width'
               answer = tline((strfind(tline,':')+1):end);
               answer = str2double(answer);
               out.PulseWidth = answer;              
           otherwise

       end
    end
elseif( contains(tline, 'Mapping') )
    tline = fgetl(fid);
    out.Positions = [];
    while( ~feof(fid) )
        tline = fgetl(fid);
        if( contains(tline,':') )
            tag = tline(1:(strfind(tline,':')-1));
            idxS = regexp(tline,'\t') + 1;
            idxE = length(tline) - regexp(tline(end:-1:1), ' ');
            if( isempty(idxE) | (idxE < idxS) )
                idxE = length(tline);
            end
            Value = str2num(tline(idxS:idxE));
            switch(tag)
                case 'Laser Power' %Deprecated
                    out.LaserPower = Value;
                case 'OptoGen Power'
                    out.LaserPower = Value;       
                case 'Time On'
                    out.TimeON = Value;
                case 'Ref X Pos'
                    out.RefX = Value;
                case 'Ref Y Pos'
                    out.RefY = Value;
                case 'MM per Pix'
                    out.MMpPix = Value;
               case 'Spacing'
                   out.Spacing = Value;     
               case 'NbReps'
                   out.NbReps = Value;   
               case 'Events Order'
                   out.EOrder = str2num(tline((strfind(tline,':')+1):end));                   
            end
        else
            out.Positions = [out.Positions; str2num(tline)];
        end
    end
end
fclose(fid);

end