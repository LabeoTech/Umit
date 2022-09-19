function TwoPtsAnalysis(FolderPath)

%%HARCODED VALUES:
% Change these variables depending on what is required as output
%
%1- OutFormat: Type of output file (txt, cvs, binary, mat, etc.)
OutFormat = 'csv'; %Options: 'txt', 'csv', 'mat', 'bin', 'xls';
%2- Output Infos: What should be outputed for further analysis (raw,
%average over repetition, average per groups, standard deviation, etc)
OutInfos = {'Raw'};
%
% GraphOut: Which graph to be generated
GraphOut = {'Raw'}; %Options: 'Raw', 'AveReps', 'AveGroup'

%%%

A = mfilename;
A = dir([A '.m']);
A.date
% Step 1: Inputs validation
disp(['Two Points Analysis of: ' FolderPath]);
if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

%Step 2: Data validation
InfoStim = Read_OptoGenParams_File(FolderPath);
AISampleRate = 10000;
AINChannels = 11;

%Step 3: Data opening
disp('Reading Analog In files.');
aiFiles = dir([FolderPath 'ai_*.bin']);
AnalogIn = [];
for indF = 1:length(aiFiles)
    fid = fopen([FolderPath aiFiles(indF).name]);
    fseek(fid,5*4,'bof');
    dat = fread(fid,inf,'double');
    dat = reshape(dat, AISampleRate, AINChannels, []);
    dat = permute(dat,[2 1 3]);
    dat = reshape(dat,11,[]);
    AnalogIn = [AnalogIn, dat];
end
clear dat indF aiFiles fid

%Step 4: Stim Events segmentation
disp('Events detection');
StimE = find((diff(AnalogIn(2,:),1,2) > 2.5) & (AnalogIn(2,2:(end))>0));
[~,index] = sort(InfoStim.EOrder);
StimOrder = StimE(index);

% IMPORTANT: the following lines determine which analog inputs must be
% considered for each paw (Channels 4, 5, 6 correspond to AI 1, 2 and 3 on
% the device. Channels 8, 9 and 10 correspond to AI 5, 6 and 7 on the device) 
ChanIdxPaw1 = 4:6;
ChanIdxPaw2 = 8:10;

disp('Signal segmentation');
Sig = zeros(2, InfoStim.NbEvnts*InfoStim.NbReps, 4500,'single');
Offset = -(50 - InfoStim.TOffset)*10;
for indS = 1:size(StimE,2)
    Sig(1, indS, :) = sqrt(sum(AnalogIn(ChanIdxPaw1, StimOrder(indS) + Offset + (0:4499)).^2,1));
    Sig(2, indS, :) = sqrt(sum(AnalogIn(ChanIdxPaw2, StimOrder(indS) + Offset + (0:4499)).^2,1));
end
Sig = bsxfun(@rdivide, Sig, mean(Sig,3));
Sig = reshape(Sig, 2, InfoStim.NbReps, InfoStim.NbEvnts,[]);

% Step 5: averaging over repetitions
mSig = squeeze(mean(Sig,2));
sSig = squeeze(std(Sig,0,2));
if( InfoStim.NbEvnts == 1 )
    mSig = permute(mSig, [1 3 2]);
    sSig = permute(sSig, [1 3 2]);
end

%Step 6: Z score for thresholding
zSig = bsxfun(@rdivide,bsxfun(@minus, mSig, mean(mSig(:,:,1:500),3)),...
    std(mSig(:,:,1:500),0,3));

%Step 7: Time vector creation:
TZeros = 50;
timeVect = linspace(-TZeros,450-TZeros, 4500); %in msec 
Results = matfile([FolderPath 'Results.mat'], 'Writable', true);
Results.timeVect = timeVect';

%Step 8: Graphs
if( any(contains(GraphOut, 'Raw')) ) %Generate a figure with all raw curves superposed:
    hf1 = figure; 
    plot(timeVect, reshape(Sig(1,:,:,:),[],4500));
    title('Raw Paw #1 (Channels 1-2-3)');
    saveas(hf1, [FolderPath 'RawDat_123.png']);
    close(hf1);
    hf2 = figure;
    plot(timeVect, reshape(Sig(2,:,:,:),[],4500));
    title('Raw Paw #2 (Channels 5-6-7)');
    saveas(hf2, [FolderPath 'RawDat_567.png']);
    close(hf2);    
end
if( any(contains(GraphOut, 'AveReps')) ) %Generate a figure with averages (events):
    hf1 = figure; 
    plot(timeVect, reshape(mean(Sig(1,:,:,:),2),[],4500));
    title('Events Average Paw #1 (Channels 1-2-3)');
    saveas(hf1, [FolderPath 'AvRepsDat_123.png']);
    close(hf1);
    hf2 = figure;
    plot(timeVect, reshape(mean(Sig(2,:,:,:),2),[],4500));
    title('Events Average Paw #2 (Channels 5-6-7)');
    saveas(hf2, [FolderPath 'AvRepsDat_567.png']);
    close(hf2);    
end
if( any(contains(GraphOut, 'AveGroup')) ) %Generate a figure with averages (group):
    NbVars = size(InfoStim.EDesc,2);
    tags = {'Laser Power Stim', 'Laser Power Pre-cond', 'InterStim delay'};
    for indV = 1:NbVars
        [~,~,idx] = unique(InfoStim.EDesc(:,indV));
        idx = cell2mat(arrayfun(@(X) (idx == X)', unique(idx),'UniformOutput', false));
        tmp = squeeze(mean(Sig(1,:,:,:),2));
        hf = figure; hold('on');
        title(tags{indV});
        for indG = 1:size(idx,1)
            plot(timeVect, mean(tmp(idx(indG,:),:),1));
        end
        if( (NbVars == 3)&(indV == 1) )
            legend(int2str(InfoStim.OptoGenPwr'))
        elseif( (NbVars == 3)&(indV == 2) )
            legend(int2str(InfoStim.PreCondPwr'))
        elseif( (NbVars == 3)&(indV == 3) )
            legend(int2str(InfoStim.IStimDelay'))
        end
        saveas(hf, [FolderPath ['Av_' tags{indV} 'Paw1_123.png']]);
        close(hf);
        tmp = squeeze(mean(Sig(2,:,:,:),2));
        hf = figure; hold('on');
        title(tags{indV});
        for indG = 1:size(idx,1)
            plot(timeVect, mean(tmp(idx(indG,:),:),1));
        end
        if( (NbVars == 3)&(indV == 1) )
            legend(int2str(InfoStim.OptoGenPwr'))
        elseif( (NbVars == 3)&(indV == 2) )
            legend(int2str(InfoStim.PreCondPwr'))
        elseif( (NbVars == 3)&(indV == 3) )
            legend(int2str(InfoStim.IStimDelay'))
        end
        saveas(hf, [FolderPath ['Av_' tags{indV} 'Paw2_567.png']]);
        close(hf);
    end
end

%Outputs:
% OutInfos options: 'Raw'
OutTable = table(timeVect','VariableNames',{'Time'});
if( any(contains(OutInfos,'Raw')) )
    tmp = array2table(reshape(permute(Sig,[4 2 1 3]),4500, []));
    Header = arrayfun(@(X) ['EvntID' int2str(floor(X/(InfoStim.NbReps + 2))+1)...
        '_Paw' int2str(ceil((mod(X-1,(InfoStim.NbReps + 2))+1)/2))...
        '_Rep' int2str(ceil((mod(X-1,(InfoStim.NbReps))+1)))],...
        1:size(tmp,2), 'UniformOutput', false);
    tmp.Properties.VariableNames = Header;
    OutTable = [OutTable, tmp];
    clear Header tmp;
end

% Saving data:
if( any(contains(OutFormat, 'txt')) )
    writetable(OutTable,[FolderPath 'RawData.txt']);
end
if( any(contains(OutFormat, 'csv')) )
    writetable(OutTable,[FolderPath 'RawData.csv']);
end
if( any(contains(OutFormat, 'xls')) )
    writetable(OutTable,[FolderPath 'RawData.xls']);
end
if( any(contains(OutFormat, 'mat')) )
    save([FolderPath 'RawData.mat'],'OutTable');
end
if( any(contains(OutFormat, 'bin')) )
    fid = fopen([FolderPath 'RawData.bin'],'w');
    fwrite(fid, table2array(OutTable), 'single');
    fclose(fid);
end



end
