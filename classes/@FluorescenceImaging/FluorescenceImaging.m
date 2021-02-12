classdef FluorescenceImaging < Modality
    % "Fluorescence Imaging" class.
    %   This class deals with the fluorescence imaging data.
    properties
        NumberOfChannels % Number of Channels in the recording.
        Fluo_channel_file
    end
    properties (SetAccess = private)
        Outputs_Ptr_file
    end
    properties (Access = {?PipelineManager})
        
    end
    
    methods
        % Constructor
        %                 function obj = FluorescenceImaging()
        %
        %
        %                 end
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function run_ImagesClassification(obj, opts)
            % RUN_IMAGESCLASSIFICATION calls the function
            % IMAGESCLASSIFICATION from the IOI library (LabeoTech).
            arguments
                obj
                opts.BinningSpatial (1,1) {mustBeNumeric} = 1
                opts.BinningTemp (1,1) {mustBeNumeric} = 1
                opts.b_SubROI (1,1) {mustBeNumericOrLogical} = false
                opts.b_IgnoreStim (1,1) {mustBeNumericOrLogical} = true
            end
            % Calls function from IOI library. Temporary for now.
            ImagesClassification(obj.RawFolder, opts.BinningSpatial, opts.BinningTemp, opts.b_SubROI, opts.b_IgnoreStim)
            cd(obj.RawFolder)
            chanList = dir('*Chan.dat'); chanList = {chanList.name};
            for i = 1:length(chanList)
                chanName = chanList{i};
                switch chanName
                    case 'fChan.dat'
                        prefix = 'Fluo';
                        MetaDataFileName = 'Data_Fluo.mat';
                    case 'gChan.dat'
                        prefix = 'Green';
                        MetaDataFileName = 'Data_green.mat';
                    case 'rChan.dat'
                        prefix = 'Red';
                        MetaDataFileName = 'Data_red.mat';
                    case 'yChan.dat'
                        prefix = 'Amber';
                        MetaDataFileName = 'Data_yellow.mat';
                end
                
                newmDfile = strrep(chanName, '.dat', '_info.mat'); % Temporary...
                [~] = movefile(MetaDataFileName, newmDfile); % Temporary...
                a =  matfile(newmDfile, 'Writable', true);
                a.fileUUID = convertStringsToChars(matlab.lang.internal.uuid);
                a.Datatype = 'single';
                a.datName = 'data';
                eval(['obj.' prefix '_channel_file = fullfile(obj.RawFolder, ''' chanName ''');'])
                
            end
            a = matfile([obj.RawFolder 'fChan_info.mat']);
            obj.SampleRateHz = a.Freq*numel(chanList);
            obj.NumberOfChannels = numel(chanList);
        end
        
        function run_GSR(obj, datFileName)
            % Applies Global Signal Regression (GSR) to imaging recording
            % data in DATFILENAME.
            % Open memMapfile:
            mmData = mapDatFile(datFileName);
            % Load Data:
            data = flipud(rot90(mmData.Data.data));
            szData = size(data);
            data = reshape(data, [], szData(3), 1);
            % Calculate GSR:
            mData = mean(data,3);
            Sig = mean(data); 
            Sig = Sig / mean(Sig);
            X = [ones(szData(3),1), Sig'];
            B = X\data';
            A = X*B;
            clear X B Sig
            data = data - A';%Centre a 0
            data = data + mData; %Centre sur valeur moyenne constante.    
            data = reshape(data,szData);
                        
            % Generate .DAT and .MAT file Paths:
            datFile = fullfile(obj.SaveFolder, 'GSR.dat');
            % Save to .DAT file and create .MAT file with metaData:
            save2Dat(datFile, data);
            % Save Paths to object's properties:
            obj.GSR_file = datFile;
        end
        
        function run_TemporalFilter(obj, datFileName, opts)
            % Applies bandpass temporal filter to Imaging recordings DATA
            % in DATFILENAME.
            %   LOW_CUTOFF and HIGH_CUTOFF are the low and high cut-off
            %   frequencies of the bandpass filter.
            arguments
                obj
                datFileName {mustBeFile}
                opts.LowCutOff {mustBeNumeric} = 1/120 % Need to find a better way to choose Default value!
                opts.HighCutOff {mustBeNumeric} = 5  % Need to find a better way to choose Default value!
            end
            Freq =  obj.SampleRateHz/obj.NumberOfChannels;
            % Open memMapfile:
            mmData = mapDatFile(datFileName);
            data = mmData.Data.data;
            % Set temporal Filter params:
            szData = size(data);
            f = fdesign.lowpass('N,F3dB', 4, opts.LowCutOff, Freq);
            lpass = design(f,'butter'); lp_sosM = lpass.sosMatrix; lp_SV = lpass.ScaleValues;
            f = fdesign.lowpass('N,F3dB', 4, opts.HighCutOff, Freq);
            hpass = design(f,'butter'); hp_sosM = hpass.sosMatrix; hp_SV = hpass.ScaleValues;
            
            for ind = 1:szData(1)
                Sig = double(squeeze(data(ind,:,:)))';
                % Apply temporal filter:
                dl = single(filtfilt(lp_sosM, lp_SV, Sig));
                dh = single(filtfilt(hp_sosM, hp_SV, Sig));
                data(ind,:,:) = (dh./dl)';                
            end
             % Generate .DAT and .MAT file Paths:
            datFile = fullfile(obj.SaveFolder, 'TempFilt.dat');
             % Save to .DAT file and create .MAT file with metaData:
            save2Dat(datFile, data);
            obj.TemporalFilter_file = datFile;
        end
        
        function run_SeedPixelCorrelation(obj, datFileName)
            % This function reduces the matrix DATA in DATFILENAME to a 64 x 64, calculates
            % a pixel-wise temporal correlation.
            % Open memMapfile:
            mmData = mapDatFile(datFileName);
            % Load data:
            data = mmData.Data.data;
            % Calculate SeedPixel Correlation:
            A = imresize3(data, [64, 64, size(data,3)]);
            B = reshape(A, [], size(A,3))';
            [CM, P] = corr(B);
            CM = reshape(CM, [64 64 64^2]);
            P = reshape(P, [64 64 64^2]);
            clear A B data
            % Create MetaData structure:
            datFile = fullfile(obj.SaveFolder, 'SeedPxCorr.dat');
            metaDat = struct('datName', {'CM', 'P'}, 'datSize', {size(CM,[1 2]), size(CM,[1 2])},...
                    'datLength', {size(CM,3) size(CM,3)}, 'Datatype', {'single', 'single'}, 'datFile', datFile);
            % Save CM and METADAT to DATFILE:
            save2Dat(datFile, CM, 'metaData', metaDat)
            % Append "P" to DATFILE:
            save2Dat(datFile, P, 'flag', '-a')
            obj.SeedPixelCorr_file = datFile;
            
        end
        
        function viz_SeedPxCorr(obj)
            %
            % TEMPORARY FUNCTION, JUST TO SEE THE RESULTS OF
            % RUN_SEEDPIXELCORRELATION....
            %
            mmData = mapDatFile(obj.SeedPixelCorr_file);
            Fr = squeeze(mmData.Data.CM(:,:,1));
            figure('Position', [474 175 560 420]);
            imagesc(Fr); axis image;
            while 1
                [x,y] = ginput(1);
                ind = sub2ind([64 64], round(y), round(x));
                Fr = squeeze(mmData.Data.CM(:,:,ind));
                imagesc(Fr); axis image;
            end
        end
        
        function externalCallTo(obj, func2Call, param)
           % Dans Param:
           % Function a appeler
           % les channels a utiliser.
           %Initialization
           Out = [];
           
           %Building needed Paths:
           %Ouvrir fichier Ptr.txt et aller chercher les bon paths.
           %Reconstruire selon les ROOTs par ex:
           %du fichier, on a: Fluo_channel_file		{Proto_ROOT}\{Subject_ROOT}\{Acq_ROOT}\fChan.dat
           Str = '{Proto_ROOT}\{Subject_ROOT}\{Acq_ROOT}\fChan.dat';
           %Remplacer les {.} selon les parents de l<objet.
           
           eval(['Out = ' func2Call '(' param ');']); 
           
           Add2TxtFile...
        end
        %         function draw_manualROI(~, s_file, nROI)
        %             % Function that plots an image from the cortex and allows the
        %             % user to manually draw a number of ROIs (NROI).
        %             %   S_FILE points to a .MAT file containing a variable FRAME with
        %             %   an image (2D array) of the cortex.
        %
        %             m_savefile = 'Logical_mask.mat';
        %             load(s_file,'frame');
        %             ax = gca;
        %             imagesc(ax, frame); colormap('gray'); axis image;
        %             a = 1;
        %             bw = zeros(size(frame,1), size(frame,2), nROI);
        %             while a <= nROI
        %                 title(ax, ['Draw ROI ' num2str(a) '/' num2str(nROI)]);
        %                 roi = drawpolygon(ax);
        %                 bw(:,:,a) = createMask(roi);
        %                 a = a + 1;
        %             end
        %             logical_mask = any(bw,3);
        %             imagesc(ax,logical_mask); axis image
        %             title(ax,'Logical Mask');
        %             save(m_savefile, 'logical_mask');
        %         end
        %
        %         function draw_manualBregmaPos(~, s_file)
        %             % Uses Matlab bult-in function GINPUT to get the X,Y position
        %             % of Bregma.
        %             m_savefile = 'Bregma_coords.mat';
        %             load(s_file,'frame');
        %             ax = gca;
        %             imagesc(ax, frame); colormap('gray'); axis image; hold on
        %             title('Click on Bregma');
        %             [x,y] = ginput(1);
        %             x = round(x); y = round(y);
        %             pos = [x y 10 10];
        %             rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', 'r');
        %             text(x, y - 15, 'Bregma', 'HorizontalAlignment',  'left'); hold off
        %             title(ax, 'Done!');
        %             bregma_coords = [x y];
        %             save(m_savefile, 'bregma_coords');
        %         end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dummyMethodForTesting(obj)
            disp(class(obj))
            disp(['This is a dummy function #1 for testing of Fluorescence IMAGING ' obj.ID '!'])
        end
        function dummyMethodForTesting2(obj)
            disp(class(obj))
            disp(['This is a dummy function #2 for testing of Fluorescence IMAGING ' obj.ID '!'])
        end
        function dummyMethodForTesting3(obj)
            disp(class(obj))
            disp(['This is a dummy function #3 for testing of Fluorescence IMAGING ' obj.ID '!'])
        end
        function dummyMethodForTesting_WithError(obj)
            %%%ERROR%%
            a = class(obj);
            cell2mat(a);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end

