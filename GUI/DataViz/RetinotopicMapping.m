classdef RetinotopicMapping < handle
    % RETINOTOPICMAPPING calculates the retinotopic maps and Visual Sign Maps
    % widefield imaging data with a visual stimulus characterize by a
    % bar that drifts over the 4 cartesian directions (0,90 180 and 270
    % deg) N times.
    
    properties
        % Structure containig the parameters for Visual Sign Map calculation.
        VSM_params = struct('phaseMapFilterSigma', 0.5,...
                            'signMapFilterSigma', 8.0,....
                            'signMapThr', 0.4,...
                            'eccMapFilterSigma', 15.0,....
                            'splitLocalMinCutStep', 5.,...
                            'closeIter', 3,...
                            'openIter', 3,...
                            'dilationIter', 15,...
                            'borderWidth', 1,...
                            'smallPatchThr', 100,...
                            'visualSpacePixelSize', 0.5,...
                            'visualSpaceCloseIter', 15,...
                            'splitOverlapThr', 1.1,...
                            'mergeOverlapThr', 0.1);                               
    end
    
    properties (SetAccess = private)
        splitMov   % cell array containing the movies split by direction.
        splitMovInfo % Structure array containing  eventID, state and shifted timestamps values.
        sr          % (num. scalar) Sample rate of the recording in Hz.    
        AzimuthMap  % Structure containing the Azimuth Phase and Amplitude maps.
        ElevationMap % Structure containing the Elevation Phase and Amplitude maps.
        VSM         % Image of the Visual Sign Map calculated from the retinotopy data.
    end
    
    methods
        function obj = RetinotopicMapping(rawData, metaData,eventInfo)
            % RETINOTOPICMAPPING Construct an instance of this class
            %   This construction method splits the raw Data and stores the data's metadata.
            
            % Save sample rate from event Info file:
            obj.sr = metaData.Freq;
            % Split Raw movie:
            obj.splitMov = cell(1,4);
            obj.splitMovInfo = repmat(struct('eventID',[],'state',[],'timestamps',[]),1,4);
            dirList = 0:90:270; % We assume that there are 4 directions recorded in ascending order!
            for i = 1:4
                obj.splitMovInfo(i).eventID = dirList(i);
                % Splits the rawMovie for each direction.
                % get interstim time:
                dirIdx = eventInfo.eventID == dirList(i);
                % Find interstimulation time:
                onIdx = find(dirIdx & eventInfo.state,2,'first');
                offIdx = find(dirIdx & ~eventInfo.state,1,'first');
                interTm = (eventInfo.timestamps(onIdx(2)) - eventInfo.timestamps(offIdx));
                
                % Find first and last states:
                idxStart = find(dirIdx & eventInfo.state,1,'first');
                idxStop  = find(dirIdx & ~eventInfo.state,1,'last');
                
                % Calculate first and last frame indices:
                frStart = round((eventInfo.timestamps(idxStart) - interTm)* metaData.Freq);
                frStop = round(eventInfo.timestamps(idxStop)* metaData.Freq);
                
                % Split the Raw Data
                obj.splitMov{i} = rawData(:,:,frStart:frStop);
                % Update "state" and "timestamps" arrays:
                obj.splitMovInfo(i).state = eventInfo.state(idxStart:idxStop);
                obj.splitMovInfo(i).timestamps = interTm + eventInfo.timestamps(idxStart:idxStop) -...
                    eventInfo.timestamps(idxStart); % Shift timestamps.
            end
        end
        
        function out = calculateFFT(obj,movIndx,varargin)
            % CALCULATEFFT will generate the FFT power and phase maps OR
            % the FFT power spectrum of a given pixel from a sweep direction 
            % of the retinotopy data. 
            % The calculation can be made using the FFT from the whole movie
            % excluding the interstim data ("classic" method) or from the normalized
            % average movie ("AllenBrain" method).
            % The "AllenBrain" method follows the procedures described in :
            % Zhuang et al.‘An Extended Retinotopic Map of Mouse Cortex’.
            % ELife 6 (6 January 2017): e18372. https://doi.org/10.7554/eLife.18372.
            %
            % Inputs:
            % movIndx(int) : Index of the movie to be analysed.
            % pixel_coords(1x2 int) : X and Y coordinates of the pixel to
            %   calculate the FFT power spectrum.            
            % FFT_frequency(int) : Stimulation frequency (e.g. number of
            %   sweeps) + 1 to get the spatial maps of the FFT amplitude and
            %   phase. By default, it gives the 1st Harmonic.
            % method (str) : "classic" (default) OR "AllenBrain". The
            %   "classic" method removes the interstim frames of the movie
            %   and calculates the FFT from the whole movie. The "AllenBrain"
            %   method, calculates the FFT from the averaged movie normalized 
            %   as DeltaR (data - median(interstim).  
            %
            % Output:
            % out(1D or 3D num. matrix) : If the pixel_coordinates is provided,
            %   the output consists of the FFT power spectrum of the pixel.
            %   Otherwise, the output is a 3-D matrix containing the
            %   amplitude and phase 2-D maps of the given FFT frequency
            %   concatenated in the 3rd dimension.
            
            p = inputParser;
            addRequired(p, 'obj')
            addRequired(p, 'movIndx', @(x) isscalar(x) && ismember(x, 1:4));
            addParameter(p, 'pixel_coords', [], @(x) length(x)== 2 || isempty(x));
            addParameter(p, 'FFT_frequency', 0,...
                @(x) isscalar(x) & x >= 0);
            addParameter(p,'method', 'classic', @(x) ismember(x, {'classic', 'AllenBrain'}))
            parse(p, obj, movIndx, varargin{:});
            movIndx = p.Results.movIndx;
            px_XY = round(p.Results.pixel_coords); % Round to integer.
            freqFr = round(p.Results.FFT_frequency); % Round to integer.
            method = p.Results.method;            
            clear p
            % Default FreqFr to the 1st Harmonic:
            if freqFr == 0
                if strcmp(method, 'classic')
                    freqFr = 1 + sum(obj.splitMovInfo(movIndx).state);
                else
                    freqFr = 2;
                end
            end                                
            %            
            if strcmp(method,'classic')
                % Remove interstim sections of the movies:
                onFr = [];
                frIndx = round(obj.splitMovInfo(movIndx).timestamps*obj.sr);
                for j = 1:2:length(frIndx)-1
                    onFr = [onFr frIndx(j):frIndx(j+1)];%#ok
                end
                mov = obj.splitMov{movIndx}(:,:,onFr);                
            else
                % Allen Brain Method
                mov = obj.averageData(movIndx);                
            end            
            % Calculate FFT:
            disp(['Calculating FFT for direction ' num2str(obj.splitMovInfo(movIndx).eventID) ' deg...'])                       
            if isempty(px_XY)
                fftData = fft(mov,[],3);
                % Generate Amplitude/Phase maps for a given Frequency "freqFr":
                out = (abs(fftData(:,:,freqFr)) * 2) / size(mov,3); % From Zhuang et al., 2017
                out = cat(3,out,mod(-1*angle(fftData(:,:,freqFr)),2*pi)); % From Zhuang et al., 2017                
            else
                fftData = fft(mov(px_XY(1), px_XY(2),:),[],3);
                % Generate the Power FFT plot for a given pixel
                out = squeeze(abs(fftData) * 2); % From Zhuang et al., 2017
            end  
            disp('Done!');
        end
        
        function averageCardinalMaps(obj, varargin)
           % This method averages the FFT data from Top-Down and Left-Right
           % to generate the Azimuth and Elevation amplitude and phase
           % maps. Here, we calculate the average of the amplitude and
           % subtract the phase maps.
           % Input:
           % method(str): "classic" (default) OR "AllenBrain". See
           % docstring of method "calculateFFT" for details.
           
           p = inputParser;
           addRequired(p, 'obj');
           addParameter(p,'method', 'classic', @(x) ismember(x, {'classic', 'AllenBrain'}))
           addParameter(p, 'FFT_frequency', 0,@(x) isscalar(x) & x >= 0);
           parse(p, obj, varargin{:});
           method = p.Results.method;
           freqFr = round(p.Results.FFT_frequency); % Round to integer.
           clear p
           
           % Average azimuth:
           disp('Averaging azimuth maps...')
           map_0 = obj.calculateFFT(1,'method',method, 'FFT_frequency', freqFr);
           map_180 = obj.calculateFFT(3,'method',method, 'FFT_frequency', freqFr);
           out = zeros(size(map_0,1), size(map_0,2), 4);
           out(:,:,1) = mean(cat(3,map_0(:,:,1), map_180(:,:,1)),3); % Average amplitude;
           out(:,:,2) = (map_0(:,:,2) - map_180(:,:,2))/2; % Subtraction of the phase;
           disp('Done!')
           % Average elevation:
           disp('Averaging elevation maps...')
           map_90 = obj.calculateFFT(2,'method',method, 'FFT_frequency', freqFr);
           map_270 = obj.calculateFFT(4,'method',method, 'FFT_frequency', freqFr);
           out(:,:,3) = mean(cat(3,map_90(:,:,1), map_270(:,:,1)),3); % Average amplitude;
           out(:,:,4) = (map_90(:,:,2) - map_270(:,:,2))/2; % Subtraction of the phase; 
           disp('Done!')
           
           obj.AzimuthMap = out(:,:,[1 2]);
           obj.ElevationMap = out(:,:,[3 4]);
        end
        
        function out = genVSM(obj,varargin)
            
            
            
            if isempty(obj.AzimuthMap)
                disp('Run "averageCardinalMaps" to create Azimuth and Elevation maps first!')
                return
            end
            % Filter Maps spatially:
            obj.applyGaussFilt;
            % Calculate visual sign map:
            obj.calc_visualSign(obj.AzimuthMap(:,:,2), obj.ElevationMap(:,:,2));
            
        end
        
        function drawROI(obj);
            
        end
        
    end
    
    methods (Access = private)
        
        function avg = averageData(obj, indx_mov)
            % This function averages the movies in "splitMov"
            % divided by the median interstim value to get a DeltaR movie:
            %       "x - median(x_interstim)"
            % Input:
            %   indx_mov (num. scalar) : index of the movie stored in
            % "splitMovInfo"
            % Output:
            %   avg: (3D movie) : average DeltaR movie.
            
            disp('Averaging movie...')
            
            % Get the data:
            data = obj.splitMov{indx_mov};
            % Infer the trial size from first sweep:
            lenTrial = round(obj.splitMovInfo(indx_mov).timestamps(2)*obj.sr);
            indx_template = (1:lenTrial) - round(obj.splitMovInfo(indx_mov).timestamps(1)*obj.sr);
            % Get the onset frame of each stimulus:
            frOnList = round(obj.splitMovInfo(indx_mov).timestamps(obj.splitMovInfo(indx_mov).state)...
                *obj.sr);
            avg = zeros([size(data,1),size(data,2), lenTrial], 'single');
            % Sum all trials
            for i = 1:length(frOnList)
                frames = indx_template + frOnList(i);
                % control for first sweep shorter than the Trial size:
                if min(frames) < 1
                    disp(['Missing frames found in first sweep of direction '...
                        num2str(obj.splitMovInfo(indx_mov).eventID)]);
                    disp('Padding missing frames with the average trial values...')
                    % Pad missing frames with the average of the trial
                    tmp_trial = data(:,:,frames(frames>0));
                    tmp_avg = mean(tmp_trial,3);
                    tmp_avg = cat(3,repmat(tmp_avg,1,1,sum(frames<1)),tmp_avg);
                    avg = avg + tmp_avg;
                else
                    avg = avg + data(:,:,frames);
                end
            end
            avg = avg/sum(obj.splitMovInfo(indx_mov).state);
            % Calculate the DeltaR
            avg = avg - median(avg(:,:,indx_template < 0), 3);
        end
        
        function applyGaussFilt(obj)
           % This method applies spatial gaussian filter to the amplitude
           % and phase components of the Azimuth and Elevation maps.
           disp('Filtering maps...')
           for i = 1:2
            obj.AzimuthMap(:,:,i) = imgaussfilt(obj.AzimuthMap(:,:,i),...
                obj.VSM_params.phaseMapFilterSigma);
            obj.ElevationMap(:,:,i) = imgaussfilt(obj.ElevationMap(:,:,i),...
                obj.VSM_params.phaseMapFilterSigma);
           end
           disp('Gaussian filter applied to Azimuth and Elevation maps');
                      
        end
        
        function calc_visualSign(obj, phaseAz,phaseEl)
            % This method calculates the visual sign maps of the Azimuth
            % and Elevation maps.
            disp('Calculating visual sign map...');
            
            % Calculate phase map gradients
            [gradAz1, gradAz2] = gradient(phaseAz);
            [gradEl1, gradEl2] = gradient(phaseEl);
            %
            gradDirAz = zeros(size(gradAz1));
            gradDirEl = zeros(size(gradEl1));
            %
            for i = 1:size(phaseAz,1)
                for j = 1:size(phaseAz,2)
                    gradDirAz(i,j) = atan2(gradAz2(i,j),gradAz1(i,j));
                    gradDirEl(i,j) = atan2(gradEl2(i,j),gradEl1(i,j));
                end
            end
            clear i j
            vsm = sin(angle(exp(1j.*gradDirAz).*exp(-1j.*gradDirEl)));
            % Filter map
            obj.VSM = imgaussfilt(vsm, obj.VSM_params.signMapFilterSigma);
                        
        end
        
        
    end
end

