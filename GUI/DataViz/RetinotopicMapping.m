classdef RetinotopicMapping
    % RETINOTOPICMAPPING calculates the retinotopic maps and Visual Sign Maps
    % widefield imaging data with a visual stimulus characterize by a
    % bar that drifts over the 4 cartesian directions (0,90 180 and 270
    % deg) N times.
    
    properties
        VSM_params  % Structure containig the parameters for Visual Sign Map calculation.
    end
    
    properties (SetAccess = private)
        splitMov   % cell array containing the movies split by direction.
        splitMovInfo % Structure array containing  eventID, state and shifted timestamps values.
        sr          % (num. scalar) Sample rate of the recording in Hz.
        powData % cell array containing the FFT amplitude movies split by direction.
        phiData % cell array containing the FFT phase movies split by direction.
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
        
        function calculateFFT(obj,varargin)
           % CALCULATEFFT will generate the FFT power and phase spectra of
           % each direction of the retinotopy data.
           % The calculation can be made using the FFT from the whole movie
           % excluding the interstim data ("classic" method) or from the normalized 
           % average movie ("AllenBrain" method). 
           % The "AllenBrain" method follows the procedures described in :
           %
           % Zhuang et al.‘An Extended Retinotopic Map of Mouse Cortex’.
           % ELife 6 (6 January 2017): e18372. https://doi.org/10.7554/eLife.18372.
           % 
           % Input:
           % method (str) : "classic" (default) OR "AllenBrain".

           p = inputParser;
           addRequired(p, 'obj')
           addOptional(p,'method', 'classic', @(x) ismember(x, {'classic', 'AllenBrain'}))
           parse(p, obj, varargin{:});
           method = p.Results.method;
           clear p
           for i = 1:4               
               if strcmp(method,'classic')
                   % Remove interstim sections of the movies:
                   onFr = [];
                   frIndx = round(obj.splitMovInfo(i).timestamps*obj.sr);
                   for j = 1:2:length(frIndx)-1
                       onFr = [onFr frIndx(j):frIndx(j+1)];
                   end
                   mov = obj.splitMov{i}(:,:,onFr);
%                     mov = obj.splitMov{i};
               else 
                   % Allen Brain Method
                   mov = obj.averageData(i);                   
               end
               
               % Calculate FFT:
               disp(['Calculating FFT for direction ' num2str(obj.splitMovInfo(i).eventID) ' deg...'])
               fftData = fft(mov,[],3); 
               clear mov
               app.powData{i} = abs(fftData);
               app.phiData{i} = mod(angle(fftData),2*pi); % From Zhuang et al., 2017
               clear fftData
           end
           disp('Done!')
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
            
             disp('Splitting data...')
            
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
            avg = avg - median(avg(:,:,indx_template<0), 3);            
        end
    end
end

