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
        FFTdata     % Structure containing the Phase and Amplitude movies for each direction.
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
                obj.splitMovInfo(i).timestamps = eventInfo.timestamps(idxStart:idxStop) -...
                    eventInfo.timestamps(idxStart); % Shift timestamps.
            end
        end
        
        function calculateFFT(method)
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

           
            
        end
        
        
        
    end
    
    methods (Access = private)
        function dataPerTrial = splitTrial(obj, idx)
            % This method splits the data in "splitMov" and reshapes it in a
            % 4D array. Missing data are imputed with the average trial
            % value.
            disp('Splitting data...')
            
            % Get the data:
            data = obj.splitMov{idx};
            % Infer the trial size in frames:
            trial_length = round((mean(diff(obj.splitMovInfo(idx).timestamps(~obj.splitMovInfo(idx).state)))) * obj.sr);
            
            
            
            
            preFr = round(sr*opts.preEventTime_sec);
            postFr = round(sr*opts.postEventTime_sec);
            centralFr = preFr + 1;
            len_trial = preFr + postFr;
            n_trial = sum(evDat.state == 1);
            timestamps = evDat.timestamps(evDat.state == 1);
            new_dims = {'E', 'Y', 'X','T'};
            [~, locB]= ismember(new_dims([2,3]), metaData.dim_names);
            % Create empty matrix:
            outData = nan([n_trial, szdat(locB), len_trial], 'single');
            
            
            
            
        end
        function averageData(obj, varargin)
            % This function averages the movies in "splitMov".
            % One can also normalize the data to obtain DeltaF/F.
            % The normalization is calculated by :
            %       "x - median(x_interstim) ./ median(x_interstim)"
            % Input:
            % normalize (bool) : Keyword argument to normalize or not the
            % averaged data.
            
            p = inputParser;
            addRequired(p, 'obj')
            addParameter(p,'normalize', false, @islogical)
            parse(p, obj, varargin{:});
            flag = p.Results.normalize;
            clear p
            
            for i = 1:4
                dat = obj.splitTrial(i)
                
            end
            
            
            
        end
    end
end

