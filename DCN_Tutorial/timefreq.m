% INPUT
data = EEG.data; % data = data matrix; channel x time/data points x trials
srate = EEG.srate; % srate = sampling rate;
% times = time vector / data time window;
% baseline = [xx xx]; time range/window to be used for baseline normalization
% min_freq = lowest frequency to analysis;
% max_freq = highest frequency to analysis;

% OUTPUT
% timefreq_data = frequency x time x channel matrix;
% times = time vector;
% frex = frequency vector;

% frequencies vector
frex = linspace(min_freq, max_freq, num_frex);
frequency = frex;

% cycles vector
cylvec = linspace(range_cycles(1),range_cycles(end), num_frex)./ (2*pi*frex);

%% wavelet parameters
wavtime = -1:1/srate:1; % length of wavelet
half_wave = (length(wavtime)-1)/2;

%% FFT parameters
nWave = length(wavtime);
nData = size(data, 2);
nConv = nWave + nData - 1;

%% Prepare data for wavelet
%Resting Data - only one trial type
%All to all connectivity
if ConnectType == 0
    %Loop through epochs and pull only epochs with appropriate number of
    %channels
    epochs2remove = [];
    for e=1:size(data,3) %begin loop through epochs
        %Only keep epoch if none of channels are NaNs
        for elec = 1:nbchan %begin loop through electrodes
            %check if the electrode is NaN within that epoch
            if sum(isnan(EEG.data(elec,:,e)))
                %add that epoch to vector
                epochs2remove = [epochs2remove e];
            else
                %do nothing
            end %end if loop for if channel within epoch is NaN
        end %end loop through electrodes
    end %end loop through epochs
    %remove all epochs with NaNs for channels of interest
    data(:,:,unique(epochs2remove))=[];
    %Seed-based connectivity
elseif ConnectType == 1
    epochs2remove=[];
    for e=1:size(data,3) %begin loop through epochs
        %Channel indices
        %concatenate Elecs4Connect and Seed to have vector of all
        %channels of interest
        Elecs4SeedBased = [Seed Elecs4Connect];
        %loop through and find indices of these electrodes
        for i=1:length(Elecs4SeedBased)
            Chan_idx (i)= find(strcmp({EEG.chanlocs.labels}, Elecs4SeedBased{i}));
        end
        %do any of these channels have NaNs in them?
        if min(sum(isnan(EEG.data(Chan_idx,:,e)))) > 0
            %if so, add that epoch to vector of epochs to be removed
            epochs2remove = [epochs2remove e];
        else
            %do nothing
        end %end if loop for if channel within epoch is NaN
    end %end loop through epochs
    %remove all epochs with NaNs for channels of interest
    data(:,:,unique(epochs2remove))=[];
    %keep only channels of interest
    data = data(Chan_idx,:,:);
end % end if all-to-all connectivity for rest

fprintf('\n\n\n*** Calculating TF for subject %d (%s) ***\n\n\n', sub, subject);

%% Run wavelet convolution
if RestorEvent ==1
    
    TrialNums.subject = subject;
    TrialNums.TrialNum = size(data,3);
    
    if size(data,3)>= mintrialnum
        timefreqs = zeros(length(frex), nData, size(data, 3));
        timefreq_strial = zeros(length(frex), nData, size(data, 3), size(data,1));
        for ch=1:size(data, 1) % Loop through all channels
            for fi=1:length(frex) % loop through all frequencies
                trial_conv = zeros(nData, size(data, 3));
                
                %% Create wavelet
                wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
                waveletX = fft(wavelet, nConv); % fft of wavelet
                waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
                
                %% Loop through all trials
                for trl=1:size(data, 3)
                    temp_data = fft(squeeze(data(ch,:,trl)), nConv);
                    
                    %% run convolution
                    temp_conv = ifft(waveletX .* temp_data);
                    temp_conv = temp_conv(half_wave+1:end-half_wave);
                    
                    %% data for phase analysis
                    trial_conv (:,trl) = temp_conv;
                    
                    %% compute power
                    trial_temppow (:,trl) = abs(temp_conv(half_wave+1:end-half_wave)).^2;
                end
                %Build dataset for phase data
                phase_data(fi,:,:) = trial_conv;
                        
                %Build dataset for TF
                timefreqs(fi,:,:) = trial_temppow;

                
            end
            %Build dataset for phase data with channels
            phase_data_strial_ch(:,:,:,ch) = phase_data;
            
            %Build dataset for TF with channels
            timefreqs_strial_ch(:,:,:,ch) = timefreqs;
        end %end loop through channels
        
    %Trial average timefrequency data
    timefreqs_data = squeeze(mean(timefreqs_strial_ch,3));
    
    else
        data2save='';
        save([save_location subject(1:end-4) '_notenoughdata.mat'],data2save)
    end %end if there's enough trials
    
    
    
    %Baseline Correction
    if BaselineCorrect == 1
        
        %% baseline time indices
        basetimeidx   = dsearchn(EEG.times', BaselineTime');
        
        Baseline = squeeze(mean(timefreqs_data(:,basetimeidx(1):basetimeidx(end),:),2));
        %loop through samples
        for t=1:size(timefreqs_data,2)
            timefreqs_baselinecorr(:,t,:) = 10*log10( squeeze(timefreqs_data(:,t,:)) ./ Baseline);
        end %end loop through frequencies
        
        if Downsample ==0
            %save out trial averaged and baseline corrected TF data for this subject for this condition
            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_baselinecorrected'];
            save(save_data, 'timefreqs_baselinecorr', 'frequency', 'time','channel_location', '-v7.3');
        elseif Downsample ==1
            
            %Downsample to 125Hz
            timefreqs_baselinecorr = timefreqs_baselinecorr(:,1:2:size(timefreqs_baselinecorr,2),:);
            
            %Downsample time variable
            ds_time = downsample(time,2);
            
            %save out trial averaged, downsampled, baseline corrected TF data for this subject for this condition
            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_baselinecorrected'];
            save(save_data, 'timefreqs_baselinecorr', 'frequency', 'ds_time','channel_location', '-v7.3');
        end
        
    elseif BaselineCorrect == 0
        
        if Downsample ==0
            %save out trial averaged and non-baseline corrected TF data
            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection'];
            save(save_data, 'timefreqs_data', 'frequency', 'time','channel_location', '-v7.3');
        elseif Downsample ==1
            %Downsample to 125Hz
            timefreqs_data = timefreqs_data(:,1:2:size(timefreqs_data,2),:);
            
            %Downsample time variable
            ds_time = downsample(time,2);
            
            %save out trial averaged and non-baseline corrected TF data
            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection'];
            save(save_data, 'timefreqs_data', 'frequency', 'ds_time','channel_location', '-v7.3');
        end
    end %end if you want baseline correction
    
    
    
    
    
elseif RestorEvent==0
    for cond=1:length(Conds)
        
        %pull out condition data
        cond_data = []; cond_EEG=[];
        cond_EEG = pop_selectevent( EEG, 'Condition',(Conds{cond}), 'deleteevents','on','deleteepochs','on','invertepochs','off');
        cond_data = eeg_checkset( cond_EEG.data );%add pulling out condition specific data
         
        %Create TrialNums structure to be saved out
        TrialNums( sub+(cond-1)+((sub-1)*(length(Conds)-1)) ).subject = subject;
        TrialNums( sub+(cond-1)+((sub-1)*(length(Conds)-1)) ).condition = Conds{cond};
        TrialNums( sub+(cond-1)+((sub-1)*(length(Conds)-1)) ).TrialNum = size(cond_data,3);
        
        %initialize matrices
        timefreqs = zeros(length(frex), nData, size(cond_data, 3));
        phase_data = zeros(length(frex), nData, size(cond_data, 3));
        phase_data_strial_ch = zeros(length(frex),  nData, size(cond_data, 3), size(cond_data,1));
        timefreqs_strial_ch = zeros(length(frex),  nData, size(cond_data, 3), size(cond_data,1));
        timefreqs_data = zeros(length(frex),  nData,  size(cond_data,1));
        
        if size(cond_data,3)< mintrialnum
            data2save='';
            save([save_location subject(1:end-4) '_notenoughdata.mat'],data2save)
            break %go onto next subject if one condition doesn't have enough trials
        else
%             if exist([save_location subject(1:end-4) '_notenoughdata.mat'],'file')==2
%                 continue %skip any condition if the participant has already been deemed to not have enough data
%             else
                for ch=1:size(cond_data, 1) % Loop through all channels
                    for fi=1:length(frex) % loop through all frequencies
                        trial_conv = zeros(nData, size(cond_data, 3));
                        trial_temppow = zeros(nData, size(cond_data, 3));
                        
                        %% Create wavelet
                        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
                        waveletX = fft(wavelet, nConv); % fft of wavelet
                        waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
                        
                        %% Loop through all trials
                        for trl=1:size(cond_data, 3)
                            temp_data = fft(squeeze(cond_data(ch,:,trl)), nConv);
                            
                            %% run convolution
                            temp_conv = ifft(waveletX .* temp_data);
                            temp_conv2 = temp_conv(half_wave+1:end-half_wave);
                            
                            %% data for phase analysis
                            trial_conv (:,trl) = temp_conv2;
                            
                            %% compute power
                            trial_temppow (:,trl) = abs(temp_conv(half_wave+1:end-half_wave)).^2;
                        end
                        
                        %Build dataset for phase data
                        phase_data(fi,:,:) = trial_conv;
                        
                        %Build dataset for TF
                        timefreqs(fi,:,:) = trial_temppow;
                        
                    end
                    
                    %Build dataset for phase data with channels
                    phase_data_strial_ch(:,:,:,ch) = phase_data;
                    
                    %Build dataset for TF with channels
                    timefreqs_strial_ch(:,:,:,ch) = timefreqs;
                    
                    
                end %end loop through channels
                
                %Trial average timefrequency data
                timefreqs_data = squeeze(mean(timefreqs_strial_ch,3));

%                     if TF_singletrial_save ==1
%                         %Downsample - NEED TO FIX THIS
%                         timefreqs_strial_ds = timefreqs_strial_ch(:,1:5:length(nData),:,:);
%                         time = time(1:5:length(nData));
% 
%                         %Save data
%                         save_data =[save_location, subject(1:end-4),'_singletrialTF_', 'condition',Conds{cond},'_',DatasetName];
%                         save(save_data, 'timefreqs_strial_ds', 'frequency', 'time', '-v7.3');
%                     else
%                         %don't save anything out
%                         timefreqs_strial_ch =[]; %blank out this large matrix
%                     end
                    
                    %Baseline Correction
                    if BaselineCorrect == 1
                        
                        %% baseline time indices
                        basetimeidx   = dsearchn(EEG.times', BaselineTime');
                            
                        Baseline = squeeze(mean(timefreqs_data(:,basetimeidx(1):basetimeidx(end),:),2));
                        %loop through samples    
                        for t=1:size(timefreqs_data,2)
                            timefreqs_baselinecorr(:,t,:) = 10*log10( squeeze(timefreqs_data(:,t,:)) ./ Baseline);
                        end %end loop through frequencies
                        
                        if Downsample ==0
                            %save out trial averaged and baseline corrected TF data for this subject for this condition    
                            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_baselinecorrected_', 'condition',Conds{cond}];
                            save(save_data, 'timefreqs_baselinecorr', 'frequency', 'time','channel_location', '-v7.3');
                        elseif Downsample ==1
                            
                            %Downsample to 125Hz
                            timefreqs_baselinecorr = timefreqs_baselinecorr(:,1:2:size(timefreqs_baselinecorr,2),:);
                            
                            %Downsample time variable
                            ds_time = downsample(time,2);
                            
                            %save out trial averaged, downsampled, baseline corrected TF data for this subject for this condition
                            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_baselinecorrected_', 'condition',Conds{cond}];
                            save(save_data, 'timefreqs_baselinecorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                        end
                            
                    elseif BaselineCorrect == 0
                        
                        if Downsample ==0
                            %save out trial averaged and non-baseline corrected TF data
                            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection_', 'condition',Conds{cond}];
                            save(save_data, 'timefreqs_data', 'frequency', 'time','channel_location', '-v7.3');
                        elseif Downsample ==1
                            %Downsample to 125Hz
                            timefreqs_data = timefreqs_data(:,1:2:size(timefreqs_data,2),:);
                            
                            %Downsample time variable
                            ds_time = downsample(time,2);
                            
                            %save out trial averaged and non-baseline corrected TF data
                            save_data =[save_location, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection_', 'condition',Conds{cond}];
                            save(save_data, 'timefreqs_data', 'frequency', 'ds_time','channel_location', '-v7.3');
                        end
                    end %end if you want baseline correction
                    
                    %Blank out power data over trials after saved to save memory space
                     timefreqs_data=[]; timefreqs_baselinecorr = [];
                    
                    %Inter-trial phase synchrony
                    if ITPS_calc == 1
                            cd(scripts_location)
                            ITPS;
                            ITPS_all=[]; %blank out matrix when done
                    else
                        %do nothing
                    end
                    
                    
                    %Inter-channel calculations
                    if ICPS_calc == 1
                            cd(scripts_location)
                            if ICPS_or_wPLI == 1 %calculate coherence
                                ICPS;
                                ICPS_all=[]; %blank out matrix when done
                            elseif ICPS_or_wPLI == 0 %calculate wPLI
                                wPLI;
                                wPLI_all=[]; %blank out matrix when done
                            end
                    else
                        %do nothing
                    end
                    
        end %end if there is enough data for this condition
    end %end loop through conditions
end %end if rest or event

