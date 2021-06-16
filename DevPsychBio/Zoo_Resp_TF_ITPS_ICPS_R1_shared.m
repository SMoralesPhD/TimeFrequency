%% Perform time frequency calculation
% This script performs the TF analyses, showed in Morales et al. (under review)
% These analyses include the computation of TF power, ITPS, and ICPS. 
% This script will save a separate file per participant.
% This scripts and less organized and are less intuitive than forthcoming
% efforts to conduct these same analyses in a more costumizable form.
% Also, for a tutorial, please see Bowers & Morales in prep. 
% The files produced by this script can be then plotted in matlab using the
% "compiling" scripts and then read into R using the ECHO_ERN_ERP_TF_ITPS_ICPS_R1.Rmd 

%%
clear
clc;

%% Adding this path for Ludwig
addpath('/data/cdl/Toolbox/eeglab13_6_5b/')

%% Locations
Data_Location = '/data/cdl/Projects/ECHO/CSD_tranformed_data/Zoo_3s/';
Save_Data = '/data/rosalind/Dropboxes/moraless/ECHO/ITPS_ICPS_redo/';

dataset_name = '_Zoo_Resp_TF_ITPS_100sub'; % '_stm_TimeFreqs'
%%

%% Subject list
cd(Data_Location)
subnum=dir('*CSD.set');
sub_list={subnum.name};
for i =1:length(sub_list)
    sub = sub_list{i};
    subject_list{i}= sub(1:15);
end

%%
% Reading in files that are included in the paper based on our exclusion
% criteria
BehData = importdata('/data/cdl/Projects/ECHO/Descriptives/ECHO_ERN_IDS_Paper_01_06_2021.csv');
subject_list = intersect(subject_list, BehData); % Deleting participants on the list
conditions = 2; % This is just for Go vs. No-Go
cutoff = 5; % Maybe the actual names to pass them for naming and stuff later


%% Events
event_markers = {'resp'}; % 'Stm+' or 'resp'

%% Change current directory
cd (Save_Data)

%% define baseline time window
baseline = [-300 -100]; 

%% frequency parameters
min_freq =  1;
max_freq = 30;
num_frex = 60;

% frequencies vector - this can be in linear or log space
frex = logspace(log10(min_freq), log10(max_freq), num_frex);

%% wavelet cycles - variable : min 3 max 10
range_cycles = [ 3 10 ];
% cycles vector
cylvec = logspace(log10(range_cycles(1)),log10(range_cycles(end)), num_frex)./ (2*pi*frex);

%% Running EEGLAB
eeglab

%% Loop through all subjects
for sub=1:length(subject_list)
    
    % load data
    subject = subject_list{sub};
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', subject);
    
    % Checking if the file already exists and if it does, skipping 
    cd(Save_Data);
    if exist([subject '_Zoo_Resp_TF_ITPS_100sub.mat'],'file') > 0
        do nothing
    else 
        TrialNums(sub).Subject = subject_list(sub);  % Creating a structure to save the number of trials for each condition
        
        
        EEG=pop_loadset('filename',[subject '_Zoo_Epoched_Matched_CSD.set'], 'filepath', Data_Location);
        EEG = eeg_checkset( EEG );
        
        %% Change sampling rate to 250Hz
        EEG = eeg_checkset( EEG );
        EEG = pop_resample( EEG, 250);
        
        %% Initialize variables for separate conditions
        EEGa = [];
        EEGb = [];
        %% Defining the conditions
        % Checking that we have enough epochs per condition before selecting epochs
        Cond1Num = length(find(strcmp({EEG.event.type},'resp') & [EEG.event.Accuracy] == 1 & [EEG.event.GoNogo] == 1 & [EEG.event.Bad] == 0));
        Cond2Num = length(find(strcmp({EEG.event.type},'resp') & [EEG.event.Accuracy] == 0 & [EEG.event.GoNogo] == 2 & [EEG.event.Bad] == 0));
        
        if Cond1Num > 0
            % Response Correct
            EEGa = pop_selectevent( EEG, 'StimType','Resp','GoNogo',1,'Accuracy',1, 'Bad', 0, 'deleteevents','on','deleteepochs','on','invertepochs','off');
            EEGa = eeg_checkset( EEGa );
            TrialNums(sub).Correct = length(EEGa.event); % Getting the number of trials
        elseif Cond1Num == 0
            TrialNums(sub).Correct = 0; % Setting the number of trials to zero
        end
        
        if Cond2Num > 0
            % Response Incorrect
            EEGb = pop_selectevent( EEG, 'StimType','Resp','GoNogo',2,'Accuracy',0, 'Bad', 0, 'deleteevents','on','deleteepochs','on','invertepochs','off');
            EEGb = eeg_checkset( EEGb );
            TrialNums(sub).Error = length(EEGb.event); % Getting the number of trials
        elseif Cond2Num == 0
            TrialNums(sub).Error = 0; % Setting the number of trials to zero
        end
        
        %% Get some data parameters
        channel_location = EEG.chanlocs;
        
        %% wavelet parameters
        srate=EEG.srate; % sampling rate
        wavtime = -.5:1/srate:.5; % length of wavelet
        half_wave = (length(wavtime)-1)/2;
        
        %% FFT parameters
        nWave = length(wavtime);
        nData = EEG.pnts;
        nConv = nWave + nData - 1;
        
        %% baseline time indices
        basetimeidx   = dsearchn(EEG.times', baseline'); % baseline indicies
        
        %% initialize output time-frequency data
        correct_tf_data = zeros( length(frex), nData, EEG.nbchan);
        error_tf_data = zeros( length(frex), nData, EEG.nbchan);
        
        %% initialize output ITPS data - freq x Data x channels
        correct_ITPS_nosubsamp_baseCorr = zeros( length(frex), nData, EEG.nbchan);
        error_ITPS_nosubsamp_baseCorr = zeros( length(frex), nData, EEG.nbchan);
        
        %blank out phase data
        correct_phase_data=[]; error_phase_data=[];
        
        %% Check if there are enough trials before convolution
        if TrialNums(sub).Correct > cutoff && TrialNums(sub).Error > cutoff
            TrialNums(sub).IncludedERP = 1;
            
            %% Run wavelet convolution
            for ch=1:EEG.nbchan % Loop through all channels
                fprintf('\n*** Processing chan %d for Subject %s ***', ch, subject);
                
                % Blanking and initializing structures that are used in the loop
                timefreqs = zeros(length(frex), nData);
                         
                correct_conv_trials=[]; error_conv_trials=[];
                ITPS_correct=[]; ITPS_error=[];
                correct_ITPS_nosub_baseCorr=[]; error_ITPS_nosub_baseCorr=[];
                correct_timefreqs=[]; error_timefreqs=[];
                correct_phase=[]; error_phase=[];
                
                for fi=1:length(frex) % loop through all frequencies
                    
                    correct_temppow = zeros(nData, EEGa.trials);
                    error_temppow = zeros(nData, EEGb.trials);
                    
                    %% Create wavelate
                    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
                    waveletX = fft(wavelet, nConv); % fft of wavelet
                    waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
                    
                    %% Loop through all trials
                    
                    for trl=1:EEGa.trials
                        if Cond1Num > cutoff
                            
                            correct_data = fft(squeeze(EEGa.data(ch,:,trl)), nConv);
                            
                            % run convolution
                            correct_conv = ifft(waveletX .* correct_data);
                            correct_conv = correct_conv(half_wave+1:end-half_wave);
                            
                            % compute power
                            correct_temppow (:,trl) = abs(correct_conv).^2;
                            
                            %save convolution
                            correct_conv_trials(:,trl) = correct_conv;
                            
                        end
                    end
                    
                    % compute ITPS without subsampling
                    ITPS_correct(fi,:) = abs(mean( exp(1i*angle(correct_conv_trials)),2));
                    %Baseline correct ITPS for this frequency, condition,channel and participant
                    correct_ITPS_nosub_baseCorr(fi,:) = ( ITPS_correct(fi,:) - mean(fi,ITPS_correct(basetimeidx(1):basetimeidx(end))));
                    
                    
                    
                    for trl=1:EEGb.trials
                        if Cond2Num > cutoff
                            
                            error_data = fft(squeeze(EEGb.data(ch,:,trl)), nConv);
                            
                            %run convolution
                            error_conv = ifft(waveletX .* error_data);
                            error_conv = error_conv(half_wave+1:end-half_wave);
                            
                            %compute power
                            error_temppow (:,trl) = abs(error_conv).^2;
                            
                            %save convolution
                            error_conv_trials(:,trl) = error_conv;
                        end
                    end
                    
                    % compute ITPS without subsampling
                    ITPS_error(fi,:) = abs(mean( exp(1i*angle(error_conv_trials)),2));
                    %Baseline correct ITPS for this frequency, condition,channel and participant
                    error_ITPS_nosub_baseCorr(fi,:) = ( ITPS_error(fi,:) - mean(fi,ITPS_error(basetimeidx(1):basetimeidx(end))));
                    
                    
                    % Power average all trials
                    correct_avg_temppow = squeeze(mean(correct_temppow, 2));
                    error_avg_temppow = squeeze(mean(error_temppow, 2));
                    
                    
                    % Baseline normalization within condition at each frequency
                    correct_timefreqs(fi,:) = 10*log10( correct_avg_temppow ./ mean(correct_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    error_timefreqs(fi,:) = 10*log10( error_avg_temppow ./ mean(error_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    
                    
                    %ICPS
                    correct_phase(fi,:,:) = correct_conv_trials;
                    error_phase(fi,:,:) = error_conv_trials;
                    
                    
                end
                
                %build dataset of baseline corrected TF with each channel
                correct_tf_data(:,:,ch) = correct_timefreqs;
                error_tf_data(:,:,ch) = error_timefreqs;
                
                %build dataset of baseline corrected ITPS with each channel
                correct_ITPS_nosubsamp_baseCorr(:,:,ch) = correct_ITPS_nosub_baseCorr;
                error_ITPS_nosubsamp_baseCorr(:,:,ch) = error_ITPS_nosub_baseCorr;
                
                
                %build dataset of phase info to be used to calculate ITPS and ICPS
                % chan x freq x time x trials
                correct_phase_data(ch,:,:,:) = correct_phase;
                error_phase_data(ch,:,:,:) = error_phase;
                
            end
            
            %% compute ITPS and ICPS
            fprintf('\n\n\n*** Doing ITPS and ICPS for subject %s ***\n\n\n', subject);
            
            %establish number of subsamples
            subsample_num = 100;
            
            %initialize temp connectivity
            Correct_phase_connectivity_temp = zeros(size(correct_phase_data,2), size(correct_phase_data,3));
            Error_phase_connectivity_temp = zeros(size(error_phase_data,2), size(error_phase_data,3));
            
            %initialize temp ITPS
            correct_ITPS_temp = zeros(size(correct_phase_data,2), size(correct_phase_data,3));
            error_ITPS_temp = zeros(size(error_phase_data,2), size(error_phase_data,3));
            
            % initializing the subsample structures. 
            correct_ICPS_subsamples = []; error_ICPS_subsamples = [];
            correct_ITPS_subsamples = []; error_ITPS_subsamples = [];
            
            correct_ICPS_subsamples_ch = []; error_ICPS_subsamples_ch = [];
            correct_ITPS_subsamples_ch = []; error_ITPS_subsamples_ch = [];
            
            % now compute connectivity from FCz seed
            for chanj=1:size(correct_phase_data, 1)
                
                % take cross-spectral density between two channels - one being the seed for ICPS
                seed = 4; % change this electrode number depending what seed you want to use - we chose 4, which is FCz
                correct_crossspecden = squeeze(correct_phase_data(seed,:,:,:) .* conj(correct_phase_data(chanj,:,:,:)));
                error_crossspecden = squeeze(error_phase_data(seed,:,:,:) .* conj(error_phase_data(chanj,:,:,:)));
                
                % Subsample trials to equate trial counts between conditions and across participants 
                for samp=1:subsample_num
                    %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                    
                    %Get indices of trials for this subsample
                    correct_subtrials = randsample(TrialNums(sub).Correct,cutoff,false);
                    error_subtrials = randsample(TrialNums(sub).Error,cutoff,false);
                    %Index into those trials and pull them out for ICPS analyses
                    correct_crossspecden_temp = correct_crossspecden(:,:,correct_subtrials);
                    error_crossspecden_temp = error_crossspecden(:,:,error_subtrials);
                    %Index into those trials and pull them out for ITPS analyses
                    correct_phase_data_temp = squeeze(correct_phase_data(chanj,:,:,correct_subtrials));
                    error_phase_data_temp = squeeze(error_phase_data(chanj,:,:,error_subtrials));
                    
                    for freq=1:size(correct_phase_data, 2)
                        fprintf('\n\n\n*** Now in frequency %d for sample %d for channel %d for subject %s ***\n\n\n', freq, samp, chanj, subject);
                        
                        % phase angle difference for ICPS (shortcut, as implemented in Cohen's 2014 book)
                        Correct_phase_connectivity_temp(freq,:) = abs(mean(exp(1i*angle(correct_crossspecden_temp(freq,:,:))),3));
                        Error_phase_connectivity_temp(freq,:) = abs(mean(exp(1i*angle(error_crossspecden_temp(freq,:,:))),3));
                        
                        % compute ITPS
                        correct_ITPS_temp(freq,:) = abs(mean( exp(1i*angle(correct_phase_data_temp(freq,:,:))),3));
                        error_ITPS_temp(freq,:) = abs(mean( exp(1i*angle(error_phase_data_temp(freq,:,:))),3));
                        
                    end
                    
                    %create matrix of subsamples for ICPS - samp x freq x time
                    correct_ICPS_subsamples(samp,:,:) = Correct_phase_connectivity_temp;
                    error_ICPS_subsamples(samp,:,:) = Error_phase_connectivity_temp;
                    
                    %create matrix of subsamples for ITPS - samp x chan x freq x time
                    correct_ITPS_subsamples(samp,:,:) = correct_ITPS_temp;
                    error_ITPS_subsamples(samp,:,:) = error_ITPS_temp;
                    
                    
                end %end loop through subsamples
                
                correct_ICPS_subsamples_ch(chanj,:,:,:) = correct_ICPS_subsamples;
                error_ICPS_subsamples_ch(chanj,:,:,:) = error_ICPS_subsamples;
                
                correct_ITPS_subsamples_ch(chanj,:,:,:) = correct_ITPS_subsamples;
                error_ITPS_subsamples_ch(chanj,:,:,:) = error_ITPS_subsamples;
                
                
            end%end loop through channels
            
            
            %Average across subsamples (second dimension) for ICPS
            correct_ICPS_data = squeeze(mean(correct_ICPS_subsamples_ch,2));
            error_ICPS_data = squeeze(mean(error_ICPS_subsamples_ch,2));
            
            %Average across subsamples (second dimension) for ITPS
            correct_ITPS_subsamp_nobaseCorr = squeeze(mean(correct_ITPS_subsamples_ch,2));
            error_ITPS_subsamp_nobaseCorr = squeeze(mean(error_ITPS_subsamples_ch,2));
            
            %Baseline Correct subsampled ITPS
            for chanj=1:size(correct_ITPS_subsamp_nobaseCorr,1)
                for fi = 1:size(correct_ITPS_subsamp_nobaseCorr,2)
                    correct_ITPS_subsamp_baseCorr(chanj,fi,:) = ( correct_ITPS_subsamp_nobaseCorr(chanj,fi,:) - mean(correct_ITPS_subsamp_nobaseCorr(chanj,fi,basetimeidx(1):basetimeidx(end))));
                    error_ITPS_subsamp_baseCorr(chanj,fi,:) = ( error_ITPS_subsamp_nobaseCorr(chanj,fi,:) - mean(error_ITPS_subsamp_nobaseCorr(chanj,fi,basetimeidx(1):basetimeidx(end))));
                end
            end
            
            %Baseline Correct subsampled ICPS
            for chanj=1:size(correct_ICPS_data,1)
                for fi = 1:size(correct_ICPS_data,2)
                    correct_icps_subsamp_baseCorr(chanj,fi,:) = ( correct_ICPS_data(chanj,fi,:) - mean(correct_ICPS_data(chanj,fi,basetimeidx(1):basetimeidx(end))));
                    error_icps_subsamp_baseCorr(chanj,fi,:) = ( error_ICPS_data(chanj,fi,:) - mean(error_ICPS_data(chanj,fi,basetimeidx(1):basetimeidx(end))));
                end
            end
            
            
            %% Save data
            time = EEG.times;
            frequency = frex;
            
            %% Downsampling
            times2save = 1:5:750;
            
            % downsample TF
            correct_tf_data_ds = correct_tf_data(:,times2save,:);
            error_tf_data_ds = error_tf_data(:,times2save,:);
            
            % downsample subsampled baseline corrected ITPS - chan x freq x time
            correct_ITPS_subsamp_baseCorr_ds = correct_ITPS_subsamp_baseCorr(:,:,times2save);
            error_ITPS_subsamp_baseCorr_ds = error_ITPS_subsamp_baseCorr(:,:,times2save);
            
            % downsample nbaseline corrected ICPS - chan x freq x time
            correct_icps_subsamp_baseCorr_ds = correct_icps_subsamp_baseCorr(:,:,times2save);
            error_icps_subsamp_baseCorr_ds = error_icps_subsamp_baseCorr(:,:,times2save);
            
            save_data=[subject, dataset_name];
            save (save_data,'time', 'times2save', 'frequency', 'correct_tf_data_ds', 'error_tf_data_ds','correct_ITPS_subsamp_baseCorr_ds','error_ITPS_subsamp_baseCorr_ds','channel_location');
            
            ICPS_dataset_name = '_Zoo_Resp_ICPS_100sub';
            save_ICPS_data = [subject, ICPS_dataset_name];
            save (save_ICPS_data, 'time','times2save', 'frequency', 'correct_icps_subsamp_baseCorr_ds','error_icps_subsamp_baseCorr_ds', 'channel_location');
            
            
        end %end loop through if they have enough trials
        
        EEG=[];
    end %end check if it already exists
end %end loop through subjects


param.date = datestr(now,'mm_dd_yyyy');
% Saving the number of trials and the people included
table_events = struct2table(TrialNums);
writetable(table_events, ['/export/data/cdl/Projects/ECHO/Descriptives/Trial_Numbers_TF_ITPS_Zoo_Resp_' param.date '.csv']);
