%% Perform Reliability of ITPS and ICPS for all trials
% This script performs the reliability of ITPS and ICPS analyses for all trials,  
% rather than by increasing numbers of trials, showed in Morales et al. (under review).
% This script will save a separate file per participant.
% The files produced by this script can then be read into R using the ECHO_ERN_ERP_TF_ITPS_ICPS_R1.Rmd 

clear all;
clc;

%% Adding this path for Ludwig
addpath('/data/cdl/Toolbox/eeglab13_6_5b/')

%% Locations
Data_Location = '/data/cdl/Projects/ECHO/CSD_tranformed_data/Zoo_3s/';
Save_Data = '/data/rosalind/Dropboxes/moraless/ECHO/Reliability/R1_subs_300_100_All/';

dataset_name = '_TimeFreqs_ITPS_ICPS'; % Dataset name to dave
%% Subject list
cd(Data_Location)
subnum=dir('*Zoo_Epoched_Matched_CSD.set'); % Looking for this data
sub_list={subnum.name};
for i =1:length(sub_list)
    sub = sub_list{i};
    subject_list{i}= sub(1:15);
end

%%
% Reading in files that are included in the paper
BehData = importdata('/data/cdl/Projects/ECHO/Descriptives/ECHO_ERN_IDS_Paper_01_06_2021.csv');
subject_list = intersect(subject_list, BehData); % Deleting participants on the list

conditions = 2; % This is just for Go vs. No-Go
cutoff = 5; % Maybe the actual names to pass them for naming and stuff later


%% Events
event_markers = {'resp'}; % 'Stm+' or 'resp'

%% Change current directory
cd (Save_Data)

%% define baseline time window
baseline = [-300 -100]; % Changed this based on reviewer comment

%% frequency parameters
% min_freq =  3;
min_freq =  1;
max_freq = 30;
num_frex = 60;

% frequencies vector
frex = logspace(log10(min_freq), log10(max_freq), num_frex);

%% wavelet cycles - variable : min 4 max 10
range_cycles = [ 3 10 ];
% cycles vector
cylvec = logspace(log10(range_cycles(1)),log10(range_cycles(end)), num_frex)./ (2*pi*frex);

%% Running EEGLAB
eeglab

%% Loop through all subjects
for sub=200:length(subject_list)
    
    % load data
    subject = subject_list{sub};
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', subject);
    
    % Check if file exists - this is helpful for running multiple instances of matlab in parallel
    if exist([subject '_TimeFreqs_ITPS_ICPS_Reliability.mat'], 'file') == 2
        % File exists.  skip....
        fprintf('\n\n\n*** Skipping subject %d (%s) ***\n\n\n', sub, subject);
    else
        TrialNums(sub).Subject = subject_list{sub};
        
        
        EEG=pop_loadset('filename',[subject '_Zoo_Epoched_Matched_CSD.set'], 'filepath', Data_Location);
        EEG = eeg_checkset( EEG );
        
        %% Change sampling rate to 50Hz
        EEG = eeg_checkset( EEG );
        EEG = pop_resample( EEG, 50);
        
        %% Initialize variables for separate conditions
        EEGa = [];
        EEGb = [];
        EEGall = [];
        %% Defining the conditions
        % Hard coding for now:
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
        
        if Cond1Num > cutoff && Cond2Num > cutoff 
            % Stimulus Novel
            EEG_all = pop_selectevent( EEG, 'StimType', 'Resp', 'Bad', 0, 'deleteevents','on','deleteepochs','on','invertepochs','off');
            EEG_all = eeg_checkset( EEG_all );
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
        basetimeidx   = dsearchn(EEG.times', baseline'); % baseline indecies
        
        %% initialize output time-frequency data
        correct_tf_data = zeros( length(frex), nData, EEG.nbchan);
        error_tf_data = zeros( length(frex), nData, EEG.nbchan);
        
        %% initialize output ITPS data - freq x Data x channels
        correct_ITPS_nosubsamp_baseCorr = zeros( length(frex), nData, EEG.nbchan);
        error_ITPS_nosubsamp_baseCorr = zeros( length(frex), nData, EEG.nbchan);
        
        %blank out phase data
        correct_phase_data=[];
        error_phase_data=[];
        
        %% Check if there are enough trials before convolution
        if TrialNums(sub).Correct > cutoff && TrialNums(sub).Error > cutoff
            TrialNums(sub).IncludedERP = 1;
            
            %% Run wavelet convolution
            for ch=1:EEG.nbchan % Loop through all channels
                if ismember(ch, [4 7 54 12 13 59 60]) % to seepd things up just doing it for chans of interest
                fprintf('\n*** Processing chan %d for Subject %s ***', ch, subject);
                
                timefreqs = zeros(length(frex), nData);
                
                correct_conv_trials=[]; error_conv_trials=[]; all_conv_trials=[];
                ITPS_correct=[]; ITPS_error=[];
                correct_ITPS_nosub_baseCorr=[]; error_ITPS_nosub_baseCorr=[];
                correct_ITPS_nosubsamp_baseCorr=[]; error_ITPS_nosubsamp_baseCorr=[];
                correct_timefreqs=[]; error_timefreqs=[];
                correct_phase=[]; error_phase=[];
                correct_timefreqs_trials = []; error_timefreqs_trials=[];
                
                for fi=1:length(frex) % loop through all frequencies
                    
                    correct_temppow = zeros(nData, EEGa.trials);
                    error_temppow = zeros(nData, EEGb.trials);
                    %             all_temppow = zeros(nData, EEGe.trials);
                    
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
                    correct_ITPS_nosub_baseCorr(fi,:) = ( ITPS_correct(fi,:) - mean(ITPS_correct(fi, basetimeidx(1):basetimeidx(end))));
                    
                    
                    
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
                    error_ITPS_nosub_baseCorr(fi,:) = ( ITPS_error(fi,:) - mean(ITPS_error(fi, basetimeidx(1):basetimeidx(end))));
                    
                    
                    for trl=1:EEG_all.trials
                        if Cond1Num > cutoff && Cond2Num > cutoff 
                            
                            all_data = fft(squeeze(EEG_all.data(ch,:,trl)), nConv);
                            
                            all_conv = ifft(waveletX .* all_data);
                            all_conv = all_conv(half_wave+1:end-half_wave);
                            
                            all_temppow (:,trl) = abs(all_conv).^2;
                            all_conv_trials(:,trl) = all_conv;
                            
                        end
                    end
                    
                    % Power average all trials
                    all_avg_temppow = squeeze(mean(all_temppow, 2));
                    
                    % Power average across trials within condition
                    correct_avg_temppow = squeeze(mean(correct_temppow, 2));
                    error_avg_temppow = squeeze(mean(error_temppow, 2));
                    
                    
                    % Baseline normalization within condition at each frequency
                    correct_timefreqs(fi,:) = 10*log10( correct_avg_temppow ./ mean(correct_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    error_timefreqs(fi,:) = 10*log10( error_avg_temppow ./ mean(error_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    
                    
                    % %Baseline normalization across condition
                    % correct_timefreqs(fi,:) = 10*log10( correct_avg_temppow ./ mean(all_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    % error_timefreqs(fi,:) = 10*log10( error_avg_temppow ./ mean(all_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    
                    %% Saving trial by trial data
                    for trl=1:size(correct_temppow,2)
                        correct_timefreqs_trials(fi,:,trl) = 10*log10( correct_temppow(:,trl) ./ mean(correct_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    end
                    for trl=1:size(error_temppow,2)
                        error_timefreqs_trials(fi,:,trl) = 10*log10( error_temppow(:,trl) ./ mean(error_avg_temppow(basetimeidx(1):basetimeidx(end))));
                    end
                    
                     
                    %ICPS
                    correct_phase(fi,:,:) = correct_conv_trials;
                    error_phase(fi,:,:) = error_conv_trials;
                    
                    
                end
                
                %build dataset of baseline corrected TF with each channel
                correct_tf_data(:,:,ch) = correct_timefreqs;
                error_tf_data(:,:,ch) = error_timefreqs;
                
                %% Trial-level data
                correct_tf_trials_data(:,:,:,ch) = correct_timefreqs_trials;
                error_tf_trials_data(:,:,:,ch) = error_timefreqs_trials;
                
                
                %             %build dataset of baseline corrected ITPS with each channel
                %             correct_ITPS_nosubsamp_baseCorr(:,:,ch) = correct_ITPS_nosub_baseCorr;
                %             error_ITPS_nosubsamp_baseCorr(:,:,ch) = error_ITPS_nosub_baseCorr;
                
                
                %build dataset of phase info to be used to calculate ITPS and ICPS
                % chan x freq x time x trials
                correct_phase_data(ch,:,:,:) = correct_phase;
                error_phase_data(ch,:,:,:) = error_phase;
                end
            end
            
            %% compute ITPS and ICPS with subsampling
            fprintf('\n\n\n*** Doing ITPS and ICPS for subject %s ***\n\n\n', subject);
            
            % Deleting all the channels that were not selected - the rest
            % should be zero
            correct_phase_data = correct_phase_data([4 7 12 13 54 59 60],:,:,:);
            error_phase_data = error_phase_data([4 7 12 13 54 59 60],:,:,:);
            
            %establish number of subsamples
            subsample_num = 100; % This was 2 at first
            
            correct_phase_connectivity_temp_1 = []; error_phase_connectivity_temp_1 = [];
            correct_phase_connectivity_temp_2 = []; error_phase_connectivity_temp_2 = [];
            correct_ITPS_temp_1 = []; error_ITPS_temp_1 = [];
            correct_ITPS_temp_2 = []; error_ITPS_temp_2 = [];
            
            correct_ICPS_subsamples = []; error_ICPS_subsamples = [];
            correct_ITPS_subsamples = []; error_ITPS_subsamples = [];
            
            correct_ICPS_subsamples_ch = []; error_ICPS_subsamples_ch = [];
            correct_ITPS_subsamples_ch = []; error_ITPS_subsamples_ch = [];
            


            % now compute connectivity from FCz seed
            
            for chanj=1:size(correct_phase_data, 1)
                
                % take cross-spectral density between two channels - one being the seed for ICPS
                seed = 1; % change this electrode number depending what seed you want to use - NOTE this is based on the subsample of channels
                correct_crossspecden = squeeze(correct_phase_data(seed,:,:,:) .* conj(correct_phase_data(chanj,:,:,:)));
                error_crossspecden = squeeze(error_phase_data(seed,:,:,:) .* conj(error_phase_data(chanj,:,:,:)));
                
                %Here, we will subsample trials to equate trial counts between
                %conditions and across participants - you can choose what
                %number of subsamples you would like to do
%                 subsample_array = [4 8 12 16 20 24 28 32]; % Create an array of subsamples to iterate through
%                 for subsample_i=1:length(subsample_array)
%                     subsample_n = subsample_array(subsample_i);
                    % Doing it for each condition separately to get different numbers of trials per condition
                    %% Correct
%                     if Cond1Num >= subsample_n 
                        for samp=1:subsample_num
                            %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                            %Get indices of trials for this subsample
                            correct_subtrials = randsample(TrialNums(sub).Correct,round(TrialNums(sub).Correct/2,0),false);
                                
                                % Divide them into half
                                correct_trials = 1:TrialNums(sub).Correct;
                                correct_subtrials_1 = correct_trials(ismember(correct_trials, correct_subtrials));
                                correct_subtrials_2 = correct_trials(~ismember(correct_trials, correct_subtrials));
                            
                            %Index into those trials and pull them out for ICPS analyses
                            correct_crossspecden_temp_1 = correct_crossspecden(:,:,correct_subtrials_1);
                            correct_crossspecden_temp_2 = correct_crossspecden(:,:,correct_subtrials_2);                            
                            
                            %Index into those trials and pull them out for ITPS analyses
                            correct_phase_data_temp_1 = squeeze(correct_phase_data(chanj,:,:,correct_subtrials_1));
                            correct_phase_data_temp_2 = squeeze(correct_phase_data(chanj,:,:,correct_subtrials_2));                   
                            
                            
                            for freq=1:size(correct_phase_data, 2)
                                %                             fprintf('\n\n\n*** Now in frequency %d for sample %d for channel %d for subject %s ***\n\n\n', freq, samp, chanj, subject);
                                
                                % phase angle difference for ICPS (shortcut, as implemented in Cohen's book)
                                correct_phase_connectivity_temp_1(freq,:) = abs(mean(exp(1i*angle(correct_crossspecden_temp_1(freq,:,:))),3));
                                correct_phase_connectivity_temp_2(freq,:) = abs(mean(exp(1i*angle(correct_crossspecden_temp_2(freq,:,:))),3));
                                
                                % compute ITPS
                                correct_ITPS_temp_1(freq,:) = abs(mean( exp(1i*angle(correct_phase_data_temp_1(freq,:,:))),3));
                                correct_ITPS_temp_2(freq,:) = abs(mean( exp(1i*angle(correct_phase_data_temp_2(freq,:,:))),3));
                            end
                            
                            
                            %create matrix of subsamples for ICPS - n_subsamples x half x samp x freq x time
                            correct_ICPS_subsamples(1,samp,:,:) = correct_phase_connectivity_temp_1;
                            correct_ICPS_subsamples(2,samp,:,:) = correct_phase_connectivity_temp_2;
                            
                            %create matrix of subsamples for ITPS - half x samp x freq x time
                            correct_ITPS_subsamples(1,samp,:,:) = correct_ITPS_temp_1;
                            correct_ITPS_subsamples(2,samp,:,:) = correct_ITPS_temp_2;
                        end
%                     end
                    
                    %% Error 
%                     if Cond2Num >= subsample_n 
                        for samp=1:subsample_num
                            %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                            %Get indices of trials for this subsample
                            error_subtrials = randsample(TrialNums(sub).Error,round(TrialNums(sub).Error/2,0),false);
                                
                                % Divide them into half
                                error_trials = 1:TrialNums(sub).Error;
                                error_subtrials_1 = error_trials(ismember(error_trials, error_subtrials));
                                error_subtrials_2 = error_trials(~ismember(error_trials, error_subtrials));
                            
                            %Index into those trials and pull them out for ICPS analyses
                            error_crossspecden_temp_1 = error_crossspecden(:,:,error_subtrials_1);
                            error_crossspecden_temp_2 = error_crossspecden(:,:,error_subtrials_2);
                            
                            
                            %Index into those trials and pull them out for ITPS analyses
                            error_phase_data_temp_1 = squeeze(error_phase_data(chanj,:,:,error_subtrials_1));
                            error_phase_data_temp_2 = squeeze(error_phase_data(chanj,:,:,error_subtrials_2));
                            
                            for freq=1:size(correct_phase_data, 2)
                                % phase angle difference for ICPS (shortcut, as implemented in Cohen's book)
                                error_phase_connectivity_temp_1(freq,:) = abs(mean(exp(1i*angle(error_crossspecden_temp_1(freq,:,:))),3));
                                error_phase_connectivity_temp_2(freq,:) = abs(mean(exp(1i*angle(error_crossspecden_temp_2(freq,:,:))),3));
                                
                                % compute ITPS
                                error_ITPS_temp_1(freq,:) = abs(mean( exp(1i*angle(error_phase_data_temp_1(freq,:,:))),3));
                                error_ITPS_temp_2(freq,:) = abs(mean( exp(1i*angle(error_phase_data_temp_2(freq,:,:))),3));
                            end
                            
                            %create matrix of subsamples for ICPS - n_subsamples x half x samp x freq x time
                            error_ICPS_subsamples(1,samp,:,:) = error_phase_connectivity_temp_1;
                            error_ICPS_subsamples(2,samp,:,:) = error_phase_connectivity_temp_2;
                            
                            %create matrix of subsamples for ITPS - n_subsamples x half x samp x freq x time
                            error_ITPS_subsamples(1,samp,:,:) = error_ITPS_temp_1;
                            error_ITPS_subsamples(2,samp,:,:) = error_ITPS_temp_2;
                        end
%                     end
%                 end %end loop through subsamples
                
                correct_ICPS_subsamples_ch(chanj,:,:,:,:) = correct_ICPS_subsamples;
                error_ICPS_subsamples_ch(chanj,:,:,:,:) = error_ICPS_subsamples;
                
                correct_ITPS_subsamples_ch(chanj,:,:,:,:) = correct_ITPS_subsamples;
                error_ITPS_subsamples_ch(chanj,:,:,:,:) = error_ITPS_subsamples;
                
            end %end loop through channels
            
            correct_ICPS_data = correct_ICPS_subsamples_ch;
            error_ICPS_data = error_ICPS_subsamples_ch;
            
            correct_ITPS_subsamp_nobaseCorr = correct_ITPS_subsamples_ch;
            error_ITPS_subsamp_nobaseCorr = error_ITPS_subsamples_ch;
            
            
            correct_ITPS_subsamp_baseCorr = []; correct_icps_subsamp_baseCorr = []; 
            %Baseline Correct subsampled ITPS
            for chanj=1:size(correct_ITPS_subsamp_nobaseCorr,1)
                for half_i=1:2
                    for samp_i=1:size(correct_ITPS_subsamp_nobaseCorr,3)
                        for fi = 1:size(correct_ITPS_subsamp_nobaseCorr,4)
                            correct_ITPS_subsamp_baseCorr(chanj,half_i,samp_i,fi,:) = ( correct_ITPS_subsamp_nobaseCorr(chanj,half_i,samp_i,fi,:) - mean(correct_ITPS_subsamp_nobaseCorr(chanj,half_i,samp_i,fi,basetimeidx(1):basetimeidx(end))));
                            correct_icps_subsamp_baseCorr(chanj,half_i,samp_i,fi,:) = ( correct_ICPS_data(chanj,half_i,samp_i,fi,:) - mean(correct_ICPS_data(chanj,half_i,samp_i,fi,basetimeidx(1):basetimeidx(end))));
                        end
                    end
                end
            end
            
            error_ITPS_subsamp_baseCorr = []; error_icps_subsamp_baseCorr = [];
            %Baseline Correct error subsampled ITPS & ICPS
            for chanj=1:size(error_ICPS_data,1)
                for half_i=1:2
                    for samp_i=1:size(correct_ITPS_subsamp_nobaseCorr,3)
                        for fi = 1:size(correct_ITPS_subsamp_nobaseCorr,4)
                            error_ITPS_subsamp_baseCorr(chanj,half_i,samp_i,fi,:) = ( error_ITPS_subsamp_nobaseCorr(chanj,half_i,samp_i,fi,:) - mean(error_ITPS_subsamp_nobaseCorr(chanj,half_i,samp_i,fi,basetimeidx(1):basetimeidx(end))));
                            error_icps_subsamp_baseCorr(chanj,half_i,samp_i,fi,:) = ( error_ICPS_data(chanj,half_i,samp_i,fi,:) - mean(error_ICPS_data(chanj,half_i,samp_i,fi,basetimeidx(1):basetimeidx(end))));
                        end
                    end
                end
            end
            %% Save Dataset
            %% Downsampling
            times2save = 1:1:150;
            
            % TF
            correct_tf_data = correct_tf_data(:,times2save,:);
            error_tf_data = error_tf_data(:,times2save,:);
            
            correct_tf_trials_data = correct_tf_trials_data(:,times2save,:,:);
            error_tf_trials_data = error_tf_trials_data(:,times2save,:,:);
            
            % ITSP - chan x trials x helf x subs x freq x time
            correct_ITPS_data = correct_ITPS_subsamp_baseCorr(:,:,:,:,times2save);
            error_ITPS_data = error_ITPS_subsamp_baseCorr(:,:,:,:,times2save);
            
            % downsample nbaseline corrected ICPS - chan x trials x helf x subs x freq x time
            correct_icps_data = correct_icps_subsamp_baseCorr(:,:,:,:,times2save);
            error_icps_data = error_icps_subsamp_baseCorr(:,:,:,:,times2save);
            
            %% Save data
            time = EEG.times(times2save);
            frequency = frex;
            
            %% Commenting this out because this is now done in a separate script rather than using trial-by-trial data
%             save_data=[subject, dataset_name];
%             save (save_data, 'time', 'frequency', 'correct_tf_data', 'error_tf_data', 'channel_location',...
%                 'correct_ITPS_data', 'error_ITPS_data', 'correct_icps_data', 'error_icps_data');
            
            
            %% For trial level data, I can select the frequencies here to save space
            % Frequency to extract
            % Theta
            freq_windows = [ 4 8];
            
            % Fine frequency indices
            for i=1:size(freq_windows,1)
                for j=1:2
                    [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
                end
            end
            
            correct_tf_trials_data_4_8 = squeeze(mean(correct_tf_trials_data(Freqs_Idx(1):Freqs_Idx(2),:,:,:),1));
            error_tf_trials_data_4_8 = squeeze(mean(error_tf_trials_data(Freqs_Idx(1):Freqs_Idx(2),:,:,:),1));
            
            % This is how we were doing it before: e.g., correct_ITPS_data_4_8 = squeeze(mean(correct_ITPS_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4));
            dims = [1 2 3 5];
            sz = size(correct_ITPS_data);
            correct_ITPS_data_4_8 = reshape(mean(correct_ITPS_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            sz = size(error_ITPS_data);
            error_ITPS_data_4_8 = reshape(mean(error_ITPS_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            
            sz = size(correct_icps_data);
            correct_icps_data_4_8 = reshape(mean(correct_icps_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            sz = size(error_icps_data);
            error_icps_data_4_8 = reshape(mean(error_icps_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            
            
            % Delta
            freq_windows = [1 4];
            % Fine frequency indices
            for i=1:size(freq_windows,1)
                for j=1:2
                    [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
                end
            end
            
            correct_tf_trials_data_1_4 = squeeze(mean(correct_tf_trials_data(Freqs_Idx(1):Freqs_Idx(2),:,:,:),1));
            error_tf_trials_data_1_4 = squeeze(mean(error_tf_trials_data(Freqs_Idx(1):Freqs_Idx(2),:,:,:),1));
            
            % This is how we were doing it before: e.g., correct_ITPS_data_4_8 = squeeze(mean(correct_ITPS_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4));
            dims = [1 2 3 5];
            sz = size(correct_ITPS_data);
            correct_ITPS_data_1_4 = reshape(mean(correct_ITPS_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            sz = size(error_ITPS_data);
            error_ITPS_data_1_4 = reshape(mean(error_ITPS_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            
            sz = size(correct_icps_data);
            correct_icps_data_1_4 = reshape(mean(correct_icps_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            sz = size(error_icps_data);
            error_icps_data_1_4 = reshape(mean(error_icps_data(:,:,:,Freqs_Idx(1):Freqs_Idx(2),:),4),sz(dims));
            
%             save_data=[subject, dataset_name, '_trial_level'];
%             save (save_data, 'time', 'frequency', 'correct_tf_trials_data_4_8', 'error_tf_trials_data_4_8', ...
%                 'correct_tf_trials_data_1_4', 'error_tf_trials_data_1_4', 'channel_location');

            %% Select time of interest to save space
            time_windows = [0 300];
            % Fine frequency indices
            for i=1:size(time_windows,1)
                for j=1:2
                    [~,Times_Idx(i,j)] = min(abs(time-time_windows(i,j)));
                end
            end
            
            dims = [1 2 3];
            sz = size(correct_ITPS_data_4_8);
            correct_ITPS_data_4_8 = reshape(mean(correct_ITPS_data_4_8(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            correct_ITPS_data_1_4 = reshape(mean(correct_ITPS_data_1_4(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            sz = size(error_ITPS_data_4_8);
            error_ITPS_data_4_8 = reshape(mean(error_ITPS_data_4_8(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            error_ITPS_data_1_4 = reshape(mean(error_ITPS_data_1_4(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            
            sz = size(correct_icps_data_4_8);
            correct_icps_data_4_8 = reshape(mean(correct_icps_data_4_8(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            correct_icps_data_1_4 = reshape(mean(correct_icps_data_1_4(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            sz = size(error_icps_data_4_8);
            error_icps_data_4_8 = reshape(mean(error_icps_data_4_8(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            error_icps_data_1_4 = reshape(mean(error_icps_data_1_4(:,:,:,Times_Idx(1):Times_Idx(2)),4),sz(dims));
            
            
%%            
            save_data=[subject, dataset_name, '_Reliability'];
            save (save_data, 'time', 'frequency', 'correct_ITPS_data_4_8', 'error_ITPS_data_4_8', ...
                'correct_ITPS_data_1_4', 'error_ITPS_data_1_4',...
                'correct_icps_data_4_8', 'error_icps_data_4_8', ...
                'correct_icps_data_1_4', 'error_icps_data_1_4',  'channel_location');
            
        end
        
        EEG=[];
        correct_tf_data = [];
        error_tf_data = [];
        correct_tf_trials_data = []; error_tf_trials_data = [];
        correct_tf_trials_data_4_8 = []; error_tf_trials_data_4_8 = [];
        correct_tf_trials_data_1_4 = []; error_tf_trials_data_1_4 = [];
        correct_ITPS_data = []; error_ITPS_data = [];
        correct_ITPS_data_4_8 = []; error_ITPS_data_4_8 = [];
        correct_ITPS_data_1_4 = []; error_ITPS_data_1_4 = [];
        correct_icps_data = []; error_icps_data = [];
        correct_icps_data_4_8 = []; error_icps_data_4_8 = [];
        correct_icps_data_1_4 = []; error_icps_data_1_4 = [];
        correct_phase_data = []; error_phase_data = [];

        
        
    end
end


