%Compute TF, ITPS, ICPS, and wPLI measures for EEG data
%Maureen Bowers 6/29/2021 based on scripts by Ranjan Debnath

clear all;
clc;

%%
%%%%% Setting paths %%%%%
%1. Data Location
data_location = '...';

%2. Save Data Location
save_location = '...';

%3. Scripts Location
scripts_location = '...';

%4. Set EEGLab path
addpath(genpath('...'));

%%
%%%%% Information about dataset and analysis procedures %%%%%%
%5. Number of channels
nbchan = 128;

%6. What kind of data? Resting State or Event-Related Data
RestorEvent = 0; %1 = rest, 0 = event

%7. What are your conditions of interest if using Event-Related Data? This
%naming convention should come from the Edit_events.m script provided.
%Note: not needed for resting state data.
Conds = {'Stim_0_0','Stim_0_1'};

%8. Minimum number of trials to analyze
mintrialnum =8; %If the participant does not have enough trials in a condition based on this cutoff, a "notenoughdata.mat" file will be saved into save_location.

%9. Would you like to baseline correct your data? NOTE: ICPS calculated over time will not be baseline corrected.
BaselineCorrect = 1; %1=Yes, 0=No

%10. What time period would you like to use to baseline correct
BaselineTime = [-400 -100];
%Put in time in ms for event-related data. For example, if you want to baseline correct from -100 to 0ms before the 
% event for evernt-related paradigms, put [-100 0]; If you would like to correct to the whole epoch for resting state data,
%put in the beginning and end time of the entire epoch (e.g., -500 to 1500 epoch, put [-500 1500])

%11. Would you like to downsample the output to 125Hz? This is done after the time-frequency computations and will have minor impacts on the resolution. 
%We recommend downsampling to reduce file size for ease of storage. All data will be initially downsampled to 250Hz. 
Downsample = 1; %0=No, 1=Yes

%12. Dataset Name - will be appended to the saved files
DatasetName = '_StudyName';

%%
%%%%% Settings for Time-Frequency Measures %%%%%

%13. Minimum and Maximum Frequency, Number of Frequency Bins, and range cycles to calculate complex Morlet wavelet decomposition
min_freq = 1;
max_freq = 30;
num_frex = 40; % number of frequency bins between minimum and maximum frequency
range_cycles = [3 10]; % wavelet cycles: min 3 max 10


%%
%%%%% Questions about phase-based measures %%%%%

%14. Would you like to calculate inter-trial phase synchrony (ITPS) in addition to TF?
ITPS_calc = 1; %(1=yes,0=no)

%15. Would you like to subsample trials? This is recommended for event-related paradigms, 
% especially when there are uneven numbers of trials in conditions. 
Subsample = 1; %1=Yes, 0=No
%How many trials to pull for each subsample? 
NumTrialsPulled = 5;
%How many times to do subsampling? We recommend doing at least 10 subsamples and to have the possiblility of using all your data
% (e.g., if you have 150 trials, do 15 subsamples of 10 trials).
NumSubsamples = 100;

%16. Would you like to calculate inter-channel phase synchrony (ICPS or wPLI) in addition to TF?
ICPS_calc = 1; %(1=yes,0=no) 

%17. Would you like to caluclate coherence or weighted phaselagidx?
ICPS_or_wPLI = 0; %(1=coherence, 0=wPLI)

%18. Inter-channel phase synchrony over trials or connectivity over time?
TimeOrTrials = 1; %0 = over time, 1 = over trials
% NOTE: over time calculations will not be able to be subsampled or
% baseline corrected

%19. Type of Connectivity to compute 
ConnectType = 1; %(0= all-to-all connectivity; 1=seed-based connectivity)
%If seed based, choose seed and which electrodes to compute connectivity:
Seed = 'E6';
Elecs4Connect = { 'E14' 'E25'}; 

% 20. Create List of subjects to loop through
subnum = dir([data_location '*.set']); % Use regex to find your files 
subject= {subnum.name};
for ii=1:length(subject)
    subject_list{ii}=subject{ii};   
end

% We recommend checking that your subject list looks like you would expect:
% subject_list

eeglab % Loading EEGLAB

%%%%%%%%%%%%%%%%%%%% COMPUTATIONS BEGIN BELOW HERE %%%%%%%%%%%%%%%%
%% loop through all subject
for sub=1:length(subject_list)
    
    % Initialize objects for this participant:
    timefreqs_data = [];
    phase_data=[];
    ITPS_all=[];
    ICPS_all=[];
    wPLI_all=[];
    EEG=[];
    
    subject = subject_list{sub};
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', sub, subject);

    % Load data
    EEG=pop_loadset('filename', [subject], 'filepath', data_location);
    EEG = pop_selectevent( EEG, 'latency','-.1 <= .1','deleteevents','on');
    
    %Downsample to 250Hz
    EEG = pop_resample(EEG, 250);
    EEG = eeg_checkset( EEG );
    
    %Keep only markers of interest
    EEG = pop_selectevent( EEG, 'Condition',Conds, 'deleteevents','on','deleteepochs','on','invertepochs','off');
    EEG = eeg_checkset( EEG );
    
    %save times
    time = EEG.times;
    
    %save channel locations
    channel_location = EEG.chanlocs;
    
    %% Compute complex morlet wavelet time frequency decomposition
    cd(scripts_location)
    timefreq;
 
end
    
%Save out table with trial numbers
TrialNums_table = struct2table(TrialNums);
writetable(TrialNums_table,[save_location 'TrialNums.csv']);

    
    

