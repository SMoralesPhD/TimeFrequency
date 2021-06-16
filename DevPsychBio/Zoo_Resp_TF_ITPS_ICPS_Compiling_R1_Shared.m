%% Directories
clear all

data_location = '/data/rosalind/Dropboxes/moraless/ECHO/ITPS_ICPS_redo/';
script_location = '/data/cdl/Projects/ECHO/Scripts/';

% Specify location to save processed data
output_dir = '/data/cdl/Projects/ECHO/CSD_TimeFreqs/Output_data/R1/';

% Load channel location file
% load('/export/data/cdl/Projects/ECHO/Scripts/channel64_location.mat');

% addpath('/export/data/cdl/Toolbox/eeglab13_6_5b/')
% eeglab
%%

% Subject list
cd(data_location)
subnum=dir('*Resp_TF_ITPS_100sub.mat');
sub_list={subnum.name};

for i =1:length(sub_list)
    sub = sub_list{i};
    % subject_list{i}= sub(1:end-11);
    subjects{i}= sub(1:end-4);
end

% % Reading in files to exclude participants based on not completing the
% % whole task or who did not have good performance
% BehData = importdata('/export/data/cdl/Projects/ECHO/Descriptives/ECHO_Zoo_Exclude.csv');
% BehData = strcat(BehData, '_Adjust_CompRem_Zoo_Referenced_Epoched'); % Adding this to match subjects
% subjects = setdiff(subjects, BehData); % Deleting participants on the list


conditions = 2; % This is just for Go vs. No-Go
cutoff = 6; % Maybe the actual names to pass them for naming and stuff later

% Current epoch for Zoo is -1000 to 2000 ms - we should prob make it
% smaller and will need to adjust the sameples below accordingly. 
pIdx =0;
% array to store data (participants x conditions x channels x samples of avg epoch data)
% This is created on the fly in the for loop below. 

Condition_Name = {'error', 'correct'};
channels = 1:60;
times2save = -1000:20:2000;

%Save Theta Frequencies
freq_windows = [4 8];
% Looping to create a file that has subj x cond x time x channels
for subj = 1:length(subjects)
    %load the file that contains all the info for that participant
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', char(subjects(subj)));
    load(char(strcat(subjects(subj), '.mat')), 'correct_tf_data_ds', 'error_tf_data_ds', ...
        'correct_ITPS_subsamp_baseCorr_ds', 'error_ITPS_subsamp_baseCorr_ds', 'frequency');
    
    subj_ICPS  = subjects{subj};
    subj_ICPS = subj_ICPS(1:end-15);
    load(char(strcat(subj_ICPS, '_ICPS_100sub.mat')), 'correct_icps_subsamp_baseCorr_ds', 'error_icps_subsamp_baseCorr_ds');

    TrialNums(subj).Subject = subjects(subj);
    
    Freqs_Idx = zeros(size (freq_windows));
    for i=1:size(freq_windows,1)
        for j=1:2
            [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
        end
    end
    
    % Downsample 
%     times2save_idx = dsearchn(time', times2save');
    
%     correct_tf_data_avg = squeeze(mean(correct_tf_data(Freqs_Idx(1):Freqs_Idx(2),times2save_idx, :),1));
%     error_tf_data_avg = squeeze(mean(error_tf_data(Freqs_Idx(1):Freqs_Idx(2),times2save_idx, :),1));
    
    correct_tf_data_avg = squeeze(mean(correct_tf_data_ds(Freqs_Idx(1):Freqs_Idx(2),:, :),1));
    error_tf_data_avg = squeeze(mean(error_tf_data_ds(Freqs_Idx(1):Freqs_Idx(2),:, :),1));
    
    correct_ITPS_data_avg = squeeze(mean(correct_ITPS_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    error_ITPS_data_avg = squeeze(mean(error_ITPS_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    
    % Changing the order of the matrix because it is flipped 
    correct_ITPS_data_avg = permute(correct_ITPS_data_avg, [2 1]);
    error_ITPS_data_avg = permute(error_ITPS_data_avg, [2 1]);
    
    correct_ICPS_data_avg = squeeze(mean(correct_icps_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    error_ICPS_data_avg = squeeze(mean(error_icps_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    
    % Changing the order of the matrix because it is flipped 
    correct_ICPS_data_avg = permute(correct_ICPS_data_avg, [2 1]);
    error_ICPS_data_avg = permute(error_ICPS_data_avg, [2 1]);
    
    
    % initialize matrices on 1st subject
    if subj==1
        Thetatf_all = zeros([ length(subjects) length(Condition_Name) size(correct_tf_data_avg) ]); % list more variables here as applicable...
        ThetaITPS_all = zeros([ length(subjects) length(Condition_Name) size(correct_ITPS_data_avg) ]); % list more variables here as applicable...
        ThetaICPS_all = zeros([ length(subjects) length(Condition_Name) size(correct_ICPS_data_avg) ]); % list more variables here as applicable...

    end
    
    % Specify conditions here
    Thetatf_all(subj,1,:,:) = correct_tf_data_avg;
    Thetatf_all(subj,2,:,:) = error_tf_data_avg;
    
    ThetaITPS_all(subj,1,:,:) = correct_ITPS_data_avg;
    ThetaITPS_all(subj,2,:,:) = error_ITPS_data_avg;
    
    ThetaICPS_all(subj,1,:,:) = correct_ICPS_data_avg;
    ThetaICPS_all(subj,2,:,:) = error_ICPS_data_avg;  
end

param.date = datestr(now,'mm_dd_yyyy');

cd(output_dir)
save(['Thetatf_all_Zoo_Resp_4_8_3s_50Hz_' param.date '.mat'], 'Thetatf_all')
save(['ThetaITPS_all_Zoo_Resp_4_8_3s_50Hz_' param.date '.mat'], 'ThetaITPS_all')
save(['ThetaICPS_all_Zoo_Resp_4_8_3s_50Hz_' param.date '.mat'], 'ThetaICPS_all')


table_events = struct2table(TrialNums);
writetable(table_events, ['/data/cdl/Projects/ECHO/Descriptives/Trial_Numbers_tf_ITPS_ICPS_all_Zoo_Resp_4_8_3s_50Hz_' param.date '.csv']);
% writetable(table_events, ['/export/data/cdl/Projects/ECHO/Descriptives/Trial_Numbers_tf_ITPS_ICPS_all_Zoo_Resp_4_8_3s_50Hz_' param.date '.csv']);


%Save Delta Frequencies
freq_windows = [1 4];
cd(data_location)

% Looping to create a file that has subj x cond x time x channels
for subj = 1:length(subjects)
    %load the file that contains all the info for that participant
    fprintf('\n\n\n*** Processing subject %s ***\n\n\n', char(subjects(subj)));
    load(char(strcat(subjects(subj), '.mat')), 'correct_tf_data_ds', 'error_tf_data_ds', ...
        'correct_ITPS_subsamp_baseCorr_ds', 'error_ITPS_subsamp_baseCorr_ds', 'frequency');
    subj_ICPS  = subjects{subj};
    subj_ICPS = subj_ICPS(1:end-15);
    load(char(strcat(subj_ICPS, '_ICPS_100sub.mat')), 'correct_icps_subsamp_baseCorr_ds', 'error_icps_subsamp_baseCorr_ds');

    TrialNums(subj).Subject = subjects(subj);
    
    Freqs_Idx = zeros(size (freq_windows));
    for i=1:size(freq_windows,1)
        for j=1:2
            [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
        end
    end
    
    % Downsample 
%     times2save_idx = dsearchn(time', times2save');
    
%     correct_tf_data_avg = squeeze(mean(correct_tf_data(Freqs_Idx(1):Freqs_Idx(2),times2save_idx, :),1));
%     error_tf_data_avg = squeeze(mean(error_tf_data(Freqs_Idx(1):Freqs_Idx(2),times2save_idx, :),1));
    
    correct_tf_data_avg = squeeze(mean(correct_tf_data_ds(Freqs_Idx(1):Freqs_Idx(2),:, :),1));
    error_tf_data_avg = squeeze(mean(error_tf_data_ds(Freqs_Idx(1):Freqs_Idx(2),:, :),1));
    
    correct_ITPS_data_avg = squeeze(mean(correct_ITPS_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    error_ITPS_data_avg = squeeze(mean(error_ITPS_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    
    % Changing the order of the matrix because it is flipped 
    correct_ITPS_data_avg = permute(correct_ITPS_data_avg, [2 1]);
    error_ITPS_data_avg = permute(error_ITPS_data_avg, [2 1]);
    
    correct_ICPS_data_avg = squeeze(mean(correct_icps_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    error_ICPS_data_avg = squeeze(mean(error_icps_subsamp_baseCorr_ds(:,Freqs_Idx(1):Freqs_Idx(2), :),2));
    
    % Changing the order of the matrix because it is flipped 
    correct_ICPS_data_avg = permute(correct_ICPS_data_avg, [2 1]);
    error_ICPS_data_avg = permute(error_ICPS_data_avg, [2 1]);
    
    
    % initialize matrices on 1st subject
    if subj==1
        Deltatf_all = zeros([ length(subjects) length(Condition_Name) size(correct_tf_data_avg) ]); % list more variables here as applicable...
        DeltaITPS_all = zeros([ length(subjects) length(Condition_Name) size(correct_ITPS_data_avg) ]); % list more variables here as applicable...
        DeltaICPS_all = zeros([ length(subjects) length(Condition_Name) size(correct_ICPS_data_avg) ]); % list more variables here as applicable...

    end
    
    % Specify conditions here
    Deltatf_all(subj,1,:,:) = correct_tf_data_avg;
    Deltatf_all(subj,2,:,:) = error_tf_data_avg;
    
    DeltaITPS_all(subj,1,:,:) = correct_ITPS_data_avg;
    DeltaITPS_all(subj,2,:,:) = error_ITPS_data_avg;
    
    DeltaICPS_all(subj,1,:,:) = correct_ICPS_data_avg;
    DeltaICPS_all(subj,2,:,:) = error_ICPS_data_avg;  
end

param.date = datestr(now,'mm_dd_yyyy');

cd(output_dir)
save(['Deltatf_all_Zoo_Resp_1_4_3s_50Hz_' param.date '.mat'], 'Deltatf_all')
save(['DeltaITPS_all_Zoo_Resp_1_4_3s_50Hz_' param.date '.mat'], 'DeltaITPS_all')
save(['DeltaICPS_all_Zoo_Resp_1_4_3s_50Hz_' param.date '.mat'], 'DeltaICPS_all')


