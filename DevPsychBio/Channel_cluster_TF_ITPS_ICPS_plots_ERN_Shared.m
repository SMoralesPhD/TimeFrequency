%% Perform time frequency calculation
clear
clc;

%% Adding this path for Ludwig
addpath('Y:/Toolbox/eeglab13_6_5b/')
eeglab
%% Locations
% data_location = 'Y:\Projects\ECHO\CSD_TimeFreqs\Zoo_3s\'; % This is the original data 
data_location = 'X:/Dropboxes/moraless/ECHO/ITPS_ICPS_redo/'; % This is the most recent
save_data = 'Y:\Projects\ECHO\';
script_location = 'Y:\Projects\ECHO\Scripts';

%% Load channel location
cd(script_location)
load('channel64_location.mat')
channel_location = channel64_location;

%%
clim=[-2 2];
% clim=[-25 25];
freq2plot=[2 15];
time2plot=[-100 500];
%% Subject list
cd(data_location)
subnum=dir('*TF_ITPS_100sub.mat');
% subnum=dir('*TimeFreqs.mat');
sub_list={subnum.name};
for i =1:length(sub_list)
    sub = sub_list{i};
    subject_list{i}= sub(1:15);
end

% subject_list = subject_list(1, 1:41);
% subject_list = subject_list(1, 22:39);

% Selecting the participants in the paper
% Reading in files that are included in the paper
BehData = importdata('Y:/Projects/ECHO/Descriptives/ECHO_ERN_IDS_Paper_01_06_2021.csv');
subject_list = intersect(subject_list, BehData); % Deleting participants on the list


%% Condition name
condition_name = {'Zoo_Resp_TF_ITPS_100sub'}; % Change this for the first and last half of the trials

%% List the group of channels
% List the group of channels
Fz = {'E4', 'E7', 'E54'};
Cz = {'E4', 'E7', 'E54', 'E51', 'E41', 'E21', 'E16'};
Pz = {'E34', 'E36', 'E33', 'E38'};

%% Find indices of the channels
for i=1:length(Fz)
Fz_idx (i)= find(strcmp({channel64_location.labels}, Fz{i}));
end

for i=1:length(Cz)
Cz_idx (i)= find(strcmp({channel64_location.labels}, Cz{i}));
end

for i=1:length(Pz)
Pz_idx (i)= find(strcmp({channel64_location.labels}, Pz{i}));
end


%% Load all observe condition data in a matrix freq x time x channel x subjects
cd(data_location)
for sub = 1:length(subject_list)
    data_condition = condition_name{1};
    subject  = subject_list{sub};
    data_file = [subject, '_', data_condition, '.mat'];
    load(data_file)
    
    % initialize matrices on 1st subject
%     if sub == 1
%         correct_tf_all = zeros([ size(correct_tf_data) length(subject_list)  ]);
%         error_tf_all = zeros([ size(error_tf_data) length(subject_list)  ]);
%         correct_itps_all = zeros([ size(correct_ITPS_data) length(subject_list)  ]);
%         error_itps_all = zeros([ size(error_ITPS_data) length(subject_list)  ]);
%     end
    correct_tf_all(:,:,:,sub) = correct_tf_data_ds;
    error_tf_all(:,:,:,sub) = error_tf_data_ds;
    correct_itps_all(:,:,:,sub) = correct_ITPS_subsamp_baseCorr_ds;
    error_itps_all(:,:,:,sub) = error_ITPS_subsamp_baseCorr_ds;
end




%%%%%%%%%%%%%%%%%%  TF Plotting Starts Here First %%%%%%%%%%%%%%%%%%
%% Make cluster
cd(save_data)

%%
clim=[-1.5 1.5];
freq2plot=[1 15];
time2plot=[-200 600];
timesforplot=time(1:5:750);

% Fz cluster for TF map
for i=1:length(Fz)
correct_tf_Fz   (:,:,i,:)   =   correct_tf_all(:,:,Fz_idx(i),:); 
error_tf_Fz   (:,:,i,:)   =   error_tf_all(:,:,Fz_idx(i),:);

end
%average over the Fz cluster for TF map
correct_tf_Fz_clustered = squeeze(mean(correct_tf_Fz, 3));
error_tf_Fz_clustered = squeeze(mean(error_tf_Fz, 3));


% Plot correct, error, and error-correct conditions
eeglab;
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_tf_Fz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Correct at FCz', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
   % cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_tf_Fz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Error at FCz', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
%figure;
subplot(3,1,3)
% contourf(time, frequency, mean(con_tf_data(:,:,find(strcmp({channel_location.labels}, {'E11'})),:), 3), 20,'linecolor','none');
contourf(timesforplot, frequency, mean(error_tf_Fz_clustered, 3) - mean(correct_tf_Fz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Difference (Error – Correct) at FCz', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);
%%
% Now doing the same with Cz
% Cz
for i=1:length(Cz)
correct_tf_Cz   (:,:,i,:)   =   correct_tf_all(:,:,Cz_idx(i),:); 
error_tf_Cz   (:,:,i,:)   =   error_tf_all(:,:,Cz_idx(i),:);

end

correct_tf_Cz_clustered = squeeze(mean(correct_tf_Cz, 3));
error_tf_Cz_clustered = squeeze(mean(error_tf_Cz, 3));


% Plot correct, error, and error-correct conditions
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_tf_Cz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Cz Correct', 'FontName','Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
    %cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_tf_Cz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Cz Error', 'FontName','Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
 subplot(3,1,3)
contourf(timesforplot, frequency, mean(error_tf_Cz_clustered, 3) - mean(correct_tf_Cz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Cz Diff', 'FontName','Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);

%%
% Now doing the same with Pz
% Pz
for i=1:length(Pz)
correct_tf_Pz   (:,:,i,:)   =   correct_tf_all(:,:,Pz_idx(i),:); 
error_tf_Pz   (:,:,i,:)   =   error_tf_all(:,:,Pz_idx(i),:);

end

correct_tf_Pz_clustered = squeeze(mean(correct_tf_Pz, 3));
error_tf_Pz_clustered = squeeze(mean(error_tf_Pz, 3));


% Plot correct, error, and error-correct conditions
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_tf_Pz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Pz Correct', 'FontName','Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
    %cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_tf_Pz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Pz Difference', 'FontName','Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
 subplot(3,1,3)
contourf(timesforplot, frequency, mean(error_tf_Pz_clustered, 3) - mean(correct_tf_Pz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Pz Diff', 'FontName','Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);

    
%% %%%%%%%%%%%%%%%%%%% TF TOPOPLOT PLOTTING %%%%%%%%%%%%%%%%%%%%%
%Find indices corresponding to what time period that you want to plot
time_windows = [0 300]; %0-300 for paper
freq_windows = [1 4];
Time_Idx=zeros(size(time_windows));
Freqs_Idx = zeros(size(freq_windows));

% Find time indices
for i=1:size(time_windows,1)
    for j=1:2
        [~,Time_Idx(i,j)] = min(abs(timesforplot-time_windows(i,j)));
    end
end

% Fine frequency indices
for i=1:size(freq_windows,1)
    for j=1:2
        [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
    end
end

%Topoplot for each condition
figure;
clim=[-1.5 1.5];

subplot(3,1,1)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot(squeeze(mean(mean(mean(correct_tf_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:,:),1), 2),4)), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title(['Correct ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

subplot(3,1,2)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot(squeeze(mean(mean(mean(error_tf_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:,:),1), 2),4)), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title(['Error ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

subplot(3,1,3)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot((squeeze(mean(mean(mean(error_tf_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:,:),1), 2),4))-(squeeze(mean(mean(mean(correct_tf_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:,:),1), 2),4)))), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title([ 'Error-Correct ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end
    
    
    
%% %%%%%%%%%%%%%%%%%%%%% ITPS MAP PLOTTING starts here %%%%%%%%%%%%%%%%   

% Make cluster
cd(save_data)
clim=[-.1,.1];
% Fz
for i=1:length(Fz)
correct_itps_Fz   (i,:,:,:)   =   correct_itps_all(Fz_idx(i),:,:,:); 
error_itps_Fz   (i,:,:,:)   =   error_itps_all(Fz_idx(i),:,:,:);

end


correct_itps_Fz_clustered = squeeze(mean(correct_itps_Fz, 1));
error_itps_Fz_clustered = squeeze(mean(error_itps_Fz, 1));


% Plot correct, error, and error-correct conditions
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_itps_Fz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Correct at FCz', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    %cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_itps_Fz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Error at FCz', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
%figure;
subplot(3,1,3)
% contourf(time, frequency, mean(con_tf_data(:,:,find(strcmp({channel_location.labels}, {'E11'})),:), 3), 20,'linecolor','none');
contourf(timesforplot, frequency, mean(error_itps_Fz_clustered, 3) - mean(correct_itps_Fz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Difference (Error – Correct) at FCz', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);
%%
% Now doing the same with Cz
% Cz
for i=1:length(Cz)
correct_itps_Cz   (:,:,i,:)   =   correct_itps_all(:,:,Cz_idx(i),:); 
error_itps_Cz   (:,:,i,:)   =   error_itps_all(:,:,Cz_idx(i),:);

end

correct_itps_Cz_clustered = squeeze(mean(correct_itps_Cz, 3));
error_itps_Cz_clustered = squeeze(mean(error_itps_Cz, 3));


% Plot error, correct, and error-correct conditions
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_itps_Cz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Cz Correct', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    %cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_itps_Cz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Cz Error', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
 subplot(3,1,3)
contourf(timesforplot, frequency, mean(error_itps_Cz_clustered, 3) - mean(correct_itps_Cz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Cz Diff', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);

%%
% Now doing the same with Pz
% Pz
for i=1:length(Pz)
correct_itps_Pz   (:,:,i,:)   =   correct_itps_all(:,:,Pz_idx(i),:); 
error_itps_Pz   (:,:,i,:)   =   error_itps_all(:,:,Pz_idx(i),:);

end

correct_itps_Pz_clustered = squeeze(mean(correct_itps_Pz, 3));
error_itps_Pz_clustered = squeeze(mean(error_itps_Pz, 3));


% Plot error, correct, and error-correct
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_itps_Pz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Pz Correct', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    %cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_itps_Pz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Pz Difference', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
 subplot(3,1,3)
contourf(timesforplot, frequency, mean(error_itps_Pz_clustered, 3) - mean(correct_itps_Pz_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Pz Diff', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);

    
%% %%%%%%%%%%%%%%%%%%%%% ITPS TOPOPLOT PLOTTING starts here %%%%%%%%%%%%%%%%   
%Find indices corresponding to what time period that you want to plot
time_windows = [0 300]; %00-300 for paper
freq_windows = [1 4];
Time_Idx=zeros(size(time_windows));
Freqs_Idx = zeros(size(freq_windows));

% Find time indices
for i=1:size(time_windows,1)
    for j=1:2
        [~,Time_Idx(i,j)] = min(abs(timesforplot-time_windows(i,j)));
    end
end

% Fine frequency indices
for i=1:size(freq_windows,1)
    for j=1:2
        [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
    end
end

%Topoplot for each condition
figure;
clim=[-.1 .1];

subplot(3,1,1)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot(squeeze(mean(mean(mean(correct_itps_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),:,Time_Idx(ti,1):Time_Idx(ti,2),:),1), 3),4)), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title(['Correct ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

subplot(3,1,2)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot(squeeze(mean(mean(mean(error_itps_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),:,Time_Idx(ti,1):Time_Idx(ti,2),:),1), 3),4)), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title(['Error ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

subplot(3,1,3)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot((squeeze(mean(mean(mean(error_itps_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),:,Time_Idx(ti,1):Time_Idx(ti,2),:),1), 3),4))-(squeeze(mean(mean(mean(correct_itps_all(Freqs_Idx(fi,1):Freqs_Idx(fi,2),:,Time_Idx(ti,1):Time_Idx(ti,2),:),1), 3),4)))), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title([ 'Error-Correct ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end
        
    
%% Load all ICPS data 
cd(data_location)
condition_name = {'Zoo_Resp_ICPS_100sub'};
for sub = 1:length(subject_list)
    data_condition = condition_name{1};
    subject  = subject_list{sub};
    data_file = [subject, '_', data_condition, '.mat'];
    load(data_file)
    
    % initialize matrices on 1st subject
%     if sub == 1
%         correct_tf_all = zeros([ size(correct_tf_data) length(subject_list)  ]);
%         error_tf_all = zeros([ size(error_tf_data) length(subject_list)  ]);
%         correct_itps_all = zeros([ size(correct_ITPS_data) length(subject_list)  ]);
%         error_itps_all = zeros([ size(error_ITPS_data) length(subject_list)  ]);
%     end
    correct_icps_all(:,:,:,sub) = correct_icps_subsamp_baseCorr_ds;
    error_icps_all(:,:,:,sub) = error_icps_subsamp_baseCorr_ds;
end

%save out all frequencies if you want
output_dir = '/export/data/cdl/Projects/ECHO/CSD_TimeFreqs/Output_data/';
cd(output_dir)
param.date = datestr(now,'mm_dd_yyyy');
save_data=['ECHO_Zoo_Resp_ICPS_100sub_' param.date];
save (save_data,'time','frequency','correct_icps_all','error_icps_all', 'channel_location');
    
  

%% %%%%%%%%%%%%%%%%%%%%% ICPS MAP PLOTTING starts here %%%%%%%%%%%%%%%%   

% Make cluster
cd(save_data)
clim=[-.05,.05];
freq2plot=[1 15];
time2plot=[-200 600];
timesforplot=time(1:5:750);
Fl = {'E59', 'E60', 'E12', 'E13'};
%% Find indices of the channels
for i=1:length(Fl)
Fl_idx (i)= find(strcmp({channel64_location.labels}, Fl{i}));
end


% Fl
for i=1:length(Fl)
correct_icps_Fl   (i,:,:,:)   =   correct_icps_all(Fl_idx(i),:,:,:); 
error_icps_Fl   (i,:,:,:)   =   error_icps_all(Fl_idx(i),:,:,:);

end


correct_icps_Fl_clustered = squeeze(mean(correct_icps_Fl, 1));
error_icps_Fl_clustered = squeeze(mean(error_icps_Fl, 1));


% Plot correct, error, and error-correct conditions
figure; 

subplot(3,1,1)
contourf(timesforplot, frequency, mean(correct_icps_Fl_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim);
     title('Correct at Fl', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    %cbar('vert',0, clim, 5);

%figure;
subplot(3,1,2)
contourf(timesforplot, frequency, mean(error_icps_Fl_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Error at Fl', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
%     cbar('vert',0, clim, 5);
    
%figure;
subplot(3,1,3)
% contourf(time, frequency, mean(con_tf_data(:,:,find(strcmp({channel_location.labels}, {'E11'})),:), 3), 20,'linecolor','none');
contourf(timesforplot, frequency, mean(error_icps_Fl_clustered, 3) - mean(correct_icps_Fl_clustered, 3), 20,'linecolor','none');
    set(gca, 'ylim', freq2plot, 'xlim', time2plot, 'clim', clim); 
    title('Difference (Error – Correct) at Fl', 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal');
    cbar('vert',0, clim, 5);


%% %%%%%%%%%%%%%%%%%%%%% ICPS TOPOPLOT PLOTTING starts here %%%%%%%%%%%%%%%%   
%Find indices corresponding to what time period that you want to plot
time_windows = [0 300];
freq_windows = [1 4];
Time_Idx=zeros(size(time_windows));
Freqs_Idx = zeros(size(freq_windows));

% Find time indices
for i=1:size(time_windows,1)
    for j=1:2
        [~,Time_Idx(i,j)] = min(abs(timesforplot-time_windows(i,j)));
    end
end

% Fine frequency indices
for i=1:size(freq_windows,1)
    for j=1:2
        [~,Freqs_Idx(i,j)] = min(abs(frequency-freq_windows(i,j)));
    end
end

%correct_icps_all_avg = squeeze(mean(correct_icps_all,4));
%error_icps_all_avg = squeeze(mean(error_icps_all,4));

figure;
subplot(3,1,1)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot(squeeze(mean(mean(mean(correct_icps_all(:,Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:),2), 3),4)), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title(['Correct ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

subplot(3,1,2)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot(squeeze(mean(mean(mean(error_icps_all(:,Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:),2), 3),4)), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title(['Error ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

subplot(3,1,3)
for ti =1:size(time_windows,1)
    for fi=1:size(freq_windows,1)
        
        %subplot(1,3,ti)
        topoplot((squeeze(mean(mean(mean(error_icps_all(:,Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:),2), 3),4))-(squeeze(mean(mean(mean(correct_icps_all(:,Freqs_Idx(fi,1):Freqs_Idx(fi,2),Time_Idx(ti,1):Time_Idx(ti,2),:),2), 3),4)))), channel_location,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title([ 'Error-Correct ' num2str(time_windows(ti,1)) '-' num2str(time_windows(ti,2)) 'ms']);
        set(gca, 'FontName','Arial', 'FontSize', 12, 'FontWeight', 'normal')
    end
end

