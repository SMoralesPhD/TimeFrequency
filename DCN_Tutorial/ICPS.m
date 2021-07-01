
% INPUT
% data = single trial wavelet decomposed signal; for resting, frequency x time x trial x channel
                                               % for event-related, condition x frequency x time x trial x channel
data=[];
data = phase_data_strial_ch;
                                               
% OUTPUT: ICPS_all
%Resting, over time, all-to-all connectivity; frequency x trials x channel x channel
%Resting, over trials, all-to-all connectivity; frequency x time x channel x channel
%Resting, over time, seed-based; frequency x trials x non-seed channel 
%Resting, over trials, seed-based; frequency x time x non-seed channel 

%Event-Related, over time, all-to-all connectivity; condition x frequency x trials x channel x channel
%Event-Related, over trials, all-to-all connectivity; condition x frequency x time x channel x channel
%Event-Related, over time, seed-based; condition x frequency x trials x non-seed channel
%Event-Related, over trials, seed-based; condition x frequency x time x non-seed channel

%% Initialize output matrices depending on type of connectivity
if TimeOrTrials == 0 %over time 
    if ConnectType == 0
        %Resting, over time, all-to-all
        ICPS_all       = zeros(size(data, 1), size(data, 3), size(data, 4), size(data, 4));
    elseif ConnectType == 1
        %Resting, over time, seed-based
        ICPS_all       = zeros(size(data, 1), size(data, 3), length(Elecs4Connect));
    end
elseif TimeOrTrials ==1 %over trials
    if ConnectType == 0
        %Resting, over trials, all-to-all
        ICPS_all       = zeros(size(data, 1), size(data, 2), size(data, 4), size(data, 4));
    elseif ConnectType == 1
        %Resting, over trials, seed-based
        ICPS_all       = zeros(size(data, 1), size(data, 2), length(Elecs4Connect));
    end
end

%% Connectivity Computations
%%Compute all-to-all connectivity
if ConnectType == 0
    if TimeOrTrials==0
        fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
        for chani=1:size(data, 4)
            for chanj=chani:size(data, 4)
                %% take cross-spectral density
                crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                %% ICPS
                for freq=1:size(data, 1)
                    ICPS_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),2))./mean(abs(crossspecden(freq,:,:)),2);
                    ICPS_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),2))./mean(abs(crossspecden(freq,:,:)),2);
                end %end loop through frequencies
            end %end second loop through channels
        end %end first loop through channels
        
    elseif TimeOrTrials==1
        %with subsampling
        if Subsample == 1
            ICPS_subsamples = zeros(NumSubsamples, size(data,1), size(data,2), size(data,4), size(data,4));
            %Begin subsampling
            for samp=1:NumSubsamples
                fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                crossspecden=[]; subtrials=[]; crossspecden_temp=[]; ICPS_temp=[];
                for chani=1:size(data, 4)
                    for chanj=chani:size(data, 4)
                        %% take cross-spectral density
                        crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                        %Get indices of trials for this subsample
                        subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                        %Index into those trials and pull them out for wPLI analyses
                        crossspecden_temp = crossspecden(:,:,subtrials);
                        for freq=1:size(data, 1)
                            ICPS_temp(freq,:,chani,chanj) = abs( mean( abs(crossspecden_temp(freq,:,:)).*sign(crossspecden_temp(freq,:,:)),3))./mean(abs(crossspecden_temp(freq,:,:)),3);
                            ICPS_temp(freq,:,chanj,chani) = abs( mean( abs(crossspecden_temp(freq,:,:)).*sign(crossspecden_temp(freq,:,:)),3))./mean(abs(crossspecden_temp(freq,:,:)),3);
                        end
                        
                    end %end second loop through channels
                end %end first loop through channels
                %create matrix of subsamples for wPLI - samp x freq x time x channel x channel
                ICPS_subsamples(samp,:,:,:,:) = ICPS_temp;
            end %end loop through subsamples
            %average over subsamples - this is final ICPS
            ICPS_all = squeeze(mean(ICPS_subsamples,1));
            
            % WITHOUT subsampling
        elseif Subsample == 0
            fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
            for chani=1:size(data, 4)
                for chanj=chani:size(data, 4)
                    %% take cross-spectral density
                    crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                    %% ICPS
                    for freq=1:size(data, 1)
                        ICPS_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),3))./mean(abs(crossspecden(freq,:,:)),3);
                        ICPS_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),3))./mean(abs(crossspecden(freq,:,:)),3);
                    end %end loop through frequencies
                end %end second loop through channels
            end %end first loop through channels
        end %end if statement for subsample
    end %end over time or over trials
    
    %Baseline Correct, Downsample, and Save Data
    if TimeOrTrials == 0 %over time does not baseline correct or downsample
        if RestorEvent==1
            %save out trial averaged and baseline corrected ICPS data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection'];
            save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
        elseif RestorEvent==0
            %save out trial averaged and baseline corrected ICPS data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection_','condition',Conds{cond}];
            save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
        end
    elseif TimeOrTrials ==1
        %Baseline Correction
        if BaselineCorrect == 1
            
            %% baseline time indices
            basetimeidx   = dsearchn(EEG.times', BaselineTime');
            
            Baseline = squeeze(mean(ICPS_all(:,basetimeidx(1):basetimeidx(end),:,:),2));
            %loop through samples
            for t=1:size(ICPS_all,2)
                ICPS_blncorr(:,t,:,:) = squeeze(ICPS_all(:,t,:,:)) - Baseline;
            end %end loop through frequencies
            
            if RestorEvent==1
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_baselinecorrected'];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_blncorr=ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_baselinecorrected'];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %if event-related, add condition name to saved file
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample==1
                    %Downsample
                    ICPS_blncorr=ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            end
        elseif BaselineCorrect == 0
            if RestorEvent==1 %Resting
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection'];
                    save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection'];
                    save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %Event-Related
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end %end if downsample
            end %end if resting or event
        end %end if you want baseline correction
    end %end if time or trials
    
    
    
    %SEED-BASED
elseif ConnectType == 1
    if TimeOrTrials == 0
        fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
            %Find index of seed electrode
            Seed_idx = find(strcmpi({EEG.chanlocs.labels},Seed));
            
            % Find indices of the non-seed channels
            for i=1:length(Elecs4Connect)
                Elecs_idx (i)= find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
            end
         %loop through channels
         for chanj=1:length(Elecs4Connect)
            % take cross-spectral density between two channels - one being the seed for ICPS
            crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
            %Initialize matrix
            ICPS_seed = zeros(size(data, 1), size(data, 3));
            %Loop through frequencies
            for freq=1:size(data, 1)
                %wPLI for nonsubsampled seed-based over trials
                ICPS_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden(freq,:,:))),2));
            end %end loop through frequencies
            ICPS_all(:,:,chanj) = ICPS_seed;
        end % end loop through channels
    elseif TimeOrTrials==1
        %WITH subsampling
        if Subsample == 1
            %initialize matrix with subsamples
            ICPS_ch = zeros(size(data,1), size(data,2),length(Elecs4Connect));
            ICPS_subsamples = zeros(NumSubsamples, size(data,1), size(data,2),length(Elecs4Connect));
            %Start subsampling
            for samp=1:NumSubsamples
                fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                crossspecden=[]; crossspecden_imag=[]; subtrials=[]; crossspecden_imag_temp=[]; 
                
                %Find index of seed electrode
                    Seed_idx = find(strcmpi({EEG.chanlocs.labels},Seed));
                    
                    % Find indices of the non-seed channels
                    for i=1:length(Elecs4Connect)
                        Elecs_idx (i)= find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                    end 
                
                %loop through channels
                for chanj=1:length(Elecs4Connect)
                    % take cross-spectral density between two channels - one being the seed for ICPS
                    crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj)))); %
                    %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                    %Get indices of trials for this subsample
                    subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                    %Index into those trials and pull them out for wPLI analyses
                    crossspecden_temp = crossspecden(:,:,subtrials);
                    %Initialize matrix
                    ICPS_seed = zeros(size(data, 1), size(data, 2));
                    %Loop through frequencies
                    for freq=1:size(data, 1)
                        ICPS_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden_temp(freq,:,:))),3));
                    end %end loop through frequencies
                    %create matrix with channelssamp x freq x time
                    ICPS_ch(:,:,chanj) = ICPS_seed;
                end %end loop through channels
                ICPS_subsamples(samp,:,:,:) = ICPS_ch;
            end%end loop through subsamples
            %average over subsamples - wPLI for subsampled seed-based connectivity
            ICPS_all = squeeze(mean(ICPS_subsamples_ch,1));
            %SEED-BASED WITHOUT subsampling
        elseif Subsample == 0
            fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
            % take cross-spectral density
            for chanj=1:length(Elecs4Connect)
                %Find index of seed electrode
                Seed_idx = find(strcmpi({EEG.chanlocs.labels},Seed));
                
                % Find indices of the non-seed channels
                for i=1:length(Elecs4Connect)
                    Elecs_idx (i)= find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                end
                
                % take cross-spectral density between two channels - one being the seed for ICPS
                crossspecden = squeeze(data(:,:,:,Seed_ids) .* conj(data(:,:,:,Elecs_idx(chanj))));
                %Initialize matrix
                ICPS_seed = zeros(size(data,1), size(data,2));
                for freq=1:size(data, 1)
                    %wPLI for nonsubsampled seed-based over trials
                    ICPS_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden(freq,:,:))),3));
                end %end loop through frequencies
                ICPS_all(:,:,chanj) = ICPS_seed;
            end % end loop through channels
        end % end if subsampling statement
    end %end if TimeOrTrials
    
    %Baseline Correct, Downsample, and Save Data
    if TimeOrTrials == 0 %over time does not baseline correct or downsample
        if RestorEvent==1
            %save out trial averaged and baseline corrected ICPS data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtime_nobaselinecorrection'];
            save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
        elseif RestorEvent==0
            %save out trial averaged and baseline corrected ICPS data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtime_nobaselinecorrection_','condition',Conds{cond}];
            save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
        end
    elseif TimeOrTrials ==1
        %Baseline Correction
        if BaselineCorrect == 1
            
            %% baseline time indices
            basetimeidx   = dsearchn(EEG.times', BaselineTime');
            
            %Initialize ICPS_blncorr
            ICPS_blncorr = zeros(size(ICPS_all));
            
            %Baseline correct
            for chanj=1:size(ICPS_all,3)
                for fi = 1:size(ICPS_all,1)
                    ICPS_blncorr(fi,:,chanj) = ( ICPS_all(fi,:,chanj) - mean(ICPS_all(fi,basetimeidx(1):basetimeidx(end),chanj)));
                end
            end

            if RestorEvent==1
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected'];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_blncorr = ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected'];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_blncorr = ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            end
        elseif BaselineCorrect == 0
            if RestorEvent==1 %resting
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection'];
                    save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection'];
                    save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %Event-Related
                if Downsample ==0
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected ICPS data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end %end if downsampling
            end %end if resting or event
        end %end if you want baseline correction
    end %end if time or trials
    
end %end if connectivity type loops

ICPS_all=[];
ICPS_blncorr=[];
Baseline=[];
crossspecden=[]; 
subtrials=[];
crossspecden_temp=[];
ICPS_temp=[];
ICPS_seed=[];
ICPS_ch=[];
ICPS_subsamples=[];


