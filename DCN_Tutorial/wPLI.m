
% INPUT
% data = single trial wavelet decomposed signal; frequency x time x trial x channel
%The input will come from the output of the strial_timefreq_decomposition_script.m
data=[];
data = phase_data_strial_ch;
                                               
% OUTPUT: wPLI_all
%over time, all-to-all connectivity; frequency x trials x channel x channel
%over trials, all-to-all connectivity; frequency x time x channel x channel
%over time, seed-based; frequency x trials x non-seed channel 
%over trials, seed-based; frequency x time x non-seed channel 

%% Initialize output matrices depending on type of connectivity
if TimeOrTrials == 0
    if ConnectType == 0
        %over time, all-to-all
        wPLI_all       = zeros(size(data, 1), size(data, 3), size(data, 4), size(data, 4));
    elseif ConnectType == 1
        %over time, seed-based
        wPLI_all       = zeros(size(data, 1), size(data, 3), length(Elecs4Connect));
    end
elseif TimeOrTrials ==1
    if ConnectType == 0
        %over trials, all-to-all
        wPLI_all       = zeros(size(data, 1), size(data, 2), size(data, 4), size(data, 4));
    elseif ConnectType == 1
        %over trials, seed-based
        wPLI_all       = zeros(size(data, 1), size(data, 2), length(Elecs4Connect));
    end
end


%% Connectivity Computations
%%Compute all-to-all connectivity
if ConnectType == 0
    if TimeOrTrials ==0
        fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) ***\n\n\n', sub, subject);
        for chani=1:size(data, 4)
            for chanj=chani:size(data, 4)
                %% take cross-spectral density
                crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                % take imaginary part of signal only
                crossspecden_imag = imag(crossspecden);
                %% weighted phase-lag index (shortcut, as implemented in Cohen's book)
                for freq=1:size(data, 1)
                    wPLI_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),2))./mean(abs(crossspecden_imag(freq,:,:)),2);
                    wPLI_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),2))./mean(abs(crossspecden_imag(freq,:,:)),2);
                end %end loop through frequencies
            end %end second loop through channels
        end %end first loop through channels
        
    elseif TimeOrTrials==1
        %with subsampling
        if Subsample == 1
            %initialize matrix with subsamples
            wPLI_subsamples = zeros(NumSubsamples, size(data,1), size(data,2), size(data,4), size(data,4));
            %Begin subsampling
            for samp=1:NumSubsamples
                fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                crossspecden=[]; crossspecden_imag=[]; subtrials=[]; crossspecden_imag_temp=[]; weighted_phaselagidx_temp=[];
                for chani=1:size(data, 4)
                    for chanj=chani:size(data, 4)
                        %% take cross-spectral density
                        crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                        % take imaginary part of signal only
                        crossspecden_imag = imag(crossspecden);
                        %Get indices of trials for this subsample
                        subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                        %Index into those trials and pull them out for wPLI analyses
                        crossspecden_imag_temp = crossspecden_imag(:,:,subtrials);
                        for freq=1:size(data, 1)
                            weighted_phaselagidx_temp(freq,:,chani,chanj) = abs( mean( abs(crossspecden_imag_temp(freq,:,:)).*sign(crossspecden_imag_temp(freq,:,:)),3))./mean(abs(crossspecden_imag_temp(freq,:,:)),3);
                            weighted_phaselagidx_temp(freq,:,chanj,chani) = abs( mean( abs(crossspecden_imag_temp(freq,:,:)).*sign(crossspecden_imag_temp(freq,:,:)),3))./mean(abs(crossspecden_imag_temp(freq,:,:)),3);
                        end %end loop through frequencies
                    end %end second loop through channels
                end %end first loop through channels
            %create matrix of subsamples for wPLI - samp x freq x time
            wPLI_subsamples(samp,:,:,:,:) = weighted_phaselagidx_temp;
            end %end loop through subsamples
            %average over subsamples - this is final wPLI
            wPLI_all = squeeze(mean(wPLI_subsamples,1));
            
            % all-to-all connectivity WITHOUT subsampling
        elseif Subsample == 0
            fprintf('\n\n\n*** Calculating wPLI for subject %d (%s)***\n\n\n', sub, subject);
            for chani=1:size(data, 4)
                for chanj=chani:size(data, 4)
                    %% take cross-spectral density
                    crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                    % take imaginary part of signal only
                    crossspecden_imag = imag(crossspecden);
                    %% weighted phase-lag index (shortcut, as implemented in Cohen's book)
                    for freq=1:size(data, 1)
                        wPLI_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),3))./mean(abs(crossspecden_imag(freq,:,:)),3);
                        wPLI_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),3))./mean(abs(crossspecden_imag(freq,:,:)),3);
                    end %end loop through frequencies
                end %end second loop through channels
            end %end first loop through channels
        end %end if subsampling
    end %end if time or trials
    
    %Baseline Correct, Downsample, and Save Data
    if TimeOrTrials == 0 %over time does not baseline correct or downsample
        if RestorEvent==1
            %save out trial averaged and baseline corrected wPLI data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection'];
            save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
        elseif RestorEvent==0
            %save out trial averaged and baseline corrected wPLI data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection_','condition',Conds{cond}];
            save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
        end
    elseif TimeOrTrials ==1
        %Baseline Correction
        if BaselineCorrect == 1
            
            %% baseline time indices
            basetimeidx   = dsearchn(EEG.times', BaselineTime');
            
            Baseline = squeeze(mean(wPLI_all(:,basetimeidx(1):basetimeidx(end),:,:),2));
            %loop through samples
            for t=1:size(ICPS_all,2)
                wPLI_blncorr(:,t,:,:) = squeeze(wPLI_all(:,t,:,:)) - Baseline;
            end %end loop through frequencies
            
            if RestorEvent==1
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    wPLI_blncorr=ICPS_baselinecorr(:,1:2:size(ITPS_baselinecorr,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %if event-related, add condition name to saved file
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample==1
                    %Downsample
                    wPLI_blncorr=wPLI_blncorr(:,1:2:size(wPLI_blncorr,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            end
        elseif BaselineCorrect == 0
            if RestorEvent==1 %Resting
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                    save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                    save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %Event-Related
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end %end if downsampling
            end %end if resting or event
        end %end if you want baseline correction
    end %end over time or trials
    
%SEED-BASED
elseif ConnectType == 1
    if TimeOrTrials==0
        fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) ***\n\n\n', sub, subject);
        
        %Find index of seed electrode
            Seed_idx = find(strcmpi({EEG.chanlocs.labels},Seed));
            
            % Find indices of the non-seed channels
            for i=1:length(Elecs4Connect)
                Elecs_idx (i)= find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
            end

        % take cross-spectral density
        for chanj=1:length(Elecs4Connect)
            % take cross-spectral density between two channels - one being the seed for ICPS
            crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
            % take imaginary part of signal only
            crossspecden_imag = imag(crossspecden);
            %initialize matrix
            weighted_phaselagidx_seed = zeros(size(data, 1), size(data, 3));
            for freq=1:size(data, 1)
                %wPLI for nonsubsampled seed-based over time
                weighted_phaselagidx_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden_imag(freq,:,:))),2));
            end %end loop through frequencies
            wPLI_all(:,:,chanj) = weighted_phaselagidx_seed;
        end % end loop through channels
    elseif TimeOrTrials==1
        %subsampling seed-based connectivity
        if Subsample == 1
            %initialize subsamples matrix
            wPLI_subsamples = zeros(NumSubsamples, size(data,1), size(data,2),length(Elecs4Connect));
            %Find index of seed electrode
            Seed_idx = find(strcmpi({EEG.chanlocs.labels},Seed));
            
            % Find indices of the non-seed channels
            for i=1:length(Elecs4Connect)
                Elecs_idx (i)= find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
            end
                    
            %Start subsampling
            for samp=1:NumSubsamples
                fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                crossspecden=[]; crossspecden_imag=[]; subtrials=[]; crossspecden_imag_temp=[];
                %Initialize matrix
                wPLI_ch = zeros(size(data,1), size(data,2),length(Elecs4Connect));
                %loop through channels
                for chanj=1:length(Elecs4Connect)
                    %Take cross spectral density
                    crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
                    % take imaginary part of signal only
                    crossspecden_imag = imag(crossspecden);
                    %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                    %Get indices of trials for this subsample
                    subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                    %Index into those trials and pull them out for wPLI analyses
                    crossspecden_imag_temp = crossspecden_imag(:,:,subtrials);
                    %Initialize matrix
                    weighted_phaselagidx_seed = zeros(size(data, 1), size(data, 2));
                    %Loop through frequencies
                    for freq=1:size(data, 1)
                        % wPLI
                        weighted_phaselagidx_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden_imag_temp(freq,:,:))),3));
                    end %end loop through frequencies
                    wPLI_ch(:,:,chanj) = weighted_phaselagidx_seed;
                end %end loop through channels
                wPLI_subsamples(samp,:,:,:) = wPLI_ch;
            end%end loop through subsamples
            %average over subsamples - wPLI for subsampled seed-based connectivity
            wPLI_all = squeeze(mean(wPLI_subsamples,1));
            
            %SEED-BASED WITHOUT subsampling
        elseif Subsample == 0
            fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) ***\n\n\n', sub, subject);
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
                % take imaginary part of signal only
                crossspecden_imag = imag(crossspecden);
                %Initialize matrix
                weighted_phaselagidx_seed = zeros(size(data, 1), size(data, 2));
                %Loop through frequencies
                for freq=1:size(data, 1)
                    %wPLI for nonsubsampled seed-based over trials
                    weighted_phaselagidx_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden_imag(freq,:,:))),3));
                end %end loop through frequencies
                wPLI_all(:,:,chanj) = weighted_phaselagidx_seed;
            end % end loop through channels
        end % end if subsampling statement
    end %end if time or trials
    
    %Baseline Correct, Downsample, and Save Data
    if TimeOrTrials == 0 %over time does not baseline correct or downsample
        if RestorEvent==1
            %save out trial averaged and baseline corrected wPLI data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection'];
            save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
        elseif RestorEvent==0
            %save out trial averaged and baseline corrected wPLI data for this subject
            save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection_','condition',Conds{cond}];
            save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
        end
    elseif TimeOrTrials ==1
        %Baseline Correction
        if BaselineCorrect == 1
            
            %% baseline time indices
            basetimeidx   = dsearchn(EEG.times', BaselineTime');
            
            %Initialize wPLI_blncorr
            wPLI_blncorr = zeros(size(wPLI_all));
            
            %Baseline Correct
           for chanj=1:size(wPLI_all,3)
                for fi = 1:size(wPLI_all,1)
                    wPLI_blncorr(fi,:,chanj) = ( wPLI_all(fi,:,chanj) - mean(wPLI_all(fi,basetimeidx(1):basetimeidx(end),chanj)));
                end
            end
            
            if RestorEvent==1
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    wPLI_blncorr=ICPS_baselinecorr(:,1:2:size(ITPS_baselinecorr,2),:,:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %if event-related, add condition name to saved file
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample==1
                    %Downsample
                    wPLI_blncorr=wPLI_blncorr(:,1:2:size(wPLI_blncorr,2),:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition',Conds{cond}];
                    save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            end
        elseif BaselineCorrect == 0
            if RestorEvent==1 %Resting
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                    save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged, downsampled, baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                    save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end
            elseif RestorEvent==0 %Event-Related
                if Downsample ==0
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                elseif Downsample ==1
                    %Downsample
                    wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:);
                    
                    %Downsample time variable
                    ds_time = downsample(time,2);
                    
                    %save out trial averaged and baseline corrected wPLI data for this subject
                    save_data =[save_location, subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition',Conds{cond}];
                    save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                end %end if downsampling
            end %end if resting or event
        end %end if you want baseline correction
    end %end if time or trials
end %end if connectivity type loops
