
% INPUT
% data = single trial wavelet decomposed signal; frequency x time x trial x channel
data=[];
data = phase_data_strial_ch;
                                               
% OUTPUT: ITPS_all
%frequency x time x channels 

%% Initialize output matrices depending on type of connectivity
ITPS_all       = zeros(size(data, 1), size(data, 2), size(data, 4));
    

%% ITPS Computations
%with subsampling
if Subsample == 1
    %initialize matrix with subsamples
    ITPS_ch = zeros(size(data,1), size(data,2), size(data,4));
    ITPS_subsamples = zeros(NumSubsamples, size(data,1), size(data,2), size(data,4));
    %Begin subsampling
    for samp=1:NumSubsamples
        fprintf('\n\n\n*** Calculating ITPS for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);    
        for chan=1:size(data, 4)
            ITPS=[];
            for freq=1:size(data, 1)
                subtrials=[]; data_temp=[];
                %Get indices of trials for this subsample
                subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                %Index into those trials and pull them out for wPLI analyses
                data_temp = squeeze(data(freq,:,subtrials,chan));
                ITPS(freq,:) = abs(mean( exp(1i*angle(data_temp)),2));
            end %end loop through frequencies
            ITPS_ch(:,:,chan) = ITPS;
        end %end loop through channels
        %create matrix of subsamples for ITPS - samp x freq x time
        ITPS_subsamples(samp,:,:,:) = ITPS_ch;
    end %end loop through subsamples
    %average over subsamples - this is final ITPS
    ITPS_all = squeeze(mean(ITPS_subsamples,1));
     
    
% WITHOUT subsampling
elseif Subsample == 0
fprintf('\n\n\n*** Calculating ITPS for subject %d (%s) ***\n\n\n', sub, subject);    
    for chan=1:size(data, 4)
        ITPS=[];
        %loop through frequencies
        for freq=1:size(data, 1)
            % pull out data needed
            data_temp=[];
            data_temp = squeeze(data(freq,:,:,chan));
            ITPS(freq,:) = abs(mean( exp(1i*angle(data_temp)),2));
        end %end loop through frequencies
        ITPS_all(:,:,chan) = ITPS;
     end %end loop through channels
end %end loop through subsampling

     %Baseline Correction
     if BaselineCorrect == 1
         
         %% baseline time indices
         basetimeidx   = dsearchn(EEG.times', BaselineTime');
         
         %Initialize ITPS_blncorr
         ITPS_blncorr = zeros(size(ITPS_all));
            
         for chanj=1:size(ITPS_all,3)
                for fi = 1:size(ITPS_all,1)
                    ITPS_blncorr(fi,:,chanj) = ( ITPS_all(fi,:,chanj) - mean(ITPS_all(fi,basetimeidx(1):basetimeidx(end),chanj)));
                end
         end
            
         if RestorEvent==1
             
             if Downsample==0
                 %save out trial averaged and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected'];
                 save(save_data, 'ITPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
             elseif Downsample==1
                 %Downsample
                 ITPS_blncorr = ITPS_blncorr(:,1:2:size(ITPS_blncorr,2),:);
                 
                 %Downsample time variable
                 ds_time = downsample(time,2);
                 
                 %save out trial averaged, downsampled, and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected'];
                 save(save_data, 'ITPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
             end
         elseif RestorEvent==0
             if Downsample ==0
                 %save out trial averaged and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected_','condition',Conds{cond}];
                 save(save_data, 'ITPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
             elseif Downsample ==1
                 %Downsample
                 ITPS_blncorr = ITPS_blncorr(:,1:2:size(ITPS_blncorr,2),:);
                 
                 %Downsample time variable
                 ds_time = downsample(time,2);
                 
                 %save out trial averaged and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected_','condition',Conds{cond}];
                 save(save_data, 'ITPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
             end
         end
         
     elseif BaselineCorrect == 0
         
         if RestorEvent==1
             if Downsample ==0
                 %save out trial averaged and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection'];
                 save(save_data, 'ITPS_all', 'frequency', 'time','channel_location', '-v7.3');
             elseif Downsample == 1
                 %Downsample
                 ITPS_all = ITPS_all(:,1:2:size(ITPS_all,2),:);
                 
                 %Downsample time variable
                 ds_time = downsample(time,2);
                 
                 %save out trial averaged and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection'];
                 save(save_data, 'ITPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
             end
         elseif RestorEvent==0 %if event-related, then save condition name in save file
             if Downsample ==0
                 %save out trial averaged and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection_','condition',Conds{cond}];
                 save(save_data, 'ITPS_all', 'frequency', 'time','channel_location', '-v7.3');
             elseif Downsample ==1
                 %Downsample
                 ITPS_all = ITPS_all(:,1:2:size(ITPS_all,2),:);
                 
                 %Downsample time variable
                 ds_time = downsample(time,2);
                 
                 %save out trial averaged, downsampled, and baseline corrected ITPS data for this subject
                 save_data =[save_location, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection_','condition',Conds{cond}];
                 save(save_data, 'ITPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
             end %end if Downsampling should be done
         end %end if resting or event
     end %end if you want baseline correction


subtrials=[];
data_temp=[];
ITPS=[];
ITPS_ch=[];
ITPS_subsamples=[];
ITPS_all=[];
ITPS_blncorr=[];
