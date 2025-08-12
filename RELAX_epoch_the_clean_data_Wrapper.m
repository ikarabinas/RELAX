%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2022) Neil Bailey

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see https://www.gnu.org/licenses/.

%% RELAX_Wrapper:
function [OutlierParticipantsToManuallyCheck,EpochRejections,RELAX_epoching_cfg] = RELAX_epoch_the_clean_data_Wrapper (RELAX_epoching_cfg)

% Load pre-processing statistics file for these participants if it already
% exists:
RELAX_epoching_cfg.OutputPath=[RELAX_epoching_cfg.CleanedPath filesep 'Epoched'];
mkdir(RELAX_epoching_cfg.OutputPath);

WarningAboutFileNumber=0;
if size(RELAX_epoching_cfg.FilesToProcess,2) > size(RELAX_epoching_cfg.files,2)
    RELAX_epoching_cfg.FilesToProcess=RELAX_epoching_cfg.FilesToProcess(1,1):size(RELAX_epoching_cfg.files,2);
    WarningAboutFileNumber=1;
end

%% Loop selected files in the directory list:
for FileNumber=RELAX_epoching_cfg.FilesToProcess(1,1:size(RELAX_epoching_cfg.FilesToProcess,2))
    
    RELAX_epoching_cfg.filename=RELAX_epoching_cfg.files{FileNumber};
    clearvars -except 'FilesWithoutConvergence' 'RELAX_epoching_cfg' 'FileNumber' 'FileName' 'Participant_IDs' 'Medianvoltageshiftwithinepoch' 'EpochRejections' 'MatCompatible_Participant_IDs' 'WarningAboutFileNumber';
    %% Load data (assuming the data is in EEGLAB .set format):

    if strcmpi(RELAX_epoching_cfg.all_data_in_1_folder_or_BIDS_format,'BIDS')
        cd(RELAX_epoching_cfg.folders{FileNumber});
    end

    EEG = pop_loadset(RELAX_epoching_cfg.filename);
    FileName = extractBefore(RELAX_epoching_cfg.files{FileNumber},".");
    Participant_IDs{1,FileNumber} = extractBefore(RELAX_epoching_cfg.files{FileNumber},".");
    MatCompatible_Participant_IDs{1,FileNumber}=strcat('ID_',Participant_IDs{1,FileNumber});
    EEG.EpochRejections.ChannelsRemaining=EEG.nbchan;

    RELAX_epoching_cfg.ms_per_sample=(1000/EEG.srate); % v2.0.1 added to enable HEOG and VEOG step function based epoch rejections

    %% Interpolate channels that were excluded back into the data:
    if strcmp(RELAX_epoching_cfg.InterpolateRejectedChannels,'yes')
        EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    end

    % Record how many channels had to be interpolated back into the data after cleaning:
    EEG.EpochRejections.TotalChannels=EEG.nbchan;
    EEG.EpochRejections.ProportionOfChannelsInterpolated=(EEG.EpochRejections.TotalChannels-EEG.EpochRejections.ChannelsRemaining)/EEG.EpochRejections.TotalChannels;

    %% Epoch Data:
    
    % If resting data, you can use this to insert triggers at a specified
    % interval and then epoch data around triggers

    if strcmp(RELAX_epoching_cfg.DataType,'Resting')==1
        EEG=eeg_regepochs(EEG,'recurrence',RELAX_epoching_cfg.restingdatatriggerinterval,'eventtype','X','extractepochs','off');      
        EEG = pop_epoch( EEG, { 'X'}, RELAX_epoching_cfg.PeriodToEpoch, 'epochinfo', 'yes');
    end
     
    % If cognitive data with triggers embedded, this epochs data from -0.5 to 1 seconds and deletes triggers other than the epoched trigger within each epoch. Example below with an emotional Go Nogo task:  
    if strcmp(RELAX_epoching_cfg.DataType,'Task')
        
        EEG = pop_epoch( EEG, RELAX_epoching_cfg.TriggersToEpoch, RELAX_epoching_cfg.PeriodToEpoch, 'epochinfo', 'yes');
        if RELAX_epoching_cfg.Reject_HEOG_movements==1 || RELAX_epoching_cfg.Reject_VEOG_movements==1
            if isfield(EEG,'Non_cleaned_electrodes')
                % v2.0.1 NWB adjusted to enable epoching of HEOG electrodes in addition to the scalp electrodes:
                EEG.Non_cleaned_electrodes = pop_epoch(EEG.Non_cleaned_electrodes, RELAX_epoching_cfg.TriggersToEpoch, RELAX_epoching_cfg.PeriodToEpoch, 'epochinfo', 'yes');
            else
                error(['You have indicated you would like to exclude epochs containing horizontal eye movements, but the initial cleaning settings do not seem to have preserved these '...
                    'electrodes in the "EEG.Non_cleaned_electrodes" variable where they are expected. Specifying the HEOG electrodes in the "RELAX_cfg.electrodes_2_keep_but_not_clean" parameter '...
                    'during cleaning may help solve this issue']);
            end
        end
        
        if strcmp(RELAX_epoching_cfg.RemoveOtherTriggers,'yes')
            % Remove triggers that are in the epoch, but aren't the trigger
            % the data is being epoched around (this can be helpful if epoching
            % data separately by condition following this script):
            EEG = pop_selectevent( EEG, 'omitlatency', '-30001<=-1','type', RELAX_epoching_cfg.TriggersToEpoch, 'deleteevents','on');
            EEG = pop_selectevent( EEG, 'omitlatency', '1<=30001','type', RELAX_epoching_cfg.TriggersToEpoch, 'deleteevents','on');
            if RELAX_epoching_cfg.Reject_HEOG_movements==1
                EEG.Non_cleaned_electrodes = pop_selectevent( EEG.Non_cleaned_electrodes, 'omitlatency', '-30001<=-1','type', RELAX_epoching_cfg.TriggersToEpoch, 'deleteevents','on');
                EEG.Non_cleaned_electrodes = pop_selectevent( EEG.Non_cleaned_electrodes, 'omitlatency', '1<=30001','type', RELAX_epoching_cfg.TriggersToEpoch, 'deleteevents','on');
            end
        end
    end
    EEG = eeg_checkset( EEG );

    %% Check for horizontal eye movements larger than the threshold in microvolts using a step function 
    % (i.e. the amplitude of the moving average of the previous 100ms compared to the following 100ms)
    % described in Luck, S. J. (2014). An introduction to the event-related potential technique. MIT press.
    if RELAX_epoching_cfg.Reject_HEOG_movements==1
        search_onset=find(EEG.times==RELAX_epoching_cfg.HEOG_rejection_period(1,1)); 
        search_offset=find(EEG.times==RELAX_epoching_cfg.HEOG_rejection_period(1,2));

        % find each relevant electrode:
        for e=1:size(EEG.Non_cleaned_electrodes.chanlocs,2)
            if strcmpi(EEG.Non_cleaned_electrodes.chanlocs(e).labels,RELAX_epoching_cfg.HEOG_electrode_labels{1,1})
                HEOG_idx1=e;
            end
            if strcmpi(EEG.Non_cleaned_electrodes.chanlocs(e).labels,RELAX_epoching_cfg.HEOG_electrode_labels{1,2})
                HEOG_idx2=e;
            end
        end
        HEOG_diff=EEG.Non_cleaned_electrodes.data(HEOG_idx1,:,:)-EEG.Non_cleaned_electrodes.data(HEOG_idx2,:,:);   
        
        for t=search_onset:search_offset
            HEOG_degree_check(1,t,:)=trimmean(HEOG_diff(1,t-round(100/RELAX_epoching_cfg.ms_per_sample):t,:),95,2)-trimmean(HEOG_diff(1,t:t+round(100/RELAX_epoching_cfg.ms_per_sample),:),95,2); % using trimmean in case the electrodes are noisy
        end
        max_HEOG_check_in_epoch=squeeze(max(abs(HEOG_degree_check),[],[1 2]));
        EEG.HEOG_based_rejections=(max_HEOG_check_in_epoch>RELAX_epoching_cfg.HEOG_rejection_threshold)';
        EEG.EpochRejections.ProportionEpochsRejectedForHEOG=mean(EEG.HEOG_based_rejections);
    end

    if RELAX_epoching_cfg.Reject_VEOG_movements==1
        search_onset=find(EEG.times==RELAX_epoching_cfg.VEOG_rejection_period(1,1)); 
        search_offset=find(EEG.times==RELAX_epoching_cfg.VEOG_rejection_period(1,2));

        % find each relevant electrode:
        for e=1:size(EEG.Non_cleaned_electrodes.chanlocs,2)
            if strcmpi(EEG.Non_cleaned_electrodes.chanlocs(e).labels,RELAX_epoching_cfg.VEOG_electrode_labels{1,1})
                VEOG_idx1=e;
            end
            if size(RELAX_epoching_cfg.VEOG_electrode_labels,2)>1
                if strcmpi(EEG.Non_cleaned_electrodes.chanlocs(e).labels,RELAX_epoching_cfg.VEOG_electrode_labels{1,2})
                    VEOG_idx2=e;
                end
            end
        end
        if RELAX_epoching_cfg.compute_VEOG_difference==1
            VEOG=EEG.Non_cleaned_electrodes.data(VEOG_idx1,:,:)-EEG.Non_cleaned_electrodes.data(VEOG_idx2,:,:); 
        else
            VEOG=EEG.Non_cleaned_electrodes.data(VEOG_idx1,:,:);
        end

        %% Check for verticle eye movements larger than the threshold in microvolts using step function:
        for t=search_onset:search_offset
            VEOG_movement_check(1,t,:)=mean(VEOG(1,t-round(100/RELAX_epoching_cfg.ms_per_sample):t,:),2)-mean(VEOG(1,t:t+round(100/RELAX_epoching_cfg.ms_per_sample),:),2);
        end
        max_VEOG_check_in_epoch=squeeze(max(abs(VEOG_movement_check),[],[1 2]));
        EEG.VEOG_based_rejections=(max_VEOG_check_in_epoch>RELAX_epoching_cfg.VEOG_rejection_threshold)';
        EEG.EpochRejections.ProportionEpochsRejectedForVEOG=mean(EEG.VEOG_based_rejections);
    end

    %% Baseline Correct Data:

    % Regression based baseline correction method (recommended):
    if strcmp(RELAX_epoching_cfg.BL_correction_method,'regression')
        
        if RELAX_epoching_cfg.NumberOfFactors==2
            [EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,'Factor_1_Level_1', RELAX_epoching_cfg.BL_correction_Factor_1_Level_1, 'Factor_2_Level_1', RELAX_epoching_cfg.BL_correction_Factor_2_Level_1); 
            % if a 2 x 2 design, this will code triggers other than 'HappyGo'/'SadGo' as -1, and Go as 1 in the first factor, and triggers other than 'HappyGo'/'HappyNogo' as -1 in the second factor
        end
        if RELAX_epoching_cfg.NumberOfFactors==1
            [EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,'Factor_1_Level_1',RELAX_epoching_cfg.BL_correction_Factor_1_Level_1); % if a 2 condition design, this will code triggers other than 'Go' as -1, and Go as 1
        end
        if RELAX_epoching_cfg.NumberOfFactors==0
            [EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg); % if only 1 stimulus condition present for each participant
        end
        disp('regression BL correction applied');
    end

    % (the benefits of this method are explained in Alday, 2019, and the specific implementation performed in Bailey et al. 2022)
    % Alday, P. M. (2019). How much baseline correction do we need in ERP research? Extended GLM model can replace baseline correction while lifting its limits. Psychophysiology, 56(12), e13451.
    % Bailey et al. (2022). Meditators probably show increased behaviour-monitoring related neural activity. 

    % Note that for designs with more than 3 conditions for a single
    % factor, linear mixed modelling to baseline correct the data is more
    % appropriate (not provided here).

    % Also note that the method provided below is effective for 
    % nice clean data with a larger number of epochs per participant,
    % but for data that is not clean or very few epochs per participant (<5 is my guess) 
    % then a method that includes all participants in the single regression 
    % and includes individual as a factor is more appropriate.

    % Traditional subtraction based baseline correction (not recommended):
    if strcmp(RELAX_epoching_cfg.BL_correction_method,'subtraction')
        EEG = pop_rmbase( EEG, [RELAX_epoching_cfg.BLperiod]);
        disp('subtraction BL correction applied (not recommended)');
    end

    if strcmp(RELAX_epoching_cfg.BL_correction_method,'none')
        disp('No BL correction applied');
    end

    %% THIS SECTION CONTAINS FUNCTIONS WHICH DETECT AND REJECT ARTIFACTS

    % Count initial epochs:
    EEG.EpochRejections.InitialEpochs=size(EEG.data,3);

    % Any one of these functions can be commented out to ignore those artifacts
    % when creating the mask:

    % This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs:
    ROIidx= 1:EEG.nbchan; 
    EEG = pop_eegthresh(EEG,1,[ROIidx],[-RELAX_epoching_cfg.reject_amp],[RELAX_epoching_cfg.reject_amp],[EEG.xmin],[EEG.xmax],1,0);
    EEG = pop_jointprob(EEG,1,[ROIidx],RELAX_epoching_cfg.SingleChannelImprobableDataThreshold,RELAX_epoching_cfg.AllChannelImprobableDataThreshold,1,0);
    EEG = pop_rejkurt(EEG,1,(1:EEG.nbchan),RELAX_epoching_cfg.SingleChannelKurtosisThreshold,RELAX_epoching_cfg.AllChannelKurtosisThreshold,1,0);
    EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
    if RELAX_epoching_cfg.Reject_HEOG_movements==1
        EEG.reject.rejglobal(EEG.HEOG_based_rejections==1)=1; % add HEOG rejections
        epochs_to_reject=EEG.reject.rejglobal;
    end
    if RELAX_epoching_cfg.Reject_VEOG_movements==1
        EEG.reject.rejglobal(EEG.VEOG_based_rejections==1)=1; % add HEOG rejections
        epochs_to_reject=EEG.reject.rejglobal;
    end
    EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal] ,0);

    if RELAX_epoching_cfg.Reject_HEOG_movements==1 || RELAX_epoching_cfg.Reject_VEOG_movements==1
        EEG.Non_cleaned_electrodes = pop_rejepoch(EEG.Non_cleaned_electrodes, [epochs_to_reject] ,0);
    end
        
    %% If you have not filtered out data below 75Hz, you could use an objective muscle slope measure to reject epochs with remaining muscle activity:   
    % Use epoched data and FFT to detect slope of log frequency log
    % power, add periods exceeding muscle threshold to mask. This method is
    % designed for use with data that includes up to 75Hz, so  is not
    % useful if frequencies below 75Hz are filtered out

    if strcmp(RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity,'yes')
        [EEG] = RELAX_Rejecting_muscle_epochs(EEG, RELAX_epoching_cfg);
    end
    
    %% CHECK MEDIAN VOLTAGE SHIFT IN EACH EPOCH FOR EACH PARTICIPANT, AND RECORD ARTIFACT REJECTION DETAILS
    % The following measures the median size of the voltage shift in each epoch
    % for each participant. Used at the end of this script to provide
    % advice on which participants to manually check as potential bad
    % data:
    if strcmp(RELAX_epoching_cfg.InterpolateRejectedChannels,'yes')
        voltageshiftwithinepoch=range(EEG.data(:,:,:),2);
        Medianvoltageshiftwithinepoch(:,FileNumber)=median(voltageshiftwithinepoch,3);
    end
    EEG.EpochRejections.EpochsRemaining=size(EEG.data,3);
    EEG.EpochRejections.ProportionOfEpochsRejected=(EEG.EpochRejections.InitialEpochs-EEG.EpochRejections.EpochsRemaining)/EEG.EpochRejections.InitialEpochs;
    EEG.EpochRejections.aFileName=FileName;

     % Order the MWF Processing statistics structure in alphabetical order:
    [~, neworder] = sort(lower(fieldnames(EEG.EpochRejections)));
    EEG.EpochRejections = orderfields(EEG.EpochRejections, neworder);
    if isfield(EEG.EpochRejections, 'NaNsForExtremeOutlierPeriods')
        EEG.EpochRejections=rmfield(EEG.EpochRejections,'NaNsForExtremeOutlierPeriods');
    end
    EpochRejections(FileNumber,:) = struct2table(EEG.EpochRejections,'AsArray',true);
    
    EEG.RELAX_settings_used_to_epoch_this_file=RELAX_epoching_cfg;
    %% Save data:
    if strcmpi(RELAX_epoching_cfg.all_data_in_1_folder_or_BIDS_format,'BIDS')
        OutputPath=[RELAX_epoching_cfg.folders{FileNumber} filesep 'Epoched'];
        mkdir(OutputPath);
    else
        OutputPath=RELAX_epoching_cfg.OutputPath;
    end

    SaveSet_epoched =[OutputPath filesep FileName '_Epoched.set'];    
    EEG = pop_saveset( EEG, SaveSet_epoched );  
end

if strcmp(RELAX_epoching_cfg.InterpolateRejectedChannels,'yes')
    %% The following checks for participants who show outlying values for the median voltage shift within each epoch:
    % The following detects outlier files in the median amount of their max-min
    % voltage shift within an epoch, after adjusting for the fact that the data
    % across all participants is likely to be positively skewed with a log transform.
    MedianvoltageshiftwithinepochLogged=log10(Medianvoltageshiftwithinepoch);
    InterQuartileRange=iqr(MedianvoltageshiftwithinepochLogged,2);
    Upper25 = prctile(MedianvoltageshiftwithinepochLogged,75,2);
    Lower25 = prctile(MedianvoltageshiftwithinepochLogged,25,2);

    % 75th% and 25th% +/- (1.5 x IQR) is the recommended outlier detection method, 
    % so this is used to recommend which participants to manually check
    % However, I find this can be a bit too sensitive upon manual inspection,
    % and that 1.75, 2, or even 2.5 can be a better threshold.
    LowerBound=size(MedianvoltageshiftwithinepochLogged,1);
    UpperBound=size(MedianvoltageshiftwithinepochLogged,1);
    for x=1:size(MedianvoltageshiftwithinepochLogged,1)
        LowerBound(x,1)=Lower25(x,1)-(2*InterQuartileRange(x,1));
        UpperBound(x,1)=Upper25(x,1)+(2*InterQuartileRange(x,1));
    end

    VoltageShiftsTooLow=MedianvoltageshiftwithinepochLogged;
    VoltageShiftsTooLow=VoltageShiftsTooLow-LowerBound;
    VoltageShiftsTooLow(0<VoltageShiftsTooLow)=0;
    CumulativeSeverityOfAmplitudesBelowThreshold=sum(VoltageShiftsTooLow,1)';

    VoltageShiftsTooHigh=MedianvoltageshiftwithinepochLogged;
    VoltageShiftsTooHigh=VoltageShiftsTooHigh-UpperBound;
    VoltageShiftsTooHigh(0>VoltageShiftsTooHigh)=0;
    CumulativeSeverityOfAmplitudesAboveThreshold=sum(VoltageShiftsTooHigh,1)';

    % Plot:
    plot(LowerBound); hold on; plot(UpperBound); 
    hold on; plot(MedianvoltageshiftwithinepochLogged); xticks([1:1:60]);xticklabels({EEG.chanlocs.labels});legend('LowerBound', 'UpperBound');

    OutlierParticipantsToManuallyCheck = table(Participant_IDs', CumulativeSeverityOfAmplitudesBelowThreshold,CumulativeSeverityOfAmplitudesAboveThreshold);

    savefileone=[RELAX_epoching_cfg.OutputPath filesep 'OutlierParticipantsToManuallyCheck'];
    save(savefileone,'OutlierParticipantsToManuallyCheck')

    LoggedMedianVoltageShiftAcrossEpochs=array2table(MedianvoltageshiftwithinepochLogged);
    
    MissingLabels=find(cellfun(@isempty,MatCompatible_Participant_IDs));
    for x=1:size(MissingLabels,2)
        MatCompatible_Participant_IDs{1,MissingLabels(x)}='NotEpoched';
    end
    LoggedMedianVoltageShiftAcrossEpochs.Properties.VariableNames =MatCompatible_Participant_IDs';
    
    for c=1:size(EEG.chanlocs,2); chanlist{c}=EEG.chanlocs(c).labels; end
    LoggedMedianVoltageShiftAcrossEpochs.Properties.RowNames =chanlist;
    savefileone=[RELAX_epoching_cfg.OutputPath filesep 'LoggedMedianVoltageShiftAcrossEpochs'];
    save(savefileone,'LoggedMedianVoltageShiftAcrossEpochs')
else
    OutlierParticipantsToManuallyCheck='Not available when you have not interpolated electrodes back into the data';
end
 
savefileone=[RELAX_epoching_cfg.OutputPath filesep 'EpochRejections'];
save(savefileone,'EpochRejections')

savefileone=[RELAX_epoching_cfg.OutputPath filesep 'RELAX_epoching_cfg'];
save(savefileone,'RELAX_epoching_cfg')

if WarningAboutFileNumber==1
    warning('You instructed RELAX to epoch more files than were in your data folder. Check all your expected files were there?');
end

clearvars -except 'OutlierParticipantsToManuallyCheck' 'RELAX_epoching_cfg' 'EpochRejections' 'EpochRejectionStats' 'FilesWithoutConvergence' 'LoggedMedianVoltageShiftAcrossEpochs';
