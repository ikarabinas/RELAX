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
function [RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingRoundOneAllParticipants, RELAXProcessingRoundTwoAllParticipants, RELAXProcessing_wICA_AllParticipants,...
        RELAXProcessing_ICA_AllParticipants, RELAXProcessingRoundThreeAllParticipants, RELAX_issues_to_check, RELAX_issues_to_check_2nd_run, RELAXProcessingExtremeRejectionsAllParticipants] = RELAX_Wrapper (RELAX_cfg)

% Load pre-processing statistics file for these participants if it already
% exists (note that this can cause errors if the number of variables
% inserted into the output table differs between participants, which can be
% caused by using different parameters in the preceding section):

tic;

RELAX_cfg.OutputPath=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep];   % use fileseparators for increased compatability 
if ~exist(RELAX_cfg.OutputPath, 'dir'); mkdir(RELAX_cfg.OutputPath); end % make dir if not present

cd(RELAX_cfg.OutputPath);
dirList=dir('*.mat');
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatisticsRoundOne.mat')==1
        load('ProcessingStatisticsRoundOne.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatisticsRoundTwo.mat')==1
        load('ProcessingStatisticsRoundTwo.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatisticsRoundThree.mat')==1
        load('ProcessingStatisticsRoundThree.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'RawMetrics.mat')==1
        load('RawMetrics.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'CleanedMetrics.mat')==1
        load('CleanedMetrics.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatistics_wICA.mat')==1
        load('ProcessingStatistics_wICA.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatistics_ICA.mat')==1
        load('ProcessingStatistics_ICA.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'RELAX_issues_to_check.mat')==1
        load('RELAX_issues_to_check.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'RELAXProcessingExtremeRejectionsAllParticipants.mat')==1
        load('RELAXProcessingExtremeRejectionsAllParticipants.mat');
    end
end

if ~isempty(RELAX_cfg.filename)
    RELAX_cfg.FilesToProcess = 1;
    RELAX_cfg.SingleFile = 1; % 1 for single file
else
    RELAX_cfg.SingleFile = 0; % 0 for multiple files
end

WarningAboutFileNumber=0;
if size(RELAX_cfg.FilesToProcess,2) > size(RELAX_cfg.files,2)
    RELAX_cfg.FilesToProcess=RELAX_cfg.FilesToProcess(1,1):size(RELAX_cfg.files,2);
    WarningAboutFileNumber=1;
end

%% Loop selected files in the directory list:
for FileNumber=RELAX_cfg.FilesToProcess(1,1:size(RELAX_cfg.FilesToProcess,2))
    
    if RELAX_cfg.SingleFile == 0
        RELAX_cfg.filename=RELAX_cfg.files{FileNumber};
    end

    clearvars -except 'RELAX_cfg' 'FileNumber' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingRoundOneAllParticipants' 'RELAXProcessingRoundTwoAllParticipants' 'RELAXProcessing_wICA_AllParticipants'...
        'RELAXProcessing_ICA_AllParticipants' 'RELAXProcessingRoundThreeAllParticipants' 'Warning' 'RELAX_issues_to_check' 'RELAX_issues_to_check_2nd_run'...
        'RELAXProcessingExtremeRejectionsAllParticipants' 'WarningAboutFileNumber';
    %% Load data (assuming the data is in EEGLAB .set format):

    %  1.1.4: fix error where PREP seems to be removed from the path after an
    % EEGLAB update:
    PrepFileLocation = which('pop_prepPipeline','-all');
    PrepFolderLocation=extractBefore(PrepFileLocation,'pop_prepPipeline.m');

    cd(RELAX_cfg.myPath);
    EEG = pop_loadset(RELAX_cfg.filename);

    FileName = extractBefore(RELAX_cfg.filename,".");
    if RELAX_cfg.SingleFile == 1 % RELAX v1.1.3 NWB added to stop RELAX crashing when trying to save due to whole folder being included twice in save file
        last_slash_pos = find(RELAX_cfg.filename == '\', 1, 'last');
        FileName = extractBetween(RELAX_cfg.filename,last_slash_pos+1,".");
        FileName = FileName{1};
    end

    EEG.RELAXProcessing.aFileName=cellstr(FileName);
    EEG.RELAXProcessingExtremeRejections.aFileName=cellstr(FileName);
    
    EEG.RELAX.Data_has_been_averagerereferenced=0;
    EEG.RELAX.Data_has_been_cleaned=0;
    RELAX_cfg.ms_per_sample=(1000/EEG.srate);

    savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_cfg'];
    save(savefileone,'RELAX_cfg')

    %% Select channels 
    if ~isempty(RELAX_cfg.caploc)
        EEG=pop_chanedit(EEG,  'lookup', RELAX_cfg.caploc);
    end
    
    %% Delete channels that are not relevant if present       
    EEG=pop_select(EEG,'nochannel',RELAX_cfg.ElectrodesToDelete);
    EEG = eeg_checkset( EEG );
    EEG.allchan=EEG.chanlocs; % take list of all included channels before any rejections

    %% Band Pass filter data: 
    if strcmp(RELAX_cfg.NotchFilterType,'Butterworth')
        % Use TESA to apply butterworth filter: 
        EEG = RELAX_filtbutter( EEG, RELAX_cfg.LineNoiseFrequency-3, RELAX_cfg.LineNoiseFrequency+3, 4, 'bandstop','acausal');
    end

    if strcmp(RELAX_cfg.LowPassFilterBeforeMWF,'no') % updated implementation, avoiding low pass filtering prior to MWF reduces chances of rank deficiencies, increasing potential values for MWF delay period 
        if strcmp(RELAX_cfg.FilterType,'Butterworth')
            EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, [], 4, 'highpass', RELAX_cfg.causal_or_acausal_filter);
        end
        if strcmp(RELAX_cfg.FilterType,'pop_eegfiltnew')
            EEG = pop_eegfiltnew(EEG,RELAX_cfg.HighPassFilter,[]);
        end
    end
    
    if strcmp(RELAX_cfg.LowPassFilterBeforeMWF,'yes') % original implementation, not recommended unless downsampling, as increases chances of rank deficiencies
        EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, RELAX_cfg.LowPassFilter, 4, 'bandpass', RELAX_cfg.causal_or_acausal_filter);
    end

    if strcmp(RELAX_cfg.DownSample,'yes')
        EEG = pop_resample(EEG,RELAX_cfg.DownSample_to_X_Hz); % downsample data (if applied, should always be applied after low pass filtering)
        RELAX_cfg.ms_per_sample=(1000/EEG.srate);
    end

    if RELAX_cfg.ms_per_sample<0.7
        warning('The sampling rate for this file is quite high. Depending on your processing power, RELAX may run slowly or even stall, especially if applying MWF cleaning');
        warning('RELAX was validated using 1000Hz sampling rates, which is still a high sample rate for most analyses. You could downsample your data by setting the relevant options in RELAX');
    end

    if strcmp(RELAX_cfg.NotchFilterType,'ZaplinePlus')
        [EEG ] = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('plotResults',0)); % requires the zapline plus plugin. Best applied after downsampling to 250Hz or 500Hz.
    end

    %% Clean flat channels and bad channels showing improbable data:
    % PREP pipeline: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
    addpath(genpath(PrepFolderLocation{1,1})); %  1.1.4: fix error where PREP seems to be removed from the path after an EEGLAB update
    noisyOut = findNoisyChannels(EEG);  
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject={};
    for x=1:size(noisyOut.noisyChannels.all,2) % loop through output of PREP's findNoisyChannels and take a record of noisy electrodes for deletion:
        PREPBasedChannelToReject{x}=EEG.chanlocs(noisyOut.noisyChannels.all(x)).labels;
        EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject = PREPBasedChannelToReject';
    end
    EEG=pop_select(EEG,'nochannel',noisyOut.noisyChannels.all); % delete noisy electrodes detected by PREP

    continuousEEG=EEG;

    [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(continuousEEG, RELAX_cfg); % Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
    [continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg); % Mark extreme periods for exclusion from MWF cleaning, and deletion before wICA cleaning

    % Use the continuous data to detect eye blinks and mark
    % these in the EEG.event as well as in the mask. The output is
    % continuous data but includes all the previous extreme period 
    % markings from the epoched data.
    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg); % use an IQR threshold method to detect and mark blinks
        if continuousEEG.RELAX.IQRmethodDetectedBlinks(1,1)==0 % If a participants doesn't show any blinks, make a note
            NoBlinksDetected{FileNumber,1}=FileName; 
            warning('No blinks were detected - if blinks are expected then you should visually inspect the file');
        end
        if RELAX_cfg.computerawmetrics==1
        [continuousEEG, epochedEEG] = RELAX_metrics_blinks(continuousEEG, epochedEEG); % record blink amplitude ratio from raw data for comparison.
        end
    end

    % Record extreme artifact rejection details for all participants in single table:
    RELAXProcessingExtremeRejectionsAllParticipants(FileNumber,:) = struct2table(epochedEEG.RELAXProcessingExtremeRejections,'AsArray',true);

    rawEEG=continuousEEG; % Take a copy of the not yet cleaned data for calculation of all cleaning SER and ARR at the end
    
    %% Mark artifacts for calculating SER and ARR, regardless of whether MWF is performed (RELAX v1.1.3 update): 
    if RELAX_cfg.computecleanedmetrics==1 && (RELAX_cfg.Do_MWF_Once==0 || RELAX_cfg.Do_MWF_Twice==0 || RELAX_cfg.Do_MWF_Thrice==0)
        [Marking_artifacts_for_SER_ARR, ~] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg); 
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
        [Marking_artifacts_for_SER_ARR] = RELAX_horizontaleye(continuousEEG, RELAX_cfg); 
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
        [Marking_artifacts_for_SER_ARR, ~] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg); % Use epoched data to add periods showing excessive drift to the mask 
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
        Marking_all_artifacts_for_SER_ARR.RELAX.NaNsForExtremeOutlierPeriods=continuousEEG.RELAX.NaNsForExtremeOutlierPeriods; 
        [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, RELAX_cfg, 'notblinks'); % If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out. 
        if isfield(continuousEEG.RELAX,'eyeblinkmask')
                Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(continuousEEG.RELAX.eyeblinkmask==1)=1; 
                [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, RELAX_cfg, 'blinks'); 
        end
        continuousEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength; 
        rawEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength; 
    end
    
    if RELAX_cfg.saveextremesrejected==1
        if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'], 'dir')
            mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'])
        end
        SaveSetExtremes_Rejected =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Extremes_Rejected', filesep FileName '_Extremes_Rejected.set'];    
        EEG = pop_saveset( rawEEG, SaveSetExtremes_Rejected ); % If desired, save data here with bad channels deleted, filtering applied, extreme outlying data periods marked
    end

    %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK AND CLEAN MUSCLE ARTIFACTS
    % Any one of these functions can be commented out to ignore those artifacts
    % when creating the mask    
    if RELAX_cfg.Do_MWF_Once==1

        % Use epoched data and FFT to detect slope of log frequency log
        % power, add periods exceeding muscle threshold to mask:
        [continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg);  
        if RELAX_cfg.computerawmetrics==1
            [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg); % record muscle contamination metrics from raw data for comparison.
        end

        EEG=continuousEEG; % Return continuousEEG to the "EEG" variable for MWF processing

        % If including eye blink cleaning in first round MWF, then insert
        % eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==1
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
        end

        % The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can
        % lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks'); % If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out.
        
        EEG.RELAX.NoiseMaskFullLengthR1=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR1=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 
  
        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        RELAX_cfg.MWFDelayPeriod=RELAX_cfg.MWFDelayPeriod_for_muscle_artifacts; 
        RELAX_cfg.MWF_delay_spacing=RELAX_cfg.MWF_delay_spacing_for_muscle_artifacts;
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);          

        EEG.RELAXProcessingRoundOne=EEG.RELAXProcessing; % Record MWF cleaning details from round 1 in EEG file          
        RELAXProcessingRoundOne=EEG.RELAXProcessingRoundOne; % Record MWF cleaning details from round 1 into file for all participants
        
        if isfield(RELAXProcessingRoundOne,'Details')
            RELAXProcessingRoundOne=rmfield(RELAXProcessingRoundOne,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundOne,'Details')
                EEG.RELAXProcessingRoundOne=rmfield(EEG.RELAXProcessingRoundOne,'Details');
            end
        end
        
        % Record processing statistics for all participants in single table:
        RELAXProcessingRoundOneAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundOne,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        % Save round 1 MWF pre-processing:
        if RELAX_cfg.saveround1==1
            if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF'], 'dir')
                mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF'])
            end
            SaveSetMWF1 =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF', filesep FileName '_MWF1.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF1 ); 
        end
    end

    %% PERFORM A SECOND ROUND OF MWF. THIS IS HELPFUL IF THE FIRST ROUND DOESN'T SUFFICIENTLY CLEAN ARTIFACTS. 

    % This has been suggested to be useful by Somers et al (2018)
    % (particularly when used in a cascading fashion). 

    % However, I can see risks. If artifact masks fall on task relevant
    % activity in both rounds of the MWF, it may be that the task relevant data
    % is just cleaned right out of the signal.
    
    if RELAX_cfg.Do_MWF_Twice==1

        EEG.RELAXProcessing.aFileName=cellstr(FileName);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        
        % If blinks weren't initially detected because they were 
        % disguised by the the muscle artifact, detect them here
        % (this happens in <1/200 cases, but is a good back up).
        if RELAX_cfg.ProbabilityDataHasNoBlinks==0
            if EEG.RELAX.IQRmethodDetectedBlinks(1,1)==0
                continuousEEG=EEG;
                [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg);
                EEG=continuousEEG;
            end
        end
        
        % If including eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask:
        if isfield(EEG.RELAX, 'eyeblinkmask')
            if RELAX_cfg.MWFRoundToCleanBlinks==2
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
            end
        end
  
        % The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can
        % lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'blinks');
        
        EEG.RELAX.NoiseMaskFullLengthR2=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        RELAX_cfg.MWFDelayPeriod=RELAX_cfg.MWFDelayPeriod_for_eye_movements; 
        RELAX_cfg.MWF_delay_spacing=RELAX_cfg.MWF_delay_spacing_for_eye_movements; % set how sparsely the delay stacking is spread
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);           
 
        EEG.RELAXProcessingRoundTwo=EEG.RELAXProcessing; % Record MWF cleaning details from round 2 in EEG file
        RELAXProcessingRoundTwo=EEG.RELAXProcessingRoundTwo; % Record MWF cleaning details from round 2 into file for all participants
        if isfield(RELAXProcessingRoundTwo,'Details')
            RELAXProcessingRoundTwo=rmfield(RELAXProcessingRoundTwo,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundTwo,'Details')
                EEG.RELAXProcessingRoundTwo=rmfield(EEG.RELAXProcessingRoundTwo,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingRoundTwoAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundTwo,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        % Save round 2 MWF pre-processing:
        if RELAX_cfg.saveround2==1
            if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF'], 'dir')
                mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF'])
            end
            SaveSetMWF2 =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF', filesep FileName '_MWF2.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF2 ); 
        end     
    end
    
    %% PERFORM A THIRD ROUND OF MWF.    
    if RELAX_cfg.Do_MWF_Thrice==1

        EEG.RELAXProcessing.aFileName=cellstr(FileName);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        % If less than 5% of data was masked as eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask in round 3:
        if isfield(EEG.RELAX,'ProportionMarkedInMWFArtifactMaskTotalR2') % NWB added to make sure function doesn't bug when trying to check this variable if it doesn't exist
            if EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2<0.05
                if isfield(EEG.RELAX, 'eyeblinkmask')
                    EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                    EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                    EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
                end
            end
        end

        % Epoch the data into 1 second epochs with a 500ms overlap. Outputs
        % both the ContinuousEEG (which has been filtered above by this
        % point) and the epoched data as EEG.
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
        
        %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK ARTIFACTS

        [continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg); % Use epoched data to add periods showing excessive drift to the mask
        
        % Use the filtered continuous data to detect horizontal eye
        % movements and mark these in the EEG.event as well as in the mask.
        % You may want to simply reject horizontal eye movements at a later
        % stage if your task requires participants to look straight ahead
        % for the entire task. Alternatively, if your task requires
        % participants to complete horizontal eye movements time locked to
        % a stimuli, this section will mark every event with these
        % horizontal eye movements as an artifact, and should not be
        % implemented.
        
        % The output is continuous data:
        [continuousEEG] = RELAX_horizontaleye(continuousEEG, RELAX_cfg);
    
        %% Return to the "EEG" variable for MWF processing:
        EEG=continuousEEG;
        
        % If including eye blink cleaning in third round MWF, then insert
        % eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==3
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
        end
 
        % The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can
        % lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks');
        
        EEG.RELAX.NoiseMaskFullLengthR3=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR3=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        RELAX_cfg.MWFDelayPeriod=RELAX_cfg.MWFDelayPeriod_for_eye_movements; 
        RELAX_cfg.MWF_delay_spacing=RELAX_cfg.MWF_delay_spacing_for_eye_movements; % set how sparsely the delay stacking is spread
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);               
        
        if isfield(EEG.RELAX, 'eyeblinkmask') % if eyeblinkmask has been created, do the following (thanks to Jane Tan for the suggested bug fix when eyeblinkmask is not created)
            EEG.RELAX=rmfield(EEG.RELAX,'eyeblinkmask'); % remove variables that are no longer necessary
        end
        
        EEG.RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 in EEG file
        RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 into file for all participants
        
        if isfield(RELAXProcessingRoundThree,'Details')
            RELAXProcessingRoundThree=rmfield(RELAXProcessingRoundThree,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundThree,'Details')
                EEG.RELAXProcessingRoundThree=rmfield(EEG.RELAXProcessingRoundThree,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingRoundThreeAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundThree,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');

        if RELAX_cfg.saveround3==1
            if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '3xMWF'], 'dir')
                mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '3xMWF'])
            end
            SaveSetMWF3 =[RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep '3xMWF', filesep FileName '_MWF3.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF3 ); 
        end         
    end
    
    %% Perform robust average re-referencing of the data, reject periods marked as extreme outliers    
    if RELAX_cfg.Do_MWF_Once==0
        EEG=continuousEEG;
    end
          
    % Reject periods that were marked as NaNs in the MWF masks because they 
    % showed extreme shift within the epoch or extremely improbable data:
    EEG = eeg_eegrej( EEG, EEG.RELAX.ExtremelyBadPeriodsForDeletion);
    
    if strcmp(RELAX_cfg.LowPassFilterBeforeMWF,'no') % if low pass filtering wasn't applied before MWF cleaning (recommended) apply it here
        if strcmp(RELAX_cfg.FilterType,'Butterworth')
            EEG = RELAX_filtbutter( EEG, [], RELAX_cfg.LowPassFilter, 4, 'lowpass', RELAX_cfg.causal_or_acausal_filter);
        end
        if strcmp(RELAX_cfg.FilterType,'pop_eegfiltnew')
            EEG = pop_eegfiltnew(EEG,[],RELAX_cfg.LowPassFilter);
        end
    end
    
    [EEG] = RELAX_average_rereference(EEG);
    EEG = eeg_checkset( EEG );  

    %% Perform wICA on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_targeted_wICA==1
        % The following cleans eye movements and muscle artifacts in the
        % independent component space by a combination of wavelet enhanced
        % ICA cleaning and targeting to restrict the cleaning to only 
        % artifact periods for eye movement components, and high pass 
        % filtering muscle components at 15Hz:
        EEG.RELAXProcessing_wICA.aFileName=cellstr(FileName);
        [EEG] = RELAX_targeted_wICA(EEG,RELAX_cfg);
        % setting 'RELAX_cfg.Report_all_wICA_info' to 1 will report proportion of ICs categorized as each category, and variance explained by ICs from each category (function is ~20s slower if this is implemented)
        EEG = eeg_checkset( EEG );
        RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_wICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    end
    
    %% Perform wICA on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        % The following performs wICA, implemented on only the components
        % marked as artifact by ICLabel.
        EEG.RELAXProcessing_wICA.aFileName=cellstr(FileName);
        [EEG,~, ~, ~, ~] = RELAX_wICA_on_ICLabel_artifacts(EEG,RELAX_cfg.ICA_method, 1, 0, EEG.srate, 5,'coif5',RELAX_cfg.Report_all_ICA_info,RELAX_cfg.ICLabel_thresholds,RELAX_cfg.Clean_other_comps); 
        % setting 'RELAX_cfg.Report_all_wICA_info' to 1 will report proportion of ICs categorized as each category, and variance explained by ICs from each category (function is ~20s slower if this is implemented)
        EEG = eeg_checkset( EEG );
        RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_wICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    end
    
    %% Perform ICA subtract on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_ICA_subtract==1
        % The following performs ICA sutraction, implemented on only the components
        % marked as artifact by ICLabel.
        EEG.RELAXProcessing_ICA.aFileName=cellstr(FileName);
        EEG = RELAX_ICA_subtract(EEG,RELAX_cfg);
        EEG = eeg_checkset( EEG );
        RELAXProcessing_ICA=EEG.RELAXProcessing_ICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_ICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_ICA,'AsArray',true);
    end
    
    EEG.RELAX.Data_has_been_cleaned=1;
    
    %% COMPUTE CLEANED METRICS:
    if RELAX_cfg.computecleanedmetrics==1    
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
        [continuousEEG, ~] = RELAX_metrics_blinks(continuousEEG, epochedEEG);
        [continuousEEG, ~] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg);

        [continuousEEG] = RELAX_metrics_final_SER_and_ARR(rawEEG, continuousEEG); % this is only a good metric for testing only the cleaning of artifacts marked for cleaning by MWF, see notes in function.

        EEG=continuousEEG;
        EEG = rmfield(EEG,'RELAXProcessing');

        if isfield(EEG,'RELAX_Metrics')
            if isfield(EEG.RELAX_Metrics, 'Cleaned')
                if isfield(EEG.RELAX_Metrics.Cleaned,'BlinkAmplitudeRatio')
                    CleanedMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio,1),FileNumber)=EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio;
                    CleanedMetrics.BlinkAmplitudeRatio(CleanedMetrics.BlinkAmplitudeRatio==0)=NaN;
                end
                if isfield(EEG.RELAX_Metrics.Cleaned,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                    CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(FileNumber)=EEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                    CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(FileNumber)=EEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
                end
                if isfield(EEG.RELAX_Metrics.Cleaned,'All_SER')
                    CleanedMetrics.All_SER(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_SER;
                    CleanedMetrics.All_ARR(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_ARR;
                end
            end
            if isfield(EEG.RELAX_Metrics, 'Raw')
                if isfield(EEG.RELAX_Metrics.Raw,'BlinkAmplitudeRatio')
                    RawMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio,1),FileNumber)=EEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio;
                    RawMetrics.BlinkAmplitudeRatio(RawMetrics.BlinkAmplitudeRatio==0)=NaN;
                end
                if isfield(EEG.RELAX_Metrics.Raw,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                    RawMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(FileNumber)=EEG.RELAX_Metrics.Raw.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                    RawMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(FileNumber)=EEG.RELAX_Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
                end
            end   
        end
    end

    %% Record warnings about potential issues:
    EEG.RELAX_issues_to_check.aFileName=cellstr(FileName);
    if size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)>RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1); % 1.1.4: fix dimension specification error
    else
        EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=0;
    end
    if (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1))...
            >=RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=...
            (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1));
    else
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=0;
    end
    if EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier>0.20
        EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier=EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier;
    else 
        EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier=0;
    end
    if isfield(EEG.RELAX, 'IQRmethodDetectedBlinks') % if IQRmethodDetectedBlinks has been created, do the following (thanks to Jane Tan for the suggested bug fix when IQRmethodDetectedBlinks is not created)
        EEG.RELAX_issues_to_check.NoBlinksDetected=(EEG.RELAX.IQRmethodDetectedBlinks==0);
    end
    if RELAX_cfg.Do_MWF_Once==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R1=isa(EEG.RELAXProcessingRoundOne.RankDeficiency,'char');
    end
    if RELAX_cfg.Do_MWF_Twice==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R2=isa(EEG.RELAXProcessingRoundTwo.RankDeficiency,'char');
    end
    if RELAX_cfg.Do_MWF_Thrice==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R3=isa(EEG.RELAXProcessingRoundThree.RankDeficiency,'char');
    end
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        if EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA>0.80
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA;
        else
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
        end
        EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA;
        EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_wICA.fastica_symm_Didnt_Converge(1,3);
    end
    if RELAX_cfg.Perform_ICA_subtract==1
        if EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA>0.80
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA;
        else
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
        end
        EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA;
        EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_ICA.fastica_symm_Didnt_Converge(1,3);
    end

    if strcmp(RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning,'yes')
        EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    end
    
    %% SAVE FILE:
    if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Cleaned_Data'], 'dir')
        mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Cleaned_Data'])
    end
    SaveSet_CleanedFile =[RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep 'Cleaned_Data', filesep FileName '_RELAX.set'];  
    EEG.RELAX_settings_used_to_clean_this_file=RELAX_cfg;
    EEG = pop_saveset( EEG, SaveSet_CleanedFile ); 
    
    % Record warnings for all participants in single table:
    try
        RELAX_issues_to_check(FileNumber,:) = struct2table(EEG.RELAX_issues_to_check,'AsArray',true);
    catch
        RELAX_issues_to_check_2nd_run(FileNumber,:) = struct2table(EEG.RELAX_issues_to_check,'AsArray',true);
        warning('The variable: "RELAX_issues_to_check" already exists and includes different settings from your current settings')
        warning('This is likely from a previous run of RELAX. Saving variable as "RELAX_issues_to_check_2nd_run" instead');
    end
    
    %% Save statistics for each participant and across participants, graph cleaning metrics:
    
    % Also set empty output variables in case these are not produced because certain
    % parameters have been switched off:

    savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAXProcessingExtremeRejectionsAllParticipants'];
    save(savefileone,'RELAXProcessingExtremeRejectionsAllParticipants')
    if RELAX_cfg.Do_MWF_Once==1
        savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundOne'];
        save(savefileone,'RELAXProcessingRoundOneAllParticipants')
    else
        RELAXProcessingRoundOneAllParticipants={};
    end
    if RELAX_cfg.Do_MWF_Twice==1
        savefiletwo=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundTwo'];
        save(savefiletwo,'RELAXProcessingRoundTwoAllParticipants')
    else
        RELAXProcessingRoundTwoAllParticipants={};
    end
    if RELAX_cfg.Do_MWF_Thrice==1
        savefilethree=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundThree'];
        save(savefilethree,'RELAXProcessingRoundThreeAllParticipants')
    else
        RELAXProcessingRoundThreeAllParticipants={};
    end
    if RELAX_cfg.Perform_wICA_on_ICLabel==1 || RELAX_cfg.Perform_targeted_wICA==1
        savefilefour=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_wICA'];
        save(savefilefour,'RELAXProcessing_wICA_AllParticipants')
    else
        RELAXProcessing_wICA_AllParticipants={}; 
    end
    if RELAX_cfg.Perform_ICA_subtract==1
        savefilefour=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_ICA'];
        save(savefilefour,'RELAXProcessing_ICA_AllParticipants')
    else
        RELAXProcessing_ICA_AllParticipants={}; 
    end
    if exist('CleanedMetrics','var')
        savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'CleanedMetrics'];
        save(savemetrics,'CleanedMetrics')
    else
        CleanedMetrics={};
    end
    if exist('RawMetrics','var')
        savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RawMetrics'];
        save(savemetrics,'RawMetrics')
    else
        RawMetrics={};
    end
    if exist('RELAX_issues_to_check','var')
        savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_issues_to_check'];
        save(savemetrics,'RELAX_issues_to_check')
    end
    if exist('RELAX_issues_to_check_2nd_run','var')
        savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_issues_to_check_2nd_run'];
        save(savemetrics,'RELAX_issues_to_check_2nd_run')
    end
    RELAX_cfg.filename=[];
    savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_cfg'];
    save(savefileone,'RELAX_cfg')    
end

set(groot, 'defaultAxesTickLabelInterpreter','none');
if RELAX_cfg.computecleanedmetrics==1
    try
        figure('Name','BlinkAmplitudeRatio','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        boxplot(CleanedMetrics.BlinkAmplitudeRatio);
        xticklabels(RELAX_cfg.files); xtickangle(90);
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
    catch
    end
    try
        figure('Name','MeanMuscleStrengthFromOnlySuperThresholdValues','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        b=bar(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues); 
        xtickangle(90); xticks([1:1:size(RELAX_cfg.files,2)]); b(1).BaseValue = RELAX_cfg.MuscleSlopeThreshold;
        xticklabels(RELAX_cfg.files); ylim([RELAX_cfg.MuscleSlopeThreshold max(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues)+1]);b.ShowBaseLine='off';
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
    catch
    end
    try
        figure('Name','ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        bar(CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel);
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        xtickangle(90); xticks([1:1:size(RELAX_cfg.files,2)]);
        xticklabels(RELAX_cfg.files);
    catch
    end
end

clearvars -except 'RELAX_cfg' 'FileNumber' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingRoundOneAllParticipants' 'RELAXProcessingRoundTwoAllParticipants' 'RELAXProcessing_wICA_AllParticipants'...
        'RELAXProcessing_ICA_AllParticipants' 'RELAXProcessingRoundThreeAllParticipants' 'Warning' 'RELAX_issues_to_check' 'RELAX_issues_to_check_2nd_run'...
        'RELAXProcessingExtremeRejectionsAllParticipants' 'WarningAboutFileNumber';

if ~exist('RELAX_issues_to_check_2nd_run','var')    
    warning('Check "RELAX_issues_to_check" to see if any issues were noted for specific files');
    RELAX_issues_to_check_2nd_run=['This variable is only here as a placeholder in case you have already run RELAX once,...' ...
        ' and had previously left the output variables in the same folder as the folder where you have saved the currently cleaned data.' ...
        'This was not an issue with the current run, so this placeholder variable was not filled'];
elseif exist('RELAX_issues_to_check_2nd_run','var')
    warning('Check "RELAX_issues_to_check_2nd_run" to see if any issues were noted for specific files');
end

if WarningAboutFileNumber==1
    warning('You instructed RELAX to clean more files than were in your data folder. Check all your expected files were there?');
end

if RELAX_cfg.ProbabilityDataHasNoBlinks<2 && sum(RELAX_issues_to_check.NoBlinksDetected)>1
    f = msgbox('RELAX did not detect any blinks for some files. Open the "RELAX_issues_to_check" struct in the workspace to check which files. We recommend visually inspecting these files to ensure there has not been an error.'...
    ,'No blinks detected for some files');    
    set(f,'Position',[500,500,450,100]);
    ah = get( f, 'CurrentAxes' );
    ch = get( ah, 'Children' );
    set( ch, 'FontSize', 12 ); %makes text bigger
end

if find(RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold>0)>0
    f = msgbox('Some files met or exceeded the electrode rejection threshold. We recommend visually inspecting the raw and cleaned files where this is the case. Open the "RELAX_issues_to_check" struct in the workspace, and check the third column. Files that exceeded the threshold will show a value above 0. Exclude files where raw data seems irretrievably noisy, or cleaned data still contains excessive noise.'...
    ,'Some files met or exceeded the electrode rejection threshold');    
    set(f,'Position',[300,300,450,150]);
    ah = get( f, 'CurrentAxes' );
    ch = get( ah, 'Children' );
    set( ch, 'FontSize', 12 ); %makes text bigger
end

toc
%% POTENTIAL IMPROVEMENTS THAT COULD BE MADE:
% 1) work out a way to threshold horizontal eye movements so the script
% catches the onset and offset, rather than just the absolute +/-2MAD on
% opposite sides of the head?
% 2) adding a requirement that the IQR blink detection method detects that
% positive amplitude shifts are biased towards frontal electrodes?
