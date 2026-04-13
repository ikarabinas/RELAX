%% Recover RELAX and GEDAI cleaning metrics from individual files
% The RELAX pipeline normally aggregates these metrics automatically, but
% this script can recover them in the case of an interrupted run and saves
% metrics under "Cleaned_Metrics_recovered.mat" in the specified data dir.

datadir = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_twICA_GEDAI-dlpfc/';

% Point to RELAX-created cleaned data folder and define savepath
cleaned_path = fullfile(datadir, 'RELAXProcessed', 'Cleaned_Data');
savepath = fullfile(datadir, 'RELAXProcessed');

% Get list of all cleaned files
all_cleaned_files = dir(fullfile(cleaned_path, '*RELAX*.set'));
fprintf('Total cleaned files found: %d\n', length(all_cleaned_files));

% Initialize struct
CleanedMetrics_recovered = struct();

for f = 1:length(all_cleaned_files)
    fprintf('Loading file %d/%d: %s\n', f, length(all_cleaned_files), all_cleaned_files(f).name);
    EEG = pop_loadset('filename', all_cleaned_files(f).name, 'filepath', cleaned_path);
    
    % Filename for unambiguous tracking
    [~, basename, ~] = fileparts(all_cleaned_files(f).name);
    CleanedMetrics_recovered(f).aFileName = basename;
    
    % RELAX_Metrics.Cleaned fields
    CleanedMetrics_recovered(f).All_SER                                                 = EEG.RELAX_Metrics.Cleaned.All_SER;
    CleanedMetrics_recovered(f).All_ARR                                                 = EEG.RELAX_Metrics.Cleaned.All_ARR;
    CleanedMetrics_recovered(f).MeanMuscleStrengthFromOnlySuperThresholdValues          = EEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues;
    CleanedMetrics_recovered(f).ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel = EEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
    
    % BlinkAmplitudeRatio may be missing for files where few or no blinks were detected
    if isfield(EEG.RELAX_Metrics.Cleaned, 'BlinkAmplitudeRatio')
        frontal_idx = [18, 19, 25, 26, 31, 32, 33, 37];
        CleanedMetrics_recovered(f).BlinkAmplitudeRatio = EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio;
        CleanedMetrics_recovered(f).MeanBAR = mean(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio);
        CleanedMetrics_recovered(f).fBAR = mean(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio(frontal_idx));
    else
        CleanedMetrics_recovered(f).BlinkAmplitudeRatio = NaN;
        CleanedMetrics_recovered(f).MeanBAR = NaN;
        CleanedMetrics_recovered(f).fBAR = NaN;
        fprintf('    WARNING: BlinkAmplitudeRatio missing in %s\n', basename);
    end
    fprintf('RELAX CleanedMetrics saved to file.\n');
    
    % GEDAI fields
    if isfield(EEG.etc, 'GEDAI')
        CleanedMetrics_recovered(f).GEDAI_SENSAI_score                = EEG.etc.GEDAI.SENSAI_score;
        CleanedMetrics_recovered(f).GEDAI_SENSAI_score_per_band       = EEG.etc.GEDAI.SENSAI_score_per_band;
        CleanedMetrics_recovered(f).GEDAI_artifact_threshold_per_band = EEG.etc.GEDAI.artifact_threshold_per_band;
        CleanedMetrics_recovered(f).GEDAI_mean_ENOVA                  = EEG.etc.GEDAI.mean_ENOVA;
        CleanedMetrics_recovered(f).GEDAI_ENOVA_per_band              = EEG.etc.GEDAI.ENOVA_per_band;
        CleanedMetrics_recovered(f).GEDAI_ENOVA_per_epoch             = EEG.etc.GEDAI.ENOVA_per_epoch;
        CleanedMetrics_recovered(f).GEDAI_total_epochs                = EEG.etc.GEDAI.total_epochs;
        CleanedMetrics_recovered(f).GEDAI_percent_rejected            = EEG.etc.GEDAI.percentage_rejected;
        
        fprintf('GEDAI metrics saved to file.\n');
    end
    
    % RELAX_issues_to_check fields
    CleanedMetrics_recovered(f).PREP_rejected_too_many_electrodes                       = EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes;
    CleanedMetrics_recovered(f).ElectrodeRejectionRecommendationsMetOrExceededThreshold = EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold;
    CleanedMetrics_recovered(f).HighProportionExcludedAsExtremeOutlier                  = EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier;
    CleanedMetrics_recovered(f).NoBlinksDetected                                        = EEG.RELAX_issues_to_check.NoBlinksDetected;
    fprintf('RELAX issues_to_check metrics saved to file.\n');
end

fprintf('Recovery complete. Total entries: %d\n', length(CleanedMetrics_recovered));
save(fullfile(savepath, 'CleanedMetrics_recovered.mat'), 'CleanedMetrics_recovered');
fprintf('Saved to CleanedMetrics_recovered.mat\n');
