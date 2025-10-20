% Concatenate sequential resting state EEG recordings and convert .mff files to .set format for use with RELAX preprocessing pipeline

N_CPU = 16;

EEG_DATA_DIR = '/athena/grosenicklab/store/tms_eeg/mdd_dlpfc';
RAW_DIR = fullfile(EEG_DATA_DIR, 'subject21_m275_dlpfc_59/m275_dlpfc_day2');
SAVE_DIR = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_cleaned/m275_dlpfc';
[~, ppt_target_day, ~] = fileparts(RAW_DIR);  % extract subject identifier from file name

if ~exist(SAVE_DIR, 'dir')
    mkdir(SAVE_DIR);
end

% Initialize EEGLAB
eeglab;

% Initialize variables to store datasets
EEG_pre1 = [];
EEG_pre2 = [];
EEG_post1 = [];
EEG_post2 = [];

% If saving data in double precision format for later use is beneficial, uncomment this line and the others related to double precision format
%pop_editoptions('option_single', 0);        % <- disable automatic conversion to single precision on saving

% Iterate over each resting state recording
raw_files = dir(fullfile(RAW_DIR, '**/*_reststate*.mff'));
for k = 1:numel(raw_files)
    file = raw_files(k);
    filepath = fullfile(file.folder, file.name);
    
    EEG = pop_mffimport(filepath);
    EEG = eeg_checkset(EEG);
    
    % Store the pre datasets
    if contains(file.name, 'reststate1_pre')
        EEG_pre1 = EEG;
        fprintf('Loaded %s reststate1_pre into memory\n', ppt_target_day);
    elseif contains(file.name, 'reststate2_pre')
        EEG_pre2 = EEG;
        fprintf('Loaded %s reststate2_pre into memory\n', ppt_target_day);
    elseif contains(file.name, 'reststate1_post')
        EEG_post1 = EEG;
        fprintf('Loaded %s reststate1_post into memory\n', ppt_target_day);
    elseif contains(file.name, 'reststate2_post')
        EEG_post2 = EEG;
        fprintf('Loaded %s reststate2_post into memory\n', ppt_target_day);
    end
end

% Concatenate and save the pre files
if ~isempty(EEG_pre1) && ~isempty(EEG_pre2)
    EEG_merged_pre = pop_mergeset(EEG_pre1, EEG_pre2);
    pre_filename = sprintf('%s_reststate_pre.set', ppt_target_day);
    EEG_merged_pre.setname = pre_filename;
    %EEG_merged_pre.data = double(EEG_merged_pre.data);  % convert to double precision
    EEG_merged_pre = eeg_checkset(EEG_merged_pre, 'makeur');  % recreate urevent structure
    EEG_merged_pre = pop_saveset(EEG_merged_pre, 'filename', pre_filename, 'filepath', SAVE_DIR);
    fprintf('Concatenated pre files saved with %d time points\n', EEG_merged_pre.pnts);
elseif ~isempty(EEG_pre1)
    pre_filename = sprintf('%s_reststate_pre.set', ppt_target_day);
    EEG_pre1.setname = pre_filename;
    %EEG_pre1.data = double(EEG_pre1.data);
    EEG_pre1 = pop_saveset(EEG_pre1, 'filename', pre_filename, 'filepath', SAVE_DIR);
    fprintf('Only 1 pre-treatment EEG file. Saved %s as .set file\n', pre_filename);
else
    warning('No pre-treatment EEG files found for participant %s', ppt_target_day);
end

% Concatenate and save post files
if ~isempty(EEG_post1) && ~isempty(EEG_post2)
    EEG_merged_post = pop_mergeset(EEG_post1, EEG_post2);
    post_filename = sprintf('%s_reststate_post.set', ppt_target_day);
    EEG_merged_post.setname = post_filename;
    %EEG_merged_post.data = double(EEG_merged_post.data);
    EEG_merged_post = eeg_checkset(EEG_merged_post, 'makeur');
    EEG_merged_post = pop_saveset(EEG_merged_post, 'filename', post_filename, 'filepath', SAVE_DIR);
    fprintf('Concatenated post files saved with %d time points\n', EEG_merged_post.pnts);
elseif ~isempty(EEG_post1)
    post_filename = sprintf('%s_reststate_post.set', ppt_target_day);
    EEG_post1.setname = post_filename;
    EEG_post1 = pop_saveset(EEG_post1, 'filename', post_filename, 'filepath', SAVE_DIR);
    fprintf('Only 1 post-treatment EEG file. Saved %s as .set file\n', post_filename);
else
    warning('No post-treatment EEG files found for participant %s', ppt_target_day);
end


eeglab redraw;
