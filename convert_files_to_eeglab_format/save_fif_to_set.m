% Save resting state EEG .fif files in .set format for use with RELAX preprocessing pipeline

N_CPU = 16;

RAW_DIR = '/athena/grosenicklab/store/tms_eeg/bart/bart_reststate_noise_by_room_55/';
SAVE_DIR = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_to_clean/bart/bart_intersession_room109';

if ~exist(SAVE_DIR, 'dir')
    mkdir(SAVE_DIR);
end

% Initialize EEGLAB
eeglab;

% Iterate over each resting state recording
raw_files = dir(fullfile(RAW_DIR, '**/*_reststate*.fif'));
for k = 1:numel(raw_files)
    file = raw_files(k);
    filepath = fullfile(file.folder, file.name);
    
    EEG = pop_fileio(filepath);

    
    EEG = eeg_checkset(EEG);              % sanity check
    
    % Convert data to double precision
    %EEG.data = double(EEG.data);
    %EEG = pop_editoptions(EEG, 'option_savetwofiles', 0, 'option_saveversion6', 1);

    % Save as .set file with suffix
    savefile = strrep(file.name, '.fif', '.set');
    EEG = pop_saveset(EEG, 'filename', savefile, 'filepath', SAVE_DIR);

    % Print the full save paths
    full_savepath = fullfile(SAVE_DIR, savefile);
    fprintf('Saving file to: %s\n', full_savepath);

end

eeglab redraw;
