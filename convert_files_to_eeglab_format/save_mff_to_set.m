% Save resting state EEG .mff files in .set format for use with RELAX preprocessing pipeline
% Run this processing in batch by calling through run_set_conversion.py

% Iterate over each resting state recording
function save_mff_to_set(raw_dir, save_dir)
    
    if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    end

    % Initialize EEGLAB
    eeglab;
    
    raw_files = dir(fullfile(raw_dir, '**/*_reststate*.mff'));
    for k = 1:numel(raw_files)
        file = raw_files(k);
        filepath = fullfile(file.folder, file.name);
        
        EEG = pop_mffimport(filepath);
        EEG = eeg_checkset(EEG);              % sanity check
        
        % Convert data to double precision
        %EEG.data = double(EEG.data);
        %EEG = pop_editoptions(EEG, 'option_savetwofiles', 0, 'option_saveversion6', 1);
        
        % Save as .set file with suffix
        savefile = strrep(file.name, '.mff', '.set');
        fprintf('%s', filepath)
        %EEG = pop_saveset(EEG, 'filename', savefile, 'filepath', save_dir, 'savemode', 'twofiles');  % uncomment if .fdt + .set saving is preferred
	EEG = pop_saveset(EEG, 'filename', savefile, 'filepath', save_dir);  % save as one .set file by default
    
        % Print the full save paths
        full_savepath = fullfile(save_dir, savefile);
        fprintf('Saving file to: %s\n', full_savepath);
    
    end
    
    eeglab redraw;
end
