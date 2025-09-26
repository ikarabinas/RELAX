
function save_ICA_topo(EEG, RELAX_cfg)
    % Save ICA topoplot images to:
    savedir = RELAX_cfg.OutputPath;
    ica_subdir = fullfile(savedir, 'RELAX_ICA_topoplots');
    if ~exist(ica_subdir, 'dir')
        mkdir(ica_subdir);
    end

    % Save in batches (e.g., 20 components per page)
    comps_per_page = 20;
    %total_comps = size(EEG.icawinv, 2);
    wanted_comps = 60;
    num_pages = ceil(wanted_comps / comps_per_page);
    display(num_pages)

    % Print file info before plotting
    fprintf('Plotting ICA components for: %s\n', EEG.setname);
    fprintf('Number of ICA components: %d\n', size(EEG.icawinv, 2));

    % Debug info
    fprintf('Size of icawinv: %s\n', mat2str(size(EEG.icawinv)));
    fprintf('Size of icaweights: %s\n', mat2str(size(EEG.icaweights)));
    fprintf('Size of icasphere: %s\n', mat2str(size(EEG.icasphere)));

    for page = 1:num_pages
        start_comp = (page-1) * comps_per_page + 1;
        end_comp = min(page * comps_per_page, wanted_comps);
        comp_range = start_comp:end_comp;
        display(start_comp);
        display(end_comp);
        display(comp_range);

        figure('visible', 'off');
        pop_topoplot(EEG, 0, comp_range, sprintf('%s - Components %d-%d', EEG.setname, start_comp, end_comp), ...
        [4 5], 0, 'electrodes', 'off');
        %pop_viewprops(EEG,0, 1:wanted_comps);

        filename = fullfile(ica_subdir, sprintf('%s_ICA_Page%d.png', EEG.setname, page));
        print(gcf, filename, '-dpng', '-r300');
        close(gcf);
    end
end
