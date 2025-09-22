% Save ICA topoplot images to:
savedir = 'home/imk2003/Documents/RELAX_ICA_topo/';

% Save in batches (e.g., 20 components per page)
comps_per_page = 20;
total_comps = size(EEG.icawinv, 2);
num_pages = ceil(total_comps / comps_per_page);

% Print file info before plotting
fprintf('Plotting ICA components for: %s\n', EEG.setname);
fprintf('Number of ICA components: %d\n', size(EEG.icawinv, 2));

% Debug info
fprintf('Size of icawinv: %s\n', mat2str(size(EEG.icawinv)));
fprintf('Size of icaweights: %s\n', mat2str(size(EEG.icaweights)));
fprintf('Size of icasphere: %s\n', mat2str(size(EEG.icasphere)));

for page = 1:num_pages
    start_comp = (page-1) * comps_per_page + 1;
    end_comp = min(page * comps_per_page, total_comps);
    comp_range = start_comp:end_comp;  % Use this variable!
    
    figure('visible', 'off');
    pop_topoplot(EEG, 0, comp_range, sprintf('%s - Components %d-%d', EEG.setname, start_comp, end_comp), ...
    [4 5], 0, 'electrodes', 'on');

    filename = sprintf('%s%s_ICA_Page%d.png', savedir, EEG.setname, page);
    print(gcf, filename, '-dpng', '-r300');
    close(gcf);
end