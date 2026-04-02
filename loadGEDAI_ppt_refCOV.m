function refCOV = loadGEDAI_ppt_refCOV(RELAX_cfg)
% Load a matching refCOV matrix based on the EEG filename in RELAX_cfg
% Inputs:
%   RELAX_cfg - struct with fields:
%       .filename         - EEG filename (e.g. 'm155_dlpfc_day3_reststate1_post_20230315_123516.set')
%       .GEDAI_refCOVPath - path to directory containing refCOV .mat files
% Outputs:
%   refCOV - matched covariance matrix

% Build prefix from first 2 underscore-delimited parts of EEG filename
% to match EEG file with ppt-specific refCOV
[~, eeg_filename_noext, ~] = fileparts(RELAX_cfg.filename);
eeg_parts = strsplit(eeg_filename_noext, '_');
eeg_prefix = strjoin(eeg_parts(1:2), '_');

% List all .mat files in refCOV directory
dirList = dir(fullfile(RELAX_cfg.GEDAI_refCOVPath, '*.mat'));
if isempty(dirList)
    error('No .mat files found in %s', RELAX_cfg.GEDAI_refCOVPath);
end

% Find a matching refCOV file. There should be one per participant.
matchIdx = find(startsWith({dirList.name}, eeg_prefix));
if isempty(matchIdx)
    error('No matching refCOV found for prefix: %s', eeg_prefix);
elseif length(matchIdx) > 1
    error('Multiple matching refCOV files found for prefix: %s', eeg_prefix);
end

% Load and return refCOV
refCOV_path = fullfile(RELAX_cfg.GEDAI_refCOVPath, dirList(matchIdx).name);
disp(['Loading custom ppt refCOV from: ' refCOV_path]);
data = load(refCOV_path);
refCOV = data.refCOV;

end
