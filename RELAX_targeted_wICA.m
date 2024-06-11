%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2024) Neil Bailey

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

%% RELAX_targeted_wICA:
% Bailey N.W., Hill A.T., Godfrey K., Perera M.P.N., Rogasch N.C., Fitzgibbon B.M., & Fitzgerald P.B. (2024). EEG is better when cleaning effectively targets artifacts. BioRxiv. https://doi.org/10.1101/2024.06.06.597688

function [EEG] = RELAX_targeted_wICA(EEG,RELAX_cfg)

    % set some defaults if not specified
    RELAX_cfg.ms_per_sample=(1000/EEG.srate); % to determine number of ms per sample (used in blink targetting)
    fastica_symm_Didnt_Converge=[0 0 0]; % to track whether fastica_symm doesn't converge
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'Clean_other_comps')==0
            RELAX_cfg.Clean_other_comps='no';
        end
        if isfield(RELAX_cfg, 'ICA_method')==0
            RELAX_cfg.ICA_method='picard';
        end
        if isfield(RELAX_cfg, 'MuscleSlopeThreshold')==0
            RELAX_cfg.MuscleSlopeThreshold=-0.31; 
        end
        if isfield(RELAX_cfg, 'Report_all_wICA_info')==0
            RELAX_cfg.Report_all_wICA_info='no'; 
        end
        if isfield(RELAX_cfg, 'ICLabel_thresholds')==0
            RELAX_cfg.ICLabel_thresholds=[0 0 0 0 0 0 0];
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.Clean_other_comps='no';
        RELAX_cfg.ICA_method='picard';
        RELAX_cfg.MuscleSlopeThreshold=-0.31; 
        RELAX_cfg.Report_all_wICA_info='no'; 
        RELAX_cfg.ICLabel_thresholds=[0 0 0 0 0 0 0];
    end

    % run ICA:
    if strcmp(RELAX_cfg.ICA_method,'extended_infomax_ICA')
        [EEG_with_ICA, ~] = pop_runica_nwb(EEG, 'extended',1,'interupt','on');
        W = EEG_with_ICA.icaweights*EEG_with_ICA.icasphere;
        A = inv(W);
        EEG_with_ICA = eeg_checkset(EEG_with_ICA, 'ica'); 
        if isempty(EEG_with_ICA.icaact)==1
            EEG_with_ICA.icaact = (EEG_with_ICA.icaweights*EEG_with_ICA.icasphere)*EEG_with_ICA.data(EEG_with_ICA.icachansind,:);      
            EEG_with_ICA.icaact = reshape( EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), EEG_with_ICA.pnts, EEG_with_ICA.trials);
        end
        Component=reshape(EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'cudaica')
        [EEG_with_ICA, ~] = pop_runica_nwb(EEG, 'cudaica', 'extended',1); 
        W = EEG_with_ICA.icaweights*EEG_with_ICA.icasphere;
        A = inv(W);
        EEG_with_ICA = eeg_checkset(EEG_with_ICA, 'ica'); 
        if isempty(EEG_with_ICA.icaact)==1
            EEG_with_ICA.icaact = (EEG_with_ICA.icaweights*EEG_with_ICA.icasphere)*EEG_with_ICA.data(EEG_with_ICA.icachansind,:);      
            EEG_with_ICA.icaact = reshape( EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), EEG_with_ICA.pnts, EEG_with_ICA.trials);
        end
        Component=reshape(EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'fastica_symm')
        % The following lines repeat fastica_symm up to 3 times in
        % the case of non-convergence, then switches to fastica_defl to
        % ensure ICA convergence (as cleaning as adversely affected by
        % non-convergence issues).
         [EEG_with_ICA, ~, NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
         fastica_symm_Didnt_Converge(1,1)=NonConvergence;
         if NonConvergence==1
             [EEG_with_ICA, ~, NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
             fastica_symm_Didnt_Converge(1,2)=NonConvergence;
         end
         if NonConvergence==1
             [EEG_with_ICA, ~, NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
             fastica_symm_Didnt_Converge(1,3)=NonConvergence;
         end
         if NonConvergence==1
             EEG_with_ICA = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'defl', 'g', 'tanh', 'stabilization', 'on');
         end
         W = EEG_with_ICA.icaweights*EEG_with_ICA.icasphere;
         A = inv(W);
         EEG_with_ICA = eeg_checkset(EEG_with_ICA, 'ica'); 
        if isempty(EEG_with_ICA.icaact)==1
            EEG_with_ICA.icaact = (EEG_with_ICA.icaweights*EEG_with_ICA.icasphere)*EEG_with_ICA.data(EEG_with_ICA.icachansind,:);      
            EEG_with_ICA.icaact = reshape( EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), EEG_with_ICA.pnts, EEG_with_ICA.trials);
        end
         Component=reshape(EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'fastica_defl')
         EEG_with_ICA = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'defl', 'g', 'tanh', 'stabilization', 'on');
         W = EEG_with_ICA.icaweights*EEG_with_ICA.icasphere;
         A = inv(W);
         EEG_with_ICA = eeg_checkset(EEG_with_ICA, 'ica'); 
        if isempty(EEG_with_ICA.icaact)==1
            EEG_with_ICA.icaact = (EEG_with_ICA.icaweights*EEG_with_ICA.icasphere)*EEG_with_ICA.data(EEG_with_ICA.icachansind,:);      
            EEG_with_ICA.icaact = reshape( EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), EEG_with_ICA.pnts, EEG_with_ICA.trials);
        end
         Component=reshape(EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'amica')
        EEG_with_ICA=EEG;
        % You'll need to install AMICA first, and in the folder that you
        % specify in the line below (with no spaces in any part of the folder or subfolders):
        % You can download AMICA via EEGLAB
        amica_file = which('runamica15');
        if ~exist(amica_file)
            disp('!! AMICA directory not found, please ensure you have AMICA installed !!');
        end
        [filepath,~,~] = fileparts(amica_file);
        cd(filepath);
        % define parameters
        numprocs = 1;       % # of nodes (default = 1)
        max_threads = 4;    % # of threads per node
        num_models = 1;     % # of models of mixture ICA
        max_iter = 2000;    % max number of learning steps
        mkdir([filepath filesep 'AMICAtmp']);
        outdir = [filepath filesep 'AMICAtmp' filesep];
        % Run AMICA:    
        [EEG_with_ICA.icaweights, EEG_with_ICA.icasphere, ~] = runamica15(EEG_with_ICA.data, 'num_chans', EEG.nbchan, 'num_models',num_models,'outdir',outdir,'numprocs', numprocs, 'max_threads', max_threads, 'max_iter',max_iter,'pcakeep', EEG.nbchan, 'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
        W = EEG_with_ICA.icaweights*EEG_with_ICA.icasphere;
        A = inv(W);
        EEG_with_ICA = eeg_checkset(EEG_with_ICA, 'ica'); 
        if isempty(EEG_with_ICA.icaact)==1
            EEG_with_ICA.icaact = (EEG_with_ICA.icaweights*EEG_with_ICA.icasphere)*EEG_with_ICA.data(EEG_with_ICA.icachansind,:);      
            EEG_with_ICA.icaact = reshape( EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), EEG_with_ICA.pnts, EEG_with_ICA.trials);
        end
        Component=reshape(EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'picard') % Run PICARD-O using default settings
        [EEG_with_ICA, ~] = pop_runica_nwb(EEG, 'picard', 'mode','ortho','tol',1e-6,'maxiter',500); % run picard
        W = EEG_with_ICA.icaweights*EEG_with_ICA.icasphere;
        A = inv(W);
        EEG_with_ICA = eeg_checkset(EEG_with_ICA, 'ica'); 
        if isempty(EEG_with_ICA.icaact)==1
            EEG_with_ICA.icaact = (EEG_with_ICA.icaweights*EEG_with_ICA.icasphere)*EEG_with_ICA.data(EEG_with_ICA.icachansind,:);      
            EEG_with_ICA.icaact = reshape( EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), EEG_with_ICA.pnts, EEG_with_ICA.trials);
        end
        Component=reshape(EEG_with_ICA.icaact, size(EEG_with_ICA.icaact,1), []);
    end
    % Use ICLabel to identify artifactual components, so that wICA can be
    % performed on them only:
    EEG_with_ICA = iclabel(EEG_with_ICA);
    IC_classifications=EEG_with_ICA.etc.ic_classification.ICLabel.classifications; % allows user to set thresholds for classification confidence before considered an artifact
    IC_classifications(IC_classifications<RELAX_cfg.ICLabel_thresholds)=0;
    [~, I]=max(IC_classifications, [], 2);
    if strcmp(RELAX_cfg.Clean_other_comps,'no')==1
        ICsMostLikelyNotBrain=(I==2 | I ==3)'; 
    elseif strcmp(RELAX_cfg.Clean_other_comps,'yes')==1
        ICsMostLikelyNotBrain=(I>1)'; 
    end
    ICsMostLikelyEye=(I==3)';

    options.muscleFreqIn=[7,70];
    options.Freq_to_compute = [1,100];

    % Calculate pwelch to enable detection of log-freq log-power slopes
    % indicative of muscle activity
    % Resize EEG.icaact if required
    if size(EEG_with_ICA.icaact,3) > 0
        eegData = reshape(EEG_with_ICA.icaact,size(EEG_with_ICA.icaact,1),[]);
    else
        eegData = EEG_with_ICA.icaact;
    end
    [pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),EEG_with_ICA.srate);
    FFTout = pxx';
    fp = fp';
    
    % Calculate FFT bins
    freq=options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
    fftBins = zeros(size(FFTout,1),size(freq,2)); %preallocate
    for a=1:size(freq,2)
        [~, index1]=min(abs(fp-((freq(1,a)-0.25))));
        [~, index2]=min(abs(fp-((freq(1,a)+0.25))));
        fftBins(:,a)=mean(FFTout(:,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
    end

    %% better muscle comp_number identification:
    comps=size(EEG_with_ICA.icaact,1);
    options.muscleFreqEx=[RELAX_cfg.LineNoiseFrequency-2 RELAX_cfg.LineNoiseFrequency+2];
    for compNum =1:comps
        % Define frequencies to include in the analysis
        if ~isempty(options.muscleFreqIn)
            [~,fin1] = min(abs(options.muscleFreqIn(1) - freq));
            [~,fin2] = min(abs(options.muscleFreqIn(2) - freq));
            freqHz = freq(1,fin1:fin2);
            freqPow = fftBins(compNum,fin1:fin2);
        else
            freqHz = freq;
            freqPow = fftBins(compNum,:);
        end
        % Define frequencies to exclude from fit
        if ~isempty(options.muscleFreqEx)
            [~,fex1] = min(abs(options.muscleFreqEx(1) - freqHz));
            [~,fex2] = min(abs(options.muscleFreqEx(2) - freqHz));
            freqHz(fex1:fex2) = [];
            freqPow(fex1:fex2) = [];
        end
        % Fit linear regression to log-log data
        p = polyfit(log(freqHz),log(freqPow),1);
        % Store the slope
        muscleRatio(compNum) = p(1);
    end
    muscle_ICs=muscleRatio>=RELAX_cfg.MuscleSlopeThreshold;
    ICsMostLikelyMuscle=(muscle_ICs==1);

    % use icablinkmetrics to double check for blink components that ICLabel
    % might have missed:
    if exist('icablinkmetrics', 'file') == 2
        EEG_all_electrodes = pop_interp(EEG, EEG.allchan, 'spherical');
        electrode_labels={EEG_all_electrodes.chanlocs.labels};
        Blink_Electrode_location=find(strcmpi(electrode_labels,RELAX_cfg.BlinkElectrodes(1,1)));
        try
            icablinkmetricsout = icablinkmetrics(EEG_with_ICA, 'ArtifactChannel', EEG_all_electrodes.data(Blink_Electrode_location,:), 'Alpha', 0.001, 'VisualizeData', 'False');
            if any(icablinkmetricsout.identifiedcomponents>0)
                ICsMostLikelyNotBrain(1,icablinkmetricsout.identifiedcomponents)=1;
                ICsMostLikelyEye(1,icablinkmetricsout.identifiedcomponents)=1;
            end
        catch
            warning('icablinkmetrics did not detect any blinks');
        end
    end

    check_padding_required = mod(size(Component,2),2^5);
    if check_padding_required ~=0
        padding = zeros(1,(2^5)-check_padding_required);
    else
        padding = [];
    end

    %% perform wavelet thresholding on eye movements (and also other components if selected), identified by ICLabel:
    disp('Using targeted approach to clean artifacts');
    for comp_number = 1:size(Component,1)
        if ICsMostLikelyNotBrain(comp_number)==1 % wavelet enhance only on artifacts identified by ICLabel
            if ~isempty(padding)
                padded_comp = [Component(comp_number,:),padding]; % pad the component with zeros if required
            else
                padded_comp = Component(comp_number,:);
            end
            [wavelet_threshold,threshold_type,~] = ddencmp('den','wv',padded_comp); % automatically obtain wavelet enhancement threshold
            if ICsMostLikelyEye(comp_number)==1
                wavelet_threshold = wavelet_threshold*2; % increase threshold for blink components based on optimal results in our informal testing
            else
                wavelet_threshold = wavelet_threshold*1;
            end
            wavelet_transform = swt(padded_comp,5,'coif5'); % apply stationary wavelet transform to each component to reduce neural contribution to component
            thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold); % remove negligible values by applying thresholding
            artifact_comp(comp_number,:) = iswt(thresholded_wavelet_transform,'coif5'); % use inverse wavelet transform to obtained the wavelet transformed component
            clear thresholded_wavelet_transform padded_comp wavelet_threshold threshold_type wavelet_transform 
        end
    end
    % pad non-artifact components with 0s in the same way that the artifact components were padded:
    if sum(ICsMostLikelyNotBrain)==0
        artifact_comp(1,:)=zeros(1,size(EEG.data,2));
        artifact_comp = [artifact_comp(:,:),padding]; % pad with zeros
    end
    for comp_number = 1:size(Component,1)
        if ICsMostLikelyNotBrain(comp_number)==0
            artifact_comp(comp_number,:)=zeros(1,size(artifact_comp,2));
        end
    end
    % remove padding
    if ~isempty(padding)
        artifact_comp = artifact_comp(:,1:end-numel(padding));
    end

    %% Restrict wICA cleaning of blink components to just blink periods:
    moving_mean_length=round(200/RELAX_cfg.ms_per_sample);
    blink_length_threshold=round(100/RELAX_cfg.ms_per_sample);
    clear M;
    for comp_number=1:size(Component,1)
        if ICsMostLikelyEye(comp_number)==1
            [z1, p1] = butter(2, [0.5 25]./(EEG.srate/2), 'bandpass');
            dataIn=Component(comp_number,:)';
            dataIn=double(dataIn);
            dataFilt1 = filtfilt(z1,p1,dataIn);
            IC_filtered = dataFilt1';
            [blink_periods,~,~]=isoutlier(IC_filtered,'median',ThresholdFactor=2);
            ix_blinkstart=find(diff(blink_periods)==1)+1;  % indices where BlinkIndexMetric goes from 0 to 1
            ix_blinkend=find(diff(blink_periods)==-1);  % indices where BlinkIndexMetric goes from 1 to 0

            [EEG_with_ICA, ~] = RELAX_blinks_IQR_method(EEG_with_ICA, EEG_with_ICA, RELAX_cfg); % use an IQR threshold method to detect and mark blinks
            blink_periods(EEG_with_ICA.RELAX.eyeblinkmask==1)=1;

            if ~isempty(ix_blinkstart)
                if ix_blinkend(1,1)<ix_blinkstart(1,1); ix_blinkend(:,1)=[]; end % if the first downshift occurs before the upshift, remove the first value in end
                if ix_blinkend(1,size(ix_blinkend,2))<ix_blinkstart(1,size(ix_blinkstart,2)); ix_blinkstart(:,size(ix_blinkstart,2))=[];end % if the last upshift occurs after the last downshift, remove the last value in start
                BlinkThresholdExceededLength=ix_blinkend-ix_blinkstart; % length of consecutive samples where blink threshold was exceeded
                BlinkRunIndex = find(BlinkThresholdExceededLength<round(blink_length_threshold/RELAX_cfg.ms_per_sample)); % find locations where blink threshold was not exceeded by more than X ms
                % find latency of the max voltage within each period where the blink
                % threshold was exceeded:
                if size(BlinkRunIndex,2)>0
                    for x=1:size(BlinkRunIndex,2)
                        o=ix_blinkstart(BlinkRunIndex(x));
                        c=ix_blinkend(BlinkRunIndex(x));
                        if c-o<round(blink_length_threshold/RELAX_cfg.ms_per_sample)
                            blink_periods(1,o:c)=0;
                        end
                    end
                end
            end
            padded_blink_periods=double(blink_periods);
            for c=flip(1:size(padded_blink_periods,2)-(moving_mean_length+1))
                if padded_blink_periods(1,c)==1 
                    padded_blink_periods(1,c:c+moving_mean_length)=1;
                end
            end
            for c=(moving_mean_length+1):size(padded_blink_periods,2)
                if padded_blink_periods(1,c)==1 
                    padded_blink_periods(1,c-moving_mean_length:c)=1;
                end
            end
            M(comp_number,:) = movmean(padded_blink_periods,[moving_mean_length moving_mean_length]);
            artifact_comp(comp_number,:)=(artifact_comp(comp_number,:).*M(comp_number,:));
        end
    end
    
    % Obtain muscle artifact for subtraction by high pass filtering data instead of wICA:
    for comp_number=1:size(Component)
        if ICsMostLikelyMuscle(comp_number)==1
            [z1, p1] = butter(2, 15./(EEG.srate/2), 'high');
            dataIn=Component(comp_number,:)';
            dataIn=double(dataIn);
            dataFilt1 = filtfilt(z1,p1,dataIn);
            artifact_comp(comp_number,:) = dataFilt1';
        end
    end

    % Remove artifact and reconstruct data:
    artifacts_in_EEG = A*artifact_comp;
    %reshape EEG signal from EEGlab format to channelsxsamples format
    Original_EEG=reshape(EEG.data, size(EEG.data,1), []);
    %subtract out wavelet artifact signal from EEG signal
    Cleaned_EEG=Original_EEG-artifacts_in_EEG;
    EEG.data = Cleaned_EEG;
        
    EEG.RELAXProcessing_wICA.fastica_symm_Didnt_Converge=fastica_symm_Didnt_Converge; % Tracks whether fastica_symm showed convergence issues (1) or not (0), and how many non-convergences. If 3 non-convergences, then fastica_defl was implemented.
    
    % Check if data might have been too short for effective ICA, using Makoto's rule
    % of thumb that ICA requires data length of ((number of channels)^2)*30
    % if data were sampled at 250 Hz (assuming that higher sampling
    % rates require the same time duration of data as low sampling rates,
    % so 1000Hz sampling rates require ((number of channels)^2)*120)
    % (https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code)
    ms_per_sample=(1000/EEG.srate);
    if ((EEG.nbchan^2)*(120/ms_per_sample))>EEG.pnts
        EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA='yes';
    else
        EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA='no';
    end
    
    if strcmp (EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA,'yes')
        warning('Data may have been shorter than recommended for effective ICA decomposition')
    end
    
    EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA=mean(ICsMostLikelyNotBrain);
    
    if strcmp(RELAX_cfg.Report_all_wICA_info,'yes')
    
        EEG.RELAXProcessing_wICA.ProportionICs_was_Brain=sum(I==1)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_wICA.ProportionICs_was_Muscle=sum(I==2)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_wICA.ProportionICs_was_Eye=sum(I==3)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_wICA.ProportionICs_was_Heart=sum(I==4)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_wICA.ProportionICs_was_LineNoise=sum(I==5)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_wICA.ProportionICs_was_ChannelNoise=sum(I==6)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_wICA.ProportionICs_was_Other=sum(I==7)/size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1);    

        ICsMostLikelyBrain=(I==1)';
        ICsMostLikelyMuscle=(I==2)';
        ICsMostLikelyEye=(I==3)';
        ICsMostLikelyHeart=(I==4)';
        ICsMostLikelyLineNoise=(I==5)';
        ICsMostLikelyChannelNoise=(I==6)';
        ICsMostLikelyOther=(I==7)';

        for x=1:size(EEG_with_ICA.etc.ic_classification.ICLabel.classifications,1)
            [~, varianceWav(x)] =compvar(EEG_with_ICA.data, EEG_with_ICA.icaact, EEG_with_ICA.icawinv, x);
        end

        BrainVariance=sum(abs(varianceWav(ICsMostLikelyBrain)));
        ArtifactVariance=sum(abs(varianceWav(~ICsMostLikelyBrain)));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_BrainICs=(BrainVariance/(BrainVariance+ArtifactVariance));

        MuscleVariance=sum(abs(varianceWav(ICsMostLikelyMuscle)));
        EyeVariance=sum(abs(varianceWav(ICsMostLikelyEye)));
        HeartVariance=sum(abs(varianceWav(ICsMostLikelyHeart)));
        LineNoiseVariance=sum(abs(varianceWav(ICsMostLikelyLineNoise)));
        ChannelNoiseVariance=sum(abs(varianceWav(ICsMostLikelyChannelNoise)));
        OtherVariance=sum(abs(varianceWav(ICsMostLikelyOther)));

        EEG.RELAXProcessing_wICA.ProportionVariance_was_MuscleICs=(MuscleVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_EyeICs=(EyeVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_HeartICs=(HeartVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_LineNoiseICs=(LineNoiseVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_ChannelNoiseICs=(ChannelNoiseVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_OtherICs=(OtherVariance/(BrainVariance+ArtifactVariance));
    
    end
end
