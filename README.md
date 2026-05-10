# RELAX-GEDAI
This fork of the Reduction of Electroencephalographic Artifacts (RELAX) incorporates the Generalized Eigenvalue De-Artifacting Instrument (GEDAI) method of lead field filtering for denoising of resting state EEG. GEDAI as implemented here can be run with with custom reference covariance matrices (refCOV) computed for participants based on their individual forward model (if anatomical MR images are available).

Specify denoising parameters and the paths to your EEG data and custom refCOV matrices in RELAX_SET_PARAMETERS_AND_RUN.m. We have had success cleaning relatively noisy files collected in an outpatient clinical setting with a combined approach: RELAX_cfg.Perform_targeted_wICA=1 and RELAX_cfg.Run_GEDAI=1. These settings will apply temporal filtering, bad channel and bad segment rejection, targeted wavelet-enhaced ICA (twICA), then follow with GEDAI.

Custom refCOV matrices can be generated from individual participant lead fields.
See ikarabinas/GEDAI-master for code and example usage of GEDAI denoising with custom participant-specific refCOVs.  

## 📧 Contact
For questions about this fork: Isabella Karabinas - imk2003@med.cornell.edu  

RELAX repository: NeilwBailey/RELAX  
GEDAI repository: neurotuning/GEDAI-master  
