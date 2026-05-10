#!/bin/bash
#SBATCH --job-name=relax-gedai
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=5-00:00:00
#SBATCH --output=relax-gedai_log_%j.txt
#SBATCH --error=relax-gedai_error_%j.txt

# Load MATLAB module
module load matlab/R2023a

# Navigate to your script directory
cd "/home/imk2003/Documents/MATLAB/eeglab/plugins/RELAX"

# Run pipeline
matlab -nodisplay -nosplash -nodesktop -r "run('RELAX_SET_PARAMETERS_AND_RUN.m'); exit"
