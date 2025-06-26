#!/bin/bash
#SBATCH --partition general
#SBATCH --ntasks 1
#SBATCH --array 1-206
#SBATCH --mem-per-cpu 5G
#SBATCH --constraint epyc128
#SBATCH --output=/home/til19015/JobOutput/tmp_output_DarkMask/slurm-%A_%a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_DarkMask/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/Analyze/DarkPixel

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX

export job njobs
matlab -nosplash -singleCompThread -r loop_DarkPixel_row\($job,$njobs\)\;exit\;
