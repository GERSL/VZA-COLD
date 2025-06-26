#!/bin/bash
#SBATCH -J Grid
#SBATCH -p general
#SBATCH --ntasks 1
#SBATCH --array 1-220
#SBATCH --nodes=1
#SBATCH --mem-per-cpu 5G
#SBATCH --constraint epyc128
#SBATCH --output=/home/til19015/JobOutput/tmp_output_grid/slurm-%A_%a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_grid/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/Analyze/Gridding

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX
export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_grid\($job,$njobs\)\;exit\;
